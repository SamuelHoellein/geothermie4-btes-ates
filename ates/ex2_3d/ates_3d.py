#!/usr/bin/env python3
"""
ATES 3D demo for OpenGeoSys 6 (Hydro-Thermal process).

Geschichtetes Reservoir:  Cap Rock (oben) | Aquifer | Cap Rock (unten)
Ein Brunnenfilter als kleine Box im Aquifer (Single-Well-Anlage).
Zyklischer Lade-/Entlade-Betrieb über CurveScaled-Parameter:
  - Pressure-Equation: Volumetric Source ±Q/V_well auf Brunnenbox
  - Temperature: Dirichlet-BC auf Brunnenbox (curve-skaliert)
Der Lateralrand des Aquifers ist ein p=0-Outlet, damit das injizierte
Wasser entweichen kann.

Alle Modellgrößen sind im CONFIG-Block einstellbar.
Aufruf:
    python ates_3d.py            # Mesh + .prj + OGS-Lauf
    python ates_3d.py --no-run   # nur Setup, kein OGS
    python ates_3d.py --no-mesh  # nur .prj (Mesh muss existieren)
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import gmsh
import numpy as np

# ======================================================================
#  CONFIG  --  hier alles anpassen
# ======================================================================
CONFIG: dict = {
    "domain": {
        "size_x_m":   200.0,    # Ausdehnung in x [m]
        "size_y_m":   200.0,    # Ausdehnung in y [m]
        "z_base_m":     0.0,    # untere Modellgrenze (z-Koordinate)
    },
    "layers": {
        "caprock_bottom_thickness_m": 30.0,
        "aquifer_thickness_m":        30.0,
        "caprock_top_thickness_m":    20.0,
    },
    "wells": {
        # Single-Well-Anlage: ein aktiver Brunnen. Der Lateralrand des
        # Aquifers wird als p=0-Outlet gesetzt, damit das am Brunnen
        # injizierte Wasser entweichen kann.
        "hot_well_xy":   ( 0.0,  0.0),     # (x, y) Lage des Brunnens
        "screen_bottom_offset_m": 5.0,     # Abstand Filterunterkante vom Aquiferboden
        "screen_top_offset_m":    5.0,     # Abstand Filteroberkante von Aquiferdecke
        "screen_dx_m":             2.0,    # x-Ausdehnung des Filtervolumens
        "screen_dy_m":             2.0,    # y-Ausdehnung des Filtervolumens
        "screen_permeability_m2":  1.0e-9, # Filterkies (>> Aquifer)
    },
    # ------------------------------------------------------------------
    # REGIONALE GRUNDWASSERSTRÖMUNG (optional)
    # ------------------------------------------------------------------
    # Eine grossräumige Hintergrund-Strömung wird durch einen linearen
    # Druckgradient auf der Lateral-Aquifer-Fläche aufgeprägt
    # (Dirichlet-BC mit p(x, y) = ρ_f · g · i · cos(α) · x + ...).
    # Die Plume wird in Hauptströmungs-Richtung verschoben.
    # Setze enable = False, um einen "ruhenden" Aquifer zu simulieren
    # (Default p = 0 als reines Outlet, kein Gradient).
    # ------------------------------------------------------------------
    "regional_gw": {
        "enable":            False,
        "gradient_m_per_m":  1.0e-3,       # hydraulischer Gradient (dimensionslos)
        "direction_deg":     0.0,          # 0° = +x, 90° = +y, ...
    },
    "mesh": {
        "size_in_well_m":       0.8,       # Elementgröße im Filtervolumen
        "size_near_wells_m":    2.0,
        "size_far_m":          12.0,
        "well_size_radius_m":  10.0,       # bis zu diesem Radius: feinmaschig
        "well_size_radius_far_m": 35.0,    # ab dem: grobmaschig
    },
    "materials": {
        # MaterialID 0 = aquifer, 1 = caprock_top, 2 = caprock_bottom
        "aquifer": {
            "permeability_m2": 1.0e-12,
            "porosity":        0.25,
            "rho_s_kg_m3":  2650.0,
            "cp_s_J_kgK":   1000.0,
            "lambda_s_W_mK":   3.0,
        },
        "caprock_top": {
            "permeability_m2": 1.0e-18,
            "porosity":        0.05,
            "rho_s_kg_m3":  2700.0,
            "cp_s_J_kgK":    900.0,
            "lambda_s_W_mK":   2.0,
        },
        "caprock_bottom": {
            "permeability_m2": 1.0e-18,
            "porosity":        0.05,
            "rho_s_kg_m3":  2700.0,
            "cp_s_J_kgK":    900.0,
            "lambda_s_W_mK":   2.0,
        },
    },
    "fluid": {
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         293.15,
        "beta_1_per_K":    0.0,        # >0 aktiviert Linear-Dichte ρ(T)=ρ_ref(1-β(T-T_ref))
        "viscosity_Pa_s":  1.0e-3,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },
    "dispersion": {
        # Mechanische thermische Dispersivität – stabilisiert advektiv dominierten
        # Transport (hoher Péclet). Realistisch für poröse Medien α_L ≈ 0.1–10 m.
        "alpha_L_m": 1.0,
        "alpha_T_m": 0.1,
    },
    "initial": {
        "T_K":  283.15,
        "p_Pa": 0.0,
    },
    "operation": {
        "mass_flow_rate_kg_s": 0.5,    # Massenstrom am Brunnen (gesamt über alle Filterpunkte)
        "T_hot_K":  353.15,            # Vorlauftemperatur Beladung, 80 °C
        # Spezifische Speicherzahl der Phasen (1/Pa). Werte > 0 sind nötig,
        # damit die Druckgleichung transient ist – sonst löst die Punktquelle
        # zu Druck-Singularitäten auf.
        "fluid_storage_1_per_Pa": 4.5e-10,    # H2O Kompressibilität
        "solid_storage_1_per_Pa": 1.0e-10,    # Korngerüst
    },
    # ------------------------------------------------------------------
    # ZYKLEN – HIER FÜR STUDIERENDE
    # ------------------------------------------------------------------
    # Ein vollständiger Zyklus besteht aus 4 aufeinander folgenden Phasen:
    #   1) charge                    – Beladung (heiß injizieren)
    #   2) storage_after_charge      – Pause/Speicherung nach Beladung
    #   3) discharge                 – Förderung (Wasser entnehmen)
    #   4) storage_after_discharge   – Pause/Speicherung nach Förderung
    #
    # Periode T_Zyklus = charge + storage_after_charge + discharge + storage_after_discharge
    # Gesamt­simulations­zeit = n_cycles * T_Zyklus
    #
    # Beispiele:
    #   - Saisonal (1 Jahr/Zyklus): 91.25 / 91.25 / 91.25 / 91.25
    #   - Sommerladung 120 d, Winterförderung 120 d, sonst Pause: 120 / 60 / 120 / 60
    #   - Phase auf 0 setzen, um sie zu deaktivieren.
    # ------------------------------------------------------------------
    "cycles": {
        "n_cycles":                     3,     # Anzahl Wiederholungen des Zyklus
        "charge_days":                 90,     # Phase 1: Beladung (Tage)
        "storage_after_charge_days":    0,     # Phase 2: Pause nach Beladung (Tage)
        "discharge_days":              90,     # Phase 3: Förderung (Tage)
        "storage_after_discharge_days": 0,     # Phase 4: Pause nach Förderung (Tage)
        "ramp_days":                    1.0,   # Sanfte Übergangsrampe zwischen Phasen (Tage)
    },
    "time": {
        "dt_seconds":           86400.0,    # 1 Tag
        "output_every_n_steps": 1,
        "gravity":              False,       # True = z-Body-Force (-9.81)
    },
    "output": {
        "prefix":    "ates_3d",
        "out_dir":   "out",
        "variables": ["T", "p", "darcy_velocity"],
    },
    "solver": {
        "linear_tol":   1.0e-12,
        "linear_iter":  10000,
        "nonlinear_iter": 20,
        "rel_tol_T":   1.0e-4,
        "rel_tol_p":   1.0e-4,
    },
}

DAY = 86400.0
G   = 9.81

# ======================================================================
#  Mesh – gmsh
# ======================================================================
def build_mesh(cfg: dict, out_dir: Path) -> Path:
    """Schichtmodell mit gmsh: 3 Schichten + 1 Brunnenbox im Aquifer."""
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"

    Lx     = cfg["domain"]["size_x_m"]
    Ly     = cfg["domain"]["size_y_m"]
    z_base = cfg["domain"]["z_base_m"]

    t_cb = cfg["layers"]["caprock_bottom_thickness_m"]
    t_aq = cfg["layers"]["aquifer_thickness_m"]
    t_ct = cfg["layers"]["caprock_top_thickness_m"]

    z_aq_bot = z_base + t_cb
    z_aq_top = z_aq_bot + t_aq
    z_top    = z_aq_top + t_ct

    hw  = cfg["wells"]["hot_well_xy"]
    sob = cfg["wells"]["screen_bottom_offset_m"]
    sot = cfg["wells"]["screen_top_offset_m"]
    dx  = cfg["wells"]["screen_dx_m"]
    dy  = cfg["wells"]["screen_dy_m"]
    z_screen_bot = z_aq_bot + sob
    z_screen_top = z_aq_top - sot
    h_screen = z_screen_top - z_screen_bot

    s_in   = cfg["mesh"]["size_in_well_m"]
    s_near = cfg["mesh"]["size_near_wells_m"]
    s_far  = cfg["mesh"]["size_far_m"]
    r_near = cfg["mesh"]["well_size_radius_m"]
    r_far  = cfg["mesh"]["well_size_radius_far_m"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ates")

    x0, y0 = -Lx / 2.0, -Ly / 2.0
    box_cb = gmsh.model.occ.addBox(x0, y0, z_base,   Lx, Ly, t_cb)
    box_aq = gmsh.model.occ.addBox(x0, y0, z_aq_bot, Lx, Ly, t_aq)
    box_ct = gmsh.model.occ.addBox(x0, y0, z_aq_top, Lx, Ly, t_ct)
    box_hw = gmsh.model.occ.addBox(hw[0] - dx / 2.0, hw[1] - dy / 2.0,
                                   z_screen_bot, dx, dy, h_screen)

    # Alles fragmentieren (konforme Schnittflächen)
    gmsh.model.occ.fragment([(3, box_cb)],
                            [(3, box_aq), (3, box_ct), (3, box_hw)])
    gmsh.model.occ.synchronize()

    # Volumen klassifizieren
    vol_aq, vol_ct, vol_cb, vol_hw = [], [], [], []
    for dim, tag in gmsh.model.getEntities(3):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (zmin + zmax)
        ext_x = xmax - xmin
        if zc > z_aq_top + 1e-6:
            vol_ct.append(tag)
        elif zc < z_aq_bot - 1e-6:
            vol_cb.append(tag)
        else:
            # innerhalb Aquifer-z – evtl. Brunnenbox?
            xc = 0.5 * (xmin + xmax)
            yc = 0.5 * (ymin + ymax)
            small = ext_x < 0.5 * Lx
            if small and abs(xc - hw[0]) < dx and abs(yc - hw[1]) < dy:
                vol_hw.append(tag)
            else:
                vol_aq.append(tag)
    if not vol_aq or not vol_hw:
        raise RuntimeError("Volumenklassifizierung fehlgeschlagen (gmsh-Fragmentierung).")

    # Top-/Bottom-Außenflächen + Lateral-Aquifer-Flächen
    surf_top, surf_bot, surf_lat_aq = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (zmin + zmax)
        # Top/Bottom: große horizontale Außenflächen
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_top) < 1e-6:
            surf_top.append(tag)
            continue
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_base) < 1e-6:
            surf_bot.append(tag)
            continue
        # Lateral-Aquifer: vertikale Außenflächen, deren z-Bereich genau die Aquiferdicke umfasst
        on_x_edge = abs(xmin - x0) < 1e-6 and abs(xmax - x0) < 1e-6
        on_x_edge_pos = abs(xmin - (x0 + Lx)) < 1e-6 and abs(xmax - (x0 + Lx)) < 1e-6
        on_y_edge = abs(ymin - y0) < 1e-6 and abs(ymax - y0) < 1e-6
        on_y_edge_pos = abs(ymin - (y0 + Ly)) < 1e-6 and abs(ymax - (y0 + Ly)) < 1e-6
        on_outer = on_x_edge or on_x_edge_pos or on_y_edge or on_y_edge_pos
        in_aquifer_z = (zmin >= z_aq_bot - 1e-6) and (zmax <= z_aq_top + 1e-6)
        if on_outer and in_aquifer_z:
            surf_lat_aq.append(tag)

    # Hüllflächen der Brunnenbox (für Distanzfeld / Neumann‑Wärme­fluss)
    surf_hw = []
    for tag in vol_hw:
        for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False):
            if d == 2:
                surf_hw.append(abs(t))
    surf_hw = sorted(set(surf_hw))

    # Physical Groups (Reihenfolge -> tag -> MaterialID nach reindex)
    gmsh.model.addPhysicalGroup(3, vol_aq, tag=1, name="aquifer")
    gmsh.model.addPhysicalGroup(3, vol_ct, tag=2, name="caprock_top")
    gmsh.model.addPhysicalGroup(3, vol_cb, tag=3, name="caprock_bottom")
    gmsh.model.addPhysicalGroup(3, vol_hw, tag=4, name="hot_well_vol")
    gmsh.model.addPhysicalGroup(2, surf_top, tag=10, name="top")
    gmsh.model.addPhysicalGroup(2, surf_bot, tag=11, name="bottom")
    gmsh.model.addPhysicalGroup(2, surf_hw,  tag=12, name="hot_well_surf")
    if surf_lat_aq:
        gmsh.model.addPhysicalGroup(2, surf_lat_aq, tag=14, name="lateral_aquifer")

    # Hülle der Brunnenbox (für Distanzfeld)
    well_surfaces: list[int] = []
    for tag in vol_hw:
        for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False):
            if d == 2:
                well_surfaces.append(abs(t))

    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", list(set(well_surfaces)))
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", s_near)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", s_far)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", r_near)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", r_far)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    # Innerhalb der Brunnenbox sehr feines Netz (Punkte der Filterbox)
    well_points: list[tuple[int, int]] = []
    for tag in vol_hw:
        for d, t in gmsh.model.getBoundary([(3, tag)], recursive=True, oriented=False):
            if d == 0:
                well_points.append((d, t))
    if well_points:
        gmsh.model.mesh.setSize(well_points, s_in)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


def _mesh_files(cfg: dict) -> dict[str, str]:
    prefix = cfg["output"]["prefix"]
    return {
        "domain":          f"{prefix}_domain.vtu",
        "top":             f"{prefix}_physical_group_top.vtu",
        "bottom":          f"{prefix}_physical_group_bottom.vtu",
        "hot_well_vol":    f"{prefix}_physical_group_hot_well_vol.vtu",
        "hot_well_surf":   f"{prefix}_physical_group_hot_well_surf.vtu",
        "lateral_aquifer": f"{prefix}_physical_group_lateral_aquifer.vtu",
    }


def convert_mesh(cfg: dict, msh_path: Path, out_dir: Path) -> dict[str, str]:
    """gmsh-Mesh -> OGS .vtu (Domäne + Subdomänen)."""
    import ogstools as ot
    prefix = cfg["output"]["prefix"]
    meshes = ot.Meshes.from_gmsh(
        filename=str(msh_path), dim=3, reindex=True, log=False
    )
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        if name == "domain":
            fname = f"{prefix}_domain.vtu"
        else:
            fname = f"{prefix}_physical_group_{name}.vtu"
        mesh.save(str(out_dir / fname), binary=True)
    return _mesh_files(cfg)


# ======================================================================
#  Zyklus-Kurven
# ======================================================================
def build_cycle_curves(cfg: dict) -> dict:
    """Zeit-stetige Kurven für Massenquelle und Brunnen-Dirichlet-T.

    Mass-Curve:  +1 Injektion, -1 Förderung, 0 Pause.
    T-Curve:     Skalierungs-Faktor für den Dirichlet-Wert am Brunnen.
      - Während Injektion = 1.0 (Brunnen liegt auf T_inj)
      - Sonst T0/T_inj   (Brunnen ruht auf Hintergrund-T0)
    """
    n     = cfg["cycles"]["n_cycles"]
    t_c   = cfg["cycles"]["charge_days"]                  * DAY
    t_sc  = cfg["cycles"]["storage_after_charge_days"]    * DAY
    t_d   = cfg["cycles"]["discharge_days"]               * DAY
    t_sd  = cfg["cycles"]["storage_after_discharge_days"] * DAY
    ramp  = max(60.0, cfg["cycles"]["ramp_days"] * DAY)

    T0    = cfg["initial"]["T_K"]
    T_hot = cfg["operation"]["T_hot_K"]
    rh = T0 / T_hot

    # Phasen: (Name, Dauer, mass, T_curve)
    phases = [
        ("charge",          t_c,  +1.0, 1.0),
        ("storage_after_c", t_sc,  0.0, rh),
        ("discharge",       t_d,  -1.0, rh),
        ("storage_after_d", t_sd,  0.0, rh),
    ]

    times = [0.0]
    v_m, v_t = [0.0], [rh]
    t_now = 0.0
    for _cyc in range(n):
        for _name, dur, m, th in phases:
            if dur <= 0.0:
                continue
            t_now += ramp
            times.append(t_now); v_m.append(m); v_t.append(th)
            hold = max(0.0, dur - ramp)
            if hold > 0.0:
                t_now += hold
                times.append(t_now); v_m.append(m); v_t.append(th)
    t_now += ramp
    times.append(t_now); v_m.append(0.0); v_t.append(rh)

    return {
        "t_total":        t_now,
        "cycle_mass_hot": (np.array(times), np.array(v_m)),
        "cycle_T_hot":    (np.array(times), np.array(v_t)),
    }


# ======================================================================
#  XML / .prj Generierung
# ======================================================================
def _se(parent: ET.Element, tag: str, text=None, **attrs) -> ET.Element:
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None:
        el.text = str(text)
    return el


def _add_const_property(parent: ET.Element, name: str, value: float) -> None:
    p = _se(parent, "property")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)


def _add_phase_aqueous(phases: ET.Element, fluid: dict, op: dict) -> None:
    ph = _se(phases, "phase")
    _se(ph, "type", "AqueousLiquid")
    props = _se(ph, "properties")

    # Dichte: optional linear in T
    p = _se(props, "property")
    _se(p, "name", "density")
    if fluid["beta_1_per_K"] > 0.0:
        _se(p, "type", "Linear")
        _se(p, "reference_value", fluid["rho_ref_kg_m3"])
        iv = _se(p, "independent_variable")
        _se(iv, "variable_name", "temperature")
        _se(iv, "reference_condition", fluid["T_ref_K"])
        _se(iv, "slope", -fluid["rho_ref_kg_m3"] * fluid["beta_1_per_K"])
    else:
        _se(p, "type", "Constant")
        _se(p, "value", fluid["rho_ref_kg_m3"])

    _add_const_property(props, "viscosity",              fluid["viscosity_Pa_s"])
    _add_const_property(props, "specific_heat_capacity", fluid["cp_J_kgK"])
    _add_const_property(props, "thermal_conductivity",   fluid["lambda_W_mK"])
    _add_const_property(props, "storage",                op["fluid_storage_1_per_Pa"])


def _add_phase_solid(phases: ET.Element, m: dict, op: dict) -> None:
    ph = _se(phases, "phase")
    _se(ph, "type", "Solid")
    props = _se(ph, "properties")
    _add_const_property(props, "density",                m["rho_s_kg_m3"])
    _add_const_property(props, "specific_heat_capacity", m["cp_s_J_kgK"])
    _add_const_property(props, "thermal_conductivity",   m["lambda_s_W_mK"])
    _add_const_property(props, "storage",                op["solid_storage_1_per_Pa"])


def _add_medium(media: ET.Element, mid: int, mat: dict, fluid: dict,
                op: dict, disp: dict) -> None:
    med = _se(media, "medium", id=mid)
    phases = _se(med, "phases")
    _add_phase_aqueous(phases, fluid, op)
    _add_phase_solid(phases, mat, op)
    props = _se(med, "properties")
    _add_const_property(props, "porosity",     mat["porosity"])
    _add_const_property(props, "permeability", mat["permeability_m2"])

    # Effektive Wärmeleitfähigkeit: Mischung aus Phasen über Porosität
    p = _se(props, "property")
    _se(p, "name", "thermal_conductivity")
    _se(p, "type", "EffectiveThermalConductivityPorosityMixing")

    _add_const_property(props, "thermal_longitudinal_dispersivity", disp["alpha_L_m"])
    _add_const_property(props, "thermal_transversal_dispersivity",  disp["alpha_T_m"])
    _add_const_property(props, "storage", 0.0)


def _curve_xml(parent: ET.Element, name: str, t: np.ndarray, v: np.ndarray) -> None:
    c = _se(parent, "curve")
    _se(c, "name", name)
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v))


def _curve_scaled_param(parent: ET.Element, name: str,
                        curve_name: str, base_param_name: str) -> None:
    p = _se(parent, "parameter")
    _se(p, "name", name)
    _se(p, "type", "CurveScaled")
    _se(p, "curve", curve_name)
    _se(p, "parameter", base_param_name)


def _const_param(parent: ET.Element, name: str, value: float) -> None:
    p = _se(parent, "parameter")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict[str, str], curves: dict) -> Path:
    prefix = cfg["output"]["prefix"]
    fluid  = cfg["fluid"]
    op     = cfg["operation"]
    init   = cfg["initial"]
    sol    = cfg["solver"]

    # Brunnenfilter-Geometrie
    h_screen = (cfg["layers"]["aquifer_thickness_m"]
                - cfg["wells"]["screen_top_offset_m"]
                - cfg["wells"]["screen_bottom_offset_m"])
    dx_w = cfg["wells"]["screen_dx_m"]
    dy_w = cfg["wells"]["screen_dy_m"]
    V_well = dx_w * dy_w * h_screen

    Q_total = op["mass_flow_rate_kg_s"]                  # kg/s gesamt am Brunnen
    q_v_mass = Q_total / V_well                          # kg/(m³·s) – Druckeq.

    root = ET.Element("OpenGeoSysProject")

    # -- Meshes
    meshes = ET.SubElement(root, "meshes")
    _mesh_keys = ["domain", "top", "bottom", "hot_well_vol", "hot_well_surf"]
    if "lateral_aquifer" in mesh_files:
        _mesh_keys.append("lateral_aquifer")
    for k in _mesh_keys:
        _se(meshes, "mesh", mesh_files[k])

    # -- Processes
    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HT")
    _se(proc, "type", "HT")
    _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables")
    _se(pv, "temperature", "T")
    _se(pv, "pressure",    "p")
    sv = _se(proc, "secondary_variables")
    _se(sv, "secondary_variable", internal_name="darcy_velocity",
        output_name="darcy_velocity")
    bf = "0 0 -9.81" if cfg["time"]["gravity"] else "0 0 0"
    _se(proc, "specific_body_force", bf)

    # -- Media (Reihenfolge entspricht den MaterialIDs nach msh2vtu reindex)
    well_mat = dict(cfg["materials"]["aquifer"])
    well_mat["permeability_m2"] = cfg["wells"]["screen_permeability_m2"]
    disp = cfg["dispersion"]
    media = _se(root, "media")
    _add_medium(media, 0, cfg["materials"]["aquifer"],        fluid, op, disp)
    _add_medium(media, 1, cfg["materials"]["caprock_top"],    fluid, op, disp)
    _add_medium(media, 2, cfg["materials"]["caprock_bottom"], fluid, op, disp)
    _add_medium(media, 3, well_mat,                           fluid, op, disp)

    # -- Time loop
    tl = _se(root, "time_loop")
    procs = _se(tl, "processes")
    p_ref = _se(procs, "process", ref="HT")
    _se(p_ref, "nonlinear_solver", "basic_picard")
    cc = _se(p_ref, "convergence_criterion")
    _se(cc, "type",      "PerComponentDeltaX")
    _se(cc, "norm_type", "NORM2")
    _se(cc, "reltols",   f"{sol['rel_tol_T']} {sol['rel_tol_p']}")
    td = _se(p_ref, "time_discretization")
    _se(td, "type", "BackwardEuler")

    # Adaptive Zeitschrittweite: IterationNumberBasedTimeStepping kann die
    # Schrittweite bei Nichtkonvergenz reduzieren (und danach wieder
    # erhöhen) — anders als FixedTimeStepping, das an harten Phasen-
    # übergängen (Flussumkehr + Dirichlet-T-Sprung am Brunnen) sonst
    # abbricht ("Time stepper cannot reduce the time step size further").
    # Normalbetrieb bleibt beim 1-Tages-Schritt (maximum_dt = dt0).
    dt0 = cfg["time"]["dt_seconds"]
    ts = _se(p_ref, "time_stepping")
    _se(ts, "type",        "IterationNumberBasedTimeStepping")
    _se(ts, "t_initial",   0.0)
    _se(ts, "t_end",       curves["t_total"])
    _se(ts, "initial_dt",  dt0)
    _se(ts, "minimum_dt",  dt0 / 64.0)   # erlaubt Reduktion an den Übergängen
    _se(ts, "maximum_dt",  dt0)          # nie gröber als der Standard-Schritt
    # Bei wenigen Nichtlinear-Iterationen wächst dt (bis maximum_dt), bei
    # vielen schrumpft es; bei Nichtkonvergenz wird der Schritt verworfen
    # und mit kleinerem dt wiederholt.
    _se(ts, "number_iterations", "1 4 8 12")
    _se(ts, "multiplier",        "1.5 1.0 0.5 0.25")

    # -- Output: feste Ausgabezeitpunkte -> gleichmäßige Snapshots,
    #    unabhängig von der jetzt variablen Schrittweite.
    out_step = dt0 * cfg["time"]["output_every_n_steps"]
    n_out    = int(np.ceil(curves["t_total"] / out_step))
    out_times = sorted(set(
        [0.0] + [min((k + 1) * out_step, curves["t_total"]) for k in range(n_out)]
    ))
    out = _se(tl, "output")
    _se(out, "type",   "VTK")
    _se(out, "prefix", prefix)
    _se(out, "fixed_output_times", " ".join(f"{t:.6f}" for t in out_times))
    _se(out, "output_iteration_results", "false")
    vars_el = _se(out, "variables")
    for v in cfg["output"]["variables"]:
        _se(vars_el, "variable", v)

    # -- Parameters
    params = _se(root, "parameters")
    _const_param(params, "T0", init["T_K"])
    _const_param(params, "p0", init["p_Pa"])

    # Regionaler GW-Druckgradient als Function-Parameter
    gw = cfg.get("regional_gw", {})
    if gw.get("enable", False):
        import math
        rho_g = cfg["fluid"]["rho_ref_kg_m3"] * 9.81   # Pa/m je m hyd. Höhe
        i     = float(gw["gradient_m_per_m"])
        alpha = math.radians(float(gw["direction_deg"]))
        gx    = -rho_g * i * math.cos(alpha)           # Druck fällt in Strömungsrichtung
        gy    = -rho_g * i * math.sin(alpha)
        p = _se(params, "parameter")
        _se(p, "name", "p_lateral_gw")
        _se(p, "type", "Function")
        _se(p, "expression", f"{init['p_Pa']} + ({gx:.6g})*x + ({gy:.6g})*y")

    # Basisamplituden + curve-skalierte Parameter
    _const_param(params, "q_mass_amp", q_v_mass)
    _const_param(params, "T_hot_amp",  op["T_hot_K"])
    _curve_scaled_param(params, "q_mass_hot", "cycle_mass_hot", "q_mass_amp")
    _curve_scaled_param(params, "T_hot_well", "cycle_T_hot",    "T_hot_amp")

    # -- Curves
    cv = _se(root, "curves")
    for name in ("cycle_mass_hot", "cycle_T_hot"):
        t, v = curves[name]
        _curve_xml(cv, name, t, v)

    # -- Process variables
    pvars = _se(root, "process_variables")

    # Temperatur
    pv_T = _se(pvars, "process_variable")
    _se(pv_T, "name", "T")
    _se(pv_T, "components", 1)
    _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs_T = _se(pv_T, "boundary_conditions")
    # Aussenränder (Top/Bottom): konstante Hintergrund-T
    for face in ("top", "bottom"):
        bc = _se(bcs_T, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files[face]).stem)
        _se(bc, "type",      "Dirichlet")
        _se(bc, "parameter", "T0")
    # Brunnentemperatur Dirichlet-BC auf Filterbox-Volumen.
    # (Hinweis: Eine 2nd-type/Neumann-Variante auf der Brunnenhülle ist
    # physikalisch sauberer, in OGS HT auf inneren Trennflächen aber
    # numerisch instabil ohne SUPG-artige Stabilisierung.)
    bc = _se(bcs_T, "boundary_condition")
    _se(bc, "mesh",      Path(mesh_files["hot_well_vol"]).stem)
    _se(bc, "type",      "Dirichlet")
    _se(bc, "parameter", "T_hot_well")

    # Druck
    pv_p = _se(pvars, "process_variable")
    _se(pv_p, "name", "p")
    _se(pv_p, "components", 1)
    _se(pv_p, "order", 1)
    _se(pv_p, "initial_condition", "p0")
    bcs_p = _se(pv_p, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs_p, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files[face]).stem)
        _se(bc, "type",      "Dirichlet")
        _se(bc, "parameter", "p0")
    # Lateral-Aquifer als Druck-Outlet, damit das am Brunnen injizierte
    # Wasser entweichen kann. Falls regional_gw.enable=True: linearer
    # Druckgradient statt konstantem p0 → Hintergrund-Strömung.
    if "lateral_aquifer" in mesh_files:
        bc = _se(bcs_p, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files["lateral_aquifer"]).stem)
        _se(bc, "type",      "Dirichlet")
        if cfg.get("regional_gw", {}).get("enable", False):
            _se(bc, "parameter", "p_lateral_gw")
        else:
            _se(bc, "parameter", "p0")
    sts_p = _se(pv_p, "source_terms")
    st = _se(sts_p, "source_term")
    _se(st, "mesh",      Path(mesh_files["hot_well_vol"]).stem)
    _se(st, "type",      "Volumetric")
    _se(st, "parameter", "q_mass_hot")

    # -- Solvers
    nls = _se(root, "nonlinear_solvers")
    n = _se(nls, "nonlinear_solver")
    _se(n, "name",          "basic_picard")
    _se(n, "type",          "Picard")
    _se(n, "max_iter",      sol["nonlinear_iter"])
    _se(n, "linear_solver", "general_linear_solver")

    lss = _se(root, "linear_solvers")
    ls = _se(lss, "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    # BiCGSTAB+ILUT mit Skalierung: iterativ, skaliert gut auf größere Netze.
    # `scaling=true` gleicht die unterschiedlichen Größenordnungen von T und p aus.
    _se(eig, "solver_type",        "BiCGSTAB")
    _se(eig, "precon_type",        "ILUT")
    _se(eig, "max_iteration_step", sol["linear_iter"])
    _se(eig, "error_tolerance",    sol["linear_tol"])
    _se(eig, "scaling",            "true")

    # Schreiben (mit XML-Deklaration und Einrückung)
    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    tree = ET.ElementTree(root)
    tree.write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


def _indent(elem: ET.Element, level: int = 0) -> None:
    pad = "\n" + level * "    "
    if len(elem):
        if not (elem.text and elem.text.strip()):
            elem.text = pad + "    "
        for child in elem:
            _indent(child, level + 1)
        if not (elem[-1].tail and elem[-1].tail.strip()):
            elem[-1].tail = pad
    if level and not (elem.tail and elem.tail.strip()):
        elem.tail = pad


# ======================================================================
#  OGS ausführen
# ======================================================================
def run_ogs(prj_path: Path) -> int:
    ogs_exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if not ogs_exe:
        print("ogs.exe nicht im PATH gefunden – bitte installieren oder --no-run nutzen.",
              file=sys.stderr)
        return 1
    cmd = [ogs_exe, str(prj_path), "-o", str(prj_path.parent)]
    print(">>", " ".join(cmd))
    return subprocess.call(cmd)


# ======================================================================
#  CLI
# ======================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description="ATES 3D OGS demo")
    ap.add_argument("--no-mesh", action="store_true", help="Mesh nicht neu erzeugen")
    ap.add_argument("--no-run",  action="store_true", help="OGS nicht ausführen")
    args = ap.parse_args()

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print("[1/3] gmsh: 3D Schichtmodell + Brunnenbox ...")
        build_mesh(CONFIG, out_dir)
        print(f"      {msh_path}")
        print("[2/3] msh2vtu: Konvertierung in OGS-Meshes ...")
        mesh_files = convert_mesh(CONFIG, msh_path, out_dir)
    else:
        mesh_files = _mesh_files(CONFIG)

    print("[3/3] OGS-Projektdatei erzeugen ...")
    curves = build_cycle_curves(CONFIG)
    prj_path = build_prj(CONFIG, out_dir, mesh_files, curves)
    print(f"      {prj_path}  (t_end = {curves['t_total']/DAY:.1f} d)")

    if args.no_run:
        return 0

    print(">>> OGS starten")
    return run_ogs(prj_path)


if __name__ == "__main__":
    sys.exit(main())

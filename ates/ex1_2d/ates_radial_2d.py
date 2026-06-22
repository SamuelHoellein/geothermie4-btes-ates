#!/usr/bin/env python3
"""
ATES 2D radialsymmetrisch — Einführungsbeispiel für OpenGeoSys 6.

Ein einzelner Brunnen sitzt auf der Symmetrieachse (r = 0). Das 3D-Problem
wird unter Annahme rotations­symmetrischer Lösung auf eine 2D-Aufgabe in
der (r, z)-Halbebene reduziert. OGS-6 berücksichtigt die Achssymmetrie
automatisch, wenn `axially_symmetric="true"` im Mesh-Eintrag steht.

VORTEIL: Bei gleicher Auflösung 100-1000x weniger Zellen als 3D
         -> Sekunden bis wenige Minuten Laufzeit.

VERWENDUNG
----------
    python ates_radial_2d.py             # Mesh + Sim + Plots
    python ates_radial_2d.py --no-run    # nur Setup
    python ates_radial_2d.py --no-mesh   # nur .prj
    python ates_radial_2d.py --no-plots  # ohne Auto-Plots

Konventionen
------------
    x = r  (Radialkoordinate, x >= 0)
    y = z  (Vertikalkoordinate)
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import gmsh
import numpy as np

# --- ogstools >=0.8 Kompatibilitäts-Shim für die alte msh2vtu-API ---
def msh2vtu(filename, output_path, output_prefix, dim, reindex=True, log_level="WARNING"):
    import ogstools as ot
    from pathlib import Path as _P
    meshes = ot.Meshes.from_gmsh(
        filename=str(filename), dim=dim, reindex=reindex, log=False
    )
    output_path = _P(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        if name == "domain":
            fname = f"{output_prefix}_domain.vtu"
        else:
            fname = f"{output_prefix}_physical_group_{name}.vtu"
        mesh.save(str(output_path / fname), binary=True)


# ======================================================================
#  CONFIG  --  alle einstellbaren Parameter
# ======================================================================
CONFIG: dict = {

    # --- Domäne ---------------------------------------------------------
    "domain": {
        "r_max_m":  500.0,    # groß genug für 30-Jahre-Drift
        "z_base_m":   0.0,
    },

    # --- Schichtdicken --------------------------------------------------
    # Cap Rock dick genug für 30 Jahre thermische Eindringtiefe (~30 m).
    "layers": {
        "caprock_bottom_thickness_m": 60.0,
        "aquifer_thickness_m":        30.0,
        "caprock_top_thickness_m":    60.0,
    },

    # --- Brunnen --------------------------------------------------------
    "well": {
        "r_well_m":                1.0,    # Brunnenradius [m]
        "screen_top_offset_m":     5.0,    # Filterfreie Zone oben [m]
        "screen_bottom_offset_m":  5.0,
        "screen_permeability_m2":  1.0e-9, # Filterkies
    },

    # --- Mesh-Auflösung -------------------------------------------------
    "mesh": {
        "size_in_well_m":         0.4,
        "size_near_well_m":       1.0,
        "size_far_m":            15.0,    # groesseres Fernfeld
        "well_size_radius_m":     5.0,
        "well_size_radius_far_m": 120.0,
    },

    # --- Materialien ----------------------------------------------------
    "materials": {
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

    # --- Fluid ----------------------------------------------------------
    "fluid": {
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         293.15,
        "beta_1_per_K":    0.0,
        "viscosity_Pa_s":  1.0e-3,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },

    # --- Dispersivität --------------------------------------------------
    "dispersion": {
        "alpha_L_m": 5.0,
        "alpha_T_m": 0.5,
    },

    # --- Anfangsbedingungen ---------------------------------------------
    "initial": {
        "T_K":  283.15,
        "p_Pa": 0.0,
        # Geothermischer Tiefen­gradient (optional)
        #   T0(z) = T_surface_K + geothermal_gradient_K_per_m · Tiefe(z)
        # Typisch: 0.03 K/m (≈ 3 K pro 100 m). Mit 0 wird T_K verwendet.
        "T_surface_K":                  283.15,
        "geothermal_gradient_K_per_m":  0.0,
    },

    # --- Betrieb --------------------------------------------------------
    "operation": {
        "mass_flow_rate_kg_s": 0.5,
        "T_hot_K":  353.15,
        "T_cold_K": 283.15,
        "fluid_storage_1_per_Pa": 4.5e-10,
        "solid_storage_1_per_Pa": 1.0e-10,
    },

    # ------------------------------------------------------------------
    # ZYKLEN
    # ------------------------------------------------------------------
    # Zwei Modi (alternativ):
    #
    # A) 4-Phasen-Zyklus (Default):
    #    Pro Zyklus 4 aufeinander folgende Phasen mit fester Dauer.
    #    Massenstrom = ±operation.mass_flow_rate_kg_s; T_inj = T_hot bzw. T_cold.
    #
    # B) Monatsprofil (überschreibt A, falls aktiviert):
    #    Liste von 12 Monats-Speicher­leistungen [W]. Positiv = laden
    #    (Wärme in den Aquifer), negativ = fördern, 0 = Stillstand.
    #    Daraus wird pro Monat der Massenstrom berechnet:
    #        ṁ_Monat = P_Monat / ( c_p,f · ( T_inj − T_0 ) )
    #    Vorlauf-Temperatur T_inj wahlweise konstant (Default: T_hot_K bei
    #    Beladung, T_cold_K bei Förderung) oder pro Monat über
    #    monthly_T_inj_K = [T_Jan, …, T_Dez].
    #    operation.mass_flow_rate_kg_s wird dann zur Referenz-/Skalierungsgröße
    #    (z. B. Maximalwert eintragen).
    #    n_cycles entspricht der Anzahl Betriebs­jahre.
    # ------------------------------------------------------------------
    "cycles": {
        "n_cycles":                          1,       # Anzahl Zyklen (Modus A) bzw. Jahre (Modus B)
        # --- Modus A: 4-Phasen-Zyklus ---
        "charge_days":                      91.25,   # Beladung (Tage)
        "storage_after_charge_days":        91.25,   # Pause nach Beladung
        "discharge_days":                   91.25,   # Förderung
        "storage_after_discharge_days":     91.25,   # Pause nach Förderung
        "ramp_days":                         7.0,    # Übergangsrampe zwischen Phasen
        # --- Modus B: Monatsprofil (auf None für Modus A) ---
        # Beispiel (Beladung Sommer, Förderung Winter):
        #   "monthly_power_W": [+50_000, +40_000, +20_000, 0, 0, 0,
        #                         0, 0, -20_000, -40_000, -50_000, -50_000],
        "monthly_power_W":                  None,
        # Optional: monatliche Vorlauf­temperatur [K]; None → Default
        # (T_hot_K bei P>0, T_cold_K bei P<0)
        "monthly_T_inj_K":                  None,
    },

    # --- Zeit + Output --------------------------------------------------
    "time": {
        "dt_seconds":            7 * 86400.0,
        "output_every_n_steps":  1,
        "gravity":               False,
    },
    "output": {
        "prefix":    "ates_radial_2d",
        "out_dir":   "out",
        "variables": ["T", "p", "darcy_velocity"],
    },
    "solver": {
        "linear_tol":      1.0e-12,
        "linear_iter":     10000,
        "nonlinear_iter":  20,
        "rel_tol_T":       1.0e-4,
        "rel_tol_p":       1.0e-4,
    },
}

DAY = 86400.0


# ======================================================================
#  1) Mesh-Generierung (2D radial)
# ======================================================================
def build_mesh(cfg: dict, out_dir: Path) -> Path:
    """2D-Mesh in (r, z)-Halbebene. x = r, y = z."""
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"

    r_max  = cfg["domain"]["r_max_m"]
    z_base = cfg["domain"]["z_base_m"]
    t_cb = cfg["layers"]["caprock_bottom_thickness_m"]
    t_aq = cfg["layers"]["aquifer_thickness_m"]
    t_ct = cfg["layers"]["caprock_top_thickness_m"]
    z_aq_bot = z_base + t_cb
    z_aq_top = z_aq_bot + t_aq
    z_top    = z_aq_top + t_ct

    r_well = cfg["well"]["r_well_m"]
    z_sc_bot = z_aq_bot + cfg["well"]["screen_bottom_offset_m"]
    z_sc_top = z_aq_top - cfg["well"]["screen_top_offset_m"]
    h_screen = z_sc_top - z_sc_bot

    m = cfg["mesh"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ates_radial_2d")

    # Rechtecke in xy-Ebene
    r_cb = gmsh.model.occ.addRectangle(0,   z_base,    0, r_max, t_cb)
    r_aq = gmsh.model.occ.addRectangle(0,   z_aq_bot,  0, r_max, t_aq)
    r_ct = gmsh.model.occ.addRectangle(0,   z_aq_top,  0, r_max, t_ct)
    r_w  = gmsh.model.occ.addRectangle(0,   z_sc_bot,  0, r_well, h_screen)

    gmsh.model.occ.fragment([(2, r_cb)],
                            [(2, r_aq), (2, r_ct), (2, r_w)])
    gmsh.model.occ.synchronize()

    # Klassifiziere Flächen anhand Schwerpunkt
    surf_aq, surf_ct, surf_cb, surf_well = [], [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.occ.getBoundingBox(dim, tag)
        yc = 0.5 * (ymin + ymax)
        small_x = xmax < r_well + 1e-3
        if yc > z_aq_top + 1e-6:
            surf_ct.append(tag)
        elif yc < z_aq_bot - 1e-6:
            surf_cb.append(tag)
        elif small_x:
            surf_well.append(tag)
        else:
            surf_aq.append(tag)
    if not all([surf_aq, surf_ct, surf_cb, surf_well]):
        raise RuntimeError("Flaechenklassifizierung fehlgeschlagen.")

    # Randkanten klassifizieren (1D)
    edge_axis, edge_top, edge_bot, edge_far = [], [], [], []
    edge_well_surf = []   # äußere Kante der Brunnenbox
    for dim, tag in gmsh.model.getEntities(1):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.occ.getBoundingBox(dim, tag)
        # Symmetrieachse: x = 0
        if abs(xmin) < 1e-6 and abs(xmax) < 1e-6:
            edge_axis.append(tag); continue
        # Top: y = z_top
        if abs(ymin - z_top) < 1e-6 and abs(ymax - z_top) < 1e-6:
            edge_top.append(tag); continue
        # Bottom: y = z_base
        if abs(ymin - z_base) < 1e-6 and abs(ymax - z_base) < 1e-6:
            edge_bot.append(tag); continue
        # Far: x = r_max
        if abs(xmin - r_max) < 1e-6 and abs(xmax - r_max) < 1e-6:
            edge_far.append(tag); continue
        # Brunnenflanke: x = r_well (vertikal) im Aquiferbereich
        if abs(xmin - r_well) < 1e-6 and abs(xmax - r_well) < 1e-6 \
                and (z_sc_bot - 1e-6) <= ymin and ymax <= (z_sc_top + 1e-6):
            edge_well_surf.append(tag); continue

    # Physical Groups (Reihenfolge -> MaterialID nach reindex)
    gmsh.model.addPhysicalGroup(2, surf_aq,   tag=1, name="aquifer")
    gmsh.model.addPhysicalGroup(2, surf_ct,   tag=2, name="caprock_top")
    gmsh.model.addPhysicalGroup(2, surf_cb,   tag=3, name="caprock_bottom")
    gmsh.model.addPhysicalGroup(2, surf_well, tag=4, name="hot_well_vol")

    gmsh.model.addPhysicalGroup(1, edge_top, tag=10, name="top")
    gmsh.model.addPhysicalGroup(1, edge_bot, tag=11, name="bottom")
    gmsh.model.addPhysicalGroup(1, edge_far, tag=12, name="far")
    if edge_axis:
        gmsh.model.addPhysicalGroup(1, edge_axis, tag=13, name="axis")
    if edge_well_surf:
        gmsh.model.addPhysicalGroup(1, edge_well_surf, tag=14, name="hot_well_surf")

    # Mesh-Verfeinerung um die Brunnenbox
    well_edges = []
    for tag in surf_well:
        for d, t in gmsh.model.getBoundary([(2, tag)], oriented=False):
            if d == 1:
                well_edges.append(abs(t))
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList", list(set(well_edges)))
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", m["size_near_well_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", m["size_far_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", m["well_size_radius_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", m["well_size_radius_far_m"])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    # Innerhalb des Brunnens noch feiner
    well_points = []
    for tag in surf_well:
        for d, t in gmsh.model.getBoundary([(2, tag)], recursive=True, oriented=False):
            if d == 0:
                well_points.append((d, t))
    if well_points:
        gmsh.model.mesh.setSize(well_points, m["size_in_well_m"])

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(2)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


def convert_mesh(cfg: dict, msh_path: Path, out_dir: Path) -> dict[str, str]:
    pass  # msh2vtu provided at module level
    prefix = cfg["output"]["prefix"]
    msh2vtu(filename=msh_path, output_path=out_dir, output_prefix=prefix,
            dim=2, reindex=True, log_level="WARNING")
    return {
        "domain":         f"{prefix}_domain.vtu",
        "top":            f"{prefix}_physical_group_top.vtu",
        "bottom":         f"{prefix}_physical_group_bottom.vtu",
        "far":            f"{prefix}_physical_group_far.vtu",
        "hot_well_vol":   f"{prefix}_physical_group_hot_well_vol.vtu",
        "hot_well_surf":  f"{prefix}_physical_group_hot_well_surf.vtu",
    }


# ======================================================================
#  2) Saisonale Kurven
# ======================================================================
def build_cycle_curves(cfg: dict) -> dict:
    cyc   = cfg["cycles"]
    n     = cyc["n_cycles"]
    ramp  = max(60.0, cyc["ramp_days"] * DAY)

    T0       = cfg["initial"]["T_K"]
    T_hot    = cfg["operation"]["T_hot_K"]
    cp_f     = cfg["fluid"]["cp_J_kgK"]
    m_nom    = cfg["operation"]["mass_flow_rate_kg_s"]    # Nenn-/Referenz-Massenstrom
    rh       = T0 / T_hot                                  # Skalierung für Dirichlet-T

    # === Monatsprofil-Modus (überschreibt 4-Phasen-Logik) ===
    monthly_P = cyc.get("monthly_power_W")
    if monthly_P is not None:
        assert len(monthly_P) == 12, "cycles.monthly_power_W muss 12 Werte enthalten."
        monthly_T = cyc.get("monthly_T_inj_K")
        if monthly_T is not None:
            assert len(monthly_T) == 12, "cycles.monthly_T_inj_K muss 12 Werte enthalten."

        if m_nom == 0:
            raise ValueError("operation.mass_flow_rate_kg_s muss > 0 sein (Referenz-Massenstrom).")
        month_dur = 365.25 / 12.0 * DAY

        times = [0.0]; v_m = [0.0]; v_T = [rh]; t_now = 0.0
        for _ in range(n):
            for m, P_month in enumerate(monthly_P):
                T_inj = monthly_T[m] if monthly_T is not None else (
                    T_hot if P_month >= 0 else cfg["operation"]["T_cold_K"])
                dT = T_inj - T0
                # Massenstrom aus Speicherleistung:  P = m_dot · cp · (T_inj − T₀)
                if abs(dT) < 1e-6 or P_month == 0:
                    m_dot = 0.0
                else:
                    m_dot = float(P_month) / (cp_f * dT)
                m_rel = m_dot / m_nom
                T_rel = T_inj / T_hot if m_dot != 0 else rh
                t_now += ramp
                times.append(t_now); v_m.append(m_rel); v_T.append(T_rel)
                hold = max(0.0, month_dur - ramp)
                if hold > 0:
                    t_now += hold
                    times.append(t_now); v_m.append(m_rel); v_T.append(T_rel)
        t_now += ramp
        times.append(t_now); v_m.append(0.0); v_T.append(rh)
        return {
            "t_total":    t_now,
            "cycle_mass": (np.array(times), np.array(v_m)),
            "cycle_T":    (np.array(times), np.array(v_T)),
        }

    # === 4-Phasen-Modus (Default) ===
    t_c   = cyc["charge_days"]                  * DAY
    t_sc  = cyc["storage_after_charge_days"]    * DAY
    t_d   = cyc["discharge_days"]               * DAY
    t_sd  = cyc["storage_after_discharge_days"] * DAY
    phases = [
        (t_c,  +1.0, 1.0),
        (t_sc,  0.0, rh ),
        (t_d,  -1.0, rh ),
        (t_sd,  0.0, rh ),
    ]
    times = [0.0]; v_m = [0.0]; v_T = [rh]; t_now = 0.0
    for _ in range(n):
        for dur, ms, ts in phases:
            if dur <= 0: continue
            t_now += ramp
            times.append(t_now); v_m.append(ms); v_T.append(ts)
            hold = max(0.0, dur - ramp)
            if hold > 0:
                t_now += hold
                times.append(t_now); v_m.append(ms); v_T.append(ts)
    t_now += ramp
    times.append(t_now); v_m.append(0.0); v_T.append(rh)
    return {
        "t_total":    t_now,
        "cycle_mass": (np.array(times), np.array(v_m)),
        "cycle_T":    (np.array(times), np.array(v_T)),
    }


# ======================================================================
#  3) OGS-Projektdatei (axisymmetrisch)
# ======================================================================
def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None: el.text = str(text)
    return el


def _const_prop(parent, name, value):
    p = _se(parent, "property"); _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)


def _add_phase_fluid(phases, fluid, op):
    ph = _se(phases, "phase"); _se(ph, "type", "AqueousLiquid")
    props = _se(ph, "properties")
    _const_prop(props, "density",                fluid["rho_ref_kg_m3"])
    _const_prop(props, "viscosity",              fluid["viscosity_Pa_s"])
    _const_prop(props, "specific_heat_capacity", fluid["cp_J_kgK"])
    _const_prop(props, "thermal_conductivity",   fluid["lambda_W_mK"])
    _const_prop(props, "storage",                op["fluid_storage_1_per_Pa"])


def _add_phase_solid(phases, mat, op):
    ph = _se(phases, "phase"); _se(ph, "type", "Solid")
    props = _se(ph, "properties")
    _const_prop(props, "density",                mat["rho_s_kg_m3"])
    _const_prop(props, "specific_heat_capacity", mat["cp_s_J_kgK"])
    _const_prop(props, "thermal_conductivity",   mat["lambda_s_W_mK"])
    _const_prop(props, "storage",                op["solid_storage_1_per_Pa"])


def _add_medium(media, mid, mat, fluid, op, disp):
    med = _se(media, "medium", id=mid)
    phases = _se(med, "phases")
    _add_phase_fluid(phases, fluid, op)
    _add_phase_solid(phases, mat, op)
    props = _se(med, "properties")
    _const_prop(props, "porosity",     mat["porosity"])
    _const_prop(props, "permeability", mat["permeability_m2"])
    p = _se(props, "property"); _se(p, "name", "thermal_conductivity")
    _se(p, "type", "EffectiveThermalConductivityPorosityMixing")
    _const_prop(props, "thermal_longitudinal_dispersivity", disp["alpha_L_m"])
    _const_prop(props, "thermal_transversal_dispersivity",  disp["alpha_T_m"])
    _const_prop(props, "storage", 0.0)


def _curve_xml(parent, name, t, v):
    c = _se(parent, "curve"); _se(c, "name", name)
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v))


def _const_param(parent, name, value):
    p = _se(parent, "parameter"); _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)


def _curve_param(parent, name, curve_name, base):
    p = _se(parent, "parameter"); _se(p, "name", name); _se(p, "type", "CurveScaled")
    _se(p, "curve", curve_name); _se(p, "parameter", base)


def _indent(elem, level=0):
    pad = "\n" + level * "    "
    if len(elem):
        if not (elem.text and elem.text.strip()): elem.text = pad + "    "
        for child in elem: _indent(child, level + 1)
        if not (elem[-1].tail and elem[-1].tail.strip()): elem[-1].tail = pad
    if level and not (elem.tail and elem.tail.strip()): elem.tail = pad


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict, curves: dict) -> Path:
    prefix = cfg["output"]["prefix"]
    fluid, op, init, sol = cfg["fluid"], cfg["operation"], cfg["initial"], cfg["solver"]
    disp = cfg["dispersion"]

    # 3D-äquivalentes Brunnenvolumen für q_v
    h_screen = (cfg["layers"]["aquifer_thickness_m"]
                - cfg["well"]["screen_top_offset_m"]
                - cfg["well"]["screen_bottom_offset_m"])
    r_well = cfg["well"]["r_well_m"]
    V_well = np.pi * r_well**2 * h_screen          # m^3 (Zylinder)
    q_v_mass = op["mass_flow_rate_kg_s"] / V_well

    root = ET.Element("OpenGeoSysProject")

    # --- Meshes mit axially_symmetric Flag
    meshes = _se(root, "meshes")
    for key in ("domain", "top", "bottom", "far", "hot_well_vol", "hot_well_surf"):
        _se(meshes, "mesh", mesh_files[key], axially_symmetric="true")

    # --- Process
    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HT"); _se(proc, "type", "HT"); _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables"); _se(pv, "temperature", "T"); _se(pv, "pressure", "p")
    sv = _se(proc, "secondary_variables")
    _se(sv, "secondary_variable", internal_name="darcy_velocity", output_name="darcy_velocity")
    _se(proc, "specific_body_force", "0 0" if not cfg["time"]["gravity"] else "0 -9.81")

    # --- Media
    well_mat = dict(cfg["materials"]["aquifer"])
    well_mat["permeability_m2"] = cfg["well"]["screen_permeability_m2"]
    media = _se(root, "media")
    _add_medium(media, 0, cfg["materials"]["aquifer"],        fluid, op, disp)
    _add_medium(media, 1, cfg["materials"]["caprock_top"],    fluid, op, disp)
    _add_medium(media, 2, cfg["materials"]["caprock_bottom"], fluid, op, disp)
    _add_medium(media, 3, well_mat,                           fluid, op, disp)

    # --- Time loop
    tl = _se(root, "time_loop")
    procs = _se(tl, "processes")
    p_ref = _se(procs, "process", ref="HT")
    _se(p_ref, "nonlinear_solver", "basic_picard")
    cc = _se(p_ref, "convergence_criterion")
    _se(cc, "type", "PerComponentDeltaX"); _se(cc, "norm_type", "NORM2")
    _se(cc, "reltols", f"{sol['rel_tol_T']} {sol['rel_tol_p']}")
    td = _se(p_ref, "time_discretization"); _se(td, "type", "BackwardEuler")
    ts = _se(p_ref, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0.0); _se(ts, "t_end", curves["t_total"])
    steps = _se(ts, "timesteps"); pair = _se(steps, "pair")
    n_steps = int(np.ceil(curves["t_total"] / cfg["time"]["dt_seconds"]))
    _se(pair, "repeat", n_steps); _se(pair, "delta_t", cfg["time"]["dt_seconds"])

    out = _se(tl, "output")
    _se(out, "type", "VTK"); _se(out, "prefix", prefix)
    out_steps = _se(out, "timesteps"); pair = _se(out_steps, "pair")
    _se(pair, "repeat", n_steps); _se(pair, "each_steps", cfg["time"]["output_every_n_steps"])
    _se(out, "output_iteration_results", "false")
    vars_el = _se(out, "variables")
    for v in cfg["output"]["variables"]:
        _se(vars_el, "variable", v)

    # --- Parameters
    params = _se(root, "parameters")
    gradient = init.get("geothermal_gradient_K_per_m", 0.0)
    if abs(gradient) > 1e-12:
        T_surf = init.get("T_surface_K", init["T_K"])
        z_tot  = (cfg["layers"]["caprock_bottom_thickness_m"]
                  + cfg["layers"]["aquifer_thickness_m"]
                  + cfg["layers"]["caprock_top_thickness_m"])
        p_el = _se(params, "parameter")
        _se(p_el, "name", "T0")
        _se(p_el, "type", "Function")
        _se(p_el, "expression", f"{T_surf} + ({gradient:.6g})*({z_tot} - y)")
    else:
        _const_param(params, "T0", init["T_K"])
    _const_param(params, "p0",   init["p_Pa"])
    _const_param(params, "q_mass_amp", q_v_mass)
    _const_param(params, "T_hot_amp",  op["T_hot_K"])
    _curve_param(params, "q_mass_well", "cycle_mass", "q_mass_amp")
    _curve_param(params, "T_well",      "cycle_T",    "T_hot_amp")

    # --- Curves
    cv = _se(root, "curves")
    for name in ("cycle_mass", "cycle_T"):
        t, v = curves[name]; _curve_xml(cv, name, t, v)

    # --- Process variables
    pvars = _se(root, "process_variables")

    pv_T = _se(pvars, "process_variable")
    _se(pv_T, "name", "T"); _se(pv_T, "components", 1); _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs = _se(pv_T, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T0")
    bc = _se(bcs, "boundary_condition")
    _se(bc, "mesh", Path(mesh_files["hot_well_vol"]).stem)
    _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T_well")

    pv_p = _se(pvars, "process_variable")
    _se(pv_p, "name", "p"); _se(pv_p, "components", 1); _se(pv_p, "order", 1)
    _se(pv_p, "initial_condition", "p0")
    bcs = _se(pv_p, "boundary_conditions")
    for face in ("top", "bottom", "far"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "p0")
    sts = _se(pv_p, "source_terms")
    st = _se(sts, "source_term")
    _se(st, "mesh", Path(mesh_files["hot_well_vol"]).stem)
    _se(st, "type", "Volumetric"); _se(st, "parameter", "q_mass_well")

    # --- Solvers
    nls = _se(root, "nonlinear_solvers")
    n = _se(nls, "nonlinear_solver")
    _se(n, "name", "basic_picard"); _se(n, "type", "Picard")
    _se(n, "max_iter", sol["nonlinear_iter"]); _se(n, "linear_solver", "general_linear_solver")
    lss = _se(root, "linear_solvers")
    ls = _se(lss, "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    _se(eig, "solver_type", "BiCGSTAB"); _se(eig, "precon_type", "ILUT")
    _se(eig, "max_iteration_step", sol["linear_iter"])
    _se(eig, "error_tolerance",    sol["linear_tol"])
    _se(eig, "scaling", "true")

    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


# ======================================================================
#  4) Run + Plots
# ======================================================================
def run_ogs(prj_path: Path) -> int:
    ogs_exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if not ogs_exe:
        print("ogs.exe nicht im PATH", file=sys.stderr); return 1
    cmd = [ogs_exe, str(prj_path), "-o", str(prj_path.parent)]
    print(">>", " ".join(cmd))
    return subprocess.call(cmd)


def make_plots(cfg: dict, out_dir: Path) -> None:
    try:
        import pyvista as pv
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("pyvista/matplotlib fehlt, Plots uebersprungen."); return

    figdir = out_dir / "figures"; figdir.mkdir(exist_ok=True)
    prefix = cfg["output"]["prefix"]
    files = sorted(out_dir.glob(f"{prefix}_ts_*_t_*.vtu"),
                   key=lambda p: int(re.search(r"_ts_(\d+)_", p.name).group(1)))
    if not files: print("Keine VTU-Ausgaben."); return

    T0 = cfg["initial"]["T_K"]
    T_hot = cfg["operation"]["T_hot_K"]
    z_aq_mid = (cfg["domain"]["z_base_m"]
                + cfg["layers"]["caprock_bottom_thickness_m"]
                + cfg["layers"]["aquifer_thickness_m"] / 2.0)
    t_cycle = sum(cfg["cycles"][k] for k in
                  ("charge_days","storage_after_charge_days",
                   "discharge_days","storage_after_discharge_days"))

    def file_for(day):
        return min(files, key=lambda p: abs(int(re.search(r"_t_(\d+)", p.name).group(1))/DAY - day))

    # 4.1 T-Feld an 5 Phasen-Zeitpunkten (xy-Schnitt = (r,z)-Halbebene)
    keypoints = [
        (1,  "Tag 1"),
        (cfg["cycles"]["charge_days"],  "Ende Beladung"),
        (cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"], "Ende Pause 1"),
        (cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"]
         + cfg["cycles"]["discharge_days"], "Ende Förderung"),
        (t_cycle, "Ende Pause 2"),
    ]
    plotter = pv.Plotter(off_screen=True, window_size=(1800, 500), shape=(1, len(keypoints)))
    pv.OFF_SCREEN = True
    for i, (day, label) in enumerate(keypoints):
        m = pv.read(file_for(day))
        plotter.subplot(0, i)
        plotter.add_mesh(m, scalars="T", cmap="coolwarm", clim=[T0, T_hot],
                         show_scalar_bar=(i == len(keypoints) - 1),
                         scalar_bar_args={"title": "T [K]"} if i == len(keypoints)-1 else None,
                         show_edges=False)
        plotter.add_text(f"{label} (Tag {day:.0f})", font_size=9, position="upper_edge")
        plotter.view_xy()
    plotter.screenshot(str(figdir / "T_field_2D.png"))
    print("  saved T_field_2D.png")

    # 4.2 T(t) an verschiedenen Radien
    probes = {
        "Brunnen (r=1 m)": (cfg["well"]["r_well_m"] - 0.1, z_aq_mid, 0),
        "r=5 m":           (5.0,  z_aq_mid, 0),
        "r=20 m":          (20.0, z_aq_mid, 0),
        "r=50 m":          (50.0, z_aq_mid, 0),
        "Cap Rock 5 m":    (5.0,  z_aq_mid + cfg["layers"]["aquifer_thickness_m"]/2 + 5, 0),
    }
    times, series = [], {k: [] for k in probes}
    for f in files:
        m = pv.read(f)
        t_d = int(re.search(r"_t_(\d+)", f.name).group(1)) / DAY
        times.append(t_d)
        pts = pv.PolyData(np.array(list(probes.values())))
        s = pts.sample(m)
        for i, k in enumerate(probes):
            series[k].append(float(s["T"][i]))
    times = np.array(times)
    fig, ax = plt.subplots(figsize=(12, 4.5))
    n_cyc = cfg["cycles"]["n_cycles"]
    for cyc in range(n_cyc):
        t0c = cyc * t_cycle
        ax.axvspan(t0c, t0c + cfg["cycles"]["charge_days"], alpha=0.08, color="red")
        ax.axvspan(t0c + cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"],
                   t0c + cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"]
                   + cfg["cycles"]["discharge_days"], alpha=0.08, color="blue")
    for k, T in series.items():
        ax.plot(times, T, lw=1.5, label=k)
    ax.axhline(T0, color="k", lw=0.6, ls=":")
    ax.axhline(T_hot, color="k", lw=0.6, ls=":")
    ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("T [K]")
    ax.set_title("ATES 2D radial — T(t) an Beobachtungspunkten")
    ax.legend(loc="upper right", fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figdir / "T_vs_time.png", dpi=130); plt.close(fig)
    print("  saved T_vs_time.png")

    # 4.3 Energiebilanz (axisymmetrisch: 2D Fläche * 2π * r für Volumen)
    domain = pv.read(out_dir / f"{prefix}_domain.vtu")
    mid = domain.cell_data["MaterialIDs"]
    cs = domain.compute_cell_sizes()
    cell_area = cs["Area"]   # 2D-Fläche
    # Centroid r-Koordinate je Zelle
    ctrs = domain.cell_centers()
    r_centroid = ctrs.points[:, 0]
    cell_vol_axi = 2 * np.pi * r_centroid * cell_area   # echtes 3D-Volumen
    a = cfg["materials"]["aquifer"]; c = cfg["materials"]["caprock_top"]
    rcp_aq = a["porosity"]*1000*4180 + (1-a["porosity"])*a["rho_s_kg_m3"]*a["cp_s_J_kgK"]
    rcp_cr = c["porosity"]*1000*4180 + (1-c["porosity"])*c["rho_s_kg_m3"]*c["cp_s_J_kgK"]
    rcp = np.where(mid == 0, rcp_aq,
          np.where(np.isin(mid, [1, 2]), rcp_cr, rcp_aq))

    E_tot = []
    for f in files:
        m = pv.read(f)
        cT = m.point_data_to_cell_data(pass_point_data=False)["T"]
        E_tot.append(float(np.sum(rcp * cell_vol_axi * (cT - T0))) / 1e9)
    E_tot = np.array(E_tot)

    fig, ax = plt.subplots(figsize=(12, 4.5))
    for cyc in range(n_cyc):
        t0c = cyc * t_cycle
        ax.axvspan(t0c, t0c + cfg["cycles"]["charge_days"], alpha=0.08, color="red", label="Beladung" if cyc==0 else None)
        ax.axvspan(t0c + cfg["cycles"]["charge_days"], t0c + cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"],
                   alpha=0.08, color="gray", label="Pause 1" if cyc==0 else None)
        t2 = t0c + cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"]
        ax.axvspan(t2, t2 + cfg["cycles"]["discharge_days"], alpha=0.08, color="blue", label="Förderung" if cyc==0 else None)
        ax.axvspan(t2 + cfg["cycles"]["discharge_days"], t0c + t_cycle, alpha=0.08, color="lightgray", label="Pause 2" if cyc==0 else None)
    ax.plot(times, E_tot, "k-", lw=2.0, label="Gespeichert (axi-symm.)")
    i_p1 = int(np.argmin(np.abs(times - (cfg["cycles"]["charge_days"]
                              + cfg["cycles"]["storage_after_charge_days"]))))
    i_f  = int(np.argmin(np.abs(times - (cfg["cycles"]["charge_days"]
                              + cfg["cycles"]["storage_after_charge_days"]
                              + cfg["cycles"]["discharge_days"]))))
    if E_tot[i_p1] > 0:
        rec = (E_tot[i_p1] - E_tot[i_f]) / E_tot[i_p1] * 100
        ax.annotate(f"Recovery Z1: {rec:.1f}%", xy=(times[i_f], E_tot[i_f]),
                    xytext=(times[i_f]*0.5, E_tot[i_p1]*0.6),
                    fontsize=10, arrowprops=dict(arrowstyle="->", color="gray"))
    ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("Energie über T0 [GJ]")
    ax.set_title("ATES 2D radial — Energiebilanz")
    ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=9, ncol=2)
    ax.set_xlim(0, times[-1]); ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(figdir / "energy_balance.png", dpi=130); plt.close(fig)
    print("  saved energy_balance.png")


# ======================================================================
#  Hauptfunktion
# ======================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description="ATES 2D radial demo")
    ap.add_argument("--no-mesh",  action="store_true")
    ap.add_argument("--no-run",   action="store_true")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print("[1/4] gmsh: 2D Radialmesh erzeugen ...")
        build_mesh(CONFIG, out_dir)
        print(f"      -> {msh_path}")
        print("[2/4] msh2vtu ...")
        mesh_files = convert_mesh(CONFIG, msh_path, out_dir)
    else:
        mesh_files = {k: f"{prefix}_physical_group_{n}.vtu" if "physical" in n else f"{prefix}_{n}.vtu"
                      for k, n in [("domain","domain"),("top","top"),("bottom","bottom"),
                                   ("far","far"),("hot_well_vol","hot_well_vol"),
                                   ("hot_well_surf","hot_well_surf")]}
        # Fix prefix
        mesh_files = {
            "domain":         f"{prefix}_domain.vtu",
            "top":            f"{prefix}_physical_group_top.vtu",
            "bottom":         f"{prefix}_physical_group_bottom.vtu",
            "far":            f"{prefix}_physical_group_far.vtu",
            "hot_well_vol":   f"{prefix}_physical_group_hot_well_vol.vtu",
            "hot_well_surf": f"{prefix}_physical_group_hot_well_surf.vtu",
        }

    print("[3/4] OGS-Projektdatei ...")
    curves = build_cycle_curves(CONFIG)
    prj_path = build_prj(CONFIG, out_dir, mesh_files, curves)
    print(f"      -> {prj_path}  (t_end = {curves['t_total']/DAY:.0f} d)")

    if args.no_run: return 0
    print("[4/4] OGS starten ...")
    rc = run_ogs(prj_path)
    if rc != 0: return rc
    if not args.no_plots:
        print("[5/4] Plots erzeugen ...")
        make_plots(CONFIG, out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
BTES 3D demo for OpenGeoSys 6 — Borehole Thermal Energy Storage.

Sondenfeld N x N im Untergrund. Je Sonde eine kleine Box als Subdomäne,
auf der eine volumetrische Wärmequelle aufgeprägt wird (positiv bei Beladung,
negativ bei Entladung). Wärmetransfer dominant über Wärmeleitung im Boden.

Geometrie- und Materialmodell sind mit der 2D-Übung (btes_radial_2d.py)
vereinheitlicht: der Untergrund ist eine freie Schichtliste (`layers`, von
oben nach unten), jede Schicht mit eigenen Materialwerten; die Sonde wird
über `borehole.depth_top_m`/`depth_bottom_m` (Tiefe unter Oberfläche)
positioniert und kann mehrere Schichten durchstoßen.

Im Gegensatz zum ATES (ates_3d.py) gibt es keinen Massenstrom — Strömungs-
permeabilität wird sehr niedrig gewählt, sodass die HT-Prozessgleichung
faktisch zu reiner Wärmeleitung degeneriert.

Aufruf:
    python btes_3d.py            # Mesh + .prj + OGS-Lauf
    python btes_3d.py --no-run   # nur Setup, kein OGS
    python btes_3d.py --no-mesh  # nur .prj erzeugen
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
#  CONFIG  --  hier alles anpassen
# ======================================================================
CONFIG: dict = {
    "domain": {
        "size_x_m":  80.0,
        "size_y_m":  80.0,
        "z_base_m":   0.0,
    },
    # ------------------------------------------------------------------
    # SCHICHTEN  (von OBEN nach UNTEN, wie in der 2D-Übung)
    # ------------------------------------------------------------------
    # Freie Liste von Bodenschichten. Jede Schicht hat eigene Material­werte.
    # Die Schicht­dicken addieren sich zur Gesamttiefe der Domain. Auf eine
    # einzelne Bodenklasse zurückzufallen: einfach nur einen Eintrag setzen.
    "layers": [
        # name              thickness_m  permeability  porosity  rho_s  cp_s  lambda_s
        {"name": "cover",    "thickness_m":  5.0,
         "permeability_m2": 1.0e-15, "porosity": 0.35,
         "rho_s_kg_m3": 1900.0, "cp_s_J_kgK": 1500.0, "lambda_s_W_mK": 1.4},
        {"name": "bedrock",  "thickness_m": 80.0,
         "permeability_m2": 1.0e-18, "porosity": 0.20,
         "rho_s_kg_m3": 2700.0, "cp_s_J_kgK":  900.0, "lambda_s_W_mK": 2.5},
        {"name": "basement", "thickness_m": 30.0,
         "permeability_m2": 1.0e-19, "porosity": 0.10,
         "rho_s_kg_m3": 2750.0, "cp_s_J_kgK":  850.0, "lambda_s_W_mK": 3.0},
    ],
    # ------------------------------------------------------------------
    # SONDEN-GEOMETRIE  (vereinheitlicht mit der 2D-Übung)
    # ------------------------------------------------------------------
    # Sondenkopf-/Sondenfuß-Tiefe gemessen VON DER OBERFLÄCHE NACH UNTEN.
    # Die Sonde kann mehrere Schichten durchstoßen. Jede Sonde wird als
    # kleine Box (borehole_dx_m × borehole_dy_m) modelliert.
    "borehole": {
        "depth_top_m":     7.0,    # Sondenkopf, Tiefe unter Oberfläche [m]
        "depth_bottom_m": 83.0,    # Sondenfuß,  Tiefe unter Oberfläche [m]
        "borehole_dx_m":   0.6,    # Sondenbox-Ausdehnung x [m]
        "borehole_dy_m":   0.6,    # Sondenbox-Ausdehnung y [m]
    },
    "field": {
        # Sondenfeld – entweder N x N Raster (n_x, n_y, spacing_m) ODER
        # eine explizite Positionsliste positions=[(x,y), ...]
        "n_x":           3,
        "n_y":           3,
        "spacing_m":     5.0,
        "positions":     None,        # bei nicht-None überschreibt es das Raster
    },
    "mesh": {
        "size_in_borehole_m":    0.4,
        "size_near_field_m":     1.5,
        "size_far_m":           15.0,
        "field_size_radius_m":   8.0,    # bis hier feinmaschig
        "field_size_radius_far_m": 30.0,
    },
    "fluid": {
        # Porenfluid (Wasser) — auch wenn k klein ist, brauchen wir die Phase
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         293.15,
        "beta_1_per_K":    0.0,
        "viscosity_Pa_s":  1.0e-3,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },
    "dispersion": {
        "alpha_L_m": 0.0,    # ohne Strömung irrelevant
        "alpha_T_m": 0.0,
    },
    "initial": {
        "T_K":  283.15,
        "p_Pa": 0.0,
        # Optionaler geothermischer Gradient (wie in der 2D-Übung). Bei
        # geothermal_gradient_K_per_m > 0 wird T0(z) tiefen­abhängig gesetzt.
        "T_surface_K":                 283.15,
        "geothermal_gradient_K_per_m": 0.0,    # 0.03 für realistischen Gradient
    },
    "operation": {
        # Heat-Source je Sonde: Leistung [W] (positiv = Beladung)
        "power_per_borehole_W":  2000.0,
        # Effektive Speicherzahlen (im HT-Prozess erforderlich)
        "fluid_storage_1_per_Pa": 4.5e-10,
        "solid_storage_1_per_Pa": 1.0e-10,
    },
    # ------------------------------------------------------------------
    # ZYKLEN – HIER FÜR STUDIERENDE
    # ------------------------------------------------------------------
    # Ein vollständiger Zyklus besteht aus 4 aufeinander folgenden Phasen:
    #   1) charge                    – Beladung (Wärme einspeisen)
    #   2) storage_after_charge      – Pause/Speicherung nach Beladung
    #   3) discharge                 – Förderung (Wärme entziehen)
    #   4) storage_after_discharge   – Pause/Speicherung nach Förderung
    #
    # Periode T_Zyklus = charge + storage_after_charge + discharge + storage_after_discharge
    # Gesamt­simulations­zeit = n_cycles * T_Zyklus
    #
    # Alternativ — Modus B (Monatsprofil): Setze cycles.monthly_power_W auf eine
    # Liste von 12 Monatsleistungen [W] (positiv = laden, negativ = fördern,
    # 0 = Stillstand). Dann wird die 4-Phasen-Logik überschrieben; jeder Monat
    # dauert 365.25/12 ≈ 30.44 d und die Sequenz wird n_cycles-mal (= Jahre)
    # wiederholt. operation.power_per_borehole_W dient dabei als Referenz­leistung
    # (Skalierung der OGS-Curve). Auf None lassen für Modus A.
    # ------------------------------------------------------------------
    "cycles": {
        "n_cycles":                        1,       # Anzahl Zyklen (A) bzw. Jahre (B)
        "charge_days":                     91.25,   # Phase 1: Beladung (Tage)
        "storage_after_charge_days":       91.25,   # Phase 2: Pause nach Beladung (Tage)
        "discharge_days":                  91.25,   # Phase 3: Förderung (Tage)
        "storage_after_discharge_days":    91.25,   # Phase 4: Pause nach Förderung (Tage)
        "ramp_days":                       7.0,     # Sanfte Übergangsrampe zwischen Phasen (Tage)
        # --- Modus B: Monatsprofil (auf None für Modus A) ---
        # Beispiel (Beladung Sommer, Förderung Winter):
        #   "monthly_power_W": [+2000, +2000, +1500, +500, 0, 0,
        #                         0, 0, -1500, -2500, -3000, -2500],
        "monthly_power_W":                 None,
    },
    "time": {
        "dt_seconds":           7 * 86400.0,
        "output_every_n_steps": 1,
        "gravity":              False,
    },
    "output": {
        "prefix":    "btes_3d",
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


def _layer_stack(cfg: dict):
    """Schichtliste von OBEN nach UNTEN (wie in cfg) in z-Grenzen umrechnen.
    Rückgabe: (layers_bottom_up, z_top), wobei jede Schicht zusätzlich
    {z_low, z_high} trägt und z_top die Oberfläche ist.
    """
    z_base = cfg["domain"].get("z_base_m", 0.0)
    bot_up = list(reversed(list(cfg["layers"])))  # bottom → top
    z = z_base
    out = []
    for L in bot_up:
        z_low = z
        z_high = z + float(L["thickness_m"])
        out.append({**L, "z_low": z_low, "z_high": z_high})
        z = z_high
    return out, z


def _borehole_positions(cfg: dict) -> list[tuple[float, float]]:
    custom = cfg["field"].get("positions")
    if custom:
        return [(float(x), float(y)) for x, y in custom]
    nx = cfg["field"]["n_x"]; ny = cfg["field"]["n_y"]; s = cfg["field"]["spacing_m"]
    x0_f = -(nx - 1) * s / 2.0
    y0_f = -(ny - 1) * s / 2.0
    return [(x0_f + ix * s, y0_f + iy * s) for ix in range(nx) for iy in range(ny)]


def _n_boreholes(cfg: dict) -> int:
    return len(_borehole_positions(cfg))


# ======================================================================
#  Mesh — gmsh
# ======================================================================
def build_mesh(cfg: dict, out_dir: Path) -> Path:
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"

    Lx = cfg["domain"]["size_x_m"]
    Ly = cfg["domain"]["size_y_m"]
    z_base = cfg["domain"]["z_base_m"]
    layers, z_top = _layer_stack(cfg)

    # Sonde: Tiefe von der Oberfläche nach unten
    dx_b = cfg["borehole"]["borehole_dx_m"]
    dy_b = cfg["borehole"]["borehole_dy_m"]
    depth_top    = cfg["borehole"]["depth_top_m"]
    depth_bottom = cfg["borehole"]["depth_bottom_m"]
    z_bh_top = z_top - depth_top
    z_bh_bot = z_top - depth_bottom
    h_bh = z_bh_top - z_bh_bot
    if h_bh <= 0:
        raise ValueError(f"borehole.depth_bottom_m ({depth_bottom}) muss > depth_top_m ({depth_top}) sein.")
    if z_bh_bot < z_base - 1e-9 or z_bh_top > z_top + 1e-9:
        raise ValueError("Sondentiefe liegt außerhalb der Schichtdomäne.")

    bh_positions = _borehole_positions(cfg)

    s_in   = cfg["mesh"]["size_in_borehole_m"]
    s_near = cfg["mesh"]["size_near_field_m"]
    s_far  = cfg["mesh"]["size_far_m"]
    r_near = cfg["mesh"]["field_size_radius_m"]
    r_far  = cfg["mesh"]["field_size_radius_far_m"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("btes")

    x0, y0 = -Lx / 2.0, -Ly / 2.0

    # Eine Box pro Schicht (von unten nach oben)
    layer_boxes = []
    for L in layers:
        b = gmsh.model.occ.addBox(x0, y0, L["z_low"], Lx, Ly, L["z_high"] - L["z_low"])
        layer_boxes.append(b)

    bh_boxes = []
    for x, y in bh_positions:
        b = gmsh.model.occ.addBox(x - dx_b / 2.0, y - dy_b / 2.0, z_bh_bot,
                                  dx_b, dy_b, h_bh)
        bh_boxes.append(b)

    gmsh.model.occ.fragment(
        [(3, layer_boxes[0])],
        [(3, b) for b in layer_boxes[1:]] + [(3, b) for b in bh_boxes],
    )
    gmsh.model.occ.synchronize()

    # Klassifikation der entstehenden Volumen:
    #  - kleine Grundfläche im Sonden-z-Bereich  -> Sondenstück (nach (x,y) gruppiert)
    #  - sonst                                    -> Bodenschicht (nach z-Mitte)
    bh_pos = np.array(bh_positions)
    vol_bh = {i: [] for i in range(len(bh_positions))}   # borehole index -> [tags]
    vol_layer = {i: [] for i in range(len(layers))}      # layer index    -> [tags]
    for dim, tag in gmsh.model.getEntities(3):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        ext = xmax - xmin
        zc = 0.5 * (zmin + zmax)
        xc = 0.5 * (xmin + xmax); yc = 0.5 * (ymin + ymax)
        small = ext < 0.3 * Lx
        if small and (z_bh_bot - 1e-3) <= zc <= (z_bh_top + 1e-3):
            d = np.hypot(bh_pos[:, 0] - xc, bh_pos[:, 1] - yc)
            vol_bh[int(np.argmin(d))].append(tag)
        else:
            for i, L in enumerate(layers):
                if L["z_low"] - 1e-6 <= zc <= L["z_high"] + 1e-6:
                    vol_layer[i].append(tag); break
    missing = [i for i, t in vol_bh.items() if not t]
    if missing:
        raise RuntimeError(f"Sonden ohne Volumen gefunden (Index {missing}).")

    # Außenflächen
    surf_top, surf_bot, surf_lat = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (zmin + zmax)
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_top) < 1e-6:
            surf_top.append(tag); continue
        if (xmax - xmin) >= 0.9 * Lx and abs(zc - z_base) < 1e-6:
            surf_bot.append(tag); continue
        on_outer = (abs(xmin - x0) < 1e-6 and abs(xmax - x0) < 1e-6) \
                   or (abs(xmin - (x0 + Lx)) < 1e-6 and abs(xmax - (x0 + Lx)) < 1e-6) \
                   or (abs(ymin - y0) < 1e-6 and abs(ymax - y0) < 1e-6) \
                   or (abs(ymin - (y0 + Ly)) < 1e-6 and abs(ymax - (y0 + Ly)) < 1e-6)
        if on_outer:
            surf_lat.append(tag)

    # Physical groups – fortlaufende Tags: erst Schichten (unten→oben),
    # dann Sonden, damit msh2vtu reindex MaterialIDs 0..(L-1) für die
    # Schichten und L..(L+N-1) für die Sonden vergibt.
    pg = 1
    for i, L in enumerate(layers):
        gmsh.model.addPhysicalGroup(3, vol_layer[i], tag=pg, name=L["name"]); pg += 1
    for i in range(len(bh_positions)):
        gmsh.model.addPhysicalGroup(3, vol_bh[i], tag=pg, name=f"bh_{i:02d}"); pg += 1
    gmsh.model.addPhysicalGroup(2, surf_top, tag=100, name="top")
    gmsh.model.addPhysicalGroup(2, surf_bot, tag=101, name="bottom")
    if surf_lat:
        gmsh.model.addPhysicalGroup(2, surf_lat, tag=102, name="lateral")

    # Mesh-Größenfeld: feiner um Sondenfeld
    all_bh_tags = [t for tags in vol_bh.values() for t in tags]
    bh_surfaces = []
    for tag in all_bh_tags:
        for d, t in gmsh.model.getBoundary([(3, tag)], oriented=False):
            if d == 2:
                bh_surfaces.append(abs(t))
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "SurfacesList", list(set(bh_surfaces)))
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", s_near)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", s_far)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", r_near)
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", r_far)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    bh_points = []
    for tag in all_bh_tags:
        for d, t in gmsh.model.getBoundary([(3, tag)], recursive=True, oriented=False):
            if d == 0:
                bh_points.append((d, t))
    if bh_points:
        gmsh.model.mesh.setSize(bh_points, s_in)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


def _mesh_files(cfg: dict) -> dict:
    prefix = cfg["output"]["prefix"]
    n_bh = _n_boreholes(cfg)
    files = {
        "domain":  f"{prefix}_domain.vtu",
        "top":     f"{prefix}_physical_group_top.vtu",
        "bottom":  f"{prefix}_physical_group_bottom.vtu",
        "lateral": f"{prefix}_physical_group_lateral.vtu",
        "_n_bh":   n_bh,
    }
    for L in cfg["layers"]:
        files[L["name"]] = f"{prefix}_physical_group_{L['name']}.vtu"
    for i in range(n_bh):
        files[f"bh_{i:02d}"] = f"{prefix}_physical_group_bh_{i:02d}.vtu"
    return files


def convert_mesh(cfg: dict, msh_path: Path, out_dir: Path) -> dict:
    prefix = cfg["output"]["prefix"]
    msh2vtu(filename=msh_path, output_path=out_dir, output_prefix=prefix,
            dim=3, reindex=True, log_level="WARNING")
    return _mesh_files(cfg)


# ======================================================================
#  Zyklen-Kurve
# ======================================================================
def build_cycle_curves(cfg: dict) -> dict:
    """Eine Kurve für alle Sonden: +1 Beladung, 0 Pause, -1 Förderung.

    Modus A (Default): 4-Phasen-Zyklus.
    Modus B (cycles.monthly_power_W ≠ None): 12 Monatsleistungen [W], skaliert
    auf operation.power_per_borehole_W (Referenzleistung).
    """
    cyc   = cfg["cycles"]
    ramp  = max(60.0, cyc["ramp_days"] * DAY)
    n     = cyc["n_cycles"]

    # === Modus B: Monatsprofil (überschreibt 4-Phasen-Logik) ===
    monthly = cyc.get("monthly_power_W")
    if monthly is not None:
        assert len(monthly) == 12, "cycles.monthly_power_W muss 12 Werte enthalten."
        P_nominal = cfg["operation"]["power_per_borehole_W"]
        if P_nominal == 0:
            raise ValueError("operation.power_per_borehole_W muss > 0 sein (Referenzleistung).")
        month_dur = 365.25 / 12.0 * DAY      # ~30.44 d
        times = [0.0]; vals = [0.0]; t_now = 0.0
        for _ in range(n):
            for P_month in monthly:
                q_rel = float(P_month) / P_nominal
                t_now += ramp; times.append(t_now); vals.append(q_rel)
                hold = max(0.0, month_dur - ramp)
                if hold > 0.0:
                    t_now += hold; times.append(t_now); vals.append(q_rel)
        t_now += ramp; times.append(t_now); vals.append(0.0)
        return {
            "t_total":  t_now,
            "cycle_q":  (np.array(times), np.array(vals)),
        }

    # === Modus A: 4-Phasen-Zyklus (Default) ===
    t_c   = cyc["charge_days"]                  * DAY
    t_sc  = cyc["storage_after_charge_days"]    * DAY
    t_d   = cyc["discharge_days"]               * DAY
    t_sd  = cyc["storage_after_discharge_days"] * DAY

    phases = [
        ("charge",          t_c,  +1.0),
        ("storage_after_c", t_sc,  0.0),
        ("discharge",       t_d,  -1.0),
        ("storage_after_d", t_sd,  0.0),
    ]

    times = [0.0]
    vals  = [0.0]
    t_now = 0.0
    for _ in range(n):
        for _name, dur, q in phases:
            if dur <= 0.0: continue
            t_now += ramp
            times.append(t_now); vals.append(q)
            hold = max(0.0, dur - ramp)
            if hold > 0.0:
                t_now += hold
                times.append(t_now); vals.append(q)
    t_now += ramp
    times.append(t_now); vals.append(0.0)
    return {
        "t_total":  t_now,
        "cycle_q":  (np.array(times), np.array(vals)),
    }


# ======================================================================
#  XML / .prj
# ======================================================================
def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None:
        el.text = str(text)
    return el


def _const_param(parent, name, value):
    p = _se(parent, "parameter")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)


def _curve_scaled_param(parent, name, curve_name, base_param_name):
    p = _se(parent, "parameter")
    _se(p, "name", name)
    _se(p, "type", "CurveScaled")
    _se(p, "curve", curve_name)
    _se(p, "parameter", base_param_name)


def _add_const_property(parent, name, value):
    p = _se(parent, "property")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)


def _add_phase_aqueous(phases, fluid, op):
    ph = _se(phases, "phase")
    _se(ph, "type", "AqueousLiquid")
    props = _se(ph, "properties")
    _add_const_property(props, "density",              fluid["rho_ref_kg_m3"])
    _add_const_property(props, "viscosity",            fluid["viscosity_Pa_s"])
    _add_const_property(props, "specific_heat_capacity", fluid["cp_J_kgK"])
    _add_const_property(props, "thermal_conductivity", fluid["lambda_W_mK"])
    _add_const_property(props, "storage",              op["fluid_storage_1_per_Pa"])


def _add_phase_solid(phases, mat, op):
    ph = _se(phases, "phase")
    _se(ph, "type", "Solid")
    props = _se(ph, "properties")
    _add_const_property(props, "density",                mat["rho_s_kg_m3"])
    _add_const_property(props, "specific_heat_capacity", mat["cp_s_J_kgK"])
    _add_const_property(props, "thermal_conductivity",   mat["lambda_s_W_mK"])
    _add_const_property(props, "storage",                op["solid_storage_1_per_Pa"])


def _add_medium(media, mid, mat, fluid, op, disp):
    med = _se(media, "medium", id=mid)
    phases = _se(med, "phases")
    _add_phase_aqueous(phases, fluid, op)
    _add_phase_solid(phases, mat, op)
    props = _se(med, "properties")
    _add_const_property(props, "porosity",     mat["porosity"])
    _add_const_property(props, "permeability", mat["permeability_m2"])
    p = _se(props, "property"); _se(p, "name", "thermal_conductivity")
    _se(p, "type", "EffectiveThermalConductivityPorosityMixing")
    _add_const_property(props, "thermal_longitudinal_dispersivity", disp["alpha_L_m"])
    _add_const_property(props, "thermal_transversal_dispersivity",  disp["alpha_T_m"])
    _add_const_property(props, "storage", 0.0)


def _curve_xml(parent, name, t, v):
    c = _se(parent, "curve")
    _se(c, "name", name)
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v))


def _indent(elem, level=0):
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


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict, curves: dict) -> Path:
    prefix = cfg["output"]["prefix"]
    fluid  = cfg["fluid"]
    op     = cfg["operation"]
    init   = cfg["initial"]
    sol    = cfg["solver"]
    disp   = cfg["dispersion"]

    # Wärmequelle pro Volumen [W/m³] je Sonde
    h_bh = cfg["borehole"]["depth_bottom_m"] - cfg["borehole"]["depth_top_m"]
    V_bh = cfg["borehole"]["borehole_dx_m"] * cfg["borehole"]["borehole_dy_m"] * h_bh
    q_v = op["power_per_borehole_W"] / V_bh

    root = ET.Element("OpenGeoSysProject")

    meshes = _se(root, "meshes")
    mesh_keys = ["domain", "top", "bottom", "lateral"]
    mesh_keys += [L["name"] for L in cfg["layers"]]
    mesh_keys += [f"bh_{i:02d}" for i in range(mesh_files["_n_bh"])]
    for key in mesh_keys:
        if mesh_files.get(key):
            _se(meshes, "mesh", mesh_files[key])

    # Process
    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HT"); _se(proc, "type", "HT"); _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables"); _se(pv, "temperature", "T"); _se(pv, "pressure", "p")
    sv = _se(proc, "secondary_variables")
    _se(sv, "secondary_variable", internal_name="darcy_velocity", output_name="darcy_velocity")
    _se(proc, "specific_body_force", "0 0 0")

    # Media: MaterialID-Reihenfolge nach msh2vtu-Reindex = erst Schichten
    # (unten→oben, passend zu build_mesh), dann Sonden.
    media = _se(root, "media")
    layers_bot_up = list(reversed(cfg["layers"]))   # passend zu mesh-Tags
    for mid, L in enumerate(layers_bot_up):
        _add_medium(media, mid, L, fluid, op, disp)
    # Sondenmaterial = mittlere Schicht (didaktische Voreinstellung; hier
    # könnte auch ein „Filterkies"-Material hinterlegt werden).
    bh_mat = layers_bot_up[len(layers_bot_up) // 2]
    for i in range(mesh_files["_n_bh"]):
        _add_medium(media, len(layers_bot_up) + i, bh_mat, fluid, op, disp)

    # Time loop
    tl = _se(root, "time_loop")
    procs = _se(tl, "processes")
    p_ref = _se(procs, "process", ref="HT")
    _se(p_ref, "nonlinear_solver", "basic_picard")
    cc = _se(p_ref, "convergence_criterion")
    _se(cc, "type", "PerComponentDeltaX")
    _se(cc, "norm_type", "NORM2")
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

    # Parameters
    params = _se(root, "parameters")
    # T0: konstant oder tiefen­abhängig (geothermischer Gradient)
    gradient = init.get("geothermal_gradient_K_per_m", 0.0)
    if abs(gradient) > 1e-12:
        T_surf = init.get("T_surface_K", init["T_K"])
        z_tot  = sum(L["thickness_m"] for L in cfg["layers"])
        p_el = _se(params, "parameter")
        _se(p_el, "name", "T0")
        _se(p_el, "type", "Function")
        # Tiefe = z_total - z;  T0 = T_surf + gradient · Tiefe
        _se(p_el, "expression", f"{T_surf} + ({gradient:.6g})*({z_tot} - z)")
    else:
        _const_param(params, "T0", init["T_K"])
    _const_param(params, "p0",   init["p_Pa"])
    _const_param(params, "q_v_amp", q_v)
    _curve_scaled_param(params, "q_v_borehole", "cycle_q", "q_v_amp")

    cv = _se(root, "curves")
    t, v = curves["cycle_q"]
    _curve_xml(cv, "cycle_q", t, v)

    # Process variables
    pvars = _se(root, "process_variables")

    # T
    pv_T = _se(pvars, "process_variable")
    _se(pv_T, "name", "T"); _se(pv_T, "components", 1); _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs = _se(pv_T, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet")
        _se(bc, "parameter", "T0")
    sts = _se(pv_T, "source_terms")
    for i in range(mesh_files["_n_bh"]):
        st = _se(sts, "source_term")
        _se(st, "mesh", Path(mesh_files[f"bh_{i:02d}"]).stem)
        _se(st, "type", "Volumetric")
        _se(st, "parameter", "q_v_borehole")

    # p (statisch, nur Dirichlet Boundary)
    pv_p = _se(pvars, "process_variable")
    _se(pv_p, "name", "p"); _se(pv_p, "components", 1); _se(pv_p, "order", 1)
    _se(pv_p, "initial_condition", "p0")
    bcs = _se(pv_p, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet")
        _se(bc, "parameter", "p0")

    # Solvers
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
#  Run + CLI
# ======================================================================
def run_ogs(prj_path: Path) -> int:
    ogs_exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if not ogs_exe:
        print("ogs.exe nicht im PATH", file=sys.stderr); return 1
    cmd = [ogs_exe, str(prj_path), "-o", str(prj_path.parent)]
    print(">>", " ".join(cmd))
    return subprocess.call(cmd)


def main() -> int:
    ap = argparse.ArgumentParser(description="BTES 3D OGS demo")
    ap.add_argument("--no-mesh", action="store_true")
    ap.add_argument("--no-run",  action="store_true")
    args = ap.parse_args()

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print("[1/3] gmsh: Sondenfeld ...")
        build_mesh(CONFIG, out_dir)
        print(f"      {msh_path}")
        print("[2/3] msh2vtu: Konvertierung ...")
        mesh_files = convert_mesh(CONFIG, msh_path, out_dir)
    else:
        mesh_files = _mesh_files(CONFIG)

    print("[3/3] OGS-Projektdatei ...")
    curves = build_cycle_curves(CONFIG)
    prj_path = build_prj(CONFIG, out_dir, mesh_files, curves)
    print(f"      {prj_path}  (t_end = {curves['t_total']/DAY:.1f} d, {mesh_files['_n_bh']} Sonden)")

    if args.no_run:
        return 0
    print(">>> OGS starten")
    return run_ogs(prj_path)


if __name__ == "__main__":
    sys.exit(main())

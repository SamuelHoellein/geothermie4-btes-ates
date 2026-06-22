#!/usr/bin/env python3
"""
BTES 2D radialsymmetrisch — Einführungsbeispiel für OpenGeoSys 6.

Eine einzelne Erdwärmesonde auf der Symmetrieachse (r = 0). Das 3D-Problem
wird unter Annahme rotations­symmetrischer Lösung auf eine 2D-Aufgabe in
der (r, z)-Halbebene reduziert. Wärmetransfer rein durch Wärmeleitung
(quasi-impermeabler Boden).

VORTEIL: Bei gleicher Auflösung 100-1000x weniger Zellen als 3D
         -> Sekunden bis wenige Minuten Laufzeit.

VERWENDUNG
----------
    python btes_radial_2d.py             # Mesh + Sim + Plots
    python btes_radial_2d.py --no-run    # nur Setup
    python btes_radial_2d.py --no-mesh   # nur .prj
    python btes_radial_2d.py --no-plots  # ohne Auto-Plots
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



CONFIG: dict = {
    "domain": {
        "r_max_m":  50.0,
        "z_base_m":  0.0,
    },
    # ------------------------------------------------------------------
    # BODENSCHICHTUNG (Mehrschicht-API)
    # ------------------------------------------------------------------
    # Liste von Bodenschichten, Reihenfolge VON OBEN NACH UNTEN (wie in
    # der Natur). Jede Schicht hat ihre eigenen Material­werte. Die
    # Schicht­dicken addieren sich zur Gesamttiefe der Domain.
    # Auf eine einzelne Bodenklasse zurückzufallen, einfach nur einen
    # Eintrag in die Liste setzen.
    # ------------------------------------------------------------------
    "layers": [
        # name              thickness_m  permeability  porosity  rho_s  cp_s  lambda_s
        {"name": "cover",      "thickness_m":  5.0,
         "permeability_m2": 1.0e-15, "porosity": 0.35,
         "rho_s_kg_m3": 1900.0, "cp_s_J_kgK": 1500.0, "lambda_s_W_mK": 1.4},
        {"name": "bedrock",    "thickness_m": 80.0,
         "permeability_m2": 1.0e-18, "porosity": 0.20,
         "rho_s_kg_m3": 2700.0, "cp_s_J_kgK":  900.0, "lambda_s_W_mK": 2.5},
        {"name": "basement",   "thickness_m": 30.0,
         "permeability_m2": 1.0e-19, "porosity": 0.10,
         "rho_s_kg_m3": 2750.0, "cp_s_J_kgK":  850.0, "lambda_s_W_mK": 3.0},
    ],
    # ------------------------------------------------------------------
    # SONDEN-GEOMETRIE
    # ------------------------------------------------------------------
    # Sondenkopf-/Sondenfuß-Tiefe gemessen VON DER OBERFLÄCHE NACH UNTEN.
    # Die Sonde kann mehrere Schichten durchstoßen — sie wird als
    # zusammenhängendes Volumen aufgebaut.
    "borehole": {
        "r_borehole_m":    0.5,    # Sondenradius [m]
        "depth_top_m":     7.0,    # Sondenkopf, Tiefe unter Oberfläche [m]
        "depth_bottom_m": 83.0,    # Sondenfuß,  Tiefe unter Oberfläche [m]
    },
    "mesh": {
        "size_in_borehole_m":         0.2,
        "size_near_borehole_m":       0.8,
        "size_far_m":                 4.0,
        "borehole_size_radius_m":     3.0,
        "borehole_size_radius_far_m": 20.0,
    },
    # `materials` wird automatisch aus `layers` abgeleitet. Hier nur
    # vorhanden, falls Skripte „materials.soil"-Verweise erwarten — wird
    # als Fallback der ersten Schicht verwendet.
    "materials": {
        "soil": {
            "permeability_m2": 1.0e-18,
            "porosity":        0.20,
            "rho_s_kg_m3":  2700.0,
            "cp_s_J_kgK":    900.0,
            "lambda_s_W_mK":   2.5,
        },
    },
    "fluid": {
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         293.15,
        "beta_1_per_K":    0.0,
        "viscosity_Pa_s":  1.0e-3,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },
    "dispersion": {
        "alpha_L_m": 0.0,
        "alpha_T_m": 0.0,
    },
    "initial": {
        # Konstante Anfangstemperatur (verwendet, wenn geothermal_gradient_K_per_m = 0):
        "T_K":  283.15,
        "p_Pa": 0.0,
        # ------------------------------------------------------------
        # Geothermischer Tiefen­gradient (optional)
        # ------------------------------------------------------------
        # Wenn ≠ 0, wird die Anfangs­temperatur tiefen­abhängig gesetzt:
        #   T0(z) = T_surface_K + geothermal_gradient_K_per_m · Tiefe(z)
        # mit Tiefe(z) = z_total − z (interne z-Koordinate, 0 am Modellboden).
        # Typischer mittel­europäischer Wert: 0.03 K/m (≈ 3 K pro 100 m).
        # ------------------------------------------------------------
        "T_surface_K":               283.15,    # T an der Oberfläche [K]
        "geothermal_gradient_K_per_m": 0.0,     # 0.03 für realistischen Gradient
    },
    "operation": {
        "power_per_borehole_W":   2000.0,
        "fluid_storage_1_per_Pa": 4.5e-10,
        "solid_storage_1_per_Pa": 1.0e-10,
    },
    # ------------------------------------------------------------------
    # ZYKLEN
    # ------------------------------------------------------------------
    # Zwei Modi (alternativ):
    #
    # A) 4-Phasen-Zyklus (Default):
    #    Pro Zyklus 4 aufeinander folgende Phasen, jede mit fester Dauer:
    #      1) charge                  – Beladung (Wärme einspeisen)
    #      2) storage_after_charge    – Pause nach Beladung
    #      3) discharge               – Förderung (Wärme entziehen)
    #      4) storage_after_discharge – Pause nach Förderung
    #    Lade-/Förder-Leistung = operation.power_per_borehole_W (mit Vorzeichen ±).
    #
    # B) Monatsprofil (überschreibt A, falls aktiviert):
    #    Eine Liste von 12 Monats­leistungen [W]. Positiv = Wärme einspeisen,
    #    negativ = Wärme entziehen, 0 = Stillstand. Jeder Monat dauert
    #    365.25/12 ≈ 30.44 d. Die Sequenz wird n_cycles-mal wiederholt
    #    (n_cycles entspricht hier "Jahren").
    #    Setze monthly_power_W = [P_Jan, P_Feb, …, P_Dez] um diesen Modus
    #    zu aktivieren; auf None setzen für Modus A.
    #    operation.power_per_borehole_W dient dann als Referenz­leistung
    #    (Skalierung in OGS-Curves). Trage z. B. den Maximalwert ein.
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
        #   "monthly_power_W": [+2000, +2000, +1500, +500, 0, 0,
        #                         0, 0, -1500, -2500, -3000, -2500],
        "monthly_power_W":                  None,
    },
    "time": {
        "dt_seconds":            7 * 86400.0,
        "output_every_n_steps":  1,
        "gravity":               False,
    },
    "output": {
        "prefix":    "btes_radial_2d",
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
    """Liefert layers von OBEN nach UNTEN (wie in cfg) plus z-Grenzen.
    Rückgabe: list of dicts {name, thickness_m, ...material..., z_low, z_high},
    sortiert von unten (z_low minimal) nach oben (höchstes z).
    Plus z_top (= Oberfläche).
    """
    z_base = cfg["domain"].get("z_base_m", 0.0)
    raw    = list(cfg["layers"])     # top → bottom in CONFIG
    bot_up = list(reversed(raw))     # bottom → top
    z = z_base
    out = []
    for L in bot_up:
        z_low  = z
        z_high = z + float(L["thickness_m"])
        out.append({**L, "z_low": z_low, "z_high": z_high})
        z = z_high
    return out, z

def build_mesh(cfg: dict, out_dir: Path) -> Path:
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"

    r_max  = cfg["domain"]["r_max_m"]
    z_base = cfg["domain"]["z_base_m"]
    layers, z_top = _layer_stack(cfg)

    # Sondenposition: Tiefe von der Oberfläche nach unten
    r_bh         = cfg["borehole"]["r_borehole_m"]
    depth_top    = cfg["borehole"]["depth_top_m"]
    depth_bottom = cfg["borehole"]["depth_bottom_m"]
    z_bh_top     = z_top - depth_top
    z_bh_bot     = z_top - depth_bottom
    h_bh         = z_bh_top - z_bh_bot
    if h_bh <= 0:
        raise ValueError(f"borehole.depth_bottom_m ({depth_bottom}) muss > depth_top_m ({depth_top}) sein.")
    if z_bh_bot < z_base or z_bh_top > z_top:
        raise ValueError("Sondentiefe liegt außerhalb der Schichtdomäne.")

    m = cfg["mesh"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("btes_radial_2d")

    # Rectangle pro Schicht (von unten nach oben)
    layer_rects = []
    for L in layers:
        tag = gmsh.model.occ.addRectangle(0, L["z_low"], 0, r_max,
                                          L["z_high"] - L["z_low"])
        layer_rects.append(tag)
    r_b = gmsh.model.occ.addRectangle(0, z_bh_bot, 0, r_bh, h_bh)

    # Fragment: alle Schicht-Rechtecke mit der Sonde
    base, others = layer_rects[0], layer_rects[1:] + [r_b]
    gmsh.model.occ.fragment([(2, base)], [(2, t) for t in others])
    gmsh.model.occ.synchronize()

    # Klassifikation der entstehenden 2D-Flächen
    surf_bh    = []
    surf_layer = {i: [] for i in range(len(layers))}   # idx wie in layers (von unten)
    for dim, tag in gmsh.model.getEntities(2):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.occ.getBoundingBox(dim, tag)
        xc = 0.5 * (xmin + xmax); yc = 0.5 * (ymin + ymax)
        ext_x = xmax - xmin
        small_x = ext_x < 2 * r_bh + 1e-3
        if small_x and (z_bh_bot - 1e-3) <= yc <= (z_bh_top + 1e-3):
            surf_bh.append(tag)
            continue
        for i, L in enumerate(layers):
            if L["z_low"] - 1e-6 <= yc <= L["z_high"] + 1e-6:
                surf_layer[i].append(tag); break
    if not surf_bh:
        raise RuntimeError("Sondenfläche nicht gefunden.")

    # Randkanten
    edge_top, edge_bot, edge_far = [], [], []
    for dim, tag in gmsh.model.getEntities(1):
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.occ.getBoundingBox(dim, tag)
        if abs(ymin - z_top) < 1e-6 and abs(ymax - z_top) < 1e-6:
            edge_top.append(tag); continue
        if abs(ymin - z_base) < 1e-6 and abs(ymax - z_base) < 1e-6:
            edge_bot.append(tag); continue
        if abs(xmin - r_max) < 1e-6 and abs(xmax - r_max) < 1e-6:
            edge_far.append(tag); continue

    # Physical Groups: Reihenfolge entspricht MaterialID nach reindex
    # (zuerst alle Bodenschichten von unten nach oben, dann Sonde)
    pg_tag = 1
    for i, L in enumerate(layers):
        gmsh.model.addPhysicalGroup(2, surf_layer[i], tag=pg_tag, name=L["name"])
        pg_tag += 1
    gmsh.model.addPhysicalGroup(2, surf_bh, tag=pg_tag, name="bh_vol")
    gmsh.model.addPhysicalGroup(1, edge_top,  tag=100, name="top")
    gmsh.model.addPhysicalGroup(1, edge_bot,  tag=101, name="bottom")
    gmsh.model.addPhysicalGroup(1, edge_far,  tag=102, name="far")

    # Für die spätere Verfeinerung benutzen wir alle Sondenkanten
    surf_soil = [t for surf_list in surf_layer.values() for t in surf_list]

    # Mesh-Verfeinerung um die Sonde
    bh_edges = []
    for tag in surf_bh:
        for d, t in gmsh.model.getBoundary([(2, tag)], oriented=False):
            if d == 1: bh_edges.append(abs(t))
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList", list(set(bh_edges)))
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin", m["size_near_borehole_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax", m["size_far_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin", m["borehole_size_radius_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax", m["borehole_size_radius_far_m"])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    bh_points = []
    for tag in surf_bh:
        for d, t in gmsh.model.getBoundary([(2, tag)], recursive=True, oriented=False):
            if d == 0: bh_points.append((d, t))
    if bh_points:
        gmsh.model.mesh.setSize(bh_points, m["size_in_borehole_m"])

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(2)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


def convert_mesh(cfg: dict, msh_path: Path, out_dir: Path) -> dict:
    prefix = cfg["output"]["prefix"]
    msh2vtu(filename=msh_path, output_path=out_dir, output_prefix=prefix,
            dim=2, reindex=True, log_level="WARNING")
    files = {
        "domain": f"{prefix}_domain.vtu",
        "top":    f"{prefix}_physical_group_top.vtu",
        "bottom": f"{prefix}_physical_group_bottom.vtu",
        "far":    f"{prefix}_physical_group_far.vtu",
        "bh_vol": f"{prefix}_physical_group_bh_vol.vtu",
    }
    for L in cfg["layers"]:
        files[L["name"]] = f"{prefix}_physical_group_{L['name']}.vtu"
    return files


def build_cycle_curves(cfg: dict) -> dict:
    cyc       = cfg["cycles"]
    n         = cyc["n_cycles"]
    ramp      = max(60.0, cyc["ramp_days"] * DAY)
    P_nominal = cfg["operation"]["power_per_borehole_W"]

    # === Monatsprofil-Modus (überschreibt 4-Phasen-Logik) ===
    monthly = cyc.get("monthly_power_W")
    if monthly is not None:
        assert len(monthly) == 12, "cycles.monthly_power_W muss 12 Werte enthalten."
        month_dur = 365.25 / 12.0 * DAY      # ~30.44 d
        if P_nominal == 0:
            raise ValueError("operation.power_per_borehole_W muss > 0 sein (Referenzleistung).")
        times = [0.0]; vals = [0.0]; t_now = 0.0
        for _ in range(n):
            for P_month in monthly:
                q_rel = float(P_month) / P_nominal
                t_now += ramp; times.append(t_now); vals.append(q_rel)
                hold = max(0.0, month_dur - ramp)
                if hold > 0:
                    t_now += hold; times.append(t_now); vals.append(q_rel)
        t_now += ramp; times.append(t_now); vals.append(0.0)
        return {"t_total": t_now, "cycle_q": (np.array(times), np.array(vals))}

    # === 4-Phasen-Modus (Default) ===
    t_c   = cyc["charge_days"]                  * DAY
    t_sc  = cyc["storage_after_charge_days"]    * DAY
    t_d   = cyc["discharge_days"]               * DAY
    t_sd  = cyc["storage_after_discharge_days"] * DAY
    phases = [(t_c, +1.0), (t_sc, 0.0), (t_d, -1.0), (t_sd, 0.0)]
    times = [0.0]; vals = [0.0]; t_now = 0.0
    for _ in range(n):
        for dur, q in phases:
            if dur <= 0: continue
            t_now += ramp; times.append(t_now); vals.append(q)
            hold = max(0.0, dur - ramp)
            if hold > 0:
                t_now += hold; times.append(t_now); vals.append(q)
    t_now += ramp; times.append(t_now); vals.append(0.0)
    return {"t_total": t_now, "cycle_q": (np.array(times), np.array(vals))}


# === prj-Helpers (identisch zu ATES 2D radial) ===
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

    # 3D-äquivalentes Sondenvolumen (Zylinder)
    h_bh = cfg["borehole"]["depth_bottom_m"] - cfg["borehole"]["depth_top_m"]
    r_bh = cfg["borehole"]["r_borehole_m"]
    V_bh = np.pi * r_bh**2 * h_bh
    q_v = op["power_per_borehole_W"] / V_bh

    root = ET.Element("OpenGeoSysProject")

    meshes = _se(root, "meshes")
    # Domain + Rand­meshes + alle Schicht-Submeshes + Sondenfilter
    mesh_keys = ["domain", "top", "bottom", "far"]
    mesh_keys += [L["name"] for L in cfg["layers"]]
    mesh_keys += ["bh_vol"]
    for key in mesh_keys:
        _se(meshes, "mesh", mesh_files[key], axially_symmetric="true")

    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HT"); _se(proc, "type", "HT"); _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables"); _se(pv, "temperature", "T"); _se(pv, "pressure", "p")
    sv = _se(proc, "secondary_variables")
    _se(sv, "secondary_variable", internal_name="darcy_velocity", output_name="darcy_velocity")
    _se(proc, "specific_body_force", "0 0")

    media = _se(root, "media")
    # MaterialID-Reihenfolge entspricht den Physical-Group-Tags nach
    # ogstools-Reindex: erst alle Bodenschichten (Reihenfolge wie in
    # CONFIG["layers"], also TOP→BOTTOM in der Liste, was nach reversed()
    # in build_mesh die Tag-Reihenfolge unten→oben ergibt), dann Sonde.
    layers_bot_up = list(reversed(cfg["layers"]))  # passend zu mesh-Tags
    for mid, L in enumerate(layers_bot_up):
        _add_medium(media, mid, L, fluid, op, disp)
    # Sondenmaterial = mittlere Schicht (Voreinstellung); didaktisch
    # könnte hier auch ein „Filterkies"-Mat hinterlegt werden.
    bh_mat = layers_bot_up[len(layers_bot_up) // 2]
    _add_medium(media, len(layers_bot_up), bh_mat, fluid, op, disp)

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

    params = _se(root, "parameters")
    # T0: konstant oder tiefen­abhängig (geothermischer Gradient)
    gradient = init.get("geothermal_gradient_K_per_m", 0.0)
    if abs(gradient) > 1e-12:
        T_surf = init.get("T_surface_K", init["T_K"])
        z_tot  = sum(L["thickness_m"] for L in cfg["layers"])
        p_el = _se(params, "parameter")
        _se(p_el, "name", "T0")
        _se(p_el, "type", "Function")
        # Tiefe = z_total - y;  T0(z) = T_surf + gradient · Tiefe
        _se(p_el, "expression", f"{T_surf} + ({gradient:.6g})*({z_tot} - y)")
    else:
        _const_param(params, "T0", init["T_K"])
    _const_param(params, "p0",   init["p_Pa"])
    _const_param(params, "q_v_amp", q_v)
    _curve_param(params, "q_v_borehole", "cycle_q", "q_v_amp")

    cv = _se(root, "curves")
    t, v = curves["cycle_q"]; _curve_xml(cv, "cycle_q", t, v)

    pvars = _se(root, "process_variables")
    pv_T = _se(pvars, "process_variable")
    _se(pv_T, "name", "T"); _se(pv_T, "components", 1); _se(pv_T, "order", 1)
    _se(pv_T, "initial_condition", "T0")
    bcs = _se(pv_T, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T0")
    sts = _se(pv_T, "source_terms")
    st = _se(sts, "source_term")
    _se(st, "mesh", Path(mesh_files["bh_vol"]).stem)
    _se(st, "type", "Volumetric"); _se(st, "parameter", "q_v_borehole")

    pv_p = _se(pvars, "process_variable")
    _se(pv_p, "name", "p"); _se(pv_p, "components", 1); _se(pv_p, "order", 1)
    _se(pv_p, "initial_condition", "p0")
    bcs = _se(pv_p, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "p0")

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
        print("pyvista/matplotlib fehlt."); return

    figdir = out_dir / "figures"; figdir.mkdir(exist_ok=True)
    prefix = cfg["output"]["prefix"]
    files = sorted(out_dir.glob(f"{prefix}_ts_*_t_*.vtu"),
                   key=lambda p: int(re.search(r"_ts_(\d+)_", p.name).group(1)))
    if not files: return

    T0 = cfg["initial"]["T_K"]
    t_cycle = sum(cfg["cycles"][k] for k in
                  ("charge_days","storage_after_charge_days",
                   "discharge_days","storage_after_discharge_days"))
    # Sondenmitte in interner z-Koordinate (Oberfläche bei z_total, Boden bei 0).
    # `layers` ist eine Liste von Schichten -> Gesamttiefe = Summe der Dicken;
    # Sondentiefen werden von der Oberfläche nach unten gemessen.
    z_total = sum(L["thickness_m"] for L in cfg["layers"])
    bh      = cfg["borehole"]
    z_mid   = z_total - 0.5 * (bh["depth_top_m"] + bh["depth_bottom_m"])
    bh_half = 0.5 * (bh["depth_bottom_m"] - bh["depth_top_m"])

    def file_for(day):
        return min(files, key=lambda p: abs(int(re.search(r"_t_(\d+)", p.name).group(1))/DAY - day))

    # Adaptive T-Range
    T_hi_glob, T_lo_glob = T0 + 0.1, T0 - 0.1
    for f in files[::max(1, len(files)//12)]:
        m = pv.read(f)
        T_hi_glob = max(T_hi_glob, float(m["T"].max()))
        T_lo_glob = min(T_lo_glob, float(m["T"].min()))
    span = max(T_hi_glob - T0, T0 - T_lo_glob, 0.5)
    T_hi, T_lo = T0 + span, T0 - span

    keypoints = [
        (1, "Tag 1"),
        (cfg["cycles"]["charge_days"], "Ende Beladung"),
        (cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"], "Ende Pause 1"),
        (cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"]
         + cfg["cycles"]["discharge_days"], "Ende Förderung"),
        (t_cycle, "Ende Pause 2"),
    ]
    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=(1800, 600), shape=(1, len(keypoints)))
    for i, (day, label) in enumerate(keypoints):
        m = pv.read(file_for(day))
        plotter.subplot(0, i)
        plotter.add_mesh(m, scalars="T", cmap="coolwarm", clim=[T_lo, T_hi],
                         show_scalar_bar=(i == len(keypoints) - 1),
                         scalar_bar_args={"title": "T [K]"} if i == len(keypoints)-1 else None)
        plotter.add_text(f"{label} (Tag {day:.0f})", font_size=9, position="upper_edge")
        plotter.view_xy()
    plotter.screenshot(str(figdir / "T_field_2D.png"))
    print("  saved T_field_2D.png")

    # T(t) Plot
    r_bh = cfg["borehole"]["r_borehole_m"]
    probes = {
        f"Sonde (r={r_bh:.1f})":    (r_bh*0.5, z_mid, 0),
        "r=2 m":                    (2.0,  z_mid, 0),
        "r=5 m":                    (5.0,  z_mid, 0),
        "r=15 m":                   (15.0, z_mid, 0),
        "Cover 2 m über BHE":       (1.0,  z_mid + bh_half + 2, 0),
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
        t2 = t0c + cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"]
        ax.axvspan(t2, t2 + cfg["cycles"]["discharge_days"], alpha=0.08, color="blue")
    for k, T in series.items():
        ax.plot(times, T, lw=1.5, label=k)
    ax.axhline(T0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("T [K]")
    ax.set_title("BTES 2D radial — T(t) an Beobachtungspunkten")
    ax.legend(loc="upper right", fontsize=9, ncol=2); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(figdir / "T_vs_time.png", dpi=130); plt.close(fig)
    print("  saved T_vs_time.png")

    # Energiebilanz
    domain = pv.read(out_dir / f"{prefix}_domain.vtu")
    mid = domain.cell_data["MaterialIDs"]
    cs = domain.compute_cell_sizes()
    cell_area = cs["Area"]
    ctrs = domain.cell_centers()
    r_centroid = ctrs.points[:, 0]
    cell_vol_axi = 2 * np.pi * r_centroid * cell_area
    soil = cfg["materials"]["soil"]
    rcp = soil["porosity"]*1000*4180 + (1 - soil["porosity"])*soil["rho_s_kg_m3"]*soil["cp_s_J_kgK"]

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
        ax.axvspan(t0c + cfg["cycles"]["charge_days"], t0c + cfg["cycles"]["charge_days"]+cfg["cycles"]["storage_after_charge_days"],
                   alpha=0.08, color="gray", label="Pause 1" if cyc==0 else None)
        t2 = t0c + cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"]
        ax.axvspan(t2, t2 + cfg["cycles"]["discharge_days"], alpha=0.08, color="blue", label="Förderung" if cyc==0 else None)
    ax.plot(times, E_tot, "k-", lw=2.0, label="Gespeichert (axi-symm.)")
    ax.axhline(0, color="k", lw=0.5)
    i_p1 = int(np.argmin(np.abs(times - (cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"]))))
    i_f  = int(np.argmin(np.abs(times - (cfg["cycles"]["charge_days"] + cfg["cycles"]["storage_after_charge_days"] + cfg["cycles"]["discharge_days"]))))
    if E_tot[i_p1] > 0:
        rec = (E_tot[i_p1] - E_tot[i_f]) / E_tot[i_p1] * 100
        ax.annotate(f"Recovery Z1: {rec:.1f}%", xy=(times[i_f], E_tot[i_f]),
                    xytext=(times[i_f]*0.5, E_tot[i_p1]*0.6),
                    fontsize=10, arrowprops=dict(arrowstyle="->", color="gray"))
    ax.set_xlabel("Zeit [Tage]"); ax.set_ylabel("Energie über T0 [GJ]")
    ax.set_title("BTES 2D radial — Energiebilanz")
    ax.grid(alpha=0.3); ax.legend(loc="upper right", fontsize=9, ncol=2)
    fig.tight_layout()
    fig.savefig(figdir / "energy_balance.png", dpi=130); plt.close(fig)
    print("  saved energy_balance.png")


def main() -> int:
    ap = argparse.ArgumentParser(description="BTES 2D radial demo")
    ap.add_argument("--no-mesh",  action="store_true")
    ap.add_argument("--no-run",   action="store_true")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print("[1/4] gmsh: 2D Radialmesh ...")
        build_mesh(CONFIG, out_dir)
        print(f"      -> {msh_path}")
        print("[2/4] msh2vtu ...")
        mesh_files = convert_mesh(CONFIG, msh_path, out_dir)
    else:
        mesh_files = {
            "domain": f"{prefix}_domain.vtu",
            "top":    f"{prefix}_physical_group_top.vtu",
            "bottom": f"{prefix}_physical_group_bottom.vtu",
            "far":    f"{prefix}_physical_group_far.vtu",
            "bh_vol": f"{prefix}_physical_group_bh_vol.vtu",
        }
        for L in CONFIG["layers"]:
            mesh_files[L["name"]] = f"{prefix}_physical_group_{L['name']}.vtu"
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

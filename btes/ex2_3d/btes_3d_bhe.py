#!/usr/bin/env python3
# coding: utf-8
"""
BTES 3D — Variante mit OpenGeoSys-Modul HEAT_TRANSPORT_BHE.

STATUS: Lauffähig. Mesh-/PRJ-Generierung und ein end-to-end-OGS-Lauf
sind mit der Default-Konfiguration (3x3-Feld, Typ 1U, ein Lade-/Entlade-
Zyklus -> 53 Zeitschritte) verifiziert. Quantitativ noch NICHT validiert:
vor produktivem Einsatz gegen die OGS-Benchmarks "BHE_1U" / "BHE_2U" /
"BHE_CXA_CXC" gegenprüfen.

Unterschiede zur Variante `btes_3d.py` (vereinfachte HT-Volumen­quelle):

- Sonden werden als **1D-Linien­elemente** im 3D-Mesh eingebettet,
  nicht als 3D-Boxen mit Wärme­quell­term.
- Prozess­typ ist `HEAT_TRANSPORT_BHE` (eigenes OGS-Modul) mit
  expliziter U-Rohr-Geometrie, Grout, Refrigerant und Strömungs­regelung.
- Studierende konfigurieren BHE-Typ (1U, 2U, CXA, CXC), Innen/Außen-
  Rohrdurchmesser, Grout-λ, Refrigerant-Massenstrom und Inlet-T-
  Steuerung über `CONFIG["bhe"]`.

Doku OGS:  https://www.opengeosys.org/docs/processes/heat-transport/bhe/
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


def msh2vtu(filename, output_path, output_prefix, dim, reindex=True, log_level="WARNING"):
    import ogstools as ot
    from pathlib import Path as _P
    meshes = ot.Meshes.from_gmsh(filename=str(filename), dim=dim,
                                  reindex=reindex, log=False)
    output_path = _P(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        fname = (f"{output_prefix}_domain.vtu" if name == "domain"
                 else f"{output_prefix}_physical_group_{name}.vtu")
        mesh.save(str(output_path / fname), binary=True)


# ======================================================================
#  CONFIG
# ======================================================================
CONFIG: dict = {
    "domain": {
        "size_x_m":  80.0,
        "size_y_m":  80.0,
        "z_base_m":   0.0,
    },
    # Schichten (von OBEN nach UNTEN) — vereinheitlicht mit btes_radial_2d.py
    # und btes_3d.py. Das Boden-Mesh ist echt geschichtet (jede Schicht eine
    # eigene MaterialID). Die eingebetteten 1D-BHE-Linien werden beim
    # Fragmentieren an den Schichtgrenzen geteilt; quert eine BHE mehrere
    # Schichten, nutzt ihr 1D-Element die mittlere Schicht als repräsentatives
    # umgebendes Material.
    "layers": [
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
    # Sonde: Tiefe von der Oberfläche nach unten (wie in den anderen Übungen)
    "borehole": {
        "depth_top_m":     7.0,           # Sondenkopf, Tiefe unter Oberfläche [m]
        "depth_bottom_m": 82.0,           # Sondenfuß,  Tiefe unter Oberfläche [m]
    },
    "field": {
        # Sondenfeld — N x N Raster ODER explizite Positionsliste
        "n_x":           3,
        "n_y":           3,
        "spacing_m":     5.0,
        "positions":     None,            # bei nicht-None überschreibt es das Raster
    },
    "mesh": {
        "size_in_borehole_m":      0.4,
        "size_near_field_m":       1.5,
        "size_far_m":             15.0,
        "field_size_radius_m":     8.0,
        "field_size_radius_far_m": 30.0,
    },
    "fluid": {
        "rho_ref_kg_m3":  1000.0,
        "T_ref_K":         293.15,
        "viscosity_Pa_s":  1.0e-3,
        "cp_J_kgK":        4180.0,
        "lambda_W_mK":     0.6,
    },
    # ------------------------------------------------------------------
    # BHE-MODELL — HIER FÜR STUDIERENDE
    # ------------------------------------------------------------------
    "bhe": {
        # Typ: "1U"  – einfach-U-Rohr
        #      "2U"  – doppel-U-Rohr (zwei parallele Schleifen)
        #      "CXA" – coaxial, Vorlauf außen
        #      "CXC" – coaxial, Vorlauf zentral
        "type": "1U",
        "borehole": {
            "diameter_m": 0.15,           # Bohrdurchmesser
        },
        "pipes": {
            "diameter_outer_m":               0.032,   # Rohr-Außendurchmesser
            "wall_thickness_m":               0.003,
            "wall_thermal_conductivity_W_mK": 0.4,
            "distance_between_pipes_m":       0.05,    # Achsabstand Vor-/Rücklauf
            "longitudinal_dispersion_length_m": 0.001,
        },
        "grout": {
            "density_kg_m3":                  2190.0,
            "porosity":                          0.0,
            "specific_heat_capacity_J_kgK":   1735.0,
            "thermal_conductivity_W_mK":         2.3,
        },
        "refrigerant": {
            "density_kg_m3":                  1052.0,  # Wasser-Glykol-Mischung
            "viscosity_Pa_s":                 0.0052,
            "specific_heat_capacity_J_kgK":   3795.0,
            "thermal_conductivity_W_mK":         0.48,
            "reference_temperature_K":         295.15,
        },
        "control": {
            # Typ: "FixedPowerConstantFlow" — feste Leistung + konstanter Flow.
            # Hinweis: Wird cycles.monthly_power_W gesetzt (Lastprofil), schaltet
            # das Skript automatisch auf das OGS-Steuerungsmodell
            # "PowerCurveConstantFlow" um (zeitabhängige Leistungskurve bei
            # konstantem Durchfluss) und ignoriert power_W.
            "type":             "FixedPowerConstantFlow",
            "power_W":           2000.0,    # + = laden, - = fördern
            "flow_rate_kg_s":    0.2,
        },
    },
    "initial": {
        "T_K":         283.15,
        "p_Pa":          0.0,
        "T_surface_K": 283.15,
        "geothermal_gradient_K_per_m": 0.0,
    },
    "cycles": {
        "n_cycles":                     1,
        "charge_days":                 91.25,
        "storage_after_charge_days":   91.25,
        "discharge_days":              91.25,
        "storage_after_discharge_days":91.25,
        "ramp_days":                    7.0,
        # --- Modus B: Monatsprofil (auf None für den 4-Phasen-Default) ---
        # Liste von 12 Monatsleistungen [W] je Sonde (+ = laden, - = fördern,
        # 0 = Stillstand). Aktiv ≠ None: dann wird über das OGS-Steuerungsmodell
        # PowerCurveConstantFlow eine zeitabhängige Leistungskurve eingespeist.
        # Jeder Monat dauert 365.25/12 ≈ 30.44 d, n_cycles = Anzahl Jahre.
        # Beispiel (Beladung Sommer, Förderung Winter):
        #   "monthly_power_W": [+2000, +2000, +1500, +500, 0, 0,
        #                         0, 0, -1500, -2500, -3000, -2500],
        "monthly_power_W":              None,
    },
    "time": {
        "dt_seconds":          7 * 86400.0,
        "output_every_n_steps": 1,
        "gravity":              False,
    },
    "output": {
        "prefix":    "btes_3d_bhe",
        "out_dir":   "out",
        "variables": ["temperature_soil"],  # plus temperature_BHE_i pro BHE
    },
    "solver": {
        "linear_tol":      1.0e-12,
        "linear_iter":     10000,
        "nonlinear_iter":  20,
    },
}

DAY = 86400.0


def build_power_curve(cfg: dict):
    """Leistungskurve [W] aus cycles.monthly_power_W für PowerCurveConstantFlow.

    Liefert (t_total, (times, vals)) mit absoluten Leistungen [W] je Sonde
    (gleiches Vorzeichen wie bhe.control.power_W: + = laden, - = fördern).
    Jeder Monat dauert 365.25/12 d, mit sanften Rampen (ramp_days) zwischen den
    Monaten. Gibt None zurück, wenn kein Monatsprofil gesetzt ist.
    """
    monthly = cfg["cycles"].get("monthly_power_W")
    if monthly is None:
        return None
    assert len(monthly) == 12, "cycles.monthly_power_W muss 12 Werte enthalten."
    n         = cfg["cycles"]["n_cycles"]
    ramp      = max(60.0, cfg["cycles"]["ramp_days"] * DAY)
    month_dur = 365.25 / 12.0 * DAY      # ~30.44 d
    times = [0.0]; vals = [0.0]; t_now = 0.0
    for _ in range(n):
        for P_month in monthly:
            P = float(P_month)
            t_now += ramp; times.append(t_now); vals.append(P)
            hold = max(0.0, month_dur - ramp)
            if hold > 0.0:
                t_now += hold; times.append(t_now); vals.append(P)
    t_now += ramp; times.append(t_now); vals.append(0.0)
    return t_now, (np.array(times), np.array(vals))


# ======================================================================
#  Hilfen
# ======================================================================
def _bhe_positions(cfg: dict) -> list[tuple[float, float]]:
    """Liste der (x, y)-Position aller BHEs."""
    fld = cfg["field"]
    if fld.get("positions"):
        return [tuple(p) for p in fld["positions"]]
    nx, ny, sp = fld["n_x"], fld["n_y"], fld["spacing_m"]
    xs = (np.arange(nx) - (nx - 1) / 2) * sp
    ys = (np.arange(ny) - (ny - 1) / 2) * sp
    return [(float(x), float(y)) for y in ys for x in xs]


def _z_for_depth(cfg: dict, depth: float) -> float:
    """Interner z-Wert für eine Tiefe (gemessen von der Oberfläche)."""
    z_base = cfg["domain"]["z_base_m"]
    z_top  = z_base + sum(float(L["thickness_m"]) for L in cfg["layers"])
    return z_top - depth


def _layer_stack(cfg: dict):
    """Schichtliste (CONFIG: oben→unten) in z-Grenzen umrechnen.
    Rückgabe: (layers_bottom_up mit z_low/z_high, z_top)."""
    z_base = cfg["domain"].get("z_base_m", 0.0)
    bot_up = list(reversed(list(cfg["layers"])))
    z = z_base
    out = []
    for L in bot_up:
        h = float(L["thickness_m"])
        out.append({**L, "z_low": z, "z_high": z + h})
        z += h
    return out, z


# ======================================================================
#  Mesh — gmsh
# ======================================================================
def build_mesh(cfg: dict, out_dir: Path) -> Path:
    msh_path = out_dir / f"{cfg['output']['prefix']}.msh"
    Lx = cfg["domain"]["size_x_m"]
    Ly = cfg["domain"]["size_y_m"]
    z_base = cfg["domain"]["z_base_m"]
    layers, z_top = _layer_stack(cfg)

    bhe_pos = _bhe_positions(cfg)
    z_bhe_top = z_top - cfg["borehole"]["depth_top_m"]
    z_bhe_bot = z_top - cfg["borehole"]["depth_bottom_m"]
    m = cfg["mesh"]
    x0, y0 = -Lx / 2.0, -Ly / 2.0

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("btes_3d_bhe")

    # Eine Box pro Schicht + eine 1D-Linie pro BHE. Das gemeinsame
    # Fragmentieren macht die Schichten konform UND imprintet die Linien
    # als interne Kanten (an den Schichtgrenzen geteilt) in die Volumen —
    # die Mesh-Knoten sitzen damit exakt auf den BHE-Achsen.
    layer_boxes = [gmsh.model.occ.addBox(x0, y0, L["z_low"], Lx, Ly,
                                         L["z_high"] - L["z_low"]) for L in layers]
    lines = []
    for x, y in bhe_pos:
        p_top = gmsh.model.occ.addPoint(x, y, z_bhe_top)
        p_bot = gmsh.model.occ.addPoint(x, y, z_bhe_bot)
        lines.append(gmsh.model.occ.addLine(p_top, p_bot))
    gmsh.model.occ.fragment([(3, b) for b in layer_boxes], [(1, l) for l in lines])
    gmsh.model.occ.synchronize()

    # Schichtvolumen nach z-Mitte zuordnen
    vol_layer = {i: [] for i in range(len(layers))}
    for dim, tag in gmsh.model.getEntities(3):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (bb[2] + bb[5])
        for i, L in enumerate(layers):
            if L["z_low"] - 1e-6 <= zc <= L["z_high"] + 1e-6:
                vol_layer[i].append(tag); break

    # BHE-Liniensegmente (an Schichtgrenzen geteilt) nach (x,y) gruppieren
    pos = np.array(bhe_pos)
    seg_bhe = {i: [] for i in range(len(bhe_pos))}
    all_segs = []
    for dim, tag in gmsh.model.getEntities(1):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        dx, dy, dz = bb[3] - bb[0], bb[4] - bb[1], bb[5] - bb[2]
        xc, yc, zc = 0.5*(bb[0]+bb[3]), 0.5*(bb[1]+bb[4]), 0.5*(bb[2]+bb[5])
        vertical = dx < 1e-6 and dy < 1e-6 and dz > 1e-9
        if vertical and (z_bhe_bot - 1e-3) <= zc <= (z_bhe_top + 1e-3):
            d = np.hypot(pos[:, 0] - xc, pos[:, 1] - yc)
            j = int(np.argmin(d))
            if d[j] < 1e-3:
                seg_bhe[j].append(tag); all_segs.append(tag)
    missing = [i for i, s in seg_bhe.items() if not s]
    if missing:
        raise RuntimeError(f"BHE ohne Liniensegment gefunden: {missing}")

    # Physical Groups: Schichten (unten→oben) -> MaterialID 0..L-1,
    # danach BHE-Linien -> MaterialID L..L+N-1
    pg = 1
    for i, L in enumerate(layers):
        gmsh.model.addPhysicalGroup(3, vol_layer[i], tag=pg, name=L["name"]); pg += 1
    for i in range(len(bhe_pos)):
        gmsh.model.addPhysicalGroup(1, seg_bhe[i], tag=pg, name=f"bhe_{i:02d}"); pg += 1

    # Domänenränder
    top_faces, bot_faces, lat_faces = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (bb[2] + bb[5])
        if (bb[3] - bb[0]) >= 0.9 * Lx and abs(zc - z_top) < 1e-6:
            top_faces.append(tag)
        elif (bb[3] - bb[0]) >= 0.9 * Lx and abs(zc - z_base) < 1e-6:
            bot_faces.append(tag)
        else:
            on_outer = (abs(bb[0]-x0) < 1e-6 and abs(bb[3]-x0) < 1e-6) \
                       or (abs(bb[0]-(x0+Lx)) < 1e-6 and abs(bb[3]-(x0+Lx)) < 1e-6) \
                       or (abs(bb[1]-y0) < 1e-6 and abs(bb[4]-y0) < 1e-6) \
                       or (abs(bb[1]-(y0+Ly)) < 1e-6 and abs(bb[4]-(y0+Ly)) < 1e-6)
            if on_outer:
                lat_faces.append(tag)
    gmsh.model.addPhysicalGroup(2, top_faces, tag=200, name="top")
    gmsh.model.addPhysicalGroup(2, bot_faces, tag=201, name="bottom")
    if lat_faces:
        gmsh.model.addPhysicalGroup(2, lat_faces, tag=202, name="lateral")

    # Verfeinerung um die BHE-Linien
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList", all_segs)
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField",  f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin",  m["size_near_field_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax",  m["size_far_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin",  m["field_size_radius_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax",  m["field_size_radius_far_m"])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)
    bhe_pts = []
    for tag in all_segs:
        for d, t in gmsh.model.getBoundary([(1, tag)], oriented=False):
            if d == 0:
                bhe_pts.append((0, t))
    if bhe_pts:
        gmsh.model.mesh.setSize(list(set(bhe_pts)), m["size_in_borehole_m"])

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


# ======================================================================
#  PRJ — HEAT_TRANSPORT_BHE
# ======================================================================
def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None: el.text = str(text)
    return el

def _const_prop(parent, name, value):
    p = _se(parent, "property"); _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)

def _indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip(): elem.text = i + "  "
        for child in elem: _indent(child, level + 1)
        if not child.tail or not child.tail.strip(): child.tail = i
    if level and (not elem.tail or not elem.tail.strip()): elem.tail = i


def _bhe_xml(parent, cfg: dict, n_bhe: int) -> None:
    """Erzeugt <borehole_heat_exchangers>-Block für N gleichartige BHEs."""
    bhe = cfg["bhe"]
    L_bhe = cfg["borehole"]["depth_bottom_m"] - cfg["borehole"]["depth_top_m"]
    bhes = _se(parent, "borehole_heat_exchangers")
    for i in range(n_bhe):
        b = _se(bhes, "borehole_heat_exchanger")
        _se(b, "type", bhe["type"])
        # Strömungs-/Temperatur-Steuerung
        ftc = _se(b, "flow_and_temperature_control")
        # OGS erwartet einen VOLUMETRISCHEN Durchfluss [m^3/s], nicht kg/s.
        # Umrechnung ueber die Refrigerant-Dichte (sonst ~1000x zu hoch ->
        # Refrigerant-DeltaT ~ 0 -> keine Waermekopplung).
        _flow_vol = bhe["control"]["flow_rate_kg_s"] / bhe["refrigerant"]["density_kg_m3"]
        if cfg["cycles"].get("monthly_power_W") is not None:
            # Modus B: zeitabhängige Leistungskurve bei konstantem Durchfluss.
            _se(ftc, "type", "PowerCurveConstantFlow")
            _se(ftc, "power_curve", "power_curve")
            _se(ftc, "flow_rate", _flow_vol)
        else:
            _se(ftc, "type", bhe["control"]["type"])
            if bhe["control"]["type"] == "FixedPowerConstantFlow":
                _se(ftc, "power",     bhe["control"]["power_W"])
                _se(ftc, "flow_rate", _flow_vol)
        # Bohrloch
        bh = _se(b, "borehole")
        _se(bh, "length",   L_bhe)
        _se(bh, "diameter", bhe["borehole"]["diameter_m"])
        # Rohre
        pipes = _se(b, "pipes")
        for pname in ("inlet", "outlet"):
            p = _se(pipes, pname)
            _se(p, "diameter",                  bhe["pipes"]["diameter_outer_m"])
            _se(p, "wall_thickness",            bhe["pipes"]["wall_thickness_m"])
            _se(p, "wall_thermal_conductivity", bhe["pipes"]["wall_thermal_conductivity_W_mK"])
        _se(pipes, "distance_between_pipes",         bhe["pipes"]["distance_between_pipes_m"])
        _se(pipes, "longitudinal_dispersion_length", bhe["pipes"]["longitudinal_dispersion_length_m"])
        # Grout
        gr = _se(b, "grout")
        _se(gr, "density",                  bhe["grout"]["density_kg_m3"])
        _se(gr, "porosity",                 bhe["grout"]["porosity"])
        _se(gr, "specific_heat_capacity",   bhe["grout"]["specific_heat_capacity_J_kgK"])
        _se(gr, "thermal_conductivity",     bhe["grout"]["thermal_conductivity_W_mK"])
        # Refrigerant
        rf = _se(b, "refrigerant")
        _se(rf, "density",                  bhe["refrigerant"]["density_kg_m3"])
        _se(rf, "viscosity",                bhe["refrigerant"]["viscosity_Pa_s"])
        _se(rf, "specific_heat_capacity",   bhe["refrigerant"]["specific_heat_capacity_J_kgK"])
        _se(rf, "thermal_conductivity",     bhe["refrigerant"]["thermal_conductivity_W_mK"])
        _se(rf, "reference_temperature",    bhe["refrigerant"]["reference_temperature_K"])


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict) -> Path:
    prefix = cfg["output"]["prefix"]
    fluid, init, sol = cfg["fluid"], cfg["initial"], cfg["solver"]
    n_bhe = len(_bhe_positions(cfg))

    # Simulationsdauer: Monatsprofil (Modus B) oder 4-Phasen-Zyklus (Modus A).
    power_curve = build_power_curve(cfg)
    if power_curve is not None:
        t_end = power_curve[0]
    else:
        t_end = cfg["cycles"]["n_cycles"] * sum(cfg["cycles"][k] for k in
            ("charge_days","storage_after_charge_days",
             "discharge_days","storage_after_discharge_days")) * DAY
    n_steps = int(t_end // cfg["time"]["dt_seconds"]) + 1

    root = ET.Element("OpenGeoSysProject")
    meshes = _se(root, "meshes")
    for key in ("domain", "top", "bottom", "lateral"):
        _se(meshes, "mesh", mesh_files[key])
    for i in range(n_bhe):
        _se(meshes, "mesh", mesh_files[f"bhe_{i:02d}"])

    # Process
    processes = _se(root, "processes")
    proc = _se(processes, "process")
    _se(proc, "name", "HeatTransportBHE")
    _se(proc, "type", "HEAT_TRANSPORT_BHE")
    _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables")
    _se(pv, "process_variable", "temperature_soil")
    for i in range(n_bhe):
        _se(pv, "process_variable", f"temperature_BHE{i+1}")
    _bhe_xml(proc, cfg, n_bhe)

    # Media — Boden + ein Material pro BHE (Identifizierung über Mesh-MaterialID)
    media = _se(root, "media")
    # MaterialIDs nach msh2vtu-Reindex: erst Schichten (unten→oben) -> 0..L-1,
    # dann BHE-Linien -> L..L+N-1. BHE-Linien queren mehrere Schichten und
    # bekommen die mittlere Schicht als repräsentatives umgebendes Material.
    layers_bot_up = list(reversed(cfg["layers"]))
    bh_mat = layers_bot_up[len(layers_bot_up) // 2]
    mats = layers_bot_up + [bh_mat] * n_bhe
    for mid, soil in enumerate(mats):
        med = _se(media, "medium", id=mid)
        phs = _se(med, "phases")
        # Aqueous-Liquid
        ph = _se(phs, "phase"); _se(ph, "type", "AqueousLiquid")
        pp = _se(ph, "properties")
        _const_prop(pp, "density",                fluid["rho_ref_kg_m3"])
        _const_prop(pp, "viscosity",              fluid["viscosity_Pa_s"])
        _const_prop(pp, "specific_heat_capacity", fluid["cp_J_kgK"])
        _const_prop(pp, "thermal_conductivity",   fluid["lambda_W_mK"])
        # Keine Grundwasserströmung -> Phasengeschwindigkeit = Nullvektor.
        # HEAT_TRANSPORT_BHE verlangt diese Eigenschaft fuer den Advektionsterm.
        _const_prop(pp, "phase_velocity",         "0 0 0")
        # Solid
        ph = _se(phs, "phase"); _se(ph, "type", "Solid")
        pp = _se(ph, "properties")
        _const_prop(pp, "density",                soil["rho_s_kg_m3"])
        _const_prop(pp, "specific_heat_capacity", soil["cp_s_J_kgK"])
        _const_prop(pp, "thermal_conductivity",   soil["lambda_s_W_mK"])
        # Mediumeigenschaften
        props = _se(med, "properties")
        _const_prop(props, "porosity",     soil["porosity"])
        _const_prop(props, "permeability", soil["permeability_m2"])
        _const_prop(props, "thermal_conductivity",
                     soil["porosity"]*fluid["lambda_W_mK"] +
                     (1-soil["porosity"])*soil["lambda_s_W_mK"])
        _const_prop(props, "storage", 0.0)

    # Time loop
    tl = _se(root, "time_loop")
    procs = _se(tl, "processes")
    pref = _se(procs, "process", ref="HeatTransportBHE")
    _se(pref, "nonlinear_solver", "basic_picard")
    conv = _se(pref, "convergence_criterion")
    _se(conv, "type", "DeltaX"); _se(conv, "norm_type", "NORM2")
    _se(conv, "reltol", sol.get("rel_tol_T", 1e-4))
    td = _se(pref, "time_discretization"); _se(td, "type", "BackwardEuler")
    ts = _se(pref, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0)
    _se(ts, "t_end", t_end)
    tsteps = _se(ts, "timesteps")
    pair = _se(tsteps, "pair")
    _se(pair, "repeat", n_steps); _se(pair, "delta_t", cfg["time"]["dt_seconds"])
    out = _se(tl, "output")
    _se(out, "type", "VTK"); _se(out, "prefix", prefix)
    out_vars = _se(out, "variables")
    _se(out_vars, "variable", "temperature_soil")
    for i in range(n_bhe):
        _se(out_vars, "variable", f"temperature_BHE{i+1}")
    out_steps = _se(out, "timesteps"); pair = _se(out_steps, "pair")
    _se(pair, "repeat", n_steps)
    _se(pair, "each_steps", cfg["time"]["output_every_n_steps"])

    # Parameters + IC + BCs
    params = _se(root, "parameters")
    gradient = init.get("geothermal_gradient_K_per_m", 0.0)
    p_T0 = _se(params, "parameter")
    _se(p_T0, "name", "T0")
    if abs(gradient) > 1e-12:
        T_surf = init.get("T_surface_K", init["T_K"])
        z_tot  = sum(float(L["thickness_m"]) for L in cfg["layers"])
        _se(p_T0, "type", "Function")
        _se(p_T0, "expression", f"{T_surf} + ({gradient:.6g})*({z_tot} - z)")
    else:
        _se(p_T0, "type", "Constant"); _se(p_T0, "value", init["T_K"])
    p_T0_bhe = _se(params, "parameter")
    _se(p_T0_bhe, "name", "T0_BHE")
    _se(p_T0_bhe, "type", "Constant")
    # Initial-T im Pipe (Vor/Rück + Grout-Knoten) — bei 1U: 4 Werte
    _se(p_T0_bhe, "value", f"{init['T_K']} {init['T_K']} {init['T_K']} {init['T_K']}")

    pvars = _se(root, "process_variables")
    pvs = _se(pvars, "process_variable")
    _se(pvs, "name", "temperature_soil")
    _se(pvs, "components", 1); _se(pvs, "order", 1)
    _se(pvs, "initial_condition", "T0")
    bcs = _se(pvs, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T0")
    for i in range(n_bhe):
        pvb = _se(pvars, "process_variable")
        _se(pvb, "name", f"temperature_BHE{i+1}")
        _se(pvb, "components", 4 if cfg["bhe"]["type"] in ("1U", "CXA", "CXC") else 8)
        _se(pvb, "order", 1)
        _se(pvb, "initial_condition", "T0_BHE")
        _se(pvb, "boundary_conditions")   # leer — über BHE-Block geregelt

    # Nonlinear solver
    nls = _se(root, "nonlinear_solvers")
    n = _se(nls, "nonlinear_solver")
    _se(n, "name", "basic_picard"); _se(n, "type", "Picard")
    _se(n, "max_iter", sol["nonlinear_iter"])
    _se(n, "linear_solver", "general_linear_solver")
    lss = _se(root, "linear_solvers")
    ls = _se(lss, "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen"); _se(eig, "solver_type", "BiCGSTAB")
    _se(eig, "precon_type", "ILUT"); _se(eig, "max_iteration_step", sol["linear_iter"])
    _se(eig, "error_tolerance", sol["linear_tol"]); _se(eig, "scaling", "true")

    # Curves — Leistungsprofil für PowerCurveConstantFlow (nur Modus B)
    if power_curve is not None:
        t, v = power_curve[1]
        curves = _se(root, "curves")
        c = _se(curves, "curve")
        _se(c, "name", "power_curve")
        _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
        _se(c, "values", " ".join(f"{x:.6e}" for x in v))

    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


# ======================================================================
#  Main
# ======================================================================
def run_ogs(prj_path: Path) -> int:
    ogs = shutil.which("ogs") or shutil.which("ogs.exe")
    if not ogs:
        print("ogs.exe nicht im PATH — Setup ist fertig, aber kein Lauf.",
              file=sys.stderr)
        return 1
    return subprocess.call([ogs, str(prj_path), "-o", str(prj_path.parent)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-mesh", action="store_true")
    ap.add_argument("--no-run",  action="store_true")
    args = ap.parse_args()

    out_dir = Path(CONFIG["output"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = CONFIG["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"
    n_bhe = len(_bhe_positions(CONFIG))

    if not args.no_mesh:
        print(f"[1/3] gmsh: 3D-Mesh mit {n_bhe} eingebetteten BHE-Linien ...")
        build_mesh(CONFIG, out_dir)
        print(f"      -> {msh_path}")
        print("[2/3] msh2vtu (gemischt-dimensional: 3D + 1D) ...")
        # OGS HEAT_TRANSPORT_BHE benötigt 3D-Soil + 1D-BHE-Elemente
        # im selben Domain-Mesh mit unterschiedlichen MaterialIDs.
        msh2vtu(filename=msh_path, output_path=out_dir,
                output_prefix=prefix, dim=[3, 1], reindex=True)

    mesh_files = {
        "domain":  f"{prefix}_domain.vtu",
        "top":     f"{prefix}_physical_group_top.vtu",
        "bottom":  f"{prefix}_physical_group_bottom.vtu",
        "lateral": f"{prefix}_physical_group_lateral.vtu",
    }
    for i in range(n_bhe):
        mesh_files[f"bhe_{i:02d}"] = f"{prefix}_physical_group_bhe_{i:02d}.vtu"

    print(f"[3/3] OGS-Projektdatei (HEAT_TRANSPORT_BHE, {n_bhe} BHEs) ...")
    prj_path = build_prj(CONFIG, out_dir, mesh_files)
    print(f"      -> {prj_path}")

    if args.no_run:
        return 0
    return run_ogs(prj_path)


if __name__ == "__main__":
    sys.exit(main())

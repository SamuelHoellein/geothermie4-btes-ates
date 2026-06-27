#!/usr/bin/env python3
# coding: utf-8
"""
Saisonaler ATES-Solarspeicher  —  OGS HEAT_CONDUCTION + Neumann-Linienquellen

Eine Doublette (Warmbrunnen + Kaltbrunnen, Abstand 120 m) im Aquifer (84–122 m).
Grundwasserfluss = 0. Warmbrunnen und Kaltbrunnen haben stets entgegengesetzte Leistung.

Ablauf:
  1. Energiebilanz aus CSVs → monatliche Leistung je Doublette [W]
  2. Mesh (gmsh): geschichtetes 3D-Gebiet + 2 × 1D-Brunnenlinien im Aquifer
  3. .msh → .vtu (ogstools)
  4. OGS-Projektdatei (.prj) schreiben
  5. OGS starten → Ergebnisse in out/<szenario>/ als *.pvd

Verwendung:
  python ates_seasonal.py
  python ates_seasonal.py --scenario after_renovation
  python ates_seasonal.py --no-mesh --no-run
"""

import argparse
import csv as _csv
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from xml.etree import ElementTree as ET

import gmsh
import numpy as np


# Nutzereinstellungen
COP                  = 4.0
T_GROUND_INITIAL_C   = 10.0
DEFAULT_SCENARIO     = "before_renovation"

WELL_SPACING_M        = 120.0   # Abstand Warm- zu Kaltbrunnen [m]
N_DOUBLETS            = 1       # Anzahl Doubletten (Energiebilanz)
AQUIFER_TOP_M         = 84.0    # Aquifer-Oberkante [m u. GOK]
AQUIFER_BOTTOM_M      = 122.0   # Aquifer-Unterkante [m u. GOK]
DOMAIN_BUFFER_M       = 100.0   # seitlicher Puffer über Brunnenfeld hinaus [m]
DOMAIN_DEPTH_BUFFER_M = 25.0    # Puffer unterhalb Aquifer [m]

POROSITY_DEFAULT = 0.02
FLUID_LAMBDA     = 0.6          # Wärmeleitfähigkeit Porenwasser [W/mK]
RAMP_DAYS        = 3

N_YEARS    = 3
DT_SECONDS = 30 * 86400.0

MESH_SIZE_NEAR_M   = 3.0
MESH_SIZE_FAR_M    = 20.0
MESH_RADIUS_NEAR_M = 15.0
MESH_RADIUS_FAR_M  = 60.0

LINEAR_TOL     = 1.0e-12
LINEAR_ITER    = 10000
NONLINEAR_ITER = 2

_ROOT    = Path(__file__).resolve().parent.parent.parent
DATA_DIR = _ROOT / "Data Input"

DAY         = 86400.0
MONTH_DAYS  = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MONTH_NAMES = ["Jan","Feb","Mär","Apr","Mai","Jun","Jul","Aug","Sep","Okt","Nov","Dez"]


# Eingangsdaten & Energiebilanz

def _read_csv(path):
    with open(path, encoding="utf-8-sig") as f:
        rows = list(_csv.DictReader(f, delimiter=";"))
    return [{k: (v or "").replace(",", ".") for k, v in r.items()} for r in rows]


def load_solar():
    rows = _read_csv(DATA_DIR / "Solarthermie.csv")
    col  = list(rows[0].keys())[1]
    return np.array([float(r[col]) for r in rows])


def load_demand(scenario):
    rows = _read_csv(DATA_DIR / "Heizwaermebedarfe.csv")
    keys = list(rows[0].keys())
    col  = keys[1] if scenario == "before_renovation" else keys[2]
    return np.array([float(r[col]) for r in rows])


def compute_monthly_powers(scenario):
    """
    Monatliche Nettoleistung [W] je Doublette (Warmbrunnen-Seite).
      WP-Entzug = Heizlast × (1 − 1/COP)
      E_net = E_solar − WP-Entzug
      P [W] = E_net [kWh] × 1000 [W/kW] / (Tage × 24 h) / N_DOUBLETS
    P > 0 = Laden (Wärme in Warmbrunnen), P < 0 = Entladen.
    Kaltbrunnen hat stets −P.
    """
    solar  = load_solar()
    hp_ext = load_demand(scenario) * (1.0 - 1.0 / COP)
    net    = solar - hp_ext
    power_W = [net[i] * 1000.0 / (MONTH_DAYS[i] * 24.0) / N_DOUBLETS
               for i in range(12)]
    return power_W, {"solar": solar, "hp_ext": hp_ext, "net": net}


def print_energy_balance(power_W, info, scenario):
    label = "vor Sanierung" if scenario == "before_renovation" else "nach Sanierung"
    s, hp, n = info["solar"], info["hp_ext"], info["net"]
    print(f"\n{'='*64}")
    print(f"  ATES Energiebilanz [{label}]")
    print(f"{'='*64}")
    print(f"  {'Monat':<5} {'Solar':>8} {'WP-Entzug':>11} {'Netto':>9} {'P/Doublette':>13}")
    print(f"  {'':5} {'[MWh]':>8} {'[MWh]':>11} {'[MWh]':>9} {'[W]':>13}")
    print("  " + "-"*49)
    for i in range(12):
        print(f"  {MONTH_NAMES[i]:<5} {s[i]/1e3:>8.1f} {hp[i]/1e3:>11.1f}"
              f" {n[i]/1e3:>9.1f} {power_W[i]:>13.0f}")
    print(f"  {'Summe':<5} {s.sum()/1e3:>8.1f} {hp.sum()/1e3:>11.1f} {n.sum()/1e3:>9.1f}")
    print(f"\n  Solardeckung WP-Entzug: {s.sum()/hp.sum()*100:.1f}%"
          f"   Max. Laden: +{max(power_W):.0f} W"
          f"   Max. Entladen: {min(power_W):.0f} W")
    print(f"{'='*64}\n")


# Geologische Schichten

def build_layers_from_csv(domain_depth_m):
    rows = _read_csv(DATA_DIR / "Ground_1.csv")
    keys = list(rows[0].keys())

    def _f(row, col, default):
        val = row[keys[col]].strip().lstrip("<>")
        try:
            return float(val)
        except ValueError:
            return default

    depths = [_f(r, 0, None) for r in rows]
    layers = []
    for i, z_top in enumerate(depths):
        if z_top is None or z_top >= domain_depth_m:
            break
        next_depth = (depths[i+1] if i+1 < len(depths) and depths[i+1] is not None
                      else domain_depth_m)
        z_bot = min(next_depth, domain_depth_m)
        h = z_bot - z_top
        if h <= 0:
            continue
        lam   = _f(rows[i], 1, 2.0)
        cv    = _f(rows[i], 3, 1.8) * 1e6
        rho_s = _f(rows[i], 4, 2.65) * 1000
        phi   = _f(rows[i], 6, 10.0) / 100
        layers.append({
            "name"       : f"layer_{i:02d}",
            "thickness_m": float(h),
            "rho_s"      : rho_s,
            "cp_s"       : cv / rho_s,
            "phi"        : phi,
            "lambda_eff" : phi * FLUID_LAMBDA + (1.0 - phi) * lam,
        })
    return layers


# Config

def build_config(scenario):
    power_W, info = compute_monthly_powers(scenario)
    print_energy_balance(power_W, info, scenario)
    layers = build_layers_from_csv(AQUIFER_BOTTOM_M + DOMAIN_DEPTH_BUFFER_M)
    return {
        "prefix" : f"ates_{scenario}",
        "domain" : {
            "Lx"    : WELL_SPACING_M + 2 * DOMAIN_BUFFER_M,
            "Ly"    : 2 * DOMAIN_BUFFER_M,
            "z_base": 0.0,
        },
        "layers" : layers,
        "aquifer": {
            "top"     : AQUIFER_TOP_M,
            "bottom"  : AQUIFER_BOTTOM_M,
            "L_active": AQUIFER_BOTTOM_M - AQUIFER_TOP_M,
        },
        "doublet": {"spacing": WELL_SPACING_M},
        "mesh"   : {"near": MESH_SIZE_NEAR_M, "far": MESH_SIZE_FAR_M,
                    "r_near": MESH_RADIUS_NEAR_M, "r_far": MESH_RADIUS_FAR_M},
        "T0_K"   : T_GROUND_INITIAL_C + 273.15,
        "cycles" : {"n": N_YEARS, "power_W": power_W, "ramp_days": RAMP_DAYS},
        "time"   : {"dt": DT_SECONDS},
        "solver" : {"tol": LINEAR_TOL, "iter_lin": LINEAR_ITER, "iter_nl": NONLINEAR_ITER},
    }


# Mesh

def msh2vtu(msh_path, out_dir, prefix):
    import ogstools as ot
    meshes = ot.Meshes.from_gmsh(filename=str(msh_path), dim=[3, 1], reindex=True, log=False)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        fname = (f"{prefix}_domain.vtu" if name == "domain"
                 else f"{prefix}_physical_group_{name}.vtu")
        mesh.save(str(out_dir / fname), binary=True)


def _layer_stack(cfg):
    z, out = cfg["domain"]["z_base"], []
    for L in reversed(cfg["layers"]):
        h = float(L["thickness_m"])
        out.append({**L, "z_low": z, "z_high": z + h})
        z += h
    return out, z


def build_mesh(cfg, out_dir):
    prefix = cfg["prefix"]
    msh_out = out_dir / f"{prefix}.msh"
    Lx, Ly  = cfg["domain"]["Lx"], cfg["domain"]["Ly"]
    z_base  = cfg["domain"]["z_base"]
    x0, y0  = -Lx / 2.0, -Ly / 2.0
    layers, z_top = _layer_stack(cfg)
    aq  = cfg["aquifer"]
    sp  = cfg["doublet"]["spacing"]
    z_well_top = z_top - aq["top"]
    z_well_bot = z_top - aq["bottom"]
    m = cfg["mesh"]

    well_pos = [(-sp / 2.0, 0.0), (sp / 2.0, 0.0)]  # warm, cold

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ates")

    boxes = [gmsh.model.occ.addBox(x0, y0, L["z_low"], Lx, Ly, L["z_high"] - L["z_low"])
             for L in layers]
    lines = []
    for x, y in well_pos:
        p1 = gmsh.model.occ.addPoint(x, y, z_well_top)
        p2 = gmsh.model.occ.addPoint(x, y, z_well_bot)
        lines.append(gmsh.model.occ.addLine(p1, p2))

    gmsh.model.occ.fragment([(3, b) for b in boxes], [(1, l) for l in lines])
    gmsh.model.occ.synchronize()

    vol_layer = {i: [] for i in range(len(layers))}
    for _, tag in gmsh.model.getEntities(3):
        bb = gmsh.model.occ.getBoundingBox(3, tag)
        zc = 0.5 * (bb[2] + bb[5])
        for i, L in enumerate(layers):
            if L["z_low"] - 1e-6 <= zc <= L["z_high"] + 1e-6:
                vol_layer[i].append(tag); break

    pos = np.array(well_pos)
    seg_well = {i: [] for i in range(2)}
    all_segs = []
    for _, tag in gmsh.model.getEntities(1):
        bb = gmsh.model.occ.getBoundingBox(1, tag)
        dx, dy, dz = bb[3]-bb[0], bb[4]-bb[1], bb[5]-bb[2]
        xc, yc, zc = 0.5*(bb[0]+bb[3]), 0.5*(bb[1]+bb[4]), 0.5*(bb[2]+bb[5])
        if dx < 1e-6 and dy < 1e-6 and dz > 1e-9 and z_well_bot-1e-3 <= zc <= z_well_top+1e-3:
            j = int(np.argmin(np.hypot(pos[:, 0]-xc, pos[:, 1]-yc)))
            if np.hypot(pos[j, 0]-xc, pos[j, 1]-yc) < 1e-3:
                seg_well[j].append(tag); all_segs.append(tag)

    if any(not s for s in seg_well.values()):
        raise RuntimeError("Mindestens ein Brunnen ohne Liniensegment.")

    pg = 1
    for i, L in enumerate(layers):
        gmsh.model.addPhysicalGroup(3, vol_layer[i], tag=pg, name=L["name"]); pg += 1
    gmsh.model.addPhysicalGroup(1, seg_well[0], tag=pg, name="well_warm"); pg += 1
    gmsh.model.addPhysicalGroup(1, seg_well[1], tag=pg, name="well_cold"); pg += 1

    top_f, bot_f = [], []
    for _, tag in gmsh.model.getEntities(2):
        bb = gmsh.model.occ.getBoundingBox(2, tag)
        zc, w = 0.5*(bb[2]+bb[5]), bb[3]-bb[0]
        if w >= 0.9 * Lx:
            if   abs(zc - z_top)  < 1e-6: top_f.append(tag)
            elif abs(zc - z_base) < 1e-6: bot_f.append(tag)
    gmsh.model.addPhysicalGroup(2, top_f, tag=200, name="top")
    gmsh.model.addPhysicalGroup(2, bot_f, tag=201, name="bottom")

    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", all_segs)
    ft = gmsh.model.mesh.field.add("Threshold")
    for k, v in [("InField", fd), ("SizeMin", m["near"]), ("SizeMax", m["far"]),
                 ("DistMin", m["r_near"]), ("DistMax", m["r_far"])]:
        gmsh.model.mesh.field.setNumber(ft, k, v)
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)

    well_pts = [(0, t) for tag in all_segs
                for dim, t in gmsh.model.getBoundary([(1, tag)], oriented=False) if dim == 0]
    if well_pts:
        gmsh.model.mesh.setSize(list(set(well_pts)), 0.5)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_out))
    gmsh.finalize()
    return msh_out


# Leistungskurve

def build_power_curve(cfg):
    """q_warm(t) [W/m]; Kaltbrunnen erhält −q über neg_unity-Parameter in OGS."""
    monthly   = cfg["cycles"]["power_W"]
    L         = cfg["aquifer"]["L_active"]
    ramp      = max(60.0, cfg["cycles"]["ramp_days"] * DAY)
    month_dur = 365.25 / 12.0 * DAY
    times, vals = [0.0], [0.0]
    t = 0.0
    for _ in range(cfg["cycles"]["n"]):
        for P in monthly:
            q = P / L
            t += ramp;                 times.append(t); vals.append(q)
            hold = max(0.0, month_dur - ramp)
            if hold > 0.0:
                t += hold;             times.append(t); vals.append(q)
    t += ramp; times.append(t); vals.append(0.0)
    return t, (np.array(times), np.array(vals))


# PRJ-Hilfsfunktionen

def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None:
        el.text = str(text)
    return el

def _prop(parent, name, value):
    p = _se(parent, "property")
    _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)

def _indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip(): elem.text = i + "  "
        for child in elem: _indent(child, level + 1)
        if not child.tail or not child.tail.strip(): child.tail = i
    if level and (not elem.tail or not elem.tail.strip()): elem.tail = i


# OGS-Projektdatei

def build_prj(cfg, out_dir, mesh_files):
    prefix  = cfg["prefix"]
    layers  = list(reversed(cfg["layers"]))  # bottom → top = MaterialID 0, 1, 2 …
    sol     = cfg["solver"]
    T0      = cfg["T0_K"]
    t_end, (t_arr, v_arr) = build_power_curve(cfg)
    n_steps = int(t_end // cfg["time"]["dt"]) + 1

    root = ET.Element("OpenGeoSysProject")

    meshes = _se(root, "meshes")
    for key in ("domain", "top", "bottom", "well_warm", "well_cold"):
        _se(meshes, "mesh", mesh_files[key])

    proc = _se(_se(root, "processes"), "process")
    _se(proc, "name", "HeatConduction")
    _se(proc, "type", "HEAT_CONDUCTION")
    _se(proc, "integration_order", 2)
    _se(_se(proc, "process_variables"), "process_variable", "temperature")

    def _add_medium(media_el, mid, soil):
        med   = _se(media_el, "medium", id=mid)
        _se(med, "phases")
        props = _se(med, "properties")
        _prop(props, "density",                soil["rho_s"])
        _prop(props, "specific_heat_capacity", soil["cp_s"])
        _prop(props, "thermal_conductivity",   soil["lambda_eff"])
        _prop(props, "porosity",               soil.get("phi", POROSITY_DEFAULT))

    media = _se(root, "media")
    for mid, soil in enumerate(layers):
        _add_medium(media, mid, soil)
    mid_soil = layers[len(layers) // 2]
    for mid in range(len(layers), len(layers) + 2):  # +2 für Warm- und Kaltbrunnen
        _add_medium(media, mid, mid_soil)

    tl  = _se(root, "time_loop")
    ptl = _se(_se(tl, "processes"), "process", ref="HeatConduction")
    _se(ptl, "nonlinear_solver", "basic_picard")
    conv = _se(ptl, "convergence_criterion")
    _se(conv, "type", "DeltaX"); _se(conv, "norm_type", "NORM2"); _se(conv, "reltol", 1e-6)
    _se(_se(ptl, "time_discretization"), "type", "BackwardEuler")
    ts = _se(ptl, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0); _se(ts, "t_end", t_end)
    pair = _se(_se(ts, "timesteps"), "pair")
    _se(pair, "repeat", n_steps); _se(pair, "delta_t", cfg["time"]["dt"])
    out_tl = _se(tl, "output")
    _se(out_tl, "type", "VTK"); _se(out_tl, "prefix", prefix)
    _se(_se(out_tl, "variables"), "variable", "temperature")
    pair2 = _se(_se(out_tl, "timesteps"), "pair")
    _se(pair2, "repeat", n_steps); _se(pair2, "each_steps", 1)

    params = _se(root, "parameters")
    def _param(name, **kw):
        p = _se(params, "parameter")
        _se(p, "name", name)
        for k, v in kw.items(): _se(p, k, v)
    _param("T0",         type="Constant",    value=T0)
    _param("unity",      type="Constant",    value="1")
    _param("neg_unity",  type="Constant",    value="-1")
    _param("warm_power", type="CurveScaled", curve="power_curve", parameter="unity")
    _param("cold_power", type="CurveScaled", curve="power_curve", parameter="neg_unity")

    pvars = _se(root, "process_variables")
    pvs   = _se(pvars, "process_variable")
    _se(pvs, "name", "temperature")
    _se(pvs, "components", 1); _se(pvs, "order", 1)
    _se(pvs, "initial_condition", "T0")
    bcs = _se(pvs, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh", Path(mesh_files[face]).stem)
        _se(bc, "type", "Dirichlet"); _se(bc, "parameter", "T0")
    sts = _se(pvs, "source_terms")
    for well, param in [("well_warm", "warm_power"), ("well_cold", "cold_power")]:
        st = _se(sts, "source_term")
        _se(st, "mesh",      Path(mesh_files[well]).stem)
        _se(st, "type",      "Volumetric")
        _se(st, "parameter", param)

    nl = _se(_se(root, "nonlinear_solvers"), "nonlinear_solver")
    _se(nl, "name", "basic_picard"); _se(nl, "type", "Picard")
    _se(nl, "max_iter", sol["iter_nl"])
    _se(nl, "linear_solver", "general_linear_solver")
    ls = _se(_se(root, "linear_solvers"), "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    _se(eig, "solver_type",        "BiCGSTAB")
    _se(eig, "precon_type",        "ILUT")
    _se(eig, "max_iteration_step", sol["iter_lin"])
    _se(eig, "error_tolerance",    sol["tol"])
    _se(eig, "scaling",            "true")

    c = _se(_se(root, "curves"), "curve")
    _se(c, "name",   "power_curve")
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t_arr))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v_arr))

    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


# OGS starten

def _find_ogs():
    ogs = shutil.which("ogs") or shutil.which("ogs.exe")
    if ogs:
        return [ogs]
    py = Path(sys.executable)
    for c in [py.parent/"ogs.exe", py.parent/"ogs",
              py.parent/"Scripts"/"ogs.exe", py.parent/"Scripts"/"ogs"]:
        if c.is_file():
            return [str(c)]
    try:
        r = subprocess.run([sys.executable, "-m", "ogs", "--version"],
                           capture_output=True, timeout=10)
        if r.returncode == 0:
            return [sys.executable, "-m", "ogs"]
    except Exception:
        pass
    return None


def run_ogs(prj_path, n_steps=None):
    cmd = _find_ogs()
    if cmd is None:
        print(f"FEHLER: ogs nicht gefunden.\n"
              f"  → conda activate ghe23\n"
              f"  → ogs {prj_path} -o {prj_path.parent}", file=sys.stderr)
        return 1

    step_re = re.compile(r"Time step #(\d+) started")
    width   = len(str(n_steps)) if n_steps else 4
    t_step1 = None

    def _fmt_eta(seconds):
        if seconds < 60:   return f"{seconds:.0f} s"
        if seconds < 3600: return f"{seconds/60:.0f} min"
        return f"{seconds/3600:.1f} h"

    proc = subprocess.Popen(
        cmd + [str(prj_path), "-o", str(prj_path.parent)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )
    for line in proc.stdout:
        m = step_re.search(line)
        if m:
            step = int(m.group(1))
            now  = time.time()
            if step == 1:
                t_step1 = now
            eta_str = ""
            if n_steps and t_step1 and step > 1:
                avg_per_step = (now - t_step1) / (step - 1)
                eta_s        = avg_per_step * (n_steps - step)
                eta_str      = f"  noch ~{_fmt_eta(eta_s)}"
            total = f"/{n_steps}" if n_steps else ""
            pct   = f" ({step/n_steps*100:.0f}%)" if n_steps else ""
            bar   = "█" * int(step/n_steps*20) + "░" * (20 - int(step/n_steps*20)) if n_steps else ""
            print(f"\r  {bar} Schritt {step:{width}}{total}{pct}{eta_str}   ",
                  end="", flush=True)
        elif any(kw in line for kw in ("warning:", "error:", "critical:")):
            print(f"\n{line.rstrip()}", flush=True)
    proc.wait()
    print()
    return proc.returncode


# Plot

def plot_results(cfg, out_dir):
    import xml.etree.ElementTree as ET
    import matplotlib.pyplot as plt
    import pyvista as pv

    pvd_path = out_dir / f"{cfg['prefix']}.pvd"
    if not pvd_path.exists():
        print(f"Kein PVD gefunden: {pvd_path}", file=sys.stderr)
        return

    entries = [
        (float(ds.get("timestep")), out_dir / ds.get("file"))
        for ds in ET.parse(pvd_path).getroot().findall(".//DataSet")
        if float(ds.get("timestep")) > 0
    ]
    if not entries:
        return

    first = pv.read(str(entries[0][1]))
    z_top = float(first.points[:, 2].max())
    aq    = cfg["aquifer"]
    sp    = cfg["doublet"]["spacing"]
    z_mid = z_top - (aq["top"] + aq["bottom"]) / 2.0

    query_points = {
        "Warmbrunnen": np.array([-sp / 2, 0.0, z_mid]),
        "Mitte (x=0)": np.array([0.0,     0.0, z_mid]),
        "Kaltbrunnen": np.array([ sp / 2, 0.0, z_mid]),
    }

    times_yr = np.array([t / (365.25 * 86400) for t, _ in entries])
    results  = {label: [] for label in query_points}

    print("Erstelle Plot ...")
    for _, vtu in entries:
        mesh = pv.read(str(vtu))
        pts  = mesh.points
        for label, query in query_points.items():
            idx = int(np.argmin(np.linalg.norm(pts - query, axis=1)))
            results[label].append(float(mesh.point_data["temperature"][idx]) - 273.15)

    fig, ax = plt.subplots(figsize=(11, 5))
    for (label, temps), color in zip(results.items(), ["tab:red", "tab:green", "tab:blue"]):
        ax.plot(times_yr, temps, label=label, color=color, linewidth=1.5)
    ax.axhline(T_GROUND_INITIAL_C, color="gray", linestyle="--", linewidth=0.8,
               label=f"Ausgangstemperatur ({T_GROUND_INITIAL_C} °C)")
    ax.axhline(0, color="red", linestyle=":", linewidth=0.8, label="0 °C")
    ax.set_xlabel("Zeit [Jahre]")
    ax.set_ylabel("Temperatur [°C]")
    ax.set_title(
        f"ATES – Aquifer-Temperatur (Tiefe {(aq['top']+aq['bottom'])/2:.0f} m) · {cfg['prefix']}")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    plt.tight_layout()

    out_png = out_dir / f"temperature_ates_{cfg['prefix']}.png"
    plt.savefig(out_png, dpi=150)
    print(f"Plot gespeichert: {out_png}")
    plt.show()


# Main

def main():
    ap = argparse.ArgumentParser(description="ATES Saisonalspeicher")
    ap.add_argument("--scenario", default=DEFAULT_SCENARIO,
                    choices=["before_renovation", "after_renovation"])
    ap.add_argument("--no-mesh", action="store_true")
    ap.add_argument("--no-run",  action="store_true")
    args = ap.parse_args()

    out_dir = Path("out") / args.scenario
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg    = build_config(args.scenario)
    prefix = cfg["prefix"]

    if not args.no_mesh:
        print(f"[1/3] Mesh: Doublette ±{WELL_SPACING_M/2:.0f} m, "
              f"Aquifer {AQUIFER_TOP_M:.0f}–{AQUIFER_BOTTOM_M:.0f} m ...")
        msh_path = build_mesh(cfg, out_dir)
        print(f"[2/3] msh → vtu ...")
        msh2vtu(msh_path, out_dir, prefix)

    mesh_files = {
        "domain":    f"{prefix}_domain.vtu",
        "top":       f"{prefix}_physical_group_top.vtu",
        "bottom":    f"{prefix}_physical_group_bottom.vtu",
        "well_warm": f"{prefix}_physical_group_well_warm.vtu",
        "well_cold": f"{prefix}_physical_group_well_cold.vtu",
    }

    print(f"[3/3] PRJ (1 Doublette, {N_YEARS} Jahre) ...")
    prj_path = build_prj(cfg, out_dir, mesh_files)
    print(f"      → {prj_path}")

    if args.no_run:
        return 0

    t_end, _ = build_power_curve(cfg)
    n_steps  = int(t_end // cfg["time"]["dt"]) + 1
    print(f"\nStarte OGS ... {n_steps} Zeitschritte, Ergebnisse in {out_dir}/")
    rc = run_ogs(prj_path, n_steps=n_steps)
    if rc == 0:
        plot_results(cfg, out_dir)
    return rc


if __name__ == "__main__":
    sys.exit(main())

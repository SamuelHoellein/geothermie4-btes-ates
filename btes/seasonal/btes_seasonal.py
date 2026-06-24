#!/usr/bin/env python3
# coding: utf-8
"""
Saisonaler BTES-Solarspeicher — OpenGeoSys HEAT_TRANSPORT_BHE

Liest Eingangsdaten aus Data Input/:
  - Ground.csv            : Schichttiefen [m], lambda [W/mK], cv [J/m³K]
  - Solarthermie.csv      : Monatlicher Solarertrag [kWh] (Gesamtanlage)
  - Heizwaermebedarfe.csv : Monatlicher Heizwärme- + WHW-Bedarf [kWh],
                            Spalte 1 = vor Sanierung, Spalte 2 = nach Sanierung

Energiebilanz je Monat:
  E_net = E_solar  -  (1 - 1/COP) * E_heizung
  E_net > 0  -> Laden   (Solar  -> BTES)
  E_net < 0  -> Entladen (BTES  -> Waermepumpe)

Verwendung:
  python btes_seasonal.py [--scenario before_renovation|after_renovation]
                          [--no-mesh] [--no-run]
"""
from __future__ import annotations

import argparse
import csv as _csv
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import gmsh
import numpy as np

# ═══════════════════════════════════════════════════════════════════════════════
#  NUTZEREINSTELLUNGEN  — alle anpassbaren Größen hier
# ═══════════════════════════════════════════════════════════════════════════════

# ── Systemparameter ────────────────────────────────────────────────────────────
COP                  = 4.0    # Jahresarbeitszahl der Wärmepumpe [-]
T_INJECTION_C        = 65.0   # Einspeis-Temperatur Solar -> BTES [°C]
T_EXTRACT_MIN_C      =  5.0   # Mindest-Entnahmetemperatur BTES -> WP [°C]
T_GROUND_INITIAL_C   = 10.0   # Anfangstemperatur Untergrund [°C]

# ── Szenario-Voreinstellung (per --scenario überschreibbar) ────────────────────
DEFAULT_SCENARIO = "before_renovation"   # oder "after_renovation"

# ── Sondenfeld — Auslegung ─────────────────────────────────────────────────────
N_BOREHOLES_TOTAL       = 100    # Gesamtzahl Sonden im Realsystem
N_BHE_X                 =   5   # Sonden in x-Richtung (simuliertes Feld)
N_BHE_Y                 =   5   # Sonden in y-Richtung (simuliertes Feld)
BOREHOLE_SPACING_M      =  8.0  # Sondenabstand [m]
BOREHOLE_DEPTH_TOP_M    =  5.0  # Sondenkopf unter Gelände [m]
BOREHOLE_DEPTH_BOTTOM_M = 105.0 # Sondenfuß  unter Gelände [m]
DOMAIN_BUFFER_M         = 25.0  # seitlicher Puffer außerhalb des Felds [m]
DOMAIN_DEPTH_BUFFER_M   = 25.0  # Tiefe unterhalb des Sondenfußes [m]

# ── Untergrund ─────────────────────────────────────────────────────────────────
# cv aus Ground.csv = bulk-Wärmekapazität (Feststoff + Porenwasser).
# Aufspaltung: cp_s = (cv - phi * rho_f * cp_f) / ((1-phi) * rho_s)
RHO_SOLID_KG_M3         = 2650.0  # Feststoffdichte aller Schichten [kg/m³]
POROSITY_DEFAULT        = 0.02    # Porosität (Festgestein, BTES)
PERMEABILITY_DEFAULT_M2 = 1.0e-18 # Permeabilität [m²]

# ── Simulation ─────────────────────────────────────────────────────────────────
N_YEARS      = 3           # Simulationsjahre (= Ladezyklen)
DT_SECONDS   = 7 * 86400.0 # Zeitschrittweite [s]  (1 Woche)
OUTPUT_EVERY = 1           # Ausgabe alle N Zeitschritte

# ── BHE-Hardware (1U-Sonde) ────────────────────────────────────────────────────
BHE_BOREHOLE_DIAMETER_M    = 0.15
BHE_PIPE_OUTER_DIAMETER_M  = 0.032
BHE_PIPE_WALL_THICKNESS_M  = 0.003
BHE_PIPE_WALL_LAMBDA       = 0.4     # W/mK
BHE_PIPE_DISTANCE_M        = 0.05   # Achsabstand Vor-/Rücklauf [m]
BHE_GROUT_DENSITY          = 2190.0  # kg/m³
BHE_GROUT_POROSITY         = 0.0
BHE_GROUT_CP               = 1735.0  # J/kgK
BHE_GROUT_LAMBDA           = 2.3     # W/mK
BHE_REFR_DENSITY           = 1052.0  # kg/m³  (Wasser-Glykol-Mischung)
BHE_REFR_VISCOSITY         = 0.0052  # Pa·s
BHE_REFR_CP                = 3795.0  # J/kgK
BHE_REFR_LAMBDA            = 0.48    # W/mK
BHE_REFR_T_REF_K           = 293.15  # K
BHE_FLOW_RATE_KG_S         = 0.2     # Massenstrom je Sonde [kg/s]

# ── Mesh ───────────────────────────────────────────────────────────────────────
MESH_SIZE_NEAR_M   = 1.5
MESH_SIZE_FAR_M    = 12.0
MESH_RADIUS_NEAR_M = 8.0
MESH_RADIUS_FAR_M  = 30.0

# ── Solver ─────────────────────────────────────────────────────────────────────
LINEAR_TOL     = 1.0e-12
LINEAR_ITER    = 10000
NONLINEAR_ITER = 20

# ── Fluid (Grundwasser im Untergrund) ──────────────────────────────────────────
FLUID_RHO    = 1000.0
FLUID_VISC   = 1.0e-3
FLUID_CP     = 4180.0
FLUID_LAMBDA = 0.6

# ═══════════════════════════════════════════════════════════════════════════════
#  PFADE  (relativ zur Projektstruktur)
# ═══════════════════════════════════════════════════════════════════════════════
_ROOT    = Path(__file__).resolve().parent.parent.parent
DATA_DIR = _ROOT / "Data Input"

DAY        = 86400.0
MONTH_DAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MONTH_NAMES = ["Jan","Feb","Mär","Apr","Mai","Jun",
               "Jul","Aug","Sep","Okt","Nov","Dez"]


# ═══════════════════════════════════════════════════════════════════════════════
#  DATENLADEN & ENERGIEBILANZ
# ═══════════════════════════════════════════════════════════════════════════════

def _read_csv(path: Path) -> list[dict]:
    """Liest semikolon-getrennte CSV mit deutschem Komma-Dezimalzeichen."""
    with open(path, encoding="utf-8-sig") as f:
        rows = list(_csv.DictReader(f, delimiter=";"))
    return [{k: v.replace(",", ".") for k, v in r.items()} for r in rows]


def load_solar() -> np.ndarray:
    """12 Monatswerte Solarertrag [kWh] aus Solarthermie.csv."""
    rows = _read_csv(DATA_DIR / "Solarthermie.csv")
    col  = list(rows[0].keys())[1]
    return np.array([float(r[col]) for r in rows])


def load_demand(scenario: str) -> np.ndarray:
    """12 Monatswerte Heizwärmebedarf [kWh] aus Heizwaermebedarfe.csv."""
    rows = _read_csv(DATA_DIR / "Heizwaermebedarfe.csv")
    keys = list(rows[0].keys())
    col  = keys[1] if scenario == "before_renovation" else keys[2]
    return np.array([float(r[col]) for r in rows])


def compute_monthly_powers(scenario: str) -> tuple[list[float], dict]:
    """
    Berechnet die monatliche Nettoleistung [W] je Sonde.

    Bilanz:
      E_net[m] = E_solar[m]  -  (1 - 1/COP) * E_demand[m]
      P_bhe[m] = E_net[m] * 1000 / (days[m] * 86400) / N_BOREHOLES_TOTAL

    Returns:
        power_W : Liste von 12 Werten [W/Sonde], + = laden, - = entladen
        info    : dict mit Energiebilanz-Größen [kWh/a]
    """
    solar   = load_solar()
    demand  = load_demand(scenario)

    hp_frac = 1.0 - 1.0 / COP          # 0.75 bei COP = 4
    hp_ext  = demand * hp_frac          # WP-Entzug aus BTES [kWh]
    net_kwh = solar - hp_ext            # Netto je Monat [kWh]

    power_W = []
    for net, days in zip(net_kwh, MONTH_DAYS):
        p_total = net * 1000.0 / (days * DAY)     # Gesamtleistung System [W]
        power_W.append(p_total / N_BOREHOLES_TOTAL)

    info = {
        "solar_kwh_a"   : float(solar.sum()),
        "demand_kwh_a"  : float(demand.sum()),
        "hp_ext_kwh_a"  : float(hp_ext.sum()),
        "net_kwh_a"     : float(net_kwh.sum()),
        "solar_monthly" : solar.tolist(),
        "demand_monthly": demand.tolist(),
        "net_monthly"   : net_kwh.tolist(),
    }
    return power_W, info


def print_energy_balance(power_W: list[float], info: dict, scenario: str) -> None:
    label = "vor Sanierung" if scenario == "before_renovation" else "nach Sanierung"
    print()
    print("=" * 66)
    print(f"  BTES Saisonalspeicher — Energiebilanz  [{label}]")
    print("=" * 66)
    print(f"  {'Monat':<6} {'Solar':>9} {'WP-Entzug':>11} {'Netto':>10} {'P/Sonde':>10}")
    print(f"  {'':6} {'[MWh]':>9} {'[MWh]':>11} {'[MWh]':>10} {'[W]':>10}")
    print("  " + "-" * 50)
    solar  = info["solar_monthly"]
    demand = info["demand_monthly"]
    net    = info["net_monthly"]
    for i in range(12):
        hp = demand[i] * (1.0 - 1.0 / COP)
        print(f"  {MONTH_NAMES[i]:<6} {solar[i]/1e3:>9.1f} {hp/1e3:>11.1f} "
              f"{net[i]/1e3:>10.1f} {power_W[i]:>10.0f}")
    print("  " + "-" * 50)
    print(f"  {'Gesamt':<6} {info['solar_kwh_a']/1e3:>9.1f} "
          f"{info['hp_ext_kwh_a']/1e3:>11.1f} {info['net_kwh_a']/1e3:>10.1f}")
    print()
    cover = info["solar_kwh_a"] / info["demand_kwh_a"] * 100
    solar_vs_hp = info["solar_kwh_a"] / info["hp_ext_kwh_a"] * 100
    max_charge   = max(power_W)
    max_discharge = min(power_W)
    print(f"  Solardeckung Heizlast    : {cover:.1f} %")
    print(f"  Solardeckung WP-Entzug   : {solar_vs_hp:.1f} %")
    print(f"  Restdefizit (Erdwärme)   : {-info['net_kwh_a']/1e3:.0f} MWh/a")
    print(f"  Max. Ladeleistung/Sonde  : +{max_charge:+.0f} W")
    print(f"  Max. Entladeleistung/Sonde: {max_discharge:.0f} W")
    print(f"  Simuliertes Feld         : {N_BHE_X}x{N_BHE_Y} = {N_BHE_X*N_BHE_Y} Sonden "
          f"(von {N_BOREHOLES_TOTAL} gesamt)")
    print()
    if info["net_kwh_a"] < 0:
        deficit = -info["net_kwh_a"]
        print(f"  HINWEIS: Jährliches BTES-Defizit = {deficit/1e3:.0f} MWh/a.")
        print(f"  Das Restdefizit wird durch natürliche Erdwärme (T_Grund) gedeckt.")
        print(f"  Bitte Mindesttemperatur T_EXTRACT_MIN = {T_EXTRACT_MIN_C}°C in den")
        print(f"  Ergebnisplots prüfen.")
    print("=" * 66)
    print()


# ═══════════════════════════════════════════════════════════════════════════════
#  SCHICHTEN AUS GROUND.CSV
# ═══════════════════════════════════════════════════════════════════════════════

def build_layers_from_csv(domain_depth_m: float) -> list[dict]:
    """
    Liest Ground.csv und gibt CONFIG-kompatible Schichten (oben -> unten) zurück.
    Schichten werden auf domain_depth_m begrenzt; unterhalb der letzten CSV-Zeile
    wird die letzte bekannte Schicht fortgeschrieben.

    Jeder Layer-Dict hat:
      name, thickness_m, permeability_m2, porosity, rho_s_kg_m3,
      cp_s_J_kgK, lambda_s_W_mK
    """
    rows   = _read_csv(DATA_DIR / "Ground.csv")
    keys   = list(rows[0].keys())
    depths  = [float(r[keys[0]]) for r in rows]
    lambdas = [float(r[keys[1]]) for r in rows]
    cvs     = [float(r[keys[2]]) for r in rows]

    phi   = POROSITY_DEFAULT
    rho_s = RHO_SOLID_KG_M3

    layers = []
    for i in range(len(depths)):
        z_top = depths[i]
        if z_top >= domain_depth_m:
            break
        z_bot = depths[i + 1] if i + 1 < len(depths) else domain_depth_m
        z_bot = min(z_bot, domain_depth_m)
        thickness = z_bot - z_top
        if thickness <= 0.0:
            continue

        lam = lambdas[i]
        cv  = cvs[i]
        # bulk cv = (1-phi)*rho_s*cp_s + phi*rho_f*cp_f  ->  solve for cp_s
        cp_s = (cv - phi * FLUID_RHO * FLUID_CP) / ((1.0 - phi) * rho_s)

        layers.append({
            "name"            : f"layer_{i:02d}",
            "thickness_m"     : float(thickness),
            "permeability_m2" : PERMEABILITY_DEFAULT_M2,
            "porosity"        : phi,
            "rho_s_kg_m3"     : rho_s,
            "cp_s_J_kgK"      : float(cp_s),
            "lambda_s_W_mK"   : float(lam),
        })
    return layers


# ═══════════════════════════════════════════════════════════════════════════════
#  CONFIG ZUSAMMENBAUEN
# ═══════════════════════════════════════════════════════════════════════════════

def build_config(scenario: str, out_dir: Path) -> dict:
    """Erstellt den vollständigen CONFIG-Dict aus Nutzereinstellungen + CSV-Daten."""
    power_W, info = compute_monthly_powers(scenario)
    print_energy_balance(power_W, info, scenario)

    domain_depth = BOREHOLE_DEPTH_BOTTOM_M + DOMAIN_DEPTH_BUFFER_M
    layers = build_layers_from_csv(domain_depth)

    # Domänengröße aus Feld + Puffer berechnen
    field_x = (N_BHE_X - 1) * BOREHOLE_SPACING_M
    field_y = (N_BHE_Y - 1) * BOREHOLE_SPACING_M
    size_x  = field_x + 2.0 * DOMAIN_BUFFER_M
    size_y  = field_y + 2.0 * DOMAIN_BUFFER_M

    prefix = f"btes_seasonal_{scenario}"

    return {
        "domain": {
            "size_x_m": size_x,
            "size_y_m": size_y,
            "z_base_m": 0.0,
        },
        "layers": layers,
        "borehole": {
            "depth_top_m"   : BOREHOLE_DEPTH_TOP_M,
            "depth_bottom_m": BOREHOLE_DEPTH_BOTTOM_M,
        },
        "field": {
            "n_x"      : N_BHE_X,
            "n_y"      : N_BHE_Y,
            "spacing_m": BOREHOLE_SPACING_M,
            "positions": None,
        },
        "mesh": {
            "size_near_field_m"     : MESH_SIZE_NEAR_M,
            "size_far_m"            : MESH_SIZE_FAR_M,
            "field_size_radius_m"   : MESH_RADIUS_NEAR_M,
            "field_size_radius_far_m": MESH_RADIUS_FAR_M,
        },
        "fluid": {
            "rho_ref_kg_m3" : FLUID_RHO,
            "T_ref_K"       : T_GROUND_INITIAL_C + 273.15,
            "viscosity_Pa_s": FLUID_VISC,
            "cp_J_kgK"      : FLUID_CP,
            "lambda_W_mK"   : FLUID_LAMBDA,
        },
        "bhe": {
            "type"    : "1U",
            "borehole": {"diameter_m": BHE_BOREHOLE_DIAMETER_M},
            "pipes": {
                "diameter_outer_m"              : BHE_PIPE_OUTER_DIAMETER_M,
                "wall_thickness_m"              : BHE_PIPE_WALL_THICKNESS_M,
                "wall_thermal_conductivity_W_mK": BHE_PIPE_WALL_LAMBDA,
                "distance_between_pipes_m"      : BHE_PIPE_DISTANCE_M,
                "longitudinal_dispersion_length_m": 0.001,
            },
            "grout": {
                "density_kg_m3"            : BHE_GROUT_DENSITY,
                "porosity"                 : BHE_GROUT_POROSITY,
                "specific_heat_capacity_J_kgK": BHE_GROUT_CP,
                "thermal_conductivity_W_mK": BHE_GROUT_LAMBDA,
            },
            "refrigerant": {
                "density_kg_m3"            : BHE_REFR_DENSITY,
                "viscosity_Pa_s"           : BHE_REFR_VISCOSITY,
                "specific_heat_capacity_J_kgK": BHE_REFR_CP,
                "thermal_conductivity_W_mK": BHE_REFR_LAMBDA,
                "reference_temperature_K"  : BHE_REFR_T_REF_K,
            },
            "control": {
                "type"         : "PowerCurveConstantFlow",
                "flow_rate_kg_s": BHE_FLOW_RATE_KG_S,
            },
        },
        "initial": {
            "T_K"                       : T_GROUND_INITIAL_C + 273.15,
            "p_Pa"                      : 0.0,
            "T_surface_K"               : T_GROUND_INITIAL_C + 273.15,
            "geothermal_gradient_K_per_m": 0.0,
        },
        "cycles": {
            "n_cycles"       : N_YEARS,
            "monthly_power_W": power_W,
            "ramp_days"      : 3.0,
        },
        "time": {
            "dt_seconds"         : DT_SECONDS,
            "output_every_n_steps": OUTPUT_EVERY,
        },
        "output": {
            "prefix"   : prefix,
            "out_dir"  : str(out_dir),
            "variables": ["temperature_soil"],
        },
        "solver": {
            "linear_tol"    : LINEAR_TOL,
            "linear_iter"   : LINEAR_ITER,
            "nonlinear_iter": NONLINEAR_ITER,
        },
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  MESH  (gmsh)  — identisch zur Logik in btes_3d_bhe.py
# ═══════════════════════════════════════════════════════════════════════════════

def msh2vtu(filename, output_path, output_prefix, dim, reindex=True):
    import ogstools as ot
    meshes = ot.Meshes.from_gmsh(filename=str(filename), dim=dim,
                                  reindex=reindex, log=False)
    out = Path(output_path)
    out.mkdir(parents=True, exist_ok=True)
    for name, mesh in meshes.items():
        fname = (f"{output_prefix}_domain.vtu" if name == "domain"
                 else f"{output_prefix}_physical_group_{name}.vtu")
        mesh.save(str(out / fname), binary=True)


def _bhe_positions(cfg: dict) -> list[tuple[float, float]]:
    fld = cfg["field"]
    if fld.get("positions"):
        return [tuple(p) for p in fld["positions"]]
    nx, ny, sp = fld["n_x"], fld["n_y"], fld["spacing_m"]
    xs = (np.arange(nx) - (nx - 1) / 2.0) * sp
    ys = (np.arange(ny) - (ny - 1) / 2.0) * sp
    return [(float(x), float(y)) for y in ys for x in xs]


def _z_for_depth(cfg: dict, depth: float) -> float:
    z_base = cfg["domain"]["z_base_m"]
    z_top  = z_base + sum(float(L["thickness_m"]) for L in cfg["layers"])
    return z_top - depth


def _layer_stack(cfg: dict):
    z_base = cfg["domain"].get("z_base_m", 0.0)
    bot_up = list(reversed(cfg["layers"]))
    z = z_base
    out = []
    for L in bot_up:
        h = float(L["thickness_m"])
        out.append({**L, "z_low": z, "z_high": z + h})
        z += h
    return out, z


def build_mesh(cfg: dict, out_dir: Path) -> Path:
    prefix   = cfg["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"
    Lx       = cfg["domain"]["size_x_m"]
    Ly       = cfg["domain"]["size_y_m"]
    z_base   = cfg["domain"]["z_base_m"]
    layers, z_top = _layer_stack(cfg)

    bhe_pos   = _bhe_positions(cfg)
    z_bhe_top = z_top - cfg["borehole"]["depth_top_m"]
    z_bhe_bot = z_top - cfg["borehole"]["depth_bottom_m"]
    m         = cfg["mesh"]
    x0, y0    = -Lx / 2.0, -Ly / 2.0

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("btes_seasonal")

    layer_boxes = [
        gmsh.model.occ.addBox(x0, y0, L["z_low"], Lx, Ly, L["z_high"] - L["z_low"])
        for L in layers
    ]
    lines = []
    for x, y in bhe_pos:
        p_top = gmsh.model.occ.addPoint(x, y, z_bhe_top)
        p_bot = gmsh.model.occ.addPoint(x, y, z_bhe_bot)
        lines.append(gmsh.model.occ.addLine(p_top, p_bot))

    gmsh.model.occ.fragment([(3, b) for b in layer_boxes],
                             [(1, l) for l in lines])
    gmsh.model.occ.synchronize()

    # Schichtvolumen nach z-Mittelpunkt zuordnen
    vol_layer = {i: [] for i in range(len(layers))}
    for dim, tag in gmsh.model.getEntities(3):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (bb[2] + bb[5])
        for i, L in enumerate(layers):
            if L["z_low"] - 1e-6 <= zc <= L["z_high"] + 1e-6:
                vol_layer[i].append(tag)
                break

    # BHE-Liniensegmente nach Position gruppieren
    pos = np.array(bhe_pos)
    seg_bhe = {i: [] for i in range(len(bhe_pos))}
    all_segs = []
    for dim, tag in gmsh.model.getEntities(1):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        dx, dy, dz = bb[3]-bb[0], bb[4]-bb[1], bb[5]-bb[2]
        xc, yc, zc = 0.5*(bb[0]+bb[3]), 0.5*(bb[1]+bb[4]), 0.5*(bb[2]+bb[5])
        if dx < 1e-6 and dy < 1e-6 and dz > 1e-9:
            if (z_bhe_bot - 1e-3) <= zc <= (z_bhe_top + 1e-3):
                d = np.hypot(pos[:, 0] - xc, pos[:, 1] - yc)
                j = int(np.argmin(d))
                if d[j] < 1e-3:
                    seg_bhe[j].append(tag)
                    all_segs.append(tag)
    missing = [i for i, s in seg_bhe.items() if not s]
    if missing:
        raise RuntimeError(f"BHE ohne Liniensegment: {missing}")

    # Physical Groups: Schichten -> MaterialID 0..L-1, BHEs -> L..L+N-1
    pg = 1
    for i, L in enumerate(layers):
        gmsh.model.addPhysicalGroup(3, vol_layer[i], tag=pg, name=L["name"])
        pg += 1
    for i in range(len(bhe_pos)):
        gmsh.model.addPhysicalGroup(1, seg_bhe[i], tag=pg, name=f"bhe_{i:02d}")
        pg += 1

    # Domänenränder
    top_faces, bot_faces, lat_faces = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (bb[2] + bb[5])
        if (bb[3]-bb[0]) >= 0.9*Lx and abs(zc - z_top) < 1e-6:
            top_faces.append(tag)
        elif (bb[3]-bb[0]) >= 0.9*Lx and abs(zc - z_base) < 1e-6:
            bot_faces.append(tag)
        else:
            on_outer = (
                (abs(bb[0]-x0)       < 1e-6 and abs(bb[3]-x0)       < 1e-6) or
                (abs(bb[0]-(x0+Lx))  < 1e-6 and abs(bb[3]-(x0+Lx))  < 1e-6) or
                (abs(bb[1]-y0)       < 1e-6 and abs(bb[4]-y0)       < 1e-6) or
                (abs(bb[1]-(y0+Ly))  < 1e-6 and abs(bb[4]-(y0+Ly))  < 1e-6)
            )
            if on_outer:
                lat_faces.append(tag)

    gmsh.model.addPhysicalGroup(2, top_faces, tag=200, name="top")
    gmsh.model.addPhysicalGroup(2, bot_faces, tag=201, name="bottom")
    if lat_faces:
        gmsh.model.addPhysicalGroup(2, lat_faces, tag=202, name="lateral")

    # Verfeinerung um BHE-Linien
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
        gmsh.model.mesh.setSize(list(set(bhe_pts)), 0.4)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints",        1)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion",            2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(str(msh_path))
    gmsh.finalize()
    return msh_path


# ═══════════════════════════════════════════════════════════════════════════════
#  LEISTUNGSKURVE  (PowerCurveConstantFlow)
# ═══════════════════════════════════════════════════════════════════════════════

def build_power_curve(cfg: dict):
    monthly = cfg["cycles"].get("monthly_power_W")
    if monthly is None:
        return None
    assert len(monthly) == 12
    n         = cfg["cycles"]["n_cycles"]
    ramp      = max(60.0, cfg["cycles"]["ramp_days"] * DAY)
    month_dur = 365.25 / 12.0 * DAY
    times = [0.0]; vals = [0.0]; t_now = 0.0
    for _ in range(n):
        for P in monthly:
            t_now += ramp;  times.append(t_now); vals.append(float(P))
            hold   = max(0.0, month_dur - ramp)
            if hold > 0.0:
                t_now += hold; times.append(t_now); vals.append(float(P))
    t_now += ramp; times.append(t_now); vals.append(0.0)
    return t_now, (np.array(times), np.array(vals))


# ═══════════════════════════════════════════════════════════════════════════════
#  PRJ  (HEAT_TRANSPORT_BHE)
# ═══════════════════════════════════════════════════════════════════════════════

def _se(parent, tag, text=None, **attrs):
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None:
        el.text = str(text)
    return el

def _const_prop(parent, name, value):
    p = _se(parent, "property")
    _se(p, "name", name); _se(p, "type", "Constant"); _se(p, "value", value)

def _indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        for child in elem:
            _indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = i
    if level and (not elem.tail or not elem.tail.strip()):
        elem.tail = i


def _bhe_xml(parent, cfg: dict, n_bhe: int) -> None:
    bhe   = cfg["bhe"]
    L_bhe = cfg["borehole"]["depth_bottom_m"] - cfg["borehole"]["depth_top_m"]
    bhes  = _se(parent, "borehole_heat_exchangers")
    for i in range(n_bhe):
        b = _se(bhes, "borehole_heat_exchanger")
        _se(b, "type", bhe["type"])
        ftc = _se(b, "flow_and_temperature_control")
        flow_vol = bhe["control"]["flow_rate_kg_s"] / bhe["refrigerant"]["density_kg_m3"]
        _se(ftc, "type",        "PowerCurveConstantFlow")
        _se(ftc, "power_curve", "power_curve")
        _se(ftc, "flow_rate",   flow_vol)
        bh = _se(b, "borehole")
        _se(bh, "length",   L_bhe)
        _se(bh, "diameter", bhe["borehole"]["diameter_m"])
        pipes = _se(b, "pipes")
        for pname in ("inlet", "outlet"):
            p = _se(pipes, pname)
            _se(p, "diameter",                  bhe["pipes"]["diameter_outer_m"])
            _se(p, "wall_thickness",            bhe["pipes"]["wall_thickness_m"])
            _se(p, "wall_thermal_conductivity", bhe["pipes"]["wall_thermal_conductivity_W_mK"])
        _se(pipes, "distance_between_pipes",         bhe["pipes"]["distance_between_pipes_m"])
        _se(pipes, "longitudinal_dispersion_length", bhe["pipes"]["longitudinal_dispersion_length_m"])
        gr = _se(b, "grout")
        _se(gr, "density",                bhe["grout"]["density_kg_m3"])
        _se(gr, "porosity",               bhe["grout"]["porosity"])
        _se(gr, "specific_heat_capacity", bhe["grout"]["specific_heat_capacity_J_kgK"])
        _se(gr, "thermal_conductivity",   bhe["grout"]["thermal_conductivity_W_mK"])
        rf = _se(b, "refrigerant")
        _se(rf, "density",                bhe["refrigerant"]["density_kg_m3"])
        _se(rf, "viscosity",              bhe["refrigerant"]["viscosity_Pa_s"])
        _se(rf, "specific_heat_capacity", bhe["refrigerant"]["specific_heat_capacity_J_kgK"])
        _se(rf, "thermal_conductivity",   bhe["refrigerant"]["thermal_conductivity_W_mK"])
        _se(rf, "reference_temperature",  bhe["refrigerant"]["reference_temperature_K"])


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict) -> Path:
    prefix    = cfg["output"]["prefix"]
    fluid     = cfg["fluid"]
    init      = cfg["initial"]
    sol       = cfg["solver"]
    n_bhe     = len(_bhe_positions(cfg))

    power_curve = build_power_curve(cfg)
    t_end   = power_curve[0]
    n_steps = int(t_end // cfg["time"]["dt_seconds"]) + 1

    root   = ET.Element("OpenGeoSysProject")
    meshes = _se(root, "meshes")
    for key in ("domain", "top", "bottom", "lateral"):
        _se(meshes, "mesh", mesh_files[key])
    for i in range(n_bhe):
        _se(meshes, "mesh", mesh_files[f"bhe_{i:02d}"])

    # Process
    processes = _se(root, "processes")
    proc      = _se(processes, "process")
    _se(proc, "name", "HeatTransportBHE")
    _se(proc, "type", "HEAT_TRANSPORT_BHE")
    _se(proc, "integration_order", 2)
    pv = _se(proc, "process_variables")
    _se(pv, "process_variable", "temperature_soil")
    for i in range(n_bhe):
        _se(pv, "process_variable", f"temperature_BHE{i+1}")
    _bhe_xml(proc, cfg, n_bhe)

    # Media
    media       = _se(root, "media")
    layers_bot  = list(reversed(cfg["layers"]))
    bh_mat      = layers_bot[len(layers_bot) // 2]
    mats        = layers_bot + [bh_mat] * n_bhe
    for mid, soil in enumerate(mats):
        med = _se(media, "medium", id=mid)
        phs = _se(med, "phases")
        ph  = _se(phs, "phase"); _se(ph, "type", "AqueousLiquid")
        pp  = _se(ph, "properties")
        _const_prop(pp, "density",                fluid["rho_ref_kg_m3"])
        _const_prop(pp, "viscosity",              fluid["viscosity_Pa_s"])
        _const_prop(pp, "specific_heat_capacity", fluid["cp_J_kgK"])
        _const_prop(pp, "thermal_conductivity",   fluid["lambda_W_mK"])
        _const_prop(pp, "phase_velocity",         "0 0 0")
        ph  = _se(phs, "phase"); _se(ph, "type", "Solid")
        pp  = _se(ph, "properties")
        _const_prop(pp, "density",                soil["rho_s_kg_m3"])
        _const_prop(pp, "specific_heat_capacity", soil["cp_s_J_kgK"])
        _const_prop(pp, "thermal_conductivity",   soil["lambda_s_W_mK"])
        props = _se(med, "properties")
        _const_prop(props, "porosity",     soil["porosity"])
        _const_prop(props, "permeability", soil["permeability_m2"])
        lam_eff = (soil["porosity"] * fluid["lambda_W_mK"] +
                   (1.0 - soil["porosity"]) * soil["lambda_s_W_mK"])
        _const_prop(props, "thermal_conductivity", lam_eff)
        _const_prop(props, "storage", 0.0)

    # Time loop
    tl    = _se(root, "time_loop")
    procs = _se(tl, "processes")
    pref  = _se(procs, "process", ref="HeatTransportBHE")
    _se(pref, "nonlinear_solver", "basic_picard")
    conv  = _se(pref, "convergence_criterion")
    _se(conv, "type", "DeltaX"); _se(conv, "norm_type", "NORM2")
    _se(conv, "reltol", sol.get("rel_tol_T", 1e-4))
    td    = _se(pref, "time_discretization"); _se(td, "type", "BackwardEuler")
    ts    = _se(pref, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0); _se(ts, "t_end", t_end)
    pair  = _se(_se(ts, "timesteps"), "pair")
    _se(pair, "repeat", n_steps); _se(pair, "delta_t", cfg["time"]["dt_seconds"])
    out   = _se(tl, "output")
    _se(out, "type", "VTK"); _se(out, "prefix", prefix)
    out_v = _se(out, "variables")
    _se(out_v, "variable", "temperature_soil")
    for i in range(n_bhe):
        _se(out_v, "variable", f"temperature_BHE{i+1}")
    pair2 = _se(_se(out, "timesteps"), "pair")
    _se(pair2, "repeat", n_steps)
    _se(pair2, "each_steps", cfg["time"]["output_every_n_steps"])

    # Parameters & IC
    params = _se(root, "parameters")
    p_T0   = _se(params, "parameter")
    _se(p_T0, "name", "T0"); _se(p_T0, "type", "Constant")
    _se(p_T0, "value", init["T_K"])
    p_T0b  = _se(params, "parameter")
    _se(p_T0b, "name", "T0_BHE"); _se(p_T0b, "type", "Constant")
    T0 = init["T_K"]
    _se(p_T0b, "value", f"{T0} {T0} {T0} {T0}")

    # Process variables
    pvars = _se(root, "process_variables")
    pvs   = _se(pvars, "process_variable")
    _se(pvs, "name", "temperature_soil")
    _se(pvs, "components", 1); _se(pvs, "order", 1)
    _se(pvs, "initial_condition", "T0")
    bcs   = _se(pvs, "boundary_conditions")
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files[face]).stem)
        _se(bc, "type",      "Dirichlet")
        _se(bc, "parameter", "T0")
    for i in range(n_bhe):
        pvb = _se(pvars, "process_variable")
        _se(pvb, "name", f"temperature_BHE{i+1}")
        _se(pvb, "components", 4)
        _se(pvb, "order", 1)
        _se(pvb, "initial_condition", "T0_BHE")
        _se(pvb, "boundary_conditions")

    # Solvers
    nls = _se(root, "nonlinear_solvers")
    n   = _se(nls, "nonlinear_solver")
    _se(n, "name", "basic_picard"); _se(n, "type", "Picard")
    _se(n, "max_iter", sol["nonlinear_iter"])
    _se(n, "linear_solver", "general_linear_solver")
    lss = _se(root, "linear_solvers")
    ls  = _se(lss, "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    _se(eig, "solver_type",        "BiCGSTAB")
    _se(eig, "precon_type",        "ILUT")
    _se(eig, "max_iteration_step", sol["linear_iter"])
    _se(eig, "error_tolerance",    sol["linear_tol"])
    _se(eig, "scaling",            "true")

    # Leistungskurve (PowerCurveConstantFlow)
    t, v   = power_curve[1]
    curves = _se(root, "curves")
    c      = _se(curves, "curve")
    _se(c, "name",   "power_curve")
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t))
    _se(c, "values", " ".join(f"{x:.6e}" for x in v))

    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


# ═══════════════════════════════════════════════════════════════════════════════
#  OGS STARTEN
# ═══════════════════════════════════════════════════════════════════════════════

def run_ogs(prj_path: Path) -> int:
    ogs_exe = shutil.which("ogs") or shutil.which("ogs.exe")
    if ogs_exe:
        cmd = [ogs_exe, str(prj_path), "-o", str(prj_path.parent)]
    else:
        # Fallback: ogs als Python-Modul aufrufen (funktioniert auch ohne
        # aktiviertes venv/conda, solange ogs im selben Python installiert ist)
        try:
            subprocess.check_call([sys.executable, "-m", "ogs", "--version"],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            cmd = [sys.executable, "-m", "ogs", str(prj_path),
                   "-o", str(prj_path.parent)]
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("FEHLER: 'ogs' weder im PATH noch als Python-Modul gefunden.\n"
                  "  -> conda activate ghe23  (oder die passende Umgebung aktivieren)\n"
                  "  -> dann das Skript erneut starten.",
                  file=sys.stderr)
            return 1
    return subprocess.call(cmd)


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="Saisonaler BTES-Solarspeicher (HEAT_TRANSPORT_BHE)")
    ap.add_argument("--scenario", default=DEFAULT_SCENARIO,
                    choices=["before_renovation", "after_renovation"],
                    help="Gebäudeszenario (Standard: %(default)s)")
    ap.add_argument("--no-mesh", action="store_true",
                    help="Mesh-Generierung überspringen (VTU muss bereits vorhanden sein)")
    ap.add_argument("--no-run",  action="store_true",
                    help="OGS-Lauf überspringen (nur Mesh + PRJ erzeugen)")
    args = ap.parse_args()

    out_dir = Path("out") / args.scenario
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[BTES Saisonalspeicher]  Szenario: {args.scenario}")

    cfg    = build_config(args.scenario, out_dir)
    prefix = cfg["output"]["prefix"]
    n_bhe  = len(_bhe_positions(cfg))

    msh_path = out_dir / f"{prefix}.msh"

    if not args.no_mesh:
        print(f"[1/3] gmsh: {N_BHE_X}x{N_BHE_Y}-Feld, {n_bhe} BHEs, "
              f"Domain {cfg['domain']['size_x_m']:.0f}x{cfg['domain']['size_y_m']:.0f} m "
              f"x {BOREHOLE_DEPTH_BOTTOM_M + DOMAIN_DEPTH_BUFFER_M:.0f} m tief ...")
        build_mesh(cfg, out_dir)
        print(f"      -> {msh_path}")
        print("[2/3] msh2vtu ...")
        msh2vtu(msh_path, out_dir, prefix, dim=[3, 1], reindex=True)

    mesh_files = {
        "domain" : f"{prefix}_domain.vtu",
        "top"    : f"{prefix}_physical_group_top.vtu",
        "bottom" : f"{prefix}_physical_group_bottom.vtu",
        "lateral": f"{prefix}_physical_group_lateral.vtu",
    }
    for i in range(n_bhe):
        mesh_files[f"bhe_{i:02d}"] = f"{prefix}_physical_group_bhe_{i:02d}.vtu"

    print(f"[3/3] OGS-Projektdatei ({n_bhe} BHEs, {N_YEARS} Jahre) ...")
    prj_path = build_prj(cfg, out_dir, mesh_files)
    print(f"      -> {prj_path}")

    if args.no_run:
        print("--no-run gesetzt: OGS nicht gestartet.")
        return 0

    print(f"\nStarte OGS ({N_YEARS * 12} Monate, dt = {DT_SECONDS/86400:.0f} Tage) ...")
    return run_ogs(prj_path)


if __name__ == "__main__":
    sys.exit(main())

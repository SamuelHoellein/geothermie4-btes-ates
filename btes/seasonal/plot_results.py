#!/usr/bin/env python3
"""
Auswertung & Plots für die BTES-3D-Übung.

Erzeugt in figures/:
  - 1_well_temperature.png       T(t) am zentralen Bohrloch (Mitteltiefe)
  - 2_field_snapshots.png        T-Feld in der xy-Mittel­ebene zu 4 Zeitpunkten
  - 3_energy_balance.png         Gespeicherte Wärmemenge & Recovery-Effizienz
  - 4_plume_extent.png           Radius der thermischen Front (ΔT > 1 K)
"""
from __future__ import annotations

import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

DELTA_T_THRESHOLD = 1.0
N_SNAPSHOTS       = 4
DAY = 86400.0


def load_config():
    """Erkennt automatisch, welche BTES-3D-Variante gelaufen ist, und lädt
    deren CONFIG. Bevorzugt die BHE-Variante (Standard-Übung mit echten
    U-Rohr-Sonden); fällt auf die vereinfachte Volumenquelle zurück."""
    here = Path(__file__).parent
    for prefix, module in (("btes_3d_bhe", "btes_3d_bhe"),
                            ("btes_3d",     "btes_3d")):
        if (here / "out" / f"{prefix}.pvd").exists():
            cfg = __import__(module).CONFIG
            return cfg
    # nichts gelaufen → vereinfachte Variante als Default importieren
    return __import__("btes_3d").CONFIG


def temperature_field(mesh) -> str:
    """Name des 3D-Temperaturfelds: `temperature_soil` (BHE) oder `T`
    (vereinfachte Variante)."""
    for name in ("temperature_soil", "T"):
        if name in mesh.point_data:
            return name
    raise KeyError("Kein Temperaturfeld (temperature_soil/T) im Mesh gefunden.")


def read_pvd(pvd_path: Path):
    tree = ET.parse(pvd_path); items = []
    base = pvd_path.parent
    for ds in tree.iter("DataSet"):
        items.append((float(ds.attrib["timestep"]), base / ds.attrib["file"]))
    items.sort(); return items


def effective_rho_cp(cfg):
    layers = cfg["layers"]
    soil = layers[len(layers) // 2]   # repräsentative (mittlere) Schicht = Sondenmaterial
    phi  = soil["porosity"]
    return (phi * cfg["fluid"]["rho_ref_kg_m3"] * cfg["fluid"]["cp_J_kgK"]
            + (1 - phi) * soil["rho_s_kg_m3"] * soil["cp_s_J_kgK"])


def borehole_positions(cfg):
    fld = cfg["field"]
    if fld.get("positions"):
        return np.array(fld["positions"])
    nx, ny, sp = fld["n_x"], fld["n_y"], fld["spacing_m"]
    xs = (np.arange(nx) - (nx - 1) / 2) * sp
    ys = (np.arange(ny) - (ny - 1) / 2) * sp
    return np.array([(x, y) for y in ys for x in xs])


def cartesian_integral(mesh, field):
    vols = mesh.compute_cell_sizes(length=False, area=False, volume=True)["Volume"]
    if field.shape[0] == mesh.n_points:
        cell_field = mesh.point_data_to_cell_data().cell_data[mesh.point_data.keys()[0]] \
            if False else None
        cell_field = np.empty(mesh.n_cells)
        for i in range(mesh.n_cells):
            pts = mesh.get_cell(i).point_ids
            cell_field[i] = field[pts].mean()
    else:
        cell_field = field
    return float(np.sum(cell_field * vols))


def main() -> int:
    CONFIG = load_config()
    out_dir = Path(CONFIG["output"]["out_dir"])
    prefix  = CONFIG["output"]["prefix"]
    pvd     = out_dir / f"{prefix}.pvd"
    if not pvd.exists():
        print(f"FEHLER: {pvd} nicht gefunden. Bitte zuerst die Simulation ausführen.")
        return 1

    fig_dir = Path(__file__).parent / "figures"
    fig_dir.mkdir(exist_ok=True)

    T0 = CONFIG["initial"]["T_K"]
    rho_cp = effective_rho_cp(CONFIG)
    positions = borehole_positions(CONFIG)

    # Domäne in btes_3d.py: z=0 unten → Oberfläche bei z_total.
    # Sondentiefe wird von der Oberfläche nach unten gemessen.
    z_total = sum(L["thickness_m"] for L in CONFIG["layers"])
    z_field_top = z_total - CONFIG["borehole"]["depth_top_m"]
    z_field_bot = z_total - CONFIG["borehole"]["depth_bottom_m"]
    z_mid = 0.5 * (z_field_top + z_field_bot)

    steps = read_pvd(pvd)
    times = np.array([t for t, _ in steps]); times_d = times / DAY

    TEMP = temperature_field(pv.read(steps[0][1]))   # "temperature_soil" (BHE) oder "T"

    # 1) T(t) zentrales BH ---------------------------------------------
    bh_center = positions[len(positions) // 2]
    probe_point = np.array([[bh_center[0], bh_center[1], z_mid]])
    probe_T = []
    for _, f in steps:
        m = pv.read(f); probe_T.append(float(pv.PolyData(probe_point).sample(m)[TEMP][0]))
    probe_T = np.array(probe_T)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, probe_T - 273.15, lw=2)
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("T zentrale Sonde, Mitteltiefe [°C]")
    ax.set_title("Temperatur im Sondenfeld")
    ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "1_well_temperature.png", dpi=130); plt.close(fig)

    # 2) Snapshots (xy-Slice in Mitteltiefe) ---------------------------
    idxs = np.linspace(0, len(steps) - 1, N_SNAPSHOTS).astype(int)
    fig, axes = plt.subplots(1, N_SNAPSHOTS, figsize=(4 * N_SNAPSHOTS, 4.5),
                              sharey=True)
    for ax, i in zip(axes, idxs):
        t, f = steps[i]
        m = pv.read(f)
        sl = m.slice(normal="z", origin=(0, 0, z_mid))
        pts = sl.points
        sc = ax.tricontourf(pts[:, 0], pts[:, 1], sl[TEMP] - 273.15,
                            levels=20, cmap="inferno")
        ax.scatter(positions[:, 0], positions[:, 1], c="cyan", s=8, marker="x")
        ax.set_xlabel("x [m]"); ax.set_title(f"t = {t/DAY:.0f} d"); ax.set_aspect("equal")
    axes[0].set_ylabel("y [m]")
    fig.colorbar(sc, ax=axes, shrink=0.7, label="T [°C]")
    fig.savefig(fig_dir / "2_field_snapshots.png", dpi=130); plt.close(fig)

    # 3) Energie -------------------------------------------------------
    energy = []
    for _, f in steps:
        m = pv.read(f); dT = np.asarray(m[TEMP]) - T0
        energy.append(rho_cp * cartesian_integral(m, dT))
    energy = np.array(energy) / 1e9

    E_max = energy.max(); E_end = energy[-1]
    recovery = (E_max - E_end) / E_max if E_max > 0 else float("nan")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, energy, lw=2)
    ax.axhline(E_max, ls="--", color="r", alpha=0.5, label=f"max. gespeichert: {E_max:.2f} GJ")
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel("ΔE [GJ]")
    ax.set_title(f"Energiebilanz – Recovery ≈ {recovery*100:.0f} %")
    ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "3_energy_balance.png", dpi=130); plt.close(fig)

    # 4) Plume-Reichweite (Radius vom Feld-Schwerpunkt) ---------------
    cx, cy = positions.mean(axis=0)
    r_front = []
    for _, f in steps:
        m = pv.read(f); pts = m.points; dT = np.asarray(m[TEMP]) - T0
        in_zone = (pts[:, 2] >= z_field_bot) & (pts[:, 2] <= z_field_top)
        mask = in_zone & (dT > DELTA_T_THRESHOLD)
        if mask.any():
            d = np.hypot(pts[mask, 0] - cx, pts[mask, 1] - cy)
            r_front.append(float(d.max()))
        else:
            r_front.append(0.0)
    r_front = np.array(r_front)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_d, r_front, lw=2)
    ax.set_xlabel("Zeit [d]"); ax.set_ylabel(f"r(ΔT > {DELTA_T_THRESHOLD} K) [m]")
    ax.set_title("Reichweite der thermischen Front (horizontaler Radius vom Feldzentrum)")
    ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(fig_dir / "4_plume_extent.png", dpi=130); plt.close(fig)

    print(f"Plots → {fig_dir}")
    print(f"  Recovery: {recovery*100:.1f} %  |  Emax: {E_max:.2f} GJ  |  max. Plume: {r_front.max():.1f} m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

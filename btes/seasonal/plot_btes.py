#!/usr/bin/env python3
# coding: utf-8
"""
Plottet die Temperatur an der mittigen Sonde über alle Zeitschritte.
Verwendung: python plot_btes.py [--scenario before_renovation|after_renovation]
"""

import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# Muss mit btes_seasonal.py übereinstimmen
N_BHE_X              = 2
N_BHE_Y              = 2
BOREHOLE_SPACING_M   = 8.0
BOREHOLE_DEPTH_TOP_M    = 5.0
BOREHOLE_DEPTH_BOTTOM_M = 159.0
T_GROUND_INITIAL_C   = 10.0


def bhe_positions():
    xs = (np.arange(N_BHE_X) - (N_BHE_X - 1) / 2.0) * BOREHOLE_SPACING_M
    ys = (np.arange(N_BHE_Y) - (N_BHE_Y - 1) / 2.0) * BOREHOLE_SPACING_M
    return [(float(x), float(y)) for y in ys for x in xs]


def nearest_temperature(mesh, query):
    """Temperatur [K] am nächsten Gitterpunkt zu query [x,y,z]."""
    idx = int(np.argmin(np.linalg.norm(mesh.points - query, axis=1)))
    return float(mesh.point_data["temperature"][idx])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", default="before_renovation",
                    choices=["before_renovation", "after_renovation"])
    args = ap.parse_args()

    out_dir  = Path("out") / args.scenario
    pvd_path = out_dir / f"btes_{args.scenario}.pvd"
    if not pvd_path.exists():
        raise FileNotFoundError(f"PVD nicht gefunden: {pvd_path}")

    # Zeitschritte aus PVD einlesen
    tree = ET.parse(pvd_path)
    entries = [
        (float(ds.get("timestep")), out_dir / ds.get("file"))
        for ds in tree.getroot().findall(".//DataSet")
        if float(ds.get("timestep")) > 0   # Initialschritt überspringen
    ]
    if not entries:
        raise RuntimeError("Keine Zeitschritte im PVD gefunden.")

    # z_top aus erstem VTU bestimmen (Surface-Niveau)
    first_mesh = pv.read(str(entries[0][1]))
    z_top = float(first_mesh.points[:, 2].max())

    # Mittlere Sonde (nächste zu Ursprung)
    positions = bhe_positions()
    cx, cy = min(positions, key=lambda p: p[0]**2 + p[1]**2)

    # Drei Tiefen entlang der Sonde
    depths = {
        "Oben (−5 m)":   z_top - BOREHOLE_DEPTH_TOP_M,
        "Mitte (−82 m)":  z_top - (BOREHOLE_DEPTH_TOP_M + BOREHOLE_DEPTH_BOTTOM_M) / 2.0,
        "Unten (−159 m)": z_top - BOREHOLE_DEPTH_BOTTOM_M,
    }

    # Temperatur je Tiefe über alle Zeitschritte extrahieren
    times_yr = np.array([t / (365.25 * 86400) for t, _ in entries])
    results  = {label: [] for label in depths}

    print(f"Lese {len(entries)} Zeitschritte ...")
    for _, vtu in entries:
        mesh = pv.read(str(vtu))
        for label, z in depths.items():
            T_K = nearest_temperature(mesh, np.array([cx, cy, z]))
            results[label].append(T_K - 273.15)

    # Plot
    fig, ax = plt.subplots(figsize=(11, 5))
    colors = ["tab:blue", "tab:orange", "tab:green"]
    for (label, temps), color in zip(results.items(), colors):
        ax.plot(times_yr, temps, label=label, color=color, linewidth=1.5)

    ax.axhline(T_GROUND_INITIAL_C, color="gray", linestyle="--",
               linewidth=0.8, label=f"Ausgangstemperatur ({T_GROUND_INITIAL_C} °C)")
    ax.axhline(0, color="red", linestyle=":", linewidth=0.8, label="0 °C (Gefriergrenze)")

    ax.set_xlabel("Zeit [Jahre]")
    ax.set_ylabel("Temperatur [°C]")
    ax.set_title(f"BTES – Temperatur an Sonde ({cx:.0f}/{cy:.0f}) · Szenario: {args.scenario}")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    plt.tight_layout()

    out_png = out_dir / f"temperature_bhe_{args.scenario}.png"
    plt.savefig(out_png, dpi=150)
    print(f"Gespeichert: {out_png}")
    plt.show()


if __name__ == "__main__":
    main()

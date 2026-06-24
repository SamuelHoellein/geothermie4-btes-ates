#!/usr/bin/env python3
# coding: utf-8
"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  Saisonaler BTES-Solarspeicher  —  OpenGeoSys HEAT_TRANSPORT_BHE            ║
╚══════════════════════════════════════════════════════════════════════════════╝

ÜBERBLICK — Was macht dieses Skript?
──────────────────────────────────────
Dieses Skript berechnet und simuliert einen saisonalen Erdwärmesondenspeicher
(BTES = Borehole Thermal Energy Storage), der mit einer Solarthermieanlage
kombiniert wird:

  Sommer:  Solarkollektoren erzeugen Wärme  →  Wärme wird im Erdreich
           gespeichert (BTES wird "geladen")
  Winter:  Wärmepumpe entnimmt Wärme aus dem Erdreich  →  Gebäude wird geheizt
           (BTES wird "entladen")

Das Skript besteht aus folgenden Schritten (in dieser Reihenfolge):

  Schritt 1  Energiebilanz berechnen
             Aus den CSV-Dateien (Solar + Heizwärmebedarf) wird für jeden
             Monat die Netto-Wärme ermittelt, die in den BTES-Sonden
             ein- bzw. ausgekoppelt werden muss.

  Schritt 2  Mesherzeugung (gmsh)
             Das Untergrundmodell wird als 3D-Finite-Elemente-Netz (Mesh)
             erstellt: Erdreich in geologischen Schichten + 25 Erdwärmesonden
             als 1D-Linien darin.

  Schritt 3  Mesh-Konvertierung  (.msh  →  .vtu)
             Das gmsh-Format wird in das VTK-Format umgewandelt, das OGS liest.

  Schritt 4  OGS-Projektdatei erstellen  (.prj)
             Alle Simulationsparameter werden in einer XML-Datei für OGS
             zusammengefasst.

  Schritt 5  OGS starten
             OpenGeoSys löst die Wärmetransportgleichungen und erzeugt
             Ergebnisdateien (.vtu/.pvd) für die Visualisierung.

ORIENTIERUNG AN DEN VORÜBUNGEN
────────────────────────────────
Dieses Skript basiert auf `btes/ex2_3d/btes_3d_bhe.py` aus dem Übungsrepository.
Die Mesh-Erzeugung (build_mesh), die PRJ-Erzeugung (build_prj) und alle
Hilfsfunktionen (_bhe_xml, _bhe_positions, _layer_stack, build_power_curve,
msh2vtu) wurden von dort übernommen und für den saisonalen Betrieb angepasst.
Der Hauptunterschied: Die Leistungskurve wird hier automatisch aus den
Eingangsdaten berechnet, statt manuell im CONFIG eingetragen zu werden.

VERWENDUNG
───────────
  python btes_seasonal.py                             # Szenario: vor Sanierung
  python btes_seasonal.py --scenario after_renovation # nach Sanierung
  python btes_seasonal.py --no-mesh                   # Mesh überspringen
  python btes_seasonal.py --no-run                    # nur Mesh + PRJ, kein OGS

ERGEBNIS
─────────
Nach dem Lauf liegen in  out/before_renovation/  (oder after_renovation/):
  *.vtu / *.pvd  — Temperaturfelder für jeden Zeitschritt
Diese können in ParaView (kostenfrei, paraview.org) visualisiert werden.
Außerdem kann plot_results.py aus btes/ex2_3d/ für einfache Plots verwendet
werden (VTU-Prefix im Skript ggf. anpassen).
"""

# ══════════════════════════════════════════════════════════════════════════════
#  IMPORTS  —  externe Bibliotheken laden
# ══════════════════════════════════════════════════════════════════════════════
# Python selbst kann nicht alles. Über "import" lädt man zusätzliche
# Bibliotheken (Pakete) nach. Diese wurden mit "pip install -r requirements.txt"
# installiert.

from __future__ import annotations  # Kompatibilität für ältere Python-Versionen

import argparse   # Verarbeitet Kommandozeilenargumente (--scenario, --no-mesh ...)
import csv as _csv  # Liest CSV-Dateien (Comma/Semicolon Separated Values)
import shutil     # Hilfsfunktionen für das Betriebssystem, z. B. "where is ogs.exe?"
import subprocess # Startet externe Programme (hier: ogs.exe) aus Python heraus
import sys        # Informationen über den laufenden Python-Prozess (Pfad, Exit-Code)
from pathlib import Path              # Plattformunabhängige Dateipfade (Windows/Linux/Mac)
from xml.etree import ElementTree as ET  # Erstellt XML-Dateien (.prj für OGS)

import gmsh    # Erzeugt Finite-Elemente-Netze (Meshes) für 3D-Geometrien
import numpy as np  # Numerische Berechnungen mit Arrays (Vektoren, Matrizen)


# ══════════════════════════════════════════════════════════════════════════════
#  NUTZEREINSTELLUNGEN  —  alle anpassbaren Größen hier ändern
# ══════════════════════════════════════════════════════════════════════════════
# Diese Konstanten steuern das gesamte Modell. Änderungen hier wirken sich
# automatisch auf Energiebilanz, Mesh, PRJ und Simulation aus.

# ── Systemparameter ───────────────────────────────────────────────────────────
COP                  = 4.0    # Jahresarbeitszahl der Wärmepumpe [-]
                               # COP = gelieferte Wärme / eingesetzte Elektrizität
                               # COP=4: 1 kWh Strom  → 4 kWh Wärme ins Gebäude
                               #        davon 3 kWh aus dem Speicher, 1 kWh Strom

T_INJECTION_C        = 65.0   # Einspeis-Temperatur Solar → BTES [°C]
                               # Wie heiß ist das Solarfluid beim Laden?
T_EXTRACT_MIN_C      =  5.0   # Mindest-Entnahmetemperatur BTES → WP [°C]
                               # Unterhalb dieser Temperatur arbeitet die WP
                               # nicht mehr effizient. Wird im Postprocessing geprüft.
T_GROUND_INITIAL_C   = 10.0   # Anfangstemperatur Untergrund [°C]
                               # Annahme: ungestörte Untergrundtemperatur

# ── Szenario-Voreinstellung (per --scenario überschreibbar) ───────────────────
DEFAULT_SCENARIO = "before_renovation"   # "before_renovation" oder "after_renovation"
                                          # Wählt die Spalte in Heizwaermebedarfe.csv

# ── Sondenfeld — Auslegung ────────────────────────────────────────────────────
N_BOREHOLES_TOTAL       = 100   # Gesamtzahl Sonden im realen System
                                 # Bestimmt die Leistung pro Sonde:
                                 # P_pro_Sonde = P_gesamt / N_BOREHOLES_TOTAL
N_BHE_X                 =   5  # Simuliertes Feld: Sonden in x-Richtung
N_BHE_Y                 =   5  # Simuliertes Feld: Sonden in y-Richtung
                                # → 5×5 = 25 Sonden werden simuliert.
                                # Diese repräsentieren das zentrale Verhalten
                                # des realen 100-Sonden-Felds.
BOREHOLE_SPACING_M      =  8.0  # Sondenabstand [m]
BOREHOLE_DEPTH_TOP_M    =  5.0  # Sondenkopf unter Gelände [m] (inaktiver Bereich oben)
BOREHOLE_DEPTH_BOTTOM_M = 105.0 # Sondenfuß unter Gelände [m] → 100 m aktive Sondenlänge
DOMAIN_BUFFER_M         = 25.0  # Seitlicher Puffer außerhalb des Felds [m]
                                 # Verhindert, dass die Randbedingungen die
                                 # Wärmeplume im Sondenfeld beeinflussen.
DOMAIN_DEPTH_BUFFER_M   = 25.0  # Tiefe unterhalb des Sondenfußes [m]

# ── Untergrund — Materialparameter ────────────────────────────────────────────
# Die gemessenen Werte in Ground.csv sind "bulk"-Größen (Gestein + Porenwasser
# zusammen). OGS benötigt aber Feststoff- und Fluidphase getrennt.
# Umrechnung: cv_bulk = (1-φ)·ρ_s·cp_s + φ·ρ_f·cp_f
#             → cp_s = (cv_bulk - φ·ρ_f·cp_f) / ((1-φ)·ρ_s)
RHO_SOLID_KG_M3         = 2650.0  # Feststoffdichte [kg/m³] (typisch für Fels/Kristallin)
POROSITY_DEFAULT        = 0.02    # Porosität φ [-] (2%, typisch für Festgestein)
PERMEABILITY_DEFAULT_M2 = 1.0e-18 # Permeabilität [m²] (extrem niedrig = kein Grundwasserfluss)
                                   # Im BTES-Modell ist Wärmeleitung dominant, kein Fließen.

# ── Simulation — Zeitsteuerung ────────────────────────────────────────────────
N_YEARS      = 3           # Anzahl zu simulierender Jahre (= Lade-/Entladezyklen)
                            # Nach ~3-5 Jahren erreicht der Speicher einen periodischen
                            # Gleichgewichtszustand (jeder Zyklus gleicht dem vorherigen).
DT_SECONDS   = 7 * 86400.0 # Zeitschrittweite [s]: 7 Tage × 86400 s/Tag = 1 Woche
                            # Kleinere Schritte → genauere Ergebnisse, aber längere Laufzeit
OUTPUT_EVERY = 1            # Ergebnisdatei (.vtu) alle N Zeitschritte speichern

# ── BHE-Hardware: 1U-Sonde (U-Rohr-Erdwärmesonde) ────────────────────────────
# Eine 1U-Sonde besteht aus: Bohrloch → Verfüllmasse (Grout) → U-Rohr mit
# Vorlauf und Rücklauf. Das Fluid (Wasser-Glykol-Gemisch) zirkuliert im U-Rohr.
BHE_BOREHOLE_DIAMETER_M    = 0.15   # Bohrloch-Durchmesser [m]
BHE_PIPE_OUTER_DIAMETER_M  = 0.032  # Rohr-Außendurchmesser [m]
BHE_PIPE_WALL_THICKNESS_M  = 0.003  # Rohrwanddicke [m]
BHE_PIPE_WALL_LAMBDA       = 0.4    # Wärmeleitfähigkeit Rohrwand [W/mK] (PE-Rohr)
BHE_PIPE_DISTANCE_M        = 0.05   # Achsabstand Vor- und Rücklaufrohr [m]
BHE_GROUT_DENSITY          = 2190.0 # Dichte Verfüllmasse [kg/m³]
BHE_GROUT_POROSITY         = 0.0    # Porosität Verfüllmasse (dicht, kein Porenraum)
BHE_GROUT_CP               = 1735.0 # Wärmekapazität Verfüllmasse [J/kgK]
BHE_GROUT_LAMBDA           = 2.3    # Wärmeleitfähigkeit Verfüllmasse [W/mK]
BHE_REFR_DENSITY           = 1052.0 # Dichte Solefluid [kg/m³] (Wasser-Glykol-Mischung)
BHE_REFR_VISCOSITY         = 0.0052 # Dynamische Viskosität Solefluid [Pa·s]
BHE_REFR_CP                = 3795.0 # Wärmekapazität Solefluid [J/kgK]
BHE_REFR_LAMBDA            = 0.48   # Wärmeleitfähigkeit Solefluid [W/mK]
BHE_REFR_T_REF_K           = 293.15 # Referenztemperatur Solefluid [K] (= 20 °C)
BHE_FLOW_RATE_KG_S         = 0.2    # Massenstrom je Sonde [kg/s]
                                     # Bestimmt die Temperaturdifferenz Vor-/Rücklauf:
                                     # ΔT = P / (ṁ · cp) → bei 5 kW: ΔT ≈ 7 K

# ── Mesh — Elementgrößen ──────────────────────────────────────────────────────
# Nahe den Sonden brauchen wir feine Elemente (große Temperaturgradienten),
# weiter weg reichen grobe Elemente (wenig ändert sich dort).
MESH_SIZE_NEAR_M   = 1.5   # Elementgröße nahe den Sonden [m]
MESH_SIZE_FAR_M    = 12.0  # Elementgröße weit weg von den Sonden [m]
MESH_RADIUS_NEAR_M = 8.0   # Radius, bis zu dem feine Elemente erzwungen werden [m]
MESH_RADIUS_FAR_M  = 30.0  # Radius, ab dem grobe Elemente verwendet werden [m]

# ── Gleichungslöser (Solver) ──────────────────────────────────────────────────
# OGS löst intern ein großes Gleichungssystem (Ax = b) iterativ.
# Diese Einstellungen steuern Genauigkeit und maximale Iterationszahl.
LINEAR_TOL     = 1.0e-12  # Toleranz für den linearen Löser (Abbruchkriterium)
LINEAR_ITER    = 10000    # Maximale Iterationen des linearen Lösers
NONLINEAR_ITER = 20       # Maximale Iterationen des nichtlinearen Lösers (Picard)

# ── Fluid im Untergrund (Porenwasser) ─────────────────────────────────────────
FLUID_RHO    = 1000.0  # Dichte Wasser [kg/m³]
FLUID_VISC   = 1.0e-3  # Viskosität Wasser [Pa·s]
FLUID_CP     = 4180.0  # Wärmekapazität Wasser [J/kgK]
FLUID_LAMBDA = 0.6     # Wärmeleitfähigkeit Wasser [W/mK]


# ══════════════════════════════════════════════════════════════════════════════
#  DATEIPFADE  —  wo liegen die Eingangsdaten?
# ══════════════════════════════════════════════════════════════════════════════
# Python-Skripte können von verschiedenen Verzeichnissen aus gestartet werden.
# Damit die CSV-Dateien immer gefunden werden, berechnen wir den absoluten Pfad
# relativ zur Position *dieses* Skripts — egal, aus welchem Ordner es gestartet wird.
#
# __file__  =  der absolute Pfad zu dieser Skriptdatei selbst
# .resolve() macht den Pfad absolut (entfernt ".." usw.)
# .parent    geht ein Verzeichnis nach oben
#
# Verzeichnisstruktur:
#   geothermie4-btes-ates/          ← _ROOT (3× .parent von dieser Datei)
#     Data Input/                   ← DATA_DIR
#       Ground.csv
#       Solarthermie.csv
#       Heizwaermebedarfe.csv
#     btes/
#       seasonal/
#         btes_seasonal.py          ← __file__ (diese Datei)
_ROOT    = Path(__file__).resolve().parent.parent.parent
DATA_DIR = _ROOT / "Data Input"

# Hilfskonstanten für Zeitrechnung
DAY        = 86400.0   # Sekunden pro Tag
MONTH_DAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # Tage je Monat
MONTH_NAMES = ["Jan", "Feb", "Mär", "Apr", "Mai", "Jun",
               "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]


# ══════════════════════════════════════════════════════════════════════════════
#  SCHRITT 1: EINGANGSDATEN LADEN & ENERGIEBILANZ BERECHNEN
# ══════════════════════════════════════════════════════════════════════════════

def _read_csv(path: Path) -> list[dict]:
    """
    Liest eine CSV-Datei und gibt eine Liste von Zeilen-Dictionaries zurück.

    WAS IST EINE CSV-DATEI?
    Eine CSV-Datei (Comma/Semicolon Separated Values) ist eine einfache
    Textdatei, in der Tabellendaten zeilenweise gespeichert sind.
    In deutschen Dateien wird ";" als Trennzeichen verwendet (statt ","),
    weil "," in Deutschland als Dezimaltrennzeichen gilt (z. B. "1,5" statt "1.5").

    WAS IST EIN DICTIONARY?
    Ein Dictionary (dict) ist eine Python-Datenstruktur, die Werte unter
    Schlüsseln speichert:  {"Monat": "Januar", "Solar": "38995"}
    Hier bekommt jede CSV-Zeile ein Dictionary: Spaltenname → Wert.

    Parameter:
        path : Dateipfad als Path-Objekt

    Gibt zurück:
        Liste von Dictionaries, eine pro Datenzeile
        Dezimalkomma ("1,5") wird zu Dezimalpunkt ("1.5") konvertiert,
        damit Python float() daraus machen kann.
    """
    with open(path, encoding="utf-8-sig") as f:
        # "utf-8-sig" verarbeitet die unsichtbare BOM-Markierung, die
        # Windows/Excel manchmal an den Dateianfang schreibt.
        rows = list(_csv.DictReader(f, delimiter=";"))
    # Deutsches Komma "," durch englischen Punkt "." ersetzen (für float())
    return [{k: v.replace(",", ".") for k, v in r.items()} for r in rows]


def load_solar() -> np.ndarray:
    """
    Lädt die 12 monatlichen Solarerträge [kWh] aus Solarthermie.csv.

    Die Datei hat zwei Spalten: Monat | Solarertrag [kWh]
    Der Solarertrag ist bereits für die Gesamtanlage (709 Module × Einzelertrag).

    Gibt zurück:
        numpy-Array mit 12 Werten, einem pro Monat (Januar bis Dezember)
    """
    rows = _read_csv(DATA_DIR / "Solarthermie.csv")
    col  = list(rows[0].keys())[1]   # zweite Spalte = Ertragswerte
    return np.array([float(r[col]) for r in rows])


def load_demand(scenario: str) -> np.ndarray:
    """
    Lädt die 12 monatlichen Heizwärmebedarfe [kWh] aus Heizwaermebedarfe.csv.

    Die Datei hat drei Spalten:
      Monat | Heizwärmebedarf (vor Sanierung) | Heizwärmebedarf (nach Sanierung)
    Der Bedarf enthält bereits Heizwärme + Warmwasser.

    Parameter:
        scenario : "before_renovation" → Spalte 2 (original)
                   "after_renovation"  → Spalte 3 (saniert, ~50 % weniger)

    Gibt zurück:
        numpy-Array mit 12 Monatswerten [kWh]
    """
    rows = _read_csv(DATA_DIR / "Heizwaermebedarfe.csv")
    keys = list(rows[0].keys())
    col  = keys[1] if scenario == "before_renovation" else keys[2]
    return np.array([float(r[col]) for r in rows])


def compute_monthly_powers(scenario: str) -> tuple[list[float], dict]:
    """
    Berechnet die monatliche Nettoleistung [W] je Erdwärmesonde.

    PHYSIKALISCHER HINTERGRUND — Energiebilanz
    ────────────────────────────────────────────
    Die Wärmepumpe deckt den Heizwärmebedarf E_demand:
      E_demand = E_aus_BTES  +  E_elektrizitaet
    Mit COP = E_demand / E_elektrizitaet folgt:
      E_aus_BTES = E_demand · (1 − 1/COP)
    Bei COP = 4:  E_aus_BTES = 0,75 · E_demand

    Die monatliche Nettowärme im BTES (Speicher):
      E_net = E_solar  −  E_aus_BTES
      E_net > 0 → Netto-Laden   (mehr Solar als Entnahme, Sommer)
      E_net < 0 → Netto-Entladen (mehr Entnahme als Solar,   Winter)

    Umrechnung Energie → Leistung je Sonde:
      P_bhe = E_net [kWh] × 1000 [Wh/kWh] / (Tage × 86400 [s/Tag])
              / N_BOREHOLES_TOTAL
    Das gibt die mittlere Leistung in Watt für den jeweiligen Monat.

    Vorzeichen-Konvention (identisch zu btes_3d_bhe.py):
      P > 0  →  Wärme wird in den Untergrund injiziert (Laden)
      P < 0  →  Wärme wird aus dem Untergrund entnommen (Entladen für WP)

    Parameter:
        scenario : "before_renovation" oder "after_renovation"

    Gibt zurück:
        power_W : Liste mit 12 Leistungswerten [W/Sonde]
        info    : Dictionary mit Jahresenergien und Monatswerten für den
                  Ausdruck der Energiebilanz-Tabelle
    """
    solar  = load_solar()         # kWh/Monat, Gesamtanlage
    demand = load_demand(scenario) # kWh/Monat, Heizwärmebedarf

    hp_frac = 1.0 - 1.0 / COP    # Anteil, den die WP aus dem BTES bezieht
                                   # Bei COP=4: hp_frac = 0,75
    hp_ext  = demand * hp_frac    # WP-Entzug aus BTES [kWh/Monat]
    net_kwh = solar - hp_ext      # Netto-Wärme im BTES [kWh/Monat]

    power_W = []
    for net, days in zip(net_kwh, MONTH_DAYS):
        p_total = net * 1000.0 / (days * DAY)      # Gesamtleistung des Systems [W]
        power_W.append(p_total / N_BOREHOLES_TOTAL) # Leistung je Sonde [W]

    # Zusammenfassung für die Ausgabe-Tabelle
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
    """
    Druckt eine formatierte Energiebilanz-Tabelle in die Konsole.
    So kann man vor dem OGS-Lauf prüfen, ob die Eingangsdaten plausibel sind.
    """
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
    cover       = info["solar_kwh_a"] / info["demand_kwh_a"] * 100
    solar_vs_hp = info["solar_kwh_a"] / info["hp_ext_kwh_a"] * 100
    print(f"  Solardeckung Heizlast     : {cover:.1f} %")
    print(f"  Solardeckung WP-Entzug    : {solar_vs_hp:.1f} %")
    print(f"  Restdefizit (Erdwärme)    : {-info['net_kwh_a']/1e3:.0f} MWh/a "
          f"(natürliche Erdwärme deckt den Rest)")
    print(f"  Max. Ladeleistung/Sonde   : +{max(power_W):.0f} W")
    print(f"  Max. Entladeleistung/Sonde: {min(power_W):.0f} W")
    print(f"  Simuliertes Feld          : {N_BHE_X}×{N_BHE_Y} = {N_BHE_X*N_BHE_Y} "
          f"Sonden (von {N_BOREHOLES_TOTAL} gesamt)")
    if info["net_kwh_a"] < 0:
        print(f"\n  HINWEIS: Jährliches BTES-Defizit = {-info['net_kwh_a']/1e3:.0f} MWh/a.")
        print(f"  → Solar deckt nur {solar_vs_hp:.0f}% des WP-Entzugs.")
        print(f"  → Restliche {100-solar_vs_hp:.0f}% kommen aus der natürlichen Erdwärme.")
        print(f"  → Bitte Mindesttemperatur {T_EXTRACT_MIN_C}°C in den Plots prüfen!")
    print("=" * 66)
    print()


# ══════════════════════════════════════════════════════════════════════════════
#  SCHICHT-MODELL AUS GROUND.CSV AUFBAUEN
# ══════════════════════════════════════════════════════════════════════════════

def build_layers_from_csv(domain_depth_m: float) -> list[dict]:
    """
    Liest Ground.csv und erzeugt die OGS-kompatible Schichtliste.

    WAS SIND "SCHICHTEN" IM OGS-MODELL?
    ─────────────────────────────────────
    Das Untergrundmodell ist in übereinanderliegende Schichten (Layers) aufgeteilt.
    Jede Schicht hat eigene Materialeigenschaften (Wärmeleitfähigkeit, Dichte usw.).
    OGS weist jedem Volumenelement des Meshes eine MaterialID zu, die bestimmt,
    welche Schichteigenschaften an diesem Ort gelten.

    Die Schichten werden von OBEN nach UNTEN angegeben (wie in Ground.csv):
      Schicht 0: 0–16 m
      Schicht 1: 16–25 m
      ...
    Die Funktion schneidet alles unterhalb von domain_depth_m ab.

    WAS LIEFERT GROUND.CSV?
    ────────────────────────
    Spalte 1: Tiefe [m] = Oberkante der Schicht
    Spalte 2: lambda [W/mK] = Wärmeleitfähigkeit
    Spalte 3: cv [J/m³K] = volumetrische Wärmekapazität des Materials

    OGS benötigt ρ_s (Dichte) und cp_s (spez. Wärmekapazität) GETRENNT.
    Die volumetrische Wärmekapazität ist definiert als:
      cv = ρ_s · cp_s
      → cp_s = cv / ρ_s

    Parameter:
        domain_depth_m : Gesamttiefe des Simulationsgebiets [m]

    Gibt zurück:
        Liste von Dictionaries, je eines pro Schicht, mit den Schlüsseln:
        name, thickness_m, permeability_m2, porosity, rho_s_kg_m3,
        cp_s_J_kgK, lambda_s_W_mK
    """
    rows    = _read_csv(DATA_DIR / "Ground.csv")
    keys    = list(rows[0].keys())
    depths  = [float(r[keys[0]]) for r in rows]  # Oberkante je Schicht [m]
    lambdas = [float(r[keys[1]]) for r in rows]  # Wärmeleitfähigkeit [W/mK]
    cvs     = [float(r[keys[2]]) for r in rows]  # volumetrische Wärmekapazität [J/m³K]

    rho_s = RHO_SOLID_KG_M3

    layers = []
    for i in range(len(depths)):
        z_top = depths[i]
        if z_top >= domain_depth_m:
            break  # Schichten unterhalb des Simulationsgebiets ignorieren
        # Unterkante: nächste Tiefe aus CSV, oder Ende des Simulationsgebiets
        z_bot     = depths[i + 1] if i + 1 < len(depths) else domain_depth_m
        z_bot     = min(z_bot, domain_depth_m)  # Auf Simulationsgebiet begrenzen
        thickness = z_bot - z_top
        if thickness <= 0.0:
            continue

        lam  = lambdas[i]
        cv   = cvs[i]
        cp_s = cv / rho_s   # spez. Wärmekapazität: cv [J/m³K] / ρ_s [kg/m³] = cp_s [J/kgK]

        layers.append({
            "name"            : f"layer_{i:02d}",
            "thickness_m"     : float(thickness),
            "permeability_m2" : PERMEABILITY_DEFAULT_M2,  # kein Grundwasserfluss im BTES
            "porosity"        : POROSITY_DEFAULT,
            "rho_s_kg_m3"     : rho_s,
            "cp_s_J_kgK"      : float(cp_s),
            "lambda_s_W_mK"   : float(lam),
        })
    return layers


# ══════════════════════════════════════════════════════════════════════════════
#  CONFIG-DICTIONARY ZUSAMMENBAUEN
# ══════════════════════════════════════════════════════════════════════════════

def build_config(scenario: str, out_dir: Path) -> dict:
    """
    Erstellt den zentralen CONFIG-Dictionary aus Nutzereinstellungen + CSV-Daten.

    WAS IST EIN CONFIG-DICTIONARY?
    ────────────────────────────────
    Ein Dictionary (dict) ist eine Python-Datenstruktur mit Schlüssel-Wert-Paaren:
      {"schluessel": wert, "anderer_schluessel": anderer_wert}
    Dicts können verschachtelt sein — ein Wert kann selbst wieder ein Dict sein.

    Der CONFIG-Dict ist die zentrale Datenschnittstelle dieses Skripts:
    Alle nachfolgenden Funktionen (build_mesh, build_prj) lesen aus diesem Dict
    und brauchen keine eigenen Argumente mehr. Das ist das gleiche Prinzip wie
    in den Vorübungen (btes_3d_bhe.py, ates_radial_2d.py).

    UNTERSCHIED ZU DEN VORÜBUNGEN:
    ───────────────────────────────
    In btes_3d_bhe.py ist der CONFIG-Dict am Skriptanfang fest eingetragen.
    Hier wird er aus den CSV-Daten BERECHNET: die Leistungskurve und die
    Schichten werden automatisch aus den Eingangsdaten erzeugt.

    Parameter:
        scenario : "before_renovation" oder "after_renovation"
        out_dir  : Ausgabeverzeichnis als Path-Objekt

    Gibt zurück:
        Vollständiger CONFIG-Dictionary für build_mesh() und build_prj()
    """
    # Energiebilanz → Monatliche Leistung je Sonde
    power_W, info = compute_monthly_powers(scenario)
    print_energy_balance(power_W, info, scenario)

    # Simulationstiefe und Schichten aus Ground.csv
    domain_depth = BOREHOLE_DEPTH_BOTTOM_M + DOMAIN_DEPTH_BUFFER_M
    layers = build_layers_from_csv(domain_depth)

    # Domänengröße: Sondenfeld + Puffer auf allen Seiten
    field_x = (N_BHE_X - 1) * BOREHOLE_SPACING_M   # Ausdehnung des Felds in x [m]
    field_y = (N_BHE_Y - 1) * BOREHOLE_SPACING_M   # Ausdehnung des Felds in y [m]
    size_x  = field_x + 2.0 * DOMAIN_BUFFER_M       # Gesamtdomäne in x [m]
    size_y  = field_y + 2.0 * DOMAIN_BUFFER_M       # Gesamtdomäne in y [m]

    prefix = f"btes_seasonal_{scenario}"  # Dateiname-Präfix für alle Ausgabedateien

    return {
        # ── Geometrie des Simulationsgebiets ──────────────────────────────────
        "domain": {
            "size_x_m": size_x,   # Breite [m]
            "size_y_m": size_y,   # Tiefe  [m] (in y-Richtung)
            "z_base_m": 0.0,      # z-Koordinate des Bodens (intern; 0 = Sohle)
        },
        # ── Geologische Schichten (aus Ground.csv) ─────────────────────────────
        "layers": layers,
        # ── Sonden-Geometrie ───────────────────────────────────────────────────
        "borehole": {
            "depth_top_m"   : BOREHOLE_DEPTH_TOP_M,
            "depth_bottom_m": BOREHOLE_DEPTH_BOTTOM_M,
        },
        # ── Sondenfeld-Layout ──────────────────────────────────────────────────
        "field": {
            "n_x"      : N_BHE_X,
            "n_y"      : N_BHE_Y,
            "spacing_m": BOREHOLE_SPACING_M,
            "positions": None,   # None = gleichmäßiges Raster (kein manuelles Layout)
        },
        # ── Mesh-Feinheit (Elementgrößen) ──────────────────────────────────────
        "mesh": {
            "size_near_field_m"      : MESH_SIZE_NEAR_M,
            "size_far_m"             : MESH_SIZE_FAR_M,
            "field_size_radius_m"    : MESH_RADIUS_NEAR_M,
            "field_size_radius_far_m": MESH_RADIUS_FAR_M,
        },
        # ── Porenfluid-Eigenschaften (Grundwasser) ─────────────────────────────
        "fluid": {
            "rho_ref_kg_m3" : FLUID_RHO,
            "T_ref_K"       : T_GROUND_INITIAL_C + 273.15,  # °C → Kelvin
            "viscosity_Pa_s": FLUID_VISC,
            "cp_J_kgK"      : FLUID_CP,
            "lambda_W_mK"   : FLUID_LAMBDA,
        },
        # ── BHE (Erdwärmesonden) ───────────────────────────────────────────────
        "bhe": {
            "type"    : "1U",   # 1U = einfaches U-Rohr (Standard in Deutschland)
            "borehole": {"diameter_m": BHE_BOREHOLE_DIAMETER_M},
            "pipes": {
                "diameter_outer_m"              : BHE_PIPE_OUTER_DIAMETER_M,
                "wall_thickness_m"              : BHE_PIPE_WALL_THICKNESS_M,
                "wall_thermal_conductivity_W_mK": BHE_PIPE_WALL_LAMBDA,
                "distance_between_pipes_m"      : BHE_PIPE_DISTANCE_M,
                "longitudinal_dispersion_length_m": 0.001,
            },
            "grout": {
                "density_kg_m3"               : BHE_GROUT_DENSITY,
                "porosity"                    : BHE_GROUT_POROSITY,
                "specific_heat_capacity_J_kgK": BHE_GROUT_CP,
                "thermal_conductivity_W_mK"   : BHE_GROUT_LAMBDA,
            },
            "refrigerant": {
                "density_kg_m3"               : BHE_REFR_DENSITY,
                "viscosity_Pa_s"              : BHE_REFR_VISCOSITY,
                "specific_heat_capacity_J_kgK": BHE_REFR_CP,
                "thermal_conductivity_W_mK"   : BHE_REFR_LAMBDA,
                "reference_temperature_K"     : BHE_REFR_T_REF_K,
            },
            "control": {
                "type"          : "Power",   # OGS 6.5.8: zeitvariable Leistung
                "flow_rate_kg_s": BHE_FLOW_RATE_KG_S,
            },
        },
        # ── Anfangsbedingungen ─────────────────────────────────────────────────
        "initial": {
            "T_K"                       : T_GROUND_INITIAL_C + 273.15,
            "p_Pa"                      : 0.0,
            "T_surface_K"               : T_GROUND_INITIAL_C + 273.15,
            "geothermal_gradient_K_per_m": 0.0,  # Auf 0 gesetzt → homogene Anfangstemperatur
        },
        # ── Betriebszyklen (monatliches Lastprofil) ────────────────────────────
        "cycles": {
            "n_cycles"       : N_YEARS,   # Anzahl Jahre = Anzahl Wiederholungen des 12-Monats-Profils
            "monthly_power_W": power_W,   # 12 Leistungswerte [W/Sonde], aus Energiebilanz
            "ramp_days"      : 3.0,       # Sanfter Übergang zwischen Monaten (3 Tage Rampe)
        },
        # ── Zeitsteuerung ──────────────────────────────────────────────────────
        "time": {
            "dt_seconds"          : DT_SECONDS,
            "output_every_n_steps": OUTPUT_EVERY,
        },
        # ── Ausgabe ────────────────────────────────────────────────────────────
        "output": {
            "prefix"   : prefix,
            "out_dir"  : str(out_dir),
            "variables": ["temperature_soil"],
        },
        # ── Gleichungslöser ────────────────────────────────────────────────────
        "solver": {
            "linear_tol"    : LINEAR_TOL,
            "linear_iter"   : LINEAR_ITER,
            "nonlinear_iter": NONLINEAR_ITER,
        },
    }


# ══════════════════════════════════════════════════════════════════════════════
#  SCHRITT 2: MESH ERZEUGEN (gmsh)
# ══════════════════════════════════════════════════════════════════════════════
#
# WAS IST EIN MESH (Finite-Elemente-Netz)?
# ─────────────────────────────────────────
# Die Finite-Elemente-Methode (FEM) teilt das Simulationsgebiet in viele kleine
# Elemente (Tetraeder oder Hexaeder in 3D) auf. In jedem Element wird die
# Differentialgleichung (Wärmetransport) vereinfacht gelöst.
# Das Netz aus diesen Elementen heißt Mesh.
#
# Feine Elemente nahe den Sonden → hohe Genauigkeit dort (große Temperaturgradienten)
# Grobe Elemente weiter weg      → schnellere Berechnung (wenig ändert sich dort)
#
# WAS MACHT GMSH?
# ───────────────
# gmsh ist eine freie Software zur Mesh-Erzeugung. Wir beschreiben die Geometrie
# (Boxen für Schichten, Linien für Sonden) und gmsh erzeugt das FE-Netz automatisch.
# Ausgabe: eine .msh-Datei mit Knoten und Elementen.
#
# ORIENTIERUNG AN DER VORÜBUNG:
# ──────────────────────────────
# build_mesh() ist fast identisch zu btes/ex2_3d/btes_3d_bhe.py.
# Der einzige Unterschied: Die Schichten kommen aus build_layers_from_csv()
# statt aus einem festen CONFIG-Dict.


def msh2vtu(filename, output_path, output_prefix, dim, reindex=True):
    """
    Konvertiert die gmsh-Meshdatei (.msh) in das VTK-Format (.vtu).

    WAS IST EINE .msh-DATEI?
    Das proprietäre Ausgabeformat von gmsh: enthält Knoten, Elemente und
    Physical Groups als Textdatei. OGS kann dieses Format NICHT direkt lesen.

    WAS IST EINE .vtu-DATEI?
    VTU (VTK Unstructured Grid) ist das Standardformat für FEM-Netze und
    Simulationsergebnisse. OGS liest und schreibt .vtu-Dateien.
    VTU-Dateien können auch in ParaView visualisiert werden.

    WAS PASSIERT BEI DER KONVERTIERUNG?
    ─────────────────────────────────────
    ogstools.Meshes.from_gmsh() liest die .msh-Datei und erzeugt:
      - *_domain.vtu        : Das gesamte 3D+1D-Netz (Bodenvolumen + Sondenlinien)
      - *_physical_group_top.vtu       : Oberfläche (für Randbedingung)
      - *_physical_group_bottom.vtu    : Boden      (für Randbedingung)
      - *_physical_group_lateral.vtu   : Seitenflächen
      - *_physical_group_bhe_XX.vtu    : je eine Datei pro Sonde

    OGS braucht Boden, Oberfläche und Seiten als separate Dateien, weil dort
    die Randbedingungen (Temperatur, Druck) gesetzt werden.

    Parameter:
        filename      : Pfad zur .msh-Datei
        output_path   : Zielordner
        output_prefix : Dateiname-Präfix (z. B. "btes_seasonal_before_renovation")
        dim           : Dimensionen, die konvertiert werden [3, 1] = 3D-Volumen + 1D-Linien
        reindex       : MaterialIDs neu nummerieren (wichtig für OGS)
    """
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
    """
    Berechnet die (x, y)-Koordinaten aller Sonden aus dem Raster-Layout.
    Das Raster ist zentriert im Koordinatenursprung (0, 0).
    Übernommen unverändert aus btes_3d_bhe.py.
    """
    fld = cfg["field"]
    if fld.get("positions"):
        return [tuple(p) for p in fld["positions"]]
    nx, ny, sp = fld["n_x"], fld["n_y"], fld["spacing_m"]
    xs = (np.arange(nx) - (nx - 1) / 2.0) * sp
    ys = (np.arange(ny) - (ny - 1) / 2.0) * sp
    return [(float(x), float(y)) for y in ys for x in xs]


def _z_for_depth(cfg: dict, depth: float) -> float:
    """
    Rechnet eine Tiefe unter Gelände [m] in eine interne z-Koordinate um.
    Intern gilt: z=0 an der Domänensohle, z=z_top an der Oberfläche.
    Die CONFIG-Schichten zählen Tiefen von der Oberfläche nach unten.
    """
    z_base = cfg["domain"]["z_base_m"]
    z_top  = z_base + sum(float(L["thickness_m"]) for L in cfg["layers"])
    return z_top - depth


def _layer_stack(cfg: dict):
    """
    Wandelt die oben-nach-unten-Schichtliste des CONFIG in eine
    unten-nach-oben-Liste mit absoluten z-Koordinaten um.
    gmsh baut die Geometrie von unten nach oben (physikalisch natürlich).
    Gibt zurück: (layers_bottom_up, z_top) — z_top = z-Koordinate der Geländeoberfläche.
    """
    z_base = cfg["domain"].get("z_base_m", 0.0)
    bot_up = list(reversed(cfg["layers"]))  # Reihenfolge umkehren: unten → oben
    z = z_base
    out = []
    for L in bot_up:
        h = float(L["thickness_m"])
        out.append({**L, "z_low": z, "z_high": z + h})
        z += h
    return out, z


def build_mesh(cfg: dict, out_dir: Path) -> Path:
    """
    Erzeugt das 3D-Finite-Elemente-Netz mit gmsh.

    GEOMETRIE:
    ───────────
    Das Modell besteht aus:
      - Mehreren übereinanderliegenden Quaderboxen (eine pro geologische Schicht)
      - 25 vertikalen 1D-Linien, die die Erdwärmesonden repräsentieren
        (eingebettet in die 3D-Schichten)

    FRAGMENT-OPERATION:
    ────────────────────
    gmsh.occ.fragment() verbindet die Schicht-Boxen und Sondenlinien so, dass
    alle Elemente konform aneinanderstoßen (keine überlappenden oder lückenhaften
    Elemente). Die Linien werden an den Schichtgrenzen in Segmente aufgeteilt.

    PHYSICAL GROUPS:
    ─────────────────
    Jeder Volumen- und Flächengruppe wird ein Name und eine ID gegeben:
      - Volumen layer_00, layer_01, ... → MaterialID 0, 1, ...
        (OGS weiß so, welche Materialeigenschaften wo gelten)
      - Linien bhe_00, bhe_01, ...      → MaterialID nach den Schichten
        (OGS identifiziert die BHE-Elemente)
      - Flächen top, bottom, lateral    → Randbedingungen

    DOMÄNENRÄNDER (Randbedingungen-Flächen):
    ─────────────────────────────────────────
    In der FEM müssen an allen Rändern des Simulationsgebiets physikalische
    Bedingungen gesetzt werden:
      - "top"     = Geländeoberfläche   → Temperatur = T0 (konstant, Wärmesenke oben)
      - "bottom"  = Domänensohle        → Temperatur = T0 (konstant, Wärme kommt von unten)
      - "lateral" = Seitenflächen       → keine Randbedingung gesetzt = adiabat
        (kein Wärmefluss durch die Seiten, weil der Puffer groß genug ist)
    Diese Flächen werden als separate .vtu-Dateien gespeichert und in der
    OGS-Projektdatei als Randbedingungen referenziert.

    MESH-VERFEINERUNG:
    ───────────────────
    gmsh verwendet eine Distance/Threshold-Feld-Kombination:
      - Nahe den Sondenachsen: Elemente ≤ MESH_SIZE_NEAR_M (fein)
      - Weit weg: Elemente ≤ MESH_SIZE_FAR_M (grob)
    Das spart Rechenzeit ohne Genauigkeitsverlust im Fernfeld.

    Übernommen aus btes_3d_bhe.py, angepasst für dynamische Schichtliste.
    """
    prefix   = cfg["output"]["prefix"]
    msh_path = out_dir / f"{prefix}.msh"
    Lx       = cfg["domain"]["size_x_m"]
    Ly       = cfg["domain"]["size_y_m"]
    z_base   = cfg["domain"]["z_base_m"]
    layers, z_top = _layer_stack(cfg)

    bhe_pos   = _bhe_positions(cfg)
    z_bhe_top = z_top - cfg["borehole"]["depth_top_m"]   # z-Koord. Sondenkopf
    z_bhe_bot = z_top - cfg["borehole"]["depth_bottom_m"] # z-Koord. Sondenfuß
    m         = cfg["mesh"]
    x0, y0    = -Lx / 2.0, -Ly / 2.0  # linke untere Ecke der Domäne (zentriert)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # gmsh-Ausgabe unterdrücken
    gmsh.model.add("btes_seasonal")

    # Eine Box (Quader) pro geologische Schicht
    layer_boxes = [
        gmsh.model.occ.addBox(x0, y0, L["z_low"], Lx, Ly, L["z_high"] - L["z_low"])
        for L in layers
    ]
    # Eine vertikale Linie pro Sonde
    lines = []
    for x, y in bhe_pos:
        p_top = gmsh.model.occ.addPoint(x, y, z_bhe_top)
        p_bot = gmsh.model.occ.addPoint(x, y, z_bhe_bot)
        lines.append(gmsh.model.occ.addLine(p_top, p_bot))

    # Schichten und Sondenlinien "verschmelzen" — konforme gemeinsame Geometrie
    gmsh.model.occ.fragment([(3, b) for b in layer_boxes],
                             [(1, l) for l in lines])
    gmsh.model.occ.synchronize()

    # Jedes 3D-Volumenelement der richtigen Schicht zuweisen (über z-Mittelpunkt)
    vol_layer = {i: [] for i in range(len(layers))}
    for dim, tag in gmsh.model.getEntities(3):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (bb[2] + bb[5])  # z-Mittelpunkt des Elements
        for i, L in enumerate(layers):
            if L["z_low"] - 1e-6 <= zc <= L["z_high"] + 1e-6:
                vol_layer[i].append(tag)
                break

    # Alle 1D-Liniensegmente der zugehörigen Sonde zuweisen
    pos = np.array(bhe_pos)
    seg_bhe  = {i: [] for i in range(len(bhe_pos))}
    all_segs = []
    for dim, tag in gmsh.model.getEntities(1):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        dx, dy, dz = bb[3]-bb[0], bb[4]-bb[1], bb[5]-bb[2]
        xc = 0.5 * (bb[0] + bb[3])
        yc = 0.5 * (bb[1] + bb[4])
        zc = 0.5 * (bb[2] + bb[5])
        # Nur vertikale Linien im Sondenbereich berücksichtigen
        if dx < 1e-6 and dy < 1e-6 and dz > 1e-9:
            if (z_bhe_bot - 1e-3) <= zc <= (z_bhe_top + 1e-3):
                # Nächste Sonde finden (minimale horizontale Distanz)
                d = np.hypot(pos[:, 0] - xc, pos[:, 1] - yc)
                j = int(np.argmin(d))
                if d[j] < 1e-3:
                    seg_bhe[j].append(tag)
                    all_segs.append(tag)
    missing = [i for i, s in seg_bhe.items() if not s]
    if missing:
        raise RuntimeError(f"BHE ohne Liniensegment: {missing}")

    # Physical Groups: Schicht-Volumen → MaterialID 0..L-1,
    # danach BHE-Linien → MaterialID L..L+N-1
    pg = 1
    for i, L in enumerate(layers):
        gmsh.model.addPhysicalGroup(3, vol_layer[i], tag=pg, name=L["name"])
        pg += 1
    for i in range(len(bhe_pos)):
        gmsh.model.addPhysicalGroup(1, seg_bhe[i], tag=pg, name=f"bhe_{i:02d}")
        pg += 1

    # Domänenränder identifizieren und als Physical Groups setzen
    # (werden als .vtu-Dateien für OGS-Randbedingungen benötigt)
    top_faces, bot_faces, lat_faces = [], [], []
    for dim, tag in gmsh.model.getEntities(2):
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        zc = 0.5 * (bb[2] + bb[5])
        if (bb[3]-bb[0]) >= 0.9*Lx and abs(zc - z_top) < 1e-6:
            top_faces.append(tag)    # Oberfläche
        elif (bb[3]-bb[0]) >= 0.9*Lx and abs(zc - z_base) < 1e-6:
            bot_faces.append(tag)    # Sohle
        else:
            on_outer = (
                (abs(bb[0]-x0)       < 1e-6 and abs(bb[3]-x0)       < 1e-6) or
                (abs(bb[0]-(x0+Lx))  < 1e-6 and abs(bb[3]-(x0+Lx))  < 1e-6) or
                (abs(bb[1]-y0)       < 1e-6 and abs(bb[4]-y0)       < 1e-6) or
                (abs(bb[1]-(y0+Ly))  < 1e-6 and abs(bb[4]-(y0+Ly))  < 1e-6)
            )
            if on_outer:
                lat_faces.append(tag)  # Seitenflächen

    gmsh.model.addPhysicalGroup(2, top_faces, tag=200, name="top")
    gmsh.model.addPhysicalGroup(2, bot_faces, tag=201, name="bottom")
    if lat_faces:
        gmsh.model.addPhysicalGroup(2, lat_faces, tag=202, name="lateral")

    # Mesh-Verfeinerung: fein nahe den Sonden, grob im Fernfeld
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "CurvesList", all_segs)
    f_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_thr, "InField",  f_dist)
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMin",  m["size_near_field_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "SizeMax",  m["size_far_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMin",  m["field_size_radius_m"])
    gmsh.model.mesh.field.setNumber(f_thr, "DistMax",  m["field_size_radius_far_m"])
    gmsh.model.mesh.field.setAsBackgroundMesh(f_thr)

    # An den Sondenknoten: kleinste Elementgröße erzwingen
    bhe_pts = []
    for tag in all_segs:
        for d, t in gmsh.model.getBoundary([(1, tag)], oriented=False):
            if d == 0:
                bhe_pts.append((0, t))
    if bhe_pts:
        gmsh.model.mesh.setSize(list(set(bhe_pts)), 0.4)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints",         1)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion",             2.2)
    gmsh.model.mesh.generate(3)  # 3D-Netz erzeugen
    gmsh.write(str(msh_path))    # .msh-Datei schreiben
    gmsh.finalize()
    return msh_path


# ══════════════════════════════════════════════════════════════════════════════
#  LEISTUNGSKURVE  (zeitabhängiges Betriebsprofil)
# ══════════════════════════════════════════════════════════════════════════════

def build_power_curve(cfg: dict):
    """
    Erzeugt die zeitabhängige Leistungskurve [W/Sonde] für OGS.

    OGS arbeitet mit einer Tabelle (Zeit, Leistung):
      [(t0, P0), (t1, P1), (t2, P2), ...]
    Zwischen den Stützpunkten interpoliert OGS linear.

    RAMPEN:
    ────────
    Zwischen zwei Monaten mit unterschiedlicher Leistung gibt es eine sanfte
    Übergangsrampe (ramp_days = 3 Tage). Sprünge (ohne Rampe) würden numerische
    Instabilitäten im Löser verursachen.

    Das Profil wird n_cycles-mal wiederholt (= N_YEARS Jahre).

    Gibt zurück:
        (t_gesamt, (times_array, values_array))
        oder None, wenn kein monatliches Profil gesetzt ist.
    """
    monthly = cfg["cycles"].get("monthly_power_W")
    if monthly is None:
        return None
    assert len(monthly) == 12, "monthly_power_W muss genau 12 Werte haben."
    n         = cfg["cycles"]["n_cycles"]        # Anzahl Wiederholungen (= Jahre)
    ramp      = max(60.0, cfg["cycles"]["ramp_days"] * DAY)  # Rampendauer [s]
    month_dur = 365.25 / 12.0 * DAY             # mittlere Monatsdauer [s] ≈ 30,44 Tage
    times = [0.0]
    vals  = [0.0]
    t_now = 0.0
    for _ in range(n):
        for P in monthly:
            t_now += ramp
            times.append(t_now); vals.append(float(P))
            hold = max(0.0, month_dur - ramp)
            if hold > 0.0:
                t_now += hold
                times.append(t_now); vals.append(float(P))
    t_now += ramp
    times.append(t_now); vals.append(0.0)
    return t_now, (np.array(times), np.array(vals))


# ══════════════════════════════════════════════════════════════════════════════
#  SCHRITT 4: OGS-PROJEKTDATEI (.prj) ERSTELLEN
# ══════════════════════════════════════════════════════════════════════════════
#
# WAS IST EINE .prj-DATEI?
# ─────────────────────────
# Eine .prj-Datei ist eine XML-Datei, die alle Informationen für die OGS-Simulation
# enthält: Meshes, Prozesstyp, Materialien, Randbedingungen, Zeitsteuerung,
# Anfangsbedingungen und Gleichungslöser.
#
# WAS IST XML?
# ─────────────
# XML (eXtensible Markup Language) ist ein Textformat für strukturierte Daten.
# Es verwendet Tags wie <name>Inhalt</name> oder verschachtelte Tags:
#   <parameter>
#     <name>T0</name>
#     <type>Constant</type>
#     <value>283.15</value>
#   </parameter>
# Python's ElementTree-Bibliothek erzeugt solche Strukturen programmatisch.
#
# WAS MACHT OGS MIT DER .prj-DATEI?
# ────────────────────────────────────
# OGS (OpenGeoSys) ist ein Finite-Elemente-Simulator für geothermische und
# hydrogeologische Prozesse. Es liest die .prj-Datei und:
#   1. Lädt das Mesh (.vtu-Dateien)
#   2. Setzt Anfangsbedingungen (Temperatur = T0 überall)
#   3. Setzt Randbedingungen (Temperatur oben/unten = T0)
#   4. Löst für jeden Zeitschritt die Wärmetransportgleichungen
#   5. Schreibt die Ergebnisse als .vtu-Dateien
#
# ORIENTIERUNG AN DER VORÜBUNG:
# ──────────────────────────────
# build_prj() und die Hilfsfunktionen _se(), _const_prop(), _indent(), _bhe_xml()
# sind fast identisch zu btes_3d_bhe.py. Der Hauptunterschied:
# Statt "PowerCurveConstantFlow" wird "Power" mit einem "CurveScaled"-Parameter
# verwendet (OGS 6.5.8 hat den alten Typ entfernt).

def _se(parent, tag, text=None, **attrs):
    """
    Hilfsfunktion: Erstellt ein XML-Kindelement unter 'parent'.
    _se steht für 'SubElement'. Kurzform für ET.SubElement() mit optionalem Text.
    """
    el = ET.SubElement(parent, tag, **{k: str(v) for k, v in attrs.items()})
    if text is not None:
        el.text = str(text)
    return el


def _const_prop(parent, name, value):
    """Hilfsfunktion: Erstellt eine <property>-Zeile mit Konstantenwert."""
    p = _se(parent, "property")
    _se(p, "name", name)
    _se(p, "type", "Constant")
    _se(p, "value", value)


def _indent(elem, level=0):
    """Hilfsfunktion: Fügt Einrückung in das XML-Element ein (schönere Ausgabe)."""
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
    """
    Erstellt den <borehole_heat_exchangers>-Block in der OGS-Projektdatei.

    Für jede der n_bhe Sonden wird ein identischer Block geschrieben mit:
      - Steuerungstyp (Power: zeitvar. Leistung + konstanter Massenstrom)
      - Bohrlochgeometrie
      - Rohrgeometrie (Vor-/Rücklauf des U-Rohrs)
      - Verfüllmasse (Grout)
      - Solefluid (Wasser-Glykol)

    POWER-STEUERUNG MIT CurveScaled (OGS 6.5.8):
    ──────────────────────────────────────────────
    In früheren OGS-Versionen gab es "PowerCurveConstantFlow": direkte Referenz
    auf eine Kurve. Ab OGS 6.5.8 wird "Power" mit einem OGS-Parameter verwendet.
    Der Parameter "bhe_power" ist vom Typ "CurveScaled" und multipliziert die
    Leistungskurve mit dem Skalierungsfaktor "unity" = 1.
    Dies ergibt die gleiche zeitvariable Leistung wie zuvor.

    Übernommen aus btes_3d_bhe.py, angepasst für OGS 6.5.8.
    """
    bhe   = cfg["bhe"]
    L_bhe = cfg["borehole"]["depth_bottom_m"] - cfg["borehole"]["depth_top_m"]
    bhes  = _se(parent, "borehole_heat_exchangers")
    for i in range(n_bhe):
        b = _se(bhes, "borehole_heat_exchanger")
        _se(b, "type", bhe["type"])   # "1U" = einfaches U-Rohr
        # Leistungssteuerung: "Power" + Massenstrom
        ftc      = _se(b, "flow_and_temperature_control")
        flow_vol = bhe["control"]["flow_rate_kg_s"] / bhe["refrigerant"]["density_kg_m3"]
        # flow_vol: OGS erwartet volumetrischen Durchfluss [m³/s], nicht Massenstrom [kg/s]
        _se(ftc, "type",      "Power")
        _se(ftc, "power",     "bhe_power")  # Verweis auf den CurveScaled-Parameter
        _se(ftc, "flow_rate", flow_vol)
        # Bohrloch-Geometrie
        bh = _se(b, "borehole")
        _se(bh, "length",   L_bhe)
        _se(bh, "diameter", bhe["borehole"]["diameter_m"])
        # Rohr-Geometrie (Vor- und Rücklauf identisch)
        pipes = _se(b, "pipes")
        for pname in ("inlet", "outlet"):
            p = _se(pipes, pname)
            _se(p, "diameter",                  bhe["pipes"]["diameter_outer_m"])
            _se(p, "wall_thickness",            bhe["pipes"]["wall_thickness_m"])
            _se(p, "wall_thermal_conductivity", bhe["pipes"]["wall_thermal_conductivity_W_mK"])
        _se(pipes, "distance_between_pipes",         bhe["pipes"]["distance_between_pipes_m"])
        _se(pipes, "longitudinal_dispersion_length", bhe["pipes"]["longitudinal_dispersion_length_m"])
        # Verfüllmasse (Grout)
        gr = _se(b, "grout")
        _se(gr, "density",                bhe["grout"]["density_kg_m3"])
        _se(gr, "porosity",               bhe["grout"]["porosity"])
        _se(gr, "specific_heat_capacity", bhe["grout"]["specific_heat_capacity_J_kgK"])
        _se(gr, "thermal_conductivity",   bhe["grout"]["thermal_conductivity_W_mK"])
        # Solefluid (Wasser-Glykol)
        rf = _se(b, "refrigerant")
        _se(rf, "density",                bhe["refrigerant"]["density_kg_m3"])
        _se(rf, "viscosity",              bhe["refrigerant"]["viscosity_Pa_s"])
        _se(rf, "specific_heat_capacity", bhe["refrigerant"]["specific_heat_capacity_J_kgK"])
        _se(rf, "thermal_conductivity",   bhe["refrigerant"]["thermal_conductivity_W_mK"])
        _se(rf, "reference_temperature",  bhe["refrigerant"]["reference_temperature_K"])


def build_prj(cfg: dict, out_dir: Path, mesh_files: dict) -> Path:
    """
    Schreibt die vollständige OGS-Projektdatei (.prj) als XML.

    AUFBAU DER .prj-DATEI (Abschnitte in Reihenfolge):
    ─────────────────────────────────────────────────────
    1. <meshes>         : Liste aller .vtu-Dateien (Domäne + Ränder + BHE-Linien)
    2. <processes>      : Prozesstyp = HEAT_TRANSPORT_BHE, BHE-Konfiguration
    3. <media>          : Materialeigenschaften je MaterialID (Schichten + BHE)
    4. <time_loop>      : Zeitsteuerung, Zeitschritte, Ausgabeintervall
    5. <parameters>     : Benannte Werte (T0, T0_BHE, bhe_power, unity)
    6. <process_variables>: Anfangs- und Randbedingungen für Temperatur
    7. <nonlinear_solvers>: Picard-Iteration (äußere Schleife)
    8. <linear_solvers>  : BiCGSTAB-Löser (innere Schleife)
    9. <curves>          : Zeitreihe der Sondenleistung [W/Sonde] vs. Zeit [s]

    MATERIALEIGENSCHAFTEN (media):
    ────────────────────────────────
    OGS unterscheidet zwei Phasen im Porenmedium:
      - AqueousLiquid: das Porenfluid (Grundwasser)
      - Solid: das Gesteins-Feststoffgerüst
    Die effektive Wärmeleitfähigkeit wird aus beiden gemischt:
      λ_eff = φ·λ_f + (1-φ)·λ_s  (Porosität × Fluid + Rest × Feststoff)

    RANDBEDINGUNGEN (boundary_conditions):
    ────────────────────────────────────────
    Dirichlet-Randbedingung: T = T0 = konstant an Ober- und Unterseite der Domäne.
    Das entspricht der physikalischen Annahme, dass weit genug ober- und unterhalb
    der Sonden die Temperatur ungestört bleibt.

    Parameter:
        cfg        : CONFIG-Dictionary (aus build_config())
        out_dir    : Ausgabeverzeichnis
        mesh_files : Dictionary {name: Dateiname} für alle .vtu-Dateien

    Gibt zurück:
        Pfad zur erzeugten .prj-Datei
    """
    prefix    = cfg["output"]["prefix"]
    fluid     = cfg["fluid"]
    init      = cfg["initial"]
    sol       = cfg["solver"]
    n_bhe     = len(_bhe_positions(cfg))

    # Leistungskurve aufbauen und Simulationsendzeit bestimmen
    power_curve = build_power_curve(cfg)
    t_end   = power_curve[0]                               # Gesamte Simulationszeit [s]
    n_steps = int(t_end // cfg["time"]["dt_seconds"]) + 1  # Anzahl Zeitschritte

    # XML-Wurzelelement
    root = ET.Element("OpenGeoSysProject")

    # ── 1. Meshes ────────────────────────────────────────────────────────────
    meshes = _se(root, "meshes")
    for key in ("domain", "top", "bottom", "lateral"):
        _se(meshes, "mesh", mesh_files[key])
    for i in range(n_bhe):
        _se(meshes, "mesh", mesh_files[f"bhe_{i:02d}"])

    # ── 2. Process ────────────────────────────────────────────────────────────
    processes = _se(root, "processes")
    proc      = _se(processes, "process")
    _se(proc, "name", "HeatTransportBHE")
    _se(proc, "type", "HEAT_TRANSPORT_BHE")  # OGS-Prozessmodul mit U-Rohr-Modell
    _se(proc, "integration_order", 2)
    # Prozessvariablen: Erdreichtemperatur + je eine BHE-Fluidtemperatur
    pv = _se(proc, "process_variables")
    _se(pv, "process_variable", "temperature_soil")
    for i in range(n_bhe):
        _se(pv, "process_variable", f"temperature_BHE{i+1}")
    _bhe_xml(proc, cfg, n_bhe)  # BHE-Konfiguration einfügen

    # ── 3. Media (Materialeigenschaften) ──────────────────────────────────────
    # MaterialID 0 = unterste Schicht, ID 1 = nächste usw. (unten → oben)
    # dann MaterialID len(layers)..len(layers)+n_bhe-1 für die BHE-Linien
    media      = _se(root, "media")
    layers_bot = list(reversed(cfg["layers"]))  # Reihenfolge unten→oben
    bh_mat     = layers_bot[len(layers_bot) // 2]  # mittlere Schicht = repräsentativ für BHE
    mats       = layers_bot + [bh_mat] * n_bhe
    for mid, soil in enumerate(mats):
        med = _se(media, "medium", id=mid)
        phs = _se(med, "phases")
        # Phase 1: Porenfluid (Wasser)
        ph = _se(phs, "phase"); _se(ph, "type", "AqueousLiquid")
        pp = _se(ph, "properties")
        _const_prop(pp, "density",                fluid["rho_ref_kg_m3"])
        _const_prop(pp, "viscosity",              fluid["viscosity_Pa_s"])
        _const_prop(pp, "specific_heat_capacity", fluid["cp_J_kgK"])
        _const_prop(pp, "thermal_conductivity",   fluid["lambda_W_mK"])
        _const_prop(pp, "phase_velocity",         "0 0 0")  # kein Grundwasserfluss
        # Phase 2: Feststoff (Gestein)
        ph = _se(phs, "phase"); _se(ph, "type", "Solid")
        pp = _se(ph, "properties")
        _const_prop(pp, "density",                soil["rho_s_kg_m3"])
        _const_prop(pp, "specific_heat_capacity", soil["cp_s_J_kgK"])
        _const_prop(pp, "thermal_conductivity",   soil["lambda_s_W_mK"])
        # Medium-Eigenschaften (Kombination beider Phasen)
        props   = _se(med, "properties")
        lam_eff = (soil["porosity"] * fluid["lambda_W_mK"] +
                   (1.0 - soil["porosity"]) * soil["lambda_s_W_mK"])
        _const_prop(props, "porosity",             soil["porosity"])
        _const_prop(props, "permeability",         soil["permeability_m2"])
        _const_prop(props, "thermal_conductivity", lam_eff)
        _const_prop(props, "storage",              0.0)  # keine Druckspeicherung

    # ── 4. Time Loop ──────────────────────────────────────────────────────────
    tl    = _se(root, "time_loop")
    procs = _se(tl, "processes")
    pref  = _se(procs, "process", ref="HeatTransportBHE")
    _se(pref, "nonlinear_solver", "basic_picard")
    conv  = _se(pref, "convergence_criterion")
    _se(conv, "type", "DeltaX"); _se(conv, "norm_type", "NORM2")
    _se(conv, "reltol", sol.get("rel_tol_T", 1e-4))
    td = _se(pref, "time_discretization"); _se(td, "type", "BackwardEuler")
    ts = _se(pref, "time_stepping")
    _se(ts, "type", "FixedTimeStepping")
    _se(ts, "t_initial", 0); _se(ts, "t_end", t_end)
    pair = _se(_se(ts, "timesteps"), "pair")
    _se(pair, "repeat", n_steps); _se(pair, "delta_t", cfg["time"]["dt_seconds"])
    # Ausgabe: .vtu-Datei nach jedem Zeitschritt
    out   = _se(tl, "output")
    _se(out, "type", "VTK"); _se(out, "prefix", prefix)
    out_v = _se(out, "variables")
    _se(out_v, "variable", "temperature_soil")  # Temperatur im Erdreich
    for i in range(n_bhe):
        _se(out_v, "variable", f"temperature_BHE{i+1}")  # Fluidtemperatur je Sonde
    pair2 = _se(_se(out, "timesteps"), "pair")
    _se(pair2, "repeat", n_steps)
    _se(pair2, "each_steps", cfg["time"]["output_every_n_steps"])

    # ── 5. Parameters ─────────────────────────────────────────────────────────
    params = _se(root, "parameters")
    # T0: Anfangstemperatur Erdreich [K]
    p_T0 = _se(params, "parameter")
    _se(p_T0, "name", "T0"); _se(p_T0, "type", "Constant")
    _se(p_T0, "value", init["T_K"])
    # T0_BHE: Anfangstemperatur BHE-Fluid (4 Werte für 1U: Vorlauf, Rücklauf, 2× Grout)
    p_T0b = _se(params, "parameter")
    _se(p_T0b, "name", "T0_BHE"); _se(p_T0b, "type", "Constant")
    T0 = init["T_K"]
    _se(p_T0b, "value", f"{T0} {T0} {T0} {T0}")
    # unity: Skalierungsfaktor = 1 (wird von bhe_power benötigt)
    p_unity = _se(params, "parameter")
    _se(p_unity, "name", "unity"); _se(p_unity, "type", "Constant")
    _se(p_unity, "value", "1")
    # bhe_power: zeitabhängige Sonden-Leistung [W/Sonde] via CurveScaled
    # CurveScaled(t) = power_curve(t) × unity = power_curve(t) × 1 = power_curve(t)
    p_pow = _se(params, "parameter")
    _se(p_pow, "name", "bhe_power"); _se(p_pow, "type", "CurveScaled")
    _se(p_pow, "curve", "power_curve"); _se(p_pow, "parameter", "unity")

    # ── 6. Process Variables (Anfangs- und Randbedingungen) ───────────────────
    pvars = _se(root, "process_variables")
    # Erdreichtemperatur
    pvs = _se(pvars, "process_variable")
    _se(pvs, "name", "temperature_soil")
    _se(pvs, "components", 1); _se(pvs, "order", 1)
    _se(pvs, "initial_condition", "T0")  # Anfangswert: überall T0
    bcs = _se(pvs, "boundary_conditions")
    # Dirichlet-Randbedingung: T = T0 an Ober- und Unterseite (fixe Temperatur)
    for face in ("top", "bottom"):
        bc = _se(bcs, "boundary_condition")
        _se(bc, "mesh",      Path(mesh_files[face]).stem)
        _se(bc, "type",      "Dirichlet")
        _se(bc, "parameter", "T0")
    # BHE-Fluidtemperaturen (für jede Sonde)
    for i in range(n_bhe):
        pvb = _se(pvars, "process_variable")
        _se(pvb, "name", f"temperature_BHE{i+1}")
        _se(pvb, "components", 4)  # 1U-Sonde: 4 Freiheitsgrade
        _se(pvb, "order", 1)
        _se(pvb, "initial_condition", "T0_BHE")
        _se(pvb, "boundary_conditions")   # leer — Steuerung über BHE-Block

    # ── 7. Nonlinear Solver (Picard-Iteration) ────────────────────────────────
    nls = _se(root, "nonlinear_solvers")
    n   = _se(nls, "nonlinear_solver")
    _se(n, "name", "basic_picard"); _se(n, "type", "Picard")
    _se(n, "max_iter", sol["nonlinear_iter"])
    _se(n, "linear_solver", "general_linear_solver")

    # ── 8. Linear Solver (BiCGSTAB) ───────────────────────────────────────────
    # BiCGSTAB = Biconjugate Gradient Stabilized: iterativer Löser für Ax=b
    # ILUT = unvollständige LU-Zerlegung als Vorkonditionierer (beschleunigt Konvergenz)
    lss = _se(root, "linear_solvers")
    ls  = _se(lss, "linear_solver")
    _se(ls, "name", "general_linear_solver")
    eig = _se(ls, "eigen")
    _se(eig, "solver_type",        "BiCGSTAB")
    _se(eig, "precon_type",        "ILUT")
    _se(eig, "max_iteration_step", sol["linear_iter"])
    _se(eig, "error_tolerance",    sol["linear_tol"])
    _se(eig, "scaling",            "true")

    # ── 9. Curves (Leistungszeitreihe) ────────────────────────────────────────
    # OGS-Kurven sind Tabellen (Zeit [s], Wert [W]) — hier die Sondenleistung.
    # "power_curve" wird vom CurveScaled-Parameter "bhe_power" referenziert.
    t, v   = power_curve[1]
    curves = _se(root, "curves")
    c      = _se(curves, "curve")
    _se(c, "name",   "power_curve")
    _se(c, "coords", " ".join(f"{x:.6e}" for x in t))   # Zeiten [s]
    _se(c, "values", " ".join(f"{x:.6e}" for x in v))   # Leistungen [W/Sonde]

    # XML-Datei einrücken und speichern
    _indent(root)
    prj_path = out_dir / f"{prefix}.prj"
    ET.ElementTree(root).write(prj_path, encoding="ISO-8859-1", xml_declaration=True)
    return prj_path


# ══════════════════════════════════════════════════════════════════════════════
#  SCHRITT 5: OGS STARTEN
# ══════════════════════════════════════════════════════════════════════════════
#
# WAS MACHT OGS?
# ───────────────
# OpenGeoSys (OGS) ist ein quelloffener Finite-Elemente-Simulator, der an der
# Universität Leipzig / Helmholtz-Zentrum UFZ entwickelt wird.
# Er löst gekoppelte Gleichungen für Wärme-, Stoff- und Strömungstransport
# im Untergrund. In unserem Fall: HEAT_TRANSPORT_BHE = Wärmetransport mit
# Erdwärmesonden (Borehole Heat Exchangers).
#
# WAS BEDEUTET subprocess?
# ─────────────────────────
# Python kann andere Programme über "subprocess" starten — ähnlich wie man
# ein Programm in der Kommandozeile aufruft. OGS ist ein separates Programm
# (ogs.exe), das Python hier für uns startet:
#   ogs.exe  projektdatei.prj  -o  ausgabeordner/


def _find_ogs() -> list[str] | None:
    """
    Sucht die ogs-Binärdatei auf dem System.

    Das Problem: "ogs" ist in der conda-Umgebung installiert, aber nicht immer
    im System-PATH. Je nach Aufruf (Terminal/VS Code) findet Python ogs.exe
    über verschiedene Wege.

    Drei Suchmethoden:
    1. System-PATH (klappt wenn conda-Umgebung aktiviert ist)
    2. Nachbar-Verzeichnis zum Python-Interpreter
       conda: python.exe liegt in <env>/,       ogs.exe in <env>/Scripts/
       venv:  python.exe liegt in <env>/Scripts/, ogs.exe auch dort
    3. python -m ogs (falls ogs als Python-Modul mit __main__ installiert)

    Gibt zurück:
        Liste mit Befehlskomponenten, z. B. ["C:/miniconda/envs/ghe23/Scripts/ogs.exe"]
        oder [sys.executable, "-m", "ogs"]
        oder None wenn nicht gefunden
    """
    # Suchmethode 1: System-PATH
    ogs = shutil.which("ogs") or shutil.which("ogs.exe")
    if ogs:
        return [ogs]

    # Suchmethode 2: Nachbarverzeichnisse des aktuellen Python-Interpreters
    py = Path(sys.executable)  # z. B. C:/miniconda/envs/ghe23/python.exe
    for candidate in [
        py.parent / "ogs.exe",
        py.parent / "ogs",
        py.parent / "Scripts" / "ogs.exe",
        py.parent / "Scripts" / "ogs",
    ]:
        if candidate.is_file():
            return [str(candidate)]

    # Suchmethode 3: python -m ogs
    try:
        r = subprocess.run([sys.executable, "-m", "ogs", "--version"],
                           capture_output=True, timeout=10)
        if r.returncode == 0:
            return [sys.executable, "-m", "ogs"]
    except Exception:
        pass

    return None  # ogs nicht gefunden


def run_ogs(prj_path: Path) -> int:
    """
    Startet OGS mit der erzeugten Projektdatei.

    Gibt den Exit-Code zurück: 0 = Erfolg, 1 = Fehler.
    """
    cmd_prefix = _find_ogs()
    if cmd_prefix is None:
        print(
            "FEHLER: 'ogs' nicht gefunden.\n"
            f"  sys.executable = {sys.executable}\n"
            "  Lösung 1: Skript aus Terminal mit aktivierter Umgebung starten:\n"
            "            conda activate ghe23\n"
            "  Lösung 2: ogs.exe manuell aufrufen:\n"
            f"            ogs {prj_path} -o {prj_path.parent}",
            file=sys.stderr,
        )
        return 1
    return subprocess.call(cmd_prefix + [str(prj_path), "-o", str(prj_path.parent)])


# ══════════════════════════════════════════════════════════════════════════════
#  HAUPTFUNKTION  —  Ablaufsteuerung
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """
    Orchestriert den gesamten Ablauf in der richtigen Reihenfolge:

      1. Kommandozeilenargumente einlesen (--scenario, --no-mesh, --no-run)
      2. Ausgabeordner anlegen
      3. CONFIG aus CSVs und Einstellungen aufbauen (+ Energiebilanz drucken)
      4. [optional] Mesh erzeugen (gmsh) + konvertieren (msh2vtu)
      5. OGS-Projektdatei schreiben (.prj)
      6. [optional] OGS starten

    KOMMANDOZEILENARGUMENTE:
    ─────────────────────────
      --scenario before_renovation   Heizwärmebedarf vor Sanierung (Standard)
      --scenario after_renovation    Heizwärmebedarf nach Sanierung
      --no-mesh                      Mesh-Schritt überspringen
                                     (nützlich wenn das Mesh bereits existiert
                                      und nur die .prj neu erzeugt werden soll)
      --no-run                       Nur Mesh + .prj erzeugen, OGS nicht starten
                                     (nützlich zum schnellen Überprüfen des Setups)

    ERGEBNIS:
    ─────────
    Alle Dateien landen in  out/<szenario>/
      *.msh                       gmsh-Meshdatei (Zwischenergebnis)
      *_domain.vtu                Gesamtes Netz (Erdreich + Sondenlinien)
      *_physical_group_*.vtu      Ränder und Sondengruppen
      *.prj                       OGS-Projektdatei (XML)
      *.pvd + *_ts_*.vtu          OGS-Ergebnisse (nach dem Lauf)

    Die .pvd/.vtu-Ergebnisdateien können in ParaView visualisiert werden:
    → Datei laden → Temperaturfeld "temperature_soil" als Colormap anzeigen.
    """
    # Kommandozeilenparser einrichten
    ap = argparse.ArgumentParser(
        description="Saisonaler BTES-Solarspeicher (HEAT_TRANSPORT_BHE)")
    ap.add_argument(
        "--scenario", default=DEFAULT_SCENARIO,
        choices=["before_renovation", "after_renovation"],
        help="Gebäudeszenario (Standard: %(default)s)")
    ap.add_argument(
        "--no-mesh", action="store_true",
        help="Mesh-Generierung überspringen (VTU muss bereits existieren)")
    ap.add_argument(
        "--no-run", action="store_true",
        help="OGS-Lauf überspringen (nur Mesh + PRJ erzeugen)")
    args = ap.parse_args()

    # Ausgabeordner anlegen (wird erstellt falls er nicht existiert)
    out_dir = Path("out") / args.scenario
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[BTES Saisonalspeicher]  Szenario: {args.scenario}")

    # CONFIG aufbauen (liest CSVs, berechnet Energiebilanz, druckt Tabelle)
    cfg    = build_config(args.scenario, out_dir)
    prefix = cfg["output"]["prefix"]
    n_bhe  = len(_bhe_positions(cfg))
    msh_path = out_dir / f"{prefix}.msh"

    # Schritt 2+3: Mesh erzeugen und konvertieren
    if not args.no_mesh:
        print(f"[1/3] gmsh: {N_BHE_X}×{N_BHE_Y}-Feld, {n_bhe} BHEs, "
              f"Domain {cfg['domain']['size_x_m']:.0f}×{cfg['domain']['size_y_m']:.0f} m "
              f"× {BOREHOLE_DEPTH_BOTTOM_M + DOMAIN_DEPTH_BUFFER_M:.0f} m tief ...")
        build_mesh(cfg, out_dir)
        print(f"      → {msh_path}")
        print("[2/3] msh2vtu (.msh → .vtu Konvertierung) ...")
        msh2vtu(msh_path, out_dir, prefix, dim=[3, 1], reindex=True)
        # dim=[3, 1]: 3D-Volumenelemente (Erdreich) + 1D-Linienelemente (Sonden)

    # Dateinamen-Dictionary für die .vtu-Dateien (wird von build_prj benötigt)
    mesh_files = {
        "domain" : f"{prefix}_domain.vtu",
        "top"    : f"{prefix}_physical_group_top.vtu",
        "bottom" : f"{prefix}_physical_group_bottom.vtu",
        "lateral": f"{prefix}_physical_group_lateral.vtu",
    }
    for i in range(n_bhe):
        mesh_files[f"bhe_{i:02d}"] = f"{prefix}_physical_group_bhe_{i:02d}.vtu"

    # Schritt 4: OGS-Projektdatei erzeugen
    print(f"[3/3] OGS-Projektdatei ({n_bhe} BHEs, {N_YEARS} Jahre) ...")
    prj_path = build_prj(cfg, out_dir, mesh_files)
    print(f"      → {prj_path}")

    if args.no_run:
        print("--no-run gesetzt: OGS nicht gestartet.")
        return 0

    # Schritt 5: OGS starten
    print(f"\nStarte OGS ({N_YEARS * 12} Monate, Δt = {DT_SECONDS/86400:.0f} Tage) ...")
    print(f"Ergebnisse werden in {out_dir}/ gespeichert.")
    print("Visualisierung: ParaView → Datei öffnen → *.pvd → temperature_soil anzeigen\n")
    return run_ogs(prj_path)


# ══════════════════════════════════════════════════════════════════════════════
# Einstiegspunkt: Wird nur ausgeführt, wenn das Skript direkt gestartet wird
# (nicht wenn es von einem anderen Skript importiert wird).
# sys.exit() gibt den Rückgabewert (0=OK, 1=Fehler) an das Betriebssystem.
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    sys.exit(main())

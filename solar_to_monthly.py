#!/usr/bin/env python3
"""
Berechnet aus einem typisierten monatlichen Kollektor­ertrag und einem
monatlichen Wärmebedarf eine 12-Werte-Liste `monthly_power_W`, die
direkt in den `CONFIG["cycles"]["monthly_power_W"]`-Block der BTES-
oder ATES-Skripte eingesetzt werden kann.

Konzept (grobes Substitut für eine ausführliche Solar­ertrag­rechnung)
---------------------------------------------------------------------
Eingang  : - typisierter monatlicher Spezifikertrag q_sol(m)
             [kWh/m²/Monat]  — hier als Konstanten­profil hinterlegt
             (Mitteleuropa, Vakuum-Röhrenkollektor, β ≈ 40°)
           - Kollektor­fläche A_koll [m²]
           - Wärme­bedarf Q_bed(m) [kWh/Monat]
Bilanz   : ΔQ(m) = A_koll · q_sol(m) − Q_bed(m)   [kWh/Monat]
           + Überschuss → in den Speicher laden
           − Defizit    → aus dem Speicher fördern
Leistung : P(m) = ΔQ(m) · 3.6·10⁶ / ( d(m) · 86400 )   [W]

Hinweis. Die hier hinterlegten q_sol-Werte sind **typische Größen­ordnungen**
für eine südausgerichtete Vakuum-Röhrenkollektoranlage in
Mitteleuropa und ersetzen eine ausführliche standort­scharfe
Auslegung. Für eine reale Anlage müsste man die Strahlungs­daten am
konkreten Standort verwenden.

Verwendung
----------
    python solar_to_monthly.py
        → Demo mit Default-Bedarf, druckt die monthly_power_W-Liste.

Als Modul:
    from solar_to_monthly import monthly_power_W
    P = monthly_power_W(A_koll_m2=30.0, demand_kWh_per_month=[...12 Werte...])
    CONFIG["cycles"]["monthly_power_W"] = P
"""
from __future__ import annotations

import sys


# ----------------------------------------------------------------------
# Typisches Monats­profil: Vakuum-Röhrenkollektor, Mitteleuropa,
# Anstellwinkel ~ 40°, Süd-Ausrichtung. Werte in kWh/m²/Monat.
# Größen­ordnung kalibriert auf ca. 940 kWh/m²/a.
# ----------------------------------------------------------------------
SOLAR_YIELD_kWh_per_m2_month = [
    31.5,   # Jan
    48.1,   # Feb
    74.8,   # Mrz
   103.0,   # Apr
   117.7,   # Mai
   116.2,   # Jun
   122.2,   # Jul
   115.4,   # Aug
    88.7,   # Sep
    65.6,   # Okt
    33.0,   # Nov
    23.3,   # Dez
]

DAYS_PER_MONTH = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
MONTHS         = ["Jan", "Feb", "Mrz", "Apr", "Mai", "Jun",
                  "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"]

# Typischer normierter Heiz­bedarfs­verlauf (Anteil am Jahres­bedarf),
# Mitteleuropa, Σ = 1. Sommer minimal, Winter maximal.
DEMAND_PROFILE_FRAC = [0.175, 0.145, 0.120, 0.080, 0.045, 0.020,
                       0.010, 0.015, 0.040, 0.095, 0.125, 0.130]


def solar_monthly_yield() -> list[float]:
    """Liefert die 12 Monats­erträge q_sol [kWh/m²/Monat]."""
    return list(SOLAR_YIELD_kWh_per_m2_month)


def monthly_power_W(A_koll_m2: float,
                    demand_kWh_per_month) -> list[float]:
    """
    Berechnet die 12-Werte-Liste `monthly_power_W` [W].

    A_koll_m2            – Kollektor­fläche [m²]
    demand_kWh_per_month – Wärme­bedarf je Monat (Liste mit 12 Werten) [kWh]
    """
    if len(demand_kWh_per_month) != 12:
        raise ValueError("demand_kWh_per_month muss exakt 12 Werte enthalten.")
    q_sol = solar_monthly_yield()
    P = []
    for m in range(12):
        delta_Q_kWh = A_koll_m2 * q_sol[m] - float(demand_kWh_per_month[m])
        seconds = DAYS_PER_MONTH[m] * 86400.0
        P.append(delta_Q_kWh * 3.6e6 / seconds)
    return P


def annual_balance_kWh(A_koll_m2: float, demand_kWh_per_month) -> float:
    """Σ_m (A · q_sol − Q_bed)  — Vorzeichen-Indikator für Auslegung."""
    return A_koll_m2 * sum(SOLAR_YIELD_kWh_per_m2_month) - float(sum(demand_kWh_per_month))


def sizing_A_koll_for_balance(demand_kWh_per_month) -> float:
    """Kollektor­fläche [m²] so, dass Jahres­bilanz exakt 0 ergibt."""
    return float(sum(demand_kWh_per_month)) / sum(SOLAR_YIELD_kWh_per_m2_month)


def demand_from_annual(Q_annual_kWh: float,
                       profile: list[float] = DEMAND_PROFILE_FRAC) -> list[float]:
    """Skaliert das normierte Lastprofil (Σ = 1) auf den Jahres­bedarf."""
    s = sum(profile)
    return [Q_annual_kWh * (p / s) for p in profile]


# Demo
A_KOLL_DEMO   = 30.0           # m²
Q_ANNUAL_DEMO = 25_000.0       # kWh/a


def _print_table(q_sol, demand, P):
    print(f"{'Monat':<5} {'q_sol':>10} {'A·q_sol':>10} {'Q_bed':>10} {'ΔQ':>10} {'P':>12}")
    print(f"{'':<5} {'[kWh/m²]':>10} {'[kWh]':>10} {'[kWh]':>10} {'[kWh]':>10} {'[W]':>12}")
    for m in range(12):
        print(f"{MONTHS[m]:<5} {q_sol[m]:>10.1f} {q_sol[m]*A_KOLL_DEMO:>10.0f}"
              f" {demand[m]:>10.0f} {q_sol[m]*A_KOLL_DEMO - demand[m]:>10.0f}"
              f" {P[m]:>12.1f}")


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    print(f"Eingaben:  A_koll = {A_KOLL_DEMO} m²,  Q_jährlich = {Q_ANNUAL_DEMO:.0f} kWh")
    print(f"q_sol-Profil:  Vakuum-Röhrenkollektor, Mitteleuropa, β ≈ 40°")
    print(f"               Jahressumme = {sum(SOLAR_YIELD_kWh_per_m2_month):.0f} kWh/m²/a")
    print()
    q_sol  = solar_monthly_yield()
    demand = demand_from_annual(Q_ANNUAL_DEMO)
    P      = monthly_power_W(A_KOLL_DEMO, demand)
    _print_table(q_sol, demand, P)

    bal = annual_balance_kWh(A_KOLL_DEMO, demand)
    A_balanced = sizing_A_koll_for_balance(demand)
    print()
    print(f"Jahres­bilanz Σ ΔQ:  {bal:+.0f} kWh"
          f"   →  {'Überschuss (Speicher überdimensioniert)' if bal > 0 else 'Defizit'}")
    print(f"Empfohlene Fläche für  Σ ΔQ = 0:  A_koll ≈ {A_balanced:.1f} m²")

    print()
    print("# In OGS-CONFIG einsetzen:")
    print("CONFIG['cycles']['monthly_power_W'] = [")
    for p in P:
        print(f"    {p:.1f},")
    print("]")

# geothermie4-btes-ates

OpenGeoSys-Übungen zu **Borehole** und **Aquifer Thermal Energy Storage**
(BTES, ATES) — Vorlesung *Geothermie 4*.

📖 **Komplettes Tutorial:** siehe [`HANDBUCH.pdf`](HANDBUCH.pdf)
(Theorie, Übungen, Aufgaben, Plot-Interpretation, ~22 Seiten).

## Schnellstart

Voraussetzungen: **Python 3.10–3.12**, Windows / Linux / macOS.

```bash
git clone https://github.com/karli-a11y/geothermie4-btes-ates.git
cd geothermie4-btes-ates
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

`ogs.exe` muss im `PATH` auffindbar sein (wird vom `ogs`-Wheel
automatisch in `Scripts/` installiert — Verzeichnis ggf. zu `PATH`
hinzufügen). Windows-Hinweise zur Microsoft-Store-Python-Installation
im Handbuch, Abschnitt 5.

## Übung starten

```bash
cd btes/ex1_2d            # oder: ates/ex1_2d, btes/ex2_3d, ates/ex2_3d
python btes_radial_2d.py  # entsprechend: ates_radial_2d.py / btes_3d.py / ates_3d.py
python plot_results.py
```

Parameter werden im `CONFIG`-Block am Anfang des Sim-Skripts editiert
(Material, Betrieb, Zyklen, Geometrie, Mesh). Ergebnisbilder landen in
`figures/`.

## Repo-Inhalt

| Pfad                          | Inhalt                                          |
|-------------------------------|-------------------------------------------------|
| `HANDBUCH.{md,pdf}`           | Tutorial (Theorie + Übungen + Aufgaben)         |
| `btes/`, `ates/`              | Übungs-Ordner mit Sim- und Plot-Skripten        |
| `solar_to_monthly.py`         | Solar­ertrag → Monatsprofil-Helfer (typisiert)  |
| `formulas/`                   | Gerenderte LaTeX-Formelgrafiken (im Handbuch)   |
| `figures_illustrations/`      | Schemazeichnungen (im Handbuch)                 |
| `requirements.txt`            | Python-Abhängigkeiten                            |

# geothermie4-btes-ates

OpenGeoSys-Übungen zu **Borehole** und **Aquifer Thermal Energy Storage**
(BTES, ATES) — Vorlesung *Geothermie 4*.

📖 **Übungsskript:** siehe [`UEBUNGSSKRIPT.pdf`](UEBUNGSSKRIPT.pdf)
(Theorie, Übungen, Aufgaben, Plot-Interpretation, ~30 Seiten).

> 💬 **Komme nicht weiter?** Wenn die Einrichtung trotz dieser Anleitung
> hakt — besonders bei den im [Troubleshooting](#7-häufige-probleme--troubleshooting)
> genannten Fehlern — meldet euch **gerne jederzeit** bei mir:
> **holler@geo.tu-darmstadt.de**. Lieber einmal kurz fragen, als eine Stunde
> mit einem PATH-Problem zu kämpfen.

---

## Inhalt

1. [Überblick: Was wird installiert?](#1-überblick-was-wird-installiert)
2. [PowerShell öffnen (Windows)](#2-powershell-öffnen-windows)
3. [Software installieren (mit Download-Links)](#3-software-installieren-mit-download-links)
4. [Python-Pakete installieren — Variante wählen](#4-python-pakete-installieren--variante-wählen)
5. [Repository herunterladen](#5-repository-herunterladen)
6. [Eine Übung rechnen & prüfen](#6-eine-übung-rechnen--prüfen)
7. [Häufige Probleme — Troubleshooting](#7-häufige-probleme--troubleshooting)
8. [Erweiterte Varianten & Repo-Inhalt](#8-erweiterte-varianten--repo-inhalt)

---

## Schnellstart (für Erfahrene)

Voraussetzungen: **Python 3.10–3.12** (nicht 3.13+!), Windows / Linux / macOS.
Wer wenig Erfahrung mit Python/Kommandozeile hat, folgt besser der
ausführlichen Anleitung ab [Abschnitt 1](#1-überblick-was-wird-installiert).

```bash
git clone https://github.com/karli-a11y/geothermie4-btes-ates.git
cd geothermie4-btes-ates

# empfohlen: isolierte Umgebung (vermeidet PATH-Probleme mit ogs.exe)
python -m venv .venv
.\.venv\Scripts\Activate.ps1        # Windows PowerShell
# (Linux/macOS: source .venv/bin/activate)

python -m pip install --upgrade pip
python -m pip install -r requirements.txt

cd btes/ex1_2d && python btes_radial_2d.py   # erste Übung (~30 s)
```

---

## 1. Überblick: Was wird installiert?

| Software | Wozu | Pflicht? |
|----------|------|----------|
| **Python 3.10 – 3.12** | Programmiersprache, in der die Skripte laufen | ✅ ja |
| **Git** | Repository herunterladen & Updates ziehen | ✅ ja (Alternative: ZIP) |
| **Python-Pakete** (`ogs`, `ogstools`, `gmsh`, …) | Simulation, Mesh, Plots — per `pip`/`conda` | ✅ ja |
| **ParaView** | 3D-Visualisierung der Ergebnisse (`.vtu`/`.pvd`) | ⭕ optional, empfohlen |
| **VS Code** | Komfortabler Editor für die `CONFIG`-Blöcke | ⭕ optional |

> ⚠️ **Python-Version ist kritisch.** Für `ogs` und `ogstools` gibt es nur
> fertige Pakete (mit den OGS-Programmen darin) für **Python 3.10, 3.11 oder
> 3.12**. Unter **3.13 / 3.14** schlägt der Lauf später mit kryptischen
> Fehlern fehl (z. B. `OGS_BIN_DIR does not exist`, Crash in
> `NodeReordering`). Bitte unbedingt eine Version aus 3.10–3.12 verwenden.

---

## 2. PowerShell öffnen (Windows)

Alle Befehle hier werden in **PowerShell** eingetippt (Windows'
Kommandozeile). So öffnest du sie:

**Variante 1 (am einfachsten):**
1. `Windows-Taste` drücken.
2. `powershell` tippen.
3. In der Trefferliste **„Windows PowerShell"** anklicken.

**Variante 2 (direkt im Projektordner):**
1. Im Datei-Explorer in den gewünschten Ordner navigieren.
2. Bei gedrückter `Shift`-Taste mit der **rechten Maustaste** in den leeren
   Bereich klicken.
3. **„PowerShell-Fenster hier öffnen"** wählen.

Ein PowerShell-Fenster zeigt eine Zeile wie `PS C:\Users\DeinName>`. Hinter
dem `>` werden Befehle eingegeben und mit `Enter` ausgeführt.

> **Mini-Spickzettel:** `cd Ordnername` wechselt in einen Ordner,
> `cd ..` geht eine Ebene hoch, `ls` (oder `dir`) listet den Inhalt.

### Erste Kontrolle: Ist Python da?

```powershell
python --version
```

- **`Python 3.10.x` / `3.11.x` / `3.12.x`** → super, weiter zu
  [Abschnitt 4](#4-python-pakete-installieren--variante-wählen).
- **`3.13` / `3.14`**, eine Fehlermeldung, oder der **Microsoft Store** öffnet
  sich → zuerst [Abschnitt 3](#3-software-installieren-mit-download-links).

---

## 3. Software installieren (mit Download-Links)

### 3.1 Python (Pflicht)

- **Download:** <https://www.python.org/downloads/release/python-3129/>
  (Python **3.12** von python.org — **nicht** aus dem Microsoft Store).
  Unter *„Files"* den **„Windows installer (64-bit)"** wählen.
- Im Installer **unten das Häkchen bei „Add python.exe to PATH" setzen**,
  dann **„Install Now"**.

> ⚠️ **Microsoft-Store-Python vermeiden.** Öffnet sich beim Tippen von
> `python` der Store, ist meist die Store-Variante aktiv. Diese legt Programme
> an einem versteckten Ort ab
> (`…\AppData\Local\Packages\PythonSoftwareFoundation.Python…\`), der nicht
> im `PATH` liegt — genau das löst später *„ogs nicht gefunden"* aus.

Danach **PowerShell neu öffnen** und prüfen: `python --version`.

### 3.2 Git (Pflicht für die „clone"-Variante)

- **Download:** <https://git-scm.com/download/win>
- Mit Standard-Einstellungen durchklicken, danach `git --version` prüfen.

*(Ohne Git geht es auch — siehe [Abschnitt 5](#5-repository-herunterladen),
ZIP-Variante.)*

### 3.3 conda — nur für Variante C (optional)

Wer lieber mit **conda** arbeitet, braucht **Miniconda**:

- **Download:** <https://www.anaconda.com/download/success> → *„Miniconda
  Installers"* → Windows 64-bit.
- Danach im Startmenü **„Anaconda Prompt"** verwenden (statt normaler
  PowerShell).

### 3.4 ParaView (optional, empfohlen)

3D-Visualisierung der `.pvd`/`.vtu`-Ergebnisse:
<https://www.paraview.org/download/>

### 3.5 VS Code (optional)

Editor für die `CONFIG`-Blöcke: <https://code.visualstudio.com/>

---

## 4. Python-Pakete installieren — Variante wählen

Hier liegt der Schlüssel gegen die *„ogs nicht im PATH"*-Probleme. **Wähle
genau eine** Variante, passend zu deinem Fall:

| Variante | Für wen? |
|----------|----------|
| **A — venv** (empfohlen) | Du hast Python von python.org installiert und willst es sauber pro Projekt halten. |
| **B — direkt mit pip** | Schnell & ohne Umgebung — nur wenn `python --version` bereits 3.10–3.12 zeigt und du es einfach willst. |
| **C — conda** | Du nutzt ohnehin Anaconda/Miniconda. |

> **Warum eine Umgebung (A/C)?** Eine virtuelle bzw. conda-Umgebung ist ein
> abgeschotteter „Kasten" nur für dieses Projekt. Solange er **aktiviert**
> ist, liegt `ogs.exe` automatisch im Suchpfad — der Befehl `ogs` wird
> gefunden, **ohne** dass du den `PATH` von Hand bearbeiten musst. Genau das
> vermeidet den häufigsten Fehler.

---

### Variante A — virtuelle Umgebung (venv) + pip  *(empfohlen)*

In **PowerShell**, im Ordner, in dem das Projekt liegen soll (ggf. erst
[Repo klonen](#5-repository-herunterladen) und hineinwechseln):

```powershell
# 1) Virtuelle Umgebung anlegen (Ordner ".venv" entsteht)
python -m venv .venv

# 2) Aktivieren  -> die Eingabezeile beginnt danach mit "(.venv)"
.\.venv\Scripts\Activate.ps1

# 3) Pakete installieren
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

> **Falls Schritt 2 mit „…Ausführung von Skripts ist deaktiviert" abbricht**,
> einmalig erlauben und Schritt 2 wiederholen:
> ```powershell
> Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
> ```

**Wichtig:** Die Umgebung in **jeder neuen PowerShell-Sitzung** erneut
aktivieren (`.\.venv\Scripts\Activate.ps1`), bevor du rechnest — erkennbar am
`(.venv)` vorne in der Zeile.

*(Linux/macOS: `source .venv/bin/activate`.)*

---

### Variante B — direkt mit pip (ohne Umgebung)

Nur empfehlenswert, wenn `python --version` schon 3.10–3.12 zeigt:

```powershell
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

> ⚠️ Hier kann es passieren, dass `ogs.exe` zwar installiert, aber nicht im
> `PATH` ist (*„ogs nicht gefunden"*). Lösungen dazu in
> [Abschnitt 7](#7-häufige-probleme--troubleshooting). Wenn dich das trifft,
> ist Variante A der unkompliziertere Weg.

---

### Variante C — conda-Umgebung

In der **Anaconda Prompt** (aus [Abschnitt 3.3](#33-conda--nur-für-variante-c-optional)):

```powershell
# 1) Umgebung mit passender Python-Version
conda create -n geothermie python=3.12

# 2) Aktivieren  -> Eingabezeile beginnt mit "(geothermie)"
conda activate geothermie

# 3) Pakete (ogs & ogstools per pip — nicht in den Standard-conda-Kanälen)
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

**Wichtig:** In **jeder neuen Sitzung** zuerst `conda activate geothermie`
(erkennbar am `(geothermie)` vorne).

---

### Installation prüfen (alle Varianten)

Bei aktiver Umgebung bzw. nach der Installation:

```powershell
ogs --version
```

Erscheint `ogs version: 6.5.x`, ist alles bereit. Erscheint *„ogs … wird
nicht erkannt"* → [Abschnitt 7](#7-häufige-probleme--troubleshooting).

---

## 5. Repository herunterladen

**Mit Git (empfohlen):**

```powershell
git clone https://github.com/karli-a11y/geothermie4-btes-ates.git
cd geothermie4-btes-ates
```

**Ohne Git (ZIP):**
1. <https://github.com/karli-a11y/geothermie4-btes-ates> öffnen.
2. Grüner Button **„Code" → „Download ZIP"**, entpacken.
3. In PowerShell mit `cd` in den entpackten Ordner wechseln.

> **Reihenfolge-Tipp (Variante A):** erst klonen, dann
> `cd geothermie4-btes-ates`, dann `python -m venv .venv` — so liegt die
> Umgebung im Projektordner.

---

## 6. Eine Übung rechnen & prüfen

Empfohlener Einstieg: die schnelle **BTES-2D-Übung** (~30 s). Bei aktiver
Umgebung, im Projektordner:

```powershell
cd btes\ex1_2d
python btes_radial_2d.py
```

Das Skript erzeugt Netz + OGS-Projektdatei, startet OGS und legt am Ende
automatisch Plots in `btes\ex1_2d\figures\` ab. Fertig, wenn unten u. a.
`saved energy_balance.png` steht. Ausführlichere Auswertegrafiken:

```powershell
python plot_results.py
```

Andere Übungen analog (jeweils zuerst per `cd` in den Ordner wechseln):

| Ordner | Skript | System | Laufzeit (Default) |
|--------|--------|--------|--------------------|
| `btes\ex1_2d` | `btes_radial_2d.py` | BTES radial 2D | ~30 s |
| `btes\ex2_3d` | `btes_3d.py` | BTES Sondenfeld 3D | ~10–30 min |
| `ates\ex1_2d` | `ates_radial_2d.py` | ATES Single-Well 2D | ~1 min |
| `ates\ex2_3d` | `ates_3d.py` | ATES Single-Well 3D | ~15–60 min |

Praktische Optionen der Sim-Skripte:

```powershell
python btes_radial_2d.py --no-run     # nur Setup (Mesh + .prj), kein OGS-Lauf
python btes_radial_2d.py --no-mesh    # Mesh wiederverwenden, nur neu rechnen
python btes_radial_2d.py --no-plots   # rechnen, aber keine Auto-Plots
```

Parameter (Material, Zyklen, Geometrie, Mesh) stehen im `CONFIG`-Block am
Anfang des jeweiligen Sim-Skripts. Details im
[`UEBUNGSSKRIPT.pdf`](UEBUNGSSKRIPT.pdf).

---

## 7. Häufige Probleme — Troubleshooting

> 💬 Wenn einer dieser Fälle bei euch auftritt und ihr nicht weiterkommt:
> einfach eine kurze Mail an **holler@geo.tu-darmstadt.de** — am besten mit
> der genauen Fehlermeldung (Screenshot reicht) und eurer Python-Version.

### „`ogs` wird nicht … erkannt" / „ogs.exe nicht im PATH"

Häufigstes Problem. `ogs.exe` ist installiert, aber nicht im Suchpfad.

1. **Umgebung aktiviert?** Steht vorne `(.venv)` bzw. `(geothermie)`?
   Falls nicht → aktivieren (`.\.venv\Scripts\Activate.ps1` bzw.
   `conda activate geothermie`). In einer aktiven Umgebung verschwindet
   dieser Fehler meist von allein.
2. **Wo liegt ogs?**
   ```powershell
   python -c "import shutil; print(shutil.which('ogs'))"
   ```
   Kommt `None` → `ogs` ist nicht installiert
   ([Abschnitt 4](#4-python-pakete-installieren--variante-wählen) bei aktiver
   Umgebung wiederholen).
3. **Notlösung ohne PATH:** OGS direkt als Modul aufrufen —
   `python -m ogs --version`. Sauberer ist aber Variante A oder C.

### `OGS_BIN_DIR does not exist` / Crash in `NodeReordering` / `posix_spawn`

**Falsche Python-Version.** `ogs`/`ogstools` haben nur Wheels für
**Python 3.10–3.12**. → Umgebung mit Python 3.12 neu aufsetzen
(`py -3.12 -m venv .venv` bzw. `conda create -n geothermie python=3.12`).

### `ModuleNotFoundError: No module named 'ogstools'` (o. ä.)

Pakete sind nicht in der aktiven Umgebung installiert → Umgebung aktivieren,
[Abschnitt 4](#4-python-pakete-installieren--variante-wählen) wiederholen.

### `Activate.ps1 kann nicht geladen werden … Ausführung von Skripts ist deaktiviert`

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```
Danach erneut aktivieren.

### Beim Tippen von `python` öffnet sich der Microsoft Store

Python von **python.org** installieren ([Abschnitt 3.1](#31-python-pflicht),
Häkchen *„Add python.exe to PATH"*), PowerShell neu öffnen. Optional die
Store-Aliase abschalten: *Einstellungen → Apps → Erweiterte
App-Einstellungen → App-Ausführungsaliase* → „python.exe"/„python3.exe" aus.

### `FileNotFoundError: … .pvd nicht gefunden` beim Plot-Skript

`plot_results.py` braucht die Simulationsergebnisse. Erst das Sim-Skript
laufen lassen, dann das Plot-Skript.

*Weitere fachliche Fehlersuche: [`UEBUNGSSKRIPT.pdf`](UEBUNGSSKRIPT.pdf),
Kapitel 5 und 10.*

---

## 8. Erweiterte Varianten & Repo-Inhalt

**Erweiterte Sim-Variante:**

- `btes/ex2_3d/btes_3d_bhe.py` — BTES 3D mit OGS-Modul
  **`HEAT_TRANSPORT_BHE`** statt vereinfachter Volumen-Quelle. Sonden als
  1D-Linien­elemente eingebettet; BHE-Typ (1U/2U/CXA/CXC), U-Rohr-Geometrie,
  Refrigerant und Steuerung über `CONFIG["bhe"]` einstellbar. **Läuft
  end-to-end durch** (Default: 3×3-Feld, Typ 1U, ein Lade-/Entlade-Zyklus).
  Quantitative Validierung gegen die OGS-Benchmarks `BHE_1U`/`BHE_2U` steht
  noch aus — Ergebnisse qualitativ interpretieren.

**Repo-Inhalt:**

| Pfad                          | Inhalt                                          |
|-------------------------------|-------------------------------------------------|
| `UEBUNGSSKRIPT.{md,pdf}`      | Übungsskript (Theorie + Übungen + Aufgaben)     |
| `btes/`, `ates/`              | Übungs-Ordner mit Sim- und Plot-Skripten        |
| `formulas/`                   | Gerenderte LaTeX-Formelgrafiken (im Übungsskript) |
| `figures_illustrations/`      | Schemazeichnungen (im Übungsskript)             |
| `requirements.txt`            | Python-Abhängigkeiten                            |

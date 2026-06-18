# Installations- & Ausführungsanleitung (Schritt für Schritt)

Diese Anleitung führt von einem frischen Windows-Rechner bis zur ersten
fertig gerechneten Übung mit Plots. Sie richtet sich an Studierende **ohne**
Vorkenntnisse in Python/Kommandozeile. Wenn etwas nicht klappt, siehe
[Abschnitt 7 — Häufige Probleme](#7-häufige-probleme-troubleshooting).

> **Wichtigster Tipp vorab:** Die mit Abstand häufigste Fehlerquelle ist
> *„`ogs.exe` / ogstools nicht im PATH gefunden"*. Wer den Empfehlungen unten
> folgt (**virtuelle Umgebung** oder **conda-Umgebung**), umgeht dieses
> Problem von vornherein — siehe [Abschnitt 4](#4-python-pakete-installieren).

---

## Inhalt

1. [Überblick: Was wird installiert?](#1-überblick-was-wird-installiert)
2. [PowerShell öffnen](#2-powershell-öffnen)
3. [Software installieren (mit Download-Links)](#3-software-installieren-mit-download-links)
4. [Python-Pakete installieren — Weg A (venv) oder Weg B (conda)](#4-python-pakete-installieren)
5. [Repository herunterladen](#5-repository-herunterladen)
6. [Eine Übung rechnen & prüfen](#6-eine-übung-rechnen--prüfen)
7. [Häufige Probleme (Troubleshooting)](#7-häufige-probleme-troubleshooting)

---

## 1. Überblick: Was wird installiert?

| Software | Wozu | Pflicht? |
|----------|------|----------|
| **Python 3.10 – 3.12** | Programmiersprache, in der die Skripte laufen | ✅ ja |
| **Git** | Repository herunterladen & Updates ziehen | ✅ ja (Alternative: ZIP) |
| **Python-Pakete** (`ogs`, `ogstools`, `gmsh`, …) | Simulation, Mesh, Plots — werden per `pip`/`conda` installiert | ✅ ja |
| **ParaView** | 3D-Visualisierung der Ergebnisse (`.vtu`/`.pvd`) | ⭕ optional, empfohlen |
| **VS Code** | Komfortabler Editor für die `CONFIG`-Blöcke | ⭕ optional |

> ⚠️ **Python-Version ist kritisch.** Für `ogs` und `ogstools` gibt es nur
> fertige Pakete (mit den OGS-Programmen darin) für **Python 3.10, 3.11
> oder 3.12**. Unter **3.13 / 3.14** schlägt der Lauf später mit kryptischen
> Fehlern fehl (z. B. `OGS_BIN_DIR does not exist`, Crash in
> `NodeReordering`). Bitte unbedingt eine Version aus dem Bereich 3.10–3.12
> verwenden.

---

## 2. PowerShell öffnen

Alle Befehle in dieser Anleitung werden in **PowerShell** eingetippt (das
ist Windows' Kommandozeile). So öffnest du sie:

**Variante 1 (am einfachsten):**
1. `Windows-Taste` drücken.
2. `powershell` tippen.
3. In der Trefferliste **„Windows PowerShell"** anklicken.

**Variante 2 (im Projektordner):**
1. Im Datei-Explorer in den Ordner navigieren, in dem das Projekt liegen
   soll (z. B. `Dokumente`).
2. Bei gedrückter `Shift`-Taste mit der **rechten Maustaste** in den leeren
   Bereich des Ordners klicken.
3. **„PowerShell-Fenster hier öffnen"** wählen — PowerShell startet direkt im
   richtigen Verzeichnis.

Ein PowerShell-Fenster zeigt eine Zeile wie
`PS C:\Users\DeinName>`. Hinter diesem `>` werden die Befehle eingegeben und
mit `Enter` ausgeführt.

> **Tipp — Verzeichnis wechseln:** Mit `cd` ("change directory") wechselt man
> in einen Ordner, z. B. `cd Dokuments`. Mit `cd ..` geht man eine Ebene
> nach oben. Mit `ls` (oder `dir`) listet man den Inhalt des aktuellen
> Ordners auf.

### Erste Kontrolle: Ist Python da?

```powershell
python --version
```

- Erscheint **`Python 3.10.x` / `3.11.x` / `3.12.x`** → super, weiter zu
  [Abschnitt 4](#4-python-pakete-installieren) (Python ist installiert).
- Erscheint **`3.13` / `3.14`** oder **eine Fehlermeldung** bzw. öffnet sich
  der **Microsoft Store** → zuerst [Abschnitt 3](#3-software-installieren-mit-download-links)
  durcharbeiten.

---

## 3. Software installieren (mit Download-Links)

### 3.1 Python (Pflicht)

- **Download:** <https://www.python.org/downloads/release/python-3129/>
  (Python **3.12** — direkt von python.org, **nicht** aus dem Microsoft Store,
  siehe Hinweis unten).
  Auf der Seite unter *„Files"* den **„Windows installer (64-bit)"** wählen.
- Im Installer **ganz unten das Häkchen bei
  „Add python.exe to PATH" setzen**, dann **„Install Now"**.

> ⚠️ **Microsoft-Store-Python vermeiden.** Tippt man unter Windows einfach
> `python` und es öffnet sich der Store, ist meist die Store-Variante aktiv.
> Diese legt Programme an einem versteckten Ort ab
> (`…\AppData\Local\Packages\PythonSoftwareFoundation.Python…\`), der nicht
> im `PATH` liegt — genau das löst später *„ogs nicht gefunden"* aus. Am
> robustesten: Installer von **python.org** verwenden **und** mit einer
> virtuellen/conda-Umgebung arbeiten ([Abschnitt 4](#4-python-pakete-installieren)).

Nach der Installation **PowerShell schließen und neu öffnen** und prüfen:

```powershell
python --version
```

### 3.2 Git (Pflicht für Variante „clone")

- **Download:** <https://git-scm.com/download/win>
- Installer mit allen Standard-Einstellungen durchklicken.
- Danach PowerShell neu öffnen und prüfen: `git --version`.

*(Ohne Git geht es auch — siehe [Abschnitt 5](#5-repository-herunterladen),
ZIP-Variante.)*

### 3.3 conda — nur für Weg B (optional)

Wer lieber mit **conda** arbeitet (oder es schon installiert hat), braucht
**Miniconda**:

- **Download:** <https://www.anaconda.com/download/success> → Abschnitt
  *„Miniconda Installers"* → Windows 64-bit.
- Nach der Installation gibt es im Startmenü einen Eintrag
  **„Anaconda Prompt"** bzw. **„Anaconda PowerShell Prompt"** — diesen für
  Weg B verwenden (statt der normalen PowerShell).

### 3.4 ParaView (optional, empfohlen)

Zur 3D-Visualisierung der Ergebnis-Dateien (`.pvd`/`.vtu`):

- **Download:** <https://www.paraview.org/download/>

### 3.5 VS Code (optional)

Komfortabler Editor, um die `CONFIG`-Blöcke der Skripte zu bearbeiten:

- **Download:** <https://code.visualstudio.com/>

---

## 4. Python-Pakete installieren

Hier liegt der Schlüssel gegen die *„ogs nicht im PATH"*-Probleme. Es gibt
zwei Wege — **wähle genau einen**:

- **Weg A — venv ("normal", nur Python + pip):** empfohlen, wenn du oben
  Python von python.org installiert hast.
- **Weg B — conda:** empfohlen, wenn du ohnehin Anaconda/Miniconda nutzt.

> **Warum überhaupt eine Umgebung?** Eine virtuelle Umgebung (venv) bzw.
> conda-Umgebung ist ein abgeschotteter „Kasten" nur für dieses Projekt.
> Solange er **aktiviert** ist, liegt `ogs.exe` automatisch im Suchpfad —
> der Befehl `ogs` wird also gefunden, **ohne** dass du den `PATH` manuell
> bearbeiten musst. Genau das vermeidet den häufigsten Fehler.

---

### Weg A — virtuelle Umgebung (venv) + pip

In **PowerShell** (zuerst in den Ordner wechseln, in dem das Projekt liegen
soll, z. B. `cd Dokuments`):

```powershell
# 1) Virtuelle Umgebung anlegen (Ordner ".venv" wird erstellt)
python -m venv .venv

# 2) Umgebung aktivieren  -> die Eingabezeile beginnt danach mit "(.venv)"
.\.venv\Scripts\Activate.ps1

# 3) pip aktualisieren und Pakete installieren
python -m pip install --upgrade pip
python -m pip install ogs ogstools gmsh meshio numpy "pyvista>=0.45" matplotlib
```

> **Falls Schritt 2 mit „… kann nicht geladen werden, da die Ausführung von
> Skripts auf diesem System deaktiviert ist" abbricht**, einmalig erlauben:
> ```powershell
> Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
> ```
> Danach Schritt 2 wiederholen.

**Wichtig:** Diese Umgebung muss in **jeder neuen PowerShell-Sitzung** erneut
aktiviert werden, bevor du eine Übung rechnest:

```powershell
.\.venv\Scripts\Activate.ps1
```

(Erkennbar am `(.venv)` vorne in der Eingabezeile.)

Wenn du das Repo schon hast, kannst du statt der langen Paketliste auch
`python -m pip install -r requirements.txt` benutzen.

---

### Weg B — conda-Umgebung

In der **Anaconda Prompt** (bzw. Anaconda PowerShell Prompt) aus
[Abschnitt 3.3](#33-conda--nur-für-weg-b-optional):

```powershell
# 1) Umgebung mit passender Python-Version anlegen
conda create -n geothermie python=3.12

# 2) Umgebung aktivieren  -> Eingabezeile beginnt mit "(geothermie)"
conda activate geothermie

# 3) Pakete installieren (ogs & ogstools kommen per pip — sie liegen
#    nicht in den Standard-conda-Kanälen)
python -m pip install --upgrade pip
python -m pip install ogs ogstools gmsh meshio numpy "pyvista>=0.45" matplotlib
```

**Wichtig:** In **jeder neuen Sitzung** zuerst aktivieren:

```powershell
conda activate geothermie
```

(Erkennbar am `(geothermie)` vorne in der Eingabezeile.)

---

### Installation prüfen (beide Wege)

Bei **aktivierter** Umgebung:

```powershell
ogs --version
```

Erscheint eine Zeile wie `ogs version: 6.5.7`, ist alles bereit. Erscheint
*„ogs … wird nicht erkannt"*, ist die Umgebung vermutlich nicht aktiviert
(kein `(.venv)`/`(geothermie)` vorne) → siehe
[Abschnitt 7](#7-häufige-probleme-troubleshooting).

---

## 5. Repository herunterladen

**Mit Git (empfohlen):** in PowerShell, im gewünschten Zielordner:

```powershell
git clone https://github.com/karli-a11y/geothermie4-btes-ates.git
cd geothermie4-btes-ates
```

**Ohne Git (ZIP):**
1. <https://github.com/karli-a11y/geothermie4-btes-ates> öffnen.
2. Grüner Button **„Code" → „Download ZIP"**.
3. ZIP entpacken, dann in PowerShell mit `cd` in den entpackten Ordner
   wechseln.

> **Reihenfolge-Tipp:** Wer Weg A (venv) nutzt, legt die Umgebung am besten
> **innerhalb** des geklonten Projektordners an (also erst klonen, dann
> `cd geothermie4-btes-ates`, dann `python -m venv .venv`).

---

## 6. Eine Übung rechnen & prüfen

Empfohlener Einstieg: die schnelle **BTES-2D-Übung** (~30 s).

Bei **aktivierter** Umgebung, im Projektordner:

```powershell
cd btes\ex1_2d
python btes_radial_2d.py
```

Das Skript erzeugt das Netz, schreibt die OGS-Projektdatei, startet OGS und
legt am Ende automatisch Plots in `btes\ex1_2d\figures\` ab. Es ist fertig,
wenn unten u. a. `saved energy_balance.png` steht.

Zusätzliche, ausführlichere Auswertegrafiken erzeugt:

```powershell
python plot_results.py
```

Andere Übungen analog (jeweils zuerst per `cd` in den Ordner wechseln):

| Ordner | Skript | System |
|--------|--------|--------|
| `btes\ex1_2d` | `btes_radial_2d.py` | BTES radial 2D (schnell) |
| `btes\ex2_3d` | `btes_3d.py` | BTES Sondenfeld 3D |
| `ates\ex1_2d` | `ates_radial_2d.py` | ATES Single-Well 2D |
| `ates\ex2_3d` | `ates_3d.py` | ATES Single-Well 3D |

Praktische Optionen für die Sim-Skripte:

```powershell
python btes_radial_2d.py --no-run     # nur Setup (Mesh + .prj), kein OGS-Lauf
python btes_radial_2d.py --no-mesh    # Mesh wiederverwenden, nur neu rechnen
python btes_radial_2d.py --no-plots   # rechnen, aber keine Auto-Plots
```

Parameter (Material, Zyklen, Geometrie …) stehen im `CONFIG`-Block am Anfang
des jeweiligen Sim-Skripts und können mit einem Editor angepasst werden.
Details dazu im Übungsskript ([`UEBUNGSSKRIPT.pdf`](UEBUNGSSKRIPT.pdf)).

---

## 7. Häufige Probleme (Troubleshooting)

### „`ogs` wird nicht als Name eines Cmdlets … erkannt" / „ogs.exe nicht im PATH"

Das ist das häufigste Problem. Ursache: `ogs.exe` wurde installiert, liegt
aber nicht im Suchpfad (`PATH`).

1. **Umgebung aktiviert?** Steht vorne in der Eingabezeile `(.venv)` bzw.
   `(geothermie)`? Falls nicht → aktivieren:
   - venv: `.\.venv\Scripts\Activate.ps1`
   - conda: `conda activate geothermie`

   In einer **aktivierten** Umgebung liegt `ogs.exe` automatisch im PATH —
   damit verschwindet dieser Fehler in aller Regel.

2. **Wo liegt ogs überhaupt?** Prüfen mit:
   ```powershell
   python -c "import shutil; print(shutil.which('ogs'))"
   ```
   - Kommt ein Pfad heraus, der `\.venv\` bzw. den conda-Umgebungsnamen
     enthält → korrekt.
   - Kommt `None`, ist `ogs` nicht installiert → Installation aus
     [Abschnitt 4](#4-python-pakete-installieren) (bei aktivierter Umgebung)
     wiederholen.

3. **Notlösung ohne PATH** (falls du partout keine Umgebung nutzen willst):
   OGS lässt sich auch direkt als Python-Modul aufrufen:
   ```powershell
   python -m ogs --version
   ```
   Sauberer ist aber Weg A oder B.

### `OGS_BIN_DIR does not exist` / Crash in `NodeReordering` / `posix_spawn`

**Falsche Python-Version.** `ogs`/`ogstools` haben nur Wheels (inkl. der
mitgelieferten Programme) für **Python 3.10–3.12**. Unter 3.13/3.14 fehlt
z. B. die `NodeReordering`-Datei.
→ Umgebung mit Python 3.12 neu aufsetzen:
- venv: neue Umgebung mit einem 3.12-Python anlegen
  (`py -3.12 -m venv .venv`, falls der Launcher `py` vorhanden ist).
- conda: `conda create -n geothermie python=3.12`.

### `ModuleNotFoundError: No module named 'ogstools'` (o. ä.)

Pakete sind nicht (in der aktiven Umgebung) installiert. Umgebung
aktivieren und Installationsschritt aus [Abschnitt 4](#4-python-pakete-installieren)
wiederholen.

### `Activate.ps1 kann nicht geladen werden … Ausführung von Skripts ist deaktiviert`

Einmalig die Ausführungs-Richtlinie lockern, dann erneut aktivieren:
```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

### Beim Tippen von `python` öffnet sich der Microsoft Store

Die Store-Weiterleitung ist aktiv und/oder es ist kein „echtes" Python im
PATH. Python von **python.org** installieren ([Abschnitt 3.1](#31-python-pflicht),
Häkchen *„Add python.exe to PATH"*), PowerShell neu öffnen. Optional die
Store-Aliase deaktivieren unter
*Einstellungen → Apps → Erweiterte App-Einstellungen → App-Ausführungsaliase*
→ „python.exe"/„python3.exe" ausschalten.

### `FileNotFoundError: … .pvd nicht gefunden` beim Plot-Skript

`plot_results.py` braucht die Ergebnisse der Simulation. Erst das Sim-Skript
laufen lassen (`python btes_radial_2d.py`), dann das Plot-Skript.

---

Weitere fachliche Hinweise und eine ausführliche Fehlersuche-Tabelle stehen
im Übungsskript [`UEBUNGSSKRIPT.pdf`](UEBUNGSSKRIPT.pdf), Kapitel 5 und 10.

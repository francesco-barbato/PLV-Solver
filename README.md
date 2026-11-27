# SolutorePLV 📐🏛️

**SolutorePLV** è un'applicazione avanzata con interfaccia grafica (GUI) e motore di calcolo analitico/simbolico per la risoluzione, verifica e dimensionamento di strutture intelaiate piane (isostatiche e una volta iperstatiche) tramite il **Principio dei Lavori Virtuali (PLV)**, sviluppata specificamente per l'ambito della **Scienza delle Costruzioni** e dell'**Ingegneria Strutturale**.

Il programma combina la potenza del calcolo algebrico esatto di **SymPy** con la ricchezza visiva di **Matplotlib** e un'interfaccia utente moderna e reattiva in **PyQt6**, producendo relazioni di calcolo analitiche dettagliate passo dopo passo ed esportabili direttamente in **PDF**.

---

## 📑 Indice

- [🪟 Come accedere alla Versione Windows 11](#-come-accedere-alla-versione-windows-11)
- [Caratteristiche Principali](#-caratteristiche-principali)
- [Architettura del Progetto](#-architettura-del-progetto)
- [Prerequisiti e Installazione](#-prerequisiti-e-installazione)
- [Modalità d'Uso](#-modalità-duso)
  - [1. Avvio dell'Interfaccia Grafica](#1-avvio-dellinterfaccia-grafica-gui)
  - [2. Utilizzo come Libreria / Script Python](#2-utilizzo-come-libreria--script-python)
- [Formato dei Dati (JSON Schema)](#-formato-dei-dati-json-schema)
- [Flusso di Risoluzione (Workflow PLV)](#-flusso-di-risoluzione-workflow-plv)
- [Compilazione Standalone (PyInstaller)](#-compilazione-standalone-pyinstaller)
- [Esempi Inclusi](#-esempi-inclusi)

---

## 🪟 Come accedere alla Versione Windows 11

Questo branch (**`Versione-Windows-11`**) mantiene esattamente la stessa struttura e lo stesso aspetto dell'applicazione, aggiungendo compatibilità completa e ottimizzazioni dedicate per **Windows 11** (isolamento dal Dark Mode di sistema, DPI/scaling 125%–175%, font ricchi di simboli matematici, supporto sia a `.` sia a `,` del tastierino numerico italiano e notazione scientifica nei campi numerici, gestione sicura dei percorsi con spazi/accenti e cartelle temporanee in `%TEMP%`).

### 1. Da GitHub (senza usare comandi Git)
1. Apri la repository: `https://github.com/francesco-barbato/PLV-Solver`.
2. In alto a sinistra sopra l'elenco dei file, apri il menu a tendina dei branch (**`main`**) e seleziona **`Versione-Windows-11`** *(oppure vai direttamente su `https://github.com/francesco-barbato/PLV-Solver/tree/Versione-Windows-11`)*.
3. Clicca su **`<> Code` → `Download ZIP`** ed estrai l'archivio in una cartella sul PC.
4. Fai doppio clic su **`Avvia_Windows11.bat`** (la prima volta installa automaticamente eventuali librerie mancanti se Python 3.9+ è già installato).

### 2. Da riga di comando (Git)
```bash
# Clonazione diretta del branch Windows 11:
git clone -b Versione-Windows-11 https://github.com/francesco-barbato/PLV-Solver.git

# Oppure, se hai già clonato il repository in locale:
git fetch origin
git checkout Versione-Windows-11
```

### 3. Guide passo-passo incluse
- **[GUIDA_WINDOWS_11.md](GUIDA_WINDOWS_11.md)** — Guida dettagliata all'esame specifica per Windows 11 (pensata anche per chi non conosce GitHub, programmazione o Scienza delle Costruzioni).
- **[GUIDA_ESAME.md](GUIDA_ESAME.md)** — Guida generale all'utilizzo e alla trascrizione del report sul foglio d'esame.

---

## 🚀 Caratteristiche Principali

### 1. Aritmetica Esatta e Calcolo Simbolico (SymPy)
- Risoluzione puramente analitica tramite **SymPy**: mantiene frazioni esatte e radicali notevoli ($\sqrt{2}, \sqrt{3}, \dots$), evitando approssimazioni numeriche e cancellazioni di cifre significative tipiche dei calcoli in virgola mobile.
- Formule analitiche continue per lo sforzo normale $N(s)$, taglio $T(s)$ e momento flettente $M(s)$ lungo l'ascissa locale $s \in [0, L]$ di ciascuna asta.

### 2. Ampia Varietà di Vincoli e Svincoli
- **Vincoli a terra**:
  - Incastro perfetto
  - Cerniera fissa
  - Carrello a scorrimento orizzontale o verticale
  - Carrello inclinato con angolo di scorrimento arbitrario
  - Doppio pendolo / Bipendolo (orizzontale, verticale o con asse inclinato)
  - Pendolo / Biella inclinata (modellata come elemento monoassiale tirante/puntone)
- **Vincoli interni**:
  - Cerniere interne su nodi (svincolo del momento flettente tra le aste collegate).
- **Cedimenti vincolari anelastici**:
  - Cedimenti impressi traslazionali ($\Delta u, \Delta v$), rotazionali ($\Delta w$) e assiali lungo l'asse del pendolo ($\Delta n$).

### 3. Azioni Esterne e Distorsioni
- Forze concentrate ($F_x, F_y$) e coppie concentrate ($M$) applicate sui nodi.
- Carichi distribuiti trasversali uniformi ($q_{\text{trans}}$) o con componenti globali cartesiane ($q_x, q_y$).
- Distorsioni termiche uniformi e differenziali (a farfalla) con coefficiente di dilatazione termica $\alpha$ e gradienti $T^+$, $T^-$.

### 4. Risoluzione tramite Principio dei Lavori Virtuali
- **Analisi cinematica preliminare**: conteggio automatico dei gradi di libertà ($3n$), gradi di vincolo e verifica del grado di iperstaticità.
- **Decomposizione automatica**:
  - **Sistema $S_0$**: sistema principale isostatico soggetto ai carichi reali, cedimenti e variazioni termiche.
  - **Sistema $S_1$**: sistema ausiliario isostatico con incognita iperstatica unitaria $X_1 = 1$.
- Calcolo analitico dei coefficienti di deformabilità ($\eta_{10}, \eta_{11}$) e dei lavori virtuali dei cedimenti e delle distorsioni termiche ($\eta_{1,\text{ced}}, \eta_{1,\text{term}}$).
- Determinazione dell'incognita iperstatica $X_1$ mediante l'equazione di congruenza elastica (equazione di Müller-Breslau):
  $$\eta_{10} + \eta_{11} X_1 = \eta_{1,\text{ced}} + \eta_{1,\text{term}}$$
- Sovrapposizione degli effetti per le sollecitazioni finali reali:
  $$M(s) = M_0(s) + X_1 \cdot M_1(s), \quad T(s) = T_0(s) + X_1 \cdot T_1(s), \quad N(s) = N_0(s) + X_1 \cdot N_1(s)$$

### 5. Controlli di Rigore Statico e Verifiche di Resistenza
- **Controllo di equilibrio nodale globale e locale**: verifica analitica automatica dell'equilibrio alla traslazione e alla rotazione di ogni singolo nodo ($\sum F_x = 0, \sum F_y = 0, \sum M = 0$).
- **Verifica delle bielle/pendoli**: identificazione dello stato di sollecitazione (tirante con trazione o puntone con compressione) e verifica a sforzo normale.
- **Verifica delle tensioni di Navier**: ricerca dei punti critici lungo le aste e calcolo di $\sigma_{\max} = \frac{|N|}{A} + \frac{|M|}{W}$, con calcolo del **Margine di Sicurezza ($MS$)**:
  $$MS = \frac{\sigma_{\text{amm}}}{\sigma_{\max}} - 1$$
- **Dimensionamento automatico della sezione**: calcolo dell'altezza $h$, della base $b$ o delle dimensioni con rapporto di forma fissato ($h/b$) per garantire un margine di sicurezza target.

### 6. Interfaccia Grafica e Reportistica (PyQt6 & Matplotlib)
- **Editor interattivo**: inserimento e modifica visuale di nodi, aste, vincoli, carichi e proprietà dei materiali.
- **Diagrammi quotati ad alta risoluzione**:
  - Schemi di corpo libero di tutte le aste con reazioni vincolari interne ed esterne.
  - Diagrammi di sollecitazione ($N, T, M$) completi di valori massimi, minimi e punti di nullo.
  - Confronto visivo tra i sistemi $S_0$ e $S_1$.
- **Reportistica HTML & Stampa PDF**: generazione di un report formattato con passaggi analitici espliciti, tabelle e grafici incorporati, con esportazione diretta in formato PDF (layout A2/A3/A4 landscape).

---

## 📂 Architettura del Progetto

```
SolutorePLV/
│
├── README.md                    # Documentazione ufficiale del progetto
├── GUIDA_WINDOWS_11.md          # Guida completa passo-passo per utenti Windows 11
├── GUIDA_ESAME.md               # Guida generale all'esame
├── Avvia_Windows11.bat          # Avvio rapido 1-clic su Windows 11
├── .gitignore                   # Configurazione esclusioni Git
├── .vscode/
│   └── settings.json            # Configurazioni dell'ambiente di sviluppo
│
├── JSON/                        # Casi studio ed esercizi d'esame reali (formato JSON)
│   ├── 060921.json
│   ├── 140722.json
│   ├── 150323.json
│   ├── 151222.json
│   ├── 160922.json
│   ├── 170522.json
│   ├── 200921.json
│   ├── 240621.json
│   ├── 261021.json
│   ├── 261026.json
│   ├── 270721.json
│   └── 280223.json
│
└── py/                          # Codice sorgente Python e risorse
    ├── ui_strutture.py          # Interfaccia Grafica Desktop (PyQt6 + Matplotlib Canvas)
    ├── solutore_plv.py          # Motore computazionale, algebra simbolica e PLV (SymPy)
    ├── SolutorePLV_Win11.spec   # Specfile PyInstaller per Windows 11 (.exe)
    ├── Crea_EXE_Windows11.bat   # Script di compilazione automatica .exe su Windows 11
    ├── SolutorePLV.spec         # Specfile per la creazione di app macOS con PyInstaller
    ├── ui_strutture.spec        # Specfile alternativo per PyInstaller
    ├── icona.ico                # Icona dell'applicazione per Windows 11
    ├── icona.icns               # Icona dell'applicazione per macOS
    └── icona.jpg                # Risorsa grafica per l'interfaccia
```

---

## 🛠️ Prerequisiti e Installazione

### Requisiti
- **Python**: versione 3.9 o superiore (consigliato Python 3.9, 3.10 o 3.11).
- **Sistema Operativo**: macOS, Linux o Windows.

### 1. Clonazione del repository (o download della cartella)
```bash
git clone <URL_REPO>
cd SolutorePLV
```

### 2. Creazione dell'ambiente virtuale (consigliato)
```bash
python3 -m venv .venv
source .venv/bin/activate       # Su macOS / Linux
# .venv\Scripts\activate       # Su Windows
```

### 3. Installazione delle dipendenze
Installare i pacchetti necessari tramite `pip`:
```bash
pip install PyQt6 matplotlib numpy sympy typing-extensions
```

---

## 🖥️ Modalità d'Uso

### 1. Avvio dell'Interfaccia Grafica (GUI)

Per avviare l'applicazione desktop:
```bash
python3 py/ui_strutture.py
```

Dall'interfaccia è possibile:
1. **Caricare una struttura esistente**: cliccare su *"Apri JSON"* e selezionare uno dei file presenti nella cartella `JSON/`.
2. **Creare o modificare la struttura**:
   - Gestire nodi (coordinate $x, y$).
   - Gestire aste (nodi di estremità).
   - Impostare vincoli a terra, cerniere interne e carichi.
   - Definire materiale ($E, \sigma_{\text{amm}}$) e dimensioni della sezione rettangolare ($b \times h$).
3. **Risolvere la struttura**: cliccare su *"Calcola Struttura (PLV)"*.
4. **Consultare il Report e i Grafici**:
   - Visualizzare l'analisi passo-passo nel tab *"Report Analitico"*.
   - Consultare i diagrammi quotati delle sollecitazioni nel tab *"Diagrammi"*.
   - Esportare l'intero report in formato PDF tramite il pulsante dedicato.

### 2. Utilizzo come Libreria / Script Python

Il motore di calcolo può essere importato e utilizzato all'interno di script o pipeline automatiche:

```python
import json
from py.solutore_plv import risolvi_struttura_plv

# 1. Carica il modello della struttura
with open("JSON/060921.json", "r", encoding="utf-8") as f:
    dati_struttura = json.load(f)

# 2. Risolvi la struttura con il PLV
# Parametri:
#   dati_json: dizionario con nodi, aste, vincoli e carichi
#   out_dir: cartella in cui salvare i diagrammi generati
#   incognita: specifica eventuale dell'incognita iperstatica (o None per automatica)
report_html, path_s0_s1, path_diagrammi = risolvi_struttura_plv(
    dati_json=dati_struttura,
    out_dir="py/temp_output",
    incognita=None
)

# 3. Utilizza i risultati
print("Calcolo completato con successo!")
print(f"Grafico S0/S1 salvato in: {path_s0_s1}")
print(f"Diagrammi M-T-N salvati in: {path_diagrammi}")
```

---

## 📄 Formato dei Dati (JSON Schema)

I file di input nella cartella `JSON/` definiscono la topologia e le azioni della struttura:

```json
{
  "nodes": [
    { "id": 1, "x": 0.0, "y": 0.0 },
    { "id": 2, "x": 0.0, "y": 2.0 },
    { "id": 3, "x": 3.0, "y": 2.0 }
  ],
  "elements": [
    { "id": 1, "n1": 1, "n2": 2 },
    { "id": 2, "n1": 2, "n2": 3 }
  ],
  "supports": [
    {
      "node_id": 1,
      "fix_u": true,
      "fix_v": true,
      "fix_w": true,
      "angle": 0.0
    },
    {
      "node_id": 3,
      "fix_v": true,
      "angle": 0.0
    }
  ],
  "internal_hinges": [
    { "node_id": 2 }
  ],
  "point_loads": [
    { "node_id": 2, "Fx": 1000.0, "Fy": 0.0, "M": 0.0 }
  ],
  "element_loads": [
    { "elem_id": 2, "tipo": "costante", "q_trans": 500.0 }
  ],
  "thermal_loads": [],
  "material": {
    "E": 2.1e11,
    "sigma_amm": 235e6,
    "MS_target": 0.0,
    "MS_operator": ">",
    "verification_type": "bending_only"
  },
  "section": {
    "b": 0.2,
    "h": 0.4,
    "design": "none"
  }
}
```

---

## 🔄 Flusso di Risoluzione (Workflow PLV)

```
[ Input Struttura (JSON o GUI) ]
              │
              ▼
[ Normalizzazione & Validazione Dati ] (Aritmetica esatta SymPy)
              │
              ▼
[ Analisi Cinematica ] (Gradi di libertà, gradi di vincolo, iperstaticità)
              │
              ▼
[ Definizione del Sistema Isostatico Principale ]
      ├──> Sistema S₀ (Carichi reali, cedimenti, variazioni termiche)
      └──> Sistema S₁ (Incognita iperstatica unitaria X₁ = 1)
              │
              ▼
[ Calcolo Analitico Caratteristiche della Sollecitazione ] (N(s), T(s), M(s))
              │
              ▼
[ Integrazione del Lavoro Virtuale (PLV) ]
      ├──> Coefficienti di deformabilità: η₁₀, η₁₁
      └──> Lavori anelastici: η₁,ced, η₁,term
              │
              ▼
[ Risoluzione Equazione di Congruenza ] ──> Calcolo incognita X₁
              │
              ▼
[ Sovrapposizione degli Effetti & Verifiche di Equilibrio ]
      ├──> Sollecitazioni finali reali per ogni asta
      ├──> Equilibrio di ciascun nodo (ΣFx = 0, ΣFy = 0, ΣM = 0)
      ├──> Verifica di resistenza (Navier) e Margine di Sicurezza (MS)
      └──> Dimensionamento ottimale della sezione (opzionale)
              │
              ▼
[ Output Grafico e Relazione Tecnica ] (Diagrammi Matplotlib + Report PDF/HTML)
```

---

## 📦 Compilazione Standalone (PyInstaller)

### Su Windows 11
È incluso lo script automatico **`py\Crea_EXE_Windows11.bat`** (basato su `py/SolutorePLV_Win11.spec`), oppure puoi eseguire da Prompt dei comandi nella cartella `py/`:

```bat
pip install pyinstaller
pyinstaller --noconfirm SolutorePLV_Win11.spec
```

L'eseguibile standalone verrà generato in `py\dist\SolutorePLV\SolutorePLV.exe`.

### Su macOS
È possibile generare il pacchetto `.app` utilizzando **PyInstaller** e la configurazione inclusa in `SolutorePLV.spec`:

```bash
# 1. Installa PyInstaller
pip install pyinstaller

# 2. Esegui la build a partire dal file .spec
pyinstaller py/SolutorePLV.spec
```

Al termine della compilazione:
- L'eseguibile o il bundle applicativo (es. `SolutorePLV.app`) sarà posizionato nella cartella `py/dist/`.
- I file temporanei di compilazione saranno posizionati nella cartella `py/build/` (ignorata da Git).

---

## 📚 Esempi Inclusi

La cartella `JSON/` raccoglie una vasta serie di esercizi d'esame universitari reali (codificati per data dell'appello, es. `140722.json` = appello del 14 Luglio 2022). Ognuno di essi copre casistiche specifiche di particolare interesse didattico:
- Portali simmetrici e asimmetrici
- Travi continue su più appoggi con sbalzi
- Strutture con aste inclinate a $45^\circ$ e carichi ortogonali
- Bielle / pendoli tiranti e puntoni
- Sistemi soggetti a carichi distribuiti uniformi e distorsioni termiche
- Cedimenti vincolari impressi

---

## 📜 Licenza

Questo progetto è distribuito per scopi didattici, accademici e di ricerca nell'ambito della meccanica delle strutture e della Scienza delle Costruzioni.

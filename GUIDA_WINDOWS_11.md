# Guida Completa Windows 11: Come Scaricare, Usare e Copiare i Dati del Solutore PLV

> Questa guida è pensata per chi usa **Windows 11** e **non conosce nulla di GitHub, di programmazione né di Scienza delle Costruzioni**. Seguendo i passaggi nell'ordine indicato potrai scaricare la versione dedicata a Windows 11, avviarla con un clic (o creare il file `.exe`), inserire i dati dell'esame senza errori di unità o virgole e trascrivere passo passo il report sul foglio.

---

## Parte 1 — Come scaricare e avviare la Versione Windows 11 (senza comandi Git)

La versione ottimizzata per Windows 11 si trova sul branch **`Versione-Windows-11`** della pagina GitHub:
`https://github.com/francesco-barbato/PLV-Solver/tree/Versione-Windows-11`

### Passo 1: Scarica lo ZIP dal branch giusto (`Versione-Windows-11`)
1. Apri il browser (Edge o Chrome) e vai su **`https://github.com/francesco-barbato/PLV-Solver`**.
2. **PASSAGGIO FONDAMENTALE:** in alto a sinistra sopra la lista dei file vedi un pulsante grigio con scritto **`main`**. Cliccalo e seleziona **`Versione-Windows-11`** (oppure apri direttamente `https://github.com/francesco-barbato/PLV-Solver/tree/Versione-Windows-11`).
3. Clicca sul pulsante verde **`<> Code`** in alto a destra.
4. Clicca su **`Download ZIP`**.
5. Apri la cartella **`Download`**, fai **clic destro** sul file `.zip` scaricato (`PLV-Solver-Versione-Windows-11.zip`) → **`Estrai tutto...`** → scegli una cartella comoda (es. `Documenti` o `Desktop`).

### Passo 2: Installa Python (solo la prima volta)
1. Vai su `https://www.python.org/downloads/` e clicca sul pulsante giallo per scaricare l'ultima versione di Python per Windows (`3.9` o successiva).
2. Apri il file `.exe` scaricato.
3. **OBBLIGATORIO — NON DIMENTICARLO:** prima di cliccare su *Install Now*, metti la spunta sulla casella in basso **`Add python.exe to PATH`**.
4. Clicca su **`Install Now`** e attendi la fine dell'installazione.

### Passo 3: Avvia il programma con un doppio clic (`Avvia_Windows11.bat`)
1. Entra nella cartella estratta `PLV-Solver-Versione-Windows-11`.
2. Fai **doppio clic** sul file **`Avvia_Windows11.bat`**.
3. **Se Windows Defender SmartScreen mostra un avviso blu** (*"PC protetto da Windows"*):
   - Clicca su **`Ulteriori informazioni`** → poi su **`Esegui comunque`** (succede normalmente per file `.bat` scaricati da Internet).
4. La prima volta il launcher controlla e installa automaticamente le librerie necessarie (`PyQt6`, `matplotlib`, `numpy`, `sympy`, `typing-extensions` — attendi 1-2 minuti con connessione Internet attiva). Le volte successive si apre **subito**.

> [!TIP]
> **Vuoi creare un vero file `SolutorePLV.exe` (opzionale)?**
> Entra nella cartella `py` e fai doppio clic su **`Crea_EXE_Windows11.bat`**: compilerà automaticamente l'eseguibile standalone in **`py\dist\SolutorePLV\SolutorePLV.exe`** pronto da copiare o avviare senza console.

*(In alternativa manuale da Prompt dei comandi `cmd`: entra nella cartella estratta, esegui `pip install PyQt6 matplotlib numpy sympy typing-extensions` la prima volta e poi `python py\ui_strutture.py`.)*

---

## Parte 2 — Come usare il programma e inserire l'esercizio d'esame

La finestra è identica in tutto e per tutto e funziona senza problemi anche se su Windows 11 hai attivato il **Tema Scuro** o lo scaling dello schermo al **125% / 150% / 175%**:
- **Colonna sinistra** (tre schede in alto):
  - **`Struttura`**: costruzione guidata dall'`1.` all'`8.` + `Riepilogo`
  - **`Materiale e sezione`**: dati del materiale, sezione assegnata e opzioni di verifica (Punto 3) / dimensionamento (Punto 4)
  - **`Progetto`**: pulsanti **`Apri JSON da File`** e **`Salva JSON`**
- **Colonna destra** (due schede in alto):
  - **`Anteprima Live`**: disegno in tempo reale (usala a ogni passo per verificare che il disegno coincida con quello dell'esame!)
  - **`Report di Calcolo`**: tutto lo svolgimento completo da ricopiare
- **Barra in basso**: menu **`Incognita iperstatica X₁:`**, pulsante verde **`▶ Risolvi struttura`** e pulsante **`🖨️ PDF`**.

> [!TIP]
> **Tastierino numerico italiano e notazione scientifica su Windows 11:**
> - Nei campi numerici puoi scrivere **indifferentemente con il punto `.` oppure con la virgola `,`** (ad esempio `1.5` oppure `1,5` dal tastierino numerico): il programma li capisce entrambi senza mai moltiplicare per 10!
> - Puoi inserire direttamente anche la notazione scientifica con la lettera `e`: ad esempio `2.1e11` (o `2,1e11`) per $2{,}1 \cdot 10^{11}$ e `1.7e-5` (o `1,7e-5`) per $1{,}7 \cdot 10^{-5}$.

---

### Tabella di conversione delle unità di misura (LA REGOLA D'ORO)

Il testo dell'esame usa spesso `kN`, `cm`, `MPa`, `GPa`, mentre il programma vuole **metri (`m`)**, **Newton (`N`)** e **Pascal (`Pa`)** (tranne i cedimenti lineari che sono già in millimetri!).

| Cosa leggi nel testo dell'esame | Unità del testo | Cosa scrivere nel programma | Esempio pratico |
| :--- | :--- | :--- | :--- |
| Lunghezze aste / coordinate nodi | `cm` | Dividi per `100` → **`m`** | `150 cm` → `1.5` (o `1,5`) |
| Dimensioni sezione `b`, `h` | `cm` o `mm` | Dividi per `100` (cm) o `1000` (mm) → **`m`** | `b = 20 cm` → `0.2` |
| Forze `F` o carichi distribuiti `q` | `kN` o `kN/m` | Moltiplica per `1000` → **`N`** o **`N/m`** | `12 kN` → `12000` |
| Coppie / Momenti concentrati `M` | `kN·m` | Moltiplica per `1000` → **`N·m`** | `5 kN·m` → `5000` |
| Modulo elastico `E` | `GPa` | Moltiplica per $10^9$ (`e9`) → **`Pa`** | `210 GPa` → `2.1e11` |
| Modulo elastico `E` | `MPa` (`N/mm²`) | Moltiplica per $10^6$ (`e6`) → **`Pa`** | `30000 MPa` → `3e10` |
| Tensione ammissibile `σ_amm` | `MPa` (`N/mm²`) | Moltiplica per $10^6$ (`e6`) → **`Pa`** | `240 MPa` → `2.4e8` |
| Cedimenti lineari (`Δx`, `Δy`, `Δn`) | `mm` o `cm` | Il campo è già in **`mm`** (`1 cm = 10 mm`) | Abbassamento `2 cm` → `-20` |
| Rotazione imposta (`Δγ`) | `gradi (°)` o `rad` | Il campo è in **`rad`** (`gradi × π / 180`) | `1°` → `0.017453` |
| Radici nel testo (es. `700√2`) | `√2 ≈ 1.414214` | Scrivi il decimale con 4-6 cifre | `700√2` → `989.949494` |

---

### Passo per passo nella scheda `Struttura` (Sezioni 1 → 8)

#### `1. Nodi` — I punti della figura
Inserisci gli spigoli, i punti dove cambiano direzione le travi, dove ci sono vincoli, cerniere interne o carichi isolati:
- Usa **`Per coordinate:`** (`X (m)` destra, `Y (m)` alto) partendo da un'origine comoda `(0, 0)`, oppure **`Da un nodo:`** (`L (m)` e `Ang (°)`) se hai lunghezza e angolo rispetto all'orizzontale.

#### `2. Aste` — I tratti dritti fra i nodi
Collega due nodi scegliendo **`Da nodo:`** e **`A nodo:`**.
- **Regola chiave:** la direzione `Da nodo -> A nodo` stabilisce dove inizia (`s = 0`) e finisce (`s = L`) l'asta, e qual è il lato **sinistro/destro** per il carico distribuito perpendicolare e per la temperatura!

#### `3. Vincoli` — Gli appoggi esterni
Abbina il disegno sul foglio al tipo esatto nel menu:

| Simbolo sul foglio | Voce nel menu | Effetto pratico |
| :--- | :--- | :--- |
| Tratteggio fisso sul muro | **`Incastro`** | Blocca X, Y e rotazione |
| Triangolo senza rotelline | **`Cerniera`** | Blocca X e Y, lascia ruotare |
| Triangolo con rotelle orizzontali | **`Carrello X (Scorre Orizz.)`** | Scorre in X, blocca verticale Y (`Ry`) |
| Triangolo con rotelle verticali | **`Carrello Y (Scorre Vert.)`** | Scorre in Y, blocca orizzontale X (`Rx`) |
| Triangolo su piano obliquo | **`Carrello Inclinato`** | Piano inclinato di `Angolo(°)` da `+x`; reazione `Rn` perpendicolare (`angolo + 90°`) |
| Bielle parallele orizzontali | **`Bipendolo (Scorre Orizz.)`** | Scorre in X, blocca Y e rotazione |
| Bielle parallele verticali | **`Bipendolo (Scorre Vert.)`** | Scorre in Y, blocca X e rotazione |
| Asta singola collegata a terra | **`Pendolo / biella inclinata (tirante o puntone)`** | `Angolo(°)` da `+x` andando **dal nodo verso terra**; lascia `Lungh. biella (m) = 0` se non data (`0 = n.i.`); `Rn > 0` = tirante |
| Doppio pendolo inclinato | **`Doppio pendolo inclinato (bipendolo)`** | Come il pendolo per l'angolo d'asse, ma blocca anche la rotazione |

#### `4. Carichi sui nodi (forze e coppie)`
- Inserisci `Fx (N)`, `Fy (N)`, `M (N·m)` oppure **`F inclinata (N)`** con `angolo (°)` misurato in senso **antiorario** da destra (`0°` = destra, `90°` = su, `180°` = sinistra, `-90°` = giù).
- Dopo l'aggiunta controlla sempre nell'**`Anteprima Live`** che la freccia dritta o curva punti nel verso esatto del foglio (se è opposta, modificala col segno `-`).

#### `5. Carichi distribuiti sulle aste`
- **`Tipo`**: `Costante`, `Triangolare (0 ➔ q)` oppure `Triangolare (q ➔ 0)`.
- **`Direzione`**:
  - `Perpendicolare all'asta`: con `q > 0` spinge sul lato **DESTRO** guardando `Da nodo -> A nodo` (su una trave orizzontale da sinistra a destra spinge **verso il basso**).
  - `Verticale (+ verso il basso)` o `Orizzontale (+ verso destra)`.
- **`Riferito a`**: `per unità di lunghezza d'asta` (standard) o `per unità di proiezione` (se specificato sul disegno per aste inclinate).
- **`Tratto`**: `Tutta l'asta`, `Prima metà (da Nodo Inizio)`, `Seconda metà (verso Nodo Fine)` oppure `Tratto personalizzato` (`da` e `a` in frazioni da `0` a `1`).

#### `6. Cerniere interne`
Seleziona il **`Nodo:`** dello snodo interno. Se arrivano **3 o più aste**, scrivi in **`Aste incernierate:`** gli ID separati da virgola (es. `2` oppure `2,3`); se arrivano 2 aste lascia vuoto. **`Carichi su asta:`** di norma resta `auto`.

#### `7. Cedimenti vincolari / rotazioni imposte`
- `Δx (mm)` (`+` verso destra), `Δy (mm, ↑ +)` (**un abbassamento verso il basso va inserito NEGATIVO, es. `-15`**!), `Δn (mm, car. inclinato)` (`+` lungo `angolo + 90°`), `Δγ (rad, antiorario +)` (in **radianti**, positivo **antiorario**).

#### `8. Variazione termica (gradiente nello spessore)`
Inserisci **`α (1/°C)`** (es. `1.7e-5`) e scegli tra **`Due fibre (sinistra / destra)`**, **`Farfalla ±ΔT (sinistra +ΔT, destra −ΔT)`** o **`Differenza totale ΔT (sinistra +ΔT/2, destra −ΔT/2)`** (sempre guardando da `Da nodo -> A nodo`: su trave da sinistra a destra, sinistra = sopra, destra = sotto).

---

### Scheda `Materiale e sezione` e salvataggio in `Progetto`
1. Imposta **`E (Pa)`** (fino a `1e15`), **`Sigma amm (Pa)`** (fino a `1e12`), **`Base b (m)`**, **`Altezza h (m)`**, **`MS verifica (punto 3)`**, **`Criterio`** (`Presso/Tensoflessione (N/A + M/W)` oppure `Solo Flessione (M/W)`), **`EI assegnato (N·m²)`** (`0` = automatico), **`Dimensionamento`** (`Nessuno`, `Incognita h`, `Incognita b`, `Rapporto h/b`), **`MS richiesto (punto 4)`** ed eventuale **`Sezione del Pendolo (biella)`**.
2. Apri la scheda **`Progetto`** e clicca su **`Salva JSON`** per salvare l'esercizio in una cartella (es. `Documenti`) prima di avviare il calcolo.

---

## Parte 3 — Avviare il calcolo e copiare i risultati sul foglio

1. Scegli una reazione nel menu **`Incognita iperstatica X₁:`** in basso e premi **`▶ Risolvi struttura`**.
2. **Se appare *"Analisi Interrotta — Struttura Labile"*:** togliere proprio quella reazione fa cadere la struttura; scegli un'altra voce dal menu **`Incognita iperstatica X₁:`** e premi di nuovo **`▶ Risolvi struttura`**.
3. Al termine, clicca in alto a destra su **`Report di Calcolo`** (e usa **`🖨️ PDF`** in basso se vuoi esportarlo in PDF).

| Sezione del `Report di Calcolo` | Cosa copiare sul foglio d'esame |
| :--- | :--- |
| **Box `Verifica Cinematica`** | Conteggio g.d.l. (`3n`), gradi di vincolo (`v`) e grado di iperstaticità `h = 1`. |
| **`1. Dati della struttura...`** | Ricontrolla tutte le conversioni, `EI`, sezione assegnata (`b`, `h`, `A`, `I`, `W`) ed eventuale parametro del Punto 4 (`h/b = ...`). |
| **`2. Scelta dell'incognita iperstatica...`** | Disegno dello schema isostatico equivalente (`X₁`). Per cerniere interne: frecce **blu** = azione sul tronco precedente, **viola** = azione uguale e opposta sul successivo. |
| **`3. Sistema S₀`** | I 3 schemi di `S₀`, tabella reazioni, equilibrio per ogni blocco **`TRATTO`**, espressioni $N_0(s), T_0(s), M_0(s)$ e diagrammi globali. |
| **`4. Sistema S₁`** | Schema scarico con `X₁ = 1`, reazioni, $N_1(s), T_1(s), M_1(s)$ e diagrammi globali. |
| **`5. Coefficienti di deformabilità (PLV)`** | Tabella degli integrali di Mohr ($\delta_{10}$, $\delta_{11}$) asta per asta. |
| **`6. Equazione di congruenza e calcolo di X₁`** | Equazione di congruenza (con eventuali cedimenti/termica) e valore finale `X₁`. |
| **`7. Sistema reale`** | Reazioni definitive, forze nelle cerniere, espressioni reali $N(s), T(s), M(s)$, punti stazionari **`s*` / `M_extr`** e diagrammi finali: **copia sempre il momento `M` dal lato mostrato dal grafico rosso (lato fibre tese)**! |
| **`8. ... (Punto 3 — Navier)`** | Tabella sezioni candidate, $\sigma_{\max}$, margine $MS$ ed esito **`VERIFICATA` / `NON VERIFICATA`** rispetto a `MS verifica (punto 3)`. |
| **`9. Dimensionamento della sezione (Punto 4)`** | Disequazione in `b` e `h` (con eventuale $X_1(b,h)$ se presenti cedimenti/ΔT), valori **`b_min`**, **`h_min`** o **`[b_min, b_max]`**, eventuale avviso di **cambio della sezione critica** e controllo $MS = MS_{\text{req}}$. |

> [!WARNING]
> **Limiti da ricordare all'esame:** il solutore risolve **telai piani 1 volta iperstatici (`h = 1`)** con **deformabilità flessionale (`EI`)** negli integrali di Mohr.

---

## Parte 4 — Checklist rapida di controllo (30 secondi prima di consegnare)

1. L'**`Anteprima Live`** corrisponde per forma, vincoli, cerniere e versi delle frecce al foglio d'esame?
2. Hai convertito `cm → m`, `kN → N`, `GPa/MPa → Pa`?
3. Un cedimento verticale verso il basso ha il **segno `-`** in `Δy (mm, ↑ +)`?
4. Se c'è temperatura `ΔT`, hai inserito **`α (1/°C)`** e verificato sinistra/destra rispetto a `Da nodo -> A nodo`?
5. Hai verificato sia `MS verifica (punto 3)` sia `MS richiesto (punto 4)`?
6. Hai salvato l'esercizio con **Progetto → Salva JSON**?

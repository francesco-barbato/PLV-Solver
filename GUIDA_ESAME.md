# Guida Completa all'Esame: Come Scaricare, Usare e Copiare i Dati del Solutore PLV

> [!TIP]
> **Stai usando Windows 11?** Su questo branch è disponibile la guida dedicata passo-passo per Windows 11 (avvio con un clic tramite `Avvia_Windows11.bat` o compilazione `.exe`, gestione virgola/punto sul tastierino e Theme/HiDPI): apri [GUIDA_WINDOWS_11.md](GUIDA_WINDOWS_11.md).

> Questa guida è pensata per chi **non conosce nulla di GitHub, di programmazione né di Scienza delle Costruzioni**. Seguendo i passaggi nell'ordine indicato riuscirai a installare il programma, inserire correttamente il disegno dell'esame senza errori di unità di misura o segni, e trascrivere passo passo i conti sul foglio.

---

## Parte 1 — Come scaricare e aprire il programma (senza usare comandi Git)

Il programma si trova su GitHub a questo indirizzo:
`https://github.com/francesco-barbato/PLV-Solver`

Su GitHub è presente il codice Python del programma (i file nella cartella `py/`) insieme alla cartella `JSON/` con alcuni esercizi già svolti di esempio.

### Metodo A — Su Mac, se hai ricevuto il file pronto `SolutorePLV.app`
Se qualcuno ti ha già passato l'applicazione compilata per macOS (`SolutorePLV.app`):
1. Fai doppio clic su `SolutorePLV.app`.
2. Se macOS mostra un avviso di sicurezza (*"Impossibile aprire perché proviene da uno sviluppatore non identificato"*):
   - Fai **clic destro** (o `Control + clic`) su `SolutorePLV.app` → **Apri** → **Apri**.
   - Se non basta, apri **Impostazioni di Sistema → Privacy e sicurezza**, scorri in basso e clicca su **Apri comunque**.

### Metodo B — Scaricare da GitHub (Windows e Mac, funziona sempre)

#### Passo 1: Scarica la cartella ZIP da GitHub
1. Apri il browser (Chrome, Safari, Edge) e vai su `https://github.com/francesco-barbato/PLV-Solver`.
2. Clicca sul pulsante verde **`<> Code`** in alto a destra della lista file.
3. Clicca su **`Download ZIP`**.
4. Trova il file scaricato (nella cartella `Download`) ed **estrailo / decomprimilo**:
   - **Su Windows**: clic destro sullo `.zip` → *Estrai tutto...* → scegli ad esempio il Desktop o Documenti.
   - **Su Mac**: doppio clic sul file `.zip` per creare automaticamente la cartella `PLV-Solver-main`.

#### Passo 2: Installa Python (solo la prima volta)
Il programma richiede **Python 3.9 o successivo**.
- **Su Windows**:
  1. Vai su `https://www.python.org/downloads/` e scarica l'ultima versione.
  2. **FONDAMENTALE**: appena apri il file di installazione, metti la spunta sulla casella in basso **`Add python.exe to PATH`** prima di cliccare su *Install Now*.
- **Su Mac**:
  - Apri l'app **Terminale** (premi `Cmd + Spazio`, scrivi `Terminale` e premi Invio).
  - Scrivi `python3 --version` e premi Invio. Se compare una versione `3.9` o superiore sei già pronto; se chiede di installare i *Command Line Developer Tools* o se non è presente, scarica Python da `https://www.python.org/downloads/`.

#### Passo 3: Installa le librerie necessarie (solo la prima volta)
Apri il **Terminale** (Mac) oppure il **Prompt dei comandi** (Windows: tasto Start, scrivi `cmd` e premi Invio), entra nella cartella appena estratta e installa le librerie richieste con il comando copia-incolla.

1. Scrivi `cd ` (con uno spazio dopo `cd`), **trascina dentro la finestra del terminale la cartella estratta** (`PLV-Solver-main`) e premi Invio.
2. Copia e incolla questo comando e premi Invio (attendi 1-2 minuti finché finisce):

```bash
pip install PyQt6 matplotlib numpy sympy typing-extensions
```

*(Su Mac, se dice `pip: command not found`, usa invece `python3 -m pip install PyQt6 matplotlib numpy sympy typing-extensions`.)*

#### Passo 4: Avvia il programma
Ogni volta che vuoi aprire il programma, dentro la cartella principale del progetto esegui:
- **Su Mac** (o Linux):
  ```bash
  python3 py/ui_strutture.py
  ```
- **Su Windows** (se `python3` non viene riconosciuto, prova `python`):
  ```bash
  python py/ui_strutture.py
  ```

> [!TIP]
> **Trucco salva-tempo prima dell'esame:** fai una prova di apertura a casa almeno un giorno prima e apri un file dalla cartella `JSON/` usando **Progetto → Apri JSON da File** per verificare che tutto funzioni senza intoppi.

---

## Parte 2 — Capire la schermata e inserire i dati dell'esame

Appena si apre la finestra trovi due colonne:
- **A sinistra**: tre schede in alto:
  - **`Struttura`** (dove costruisci il disegno dell'esercizio, dall'1 all'8)
  - **`Materiale e sezione`** (dove metti i dati numerici del materiale e delle dimensioni della trave)
  - **`Progetto`** (dove trovi **`Apri JSON da File`** e **`Salva JSON`** per salvare o ricaricare il tuo esercizio)
- **A destra**: due schede in alto:
  - **`Anteprima Live`** (il disegno che si aggiorna **in diretta** mentre inserisci i dati: usala sempre per controllare a vista che il disegno assomigli al 100% a quello del foglio d'esame!)
  - **`Report di Calcolo`** (qui appariranno tutti i calcoli da copiare)
- **In basso**: la barra con **`Incognita iperstatica X₁:`**, **`▶ Risolvi struttura`** e **`🖨️ PDF`**.

---

### Tabella di conversione delle unità di misura (LA REGOLA D'ORO)

Quasi tutti gli errori all'esame succedono perché il testo usa `kN`, `cm`, `MPa`, `GPa`, mentre il programma chiede **metri (`m`)**, **Newton (`N`)** e **Pascal (`Pa`)** (tranne i cedimenti che sono già in millimetri!).

| Cosa leggi nel testo dell'esame | Unità del testo | Cosa devi scrivere nel programma | Esempio pratico |
| :--- | :--- | :--- | :--- |
| Lunghezze di aste / nodi | `cm` | Dividi per `100` → **`m`** | `150 cm` → `1.5` |
| Dimensioni sezione `b`, `h` | `cm` o `mm` | Dividi per `100` (se cm) o `1000` (se mm) → **`m`** | `b = 20 cm` → `0.2` |
| Forze concentrate o carichi `F`, `q` | `kN` o `kN/m` | Moltiplica per `1000` → **`N`** o **`N/m`** | `12 kN` → `12000` |
| Coppie / Momenti concentrati `M` | `kN·m` | Moltiplica per `1000` → **`N·m`** | `5 kN·m` → `5000` |
| Modulo elastico `E` | `GPa` | Aggiungi 9 zeri (oppure `e9`) → **`Pa`** | `210 GPa` → `2.1e11` (oppure `210000000000`) |
| Modulo elastico `E` | `MPa` (o `N/mm²`) | Aggiungi 6 zeri (oppure `e6`) → **`Pa`** | `30000 MPa` → `3e10` |
| Tensione ammissibile `σ_amm` | `MPa` (o `N/mm²`) | Aggiungi 6 zeri (oppure `e6`) → **`Pa`** | `240 MPa` → `2.4e8` (oppure `240000000`) |
| Cedimenti lineari (`Δx`, `Δy`, `Δn`) | `mm` o `cm` | Il campo vuole già **`mm`** (`1 cm = 10 mm`) | Abbassamento `2 cm` → `-20` |
| Rotazione imposta (`Δγ`) | `gradi (°)` o `rad` | Il campo vuole **`rad`** (`gradi × π / 180`) | `1°` → `0.017453` |
| Radici nel testo (es. `700√2`) | `√2 ≈ 1.414214` | Scrivi il valore decimale con almeno 4-6 cifre | `700√2` → `989.949494` |

> [!TIP]
> Nei campi numerici puoi scrivere direttamente la notazione scientifica con la lettera `e`: ad esempio `2.1e11` significa $2{,}1 \cdot 10^{11}$ e `1.7e-5` significa $1{,}7 \cdot 10^{-5}$.

---

### Passo per passo nella scheda `Struttura` (Sezioni 1 → 8)

#### `1. Nodi` — I punti della figura
I **nodi** sono semplicemente i punti caratteristici del disegno: gli spigoli, i punti dove la trave cambia direzione, i punti dove ci sono gli appoggi o le forze/cerniere.
1. Scegli un punto comodo come origine `(0, 0)` (di solito l'angolo o l'appoggio più in basso a sinistra).
2. Usa **`Per coordinate:`** (`X (m)` verso destra, `Y (m)` verso l'alto) oppure **`Da un nodo:`** (`L (m)` e `Ang (°)`) se sul disegno hai una trave inclinata di cui sai lunghezza e angolo rispetto all'orizzontale.
3. Controlla nell'**`Anteprima Live`** che i pallini dei nodi siano nella posizione giusta.

#### `2. Aste` — I tratti dritti che collegano i nodi
Le **aste** sono i segmenti dritti della struttura.
1. Scegli **`Da nodo:`** e **`A nodo:`** per collegare due nodi.
2. **ATTENZIONE ALL'ORDINE!** L'ordine in cui selezioni `Da nodo -> A nodo` stabilisce:
   - Dove inizia (`s = 0`) e dove finisce (`s = L`) l'asta nei calcoli;
   - Qual è il lato **"destro"** e il lato **"sinistro"** dell'asta (immagina di camminare in piedi lungo l'asta partendo da `Da nodo` guardando verso `A nodo`).

> [!IMPORTANT]
> Se una trave è lunga 4 metri e a metà c'è un cambio di direzione, un appoggio o una cerniera interna, inserisci un nodo nel punto intermedio e crea due aste separate. Anche per forze isolate applicate a metà asta serve un nodo intermedio; per un carico distribuito solo su una parte dell'asta puoi invece usare il menu **`Tratto`** nella sezione `5`.

#### `3. Vincoli` — Gli appoggi fissati al muro o al pavimento
Un **vincolo** è un simbolo disegnato sul nodo che tiene ferma la struttura (un triangolino con le rotelle, un tratteggio da muro, una biella, ecc.). Scegli il nodo e abbina il simbolo sul foglio d'esame alla voce esatta nel menu:

| Simbolo / Descrizione sul foglio | Nome da scegliere nel programma | Che cosa fa (in parole povere) |
| :--- | :--- | :--- |
| Tratteggio fisso direttamente sul muro (senza cerchietto) | **`Incastro`** | Blocca tutto: non sposta in X, non sposta in Y, non ruota |
| Triangolo fissato a terra (senza rotelline sotto) | **`Cerniera`** | Non si sposta né a destra né in su, ma può ruotare |
| Triangolo con **rotelline orizzontali** (scorre a destra/sinistra) | **`Carrello X (Scorre Orizz.)`** | Scorre in orizzontale, blocca il movimento verticale (reazione `Ry`) |
| Triangolo con **rotelline verticali** (su un muro verticale) | **`Carrello Y (Scorre Vert.)`** | Scorre in verticale, blocca il movimento orizzontale (reazione `Rx`) |
| Triangolo con rotelline **su un piano obliquo** | **`Carrello Inclinato`** | Scorre lungo il piano inclinato di `Angolo(°)` rispetto all'asse orizzontale `+x`; reazione `Rn` perpendicolare al piano (`angolo + 90°`) |
| Due sbarrette parallele orizzontali (scorre orizzontalmente ma non ruota) | **`Bipendolo (Scorre Orizz.)`** | Scorre in orizzontale X, blocca verticale Y e rotazione |
| Due sbarrette parallele verticali (scorre verticalmente ma non ruota) | **`Bipendolo (Scorre Vert.)`** | Scorre in verticale Y, blocca orizzontale X e rotazione |
| Un'asta singola collegata a terra con due cerniere alle estremità | **`Pendolo / biella inclinata (tirante o puntone)`** | **`Angolo(°):`** direzione dell'asta misurata da `+x` andando **dal nodo della struttura verso terra**. Lascia **`Lungh. biella (m):` = `0`** se il testo non la dà (`0 = n.i.`). Nel risultato, `Rn > 0` significa tirante (in trazione) |
| Doppio pendolo inclinato su direzione obliqua | **`Doppio pendolo inclinato (bipendolo)`** | Come il pendolo per l'angolo d'asse (dal nodo verso terra), ma blocca anche la rotazione |

#### `4. Carichi sui nodi (forze e coppie)` — Frecce applicate in un punto
Qui inserisci le forze isolate (`Fx`, `Fy`), le forze oblique (`F inclinata` + `angolo (°)`) e i momenti/coppie concentrate (`M (N·m)` — le frecce curve disegnate su un nodo):
- Per **`F inclinata (N)`**, `angolo (°)` si misura **in senso antiorario** partendo da destra (`+x`):
  - `0°` = freccia verso **destra**
  - `90°` = freccia verso **l'alto**
  - `180°` = freccia verso **sinistra**
  - `-90°` (o `270°`) = freccia verso **il basso**
- Per **`M (N·m)`** (e anche per **`Fx`**, **`Fy`**), dopo averlo aggiunto **guarda sempre nell'`Anteprima Live`** se la freccia dritta o curva sul nodo gira/punta nello stesso verso del disegno d'esame: se punta al contrario, toglila e reinseriscila col segno opposto (`-`).

#### `5. Carichi distribuiti sulle aste` — Frecce "a pettine" o "a triangolo" lungo un'asta
Scegli l'asta e compila i campi:
- **`Tipo`**:
  - **`Costante`**: rettangolo di freccette tutte uguali (`q`)
  - **`Triangolare (0 ➔ q)`**: parte da `0` su `Da nodo` e arriva a `q` su `A nodo`
  - **`Triangolare (q ➔ 0)`**: parte da `q` su `Da nodo` e scende a `0` su `A nodo`
- **`Direzione`**:
  - **`Perpendicolare all'asta`**: con `q > 0` spinge verso il lato **DESTRO** di chi cammina da `Da nodo` a `A nodo`. Esempio pratico: su una trave orizzontale disegnata da sinistra a destra, **`q > 0` spinge verso il basso**. Su un pilastro verticale disegnato dal basso verso l'alto, **`q > 0` spinge verso destra**.
  - **`Verticale (+ verso il basso)`**: sempre verticale verso il basso se `q > 0` (indipendentemente da come è inclinata l'asta).
  - **`Orizzontale (+ verso destra)`**: sempre orizzontale verso destra se `q > 0`.
- **`Riferito a`**: usa **`per unità di lunghezza d'asta`** nei casi normali; usa **`per unità di proiezione`** solo se su un'asta inclinata il testo specifica che il carico è riferito alla proiezione orizzontale/verticale.
- **`Tratto`**: puoi applicarlo su **`Tutta l'asta`**, **`Prima metà (da Nodo Inizio)`**, **`Seconda metà (verso Nodo Fine)`**, oppure **`Tratto personalizzato`** impostando **`da`** e **`a`** come frazioni tra `0` e `1` (es. da `0.25` a `0.75`).

*Controlla SEMPRE la direzione delle freccette nell'`Anteprima Live`!*

#### `6. Cerniere interne` — Cerchietti interni che uniscono più aste
Una **cerniera interna** è un pallino/cerchietto dentro la struttura (non sul muro) dove due o più aste si toccano potendo ruotare liberamente una rispetto all'altra.
1. Seleziona il **`Nodo:`** dove si trova la cerniera.
2. **`Aste incernierate:`**
   - Se in quel nodo arrivano solo **2 aste**, puoi lasciare il campo **vuoto** (fa in automatico).
   - Se nel nodo arrivano **3 o più aste**, **DEVI scrivere** gli ID delle aste effettivamente incernierate separati da virgola (es. `2` oppure `2,3`).
3. **`Carichi su asta:`** lascia `auto` di default, oppure indica l'ID dell'asta che "regge" eventuali forze o vincoli applicati proprio su quel nodo.

#### `7. Cedimenti vincolari / rotazioni imposte` — Appoggi che cedono o ruotano
Se il testo dice che un appoggio "cede" / "si sposta" / "ruota di un angolo":
- **`Δx (mm)`**: positivo se l'appoggio si sposta verso **destra** (`+x`).
- **`Δy (mm, ↑ +)`**: positivo verso l'**alto** (`+y`). **ATTENZIONE:** quasi sempre all'esame un cedimento verticale è un *abbassamento* verso il basso, quindi va scritto **con il segno MENO** (esempio: *"cedimento verticale verso il basso di 15 mm"* → scrivi **`-15`**!).
- **`Δn (mm, car. inclinato)`**: cedimento per il carrello inclinato, positivo nel verso `angolo + 90°`.
- **`Δγ (rad, antiorario +)`**: rotazione imposta in **radianti**, positiva in senso **antiorario** (se il testo dà un cedimento rotazionale orario, metti il segno **`-`**!).

#### `8. Variazione termica (gradiente nello spessore)` — Riscaldamento/raffreddamento (`ΔT`)
Se il testo indica temperature diverse sulle due facce di un'asta (es. `+ΔT` sopra e `-ΔT` sotto, oppure un salto termico nello spessore):
1. Prima di tutto inserisci il coefficiente di dilatazione termica **`α (1/°C)`** (es. `1.7e-5`).
2. Ricorda che **"sinistra"** e **"destra"** sono sempre rispetto a chi cammina lungo l'asta da **`Da nodo -> A nodo`** (su un'asta orizzontale da sinistra a destra: **sinistra = lato superiore**, **destra = lato inferiore**).
3. Scegli la modalità che corrisponde ai dati del testo:
   - **`Due fibre (sinistra / destra)`**: quando il testo ti dà i due numeri separati delle due facce.
   - **`Farfalla ±ΔT (sinistra +ΔT, destra −ΔT)`**: quando hai un diagramma termico simmetrico indicato con `+ΔT` su una faccia e `−ΔT` sull'altra.
   - **`Differenza totale ΔT (sinistra +ΔT/2, destra −ΔT/2)`**: quando il testo ti dà solo il salto termico complessivo $\Delta T = T_{\text{sin}} - T_{\text{des}}$.

---

### Passo nella scheda `Materiale e sezione`

Qui imposti i dati numerici della trave e scegli cosa ti chiede il testo ai Punti 3 e 4:
1. **`E (Pa)`**: modulo elastico in Pascal (puoi arrivare fino a `1e15`).
2. **`Sigma amm (Pa)`**: tensione ammissibile del materiale in Pascal (puoi arrivare fino a `1e12`).
3. **`Base b (m)`** e **`Altezza h (m)`**: misure della sezione **assegnata** dal testo (in metri!).
4. **`MS verifica (punto 3)`**: il margine di sicurezza richiesto per la verifica iniziale (se non richiesto, controlla il valore inserito).
5. **`Criterio`**:
   - **`Presso/Tensoflessione (N/A + M/W)`**: formula completa (tiene conto sia della forza assiale `N` sia del momento flettente `M`, standard per i telai piani).
   - **`Solo Flessione (M/W)`**: usa solo `M/W`, solo se il testo dice esplicitamente di trascurare lo sforzo normale `N`.
6. **`EI assegnato (N·m²)`**: lascia **`0`** se vuoi che il programma calcoli da solo $EI = E \cdot \frac{b \cdot h^3}{12}$; inserisci un numero solo se il testo ti dà direttamente il prodotto `EI`.
7. **`Dimensionamento`** (per il Punto 4):
   - **`Nessuno (solo verifica)`**: fa solo la verifica (Punto 3).
   - **`Incognita h (b assegnata)`**: fissa `b` e calcola `h` minima.
   - **`Incognita b (h assegnata)`**: fissa `h` e calcola `b` minima.
   - **`Rapporto h/b assegnato (incognita b)`**: imposta il rapporto scelto in **`Rapporto h/b`** (es. `2` se `h = 2b`) e trova `b` (e quindi `h`).
8. **`MS richiesto (punto 4)`**: il margine di sicurezza desiderato nel dimensionamento del Punto 4.
9. **`Sezione del Pendolo (biella)`**: da compilare solo se è presente una biella/pendolo con una propria sezione specificata.

> [!TIP]
> Appena hai finito di inserire tutti i dati, apri la scheda **`Progetto`** e clicca su **`Salva JSON`** per salvare il file dell'esercizio in una tua cartella! Così se sbagli qualcosa o vuoi riprovare un altro valore, non devi riscrivere tutto da capo.

---

## Parte 3 — Avviare il calcolo, risolvere blocchi e copiare i dati sul foglio

### 1. Scegliere l'incognita `X₁` e premere Risolvi
1. In basso, apri il menu **`Incognita iperstatica X₁:`** e scegli una delle reazioni proposte (una vale l'altra dal punto di vista del risultato finale, ma alcune rendono i conti sul foglio molto più corti: in genere conviene togliere il vincolo più semplice come un carrello o un pendolo).
2. Premi **`▶ Risolvi struttura`**.
3. **Se compare l'errore *"Analisi Interrotta — Struttura Labile"*:** non è un errore del programma! Significa semplicemente che togliendo proprio *quella* reazione la struttura cadrebbe (diventerebbe *labile*). **Soluzione:** scegli una voce diversa nel menu **`Incognita iperstatica X₁:`** in basso e premi di nuovo **`▶ Risolvi struttura`**.
4. Quando il calcolo finisce senza errori, clicca in alto a destra sulla scheda **`Report di Calcolo`** per vedere tutto lo svolgimento completo. Se vuoi conservarlo o stamparlo, clicca in basso su **`🖨️ PDF`**.

---

### 2. Mappa esatta: cosa copiare da ogni sezione del `Report di Calcolo`

Il report è già scritto esattamente nello stesso ordine di uno svolgimento d'esame universitario:

| Sezione nel `Report di Calcolo` | Cosa trovi dentro e come copiarla sul foglio d'esame |
| :--- | :--- |
| **Box `Verifica Cinematica`** | Copia il conteggio dei gradi di libertà (`3n` per i tratti rigidi), la somma dei vincoli (`v`), l'analisi dei vincoli esterni e interni e il risultato (`h = 1`, struttura 1 volta iperstatica). |
| **`1. Dati della struttura...`** | Ricontrolla subito che unità, geometria, carichi, eventuali cedimenti (`Δx`, `Δy`, `Δn`, `Δγ`), effetti termici (`ΔT`, `α`), dati della sezione assegnata (`b`, `h`, `A`, `I`, `W`, `EI`, `E`, `σ_amm`) e parametro del Punto 4 (`h/b = ...`) combacino con il testo dell'esame. |
| **`2. Scelta dell'incognita iperstatica...`** | Copia lo schema isostatico equivalente con l'incognita `X₁` messa in evidenza. **Fai attenzione ai colori se `X₁` è su una cerniera interna:** nel disegno, le frecce **blu** rappresentano l'azione sul tronco precedente, mentre le frecce **viola** sono uguali e opposte sul tronco successivo. |
| **`3. Sistema S₀`** | Copia:<br>• I **3 schemi** del sistema `S₀` (con i soli carichi esterni e `X₁ = 0`);<br>• La tabella delle **reazioni vincolari**;<br>• L'**equilibrio corpo per corpo** (ogni blocco **`TRATTO`** con le azioni alle estremità e le equazioni $\sum F_x = 0$, $\sum F_y = 0$, $\sum M = 0$);<br>• Le espressioni di $N_0(s)$, $T_0(s)$, $M_0(s)$ tratto per tratto;<br>• I grafici/diagrammi globali $N_0$, $T_0$, $M_0$. |
| **`4. Sistema S₁`** | Copia lo schema scarico con applicata **solo `X₁ = 1`**, le relative reazioni, le funzioni $N_1(s)$, $T_1(s)$, $M_1(s)$ e i relativi diagrammi globali $N_1$, $T_1$, $M_1$. |
| **`5. Coefficienti di deformabilità (PLV)`** | Copia la tabella/svolgimento degli **integrali di Mohr** asta per asta per calcolare $\delta_{10}$ e $\delta_{11}$ (con i segni e il denominatore `EI` già messi in ordine). |
| **`6. Equazione di congruenza e calcolo di X₁`** | Copia l'equazione di congruenza $\delta_{10} + \delta_{11} X_1 + \text{eventuali termini dovuti a cedimenti o }\Delta T = 0$ e il valore numerico finale di **`X₁`**. Se $X_1$ esce negativo sul foglio, significa semplicemente che il verso reale è opposto a quello ipotizzato nello schema `S₁`. |
| **`7. Sistema reale`** | Qui c'è il risultato finale completo della struttura:<br>• **Reazioni vincolari reali** e **forze nelle cerniere interne**;<br>• Espressioni reali $N(s)$, $T(s)$, $M(s)$ per ogni asta (con la sovrapposizione $N = N_0 + X_1 N_1$, ecc.);<br>• **Punti stazionari / estremi (`s*` e `M_extr`)**: controlla sempre i valori di massimo/minimo riportati per i tratti con andamento parabolico del momento!<br>• **Diagrammi finali (`N`, `T`, `M`)**: **REGOLA FONDAMENTALE DI SCIENZA DELLE COSTRUZIONI:** il diagramma del momento **`M`** è disegnato **dal lato delle fibre tese** (cioè dalla parte della trave che "tira" / si allunga). Quando copi il diagramma sul foglio, **copia la forma esattamente dal lato mostrato dal grafico del report** (sopra o sotto la linea della trave)! |
| **`8. ... (Punto 3 — Navier)`** | Anche se stai facendo il dimensionamento del Punto 4, questa sezione mostra sempre tutta la verifica della sezione iniziale assegnata:<br>• Tabella delle sezioni candidate più sollecitate;<br>• Formula della tensione massima $\sigma_{\max} = \frac{\lvert N \rvert}{A} + \frac{\lvert M \rvert}{W}$ (oppure solo $\frac{\lvert M \rvert}{W}$ se hai scelto Solo Flessione);<br>• Margine di sicurezza $MS = \frac{\sigma_{\text{amm}}}{\sigma_{\max}} - 1$;<br>• Esito finale **`VERIFICATA`** oppure **`NON VERIFICATA`** rispetto a `MS verifica (punto 3)`. |
| **`9. Dimensionamento della sezione (Punto 4)`** | Appare subito dopo il Punto 3 quando è attiva una modalità di dimensionamento:<br>• Copia la disequazione di progetto impostata in funzione di `b` e `h` (nota: se hai **cedimenti** o **variazioni termiche**, il report mostra esplicitamente anche la dipendenza $X_1(b,h)$ perché cambiando `b` e `h` cambia `EI` e quindi cambiano le sollecitazioni!);<br>• Copia la soluzione numerica (**`b_min`**, **`h_min`** oppure l'intervallo **`[b_min, b_max]`**);<br>• **Attenzione all'avviso critico:** se cambia la sezione più sollecitata al variare delle dimensioni, il report segnala esplicitamente il **cambio della sezione critica** — riportalo sul foglio insieme al controllo finale $MS = MS_{\text{req}}$. |

> [!TIP]
> **Come ricopiare i numeri senza fare errori:** nei passaggi intermedi (es. $\delta_{10}$ e $\delta_{11}$) tieni almeno 3-4 cifre significative per evitare errori di arrotondamento quando fai la verifica con la calcolatrice; per i risultati finali segui le cifre richieste dal docente.

> [!WARNING]
> **Limiti del programma da ricordare sempre all'esame:**
> 1. Risolve **strutture piane** con **1 solo grado di iperstaticità (`h = 1`)** (non è adatto a strutture isostatiche `h = 0` né due o più volte iperstatiche `h ≥ 2`).
> 2. Gli integrali di Mohr in `5.` considerano la **deformabilità flessionale (`EI`)** (non includono eventuale deformabilità assiale `EA` o da taglio `GA` delle travi).

---

## Parte 4 — Checklist rapida di controllo (30 secondi prima di copiare)

1. L'**`Anteprima Live`** mostra esattamente la stessa forma, gli stessi vincoli, le stesse cerniere e le stesse frecce nel verso giusto (su/giù/destra/sinistra) del foglio d'esame?
2. Hai convertito **`cm → m`**, **`kN → N`**, **`GPa/MPa → Pa`**?
3. Se c'è un cedimento verticale verso il basso, hai messo il **segno `-`** a **`Δy (mm, ↑ +)`**?
4. Se c'è variazione termica, hai inserito **`α (1/°C)`** e verificato lato sinistro/destro guardando nel verso **`Da nodo -> A nodo`**?
5. Hai controllato che **`MS verifica (punto 3)`** e **`MS richiesto (punto 4)`** coincidano con quelli richiesti dalle rispettive domande del compito?
6. Hai salvato il file con **Progetto → Salva JSON**?

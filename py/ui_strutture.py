import sys
import os
import json
import traceback

# --- DEFINIZIONE DELLA CARTELLA BASE SICURA ---
if getattr(sys, 'frozen', False):
    base_dir = os.path.dirname(sys.executable)
else:
    base_dir = os.path.dirname(os.path.abspath(__file__))

# Imposta la cartella di lavoro corrente
os.chdir(base_dir)

# --- COMPATIBILITÀ PYTHON 3.9 PER 'Self' ---
import typing
try:
    from typing_extensions import Self
    if not hasattr(typing, 'Self'):
        setattr(typing, 'Self', Self)
except ImportError:
    pass

from PyQt6.QtCore import Qt, QThread, pyqtSignal, QUrl, QMarginsF
from PyQt6.QtGui import QFont, QIcon, QPageLayout, QPageSize
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QPushButton, QTextEdit, QLabel, QFileDialog, QMessageBox, 
    QComboBox, QSplitter, QTabWidget, QTextBrowser
)
# IMPORTANTE: Aggiunto modulo per la stampa PDF
from PyQt6.QtPrintSupport import QPrinter

try:
    from solutore_plv import risolvi_struttura_plv, normalizza_dati, valida_dati
except Exception as e:
    print(f"ERRORE CRITICO NELL'IMPORT DI solutore_plv: {e}")
    normalizza_dati = valida_dati = None

import math
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QPushButton, QTextEdit, QLabel, QFileDialog, QMessageBox, 
    QComboBox, QSplitter, QTabWidget, QTextBrowser,
    QGroupBox, QFormLayout, QDoubleSpinBox, QListWidget, QScrollArea, QGridLayout, QLineEdit,
    QDialog, QDialogButtonBox, QListWidgetItem, QToolButton, QFrame, QSizePolicy
)
import copy
import html
import re
import shutil
import time

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.markers import MarkerStyle
from matplotlib.patches import Polygon as PoligonoMpl, Circle as CerchioMpl

# Tipi di vincolo a terra selezionabili (l'ordine e' quello dei menu a tendina; _costruisci_vincolo/_tipo_da_vincolo li traducono)
TIPI_VINCOLO = [
    "Incastro",
    "Cerniera",
    "Carrello X (Scorre Orizz.)",
    "Carrello Y (Scorre Vert.)",
    "Carrello Inclinato",
    "Bipendolo (Scorre Orizz.)",
    "Bipendolo (Scorre Vert.)",
    "Pendolo / biella inclinata (tirante o puntone)",
    "Doppio pendolo inclinato (bipendolo)",
]
IDX_PENDOLO, IDX_BIPENDOLO_INCL = 7, 8


def costruisci_vincolo(idx_tipo, nodo, ang, lung=0.0):
    """Dizionario 'supports' per il tipo scelto. Per pendolo e doppio pendolo `ang` e' la direzione dell'ASSE delle bielle
    (da +x, antioraria, orientata dal nodo verso il vincolo a terra); per il carrello inclinato e' l'angolo della superficie."""
    s = {"node_id": nodo, "fix_u": False, "fix_v": False, "fix_w": False, "angle": ang}
    if idx_tipo == 0:   s.update({"fix_u": True, "fix_v": True, "fix_w": True})
    elif idx_tipo == 1: s.update({"fix_u": True, "fix_v": True})
    elif idx_tipo == 2: s.update({"fix_v": True})
    elif idx_tipo == 3: s.update({"fix_u": True})
    elif idx_tipo == 4: s.update({"rotated_roller": True, "angle": ang})
    elif idx_tipo == 5: s.update({"fix_v": True, "fix_w": True})
    elif idx_tipo == 6: s.update({"fix_u": True, "fix_w": True})
    elif idx_tipo in (IDX_PENDOLO, IDX_BIPENDOLO_INCL):
        s.pop("angle", None)
        s.update({"axis_angle": ang, "rotated_roller": True, "angle": ang - 90})
        s["pendolo" if idx_tipo == IDX_PENDOLO else "bipendolo"] = True
        if idx_tipo == IDX_BIPENDOLO_INCL: s["fix_w"] = True
        if lung and lung > 0:
            s["pend_dx"] = float(f"{lung * math.cos(math.radians(ang)):.12g}")
            s["pend_dy"] = float(f"{lung * math.sin(math.radians(ang)):.12g}")
    return s


def tipo_da_vincolo(s):
    """Indice in TIPI_VINCOLO di un vincolo esistente (None se non e' riconducibile a un tipo del form)."""
    if s.get("pendolo"): return IDX_PENDOLO
    if s.get("bipendolo"): return IDX_BIPENDOLO_INCL
    u, v, w = bool(s.get("fix_u")), bool(s.get("fix_v")), bool(s.get("fix_w"))
    if s.get("rotated_roller"): return 4 if not w else None
    return {(True, True, True): 0, (True, True, False): 1, (False, True, False): 2, (True, False, False): 3,
            (False, True, True): 5, (True, False, True): 6}.get((u, v, w))


def angolo_vincolo(s):
    """Angolo mostrato/modificato nel form: asse della biella per pendolo/bipendolo, angolo del carrello altrimenti."""
    if s.get("pendolo") or s.get("bipendolo"):
        return float(s.get("axis_angle", float(s.get("angle", 0)) + 90))
    return float(s.get("angle", 0))


class Sezione(QWidget):
    """Riquadro comprimibile: l'intestazione e' un pulsante che apre/chiude il contenuto.
    Espone setLayout() come un QGroupBox, quindi il codice dei form resta invariato."""

    def __init__(self, titolo, aperta=False, parent=None):
        super().__init__(parent)
        self._titolo = titolo
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        self._head = QToolButton()
        self._head.setObjectName("sezHeader")
        self._head.setCheckable(True)
        self._head.setChecked(aperta)
        self._head.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self._head.setStyleSheet("text-align: left;")
        self._head.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._head.toggled.connect(self._aggiorna)
        self._body = QFrame()
        self._body.setObjectName("sezBody")
        outer.addWidget(self._head)
        outer.addWidget(self._body)
        self._aggiorna(aperta)

    def setLayout(self, layout):
        layout.setContentsMargins(8, 8, 8, 8)
        self._body.setLayout(layout)

    def _aggiorna(self, aperta):
        self._head.setText(("▾  " if aperta else "▸  ") + self._titolo)
        self._body.setVisible(aperta)

    def apri(self, aperta=True):
        self._head.setChecked(aperta)


class SolverThread(QThread):
    finished_signal = pyqtSignal(str, str, str)
    error_signal = pyqtSignal(str)

    def __init__(self, dati_json, out_dir, incognita):
        super().__init__()
        self.dati_json = dati_json
        self.out_dir = out_dir
        self.incognita = incognita

    def run(self):
        try:
            report_html, img_s1, img_diag = risolvi_struttura_plv(self.dati_json, self.out_dir, self.incognita)
            self.finished_signal.emit(report_html, img_s1, img_diag)
        except Exception as e:
            error_details = traceback.format_exc()
            self.error_signal.emit(error_details)

class UIStrutture(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Risolutore PLV - Metodo delle forze")
        # Finestra dimensionata sullo schermo (92% x 90%), centrata: il pannello dati ha subito spazio e la figura resta visibile
        _geo = QApplication.primaryScreen().availableGeometry()
        _w, _h = int(_geo.width() * 0.92), int(_geo.height() * 0.90)
        self.resize(_w, _h)
        self.move(_geo.x() + (_geo.width() - _w) // 2, _geo.y() + (_geo.height() - _h) // 2)
        self.setMinimumSize(1000, 600)
        self.output_dir = os.path.abspath("temp_output")
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Struttura dati in memoria per il Costruttore Interattivo
        self.builder_data = {
            "nodes": [],
            "elements": [],
            "supports": [],
            "point_loads": [],
            "element_loads": [],
            "thermal_loads": [],
            "material": {"E": 2.1e11, "sigma_amm": 235e6},
            "section": {"b": 0.2, "h": 0.4}
        }
        self.node_id_counter = 1
        self.elem_id_counter = 1
        self._sync_in_corso = False   # True mentre JSON e pannello si aggiornano a vicenda (evita cicli)
        self._undo = []               # istantanee JSON per "Annulla ultima modifica"
        self._sel = None              # (categoria, indice) dell'elemento selezionato nel riepilogo (evidenziato in anteprima)
        self._run_id = 0              # ogni calcolo scrive in una cartella nuova: niente immagini vecchie dalla cache
        
        self.setup_ui()
        self._collega_segnali_proprieta()
        self.setStyleSheet(self.get_stylesheet())

    def get_stylesheet(self):
        return """
            /* Sfondo generale chiaro per evitare interferenze con la Dark Mode */
            QMainWindow, QSplitter { background-color: #f8f9fa; }
            
            /* Forza testo nero e sfondo bianco di base per i contenitori */
            QWidget { color: #000000; }
            QTabWidget::pane, QScrollArea, QScrollArea > QWidget > QWidget { 
                background-color: #ffffff; 
                border: 1px solid #ced4da; 
                border-radius: 4px; 
            }
            
            /* Stile Pulsanti */
            QPushButton { background-color: #0d6efd; color: #ffffff; border: none; padding: 7px 12px; border-radius: 6px; font-weight: bold; font-size: 13px;}
            QPushButton:hover { background-color: #0b5ed7; }
            QPushButton:disabled { background-color: #cccccc; color: #666666; }
            
            #btnSolve { background-color: #198754; font-size: 15px; padding: 11px; color: #ffffff; }
            #btnSolve:hover { background-color: #157347; }
            
            #btnPdf { background-color: #6c757d; font-size: 13px; padding: 11px; color: #ffffff; }
            #btnPdf:hover { background-color: #5c636a; }
            
            #btnAggiungi { background-color: #198754; padding: 5px 10px; font-size: 13px;}
            #btnSecondario { background-color: #e9ecef; color: #212529; font-weight: normal; }
            #btnSecondario:hover { background-color: #dee2e6; }
            #btnPericolo { background-color: #ffffff; color: #dc3545; border: 1px solid #dc3545; }
            #btnPericolo:hover { background-color: #f8d7da; }

            /* Sezioni comprimibili del tab Struttura */
            QToolButton#sezHeader { background-color: #e9ecef; color: #212529; border: 1px solid #ced4da; border-radius: 5px;
                                    padding: 7px 10px; text-align: left; font-weight: bold; font-size: 13px; }
            QToolButton#sezHeader:checked { background-color: #dbe7ff; color: #0a58ca; border-color: #9ec5fe;
                                            border-bottom-left-radius: 0; border-bottom-right-radius: 0; }
            QFrame#sezBody { background-color: #fdfdfd; border: 1px solid #ced4da; border-top: none;
                             border-bottom-left-radius: 5px; border-bottom-right-radius: 5px; }
            #barraCalcolo { background-color: #f1f3f5; border: 1px solid #ced4da; border-radius: 6px; }
            
            /* Stile Campi di testo e Tabelle */
            QTextEdit, QTextBrowser, QListWidget, QLineEdit { background-color: #ffffff; border: 1px solid #ced4da; border-radius: 4px; padding: 10px; color: #000000; }
            QLabel { font-weight: bold; color: #000000; background-color: transparent; }
            
            /* Stile dei riquadri (Box) della schermata Disegna */
            QGroupBox { font-weight: bold; border: 1px solid #ced4da; border-radius: 6px; margin-top: 15px; padding-top: 20px; background-color: #fdfdfd;}
            QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px 0 5px; background-color: transparent; }
            
            QComboBox, QDoubleSpinBox { border: 1px solid #ced4da; border-radius: 4px; padding: 6px; background-color: #ffffff; color: #000000; }
            /* 🟢 REGOLA SALVA-VITA PER I MENU A TENDINA (Lista a discesa) */
            QComboBox QAbstractItemView {
                background-color: #ffffff;
                color: #000000;
                selection-background-color: #0d6efd;
                selection-color: #ffffff;
            }
            
            /* Stile Schede (Tab) */
            QTabBar::tab { background: #e9ecef; color: #000000; padding: 8px 15px; border-top-left-radius: 4px; border-top-right-radius: 4px;}
            QTabBar::tab:selected { background: #ffffff; font-weight: bold; color: #0d6efd; border-bottom: 2px solid #ffffff;}
        """

    def setup_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)
        
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)   # nessun pannello puo' sparire trascinando la maniglia
        splitter.setHandleWidth(8)
        main_layout.addWidget(splitter)
        self.splitter = splitter

        left_panel = QWidget()
        left_panel.setMinimumWidth(560)
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 10, 0)
        
        self.tabs = QTabWidget()
        self.tab_disegna = QWidget()
        self.tab_props = QWidget()
        self.tab_json = QWidget()
        
        self.tabs.addTab(self.tab_disegna, "✏️ Struttura")
        self.tabs.addTab(self.tab_props, "🔩 Materiale e sezione")
        self.tabs.addTab(self.tab_json, "📄 JSON")
        
        self.setup_disegna_tab() # Inizializza la UI (tab Struttura + tab Materiale e sezione)
        
        # --- TAB JSON ---
        tab_json_layout = QVBoxLayout(self.tab_json)
        self.btn_load = QPushButton("Apri JSON da File")
        self.btn_load.clicked.connect(self.load_json)
        self.text_json = QTextEdit()
        self.text_json.setFont(QFont("Courier", 11))
        self.text_json.textChanged.connect(self.aggiorna_anteprima_json)
        tab_json_layout.addWidget(self.btn_load)
        tab_json_layout.addWidget(self.text_json)

        left_layout.addWidget(self.tabs, 1)

        # --- BARRA DI CALCOLO: sempre visibile sotto i tab (incognita iperstatica + Risolvi + PDF) ---
        barra = QFrame()
        barra.setObjectName("barraCalcolo")
        barra_layout = QVBoxLayout(barra)
        barra_layout.setContentsMargins(8, 8, 8, 8)
        barra_layout.setSpacing(6)
        riga_inc = QHBoxLayout()
        riga_inc.addWidget(QLabel("Incognita iperstatica X₁:"))
        self.combo_incognita = QComboBox()
        self.combo_incognita.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.combo_incognita.setToolTip("Reazione da liberare nel sistema principale. Si aggiorna da sola con i vincoli inseriti.")
        riga_inc.addWidget(self.combo_incognita, 1)
        barra_layout.addLayout(riga_inc)
        riga_az = QHBoxLayout()
        self.btn_solve = QPushButton("▶ Risolvi struttura")
        self.btn_solve.setObjectName("btnSolve")
        self.btn_solve.clicked.connect(self.run_solver)
        self.btn_pdf = QPushButton("🖨️ PDF")
        self.btn_pdf.setObjectName("btnPdf")
        self.btn_pdf.setToolTip("Esporta il report di calcolo in PDF (disponibile dopo il calcolo)")
        self.btn_pdf.clicked.connect(self.esporta_pdf)
        self.btn_pdf.setEnabled(False)
        self.btn_solve.setMinimumHeight(42)
        self.btn_pdf.setMinimumHeight(42)
        riga_az.addWidget(self.btn_solve, 3)
        riga_az.addWidget(self.btn_pdf, 1)
        barra_layout.addLayout(riga_az)
        left_layout.addWidget(barra)
        
        # --- SOSTITUISCI IL VECCHIO right_panel CON QUESTO ---
        right_panel = QWidget()
        right_panel.setMinimumWidth(380)   # la figura resta sempre visibile, anche con il pannello dati molto largo
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(10, 0, 0, 0)
        
        self.right_tabs = QTabWidget()
        
        # 🟢 NUOVA SCHEDA A DESTRA: Anteprima Live
        self.tab_preview = QWidget()
        preview_layout = QVBoxLayout(self.tab_preview)
        
        # Creiamo il canvas grafico di Matplotlib per il disegno in tempo reale
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        preview_layout.addWidget(self.canvas)
        self.right_tabs.addTab(self.tab_preview, "👀 Anteprima Live")
        
        # 🟢 SECONDA SCHEDA: Il Report HTML di Calcolo
        self.tab_report = QWidget()
        report_layout = QVBoxLayout(self.tab_report)
        self.html_view = QTextBrowser()
        self.html_view.setOpenExternalLinks(False)
        self.html_view.document().setDefaultStyleSheet("""
            body { font-family: Arial, sans-serif; font-size: 14px; color: #000000; background-color: #ffffff; line-height: 1.5; }
            h1 { color: #0d6efd; border-bottom: 2px solid #0d6efd; padding-bottom: 5px; font-size: 20px; }
            h2 { color: #dc3545; background-color: #f8d7da; padding: 8px; border-left: 4px solid #dc3545; margin-top: 20px; font-size: 16px;}
            table { width: 100%; border-collapse: collapse; margin-top: 10px; margin-bottom: 15px; background-color: #ffffff; }
            th { background-color: #e9ecef; border: 1px solid #ced4da; padding: 8px; font-weight: bold; color: #000000; font-size: 13px; }
            td { border: 1px solid #ced4da; padding: 8px; text-align: center; color: #000000; font-size: 13px; }
            .formula { font-family: 'Courier New', Courier, monospace; font-size: 14px; font-weight: bold; color: #0b5ed7; background-color: #f1f3f5; padding: 2px 6px; border-radius: 4px; border: 1px solid #dee2e6;}
        """)
        self.html_view.setHtml("<h2>Pronto</h2><p>Costruisci la struttura per vederla nell'Anteprima Live.</p>")
        report_layout.addWidget(self.html_view)
        
        self.right_tabs.addTab(self.tab_report, "📄 Report di Calcolo")
        right_layout.addWidget(self.right_tabs)

        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        # Allo stesso avvio: ~58% ai dati, ~42% alla figura. Se si allarga la finestra, la figura (stretch 1) cresce di piu'.
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        _tot = max(self.width() - 40, 1000)
        splitter.setSizes([int(_tot * 0.58), int(_tot * 0.42)])
        
        # Facciamo un primo disegno vuoto
        self.update_live_preview()

    def setup_disegna_tab(self):
        """Inizializza la UI tabellare/form per la costruzione della struttura."""
        layout = QVBoxLayout(self.tab_disegna)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)
        scroll_layout.setSpacing(6)

        # Tab "Materiale e sezione": stessa struttura (scroll + colonna), riempita in fondo a questo metodo
        layout_p = QVBoxLayout(self.tab_props)
        props_scroll = QScrollArea()
        props_scroll.setWidgetResizable(True)
        props_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        props_content = QWidget()
        props_layout = QVBoxLayout(props_content)
        props_scroll.setWidget(props_content)
        layout_p.addWidget(props_scroll)

        # === GRUPPO 1: AGGIUNGI NODI ===
        grp_nodi = Sezione("1. Nodi", aperta=True)
        layout_nodi = QVBoxLayout()
        form_cart = QGridLayout()
        self.spin_x = QDoubleSpinBox(); self.spin_x.setRange(-9999, 9999); self.spin_x.setDecimals(3)
        self.spin_y = QDoubleSpinBox(); self.spin_y.setRange(-9999, 9999); self.spin_y.setDecimals(3)
        btn_add_cart = QPushButton("+ Nodo")
        btn_add_cart.setObjectName("btnAggiungi")
        btn_add_cart.clicked.connect(self.add_nodo_cartesiano)
        form_cart.addWidget(QLabel("Per coordinate:"), 0, 0)
        form_cart.addWidget(QLabel("X (m):"), 0, 1); form_cart.addWidget(self.spin_x, 0, 2)
        form_cart.addWidget(QLabel("Y (m):"), 0, 3); form_cart.addWidget(self.spin_y, 0, 4)
        form_cart.addWidget(btn_add_cart, 0, 5)
        form_cart.setColumnStretch(2, 1); form_cart.setColumnStretch(4, 1)

        form_pol = QGridLayout()
        self.combo_rif_nodo = QComboBox() 
        self.spin_L = QDoubleSpinBox(); self.spin_L.setRange(0, 9999); self.spin_L.setDecimals(3)
        self.spin_ang = QDoubleSpinBox(); self.spin_ang.setRange(-360, 360); self.spin_ang.setDecimals(1)
        btn_add_pol = QPushButton("+ Nodo")
        btn_add_pol.setObjectName("btnAggiungi")
        btn_add_pol.clicked.connect(self.add_nodo_polare)
        form_pol.addWidget(QLabel("Da un nodo:"), 0, 0)
        form_pol.addWidget(self.combo_rif_nodo, 0, 1, 1, 4)
        form_pol.addWidget(QLabel("L (m):"), 1, 1); form_pol.addWidget(self.spin_L, 1, 2)
        form_pol.addWidget(QLabel("Ang (°):"), 1, 3); form_pol.addWidget(self.spin_ang, 1, 4)
        form_pol.addWidget(btn_add_pol, 1, 5)
        form_pol.setColumnStretch(2, 1); form_pol.setColumnStretch(4, 1)

        layout_nodi.addLayout(form_cart)
        layout_nodi.addLayout(form_pol)
        grp_nodi.setLayout(layout_nodi)
        scroll_layout.addWidget(grp_nodi)

        # === GRUPPO 2: AGGIUNGI ASTE ===
        grp_aste = Sezione("2. Aste", aperta=True)
        layout_aste = QGridLayout()
        self.combo_n1 = QComboBox()
        self.combo_n2 = QComboBox()
        btn_add_asta = QPushButton("+ Asta")
        btn_add_asta.setObjectName("btnAggiungi")
        btn_add_asta.clicked.connect(self.add_asta)
        layout_aste.addWidget(QLabel("Da nodo:"), 0, 0); layout_aste.addWidget(self.combo_n1, 0, 1)
        layout_aste.addWidget(QLabel("A nodo:"), 0, 2); layout_aste.addWidget(self.combo_n2, 0, 3)
        layout_aste.addWidget(btn_add_asta, 0, 4)
        layout_aste.setColumnStretch(1, 1); layout_aste.setColumnStretch(3, 1)
        grp_aste.setLayout(layout_aste)
        scroll_layout.addWidget(grp_aste)

        # === GRUPPO 3: AGGIUNGI VINCOLI ===
        grp_vincoli = Sezione("3. Vincoli", aperta=True)
        layout_vinc = QGridLayout()
        self.combo_nodo_vinc = QComboBox()
        self.combo_tipo_vinc = QComboBox()
        self.combo_tipo_vinc.addItems(TIPI_VINCOLO)
        self.combo_tipo_vinc.setToolTip("Pendolo = biella cerniera-cerniera inclinata: reazione solo lungo il suo asse "
                                        "(il solutore dice se lavora da TIRANTE o da PUNTONE e ne dimensiona la sezione).\n"
                                        "Doppio pendolo inclinato (bipendolo) = due bielle parallele: reazione lungo l'asse + momento.\n"
                                        "Piu' pendoli non paralleli sullo stesso nodo = doppio pendolo a bielle concorrenti.")
        self.lbl_ang_vinc = QLabel("Angolo(°):")
        self.spin_ang_vinc = QDoubleSpinBox()
        self.spin_ang_vinc.setRange(-360, 360)
        self.spin_ang_vinc.setToolTip("Gradi rispetto all'orizzontale")
        self.lbl_lung_vinc = QLabel("Lungh. biella (m):")
        self.spin_lung_vinc = QDoubleSpinBox()
        self.spin_lung_vinc.setRange(0, 9999); self.spin_lung_vinc.setDecimals(3)
        self.spin_lung_vinc.setToolTip("Facoltativa (0 = non indicata): serve solo per il disegno e per la verifica a carico di punta.")
        btn_add_vinc = QPushButton("+ Vincolo")
        btn_add_vinc.setObjectName("btnAggiungi")
        btn_add_vinc.clicked.connect(self.add_vincolo)
        self.combo_tipo_vinc.currentIndexChanged.connect(self._aggiorna_etichette_vincolo)

        layout_vinc.addWidget(QLabel("Nodo:"), 0, 0)
        layout_vinc.addWidget(self.combo_nodo_vinc, 0, 1)
        layout_vinc.addWidget(QLabel("Tipo:"), 0, 2)
        layout_vinc.addWidget(self.combo_tipo_vinc, 0, 3)
        layout_vinc.addWidget(self.lbl_ang_vinc, 1, 0)
        layout_vinc.addWidget(self.spin_ang_vinc, 1, 1)
        layout_vinc.addWidget(self.lbl_lung_vinc, 1, 2)
        layout_vinc.addWidget(self.spin_lung_vinc, 1, 3)
        layout_vinc.addWidget(btn_add_vinc, 2, 3)
        grp_vincoli.setLayout(layout_vinc)
        scroll_layout.addWidget(grp_vincoli)
        self._aggiorna_etichette_vincolo()

        # === GRUPPO 4: AGGIUNGI CARICHI CONCENTRATI ===
        grp_carichi = Sezione("4. Carichi sui nodi (forze e coppie)", aperta=True)
        layout_car = QGridLayout()
        self.combo_nodo_carico = QComboBox()
        self.spin_Fx = QDoubleSpinBox(); self.spin_Fx.setRange(-999999, 999999); self.spin_Fx.setDecimals(2)
        self.spin_Fy = QDoubleSpinBox(); self.spin_Fy.setRange(-999999, 999999); self.spin_Fy.setDecimals(2)
        self.spin_M = QDoubleSpinBox(); self.spin_M.setRange(-999999, 999999); self.spin_M.setDecimals(2)
        # Forza inclinata: modulo + angolo, scomposti in Fx, Fy all'inserimento (es. 700*sqrt2 N a 45 gradi)
        self.spin_Finc = QDoubleSpinBox(); self.spin_Finc.setRange(-999999, 999999); self.spin_Finc.setDecimals(6)
        self.spin_Finc_ang = QDoubleSpinBox(); self.spin_Finc_ang.setRange(-360, 360); self.spin_Finc_ang.setDecimals(2)
        self.spin_Finc.setToolTip("Modulo di una forza inclinata (opzionale): viene scomposto in Fx e Fy e sommato alle componenti indicate.\n"
                                  "Per 700·√2 N scrivi 989.949494 (6 decimali): le componenti vengono arrotondate al 1e-4 N (→ 700, 700).")
        self.spin_Finc_ang.setToolTip("Angolo della forza inclinata in gradi, misurato da +x in senso antiorario\n"
                                      "(0 = verso destra, 90 = verso l'alto, 180 = verso sinistra, -90 = verso il basso).")
        btn_add_carico = QPushButton("+ Carico")
        btn_add_carico.setObjectName("btnAggiungi")
        btn_add_carico.clicked.connect(self.add_carico)
        
        # Riga 1: nodo + componenti; riga 2: forza inclinata (modulo, angolo) + pulsante
        layout_car.addWidget(QLabel("Nodo:"), 0, 0);  layout_car.addWidget(self.combo_nodo_carico, 0, 1)
        layout_car.addWidget(QLabel("M (N·m):"), 0, 2); layout_car.addWidget(self.spin_M, 0, 3)
        layout_car.addWidget(QLabel("Fx (N):"), 1, 0); layout_car.addWidget(self.spin_Fx, 1, 1)
        layout_car.addWidget(QLabel("Fy (N):"), 1, 2); layout_car.addWidget(self.spin_Fy, 1, 3)
        layout_car.addWidget(QLabel("F inclinata (N):"), 2, 0); layout_car.addWidget(self.spin_Finc, 2, 1)
        layout_car.addWidget(QLabel("angolo (°):"), 2, 2); layout_car.addWidget(self.spin_Finc_ang, 2, 3)
        layout_car.addWidget(btn_add_carico, 3, 3)
        grp_carichi.setLayout(layout_car)
        scroll_layout.addWidget(grp_carichi)

        # === GRUPPO 5: CARICHI DISTRIBUITI (q) ===
        grp_dist = Sezione("5. Carichi distribuiti sulle aste")
        layout_dist = QGridLayout()
        self.combo_asta_carico = QComboBox()

        self.combo_tipo_q = QComboBox()
        self.combo_tipo_q.addItems(["Costante", "Triangolare (0 ➔ q)", "Triangolare (q ➔ 0)"])

        self.combo_dir_q = QComboBox()
        self.combo_dir_q.addItems(["Perpendicolare all'asta", "Verticale (+ verso il basso)", "Orizzontale (+ verso destra)"])
        self.combo_dir_q.setToolTip("Perpendicolare: q>0 verso la normale DESTRA dell'asta (guardando da Nodo Inizio a Nodo Fine;\n"
                                    "per un'asta orizzontale da sinistra a destra = verso il basso).\n"
                                    "Verticale/Orizzontale: componente globale del carico.")
        self.combo_rif_q = QComboBox()
        self.combo_rif_q.addItems(["per unità di lunghezza d'asta", "per unità di proiezione"])
        self.combo_rif_q.setToolTip("Solo per carichi verticali/orizzontali. 'Proiezione' = N/m riferiti alla lunghezza proiettata\n"
                                    "ortogonalmente alla direzione del carico (es. carico verticale per metro di luce orizzontale).")

        self.spin_q = QDoubleSpinBox(); self.spin_q.setRange(-999999, 999999); self.spin_q.setDecimals(2)
        btn_add_dist = QPushButton("+ Carico Distr.")
        btn_add_dist.setObjectName("btnAggiungi")
        btn_add_dist.clicked.connect(self.add_carico_distribuito)

        layout_dist.addWidget(QLabel("Asta:"), 0, 0); layout_dist.addWidget(self.combo_asta_carico, 0, 1)
        layout_dist.addWidget(QLabel("Tipo:"), 0, 2); layout_dist.addWidget(self.combo_tipo_q, 0, 3)
        layout_dist.addWidget(QLabel("Direzione:"), 1, 0); layout_dist.addWidget(self.combo_dir_q, 1, 1)
        layout_dist.addWidget(QLabel("Riferito a:"), 1, 2); layout_dist.addWidget(self.combo_rif_q, 1, 3)
        layout_dist.addWidget(QLabel("q max:"), 2, 0); layout_dist.addWidget(self.spin_q, 2, 1)

        # Tratto caricato: tutta l'asta, una meta' oppure un tratto qualsiasi (frazioni della lunghezza, 0 = Nodo Inizio)
        self.combo_tratto_q = QComboBox()
        self.combo_tratto_q.addItems(["Tutta l'asta", "Prima metà (da Nodo Inizio)", "Seconda metà (verso Nodo Fine)", "Tratto personalizzato"])
        self.combo_tratto_q.setToolTip("Per caricare solo una parte dell'asta questa viene divisa con un nodo intermedio\n"
                                       "(nessuna cerniera: la continuità è mantenuta) e il carico va sul tratto scelto.")
        self.spin_t0 = QDoubleSpinBox(); self.spin_t0.setRange(0, 1); self.spin_t0.setDecimals(4); self.spin_t0.setSingleStep(0.1); self.spin_t0.setValue(0.0)
        self.spin_t1 = QDoubleSpinBox(); self.spin_t1.setRange(0, 1); self.spin_t1.setDecimals(4); self.spin_t1.setSingleStep(0.1); self.spin_t1.setValue(1.0)
        self.spin_t0.setToolTip("Inizio del tratto, come frazione della lunghezza (0 = Nodo Inizio, 1 = Nodo Fine)")
        self.spin_t1.setToolTip("Fine del tratto, come frazione della lunghezza")
        self.lbl_t0, self.lbl_t1 = QLabel("da:"), QLabel("a:")
        riga_tratto = QHBoxLayout()
        for w_ in (self.combo_tratto_q, self.lbl_t0, self.spin_t0, self.lbl_t1, self.spin_t1):
            riga_tratto.addWidget(w_)
        layout_dist.addWidget(QLabel("Tratto:"), 3, 0); layout_dist.addLayout(riga_tratto, 3, 1, 1, 3)
        self.combo_tratto_q.currentIndexChanged.connect(self._aggiorna_tratto_q)
        self._aggiorna_tratto_q()
        layout_dist.addWidget(btn_add_dist, 4, 3)
        grp_dist.setLayout(layout_dist)
        scroll_layout.addWidget(grp_dist)

        # === GRUPPO 6: CERNIERE INTERNE ===
        grp_cerniera_int = Sezione("6. Cerniere interne")
        layout_cern_int = QGridLayout()
        self.combo_nodo_cerniera_int = QComboBox()
        
        btn_add_cerniera_int = QPushButton("+ Cerniera Interna")
        btn_add_cerniera_int.setObjectName("btnAggiungi")
        btn_add_cerniera_int.clicked.connect(self.add_cerniera_interna)
        
        self.edit_elems_cerniera = QLineEdit()
        self.edit_elems_cerniera.setPlaceholderText("id aste incernierate, es. 2 oppure 2,3 (vuoto = automatico)")
        self.edit_elems_cerniera.setToolTip("Obbligatorio se nel nodo convergono 3 o più aste: indica le aste che sono\n"
                                            "realmente incernierate in quel nodo (le altre restano incastrate fra loro).\n"
                                            "Se le indichi tutte, la cerniera è comune a tutti i corpi.")
        layout_cern_int.addWidget(QLabel("Nodo:"), 0, 0)
        layout_cern_int.addWidget(self.combo_nodo_cerniera_int, 0, 1, 1, 3)
        self.edit_carichi_su = QLineEdit()
        self.edit_carichi_su.setPlaceholderText("auto")
        self.edit_carichi_su.setToolTip("Corpo a cui appartengono coppie concentrate, vincoli (appoggi, incastri) e forze applicati\n"
                                        "nel nodo della cerniera: indica l'id di una delle aste incernierate (il corpo è quello di quell'asta).\n"
                                        "Vuoto = il corpo collegato rigidamente al nodo (o l'ultima asta, se sono tutte incernierate).")
        layout_cern_int.addWidget(QLabel("Aste incernierate:"), 1, 0)
        layout_cern_int.addWidget(self.edit_elems_cerniera, 1, 1, 1, 3)
        layout_cern_int.addWidget(QLabel("Carichi su asta:"), 2, 0)
        layout_cern_int.addWidget(self.edit_carichi_su, 2, 1)
        layout_cern_int.addWidget(btn_add_cerniera_int, 2, 3)
        grp_cerniera_int.setLayout(layout_cern_int)
        scroll_layout.addWidget(grp_cerniera_int)

        # === GRUPPO 7: CEDIMENTI VINCOLARI E ROTAZIONI IMPOSTE ===
        grp_ced = Sezione("7. Cedimenti vincolari / rotazioni imposte")
        layout_ced = QGridLayout()
        self.combo_nodo_ced = QComboBox()
        self.spin_ced_x = QDoubleSpinBox(); self.spin_ced_x.setRange(-9999, 9999); self.spin_ced_x.setDecimals(3)
        self.spin_ced_y = QDoubleSpinBox(); self.spin_ced_y.setRange(-9999, 9999); self.spin_ced_y.setDecimals(3)
        self.spin_ced_n = QDoubleSpinBox(); self.spin_ced_n.setRange(-9999, 9999); self.spin_ced_n.setDecimals(3)
        self.spin_ced_w = QDoubleSpinBox(); self.spin_ced_w.setRange(-10, 10); self.spin_ced_w.setDecimals(6)
        self.spin_ced_x.setToolTip("mm, positivo verso +x")
        self.spin_ced_y.setToolTip("mm, positivo verso +y (cedimento verso il basso = valore negativo)")
        self.spin_ced_n.setToolTip("mm, positivo lungo la normale del carrello inclinato (angolo+90°)")
        self.spin_ced_w.setToolTip("rad, positiva antioraria")
        btn_add_ced = QPushButton("+ Cedimento")
        btn_add_ced.setObjectName("btnAggiungi")
        btn_add_ced.clicked.connect(self.add_cedimento)
        layout_ced.addWidget(QLabel("Nodo (con vincolo):"), 0, 0); layout_ced.addWidget(self.combo_nodo_ced, 0, 1)
        layout_ced.addWidget(QLabel("Δx (mm):"), 1, 0); layout_ced.addWidget(self.spin_ced_x, 1, 1)
        layout_ced.addWidget(QLabel("Δy (mm, ↑ +):"), 1, 2); layout_ced.addWidget(self.spin_ced_y, 1, 3)
        layout_ced.addWidget(QLabel("Δn (mm, car. inclinato):"), 2, 0); layout_ced.addWidget(self.spin_ced_n, 2, 1)
        layout_ced.addWidget(QLabel("Δφ (rad, antiorario +):"), 2, 2); layout_ced.addWidget(self.spin_ced_w, 2, 3)
        layout_ced.addWidget(btn_add_ced, 3, 3)
        grp_ced.setLayout(layout_ced)
        scroll_layout.addWidget(grp_ced)

        # === GRUPPO 8: VARIAZIONI TERMICHE (GRADIENTE) ===
        grp_term = Sezione("8. Variazione termica (gradiente nello spessore)")
        layout_term = QGridLayout()
        self.combo_asta_term = QComboBox()
        self.spin_T_pos = QDoubleSpinBox(); self.spin_T_pos.setRange(-999, 999); self.spin_T_pos.setDecimals(2)
        self.spin_T_neg = QDoubleSpinBox(); self.spin_T_neg.setRange(-999, 999); self.spin_T_neg.setDecimals(2)
        self.spin_T_pos.setToolTip("Variazione di temperatura della fibra a SINISTRA dell'asta (guardando da Nodo Inizio a Nodo Fine)")
        self.spin_T_neg.setToolTip("Variazione di temperatura della fibra a DESTRA dell'asta (guardando da Nodo Inizio a Nodo Fine)")
        self.combo_term_modo = QComboBox()
        self.combo_term_modo.addItems(["Due fibre (sinistra / destra)",
                                       "Farfalla ±ΔT (sinistra +ΔT, destra −ΔT)",
                                       "Differenza totale ΔT (sinistra +ΔT/2, destra −ΔT/2)"])
        self.combo_term_modo.setToolTip("Come leggere il dato della prova:\n"
                                        "• Due fibre: si assegnano le due temperature (T sinistra, T destra), qualsiasi valore.\n"
                                        "• Farfalla ±ΔT: «ΔT = 35 °C» significa +35 su una fibra e −35 sull'altra: differenza fra le fibre = 2ΔT. Usa solo il campo sinistra.\n"
                                        "• Differenza totale ΔT: «ΔT = 35 °C» è la differenza fra le due fibre (±17.5). Usa solo il campo sinistra.\n"
                                        "La curvatura termica è κ = α·(T_sx − T_dx)/h: i tre modi differiscono di un fattore 2 o più.")
        self.combo_term_modo.currentIndexChanged.connect(self._aggiorna_modo_termico)
        self.spin_alpha = QDoubleSpinBox(); self.spin_alpha.setDecimals(8); self.spin_alpha.setRange(0, 1)
        self.spin_alpha.setSingleStep(1e-6)
        self.spin_alpha.setSpecialValueText("— da inserire —")   # 0 = non impostato: va scritto solo se c'e' un carico termico
        self.spin_alpha.setValue(0)
        self.spin_alpha.setToolTip("Coefficiente di dilatazione termica del testo della prova (es. 1.7e-5 1/°C).\n"
                                   "Va inserito solo nelle prove con variazioni termiche.")
        btn_add_term = QPushButton("+ Variazione Termica")
        btn_add_term.setObjectName("btnAggiungi")
        btn_add_term.clicked.connect(self.add_termica)
        layout_term.addWidget(QLabel("Asta:"), 0, 0); layout_term.addWidget(self.combo_asta_term, 0, 1)
        layout_term.addWidget(QLabel("Modo:"), 1, 0); layout_term.addWidget(self.combo_term_modo, 1, 1, 1, 3)
        self.lbl_T_pos, self.lbl_T_neg = QLabel("ΔT fibra sinistra (°C):"), QLabel("ΔT fibra destra (°C):")
        layout_term.addWidget(self.lbl_T_pos, 2, 0); layout_term.addWidget(self.spin_T_pos, 2, 1)
        layout_term.addWidget(self.lbl_T_neg, 2, 2); layout_term.addWidget(self.spin_T_neg, 2, 3)
        layout_term.addWidget(QLabel("α (1/°C):"), 3, 0); layout_term.addWidget(self.spin_alpha, 3, 1)
        layout_term.addWidget(btn_add_term, 3, 3)
        self._aggiorna_modo_termico()
        grp_term.setLayout(layout_term)
        scroll_layout.addWidget(grp_term)

        # === PROPRIETÀ MATERIALE E SEZIONE ===
        group_props = QGroupBox("Proprietà Materiale e Sezione (Rettangolare)")
        grid_props = QGridLayout()
        
        # Modulo di Young E
        grid_props.addWidget(QLabel("E (Pa):"), 0, 0)
        self.spin_E = QDoubleSpinBox()
        self.spin_E.setRange(1e9, 1e12)
        self.spin_E.setValue(2.1e11)  # Default Acciaio
        self.spin_E.setDecimals(1)
        grid_props.addWidget(self.spin_E, 0, 1)
        
        # Tensione ammissibile sigma_amm
        grid_props.addWidget(QLabel("Sigma amm (Pa):"), 1, 0)
        self.spin_sigma = QDoubleSpinBox()
        self.spin_sigma.setRange(1e6, 1e9)
        self.spin_sigma.setValue(235e6)  # Default S235
        grid_props.addWidget(self.spin_sigma, 1, 1)
        
        # Base b
        grid_props.addWidget(QLabel("Base b (m):"), 2, 0)
        self.spin_b = QDoubleSpinBox()
        self.spin_b.setRange(0.001, 10.0)
        self.spin_b.setDecimals(4)
        self.spin_b.setValue(0.2)
        self.spin_b.setSingleStep(0.05)
        grid_props.addWidget(self.spin_b, 2, 1)
        
        # Altezza h
        grid_props.addWidget(QLabel("Altezza h (m):"), 3, 0)
        self.spin_h = QDoubleSpinBox()
        self.spin_h.setRange(0.001, 10.0)
        self.spin_h.setDecimals(4)
        self.spin_h.setValue(0.4)
        self.spin_h.setSingleStep(0.05)
        grid_props.addWidget(self.spin_h, 3, 1)

        # Margine di Sicurezza Desiderato con Operatore Relazionale
        grid_props.addWidget(QLabel("MS verifica (punto 3):"), 4, 0)
        
        layout_ms = QHBoxLayout()
        self.combo_ms_op = QComboBox()
        self.combo_ms_op.addItems(["≥ (Maggiore o uguale)", "> (Maggiore di)", "= (Pari a)",
                                   "≤ (Minore o uguale)", "< (Minore di)"])

        self.spin_MS = QDoubleSpinBox()
        self.spin_MS.setRange(-1.0, 10.0)   # MS < 0 e' una richiesta presente nelle prove
        self.spin_MS.setValue(0.0)
        self.spin_MS.setSingleStep(0.5)
        
        layout_ms.addWidget(self.combo_ms_op)
        layout_ms.addWidget(self.spin_MS)
        grid_props.addLayout(layout_ms, 4, 1)

        # --- NUOVO: Criterio di Verifica ---
        grid_props.addWidget(QLabel("Criterio:"), 5, 0)
        self.combo_ver_type = QComboBox()
        self.combo_ver_type.addItems(["Presso/Tensoflessione (N/A + M/W)", "Solo Flessione (M/W)"])
        grid_props.addWidget(self.combo_ver_type, 5, 1)
        # --- EI assegnato e dimensionamento ---
        grid_props.addWidget(QLabel("EI assegnato (N·m²):"), 6, 0)
        self.spin_EI = QDoubleSpinBox()
        self.spin_EI.setRange(0, 1e15)
        self.spin_EI.setDecimals(1)
        self.spin_EI.setSpecialValueText("— calcola da E·b·h³/12 —")
        self.spin_EI.setValue(0)
        self.spin_EI.setToolTip("Se la prova fornisce EI (o una rigidezza assegnata) inseriscilo qui: ha la precedenza su E, b, h\nnel calcolo di X1 con cedimenti/distorsioni.")
        grid_props.addWidget(self.spin_EI, 6, 1)

        grid_props.addWidget(QLabel("Dimensionamento:"), 7, 0)
        self.combo_design = QComboBox()
        self.combo_design.addItems(["Nessuno (solo verifica)", "Incognita h (b assegnata)", "Incognita b (h assegnata)", "Rapporto h/b assegnato (incognita b)"])
        grid_props.addWidget(self.combo_design, 7, 1)

        self.lbl_ratio = QLabel("Rapporto h/b:")
        grid_props.addWidget(self.lbl_ratio, 8, 0)
        self.spin_ratio = QDoubleSpinBox()
        self.spin_ratio.setRange(0.1, 20.0); self.spin_ratio.setDecimals(3); self.spin_ratio.setValue(2.0)
        grid_props.addWidget(self.spin_ratio, 8, 1)
        # MS della richiesta di progetto (punto 4 / pendolo): separato dall'MS di verifica del punto 3
        self.lbl_ms_design = QLabel("MS richiesto (punto 4):")
        grid_props.addWidget(self.lbl_ms_design, 9, 0)
        self.w_ms_design = QWidget()
        layout_ms_d = QHBoxLayout(self.w_ms_design)
        layout_ms_d.setContentsMargins(0, 0, 0, 0)
        self.combo_ms_op_design = QComboBox()
        self.combo_ms_op_design.addItems(["≥ (Maggiore o uguale)", "> (Maggiore di)", "= (Pari a)",
                                          "≤ (Minore o uguale)", "< (Minore di)"])
        self.spin_MS_design = QDoubleSpinBox()
        self.spin_MS_design.setRange(-1.0, 10.0)
        self.spin_MS_design.setValue(0.0)
        self.spin_MS_design.setSingleStep(0.5)
        self.w_ms_design.setToolTip("MS richiesto per la relazione fra i parametri di sezione (dimensionamento) e per il pendolo.\n"
                                    "Indipendente dall'MS del punto 3 (verifica della sezione maggiormente sollecitata).")
        layout_ms_d.addWidget(self.combo_ms_op_design)
        layout_ms_d.addWidget(self.spin_MS_design)
        grid_props.addWidget(self.w_ms_design, 9, 1)
        grid_props.addWidget(QLabel("Le modifiche si applicano subito al JSON: non serve salvare."), 10, 0, 1, 2)
        grid_props.setColumnStretch(1, 1)

        group_props.setLayout(grid_props)
        props_layout.addWidget(group_props)
        # === SEZIONE DEL PENDOLO (BIELLA): usata solo se esiste almeno un pendolo / doppio pendolo ===
        grp_pend = QGroupBox("Sezione del Pendolo (biella) - tirante o puntone")
        grid_pend = QGridLayout()
        self.combo_pend_mode = QComboBox()
        self.combo_pend_mode.addItems(["Dimensiona (usa MS richiesto)", "Verifica (sezione assegnata)"])
        self.combo_pend_shape = QComboBox()
        self.combo_pend_shape.addItems(["Solo area A", "Rettangolare b x h", "Circolare piena (d)"])
        self.spin_pend_A = QDoubleSpinBox(); self.spin_pend_A.setRange(0, 100); self.spin_pend_A.setDecimals(8)
        self.spin_pend_b = QDoubleSpinBox(); self.spin_pend_b.setRange(0, 10); self.spin_pend_b.setDecimals(5)
        self.spin_pend_h = QDoubleSpinBox(); self.spin_pend_h.setRange(0, 10); self.spin_pend_h.setDecimals(5)
        self.spin_pend_d = QDoubleSpinBox(); self.spin_pend_d.setRange(0, 10); self.spin_pend_d.setDecimals(5)
        self.spin_pend_ratio = QDoubleSpinBox(); self.spin_pend_ratio.setRange(0.1, 20); self.spin_pend_ratio.setDecimals(3)
        self.spin_pend_ratio.setValue(1.0)
        self.spin_pend_ratio.setToolTip("Solo dimensionamento rettangolare: h = rapporto * b")
        grid_pend.addWidget(QLabel("Modo:"), 0, 0); grid_pend.addWidget(self.combo_pend_mode, 0, 1)
        grid_pend.addWidget(QLabel("Forma:"), 0, 2); grid_pend.addWidget(self.combo_pend_shape, 0, 3)
        grid_pend.addWidget(QLabel("A (m²):"), 1, 0); grid_pend.addWidget(self.spin_pend_A, 1, 1)
        grid_pend.addWidget(QLabel("d (m):"), 1, 2); grid_pend.addWidget(self.spin_pend_d, 1, 3)
        grid_pend.addWidget(QLabel("b (m):"), 2, 0); grid_pend.addWidget(self.spin_pend_b, 2, 1)
        grid_pend.addWidget(QLabel("h (m):"), 2, 2); grid_pend.addWidget(self.spin_pend_h, 2, 3)
        grid_pend.addWidget(QLabel("h/b:"), 3, 0); grid_pend.addWidget(self.spin_pend_ratio, 3, 1)
        grid_pend.addWidget(QLabel("MS richiesto e σ_amm: stessi valori del pannello Proprietà."), 3, 2, 1, 2)
        grp_pend.setLayout(grid_pend)
        self.grp_pend = grp_pend
        props_layout.addWidget(grp_pend)
        props_layout.addStretch()

        # === RIEPILOGO: ogni riga e' un elemento reale della struttura, selezionabile, modificabile ed eliminabile ===
        scroll_layout.addWidget(QLabel("Riepilogo (doppio clic su una riga = modifica):"))
        self.lbl_avviso = QLabel("")
        self.lbl_avviso.setWordWrap(True)
        self.lbl_avviso.setVisible(False)
        scroll_layout.addWidget(self.lbl_avviso)
        self.list_riepilogo = QListWidget()
        self.list_riepilogo.setMinimumHeight(180)
        self.list_riepilogo.itemDoubleClicked.connect(lambda _it: self.modifica_selezionato())
        self.list_riepilogo.currentItemChanged.connect(self._on_selezione_riepilogo)
        scroll_layout.addWidget(self.list_riepilogo)

        riga_btn = QHBoxLayout()
        self.btn_modifica = QPushButton("✏️ Modifica selezionato")
        self.btn_modifica.clicked.connect(self.modifica_selezionato)
        self.btn_elimina = QPushButton("🗑️ Elimina selezionato")
        self.btn_elimina.clicked.connect(self.elimina_selezionato)
        self.btn_annulla = QPushButton("↩ Annulla ultima modifica")
        self.btn_annulla.clicked.connect(self.annulla_ultima)
        self.btn_annulla.setEnabled(False)
        btn_reset = QPushButton("Svuota tutto")
        btn_reset.clicked.connect(self.reset_builder)
        self.btn_modifica.setText("✏️ Modifica")
        self.btn_elimina.setText("🗑️ Elimina")
        self.btn_annulla.setText("↩ Annulla")
        self.btn_modifica.setObjectName("btnSecondario")
        self.btn_annulla.setObjectName("btnSecondario")
        self.btn_elimina.setObjectName("btnPericolo")
        btn_reset.setObjectName("btnPericolo")
        for b_ in (self.btn_modifica, self.btn_elimina, self.btn_annulla, btn_reset):
            riga_btn.addWidget(b_)
        scroll_layout.addLayout(riga_btn)
        scroll_layout.addStretch()
        self._aggiorna_visibilita_props()

        scroll.setWidget(scroll_content)
        layout.addWidget(scroll)

    def _proprieta_da_ui(self, dati):
        """Copia nel dizionario `dati` i valori correnti di materiale, sezione, MS, criterio, EI e dimensionamento."""
        mat = dati.setdefault("material", {})
        sec = dati.setdefault("section", {})
        mat["E"] = self.spin_E.value()
        mat["sigma_amm"] = self.spin_sigma.value()
        sec["b"] = self.spin_b.value()
        sec["h"] = self.spin_h.value()
        # alpha solo se impostato: cosi' un carico termico senza alpha viene segnalato dal solver invece di usare un valore a caso
        if self.spin_alpha.value() > 0:
            mat["alpha"] = self.spin_alpha.value()
        else:
            mat.pop("alpha", None)
        if self.spin_EI.value() > 0:
            mat["EI"] = self.spin_EI.value()
        else:
            mat.pop("EI", None)

        mat["MS_target"] = self.spin_MS.value()
        op_text = self.combo_ms_op.currentText().strip()
        # Si guarda il primo carattere: "≥ (Maggiore o uguale)" ecc. non vanno confusi con "="
        mat["MS_operator"] = {"≥": ">=", ">": ">", "=": "=", "≤": "<=", "<": "<"}.get(op_text[:1], ">=")
        mat["MS_design_target"] = self.spin_MS_design.value()
        op_d = self.combo_ms_op_design.currentText().strip()
        mat["MS_design_operator"] = {"≥": ">=", ">": ">", "=": "=", "≤": "<=", "<": "<"}.get(op_d[:1], ">=")
        mat["verification_type"] = "bending_only" if "Solo Flessione" in self.combo_ver_type.currentText() else "combined"

        sec["design"] = ["none", "h", "b", "ratio"][self.combo_design.currentIndex()]
        if sec["design"] == "ratio":
            sec["ratio"] = self.spin_ratio.value()
        else:
            sec.pop("ratio", None)
        self._pendolo_da_ui(dati)
        return dati

    def _aggiorna_visibilita_props(self, _i=None):
        """Mostra solo cio' che serve: la sezione del pendolo se esiste un pendolo, h/b solo col dimensionamento a rapporto."""
        ha_pend = any(s.get("pendolo") or s.get("bipendolo") for s in self.builder_data.get("supports", []))
        self.grp_pend.setVisible(ha_pend)
        rapporto = self.combo_design.currentIndex() == 3
        self.lbl_ratio.setVisible(rapporto)
        self.spin_ratio.setVisible(rapporto)
        # l'MS di progetto serve col dimensionamento della sezione o con un pendolo da dimensionare
        serve_ms_d = self.combo_design.currentIndex() != 0 or ha_pend
        self.lbl_ms_design.setVisible(serve_ms_d)
        self.w_ms_design.setVisible(serve_ms_d)

    def _pendolo_da_ui(self, dati):
        """Scrive 'pendolo_section' solo se la struttura ha almeno un pendolo / doppio pendolo inclinato."""
        if not any(s.get("pendolo") or s.get("bipendolo") for s in dati.get("supports", [])):
            dati.pop("pendolo_section", None)
            return dati
        forma = ["area", "rettangolare", "circolare"][self.combo_pend_shape.currentIndex()]
        modo = "verifica" if self.combo_pend_mode.currentIndex() == 1 else "design"
        ps = {"shape": forma, "mode": modo}
        if modo == "verifica":
            if forma == "area": ps["A"] = self.spin_pend_A.value()
            elif forma == "rettangolare": ps.update({"b": self.spin_pend_b.value(), "h": self.spin_pend_h.value()})
            else: ps["d"] = self.spin_pend_d.value()
        elif forma == "rettangolare":
            ps["ratio"] = self.spin_pend_ratio.value()
        dati["pendolo_section"] = ps
        return dati

    def salva_proprieta_materiale(self):
        """Salva E, sigma_amm, b, h, EI, MS_target, criterio e dimensionamento dentro builder_data."""
        if not hasattr(self, 'builder_data') or self._sync_in_corso:
            return
        self._proprieta_da_ui(self.builder_data)
        self.update_json_from_builder()

    def _collega_segnali_proprieta(self):
        """Ogni modifica del pannello proprieta' aggiorna subito il JSON: il JSON e' la sola fonte di verita'
        usata dal solver, quindi pannello e JSON non possono piu' divergere."""
        for w in (self.spin_E, self.spin_sigma, self.spin_b, self.spin_h, self.spin_MS, self.spin_MS_design,
                  self.spin_EI, self.spin_ratio, self.spin_alpha,
                  self.spin_pend_A, self.spin_pend_b, self.spin_pend_h, self.spin_pend_d, self.spin_pend_ratio):
            w.valueChanged.connect(lambda _v=None: self.salva_proprieta_materiale())
        for w in (self.combo_ms_op, self.combo_ms_op_design, self.combo_ver_type, self.combo_design, self.combo_pend_mode, self.combo_pend_shape):
            w.currentIndexChanged.connect(lambda _i=None: self.salva_proprieta_materiale())
        self.combo_design.currentIndexChanged.connect(self._aggiorna_visibilita_props)

    # ------------------------------------------------------------------
    # UTILITA' COMUNI: annulla, avvisi, validazione, aggiornamento coordinato di tutte le viste
    # ------------------------------------------------------------------
    def _snapshot(self):
        """Salva lo stato corrente prima di una modifica: 'Annulla ultima modifica' lo ripristina (ultime 30)."""
        try:
            self._undo.append(json.dumps(self.builder_data))
            del self._undo[:-30]
            self.btn_annulla.setEnabled(True)
        except Exception:
            pass

    def annulla_ultima(self):
        if not self._undo:
            return
        try:
            dati = json.loads(self._undo.pop())
        except Exception:
            return
        self.btn_annulla.setEnabled(bool(self._undo))
        self._sel = None
        self.sync_builder_da_json(dati)
        self._avviso("Ultima modifica annullata.", errore=False)

    def _avviso(self, msg, errore=True):
        """Messaggio visibile sopra il riepilogo (arancione = attenzione, verde = informazione)."""
        if not msg:
            self.lbl_avviso.setVisible(False)
            return
        colore, sfondo = ("#842029", "#f8d7da") if errore else ("#0f5132", "#d1e7dd")
        self.lbl_avviso.setStyleSheet(f"color:{colore}; background-color:{sfondo}; padding:6px; border-radius:4px;")
        self.lbl_avviso.setText(msg)
        self.lbl_avviso.setVisible(True)

    def _errori_dati(self, dati):
        """Stessi controlli del solutore (id duplicati, nodi coincidenti, vincoli sovrapposti, riferimenti inesistenti...)."""
        if normalizza_dati is None or valida_dati is None:
            return []
        try:
            return [html.unescape(re.sub(r"<[^>]+>", "", m)) for m in valida_dati(normalizza_dati(copy.deepcopy(dati)))]
        except Exception as e:
            return [f"Dati incompleti o non numerici: {e}"]

    def _mostra_errori(self):
        if not self.builder_data.get("nodes"):
            self._avviso("")
            return
        errs = self._errori_dati(self.builder_data)
        self._avviso("⚠ " + "  |  ".join(errs) if errs else "")

    def _dopo_modifica(self):
        """Allinea in un colpo solo tendine, riepilogo, JSON, controlli e anteprima."""
        self._pendolo_da_ui(self.builder_data)      # la sezione del pendolo compare/scompare con i pendoli
        self._aggiorna_visibilita_props()
        self.update_dropdowns()
        self.refresh_riepilogo()
        self.update_json_from_builder()
        self._mostra_errori()
        self.update_live_preview()

    # ------------------------------------------------------------------
    # RIEPILOGO: costruito dai dati, ogni riga e' un elemento reale (categoria, indice)
    # ------------------------------------------------------------------
    @staticmethod
    def _g(v):
        try:
            return f"{float(v):g}"
        except Exception:
            return str(v)

    def _desc_vincolo(self, s):
        idx = tipo_da_vincolo(s)
        nome = TIPI_VINCOLO[idx] if idx is not None else "Vincolo (da JSON)"
        extra = ""
        if idx == 4: extra = f" [angolo {self._g(s.get('angle', 0))}°]"
        if idx in (IDX_PENDOLO, IDX_BIPENDOLO_INCL): extra = f" [asse {self._g(angolo_vincolo(s))}°]"
        return f"Vincolo N{s.get('node_id')}: {nome}{extra}"

    def _righe_riepilogo(self):
        d, g = self.builder_data, self._g
        righe = []
        for i, n in enumerate(d.get("nodes", [])):
            righe.append((("nodo", i), f"Nodo {n['id']}: X={g(n['x'])}, Y={g(n['y'])}"))
        for i, e in enumerate(d.get("elements", [])):
            righe.append((("asta", i), f"Asta {e['id']}: Nodo {e['n1']} → Nodo {e['n2']}"))
        for i, s in enumerate(d.get("supports", [])):
            righe.append((("vincolo", i), self._desc_vincolo(s)))
        for i, s in enumerate(d.get("supports", [])):
            parti = []
            if s.get("delta_u"): parti.append(f"Δx={g(float(s['delta_u']) * 1000)} mm")
            if s.get("delta_v"): parti.append(f"Δy={g(float(s['delta_v']) * 1000)} mm")
            if s.get("delta_n"): parti.append(f"Δn={g(float(s['delta_n']) * 1000)} mm")
            if s.get("delta_w"): parti.append(f"Δφ={g(s['delta_w'])} rad")
            if parti:
                righe.append((("cedimento", i), f"Cedimento N{s.get('node_id')}: " + ", ".join(parti)))
        for i, p in enumerate(d.get("point_loads", [])):
            righe.append((("carico_nodo", i), f"Carico su N{p['node_id']}: Fx={g(p.get('Fx', 0))}, Fy={g(p.get('Fy', 0))}, M={g(p.get('M', 0))}"))
        for i, q in enumerate(d.get("element_loads", [])):
            parti_q = []
            if q.get("q_trans", 0) != 0: parti_q.append(f"q⊥ = {g(q['q_trans'])}")
            if q.get("q_x", 0) != 0: parti_q.append(f"qx = {g(q['q_x'])}")
            if q.get("q_y", 0) != 0: parti_q.append(f"qy = {g(q['q_y'])}")
            righe.append((("carico_distr", i), f"Carico distribuito su Asta {q['elem_id']}: {', '.join(parti_q) or 'q = 0'} N/m "
                                              f"({q.get('tipo', 'costante')}, {q.get('riferimento', 'lunghezza')})"))
        for i, t in enumerate(d.get("thermal_loads", [])):
            righe.append((("termica", i), f"ΔT su Asta {t['elem_id']}: sinistra {g(t.get('T_pos', 0))} °C / destra {g(t.get('T_neg', 0))} °C"))
        for i, c in enumerate(d.get("internal_hinges", [])):
            righe.append((("cerniera", i), f"Cerniera Interna al Nodo N{c['node_id']}" + (f" (aste incernierate: {c['elems']})" if c.get("elems") else "")
                          + (f" [carichi del nodo sull'asta {c['carichi_su']}]" if c.get("carichi_su") not in (None, "") else "")))
        return righe

    def refresh_riepilogo(self):
        """Ricostruisce l'elenco dai dati e ripristina la selezione. Non genera segnali di selezione."""
        self.list_riepilogo.blockSignals(True)
        try:
            self.list_riepilogo.clear()
            riga_sel = -1
            for k, (chiave, testo) in enumerate(self._righe_riepilogo()):
                it = QListWidgetItem(testo)
                it.setData(Qt.ItemDataRole.UserRole, chiave)
                self.list_riepilogo.addItem(it)
                if chiave == self._sel:
                    riga_sel = k
            if riga_sel >= 0:
                self.list_riepilogo.setCurrentRow(riga_sel)
            else:
                self._sel = None
        finally:
            self.list_riepilogo.blockSignals(False)

    def _on_selezione_riepilogo(self, corrente, _precedente):
        self._sel = corrente.data(Qt.ItemDataRole.UserRole) if corrente is not None else None
        self.update_live_preview()      # evidenzia l'elemento selezionato

    def _aggiorna_etichette_vincolo(self, _i=None):
        idx = self.combo_tipo_vinc.currentIndex()
        pend = idx in (IDX_PENDOLO, IDX_BIPENDOLO_INCL)
        self.lbl_ang_vinc.setText("Asse biella (°):" if pend else "Angolo(°):")
        self.spin_ang_vinc.setToolTip(
            "Direzione dell'ASSE della biella, da +x in senso antiorario, orientata dal NODO verso il vincolo a terra.\n"
            "Es. biella che va dal nodo verso il basso-destra a 45° => -45°." if pend
            else "Gradi rispetto all'orizzontale (solo per il carrello inclinato)")
        self.lbl_lung_vinc.setVisible(pend)
        self.spin_lung_vinc.setVisible(pend)

    def update_dropdowns(self):
        """Aggiorna i menu a tendina con i nodi e le aste inserite."""
        combos_nodo = (self.combo_rif_nodo, self.combo_n1, self.combo_n2, self.combo_nodo_vinc, self.combo_nodo_carico,
                       self.combo_nodo_cerniera_int, self.combo_nodo_ced)
        combos_asta = (self.combo_asta_carico, self.combo_asta_term)
        memo = {id(c): c.currentData() for c in combos_nodo + combos_asta}      # non perdere la scelta a ogni modifica
        for c in combos_nodo + combos_asta:
            c.blockSignals(True); c.clear()
        for node in self.builder_data["nodes"]:
            label = f"N{node['id']} ({self._g(node['x'])}, {self._g(node['y'])})"
            for c in combos_nodo:
                c.addItem(label, node['id'])
        for el in self.builder_data["elements"]:
            label_asta = f"Asta {el['id']} (Da N{el['n1']} a N{el['n2']})"
            for c in combos_asta:
                c.addItem(label_asta, el['id'])
        for c in combos_nodo + combos_asta:
            k = c.findData(memo.get(id(c)))
            if k >= 0:
                c.setCurrentIndex(k)
            c.blockSignals(False)

    def update_json_from_builder(self):
        """Traduce il dizionario python in testo JSON nella tab apposita."""
        json_text = json.dumps(self.builder_data, indent=2)
        prev = self._sync_in_corso
        self._sync_in_corso = True      # scrittura programmatica: non deve rientrare come "modifica manuale"
        try:
            self.text_json.setPlainText(json_text)
        finally:
            self._sync_in_corso = prev

    # ------------------------------------------------------------------
    # AGGIUNTA ELEMENTI
    # ------------------------------------------------------------------
    def add_nodo_cartesiano(self):
        self._snapshot()
        self.builder_data["nodes"].append({"id": self.node_id_counter, "x": self.spin_x.value(), "y": self.spin_y.value()})
        self.node_id_counter += 1
        self._dopo_modifica()

    def add_nodo_polare(self):
        ref_id = self.combo_rif_nodo.currentData()
        if not ref_id: return
        ref_node = next((n for n in self.builder_data["nodes"] if n["id"] == ref_id), None)
        if not ref_node: return
        L = self.spin_L.value()
        ang_rad = math.radians(self.spin_ang.value())
        x = ref_node["x"] + L * math.cos(ang_rad)
        y = ref_node["y"] + L * math.sin(ang_rad)
        # NIENTE round(..., 4): si elimina solo il rumore di macchina (~1e-16) mantenendo
        # 12 cifre significative; il solver recupera i valori esatti (es. 2*sqrt(2)).
        x = 0.0 if abs(x) < 1e-12 else float(f"{x:.12g}")
        y = 0.0 if abs(y) < 1e-12 else float(f"{y:.12g}")
        self._snapshot()
        self.builder_data["nodes"].append({"id": self.node_id_counter, "x": x, "y": y})
        self.node_id_counter += 1
        self._dopo_modifica()

    def add_asta(self):
        n1 = self.combo_n1.currentData()
        n2 = self.combo_n2.currentData()
        if not n1 or not n2 or n1 == n2:
            self._avviso("Scegli due nodi diversi per l'asta.")
            return
        self._snapshot()
        self.builder_data["elements"].append({"id": self.elem_id_counter, "n1": n1, "n2": n2})
        self.elem_id_counter += 1
        self._dopo_modifica()

    def add_vincolo(self):
        n_id = self.combo_nodo_vinc.currentData()
        if not n_id: return
        support = costruisci_vincolo(self.combo_tipo_vinc.currentIndex(), n_id, self.spin_ang_vinc.value(), self.spin_lung_vinc.value())
        # Si rifiuta solo un vincolo che crea NUOVI errori (es. due vincoli sulla stessa traslazione o due bielle parallele).
        # Piu' vincoli su uno stesso nodo sono ammessi: due pendoli non paralleli formano un doppio pendolo.
        prova = copy.deepcopy(self.builder_data)
        prova["supports"].append(support)
        nuovi = [e for e in self._errori_dati(prova) if e not in self._errori_dati(self.builder_data)]
        if nuovi:
            self._avviso("Vincolo NON aggiunto: " + "  |  ".join(nuovi))
            return
        self._snapshot()
        self.builder_data["supports"].append(support)
        self._dopo_modifica()

    def add_carico(self):
        """Salva un carico concentrato su un nodo"""
        n_id = self.combo_nodo_carico.currentData()
        if not n_id: return
        fx, fy, m = self.spin_Fx.value(), self.spin_Fy.value(), self.spin_M.value()
        f_inc = self.spin_Finc.value()
        if f_inc != 0:
            ang_ = math.radians(self.spin_Finc_ang.value())
            fx += f_inc * math.cos(ang_)
            fy += f_inc * math.sin(ang_)
            # ripulisce il rumore numerico delle direzioni ortogonali (cos 90 = 6e-17 ...)
            fx, fy = round(fx, 4) + 0.0, round(fy, 4) + 0.0
        if fx == 0 and fy == 0 and m == 0: return
        self._snapshot()
        self.builder_data["point_loads"].append({"node_id": n_id, "Fx": fx, "Fy": fy, "M": m})
        self.spin_Finc.setValue(0.0)
        self._dopo_modifica()

    @staticmethod
    def _applica_cedimenti(supports_nodo, dx, dy, dn, dw):
        """Assegna i cedimenti ai vincoli del nodo, solo sulle componenti realmente vincolate. Ritorna le componenti applicate."""
        fatti = set()
        for sup in supports_nodo:
            if sup.get("rotated_roller"):
                if dn != 0 and "n" not in fatti: sup["delta_n"] = dn; fatti.add("n")
            else:
                if dx != 0 and sup.get("fix_u") and "x" not in fatti: sup["delta_u"] = dx; fatti.add("x")
                if dy != 0 and sup.get("fix_v") and "y" not in fatti: sup["delta_v"] = dy; fatti.add("y")
            if dw != 0 and sup.get("fix_w") and "w" not in fatti: sup["delta_w"] = dw; fatti.add("w")
        return fatti

    def add_cedimento(self):
        """Assegna cedimenti (mm) / rotazione imposta (rad) ai vincoli gia' esistenti del nodo"""
        n_id = self.combo_nodo_ced.currentData()
        if not n_id: return
        sups = [v for v in self.builder_data["supports"] if v["node_id"] == n_id]
        if not sups:
            self._avviso(f"Nessun vincolo sul Nodo N{n_id}: inserisci prima il vincolo (sez. 3).")
            return
        dx, dy, dn, dw = (self.spin_ced_x.value() / 1000.0, self.spin_ced_y.value() / 1000.0,
                          self.spin_ced_n.value() / 1000.0, self.spin_ced_w.value())
        if dx == 0 and dy == 0 and dn == 0 and dw == 0: return
        self._snapshot()
        fatti = self._applica_cedimenti(sups, dx, dy, dn, dw)
        richiesti = {k for k, v in (("x", dx), ("y", dy), ("n", dn), ("w", dw)) if v != 0}
        if not fatti:
            self._undo.pop()
            self.btn_annulla.setEnabled(bool(self._undo))
            self._avviso(f"Nessun cedimento assegnato: il vincolo del Nodo N{n_id} non blocca le componenti richieste "
                         "(un cedimento ha senso solo nelle direzioni vincolate).")
            return
        self._dopo_modifica()
        if richiesti - fatti:
            self._avviso("Alcune componenti sono state ignorate perche' non vincolate nel nodo: "
                         + ", ".join(sorted(richiesti - fatti)))

    def _aggiorna_modo_termico(self):
        """Nei modi 'farfalla' e 'differenza totale' il dato e' uno solo (campo sinistra); il campo destra e' dedotto."""
        modo = self.combo_term_modo.currentIndex()
        self.spin_T_neg.setEnabled(modo == 0)
        self.lbl_T_pos.setText("ΔT fibra sinistra (°C):" if modo == 0 else "ΔT (°C):")
        self.lbl_T_neg.setText("ΔT fibra destra (°C):" if modo == 0 else "(dedotta)")

    def add_termica(self):
        """Salva una variazione termica a gradiente su un'asta"""
        el_id = self.combo_asta_term.currentData()
        if not el_id: return
        modo = self.combo_term_modo.currentIndex()
        v_ = self.spin_T_pos.value()
        if modo == 1:   tp, tn = v_, -v_              # farfalla +DT / -DT
        elif modo == 2: tp, tn = v_ / 2.0, -v_ / 2.0  # DT = differenza fra le due fibre
        else:           tp, tn = v_, self.spin_T_neg.value()
        if tp == 0 and tn == 0: return
        if self.spin_alpha.value() <= 0:
            QMessageBox.warning(self, "Coefficiente α mancante",
                                "Inserisci il coefficiente di dilatazione termica α (campo «α (1/°C)») prima di aggiungere la variazione termica.")
            return
        self._snapshot()
        self.builder_data.setdefault("thermal_loads", []).append({"elem_id": el_id, "T_pos": tp, "T_neg": tn})
        self.builder_data.setdefault("material", {})["alpha"] = self.spin_alpha.value()
        self._dopo_modifica()

    @staticmethod
    def _testo_a_elenco_aste(testo):
        return [int(t) for t in testo.replace(";", ",").replace(" ", ",").split(",") if t.strip()]

    @staticmethod
    def _valida_cerniera(dati, n_id, elems, carichi_su=None):
        """Messaggio d'errore oppure None. Stesse regole per inserimento e modifica."""
        incidenti = [int(e["id"]) for e in dati["elements"] if n_id in (e["n1"], e["n2"])]
        if carichi_su is not None:
            incernierate = list(elems) if elems else incidenti[:1]
            if int(carichi_su) not in incernierate:
                return (f"'Carichi su asta' = {carichi_su} deve essere una delle aste incernierate nel Nodo N{n_id} "
                        f"({incernierate}): è l'asta il cui corpo porta coppie, vincoli e forze del nodo")
        if any(i not in incidenti for i in elems):
            return f"Le aste incernierate devono convergere nel Nodo N{n_id} (aste del nodo: {incidenti})"
        if len(incidenti) < 2:
            return f"Nel Nodo N{n_id} convergono meno di 2 aste: una cerniera interna richiede almeno 2 aste"
        if len(incidenti) > 2 and not elems:
            return f"Nel Nodo N{n_id} convergono {len(incidenti)} aste {incidenti}: indica quali sono incernierate"
        return None

    def add_cerniera_interna(self):
        """Salva una cerniera interna su un nodo"""
        n_id = self.combo_nodo_cerniera_int.currentData()
        if not n_id: return
        try:
            elems = self._testo_a_elenco_aste(self.edit_elems_cerniera.text())
        except ValueError:
            self._avviso("Aste incernierate: usa solo numeri separati da virgola")
            return
        try:
            su_ = self._testo_a_elenco_aste(self.edit_carichi_su.text())
        except ValueError:
            self._avviso("Carichi su asta: usa un solo numero (id dell'asta)")
            return
        if len(su_) > 1:
            self._avviso("Carichi su asta: indica una sola asta")
            return
        carichi_su = su_[0] if su_ else None
        errore = self._valida_cerniera(self.builder_data, n_id, elems, carichi_su)
        if errore:
            self._avviso(errore)
            return
        cerniera = {"node_id": n_id}
        if elems:
            cerniera["elems"] = elems
        if carichi_su is not None:
            cerniera["carichi_su"] = carichi_su
        self._snapshot()
        self.builder_data.setdefault("internal_hinges", []).append(cerniera)
        self.edit_elems_cerniera.clear()
        self.edit_carichi_su.clear()
        self._dopo_modifica()

    def reset_builder(self):
        if self.builder_data.get("nodes") and QMessageBox.question(
                self, "Svuota tutto", "Cancellare l'intera struttura? (puoi annullare con '↩ Annulla ultima modifica')",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No) != QMessageBox.StandardButton.Yes:
            return
        self._snapshot()
        self.builder_data = {
            "nodes": [], "elements": [], "supports": [],
            "point_loads": [], "element_loads": [], "thermal_loads": [],
            "material": {"E": 2.1e11, "sigma_amm": 235e6}, "section": {"b": 0.2, "h": 0.4}
        }
        self.node_id_counter = 1
        self.elem_id_counter = 1
        self._sel = None
        self._dopo_modifica()

    # ------------------------------------------------------------------
    # MODIFICA ED ELIMINAZIONE DI SINGOLI ELEMENTI (senza rifare la struttura)
    # ------------------------------------------------------------------
    @staticmethod
    def _spin(valore, lo=-9999.0, hi=9999.0, dec=3):
        w = QDoubleSpinBox(); w.setRange(lo, hi); w.setDecimals(dec); w.setValue(float(valore))
        return w

    def _combo_nodi(self, dati, nodo_id):
        c = QComboBox()
        for n in dati["nodes"]:
            c.addItem(f"N{n['id']} ({self._g(n['x'])}, {self._g(n['y'])})", n["id"])
        c.setCurrentIndex(max(0, c.findData(nodo_id)))
        return c

    def _combo_aste(self, dati, asta_id):
        c = QComboBox()
        for e in dati["elements"]:
            c.addItem(f"Asta {e['id']} (N{e['n1']} → N{e['n2']})", e["id"])
        c.setCurrentIndex(max(0, c.findData(asta_id)))
        return c

    def _esegui_dialogo(self, titolo, righe):
        dlg = QDialog(self)
        dlg.setWindowTitle(titolo)
        form = QFormLayout(dlg)
        for etichetta, widget in righe:
            form.addRow(etichetta, widget)
        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        form.addRow(bb)
        # IMPORTANTE: il chiamante legge i valori dai widget DOPO la chiusura del dialogo. Se il dialogo venisse
        # distrutto all'uscita da questa funzione, i widget sarebbero eliminati (RuntimeError -> crash del programma).
        # Lo si tiene vivo in self finche' non ne viene aperto un altro.
        vecchio = getattr(self, "_dlg_vivo", None)
        self._dlg_vivo = dlg
        if vecchio is not None:
            vecchio.deleteLater()
        return dlg.exec() == QDialog.DialogCode.Accepted

    def modifica_selezionato(self):
        it = self.list_riepilogo.currentItem()
        if it is None:
            self._avviso("Seleziona prima una riga del riepilogo da modificare.")
            return
        cat, i = it.data(Qt.ItemDataRole.UserRole)
        nuovo = copy.deepcopy(self.builder_data)
        try:
            ok = getattr(self, f"_edit_{cat}")(nuovo, i)
        except Exception:
            # un errore in un dialogo non deve chiudere il programma
            self._avviso("Modifica non riuscita: " + traceback.format_exc().strip().splitlines()[-1])
            return
        if ok:
            self._snapshot()
            self.builder_data = nuovo
            self._dopo_modifica()

    def _edit_nodo(self, d, i):
        n = d["nodes"][i]
        sx, sy = self._spin(n["x"]), self._spin(n["y"])
        if not self._esegui_dialogo(f"Modifica Nodo {n['id']}", [("X (m):", sx), ("Y (m):", sy)]):
            return False
        n["x"], n["y"] = sx.value(), sy.value()
        return True

    def _edit_asta(self, d, i):
        e = d["elements"][i]
        c1, c2 = self._combo_nodi(d, e["n1"]), self._combo_nodi(d, e["n2"])
        if not self._esegui_dialogo(f"Modifica Asta {e['id']}", [("Nodo inizio:", c1), ("Nodo fine:", c2)]):
            return False
        if c1.currentData() == c2.currentData():
            self._avviso("Un'asta deve collegare due nodi diversi: modifica annullata.")
            return False
        e["n1"], e["n2"] = c1.currentData(), c2.currentData()
        return True

    def _edit_vincolo(self, d, i):
        s = d["supports"][i]
        idx0 = tipo_da_vincolo(s)
        if idx0 is None:
            self._avviso("Questo vincolo non e' riconducibile a un tipo del form: modificalo nel tab JSON.")
            return False
        lung0 = math.hypot(float(s["pend_dx"]), float(s["pend_dy"])) if s.get("pend_dx") is not None and s.get("pend_dy") is not None else 0.0
        cn = self._combo_nodi(d, s["node_id"])
        ct = QComboBox(); ct.addItems(TIPI_VINCOLO); ct.setCurrentIndex(idx0)
        sa = self._spin(angolo_vincolo(s), -360, 360, 3)
        sl = self._spin(lung0, 0, 9999, 3)
        sdx, sdy = self._spin(float(s.get("delta_u", 0)) * 1000), self._spin(float(s.get("delta_v", 0)) * 1000)
        sdn, sdw = self._spin(float(s.get("delta_n", 0)) * 1000), self._spin(float(s.get("delta_w", 0)), -10, 10, 6)
        righe = [("Nodo:", cn), ("Tipo:", ct), ("Angolo / asse biella (°):", sa), ("Lunghezza biella (m, 0 = n.i.):", sl),
                 ("Cedimento Δx (mm):", sdx), ("Cedimento Δy (mm, ↑ +):", sdy), ("Cedimento Δn (mm):", sdn), ("Rotazione Δφ (rad):", sdw)]
        if not self._esegui_dialogo(f"Modifica vincolo (N{s['node_id']})", righe):
            return False
        nuovo = costruisci_vincolo(ct.currentIndex(), cn.currentData(), sa.value(), sl.value())
        self._applica_cedimenti([nuovo], sdx.value() / 1000.0, sdy.value() / 1000.0, sdn.value() / 1000.0, sdw.value())
        d["supports"][i] = nuovo
        return True

    def _edit_cedimento(self, d, i):
        return self._edit_vincolo(d, i)

    def _edit_carico_nodo(self, d, i):
        p = d["point_loads"][i]
        cn = self._combo_nodi(d, p["node_id"])
        sfx, sfy, sm = (self._spin(p.get("Fx", 0), -999999, 999999, 2), self._spin(p.get("Fy", 0), -999999, 999999, 2),
                        self._spin(p.get("M", 0), -999999, 999999, 2))
        if not self._esegui_dialogo("Modifica carico concentrato", [("Nodo:", cn), ("Fx (N):", sfx), ("Fy (N):", sfy), ("M (N·m):", sm)]):
            return False
        p.update({"node_id": cn.currentData(), "Fx": sfx.value(), "Fy": sfy.value(), "M": sm.value()})
        return True

    @staticmethod
    def costruisci_carico_distribuito(el_id, idx_tipo, idx_dir, idx_rif, q):
        tipo_q = ["costante", "triangolare_crescente", "triangolare_decrescente"][idx_tipo]
        carico = {"elem_id": el_id, "tipo": tipo_q}
        if idx_dir == 0: carico["q_trans"] = q
        elif idx_dir == 1: carico["q_y"] = -q          # + verso il basso
        else: carico["q_x"] = q                        # + verso destra
        if idx_dir != 0 and idx_rif == 1:
            carico["riferimento"] = "proiezione"
        return carico

    def _edit_carico_distr(self, d, i):
        q = d["element_loads"][i]
        comp = [k for k in ("q_trans", "q_x", "q_y") if q.get(k, 0) != 0]
        if len(comp) > 1:
            self._avviso("Questo carico ha piu' componenti (q⊥, qx, qy): modificalo nel tab JSON oppure eliminalo e reinseriscilo.")
            return False
        idx_dir, val = (1, -q["q_y"]) if "q_y" in comp else ((2, q["q_x"]) if "q_x" in comp else (0, q.get("q_trans", 0)))
        ca = self._combo_aste(d, q["elem_id"])
        ct = QComboBox(); ct.addItems(["Costante", "Triangolare (0 ➔ q)", "Triangolare (q ➔ 0)"])
        ct.setCurrentIndex({"triangolare_crescente": 1, "triangolare_decrescente": 2}.get(q.get("tipo", "costante"), 0))
        cd = QComboBox(); cd.addItems(["Perpendicolare all'asta", "Verticale (+ verso il basso)", "Orizzontale (+ verso destra)"])
        cd.setCurrentIndex(idx_dir)
        cr = QComboBox(); cr.addItems(["per unità di lunghezza d'asta", "per unità di proiezione"])
        cr.setCurrentIndex(1 if q.get("riferimento") == "proiezione" else 0)
        sq = self._spin(val, -999999, 999999, 2)
        if not self._esegui_dialogo("Modifica carico distribuito", [("Asta:", ca), ("Tipo:", ct), ("Direzione:", cd), ("Riferito a:", cr), ("q max (N/m):", sq)]):
            return False
        if sq.value() == 0:
            self._avviso("q = 0: per rimuovere il carico usa 'Elimina selezionato'.")
            return False
        d["element_loads"][i] = self.costruisci_carico_distribuito(ca.currentData(), ct.currentIndex(), cd.currentIndex(), cr.currentIndex(), sq.value())
        return True

    def _edit_termica(self, d, i):
        t = d["thermal_loads"][i]
        ca = self._combo_aste(d, t["elem_id"])
        sp_, sn_ = self._spin(t.get("T_pos", 0), -999, 999, 2), self._spin(t.get("T_neg", 0), -999, 999, 2)
        if not self._esegui_dialogo("Modifica variazione termica", [("Asta:", ca), ("ΔT fibra sinistra (°C):", sp_), ("ΔT fibra destra (°C):", sn_)]):
            return False
        t.update({"elem_id": ca.currentData(), "T_pos": sp_.value(), "T_neg": sn_.value()})
        return True

    def _edit_cerniera(self, d, i):
        c = d["internal_hinges"][i]
        cn = self._combo_nodi(d, c["node_id"])
        le = QLineEdit(", ".join(str(x) for x in c.get("elems", [])))
        le.setPlaceholderText("id aste incernierate (vuoto = automatico)")
        ls = QLineEdit("" if c.get("carichi_su") in (None, "") else str(c["carichi_su"]))
        ls.setPlaceholderText("id asta (vuoto = automatico)")
        if not self._esegui_dialogo("Modifica cerniera interna", [("Nodo:", cn), ("Aste incernierate:", le), ("Carichi del nodo su asta:", ls)]):
            return False
        try:
            elems = self._testo_a_elenco_aste(le.text())
            su_ = self._testo_a_elenco_aste(ls.text())
        except ValueError:
            self._avviso("Aste incernierate / carichi su asta: usa solo numeri")
            return False
        if len(su_) > 1:
            self._avviso("Carichi su asta: indica una sola asta")
            return False
        carichi_su = su_[0] if su_ else None
        errore = self._valida_cerniera(d, cn.currentData(), elems, carichi_su)
        if errore:
            self._avviso(errore)
            return False
        c.pop("elems", None)
        c.pop("carichi_su", None)
        c["node_id"] = cn.currentData()
        if elems:
            c["elems"] = elems
        if carichi_su is not None:
            c["carichi_su"] = carichi_su
        return True

    # --- eliminazione, con rimozione guidata di cio' che dipende dall'elemento ---
    _CHIAVE_LISTA = {"nodo": "nodes", "asta": "elements", "vincolo": "supports", "carico_nodo": "point_loads",
                     "carico_distr": "element_loads", "termica": "thermal_loads", "cerniera": "internal_hinges"}

    @staticmethod
    def _pulisci_riferimenti(d):
        """Rimuove cio' che non puo' esistere senza l'elemento eliminato. Ritorna l'elenco leggibile di cio' che e' stato tolto."""
        tolti = []
        ids_n = {int(n["id"]) for n in d["nodes"]}
        el_ok = [e for e in d["elements"] if int(e["n1"]) in ids_n and int(e["n2"]) in ids_n]
        if len(el_ok) != len(d["elements"]): tolti.append(f"{len(d['elements']) - len(el_ok)} asta/e collegata/e al nodo")
        d["elements"] = el_ok
        ids_e = {int(e["id"]) for e in el_ok}
        for chiave, campo, ids, nome in (("supports", "node_id", ids_n, "vincolo/i"), ("point_loads", "node_id", ids_n, "carico/chi sul nodo"),
                                         ("element_loads", "elem_id", ids_e, "carico/chi distribuito/i"),
                                         ("thermal_loads", "elem_id", ids_e, "variazione/i termica/he")):
            vecchi = d.get(chiave, [])
            nuovi = [x for x in vecchi if int(x[campo]) in ids]
            if len(nuovi) != len(vecchi): tolti.append(f"{len(vecchi) - len(nuovi)} {nome}")
            d[chiave] = nuovi
        cern_ok = []
        for c in d.get("internal_hinges", []):
            inc = [int(e["id"]) for e in el_ok if int(c["node_id"]) in (int(e["n1"]), int(e["n2"]))]
            if int(c["node_id"]) not in ids_n or len(inc) < 2:
                tolti.append(f"cerniera interna al nodo {c['node_id']}")
                continue
            if c.get("elems"):
                el_c = [x for x in c["elems"] if int(x) in inc]
                if not el_c:
                    tolti.append(f"cerniera interna al nodo {c['node_id']}")
                    continue
                c["elems"] = el_c
            if c.get("carichi_su") not in (None, "") and int(c["carichi_su"]) not in (inc if not c.get("elems") else [int(x) for x in c["elems"]]):
                c.pop("carichi_su", None)
            cern_ok.append(c)
        d["internal_hinges"] = cern_ok
        return tolti

    def elimina_selezionato(self):
        it = self.list_riepilogo.currentItem()
        if it is None:
            self._avviso("Seleziona prima una riga del riepilogo da eliminare.")
            return
        cat, i = it.data(Qt.ItemDataRole.UserRole)
        nuovo = copy.deepcopy(self.builder_data)
        if cat == "cedimento":
            for k in ("delta_u", "delta_v", "delta_w", "delta_n"):
                nuovo["supports"][i].pop(k, None)
            tolti = []
        else:
            nuovo[self._CHIAVE_LISTA[cat]].pop(i)
            tolti = self._pulisci_riferimenti(nuovo)
        testo = f"Eliminare: {it.text()} ?"
        if tolti:
            testo += "\n\nVerranno rimossi anche (dipendono dall'elemento):\n - " + "\n - ".join(tolti)
        if QMessageBox.question(self, "Elimina elemento", testo,
                                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No) != QMessageBox.StandardButton.Yes:
            return
        self._snapshot()
        self.builder_data = nuovo
        self._sel = None
        self._dopo_modifica()

    # --- LE RESTANTI FUNZIONI (load_json, aggiorna_anteprima_json, run_solver, esporta_pdf) 
    # RIMANGONO IDENTICHE A QUELLE DEL TUO CODICE ORIGINALE ---
    
    def load_json(self):
        file_name, _ = QFileDialog.getOpenFileName(self, "Apri JSON", "", "JSON Files (*.json)")
        if file_name:
            with open(file_name, 'r') as f:
                contenuto = f.read()
            try:
                dati = json.loads(contenuto)
            except Exception as e:
                QMessageBox.critical(self, "Errore JSON", f"JSON non valido:\n{e}")
                return
            self.sync_builder_da_json(dati)
            self.tabs.setCurrentWidget(self.tab_disegna)

    def sync_builder_da_json(self, dati, riscrivi_json=True):
        """Porta un JSON (da file o scritto a mano) dentro builder_data, pannello, contatori, riepilogo e anteprima,
        cosi' aggiungere elementi dopo non sovrascrive i dati caricati. riscrivi_json=False: il testo del tab JSON
        e' quello che l'utente sta digitando e non va rigenerato."""
        prev = self._sync_in_corso
        self._sync_in_corso = True      # i setValue sul pannello non devono riscrivere il JSON
        try:
            self._sync_builder_core(dati, riscrivi_json)
        finally:
            self._sync_in_corso = prev

    def _sync_builder_core(self, dati, riscrivi_json):
        for chiave in ("nodes", "elements", "supports", "point_loads", "element_loads", "thermal_loads"):
            dati.setdefault(chiave, [])
        dati.setdefault("material", {"E": 2.1e11, "sigma_amm": 235e6})
        dati.setdefault("section", {"b": 0.2, "h": 0.4})
        self.builder_data = dati
        self.node_id_counter = max([int(n["id"]) for n in dati["nodes"]], default=0) + 1
        self.elem_id_counter = max([int(e["id"]) for e in dati["elements"]], default=0) + 1

        self._sel = None
        mat, sec = dati["material"], dati["section"]
        if "E" in mat: self.spin_E.setValue(float(mat["E"]))
        if "sigma_amm" in mat: self.spin_sigma.setValue(float(mat["sigma_amm"]))
        if "b" in sec: self.spin_b.setValue(float(sec["b"]))
        if "h" in sec: self.spin_h.setValue(float(sec["h"]))
        self.spin_alpha.setValue(float(mat["alpha"]) if "alpha" in mat else 0)
        if "MS_target" in mat: self.spin_MS.setValue(float(mat["MS_target"]))
        _ops = {">=": 0, ">": 1, "=": 2, "<=": 3, "<": 4}
        self.combo_ms_op.setCurrentIndex(_ops.get(mat.get("MS_operator", ">="), 0))
        # vecchi JSON senza campi di progetto: si parte dai valori di verifica (comportamento precedente)
        self.spin_MS_design.setValue(float(mat.get("MS_design_target", mat.get("MS_target", 0.0))))
        self.combo_ms_op_design.setCurrentIndex(_ops.get(mat.get("MS_design_operator", mat.get("MS_operator", ">=")), 0))
        self.combo_ver_type.setCurrentIndex(1 if mat.get("verification_type") == "bending_only" else 0)
        self.spin_EI.setValue(float(mat["EI"]) if "EI" in mat else 0)
        self.combo_design.setCurrentIndex({"none": 0, "h": 1, "b": 2, "ratio": 3}.get(sec.get("design", "none"), 0))
        if "ratio" in sec: self.spin_ratio.setValue(float(sec["ratio"]))
        ps = dati.get("pendolo_section") or {}
        self.combo_pend_mode.setCurrentIndex(1 if ps.get("mode") == "verifica" else 0)
        self.combo_pend_shape.setCurrentIndex({"rettangolare": 1, "circolare": 2}.get(ps.get("shape"), 0))
        for chiave, w in (("A", self.spin_pend_A), ("b", self.spin_pend_b), ("h", self.spin_pend_h),
                          ("d", self.spin_pend_d), ("ratio", self.spin_pend_ratio)):
            if chiave in ps:
                w.setValue(float(ps[chiave]))

        self._aggiorna_visibilita_props()
        self.update_dropdowns()
        self.refresh_riepilogo()
        if riscrivi_json:
            self.update_json_from_builder()
        self._mostra_errori()
        self.update_live_preview()

    def add_carico_distribuito(self):
        """Salva un carico distribuito su un'asta (perpendicolare, verticale o orizzontale)"""
        el_id = self.combo_asta_carico.currentData()
        if not el_id: return
        q_val = self.spin_q.value()
        if q_val == 0: return
        idx_tratto = self.combo_tratto_q.currentIndex()
        t0, t1 = {0: (0.0, 1.0), 1: (0.0, 0.5), 2: (0.5, 1.0)}.get(idx_tratto, (self.spin_t0.value(), self.spin_t1.value()))
        if not (0.0 <= t0 < t1 <= 1.0):
            self._avviso("Tratto non valido: servono 0 ≤ da < a ≤ 1 (frazioni della lunghezza dell'asta).")
            return
        nuovo = copy.deepcopy(self.builder_data)
        try:
            id_tratto = self._ritaglia_tratto(nuovo, el_id, t0, t1)
        except ValueError as e_:
            self._avviso(str(e_))
            return
        carico = self.costruisci_carico_distribuito(id_tratto, self.combo_tipo_q.currentIndex(), self.combo_dir_q.currentIndex(),
                                                    self.combo_rif_q.currentIndex(), q_val)
        nuovo["element_loads"].append(carico)
        self._snapshot()
        self.builder_data = nuovo
        self._dopo_modifica()

    def _aggiorna_tratto_q(self, _i=None):
        pers = self.combo_tratto_q.currentIndex() == 3
        for w_ in (self.lbl_t0, self.spin_t0, self.lbl_t1, self.spin_t1):
            w_.setVisible(pers)

    def _spezza_asta(self, d, el_id, t):
        """Divide l'asta el_id al punto a frazione t (0<t<1) con un nodo intermedio. La prima parte mantiene l'id originale
        (Nodo Inizio -> nuovo nodo), la seconda ha un id nuovo (nuovo nodo -> Nodo Fine). Carichi costanti e variazioni
        termiche vengono replicati sulla seconda parte; le cerniere nel Nodo Fine passano alla seconda parte.
        Ritorna l'id della seconda parte."""
        el = next(e for e in d["elements"] if e["id"] == el_id)
        for q in d.get("element_loads", []):
            if q["elem_id"] == el_id and q.get("tipo", "costante") != "costante":
                raise ValueError(f"L'asta {el_id} ha già un carico triangolare: non può essere divisa senza alterarlo. "
                                 "Elimina quel carico, dividi l'asta e reinseriscilo sul tratto voluto.")
        nodi = {n["id"]: n for n in d["nodes"]}
        n1, n2 = nodi[el["n1"]], nodi[el["n2"]]
        x = float(f"{n1['x'] + t * (n2['x'] - n1['x']):.12g}")
        y = float(f"{n1['y'] + t * (n2['y'] - n1['y']):.12g}")
        nid = max([int(n["id"]) for n in d["nodes"]] + [self.node_id_counter - 1]) + 1
        eid = max([int(e["id"]) for e in d["elements"]] + [self.elem_id_counter - 1]) + 1
        d["nodes"].append({"id": nid, "x": x, "y": y})
        vecchio_n2 = el["n2"]
        el["n2"] = nid
        pos = d["elements"].index(el)
        d["elements"].insert(pos + 1, {"id": eid, "n1": nid, "n2": vecchio_n2})
        for chiave in ("element_loads", "thermal_loads"):
            copie = [dict(c, elem_id=eid) for c in d.get(chiave, []) if c["elem_id"] == el_id]
            d.setdefault(chiave, []).extend(copie)
        for c in d.get("internal_hinges", []):
            if c.get("elems") and int(c["node_id"]) == int(vecchio_n2):
                c["elems"] = [eid if int(x_) == int(el_id) else x_ for x_ in c["elems"]]
            if c.get("carichi_su") not in (None, "") and int(c["node_id"]) == int(vecchio_n2) and int(c["carichi_su"]) == int(el_id):
                c["carichi_su"] = eid
        self.node_id_counter = nid + 1
        self.elem_id_counter = eid + 1
        return eid

    def _ritaglia_tratto(self, d, el_id, t0, t1, eps=1e-9):
        """Rende il tratto [t0, t1] (frazioni della lunghezza di el_id) un'asta a se' stante; ritorna il suo id.
        Se il tratto e' tutta l'asta non cambia nulla."""
        id_corr, lo = el_id, 0.0       # asta che contiene la parte ancora da trattare, dalla frazione `lo` a 1
        id_tratto = el_id
        if t0 > eps:
            id_corr = self._spezza_asta(d, el_id, t0)
            lo = t0
            id_tratto = id_corr
        if t1 < 1 - eps:
            rel = (t1 - lo) / (1.0 - lo)
            self._spezza_asta(d, id_corr, rel)
            id_tratto = id_corr
        return id_tratto

    def aggiorna_anteprima_json(self):
        try:
            dati = json.loads(self.text_json.toPlainText())
            scelta_precedente = self.combo_incognita.currentData()   # per non perdere la selezione a ogni tasto
            self.combo_incognita.clear()

            # 1. Incognite esterne (Vincoli a terra). 'sup_idx' distingue piu' pendoli sullo stesso nodo.
            for i_s, v in enumerate(dati.get('supports', [])):
                n_id = v['node_id']
                if v.get('fix_u'): self.combo_incognita.addItem(f"Reazione Orizz. (Fx) al Nodo {n_id}", {'node': n_id, 'type': 'Fx'})
                if v.get('fix_v'): self.combo_incognita.addItem(f"Reazione Vert. (Fy) al Nodo {n_id}", {'node': n_id, 'type': 'Fy'})
                if v.get('fix_w'): self.combo_incognita.addItem(f"Momento Reattivo (M) al Nodo {n_id}", {'node': n_id, 'type': 'M'})
                if v.get('rotated_roller'):
                    if v.get('pendolo') or v.get('bipendolo'):
                        ang_n = float(v['axis_angle']) - 90 if 'axis_angle' in v else float(v.get('angle', 0))
                        nome_r = "Reazione del pendolo (Rn)" if v.get('pendolo') else "Reazione del doppio pendolo inclinato (Rn)"
                        etic = f"{nome_r}, asse {ang_n + 90:g}°, al Nodo {n_id}"
                    else:
                        ang_n = float(v.get('angle', 0))
                        etic = f"Reazione Normale (Rn) al Nodo {n_id}"
                    self.combo_incognita.addItem(etic, {'node': n_id, 'type': 'Fn', 'angle': ang_n, 'sup_idx': i_s})
                
            # Nota: l'incognita "M_int" (momento in cerniera) e' stata tolta. In una cerniera interna il
            # momento e' gia' nullo per definizione, quindi non puo' essere una incognita iperstatica.

            # Ripristina la selezione precedente se esiste ancora
            for i in range(self.combo_incognita.count()):
                if self.combo_incognita.itemData(i) == scelta_precedente:
                    self.combo_incognita.setCurrentIndex(i)
                    break

        except Exception:
            self.combo_incognita.clear()
            return

        # Modifica manuale del testo: allinea builder, pannello, riepilogo e anteprima al JSON digitato.
        # (Durante la digitazione il JSON puo' essere valido ma incompleto: gli errori si ignorano, il testo resta com'e'.)
        if not self._sync_in_corso and isinstance(dati, dict):
            try:
                self.sync_builder_da_json(dati, riscrivi_json=False)
            except Exception:
                pass

    def run_solver(self):
        try: 
            parsed_data = json.loads(self.text_json.toPlainText())
        except Exception as e: 
            QMessageBox.critical(self, "Errore JSON", f"JSON non valido:\n{e}")
            return
        # Il JSON e' la fonte di verita' (pannello e JSON sono sincronizzati in entrambi i versi).
        # Il pannello compila solo le chiavi mancanti: niente piu' sovrascritture silenziose dei valori scritti a mano.
        if not isinstance(parsed_data, dict):
            QMessageBox.critical(self, "Errore JSON", "Il JSON deve essere un oggetto (graffe { }).")
            return
        da_pannello = self._proprieta_da_ui({})
        for gruppo in ("material", "section"):
            dest = parsed_data.setdefault(gruppo, {})
            for k, v in da_pannello.get(gruppo, {}).items():
                dest.setdefault(k, v)
        if "pendolo_section" not in parsed_data:
            self._pendolo_da_ui(parsed_data)
        
        incognita = self.combo_incognita.currentData()
        if not incognita: 
            QMessageBox.warning(self, "Attenzione", "Seleziona un'incognita valida.")
            return

        self.btn_solve.setText("⏳ Calcolo in corso...")
        self.btn_solve.setEnabled(False)
        self.btn_pdf.setEnabled(False)
        self.html_view.setHtml("<h2>Elaborazione in corso...</h2>")
        
        # Ogni calcolo scrive in una cartella nuova: QTextBrowser tiene in cache le immagini per nome/percorso e, con nomi
        # fissi (s0_libero_1.png...), mostrerebbe i grafici del calcolo precedente. Le cartelle vecchie vengono eliminate.
        self._run_id += 1
        self.run_dir = os.path.join(self.output_dir, f"run_{int(time.time())}_{self._run_id}")
        os.makedirs(self.run_dir, exist_ok=True)
        for nome in os.listdir(self.output_dir):
            if nome.startswith("run_") and os.path.join(self.output_dir, nome) != self.run_dir:
                shutil.rmtree(os.path.join(self.output_dir, nome), ignore_errors=True)
        self.thread = SolverThread(parsed_data, self.run_dir, incognita)
        self.thread.finished_signal.connect(self.on_solver_finished)
        self.thread.error_signal.connect(self.on_solver_error)
        self.thread.start()

    def on_solver_finished(self, test_report, img_s1, img_diag):
        self.btn_solve.setText("Risolvi Struttura")
        self.btn_solve.setEnabled(True)

        # Errore bloccante / limite del solutore: il backend restituisce solo il report (box rosso), senza immagini
        if not img_s1 or not img_diag:
            self.html_view.setHtml(f"<h1>Svolgimento Metodo delle Forze (PLV)</h1>{test_report}")
            self.btn_pdf.setEnabled(False)
            self.right_tabs.setCurrentWidget(self.tab_report)
            return

        run_dir = getattr(self, "run_dir", self.output_dir)
        url_s1 = QUrl.fromLocalFile(os.path.join(run_dir, img_s1)).toString()
        url_diag = QUrl.fromLocalFile(os.path.join(run_dir, img_diag)).toString()
        
        html_finale = f"""
        <h1>Svolgimento Metodo delle Forze (PLV)</h1>
        <h2>1. Schema Strutturale (S1)</h2>
        <div style="text-align: center;"><img src="{url_s1}" width="700"></div>
        {test_report}
        <h2>11. Diagrammi delle Sollecitazioni Finali (M, T, N)</h2>
        <div style="text-align: center;"><img src="{url_diag}" width="750"></div>
        """
        self.html_view.setHtml(html_finale)
        self.btn_pdf.setEnabled(True)

    def on_solver_error(self, error_msg):
        self.btn_solve.setText("Risolvi Struttura")
        self.btn_solve.setEnabled(True)
        self.btn_pdf.setEnabled(False)
        QMessageBox.critical(self, "Errore di Calcolo nel Backend", f"<pre>{error_msg}</pre>")
        self.html_view.setHtml(f"<h2 style='color:red;'>Errore durante l'esecuzione:</h2><pre>{error_msg}</pre>")

    def _disegna_selezione(self, nodes):
        """Alone giallo sull'elemento selezionato: cosi' si vede subito quale riga si sta per modificare o eliminare."""
        if not self._sel:
            return
        cat, i = self._sel
        d = self.builder_data
        try:
            nodo = None
            asta = None
            if cat == "nodo": nodo = d["nodes"][i]["id"]
            elif cat in ("vincolo", "cedimento"): nodo = d["supports"][i]["node_id"]
            elif cat == "carico_nodo": nodo = d["point_loads"][i]["node_id"]
            elif cat == "cerniera": nodo = d["internal_hinges"][i]["node_id"]
            elif cat == "asta": asta = d["elements"][i]["id"]
            elif cat == "carico_distr": asta = d["element_loads"][i]["elem_id"]
            elif cat == "termica": asta = d["thermal_loads"][i]["elem_id"]
            if nodo in nodes:
                self.ax.plot(*nodes[nodo], marker='o', markersize=34, markerfacecolor='none',
                             markeredgecolor='#ffc107', markeredgewidth=4, zorder=1)
            if asta is not None:
                el = next((e for e in d["elements"] if e["id"] == asta), None)
                if el and el["n1"] in nodes and el["n2"] in nodes:
                    self.ax.plot([nodes[el["n1"]][0], nodes[el["n2"]][0]], [nodes[el["n1"]][1], nodes[el["n2"]][1]],
                                 color='#ffc107', linewidth=11, alpha=0.55, solid_capstyle='round', zorder=1)
        except (IndexError, KeyError):
            pass

    # ------------------------------------------------------------------
    # SIMBOLI DEI VINCOLI "COME NELLE PROVE": ogni tipo ha la sua figura, orientata dalla parte opposta alla struttura
    # ------------------------------------------------------------------
    @staticmethod
    def _sistema_vincolo(x, y, d):
        """Sistema locale del simbolo: P(a, b) = punto a distanza `a` dal nodo lungo d (verso la terra) e `b` lungo la normale."""
        t = (-d[1], d[0])
        return lambda a, b: (x + d[0] * a + t[0] * b, y + d[1] * a + t[1] * b)

    def _terra(self, P, a, mezza, u, colore, lw=2.5):
        """Linea di terra a distanza `a` dal nodo, larga 2*mezza, con il tratteggio obliquo verso l'esterno."""
        self.ax.plot(*zip(P(a, -mezza), P(a, mezza)), color=colore, linewidth=lw, zorder=3, solid_capstyle='butt')
        n = 6
        for i in range(n + 1):
            b = -mezza + 2 * mezza * i / n
            self.ax.plot(*zip(P(a, b), P(a + 0.25 * u, b - 0.25 * u)), color=colore, linewidth=1.2, zorder=3)

    def _verso_struttura(self, nid, nodes):
        """Somma dei versori dal nodo verso le aste che vi convergono (None se non c'e' nessuna asta)."""
        sx = sy = 0.0
        x0, y0 = nodes[nid]
        for el in self.builder_data.get('elements', []):
            if nid not in (el['n1'], el['n2']):
                continue
            altro = el['n2'] if el['n1'] == nid else el['n1']
            if altro not in nodes:
                continue
            dx, dy = nodes[altro][0] - x0, nodes[altro][1] - y0
            L = math.hypot(dx, dy)
            if L > 0:
                sx += dx / L
                sy += dy / L
        return (sx, sy) if math.hypot(sx, sy) > 1e-9 else None

    def _dir_terra(self, nid, nodes, modo, ang_superficie=0.0):
        """Versore dal nodo verso la terra, sempre dalla parte opposta alle aste.
        modo: 'libero' (si sceglie l'asse dominante: incastro, cerniera), 'verticale' / 'orizzontale' (asse imposto dal tipo di
        vincolo; con aste dai due lati vince il basso / la sinistra), 'normale' (carrello inclinato: normale alla superficie)."""
        m = self._verso_struttura(nid, nodes)
        mx, my = (-m[0], -m[1]) if m else (0.0, -1.0)
        eps = 1e-9
        vert = (0.0, 1.0 if my > eps else -1.0)
        oriz = (1.0 if mx > eps else -1.0, 0.0)
        if modo == 'verticale':
            return vert
        if modo == 'orizzontale':
            return oriz
        if modo == 'normale':
            a = math.radians(ang_superficie + 90.0)
            nx, ny = math.cos(a), math.sin(a)
            dot = nx * mx + ny * my
            if dot > eps: return (nx, ny)
            if dot < -eps: return (-nx, -ny)
            if ny < -eps: return (nx, ny)
            if ny > eps: return (-nx, -ny)
            return (nx, ny) if nx < 0 else (-nx, -ny)
        return oriz if abs(mx) > abs(my) + eps else vert

    def _disegna_vincolo(self, sup, x, y, u, nodes):
        """Disegna il simbolo del vincolo a terra: incastro (muro tratteggiato), cerniera (triangolo), carrello (triangolo su
        rotelle), bipendolo (due bielle parallele), carrello inclinato (ruotato con la superficie)."""
        nid = sup['node_id']
        colore = '#dc3545'
        fu, fv, fw = bool(sup.get('fix_u')), bool(sup.get('fix_v')), bool(sup.get('fix_w'))

        def triangolo(P, h, w):
            self.ax.add_patch(PoligonoMpl([P(0, 0), P(h, -w), P(h, w)], closed=True, facecolor='white',
                                          edgecolor=colore, linewidth=2, zorder=3))

        if sup.get('rotated_roller') or (fv and not fu and not fw) or (fu and not fv and not fw):
            # CARRELLO: triangolo + due rotelle + terra. Inclinato: la terra e' parallela alla superficie (normale = angle + 90)
            if sup.get('rotated_roller'):
                d = self._dir_terra(nid, nodes, 'normale', float(sup.get('angle', 0)))
            else:
                d = self._dir_terra(nid, nodes, 'verticale' if fv else 'orizzontale')
            P = self._sistema_vincolo(x, y, d)
            h, w, r = 0.50 * u, 0.38 * u, 0.10 * u
            triangolo(P, h, w)
            for b in (-0.5 * w, 0.5 * w):
                cx, cy = P(h + r, b)
                self.ax.add_patch(CerchioMpl((cx, cy), r, facecolor='white', edgecolor=colore, linewidth=1.8, zorder=3))
            self._terra(P, h + 2 * r, 1.3 * w, u, colore)
        elif fu and fv and fw:
            # INCASTRO: muro tratteggiato sul lato opposto alle aste
            d = self._dir_terra(nid, nodes, 'libero')
            P = self._sistema_vincolo(x, y, d)
            self._terra(P, 0.0, 0.65 * u, u, colore, lw=4)
        elif fu and fv:
            # CERNIERA: triangolo con il vertice nel nodo e la base a terra
            d = self._dir_terra(nid, nodes, 'libero')
            P = self._sistema_vincolo(x, y, d)
            h, w = 0.55 * u, 0.42 * u
            triangolo(P, h, w)
            self._terra(P, h, 1.3 * w, u, colore)
        elif (fv and fw) or (fu and fw):
            # BIPENDOLO (glifo): due bielle parallele con cerniere agli estremi. Scorre in orizzontale (blocca Y e rotazione)
            # => bielle verticali; scorre in verticale (blocca X e rotazione) => bielle orizzontali.
            d = self._dir_terra(nid, nodes, 'verticale' if fv else 'orizzontale')
            P = self._sistema_vincolo(x, y, d)
            Lr, off = 0.95 * u, 0.30 * u
            self.ax.plot(*zip(P(0, -off), P(0, off)), color=colore, linewidth=3, zorder=3)      # attacco rigido alla struttura
            for b in (-off, off):
                self.ax.plot(*zip(P(0, b), P(Lr, b)), color=colore, linewidth=2.5, zorder=3)
                for a_ in (0.0, Lr):
                    px, py = P(a_, b)
                    self.ax.plot([px], [py], marker='o', markersize=6, markerfacecolor='white', markeredgecolor=colore, zorder=4)
            self._terra(P, Lr, 0.6 * u, u, colore)
        else:
            # vincolo non classificabile (es. solo rotazione bloccata): segnaposto quadrato
            self.ax.plot(x, y, marker='s', color=colore, markersize=14, markerfacecolor='none', zorder=3)

    def _disegna_pendolo(self, sup, x, y, u=0.12):
        """Pendolo = biella cerniera-cerniera lungo `axis_angle` (nodo -> vincolo a terra); doppio pendolo = due bielle parallele."""
        ang = math.radians(angolo_vincolo(sup))
        ux, uy = math.cos(ang), math.sin(ang)
        nx, ny = -uy, ux
        if sup.get("pend_dx") is not None and sup.get("pend_dy") is not None:
            L = math.hypot(float(sup["pend_dx"]), float(sup["pend_dy"]))
        else:
            L = 0.9
        L = max(L, 0.3)
        colore = '#6f42c1'
        offs = (-0.12, 0.12) if sup.get("bipendolo") else (0.0,)
        for o in offs:
            xa, ya, xb, yb = x + nx * o, y + ny * o, x + nx * o + ux * L, y + ny * o + uy * L
            self.ax.plot([xa, xb], [ya, yb], color=colore, linewidth=3, zorder=3)
            self.ax.plot([xa, xb], [ya, yb], marker='o', linestyle='none', markerfacecolor='white',
                         markeredgecolor=colore, markersize=7, zorder=4)
        # terra tratteggiata all'estremo della/e biella/e
        self._terra(self._sistema_vincolo(x, y, (ux, uy)), L, 0.4 * u if sup.get("bipendolo") else 0.3 * u, u, '#343a40')
        etichetta = ("2P " if sup.get("bipendolo") else "P ") + f"{angolo_vincolo(sup):g}°"
        self.ax.text(x + ux * L * 0.5 + nx * 0.25, y + uy * L * 0.5 + ny * 0.25, etichetta, color=colore, fontsize=9,
                     weight='bold', ha='center', va='center', bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7), zorder=6)

    def update_live_preview(self):
        """Disegna la struttura sul canvas di Matplotlib in tempo reale."""
        self.ax.clear()
        self.ax.set_title("Costruzione in tempo reale", fontsize=12, weight='bold', color='#495057')
        
        # Prepariamo un dizionario per trovare rapidamente le coordinate dei nodi
        nodes = {n['id']: (n['x'], n['y']) for n in self.builder_data.get('nodes', [])}
        
        if not nodes:
            self.ax.text(0.5, 0.5, "Aggiungi un nodo per iniziare", 
                         ha='center', va='center', color='gray', fontsize=12)
            self.ax.set_xlim(0, 1)
            self.ax.set_ylim(0, 1)
            self.ax.axis('off')
            self.canvas.draw()
            return
            
        # Scala di disegno: dipende SOLO dalla struttura (nodi e bielle), non dai carichi. Cosi' frecce, carichi distribuiti,
        # momenti e cedimenti hanno dimensioni proporzionate e aggiungere un carico non cambia lo zoom dell'anteprima.
        _estremi = list(nodes.values())
        for _s in self.builder_data.get('supports', []):
            if (_s.get('pendolo') or _s.get('bipendolo')) and _s['node_id'] in nodes:
                _a = math.radians(angolo_vincolo(_s))
                _Lp = max(math.hypot(float(_s['pend_dx']), float(_s['pend_dy'])) if _s.get('pend_dx') is not None and _s.get('pend_dy') is not None else 0.9, 0.3)
                _px, _py = nodes[_s['node_id']]
                _estremi.append((_px + _Lp * math.cos(_a), _py + _Lp * math.sin(_a)))
        _xs = [p[0] for p in _estremi]
        _ys = [p[1] for p in _estremi]
        _x0, _x1, _y0, _y1 = min(_xs), max(_xs), min(_ys), max(_ys)
        span = max(_x1 - _x0, _y1 - _y0, 1.0)
        u = 0.12 * span

        # 1. Disegna le aste (linee nere doppie)
        for el in self.builder_data.get('elements', []):
            n1, n2 = el['n1'], el['n2']
            if n1 in nodes and n2 in nodes:
                x_vals = [nodes[n1][0], nodes[n2][0]]
                y_vals = [nodes[n1][1], nodes[n2][1]]
                self.ax.plot(x_vals, y_vals, color='#212529', linewidth=3, zorder=2)
        
        # 1-bis. Evidenzia l'elemento selezionato nel riepilogo (alone giallo su nodo o asta)
        self._disegna_selezione(nodes)

        # 2. Disegna i vincoli (triangoli/quadrati rossi; pendoli = bielle viola lungo il loro asse)
        for sup in self.builder_data.get('supports', []):
            nid = sup['node_id']
            if nid in nodes:
                x, y = nodes[nid]
                if sup.get('pendolo') or sup.get('bipendolo'):
                    self._disegna_pendolo(sup, x, y, u)
                    continue
                self._disegna_vincolo(sup, x, y, u, nodes)
        
        # 3. Disegna i nodi (cerchietti bianchi con bordo blu)
        for nid, (x, y) in nodes.items():
            self.ax.plot(x, y, marker='o', color='white', markeredgecolor='#0d6efd', 
                         markeredgewidth=2, markersize=16, zorder=4)
            self.ax.text(x, y, str(nid), color='#000000', fontsize=9, 
                         ha='center', va='center', weight='bold', zorder=5)
            
        # 3-bis. Cerniere interne: anello arancione attorno al nodo
        for cern in self.builder_data.get('internal_hinges', []):
            if cern['node_id'] in nodes:
                hx, hy = nodes[cern['node_id']]
                self.ax.plot(hx, hy, marker='o', markersize=26, markerfacecolor='none',
                             markeredgecolor='#fd7e14', markeredgewidth=2.5, zorder=4)

        # 4. Disegna i carichi concentrati (frecce magenta)
        for load in self.builder_data.get('point_loads', []):
            nid = load['node_id']
            if nid in nodes:
                x, y = nodes[nid]
                fx, fy, mom = load.get('Fx', 0), load.get('Fy', 0), load.get('M', 0)
                if fx != 0:
                    dx = -u if fx > 0 else u
                    self.ax.annotate('', xy=(x, y), xytext=(x + dx, y),
                                     arrowprops=dict(facecolor='magenta', edgecolor='magenta', width=2.5, headwidth=9), zorder=6)
                if fy != 0:
                    dy = -u if fy > 0 else u
                    self.ax.annotate('', xy=(x, y), xytext=(x, y + dy),
                                     arrowprops=dict(facecolor='magenta', edgecolor='magenta', width=2.5, headwidth=9), zorder=6)
                if mom != 0:
                    # freccia circolare: antioraria se M > 0, oraria se M < 0 (stessa convenzione del solutore)
                    r = 0.55 * u
                    t0, t1 = math.radians(40), math.radians(320)
                    th = [t0 + (t1 - t0) * k / 40 for k in range(41)]
                    if mom < 0:
                        th = th[::-1]
                    ax_ = [x + r * math.cos(t) for t in th]
                    ay_ = [y + r * math.sin(t) for t in th]
                    self.ax.plot(ax_[:-3], ay_[:-3], color='magenta', linewidth=2.5, zorder=6)
                    self.ax.annotate('', xy=(ax_[-1], ay_[-1]), xytext=(ax_[-4], ay_[-4]),
                                     arrowprops=dict(arrowstyle='-|>', color='magenta', lw=2.5, mutation_scale=16), zorder=6)
                    self.ax.text(x + r * 1.15, y + r * 1.15, f"M={mom:g}", color='magenta', fontsize=9, weight='bold',
                                 ha='left', va='bottom', bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7), zorder=7)
                    
        # 5. Disegna carichi distribuiti (Poligono geometrico + frecce interne)
        for q_load in self.builder_data.get('element_loads', []):
            el_id = q_load['elem_id']
            tipo_q = q_load.get('tipo', 'costante')
            q_t, q_xx, q_yy = q_load.get('q_trans', 0), q_load.get('q_x', 0), q_load.get('q_y', 0)
            if q_t != 0 or q_xx != 0 or q_yy != 0:
                asta = next((e for e in self.builder_data.get('elements', []) if e['id'] == el_id), None)
                if asta and asta['n1'] in nodes and asta['n2'] in nodes:
                    x1, y1 = nodes[asta['n1']]
                    x2, y2 = nodes[asta['n2']]

                    dx_asta = x2 - x1
                    dy_asta = y2 - y1
                    L = (dx_asta**2 + dy_asta**2)**0.5
                    if L == 0: continue

                    ux, uy = dx_asta / L, dy_asta / L
                    # Direzione REALE del carico, con la stessa convenzione del solutore (q_trans>0 verso la normale destra)
                    wx = q_t * uy + q_xx
                    wy = -q_t * ux + q_yy
                    wmod = (wx**2 + wy**2)**0.5
                    if wmod == 0: continue

                    plot_scale = u
                    offset_x = wx / wmod * plot_scale
                    offset_y = wy / wmod * plot_scale

                    if tipo_q == 'costante':
                        poly_x = [x1, x1 - offset_x, x2 - offset_x, x2, x1]
                        poly_y = [y1, y1 - offset_y, y2 - offset_y, y2, y1]
                    elif tipo_q == 'triangolare_crescente':
                        poly_x = [x1, x2 - offset_x, x2, x1]
                        poly_y = [y1, y2 - offset_y, y2, y1]
                    else:
                        poly_x = [x1, x1 - offset_x, x2, x1]
                        poly_y = [y1, y1 - offset_y, y2, y1]

                    self.ax.fill(poly_x, poly_y, color='deepskyblue', alpha=0.15, zorder=5)
                    self.ax.plot(poly_x, poly_y, color='deepskyblue', linewidth=2, zorder=6)

                    num_arrows = 9
                    for i in range(num_arrows):
                        t = i / (num_arrows - 1)
                        px = x1 + t * dx_asta
                        py = y1 + t * dy_asta

                        if tipo_q == 'costante': intensity = 1.0
                        elif tipo_q == 'triangolare_crescente': intensity = t
                        else: intensity = 1.0 - t

                        tail_x = px - offset_x * intensity
                        tail_y = py - offset_y * intensity

                        if intensity > 0.05:
                            self.ax.annotate('', xy=(px, py), xytext=(tail_x, tail_y),
                                             arrowprops=dict(facecolor='deepskyblue', edgecolor='deepskyblue',
                                                             width=1.5, headwidth=6), zorder=6)

                    xm, ym = x1 + 0.5 * dx_asta, y1 + 0.5 * dy_asta
                    text_offset_x = - offset_x * 1.4
                    text_offset_y = - offset_y * 1.4

                    simbolo = "q" if tipo_q == 'costante' else ("q(s) ◢" if tipo_q == 'triangolare_crescente' else "q(s) ◣")
                    parti = []
                    if q_t != 0: parti.append(f"q⊥ = {q_t:g}")
                    if q_xx != 0: parti.append(f"qx = {q_xx:g}")
                    if q_yy != 0: parti.append(f"qy = {q_yy:g}")
                    self.ax.text(xm + text_offset_x, ym + text_offset_y, f"{simbolo}: " + ", ".join(parti),
                                 color='darkblue', fontsize=10, ha='center', va='center', weight='bold',
                                 bbox=dict(boxstyle='round,pad=0.2', fc='lightcyan', ec='deepskyblue', alpha=0.9), zorder=7)

        # 6. DISEGNA CEDIMENTI VINCOLARI E ROTAZIONI IMPOSTE
        for sup in self.builder_data.get('supports', []):
            nid = sup['node_id']
            if nid in nodes:
                x, y = nodes[nid]
                
                dx = sup.get('delta_u', 0)
                dy = sup.get('delta_v', 0)
                dn = sup.get('delta_n', 0)
                dw = sup.get('delta_w', 0)
                
                color_ced = 'darkorange'
                arrow_len = u
                
                # Cedimento Orizzontale (Δx)
                if dx != 0:
                    dir_x = 1 if dx > 0 else -1
                    self.ax.annotate('', xy=(x + dir_x * arrow_len, y - 0.3), xytext=(x, y - 0.3),
                                     arrowprops=dict(edgecolor=color_ced, arrowstyle='->', lw=2.5, linestyle='--'), zorder=6)
                    self.ax.text(x + dir_x * arrow_len * 0.5, y - 0.5, f"Δx={dx*1000:g}mm", 
                                 color=color_ced, fontsize=9, ha='center', weight='bold',
                                 bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))
                                 
                # Cedimento Verticale (Δy)
                if dy != 0:
                    dir_y = 1 if dy > 0 else -1
                    self.ax.annotate('', xy=(x + 0.3, y + dir_y * arrow_len), xytext=(x + 0.3, y),
                                     arrowprops=dict(edgecolor=color_ced, arrowstyle='->', lw=2.5, linestyle='--'), zorder=6)
                    self.ax.text(x + 0.45, y + dir_y * arrow_len * 0.5, f"Δy={dy*1000:g}mm", 
                                 color=color_ced, fontsize=9, ha='left', va='center', weight='bold',
                                 bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))
                                 
                # Cedimento Inclinato (Δn)
                if dn != 0:
                    ang_rad = math.radians(sup.get('angle', 0) + 90)
                    dir_nx = math.cos(ang_rad) * (1 if dn > 0 else -1)
                    dir_ny = math.sin(ang_rad) * (1 if dn > 0 else -1)
                    self.ax.annotate('', xy=(x + dir_nx * arrow_len, y + dir_ny * arrow_len), xytext=(x, y),
                                     arrowprops=dict(edgecolor=color_ced, arrowstyle='->', lw=2.5, linestyle='--'), zorder=6)
                    self.ax.text(x + dir_nx * arrow_len * 1.2, y + dir_ny * arrow_len * 1.2, f"Δn={dn*1000:g}mm", 
                                 color=color_ced, fontsize=9, ha='center', va='center', weight='bold',
                                 bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))
                                 
                # Rotazione Imposta (Δφ)
                if dw != 0:
                    rad_curve = 0.4 if dw > 0 else -0.4
                    self.ax.annotate('', xy=(x - 0.4, y + 0.4), xytext=(x + 0.4, y + 0.4),
                                     arrowprops=dict(arrowstyle="->", color=color_ced, lw=2.5, linestyle='--',
                                                     connectionstyle=f"arc3,rad={rad_curve}"), zorder=6)
                    self.ax.text(x, y + 0.7, f"Δφ={dw:g}rad", 
                                 color=color_ced, fontsize=9, ha='center', weight='bold',
                                 bbox=dict(boxstyle='round,pad=0.1', fc='white', ec=color_ced, alpha=0.8))

        # 7. DISEGNA VARIAZIONI TERMICHE (Diagramma a farfalla)
        for t_load in self.builder_data.get('thermal_loads', []):
            el_id = t_load['elem_id']
            t_pos = t_load.get('T_pos', 0)
            t_neg = t_load.get('T_neg', 0)
            
            asta = next((e for e in self.builder_data.get('elements', []) if e['id'] == el_id), None)
            if asta and asta['n1'] in nodes and asta['n2'] in nodes:
                x1, y1 = nodes[asta['n1']]
                x2, y2 = nodes[asta['n2']]
                
                # Troviamo il punto medio dell'asta
                xm, ym = (x1 + x2) / 2, (y1 + y2) / 2
                dx, dy = x2 - x1, y2 - y1
                L = (dx**2 + dy**2)**0.5
                if L == 0: continue
                
                # Versori asse e normale
                ux, uy = dx / L, dy / L
                nx, ny = -uy, ux
                
                # Dimensione della "farfalla"
                h_t = 0.75 * u  # Altezza lungo la normale
                w_t = 0.5 * u   # Larghezza lungo l'asse dell'asta
                
                # Triangolo Fibra Sinistra (+T)
                if t_pos != 0:
                    p_pos_x, p_pos_y = xm + nx * h_t, ym + ny * h_t
                    p_base_pos_x, p_base_pos_y = xm + ux * w_t, ym + uy * w_t
                    
                    poly_pos_x = [xm, p_pos_x, p_base_pos_x, xm]
                    poly_pos_y = [ym, p_pos_y, p_base_pos_y, ym]
                    
                    self.ax.fill(poly_pos_x, poly_pos_y, color='tomato', alpha=0.4, zorder=5)
                    self.ax.plot(poly_pos_x, poly_pos_y, color='red', lw=1.5, zorder=6)
                    self.ax.text(p_pos_x + nx*0.2, p_pos_y + ny*0.2, f"{t_pos:g}°C", 
                                 color='red', fontsize=9, weight='bold', ha='center', va='center',
                                 bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))
                                 
                # Triangolo Fibra Destra (-T)
                if t_neg != 0:
                    p_neg_x, p_neg_y = xm - nx * h_t, ym - ny * h_t
                    p_base_neg_x, p_base_neg_y = xm - ux * w_t, ym - uy * w_t
                    
                    poly_neg_x = [xm, p_neg_x, p_base_neg_x, xm]
                    poly_neg_y = [ym, p_neg_y, p_base_neg_y, ym]
                    
                    self.ax.fill(poly_neg_x, poly_neg_y, color='cornflowerblue', alpha=0.4, zorder=5)
                    self.ax.plot(poly_neg_x, poly_neg_y, color='blue', lw=1.5, zorder=6)
                    self.ax.text(p_neg_x - nx*0.2, p_neg_y - ny*0.2, f"{t_neg:g}°C", 
                                 color='blue', fontsize=9, weight='bold', ha='center', va='center',
                                 bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.7))
            
        # Limiti fissati dalla sola struttura (+ margine per ospitare carichi, momenti e cedimenti): i carichi non
        # entrano nel calcolo dello zoom, quindi la figura non si rimpicciolisce quando se ne aggiungono.
        marg = 0.30 * span
        self.ax.set_xlim(_x0 - marg, _x1 + marg)
        self.ax.set_ylim(_y0 - marg, _y1 + marg)
        self.ax.set_aspect('equal', adjustable='datalim')
        self.ax.grid(True, linestyle='--', alpha=0.5)
        self.canvas.draw()

    def esporta_pdf(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Salva Report PDF", "Report_Struttura.pdf", "PDF Files (*.pdf)")
        if not file_path: return
            
        try:
            printer = QPrinter(QPrinter.PrinterMode.ScreenResolution)
            printer.setOutputFormat(QPrinter.OutputFormat.PdfFormat)
            printer.setOutputFileName(file_path)
            
            layout = QPageLayout()
            layout.setPageSize(QPageSize(QPageSize.PageSizeId.A2))
            layout.setOrientation(QPageLayout.Orientation.Landscape)
            layout.setMargins(QMarginsF(15, 15, 15, 15))
            printer.setPageLayout(layout)
            
            temp_doc = self.html_view.document().clone()
            temp_doc.setPageSize(printer.pageRect(QPrinter.Unit.Point).size())
            temp_doc.print(printer)
            
            QMessageBox.information(self, "Successo", "Report esportato in PDF con successo!")
        except Exception as e:
            QMessageBox.critical(self, "Errore", f"Errore durante il salvataggio del PDF:\n{str(e)}")

def _gestore_eccezioni(tipo, valore, tb):
    """PyQt6 termina il processo se uno slot solleva un'eccezione non gestita: la si mostra invece in una finestra."""
    testo = "".join(traceback.format_exception(tipo, valore, tb))
    print(testo, file=sys.stderr)
    try:
        QMessageBox.critical(None, "Errore inatteso", f"<pre>{html.escape(testo)}</pre>")
    except Exception:
        pass


if __name__ == "__main__":
    sys.excepthook = _gestore_eccezioni
    app = QApplication(sys.argv)
    window = UIStrutture()
    window.show()
    sys.exit(app.exec())
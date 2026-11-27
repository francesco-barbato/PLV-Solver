import os
import sys
import json
import math
import struct
import threading
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import url2pathname

import sympy as sp
import numpy as np

_MONO_CSS = 'Consolas, "Courier New", Courier, monospace' if sys.platform == "win32" else '"Courier New", Courier, monospace'

def _to_file_uri(img_path):
    """Converte un percorso locale (anche Windows con spazi, accenti o drive C:\\) in URI file:// valido."""
    if not img_path:
        return ""
    s = str(img_path)
    if s.startswith("file://"):
        return s
    return Path(os.path.abspath(s)).as_uri()

def _risolvi_percorso_img(img_ref):
    """Restituisce il percorso locale OS da path o URI file://."""
    s = str(img_ref or "")
    if s.startswith("file://"):
        parsed = urlparse(s)
        path = url2pathname(parsed.path)
        if parsed.netloc:
            return f"//{parsed.netloc}{path}" if os.name != "nt" else f"\\\\{parsed.netloc}{path}"
        return path
    return s if os.path.isabs(s) else ("/" + s.lstrip("/"))

def _salva_figura_png(fig, out_dir, filename, **save_kw):
    """Salva la figura Matplotlib e la chiude sempre in sicurezza; restituisce il path normalizzato."""
    _plt, _ = _get_mpl()
    filepath = os.path.abspath(os.path.join(out_dir, filename))
    try:
        save_kw.setdefault("bbox_inches", "tight")
        fig.savefig(filepath, **save_kw)
    finally:
        _plt.close(fig)
    return filepath.replace("\\", "/")

plt = None
patches = None
_MPL_LOCK = threading.Lock()

def _get_mpl():
    """Import lazy di matplotlib (pyplot backend Agg + patches) solo al momento del disegno."""
    global plt, patches
    if plt is not None and patches is not None:
        return plt, patches
    with _MPL_LOCK:
        if plt is None or patches is None:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as _plt
            import matplotlib.patches as _patches
            plt = _plt
            patches = _patches
        return plt, patches

# =====================================================================
# FASE 1 - ARITMETICA ESATTA, PULIZIA RUMORE, ERRORI BLOCCANTI
# =====================================================================
TOL_RUMORE = 1e-10
_COSTANTI = None

ULTIMI_RISULTATI = {}   # riempito a ogni calcolo: reazioni finali, X1, MS, pendoli, residui dei controlli

class ErroreSolutore(Exception):
    """Errore bloccante: il messaggio contiene il report HTML parziale + box d'errore."""
    pass

def ex(v, default=0):
    """Converte un numero (int/float/str/SymPy) in valore ESATTO SymPy.
    4.0 -> 4 ; 0.5 -> 1/2 ; 2.8284271247461903 -> 2*sqrt(2) ; 1.7320508 -> sqrt(3)."""
    global _COSTANTI
    if v is None:
        v = default
    if isinstance(v, sp.Basic):
        return v
    if isinstance(v, (bool, np.bool_)):
        return sp.Integer(int(v))
    if isinstance(v, (int, np.integer)):
        return sp.Integer(int(v))
    v = float(v)
    r6 = round(v, 6)
    if abs(v - r6) < 1e-12:                     # decimale "pulito" (anche con rumore ~1e-16)
        return sp.Rational(repr(r6)) if r6 == r6 else sp.Integer(0)
    try:                                        # radicali notevoli (45deg, 30/60deg, arctan 2...)
        if _COSTANTI is None:
            _COSTANTI = [sp.sqrt(2), sp.sqrt(3), sp.sqrt(5)]
        r = sp.nsimplify(v, _COSTANTI, tolerance=TOL_RUMORE)
        if abs(float(r) - v) <= 1e-9 * max(1.0, abs(v)) and sp.count_ops(r) <= 8:
            return r
    except Exception:
        pass
    return sp.Rational(repr(v))

def pulisci_espr(expr, tol=TOL_RUMORE):
    """Rimuove il rumore numerico: Float < tol -> 0, altri Float -> numeri esatti."""
    expr = sp.sympify(expr)
    sost = {}
    for f in expr.atoms(sp.Float):
        sost[f] = sp.Integer(0) if abs(f) < tol else ex(float(f))
    return expr.xreplace(sost) if sost else expr

def hyp(dx, dy):
    """Lunghezza ESATTA sqrt(dx^2+dy^2) (sostituisce math.hypot)."""
    return sp.sqrt(sp.expand(ex(dx) ** 2 + ex(dy) ** 2))

def cs_ang(angolo_deg):
    """(cos, sin) ESATTI di un angolo in gradi (sostituisce math.cos/sin(math.radians))."""
    ar = ex(angolo_deg) * sp.pi / 180
    c, s = sp.cos(ar), sp.sin(ar)
    if c.has(sp.cos) or c.has(sp.sin):
        c = ex(float(c.evalf(30)))
    if s.has(sp.cos) or s.has(sp.sin):
        s = ex(float(s.evalf(30)))
    return c, s

def normalizza_dati(d):
    """Copia del JSON con tutti i valori numerici convertiti in esatti SymPy."""
    out = dict(d)
    out['nodes'] = [{**n, 'id': int(n['id']), 'x': ex(n['x']), 'y': ex(n['y'])} for n in d.get('nodes', [])]
    out['elements'] = [{**e, 'id': int(e['id']), 'n1': int(e['n1']), 'n2': int(e['n2'])} for e in d.get('elements', [])]
    sups = []
    for s in d.get('supports', []):
        s2 = dict(s); s2['node_id'] = int(s['node_id'])
        # Pendolo (biella) e doppio pendolo (bipendolo) inclinati: 'axis_angle' = direzione dell'asse delle bielle,
        # misurata da +x e orientata dal nodo verso il vincolo a terra. La reazione Rn agisce lungo l'asse
        # (Rn > 0 = la struttura e' tirata verso il vincolo a terra = TIRANTE), quindi equivale a un carrello
        # inclinato la cui normale (angle + 90) coincide con l'asse: angle = axis_angle - 90.
        if s2.get('pendolo') or s2.get('bipendolo'):
            ax_ang = ex(s2.get('axis_angle', ex(s2.get('angle', 0)) + 90))
            s2['axis_angle'] = ax_ang
            s2['rotated_roller'] = True
            s2['angle'] = ax_ang - 90
            if s2.get('bipendolo'):
                s2['fix_w'] = True
        if 'angle' in s2: s2['angle'] = ex(s2['angle'])
        for k in ('delta_u', 'delta_v', 'delta_w', 'delta_n'):
            if k in s2: s2[k] = ex(s2[k])
        sups.append(s2)
    out['supports'] = sups
    out['internal_hinges'] = [{**c, 'node_id': int(c['node_id'])} for c in d.get('internal_hinges', [])]
    out['point_loads'] = [{**p, 'node_id': int(p['node_id']), 'Fx': ex(p.get('Fx', 0)), 'Fy': ex(p.get('Fy', 0)), 'M': ex(p.get('M', 0))} for p in d.get('point_loads', [])]
    out['element_loads'] = [{**q, 'elem_id': int(q['elem_id']), 'q_trans': ex(q.get('q_trans', 0)),
                             'q_x': ex(q.get('q_x', 0)), 'q_y': ex(q.get('q_y', 0))} for q in d.get('element_loads', [])]
    th = []
    for t in d.get('thermal_loads', []):
        t2 = {**t, 'elem_id': int(t['elem_id']), 'T_pos': ex(t.get('T_pos', 0)), 'T_neg': ex(t.get('T_neg', 0))}
        if 'alpha' in t2: t2['alpha'] = ex(t2['alpha'])
        th.append(t2)
    out['thermal_loads'] = th
    mat = dict(d.get('material', {}))
    mat.setdefault('E', 2.1e11); mat.setdefault('sigma_amm', 235e6)
    for k in ('E', 'sigma_amm', 'alpha', 'MS_target', 'MS_design_target', 'EI'):
        if k in mat: mat[k] = ex(mat[k])
    out['material'] = mat
    sec = dict(d.get('section', {}))
    sec.setdefault('b', 0.2); sec.setdefault('h', 0.4)      # un solo default coerente (prima: 0.2/0.4 vs 0.12/0.15)
    for k in ('b', 'h'):
        sec[k] = ex(sec[k])
    sec.setdefault('design', 'none')          # 'none' | 'h' (b fissata) | 'b' (h fissata) | 'ratio' (h = ratio * b)
    if 'ratio' in sec: sec['ratio'] = ex(sec['ratio'])
    out['section'] = sec
    return out

def valida_dati(d):
    """Controlli di coerenza sui dati NORMALIZZATI. Restituisce l'elenco dei messaggi d'errore (vuoto = tutto ok).
    Evita risultati falsi dovuti a errori di inserimento: riferimenti a nodi/aste inesistenti, id duplicati,
    aste di lunghezza nulla, vincoli sovrapposti sullo stesso nodo."""
    err = []
    ids_n = [n['id'] for n in d.get('nodes', [])]
    ids_e = [e['id'] for e in d.get('elements', [])]
    nodi_set, el_set = set(ids_n), set(ids_e)
    if len(ids_n) != len(nodi_set): err.append("Ci sono nodi con lo stesso id.")
    if len(ids_e) != len(el_set): err.append("Ci sono aste con lo stesso id.")
    if not ids_n: err.append("Nessun nodo inserito.")
    if not ids_e: err.append("Nessuna asta inserita.")
    coord = {n['id']: (n['x'], n['y']) for n in d.get('nodes', [])}
    for e in d.get('elements', []):
        if e['n1'] not in nodi_set or e['n2'] not in nodi_set:
            err.append(f"L'asta {e['id']} usa un nodo che non esiste (N{e['n1']} &rarr; N{e['n2']}).")
        elif e['n1'] == e['n2']:
            err.append(f"L'asta {e['id']} ha lo stesso nodo come inizio e fine.")
        elif coord[e['n1']] == coord[e['n2']]:
            err.append(f"L'asta {e['id']} ha lunghezza nulla (i nodi {e['n1']} e {e['n2']} coincidono).")
    for i_a in range(len(ids_n)):
        for i_b in range(i_a + 1, len(ids_n)):
            if coord[ids_n[i_a]] == coord[ids_n[i_b]]:
                err.append(f"I nodi {ids_n[i_a]} e {ids_n[i_b]} hanno le stesse coordinate.")
    fix_per_nodo = {}
    asse_per_nodo = {}
    for s in d.get('supports', []):
        nid = s['node_id']
        if nid not in nodi_set:
            err.append(f"Un vincolo e' assegnato al nodo {nid}, che non esiste.")
            continue
        for chiave, nome in (('fix_u', 'traslazione X'), ('fix_v', 'traslazione Y'), ('fix_w', 'rotazione')):
            if s.get(chiave):
                if (nid, chiave) in fix_per_nodo:
                    err.append(f"Il nodo {nid} ha due vincoli che bloccano la stessa {nome}: vincolo sovrapposto, eliminane uno.")
                fix_per_nodo[(nid, chiave)] = True
        if s.get('rotated_roller'):
            ang_asse = float(s['angle']) + 90.0
            for altro in asse_per_nodo.get(nid, []):
                diff = (ang_asse - altro) % 180.0
                if diff < 1e-6 or abs(diff - 180.0) < 1e-6:
                    err.append(f"Il nodo {nid} ha due vincoli inclinati con la stessa direzione di reazione: sono sovrapposti.")
            asse_per_nodo.setdefault(nid, []).append(ang_asse)
    for p in d.get('point_loads', []):
        if p['node_id'] not in nodi_set: err.append(f"Un carico concentrato e' assegnato al nodo {p['node_id']}, che non esiste.")
    for q in d.get('element_loads', []):
        if q['elem_id'] not in el_set: err.append(f"Un carico distribuito e' assegnato all'asta {q['elem_id']}, che non esiste.")
    for tl in d.get('thermal_loads', []):
        if tl['elem_id'] not in el_set: err.append(f"Una variazione termica e' assegnata all'asta {tl['elem_id']}, che non esiste.")
    visti_cern = set()
    for c in d.get('internal_hinges', []):
        if c['node_id'] not in nodi_set: err.append(f"Una cerniera interna e' assegnata al nodo {c['node_id']}, che non esiste.")
        elif c['node_id'] in visti_cern:
            err.append(f"Il nodo {c['node_id']} ha piu' di una cerniera interna: ne basta una (indica nel campo elems tutte le aste incernierate), "
                       "altrimenti il numero di corpi e il grado di iperstaticita' risultano sbagliati.")
        visti_cern.add(c['node_id'])
    # connettivita': tutti i nodi devono appartenere alla stessa struttura (un nodo isolato o due strutture separate
    # falserebbero il conteggio dei corpi e quindi il grado di iperstaticita')
    if ids_n and not any(e_.startswith("L'asta") and 'non esiste' in e_ for e_ in err):
        padre = {n_: n_ for n_ in ids_n}
        def _rad(a_):
            while padre[a_] != a_:
                padre[a_] = padre[padre[a_]]
                a_ = padre[a_]
            return a_
        for e in d.get('elements', []):
            if e['n1'] in padre and e['n2'] in padre:
                padre[_rad(e['n1'])] = _rad(e['n2'])
        isolati = [n_ for n_ in ids_n if not any(n_ in (e['n1'], e['n2']) for e in d.get('elements', []))]
        if isolati and len(ids_n) > 1:
            err.append(f"Il nodo {', '.join(str(n_) for n_ in isolati)} non e' collegato a nessuna asta: eliminalo o collegalo.")
        elif len({_rad(n_) for n_ in ids_n}) > 1:
            err.append("La struttura e' formata da parti non collegate fra loro: il solutore tratta una sola struttura connessa.")
    # distorsione termica senza coefficiente di dilatazione: l'effetto sparirebbe in silenzio (alpha = 0)
    mat_ = d.get('material', {})
    for tl in d.get('thermal_loads', []):
        if tl['elem_id'] in el_set and (tl.get('T_pos', 0) != 0 or tl.get('T_neg', 0) != 0) \
                and 'alpha' not in tl and 'alpha' not in mat_:
            err.append("Sono presenti variazioni termiche ma manca il coefficiente di dilatazione: aggiungi \"alpha\" in material (es. 1.7e-5).")
            break
    return err

def _simboli_op(op):
    """(simbolo MS, simbolo sulla tensione) per l'operatore del margine di sicurezza."""
    return {">=": ("&ge;", "&le;"), ">": ("&gt;", "&lt;"), "<=": ("&le;", "&ge;"), "<": ("&lt;", "&gt;")}.get(op, ("=", "="))

def _verifica_op(MS, op, target, tol=1e-9):
    """True se MS soddisfa la richiesta 'MS op target' (=, >, >=, <, <=)."""
    if MS == float('inf'):
        return op in (">=", ">")
    if op == ">=": return MS >= target - tol
    if op == ">": return MS > target + tol
    if op == "<=": return MS <= target + tol
    if op == "<": return MS < target - tol
    return math.isclose(MS, target, abs_tol=1e-2)

def _box_errore(titolo, msg):
    return ("<div style='background-color: #f8d7da; padding: 15px; border-left: 5px solid #dc3545; margin-bottom: 20px;'>"
            f"<h4 style='color: #721c24; margin-top: 0;'>🛑 {titolo}</h4>"
            f"<p style='color: #721c24;'>{msg}</p></div>")
def _box_avviso(titolo, msg):
    return ("<div style='background-color: #fff3cd; padding: 15px; border-left: 5px solid #ffc107; margin-bottom: 20px;'>"
            f"<h4 style='color: #856404; margin-top: 0;'>⚠️ {titolo}</h4>"
            f"<p style='color: #856404;'>{msg}</p></div>")

def risolvi_sistema_lineare(equazioni, incognite, nome, report):
    """Risolve il sistema lineare in modo ESATTO. Se non esiste soluzione unica
    (sistema vuoto, labile, indeterminato, incompatibile) solleva ErroreSolutore."""
    def stop(titolo, msg):
        report.append(_box_errore(titolo, msg))
        raise ErroreSolutore("\n".join(report))
    if not incognite:
        stop(f"Sistema {nome} impossibile", "Non ci sono reazioni vincolari incognite: la struttura non ha vincoli utili.")
    A, b_ = sp.linear_eq_to_matrix([sp.sympify(e) for e in equazioni], incognite)
    rA = A.rank()
    if rA < len(incognite):
        stop(f"Sistema {nome} indeterminato (struttura labile o vincoli mal posti)",
             f"Le equazioni di equilibrio non determinano tutte le {len(incognite)} reazioni (rango {rA}). "
             "Cause tipiche: tre vincoli concorrenti o paralleli, cerniere interne allineate con i vincoli, incognita iperstatica scelta male.")
    if A.row_join(b_).rank() > rA:
        stop(f"Sistema {nome} incompatibile",
             "Le equazioni di equilibrio sono in contraddizione: la struttura è labile in una sua parte oppure una cerniera interna "
             "non è stata modellata correttamente. Nessun risultato numerico viene prodotto.")
    sol = sp.linsolve((A, b_), *incognite)
    if not sol:
        stop(f"Sistema {nome} senza soluzione", "Il sistema lineare di equilibrio non ammette soluzione.")
    vals = list(sol)[0]
    risultato = {v: sp.nsimplify(sp.simplify(x)) if not x.free_symbols else sp.simplify(x) for v, x in zip(incognite, vals)}
    if not risultato or any(x.free_symbols for x in risultato.values()):
        stop(f"Sistema {nome} non univocamente risolvibile", "Alcune reazioni restano espresse in funzione di altre incognite.")
    return risultato

def _trim_zeros_num_str(s):
    """Elimina zeri decimali finali inutili (es. '1200.00' -> '1200', '2.050' -> '2.05')."""
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return "0" if s in ("", "-", "-0", "&minus;", "&minus;0") else s


def _fmt_sci(val, sig=3, plain=False):
    """Notazione scientifica pulita con `sig` cifre significative (HTML o Unicode plain)."""
    v = float(val)
    if v == 0.0 or abs(v) < 1e-12:
        return "0"
    sig_dec = max(0, int(sig) - 1)
    mant, exp_s = f"{abs(v):.{sig_dec}e}".split("e")
    mant = _trim_zeros_num_str(mant)
    exp_i = int(exp_s)
    if exp_i == 0:
        sign_p = ("-" if plain else "&minus;") if v < 0 else ""
        return f"{sign_p}{mant}"
    if plain:
        sign = "-" if v < 0 else ""
        sup_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
        return f"{sign}{mant} × 10{str(exp_i).translate(sup_map)}"
    sign = "&minus;" if v < 0 else ""
    exp_str = f"&minus;{abs(exp_i)}" if exp_i < 0 else str(exp_i)
    return f"{sign}{mant} &times; 10<sup>{exp_str}</sup>"


def fmt_smart(val, min_sig=None, fixed_min=0.001, dec_large=2, sci_sig=3, plain=False):
    """Formatta numeri adattando la precisione all'ordine di grandezza senza troncare valori piccoli."""
    v = float(val)
    if v == 0.0 or abs(v) < 1e-12:
        return "0"
    av = abs(v)
    if av < fixed_min:
        sig_use = min_sig if min_sig is not None else sci_sig
        return _fmt_sci(v, sig=sig_use, plain=plain)
    if min_sig is not None:
        e10 = math.floor(math.log10(av))
        dec = max(dec_large if av >= 1.0 else 4, min_sig - 1 - e10)
    elif av >= 1.0:
        dec = dec_large
    else:
        dec = 4
    return _trim_zeros_num_str(f"{v:.{dec}f}")


def _fmt_sci_dispatch(val, dec=2, sig=None, plain=False):
    """Formatta in notazione scientifica (o smart per ordini intermedi) condiviso tra HTML e plain."""
    v = float(val)
    if v == 0.0 or abs(v) < 1e-12:
        return "0"
    if 0.001 <= abs(v) < 10000:
        return fmt_smart(v, min_sig=sig, dec_large=dec, plain=plain)
    sig_use = sig if sig is not None else max(1, int(dec) + 1)
    return _fmt_sci(v, sig=sig_use, plain=plain)


def fmt_sci_plain(val, dec=2, sig=None):
    """Notazione scientifica Unicode plain text (per Matplotlib) senza tag HTML."""
    return _fmt_sci_dispatch(val, dec=dec, sig=sig, plain=True)


def fmt_num(val, dec=2, plain=False):
    """Formatta un numero reale in stringa pulita evitando troncamenti su valori piccoli."""
    v = float(val)
    if v == 0.0 or abs(v) < 1e-12:
        return "0"
    if abs(v) >= 0.5 and abs(v - round(v)) < 1e-7:
        return str(int(round(v)))
    if abs(v) < 10 ** (-dec + 1):
        return fmt_smart(v, dec_large=dec, plain=plain)
    s = f"{v:.{dec}f}"
    return "0" if s.lstrip("-") in ("0", "0.0", "0.00", "0.000", "0.0000") else s


def fmt_sci_html(val, dec=2, sig=None):
    """Notazione scientifica HTML (es. 2.05 × 10¹¹) al posto di 2.05e+11; numeri medi in formato smart."""
    return _fmt_sci_dispatch(val, dec=dec, sig=sig, plain=False)

def sym_reaz_html(name):
    """Converte nomi di reazione (H_8, V_3, M_1, Rn_7, Hc_5, Vc_5, X_1) in pedici HTML senza underscore."""
    s = str(name)
    if s in ("X1", "X_1"):
        return "X<sub>1</sub>"
    parts = s.split("_")
    if len(parts) == 2:
        if parts[0] in ("Hc", "Vc"):
            return f"{parts[0][0]}<sub>c{parts[1]}</sub>"
        pref = "R" if parts[0] == "Rn" else parts[0]
        sub = f"n{parts[1]}" if parts[0] == "Rn" else parts[1]
        return f"{pref}<sub>{sub}</sub>"
    if len(parts) == 3:
        if parts[0] in ("Hc", "Vc"):
            return f"{parts[0][0]}<sub>c{parts[1]}.{parts[2]}</sub>"
        if parts[0] == "Rn":
            return f"R<sub>n{parts[1]}.{parts[2]}</sub>"
    return s.replace("_", "")

def sym_reaz_mpl(name):
    """Converte nomi di reazione in mathtext Matplotlib senza underscore."""
    s = str(name)
    if s in ("X1", "X_1"):
        return r"$X_{1}$"
    parts = s.split("_")
    if len(parts) == 2:
        if parts[0] in ("Hc", "Vc"):
            return rf"${parts[0][0]}_{{c{parts[1]}}}$"
        pref = "R" if parts[0] == "Rn" else parts[0]
        sub = f"n{parts[1]}" if parts[0] == "Rn" else parts[1]
        return rf"${pref}_{{{sub}}}$"
    if len(parts) == 3:
        if parts[0] in ("Hc", "Vc"):
            return rf"${parts[0][0]}_{{c{parts[1]}.{parts[2]}}}$"
        if parts[0] == "Rn":
            return rf"$R_{{n{parts[1]}.{parts[2]}}}$"
    return s.replace("_", "")

def _frac_html(num, den, border_col="#0b5ed7"):
    """Frazione verticale HTML compatibile con QTextBrowser/PDF (senza bordi ereditati)."""
    return (
        "<table class='frac-table' cellpadding='0' cellspacing='0' "
        "style='margin:0; border:none; width:auto; background:transparent; display:inline-table; vertical-align:middle;'>"
        f"<tr><td style='border:none; border-bottom:1px solid {border_col}; text-align:center; padding:0 3px; background:transparent; font-size:12px;'>{num}</td></tr>"
        f"<tr><td style='border:none; text-align:center; padding:0 3px; background:transparent; font-size:12px;'>{den}</td></tr>"
        "</table>"
    )

def _eq_table_html(*elements):
    """Riga di equazione orizzontale senza bordi, sicura per frazioni impilate in QTextBrowser."""
    cells = "".join(
        f"<td valign='middle' style='border:none; padding:0 3px; background:transparent; width:auto; text-align:left;'>{el}</td>"
        for el in elements
    )
    return (
        "<table class='eq-table' cellpadding='0' cellspacing='0' "
        f"style='margin:2px 0; border:none; width:auto; background:transparent; font-family:{_MONO_CSS}; font-size:14px; font-weight:bold; color:#0b5ed7;'>"
        f"<tr>{cells}</tr></table>"
    )

def _format_coeff_num(c_expr, decimale=False, stacked_frac=False):
    """Restituisce (abs_str, is_neg, is_unit, is_stacked) per un fattore numerico/radicale SymPy."""
    c_expr = sp.radsimp(sp.simplify(c_expr))
    val = float(sp.N(c_expr))
    if abs(val) < 1e-10:
        return "0", False, False, False
    is_neg = val < 0
    c_pos = -c_expr if is_neg else c_expr

    if decimale:
        s_abs = fmt_num(abs(val), 2)
        if s_abs == "0":
            return "0", False, False, False
        return s_abs, is_neg, (s_abs == "1"), False

    if c_pos.is_Integer:
        iv = int(c_pos)
        return str(iv), is_neg, (iv == 1), False

    if c_pos.is_Rational:
        p, q = int(c_pos.p), int(c_pos.q)
        if q > 10 or p >= 1000:
            s_abs = fmt_num(abs(val), 2)
            return s_abs, is_neg, (s_abs == "1"), False
        if stacked_frac:
            return _frac_html(str(p), str(q)), is_neg, False, True
        return f"<sup>{p}</sup>&frasl;<sub>{q}</sub>", is_neg, False, False

    r_num, r_rad = sp.Integer(1), None
    fattori = c_pos.as_ordered_factors()
    rad_ok = True
    for f_ in fattori:
        if f_.is_Rational:
            r_num *= f_
        elif isinstance(f_, sp.Pow) and f_.exp == sp.Rational(1, 2) and f_.base.is_Integer and r_rad is None:
            r_rad = int(f_.base)
        else:
            rad_ok = False
            break
    if rad_ok and r_rad is not None and r_num.is_Rational:
        p, q = int(r_num.p), int(r_num.q)
        if q <= 10 and p < 10000:
            rad_html = f"&radic;{r_rad}"
            num_str = f"{'' if p == 1 else p}{rad_html}"
            if q == 1:
                return num_str, is_neg, False, False
            if stacked_frac:
                return _frac_html(num_str, str(q)), is_neg, False, True
            return f"<sup>{num_str}</sup>&frasl;<sub>{q}</sub>", is_neg, False, False

    s_abs = fmt_num(abs(val), 2)
    return s_abs, is_neg, (s_abs == "1"), False

def pulisci_formula(expr, decimale=False, stacked_frac=False):
    """Formatta espressioni SymPy in HTML matematico pulito (niente *, ^, _, √(k) o frazioni giganti)."""
    e = sp.expand(pulisci_espr(expr))
    if e == 0:
        return "0"

    s_sym = sp.Symbol("s")
    if decimale and not e.free_symbols:
        return fmt_num(float(sp.N(e)), 2).replace("-", "&minus;")

    if decimale and e.free_symbols <= {s_sym}:
        try:
            poly = sp.Poly(e, s_sym)
            coeffs = poly.all_coeffs()
            deg_max = len(coeffs) - 1
            parti = []
            for idx_c, c_ in enumerate(coeffs):
                deg = deg_max - idx_c
                val = float(sp.N(c_))
                if abs(val) < 1e-9:
                    continue
                s_abs = fmt_num(abs(val), 2)
                if s_abs == "0":
                    continue
                neg = val < 0
                if deg == 0:
                    term = s_abs
                elif deg == 1:
                    term = "s" if s_abs == "1" else f"{s_abs}s"
                else:
                    term = f"s<sup>{deg}</sup>" if s_abs == "1" else f"{s_abs}s<sup>{deg}</sup>"
                if not parti:
                    parti.append(f"&minus;{term}" if neg else term)
                else:
                    parti.append(f" &minus; {term}" if neg else f" + {term}")
            return "".join(parti) if parti else "0"
        except Exception:
            pass

    termini = list(sp.Add.make_args(e))
    has_reaz = any(str(sym) != "s" for sym in e.free_symbols)

    def _chiave_ord(t_):
        syms = [sym for sym in t_.free_symbols]
        if not syms:
            return (2, 0, "")
        if not has_reaz and s_sym in syms:
            return (0, -int(sp.degree(t_, s_sym)), "")
        nomi = sorted(str(sym) for sym in syms)
        m_num = 0
        for p_ in nomi[0].split("_")[1:]:
            if p_.isdigit():
                m_num = int(p_)
                break
        return (0, m_num, nomi[0])

    termini.sort(key=_chiave_ord)

    parti_out = []
    any_stacked = False
    for t_ in termini:
        f_list = sp.Mul.make_args(t_)
        num_factors = [f_ for f_ in f_list if not f_.free_symbols]
        sym_factors = [f_ for f_ in f_list if f_.free_symbols]
        c_num = sp.Mul(*num_factors) if num_factors else sp.Integer(1)
        sym_part = sp.Mul(*sym_factors) if sym_factors else sp.Integer(1)
        s_abs, neg, is_unit, is_st = _format_coeff_num(c_num, decimale=decimale, stacked_frac=stacked_frac)
        if s_abs == "0":
            continue
        if is_st:
            any_stacked = True

        var_str = ""
        if sym_part != 1:
            vp = []
            for f_ in sym_part.as_ordered_factors():
                if isinstance(f_, sp.Pow):
                    base_s = sym_reaz_html(str(f_.base)) if str(f_.base) != "s" else "s"
                    vp.append(f"{base_s}<sup>{f_.exp}</sup>")
                else:
                    vp.append(sym_reaz_html(str(f_)) if str(f_) != "s" else "s")
            var_str = "".join(vp)

        if not parti_out:
            segno = "&minus;" if neg else ""
        else:
            segno = " &minus; " if neg else " + "

        coeff_str = "" if (var_str and is_unit) else s_abs
        parti_out.append((segno, coeff_str, var_str, is_st))

    if not parti_out:
        return "0"

    if stacked_frac and any_stacked:
        celle = []
        for idx_p, (seg, c_str, v_str, is_st) in enumerate(parti_out):
            if seg.strip():
                celle.append(seg if idx_p > 0 else "&minus;")
            if is_st:
                celle.append(c_str)
                if v_str:
                    celle.append(v_str)
            else:
                sep = " " if (c_str and v_str and "&radic;" in c_str) else ""
                celle.append(f"{c_str}{sep}{v_str}")
        return celle

    out_str = []
    for seg, c_str, v_str, _ in parti_out:
        sep = " " if (c_str and v_str and ("&radic;" in c_str or "</sup>" in c_str)) else ""
        out_str.append(f"{seg}{c_str}{sep}{v_str}")
    return "".join(out_str)

def formatta_equazione_html(lhs_expr, rhs_str="0"):
    """Formatta un'equazione in HTML pulito, usando una tabella orizzontale se contiene frazioni impilate."""
    res = pulisci_formula(lhs_expr, decimale=False, stacked_frac=True)
    if isinstance(res, list):
        return _eq_table_html(*res, f"= {rhs_str}")
    return f"<span class='formula'>{res} = {rhs_str}</span>"

def plot_sistema_S0(nodi, elementi, reazioni, valori_calcolati, filename, out_dir, tipo_grafico="versi", con_unita=False, cerniere_plot=None):
    plt, patches = _get_mpl()
    import re as _re
    nodi = {k: (float(v[0]), float(v[1])) for k, v in nodi.items()}   # solo per il disegno
    fig, ax = plt.subplots(figsize=(10.2, 3.9), dpi=140)

    xs_n = [c[0] for c in nodi.values()]
    ys_n = [c[1] for c in nodi.values()]
    span_x = max(xs_n) - min(xs_n)
    span_y = max(ys_n) - min(ys_n)
    span_g = max(span_x, span_y, 1.0)
    bar_lens = [math.hypot(nodi[int(e['n2'])][0] - nodi[int(e['n1'])][0],
                           nodi[int(e['n2'])][1] - nodi[int(e['n1'])][1]) for e in elementi]
    L_arr_ext = min(0.145 * span_g, max(0.095 * span_g, max(bar_lens + [0.5]) * 0.16))
    off_s = 0.075 * span_g

    bar_segs = [(nodi[int(e['n1'])], nodi[int(e['n2'])]) for e in elementi]
    node_pts = list(nodi.values())
    placed_boxes = []   # (cx, cy, hw, hh)
    drawn_arrows = []   # ((x1, y1), (x2, y2))
    all_x = list(xs_n)
    all_y = list(ys_n)

    def _dist_pt_seg(px, py, ax_, ay_, bx_, by_):
        dx_, dy_ = bx_ - ax_, by_ - ay_
        l2 = dx_ * dx_ + dy_ * dy_
        if l2 < 1e-12:
            return math.hypot(px - ax_, py - ay_)
        t = max(0.0, min(1.0, ((px - ax_) * dx_ + (py - ay_) * dy_) / l2))
        return math.hypot(px - (ax_ + t * dx_), py - (ay_ + t * dy_))

    def _est_box(txt):
        clean = _re.sub(r'[\$_{}\\]', '', txt)
        n_ch = max(2, len(clean))
        hw = (0.0205 * n_ch + 0.028) * span_g
        hh = 0.048 * span_g
        return hw, hh

    def _score_box(cx, cy, hw, hh, link_x, link_y):
        pen = 0.0
        # 1. Collisione con altri riquadri già piazzati (con piccolo margine di sicurezza)
        pad_b = 0.014 * span_g
        for lx, ly, lhw, lhh in placed_boxes:
            ox = max(0.0, (hw + lhw + pad_b) - abs(cx - lx))
            oy = max(0.0, (hh + lhh + pad_b) - abs(cy - ly))
            if ox > 0.0 and oy > 0.0:
                pen += 8000.0 * (ox / span_g) * (oy / span_g) + 120.0
        # 2. Collisione con i nodi
        r_node = 0.072 * span_g
        for nx0, ny0 in node_pts:
            ox = max(0.0, (hw + r_node) - abs(cx - nx0))
            oy = max(0.0, (hh + r_node) - abs(cy - ny0))
            if ox > 0.0 and oy > 0.0:
                pen += 4500.0 * (ox / span_g) * (oy / span_g) + 80.0
        # 3. Collisione con aste strutturali (griglia 5x3 su tutto l'ingombro del riquadro)
        sample_pts = [
            (cx + sx * hw * 0.98, cy + sy * hh * 0.92)
            for sx in (-1.0, -0.5, 0.0, 0.5, 1.0)
            for sy in (-1.0, 0.0, 1.0)
        ]
        bar_clear = max(hh * 0.65, 0.036 * span_g)
        for (x1, y1), (x2, y2) in bar_segs:
            d_min = min(_dist_pt_seg(px_, py_, x1, y1, x2, y2) for px_, py_ in sample_pts)
            if d_min < bar_clear:
                pen += 2200.0 * ((bar_clear - d_min) / span_g) + 60.0
        # 4. Collisione con frecce già disegnate
        arr_clear = max(hh * 0.50, 0.028 * span_g)
        for (ax1, ay1), (ax2, ay2) in drawn_arrows:
            d_min = min(_dist_pt_seg(px_, py_, ax1, ay1, ax2, ay2) for px_, py_ in sample_pts)
            if d_min < arr_clear:
                pen += 2000.0 * ((arr_clear - d_min) / span_g) + 55.0
        # 5. Evita se possibile che il richiamo tratteggiato tagli i nodi
        for nx0, ny0 in node_pts:
            if math.hypot(nx0 - link_x, ny0 - link_y) > 0.05 * span_g:
                dn_ = _dist_pt_seg(nx0, ny0, link_x, link_y, cx, cy)
                if dn_ < 0.045 * span_g:
                    pen += 120.0
        # Leggera preferenza per rimanere vicino al proprio punto di aggancio
        pen += 8.0 * math.hypot(cx - link_x, cy - link_y) / span_g
        return pen

    def _place_label(pref_x, pref_y, link_x, link_y, label_txt, col, box_style, fontsize=8.8):
        hw, hh = _est_box(label_txt)
        ux_p, uy_p = pref_x - link_x, pref_y - link_y
        lp = math.hypot(ux_p, uy_p)
        if lp > 1e-9:
            ux_p, uy_p = ux_p / lp, uy_p / lp
        else:
            ux_p, uy_p = 0.0, 1.0
        nx_p, ny_p = -uy_p, ux_p

        cands = [(pref_x, pref_y)]
        # Offset direzionali intorno all'uscita della freccia
        for dist_k in (0.04, 0.09, 0.14, 0.20, 0.27):
            rk = dist_k * span_g
            cands.append((pref_x + rk * ux_p, pref_y + rk * uy_p))
            for s_side in (1.0, -1.0):
                cands.append((pref_x + rk * nx_p * s_side, pref_y + rk * ny_p * s_side))
                cands.append((pref_x + rk * 0.65 * (ux_p + nx_p * s_side), pref_y + rk * 0.65 * (uy_p + ny_p * s_side)))
        # Anello completo attorno a link per casi molto affollati
        for r_k in (0.14, 0.21, 0.29, 0.37):
            rk = r_k * span_g
            for deg in range(0, 360, 30):
                rad = math.radians(deg)
                cands.append((link_x + rk * math.cos(rad), link_y + rk * math.sin(rad)))

        best_x, best_y = cands[0]
        best_s = _score_box(best_x, best_y, hw, hh, link_x, link_y)
        for cx_, cy_ in cands[1:]:
            sc_ = _score_box(cx_, cy_, hw, hh, link_x, link_y)
            if sc_ < best_s:
                best_s = sc_
                best_x, best_y = cx_, cy_

        placed_boxes.append((best_x, best_y, hw + 0.008 * span_g, hh + 0.008 * span_g))
        all_x.extend([best_x - hw, best_x + hw, link_x])
        all_y.extend([best_y - hh, best_y + hh, link_y])

        # Piccolo segmento di richiamo tratteggiato solo se il box e' distante dal punto di attacco
        dx_r = max(0.0, abs(link_x - best_x) - hw * 0.9)
        dy_r = max(0.0, abs(link_y - best_y) - hh * 0.9)
        if math.hypot(dx_r, dy_r) > 0.025 * span_g:
            ax.plot([link_x, best_x], [link_y, best_y], color=col, lw=0.9, ls=':', alpha=0.75, zorder=5)

        ax.text(best_x, best_y, label_txt, color=col, weight='bold', fontsize=fontsize,
                bbox=box_style, ha='center', va='center', zorder=7)
        return best_x, best_y

    def _draw_fixed_arrow(tip_x, tip_y, dx_u, dy_u, arr_len, col, lw_shaft=2.2, hw_arr=8.5):
        tail_x = tip_x - dx_u * arr_len
        tail_y = tip_y - dy_u * arr_len
        ax.annotate('', xy=(tip_x, tip_y), xytext=(tail_x, tail_y),
                    arrowprops=dict(facecolor=col, edgecolor=col, width=lw_shaft, headwidth=hw_arr, shrink=0),
                    zorder=6)
        drawn_arrows.append(((tail_x, tail_y), (tip_x, tip_y)))
        all_x.extend([tip_x, tail_x])
        all_y.extend([tip_y, tail_y])
        return tail_x, tail_y

    # 1. Disegna le travi e gli assi monoassiali (ascissa s)
    for el in elementi:
        n1, n2 = nodi[int(el['n1'])], nodi[int(el['n2'])]
        ax.plot([n1[0], n2[0]], [n1[1], n2[1]], 'k-', linewidth=3.2, zorder=1)
        if tipo_grafico == "monoassiali":
            dx, dy = n2[0] - n1[0], n2[1] - n1[1]
            L = math.hypot(dx, dy)
            if L > 0:
                ux, uy = dx / L, dy / L
                nx, ny = -uy, ux
                ax.annotate('', xy=(n1[0] + ux * L * 0.44, n1[1] + uy * L * 0.44),
                            xytext=(n1[0] + ux * L * 0.08, n1[1] + uy * L * 0.08),
                            arrowprops=dict(edgecolor='purple', arrowstyle='->', lw=2.2), zorder=4)
                sx_c = n1[0] + ux * L * 0.26 - nx * off_s * 0.55
                sy_c = n1[1] + uy * L * 0.26 - ny * off_s * 0.55
                placed_boxes.append((sx_c, sy_c, 0.042 * span_g, 0.042 * span_g))
                ax.text(sx_c, sy_c, 's', color='purple', weight='bold', fontsize=10.5,
                        ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85), zorder=5)

    # 2. Disegna i nodi (cerchi + ID; eventuale anello arancione sulle cerniere interne)
    cern_nodi = {c['nodo'] for c in (cerniere_plot or [])}
    for nid, coords in nodi.items():
        if nid in cern_nodi:
            ax.scatter([coords[0]], [coords[1]], c='none', s=270, zorder=2, edgecolors='#fd7e14', linewidths=2.2)
        ax.scatter(coords[0], coords[1], c='white', s=150, zorder=3, edgecolors='black', linewidth=2.2)
        ax.text(coords[0], coords[1], f"{nid}", fontsize=10.5, weight='bold', ha='center', va='center', color='black', zorder=4)

    # 3. Reazioni vincolari esterne (direzione fisica invariante + riquadro anti-collisione)
    for r in reazioni:
        xn, yn = nodi[r['nodo']]
        sym_lbl = sym_reaz_mpl(str(r['sym']))
        if tipo_grafico == "versi":
            label = sym_lbl
            colore_freccia = 'darkorange'
        else:
            valore = valori_calcolati.get(r['sym'], 0)
            try:
                num_val = float(valore)
                valore_str = fmt_num(num_val, 2, plain=True)
                colore_freccia = 'red' if num_val < 0 else 'forestgreen'
            except Exception:
                valore_str = str(valore)
                colore_freccia = 'forestgreen'
            unita_str = (" N·m" if r['tipo'] == 'M' else " N") if con_unita else ""
            label = f"{sym_lbl} = {valore_str}{unita_str}"

        box_style = dict(boxstyle='round,pad=0.24', fc='white', ec=colore_freccia, alpha=0.95, lw=1.4)

        if r['tipo'] == 'Fx':
            tx, ty = _draw_fixed_arrow(xn, yn, 1.0, 0.0, L_arr_ext, colore_freccia)
            _place_label(tx - 0.06 * span_g, ty, tx, ty, label, colore_freccia, box_style, fontsize=9.2)
        elif r['tipo'] == 'Fy':
            tx, ty = _draw_fixed_arrow(xn, yn, 0.0, 1.0, L_arr_ext, colore_freccia)
            _place_label(tx, ty - 0.06 * span_g, tx, ty, label, colore_freccia, box_style, fontsize=9.2)
        elif r['tipo'] == 'M':
            r_arc = 0.062 * span_g
            arc = patches.Arc((xn, yn), width=2 * r_arc, height=2 * r_arc, angle=0,
                              theta1=30, theta2=295, lw=2.3, color=colore_freccia, zorder=5)
            ax.add_patch(arc)
            th_end = math.radians(295)
            ex_a, ey_a = xn + r_arc * math.cos(th_end), yn + r_arc * math.sin(th_end)
            ax.annotate('', xy=(ex_a + 0.018 * span_g, ey_a + 0.008 * span_g), xytext=(ex_a - 0.01 * span_g, ey_a - 0.005 * span_g),
                        arrowprops=dict(facecolor=colore_freccia, edgecolor=colore_freccia, width=1.8, headwidth=7, shrink=0), zorder=6)
            lx_a, ly_a = xn + r_arc * 0.707, yn + r_arc * 0.707
            _place_label(xn + 0.14 * span_g, yn + 0.12 * span_g, lx_a, ly_a, label, colore_freccia, box_style, fontsize=9.2)
        elif r['tipo'] == 'Fn':
            _ang = math.radians(float(r.get('angle', 0)) + 90)
            _dx, _dy = math.cos(_ang), math.sin(_ang)
            if r.get('pendolo'):
                xa, ya = xn + L_arr_ext * 1.1 * _dx, yn + L_arr_ext * 1.1 * _dy
                ax.plot([xn, xa], [yn, ya], color='dimgray', lw=3, zorder=1, solid_capstyle='round')
                ax.plot([xa], [ya], marker='o', markersize=9, markerfacecolor='white', markeredgecolor='dimgray', zorder=3)
                nx_, ny_ = -_dy, _dx
                for k_ in (-2, -1, 0, 1, 2):
                    bx_, by_ = xa + nx_ * 0.10 * L_arr_ext * k_, ya + ny_ * 0.10 * L_arr_ext * k_
                    ax.plot([bx_, bx_ + _dx * 0.16 * L_arr_ext + nx_ * 0.05 * L_arr_ext],
                            [by_, by_ + _dy * 0.16 * L_arr_ext + ny_ * 0.05 * L_arr_ext], color='dimgray', lw=1, zorder=1)
                node_pts.append((xa, ya))
            tx, ty = _draw_fixed_arrow(xn, yn, _dx, _dy, L_arr_ext, colore_freccia)
            _place_label(tx - _dx * 0.06 * span_g, ty - _dy * 0.06 * span_g, tx, ty, label, colore_freccia, box_style, fontsize=9.2)

    # 4. Reazioni interne alle cerniere (legenda bicolore: Blu = tratto A, Viola = tratto B opposto)
    for ch in (cerniere_plot or []):
        xc, yc = nodi[ch['nodo']]
        ux_d, uy_d = ch['dir']
        col_c = ch['colore']
        sgn_c = ch['sign']  # +1 sul tratto precedente/ramo (colore A), -1 sul tratto successivo/carrier (colore B)
        nx_c, ny_c = -uy_d, ux_d

        # Trova la lunghezza effettiva dell'asta del tronco uscente dalla cerniera
        len_tr = span_g * 0.5
        for el_ in elementi:
            a1, a2 = int(el_['n1']), int(el_['n2'])
            if ch['nodo'] not in (a1, a2):
                continue
            alt = a2 if a1 == ch['nodo'] else a1
            vx_, vy_ = nodi[alt][0] - xc, nodi[alt][1] - yc
            lv_ = math.hypot(vx_, vy_)
            if lv_ > 1e-9 and (vx_ * ux_d + vy_ * uy_d) / lv_ > 0.85:
                len_tr = lv_
                break

        anc_dist = min(0.145 * span_g, max(0.08 * span_g, 0.40 * len_tr))
        arr_c = min(0.115 * span_g, max(0.075 * span_g, 0.28 * len_tr))
        anc_x = xc + ux_d * anc_dist
        anc_y = yc + uy_d * anc_dist
        box_c = dict(boxstyle='round,pad=0.22', fc='white', ec=col_c, alpha=0.95, lw=1.3)

        for sym_k, is_horiz in ((ch['sym_H'], True), (ch['sym_V'], False)):
            sym_m = sym_reaz_mpl(str(sym_k))
            if tipo_grafico == "versi":
                dir_x = float(sgn_c) if is_horiz else 0.0
                dir_y = 0.0 if is_horiz else float(sgn_c)
                lbl_c = sym_m if sgn_c > 0 else ("-" + sym_m)
                fs_c = 8.6
            else:
                val_raw = float(valori_calcolati.get(sym_k, 0))
                val_eff = sgn_c * val_raw
                unita_c = " N" if con_unita else ""
                if abs(val_eff) < 1e-9:
                    continue
                dir_phys = 1.0 if val_eff >= 0 else -1.0
                dir_x = dir_phys if is_horiz else 0.0
                dir_y = 0.0 if is_horiz else dir_phys
                lbl_c = f"{sym_m}={fmt_num(abs(val_eff), 2, plain=True)}{unita_c}"
                fs_c = 8.4

            # Mantieni sempre la direzione esatta (dir_x, dir_y) senza attraversare il nodo cerniera
            par_dot = abs(dir_x * ux_d + dir_y * uy_d)
            if par_dot >= 0.707:
                # Forza quasi collineare all'asta: sfalsa leggermente dal tratto nero e resta lato tronco
                perp_side = 1.0 if sgn_c > 0 else -1.0
                cx0 = anc_x + nx_c * (0.032 * span_g * perp_side)
                cy0 = anc_y + ny_c * (0.032 * span_g * perp_side)
                tx_test = cx0 - 0.5 * arr_c * dir_x
                ty_test = cy0 - 0.5 * arr_c * dir_y
                ptx_test = cx0 + 0.5 * arr_c * dir_x
                pty_test = cy0 + 0.5 * arr_c * dir_y
                min_proj = min((tx_test - xc) * ux_d + (ty_test - yc) * uy_d,
                               (ptx_test - xc) * ux_d + (pty_test - yc) * uy_d)
                min_clear = 0.075 * span_g
                if min_proj < min_clear:
                    shift_u = min_clear - min_proj
                    cx0 += ux_d * shift_u
                    cy0 += uy_d * shift_u
                tip_x, tip_y = cx0 + 0.5 * arr_c * dir_x, cy0 + 0.5 * arr_c * dir_y
                pref_x = cx0 + nx_c * (0.075 * span_g * perp_side)
                pref_y = cy0 + ny_c * (0.075 * span_g * perp_side)
            else:
                # Forza perpendicolare all'asta: punta esattamente sull'ancora dell'asta
                tip_x, tip_y = anc_x, anc_y
                pref_x = tip_x - dir_x * (arr_c + 0.05 * span_g)
                pref_y = tip_y - dir_y * (arr_c + 0.05 * span_g)

            tail_x, tail_y = _draw_fixed_arrow(tip_x, tip_y, dir_x, dir_y, arr_c, col_c, lw_shaft=1.9, hw_arr=7.2)
            link_cx, link_cy = 0.5 * (tail_x + tip_x) if par_dot >= 0.707 else tail_x, 0.5 * (tail_y + tip_y) if par_dot >= 0.707 else tail_y
            _place_label(pref_x, pref_y, link_cx, link_cy, lbl_c, col_c, box_c, fontsize=fs_c)

    pad_x = 0.08 * span_g
    pad_y = 0.08 * span_g
    ax.set_xlim(min(all_x) - pad_x, max(all_x) + pad_x)
    ax.set_ylim(min(all_y) - pad_y, max(all_y) + pad_y)
    ax.margins(x=0.12, y=0.22)
    ax.set_aspect('equal')
    ax.axis('off')
    return _salva_figura_png(fig, out_dir, filename, pad_inches=0.16, dpi=140)

def plot_diagrammi_singola_asta(L, expr_M, expr_T, expr_N, nome_asta, filename, out_dir, nome_sistema="S₀"):
    """Genera i 3 grafici (M, T, N) affiancati per una singola asta (Versione Ant-Crash Termico)."""
    plt, _ = _get_mpl()
    fig, axes = plt.subplots(1, 3, figsize=(10, 2.5), dpi=120)
    fig.suptitle(f"Diagrammi Sistema {nome_sistema} - Asta {nome_asta}", fontsize=11, weight='bold')
    ped = "₁" if ("1" in str(nome_sistema) or "₁" in str(nome_sistema)) else "₀"
    unita_m = "" if ped == "₁" else " [N·m]"
    unita_f = "" if ped == "₁" else " [N]"
    
    s_sym = sp.Symbol("s")
    s_vals = np.linspace(0, float(L), 100)
    
    def eval_sym(expr):
        try:
            expr = sp.sympify(expr)
            
            # 1. FILTRO TERMICO: Azzera la polvere numerica (valori < 10^-10)
            expr = expr.xreplace({n: 0.0 for n in expr.atoms(sp.Float) if abs(n) < 1e-10})
            
            # 2. Pone a 0 eventuali simboli residui non calcolati
            simboli_extra = [sym for sym in expr.free_symbols if str(sym) != 's']
            for sym in simboli_extra:
                expr = expr.subs(sym, 0.0)
                
            if len(expr.free_symbols) == 0:
                return np.full_like(s_vals, float(expr.evalf()))
                
            # 3. Trasforma in funzione Numpy
            func = sp.lambdify(s_sym, expr, 'numpy')
            y = func(s_vals)
            
            if np.isscalar(y):
                return np.full_like(s_vals, float(y))
                
            # 4. Forza la rimozione di eventuali residui complessi/immaginari nascosti
            return np.real(np.asarray(y, dtype=complex)).astype(float)
            
        except Exception:
            # 5. Fallback estremo punto per punto
            try:
                return np.array([float(sp.re(expr.subs(s_sym, v).evalf())) for v in s_vals], dtype=float)
            except:
                return np.zeros_like(s_vals)

    M_vals = eval_sym(expr_M)
    T_vals = eval_sym(expr_T)
    N_vals = eval_sym(expr_N)
    
    # Momento (rosso): M positivo = fibre tese a sinistra dell'asta (n1->n2).
    axes[0].plot(s_vals, M_vals, color='crimson', lw=2)
    axes[0].fill_between(s_vals, 0, M_vals, color='crimson', alpha=0.2)
    axes[0].axhline(0, color='black', lw=1)
    axes[0].set_title(f"Momento M{ped}(s){unita_m}", fontsize=10)
    
    # Taglio (verde)
    axes[1].plot(s_vals, T_vals, color='forestgreen', lw=2)
    axes[1].fill_between(s_vals, 0, T_vals, color='forestgreen', alpha=0.2)
    axes[1].axhline(0, color='black', lw=1)
    axes[1].set_title(f"Taglio T{ped}(s){unita_f}", fontsize=10)
    
    # Sforzo Normale (blu)
    axes[2].plot(s_vals, N_vals, color='dodgerblue', lw=2)
    axes[2].fill_between(s_vals, 0, N_vals, color='dodgerblue', alpha=0.2)
    axes[2].axhline(0, color='black', lw=1)
    axes[2].set_title(f"Sforzo Normale N{ped}(s){unita_f}", fontsize=10)
    
    for ax in axes:
        ax.set_xlabel("ascissa s [m]", fontsize=8)
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.margins(x=0.06, y=0.18)

    try:
        plt.tight_layout()
        return _salva_figura_png(fig, out_dir, filename, pad_inches=0.15)
    except Exception:
        plt.close(fig)
        raise

def _dist_punto_segmento(px, py, x1, y1, x2, y2):
    """Distanza minima del punto (px, py) dal segmento (x1,y1)-(x2,y2)."""
    dx, dy = x2 - x1, y2 - y1
    l2 = dx * dx + dy * dy
    if l2 <= 1e-18:
        return math.hypot(px - x1, py - y1)
    t = max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / l2))
    return math.hypot(px - (x1 + t * dx), py - (y1 + t * dy))

def _seg_interseca(ax_, ay_, bx_, by_, cx_, cy_, dx_, dy_, eps=1e-9):
    """True se il segmento AB interseca propriamente il segmento CD."""
    def _ccw(p1x, p1y, p2x, p2y, p3x, p3y):
        return (p3y - p1y) * (p2x - p1x) - (p2y - p1y) * (p3x - p1x)
    d1 = _ccw(ax_, ay_, bx_, by_, cx_, cy_)
    d2 = _ccw(ax_, ay_, bx_, by_, dx_, dy_)
    d3 = _ccw(cx_, cy_, dx_, dy_, ax_, ay_)
    d4 = _ccw(cx_, cy_, dx_, dy_, bx_, by_)
    if ((d1 > eps and d2 < -eps) or (d1 < -eps and d2 > eps)) and ((d3 > eps and d4 < -eps) or (d3 < -eps and d4 > eps)):
        return True
    return False

def _disegna_diagramma_su_ax(ax, nodi_f, elementi, expr_dict, titolo, col, show_labels=True):
    """Disegna un diagramma globale (M, T o N) su un Axes Matplotlib con scala unica e etichette anti-sovrapposizione."""
    from matplotlib.path import Path as MplPath
    ax.set_title(titolo, fontsize=12, weight='bold', pad=12)
    x_sym = sp.Symbol("s")
    xs_all = [c[0] for c in nodi_f.values()]
    ys_all = [c[1] for c in nodi_f.values()]
    max_span = max(max(xs_all) - min(xs_all), max(ys_all) - min(ys_all), 1.0)
    lunghezze = [math.hypot(nodi_f[int(e['n2'])][0] - nodi_f[int(e['n1'])][0],
                            nodi_f[int(e['n2'])][1] - nodi_f[int(e['n1'])][1]) for e in elementi]
    min_L = min([l for l in lunghezze if l > 1e-9], default=max_span)

    # 1. Pre-valutazione per scala uniforme su tutte le aste
    eval_bars = []
    global_max = 1e-5
    for el in elementi:
        eid = el['id']
        n1, n2 = int(el['n1']), int(el['n2'])
        x1, y1 = nodi_f[n1]
        x2, y2 = nodi_f[n2]
        L = math.hypot(x2 - x1, y2 - y1)
        if L <= 1e-9:
            continue
        x_vals = np.linspace(0.0, L, 121)
        expr = sp.sympify(expr_dict[eid])
        expr = expr.xreplace({n: 0.0 for n in expr.atoms(sp.Float) if abs(n) < 1e-10})
        try:
            func = sp.lambdify(x_sym, expr, 'numpy')
            y_vals = func(x_vals)
            if np.isscalar(y_vals):
                y_vals = np.full(len(x_vals), float(y_vals))
            y_vals = np.real(np.asarray(y_vals, dtype=complex)).astype(float)
        except Exception:
            y_vals = np.array([float(sp.re(expr.subs(x_sym, v).evalf())) for v in x_vals])
        global_max = max(global_max, float(np.max(np.abs(y_vals))))
        eval_bars.append((el, L, x1, y1, x2, y2, x_vals, y_vals, expr))

    # Scala globale compatta per evitare sovrapposizioni nei nodi d'angolo tra aste corte e lunghe
    amp = min(max_span * 0.105, max(min_L * 0.34, max_span * 0.052))
    scale_f = amp / global_max

    # 2. Scheletro struttura e numeri nodi
    segmenti = []
    for el in elementi:
        n1, n2 = int(el['n1']), int(el['n2'])
        x1, y1 = nodi_f[n1]
        x2, y2 = nodi_f[n2]
        ax.plot([x1, x2], [y1, y2], color='#868e96', linestyle='--', linewidth=1.6, zorder=1)
        segmenti.append((x1, y1, x2, y2))

    for nid, (nx, ny) in nodi_f.items():
        ax.scatter([nx], [ny], c='white', s=70, edgecolors='#495057', linewidths=1.2, zorder=4)
        ax.text(nx, ny, str(nid), fontsize=7, weight='bold', color='#495057', ha='center', va='center', zorder=5)

    # 3. Pass 1: disegno di tutti i diagrammi e raccolta geometrie complete + etichette candidate
    soglia_lab = max(1e-3, 1e-3 * global_max)
    etichette_fatte = set()
    all_px = list(xs_all)
    all_py = list(ys_all)
    diag_segs = []      # [(eid, x1, y1, x2, y2)] curve del diagramma
    fill_paths = []     # [MplPath] poligoni riempiti (non nulli)
    pending_labels = []

    for el, L, x1, y1, x2, y2, x_vals, y_vals, expr in eval_bars:
        eid = el['id']
        dx, dy = x2 - x1, y2 - y1
        ux, uy = -dy / L, dx / L
        tx, ty = dx / L, dy / L
        bx = x1 + tx * x_vals
        by = y1 + ty * x_vals
        px = bx + ux * y_vals * scale_f
        py = by + uy * y_vals * scale_f

        all_px.extend(px.tolist())
        all_py.extend(py.tolist())

        ax.fill(np.append(bx, px[::-1]), np.append(by, py[::-1]), col, alpha=0.22, zorder=2)
        ax.plot(px, py, col, linewidth=2.0, zorder=3)

        if float(np.max(np.abs(y_vals))) > soglia_lab:
            step_p = max(1, len(x_vals) // 16)
            idxs = list(range(0, len(x_vals), step_p))
            if idxs[-1] != len(x_vals) - 1:
                idxs.append(len(x_vals) - 1)
            poly_xy = np.vstack([
                np.column_stack([bx[idxs], by[idxs]]),
                np.column_stack([px[idxs[::-1]], py[idxs[::-1]]]),
                [[bx[0], by[0]]]
            ])
            fill_paths.append(MplPath(poly_xy))
            for k_s in range(len(idxs) - 1):
                i0, i1 = idxs[k_s], idxs[k_s + 1]
                diag_segs.append((eid, float(px[i0]), float(py[i0]), float(px[i1]), float(py[i1])))

        if not show_labels:
            continue

        coef_k = _coeff_poly(expr, x_sym)
        chiavi = [(0.0, 'estremo'), (L, 'estremo')]
        if coef_k is not None:
            c_k = np.trim_zeros(np.array(coef_k, dtype=float), 'f')
            if len(c_k) > 1:
                for r_ in _radici_in(list(np.polyder(c_k)), L):
                    if 1e-9 * L < r_ < L * (1 - 1e-9):
                        chiavi.append((float(r_), 'estremo interno'))
                for r_ in _radici_in(list(c_k), L):
                    if 1e-9 * L < r_ < L * (1 - 1e-9):
                        chiavi.append((float(r_), 'zero'))
            val_k = lambda t_, cc=c_k: float(np.polyval(cc, t_)) if len(cc) else 0.0
        else:
            chiavi.append((float(x_vals[int(np.argmax(np.abs(y_vals)))]), 'estremo interno'))
            val_k = lambda t_, xv=x_vals, yv=y_vals: float(np.interp(t_, xv, yv))

        for t_k, tipo_k in chiavi:
            v_k = val_k(t_k)
            bxk = x1 + tx * t_k
            byk = y1 + ty * t_k
            pxk, pyk = bxk + ux * v_k * scale_f, byk + uy * v_k * scale_f
            if tipo_k == 'zero':
                ax.plot([bxk], [byk], marker='o', markersize=4.5, markerfacecolor='white', markeredgecolor=col, zorder=5)
                continue
            if abs(v_k) <= soglia_lab:
                continue
            chiave_lab = (round(pxk, 3), round(pyk, 3), round(v_k, 1))
            if chiave_lab in etichette_fatte:
                continue
            etichette_fatte.add(chiave_lab)
            ax.plot([pxk], [pyk], marker='o', markersize=3.5, color=col, zorder=5)
            pending_labels.append((eid, pxk, pyk, tx, ty, ux, uy, v_k, L, t_k))

    # 4. Pass 2: posizionamento globale collision-free (con conoscenza di tutte le curve e fill)
    placed_boxes = []   # [(cx, cy, hw, hh)]
    placed_leaders = [] # [(pxk, pyk, cx, cy)]
    bar_pad = max_span * 0.024
    r_cands = (max_span * 0.080, max_span * 0.125, max_span * 0.175)
    angoli_deg = (0, 35, -35, 65, -65, 95, -95, 130, -130, 160, -160, 180)

    for eid, pxk, pyk, tx, ty, ux, uy, v_k, L, t_k in pending_labels:
        lbl_txt = fmt_num(v_k, 2, plain=True)
        hw = max_span * (0.018 + 0.0115 * len(lbl_txt))
        hh = max_span * 0.026
        sg_ = 1.0 if v_k >= 0 else -1.0
        nx_out, ny_out = ux * sg_, uy * sg_
        base_ang = math.atan2(ny_out, nx_out)

        # Direzione tangente verso l'esterno del nodo d'estremita' (aiuta nei nodi d'angolo)
        if t_k < 0.25 * L:
            out_tx, out_ty = -tx, -ty
        elif t_k > 0.75 * L:
            out_tx, out_ty = tx, ty
        else:
            out_tx, out_ty = 0.0, 0.0

        best_score = float('inf')
        best_c = (pxk + nx_out * r_cands[0], pyk + ny_out * r_cands[0])

        for r_idx, rad in enumerate(r_cands):
            for a_idx, da_deg in enumerate(angoli_deg):
                ang = base_ang + math.radians(da_deg)
                dx_c, dy_c = math.cos(ang), math.sin(ang)
                cx = pxk + dx_c * rad
                cy = pyk + dy_c * rad

                # Preferenza: normale esterna o allontanamento verso l'esterno sui nodi terminali
                score = a_idx * 0.08 + r_idx * 1.1
                if out_tx or out_ty:
                    dot_out = dx_c * out_tx + dy_c * out_ty
                    if dot_out < -0.15:
                        score += 1.5

                # 1) Collisione box-box con etichette gia' piazzate
                for lx, ly, lhw, lhh in placed_boxes:
                    ox_ = (hw + lhw) * 1.15 - abs(cx - lx)
                    oy_ = (hh + lhh) * 1.25 - abs(cy - ly)
                    if ox_ > 0 and oy_ > 0:
                        score += 4000.0 + (ox_ + oy_) * 2000.0

                # 2) Campioni del box contro fill polygons, aste e curve diagramma
                samp_pts = [
                    (cx, cy),
                    (cx - hw * 0.8, cy), (cx + hw * 0.8, cy),
                    (cx, cy - hh * 0.8), (cx, cy + hh * 0.8),
                    (cx - hw * 0.8, cy - hh * 0.8), (cx + hw * 0.8, cy + hh * 0.8),
                    (cx - hw * 0.8, cy + hh * 0.8), (cx + hw * 0.8, cy - hh * 0.8),
                ]
                for fpoly in fill_paths:
                    hits_in = fpoly.contains_points(samp_pts)
                    if np.any(hits_in):
                        score += 220.0 * float(np.sum(hits_in))

                for sx1, sy1, sx2, sy2 in segmenti:
                    d_min_bar = min(_dist_punto_segmento(px_s, py_s, sx1, sy1, sx2, sy2) for px_s, py_s in samp_pts)
                    if d_min_bar < bar_pad:
                        score += (bar_pad - d_min_bar) / max(bar_pad, 1e-6) * 260.0 + 180.0

                for deid, dx1, dy1, dx2, dy2 in diag_segs:
                    d_min_curv = min(_dist_punto_segmento(px_s, py_s, dx1, dy1, dx2, dy2) for px_s, py_s in samp_pts)
                    if d_min_curv < bar_pad:
                        score += (bar_pad - d_min_curv) / max(bar_pad, 1e-6) * 240.0 + 160.0

                for nx_c, ny_c in nodi_f.values():
                    d_nd = min(math.hypot(px_s - nx_c, py_s - ny_c) for px_s, py_s in samp_pts)
                    if d_nd < bar_pad * 1.3:
                        score += (bar_pad * 1.3 - d_nd) / max(bar_pad, 1e-6) * 240.0 + 160.0

                # 3) Incrocio della leader line con aste, curve d'altre aste o altre leader
                lp1x = pxk + (cx - pxk) * 0.08
                lp1y = pyk + (cy - pyk) * 0.08
                for sx1, sy1, sx2, sy2 in segmenti:
                    if _seg_interseca(lp1x, lp1y, cx, cy, sx1, sy1, sx2, sy2):
                        score += 420.0
                for deid, dx1, dy1, dx2, dy2 in diag_segs:
                    if _seg_interseca(lp1x, lp1y, cx, cy, dx1, dy1, dx2, dy2):
                        score += 380.0 if deid != eid else 250.0
                for lpx1, lpy1, lpx2, lpy2 in placed_leaders:
                    if _seg_interseca(lp1x, lp1y, cx, cy, lpx1, lpy1, lpx2, lpy2):
                        score += 650.0

                if score < best_score:
                    best_score = score
                    best_c = (cx, cy)
                    best_box = (cx, cy, hw, hh)

        placed_boxes.append(best_box)
        placed_leaders.append((pxk, pyk, best_c[0], best_c[1]))
        all_px.extend([best_c[0] - hw, best_c[0] + hw])
        all_py.extend([best_c[1] - hh, best_c[1] + hh])
        ax.annotate(
            lbl_txt,
            xy=(pxk, pyk),
            xytext=best_c,
            textcoords='data',
            arrowprops=dict(arrowstyle='-', color=col, lw=0.85, alpha=0.65),
            color=col,
            weight='bold',
            fontsize=8,
            ha='center',
            va='center',
            bbox=dict(boxstyle='round,pad=0.16', fc='white', ec=col, alpha=0.92),
            zorder=6,
        )

    pad_x = max_span * 0.06
    pad_y = max_span * 0.06
    ax.set_xlim(min(all_px) - pad_x, max(all_px) + pad_x)
    ax.set_ylim(min(all_py) - pad_y, max(all_py) + pad_y)
    ax.margins(x=0.12, y=0.22)
    ax.set_aspect('equal')
    ax.axis('off')

def plot_momenti_S0_S1(nodi, elementi, M0_dict, M1_dict, grado_iperstaticita, out_dir, filename="plv_momenti_S0_S1.png"):
    """Plotta i 2 diagrammi globali di momento flettente M0 e M1 sul telaio intero a supporto degli integrali di Mohr."""
    plt, _ = _get_mpl()
    nodi_f = {k: (float(v[0]), float(v[1])) for k, v in nodi.items()}
    if grado_iperstaticita > 0:
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), dpi=130)
        _disegna_diagramma_su_ax(axes[0], nodi_f, elementi, M0_dict, "Momento Flettente M₀(s) [N·m] (Sistema S₀)\n(lato fibre tese)", "crimson")
        _disegna_diagramma_su_ax(axes[1], nodi_f, elementi, M1_dict, "Momento Flettente M₁(s) (Sistema S₁, X₁ = 1)\n(lato fibre tese)", "darkorange")
    else:
        fig, ax = plt.subplots(1, 1, figsize=(7.5, 5.0), dpi=130)
        _disegna_diagramma_su_ax(ax, nodi_f, elementi, M0_dict, "Momento Flettente M₀(s) [N·m] (Sistema Isostatico)\n(lato fibre tese)", "crimson")
    try:
        plt.tight_layout()
        return _salva_figura_png(fig, out_dir, filename, pad_inches=0.16, dpi=130)
    except Exception:
        plt.close(fig)
        raise

def _scala_dim_png(img_path, max_w, max_h=None):
    """Calcola attributi width/height proporzionali leggendo l'header PNG per stare nella pagina."""
    try:
        p_abs = _risolvi_percorso_img(img_path)
        with open(p_abs, "rb") as fp:
            hdr = fp.read(24)
        if len(hdr) >= 24 and hdr[:8] == b'\x89PNG\r\n\x1a\n':
            w_px, h_px = struct.unpack('>II', hdr[16:24])
            if w_px > 0 and h_px > 0:
                sc = min(float(max_w) / w_px, float(max_h) / h_px) if max_h else (float(max_w) / w_px)
                return max(80, int(round(w_px * sc))), max(50, int(round(h_px * sc)))
    except Exception:
        pass
    return int(max_w), None

def _tabella_html(headers, rows, keep=True, compact=False):
    """Genera una tabella HTML con <thead> esplicito (ripetizione riga testata) e class anti-spezzamento."""
    cl = []
    if keep: cl.append("keep-table")
    if compact: cl.append("compact-table")
    cls_attr = f" class='{' '.join(cl)}'" if cl else ""
    ths = "".join(f"<th>{h}</th>" for h in headers)
    trs = "".join(rows)
    return f"<table{cls_attr}><thead><tr>{ths}</tr></thead><tbody>{trs}</tbody></table>"

def _img_block_html(title_html, img_path, max_w, max_h=None, page_break_before=False):
    """Blocco titolo+immagine centrato, dimensionato proporzionalmente e indivisibile."""
    pb = " page-break-before: always;" if page_break_before else ""
    w_o, h_o = _scala_dim_png(img_path, max_w, max_h)
    dim_attr = f"width='{w_o}' height='{h_o}'" if h_o else f"width='{w_o}'"
    uri = _to_file_uri(img_path)
    return (
        f"<div class='keep-together' align='center' style='text-align:center; page-break-inside: avoid; break-inside: avoid;{pb} margin: 4px 0 8px 0;'>"
        f"<div class='keep-next' align='center' style='text-align:center; page-break-after: avoid; break-after: avoid; margin-bottom: 4px; line-height: 1.2;'>{title_html}</div>"
        f"<div align='center' style='text-align:center; line-height: 1.0; font-size: 1px; margin: 0; padding: 0;'><img src='{uri}' {dim_attr}></div>"
        "</div>"
    )

def _salva_trio_diagrammi(nodi_f, elementi, M_dict, T_dict, N_dict, tit_M, tit_T, tit_N, filename, out_dir):
    """Genera e salva i 3 diagrammi globali impilati (M, T, N) su un'unica pagina compatta."""
    plt, _ = _get_mpl()
    fig, axes = plt.subplots(3, 1, figsize=(9.6, 8.6), dpi=130)
    try:
        _disegna_diagramma_su_ax(axes[0], nodi_f, elementi, M_dict, tit_M, "crimson")
        _disegna_diagramma_su_ax(axes[1], nodi_f, elementi, T_dict, tit_T, "forestgreen")
        _disegna_diagramma_su_ax(axes[2], nodi_f, elementi, N_dict, tit_N, "dodgerblue")
        fig.tight_layout(h_pad=1.0)
        return _salva_figura_png(fig, out_dir, filename, pad_inches=0.16, dpi=130)
    except Exception:
        plt.close(fig)
        raise

def plot_schema_generale(nodi, elementi, supports_n, internal_hinges, point_loads, element_loads, thermal_loads, out_dir, filename="plv_S0_S1.png"):
    """Disegna lo schema reale completo della struttura (aste, nodi, vincoli, cerniere, carichi)."""
    plt, patches = _get_mpl()
    nodi_f = {k: (float(v[0]), float(v[1])) for k, v in nodi.items()}
    xs = [c[0] for c in nodi_f.values()]
    ys = [c[1] for c in nodi_f.values()]
    span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
    u = 0.12 * span

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=140)

    def _verso_str(nid):
        sx = sy = 0.0
        x0, y0 = nodi_f[nid]
        for el in elementi:
            n1, n2 = int(el['n1']), int(el['n2'])
            if nid not in (n1, n2):
                continue
            altro = n2 if n1 == nid else n1
            dx, dy = nodi_f[altro][0] - x0, nodi_f[altro][1] - y0
            L = math.hypot(dx, dy)
            if L > 0:
                sx += dx / L; sy += dy / L
        return (sx, sy) if math.hypot(sx, sy) > 1e-9 else None

    def _dir_terra(nid, modo, ang_sup=0.0):
        m = _verso_str(nid)
        mx, my = (-m[0], -m[1]) if m else (0.0, -1.0)
        eps = 1e-9
        vert = (0.0, 1.0 if my > eps else -1.0)
        oriz = (1.0 if mx > eps else -1.0, 0.0)
        if modo == 'verticale': return vert
        if modo == 'orizzontale': return oriz
        if modo == 'normale':
            a = math.radians(ang_sup + 90.0)
            nx, ny = math.cos(a), math.sin(a)
            dot = nx * mx + ny * my
            if dot > eps: return (nx, ny)
            if dot < -eps: return (-nx, -ny)
            if ny < -eps: return (nx, ny)
            if ny > eps: return (-nx, -ny)
            return (nx, ny) if nx < 0 else (-nx, -ny)
        return oriz if abs(mx) > abs(my) + eps else vert

    def _sist(x, y, d):
        t = (-d[1], d[0])
        return lambda a, b: (x + d[0] * a + t[0] * b, y + d[1] * a + t[1] * b)

    def _terra(P, a, mezza, col, lw=2.2):
        ax.plot(*zip(P(a, -mezza), P(a, mezza)), color=col, lw=lw, zorder=3)
        for i in range(7):
            b_ = -mezza + 2 * mezza * i / 6
            ax.plot(*zip(P(a, b_), P(a + 0.22 * u, b_ - 0.22 * u)), color=col, lw=1.1, zorder=3)

    # 1. Aste e badge identificativo asta
    for el in elementi:
        n1, n2 = int(el['n1']), int(el['n2'])
        x1, y1 = nodi_f[n1]; x2, y2 = nodi_f[n2]
        ax.plot([x1, x2], [y1, y2], color='#212529', lw=3.2, zorder=2)
        xm, ym = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        L = math.hypot(x2 - x1, y2 - y1)
        if L > 0:
            nx, ny = -(y2 - y1) / L, (x2 - x1) / L
            ax.text(xm - nx * 0.22 * u, ym - ny * 0.22 * u, f"A{el['id']}", fontsize=8, weight='bold',
                    color='#495057', ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.15', fc='#f8f9fa', ec='#ced4da', alpha=0.9), zorder=3)

    # 2. Vincoli a terra
    for sup in supports_n:
        nid = int(sup['node_id'])
        if nid not in nodi_f: continue
        x, y = nodi_f[nid]
        col = '#dc3545'
        fu, fv, fw = bool(sup.get('fix_u')), bool(sup.get('fix_v')), bool(sup.get('fix_w'))
        if sup.get('pendolo') or sup.get('bipendolo'):
            ang_p = math.radians(float(sup.get('axis_angle', float(sup.get('angle', 0)) + 90)))
            ux, uy = math.cos(ang_p), math.sin(ang_p)
            nx, ny = -uy, ux
            Lp = max(math.hypot(float(sup['pend_dx']), float(sup['pend_dy'])) if sup.get('pend_dx') is not None and sup.get('pend_dy') is not None else 0.9, 0.35)
            offs = (-0.12 * u, 0.12 * u) if sup.get('bipendolo') else (0.0,)
            for o in offs:
                ax.plot([x + nx * o, x + nx * o + ux * Lp], [y + ny * o, y + ny * o + uy * Lp], color='#6f42c1', lw=2.8, zorder=3)
                ax.plot([x + nx * o, x + nx * o + ux * Lp], [y + ny * o, y + ny * o + uy * Lp], 'o', mfc='white', mec='#6f42c1', ms=6, zorder=4)
            _terra(_sist(x, y, (ux, uy)), Lp, 0.4 * u, '#343a40')
        elif sup.get('rotated_roller') or (fv and not fu and not fw) or (fu and not fv and not fw):
            d = _dir_terra(nid, 'normale', float(sup.get('angle', 0))) if sup.get('rotated_roller') else _dir_terra(nid, 'verticale' if fv else 'orizzontale')
            P = _sist(x, y, d)
            h_t, w_t, r_c = 0.50 * u, 0.38 * u, 0.10 * u
            ax.add_patch(patches.Polygon([P(0, 0), P(h_t, -w_t), P(h_t, w_t)], closed=True, facecolor='white', edgecolor=col, lw=1.8, zorder=3))
            for b_ in (-0.5 * w_t, 0.5 * w_t):
                ax.add_patch(patches.Circle(P(h_t + r_c, b_), r_c, facecolor='white', edgecolor=col, lw=1.5, zorder=3))
            _terra(P, h_t + 2 * r_c, 1.3 * w_t, col)
        elif fu and fv and fw:
            P = _sist(x, y, _dir_terra(nid, 'libero'))
            _terra(P, 0.0, 0.65 * u, col, lw=3.8)
        elif fu and fv:
            P = _sist(x, y, _dir_terra(nid, 'libero'))
            h_t, w_t = 0.55 * u, 0.42 * u
            ax.add_patch(patches.Polygon([P(0, 0), P(h_t, -w_t), P(h_t, w_t)], closed=True, facecolor='white', edgecolor=col, lw=1.8, zorder=3))
            _terra(P, h_t, 1.3 * w_t, col)
        elif (fv and fw) or (fu and fw):
            P = _sist(x, y, _dir_terra(nid, 'verticale' if fv else 'orizzontale'))
            Lr, off = 0.90 * u, 0.28 * u
            ax.plot(*zip(P(0, -off), P(0, off)), color=col, lw=2.8, zorder=3)
            for b_ in (-off, off):
                ax.plot(*zip(P(0, b_), P(Lr, b_)), color=col, lw=2.2, zorder=3)
                for a_ in (0.0, Lr):
                    ax.plot(*P(a_, b_), marker='o', ms=5, mfc='white', mec=col, zorder=4)
            _terra(P, Lr, 0.6 * u, col)

    # 3. Cerniere interne
    for cern in internal_hinges:
        nid = int(cern['node_id'])
        if nid in nodi_f:
            ax.plot(*nodi_f[nid], marker='o', markersize=20, mfc='none', mec='#fd7e14', mew=2.4, zorder=4)

    # 4. Carichi concentrati
    for load in point_loads:
        nid = int(load['node_id'])
        if nid not in nodi_f: continue
        x, y = nodi_f[nid]
        fx, fy, mom = float(load.get('Fx', 0)), float(load.get('Fy', 0)), float(load.get('M', 0))
        if abs(fx) > 1e-9:
            dx_a = -u if fx > 0 else u
            ax.annotate(f"{fmt_num(fx, 2)} N", xy=(x, y), xytext=(x + dx_a * 1.15, y),
                        arrowprops=dict(facecolor='magenta', edgecolor='magenta', width=2.2, headwidth=8),
                        color='magenta', fontsize=8.5, weight='bold', va='center',
                        ha='right' if dx_a < 0 else 'left', zorder=6)
        if abs(fy) > 1e-9:
            dy_a = -u if fy > 0 else u
            ax.annotate(f"{fmt_num(fy, 2)} N", xy=(x, y), xytext=(x, y + dy_a * 1.15),
                        arrowprops=dict(facecolor='magenta', edgecolor='magenta', width=2.2, headwidth=8),
                        color='magenta', fontsize=8.5, weight='bold', ha='center',
                        va='top' if dy_a < 0 else 'bottom', zorder=6)
        if abs(mom) > 1e-9:
            r = 0.50 * u
            t0, t1 = math.radians(40), math.radians(320)
            th = [t0 + (t1 - t0) * k / 32 for k in range(33)]
            if mom < 0: th = th[::-1]
            ax.plot([x + r * math.cos(t) for t in th[:-3]], [y + r * math.sin(t) for t in th[:-3]], color='magenta', lw=2.2, zorder=6)
            ax.annotate('', xy=(x + r * math.cos(th[-1]), y + r * math.sin(th[-1])),
                        xytext=(x + r * math.cos(th[-4]), y + r * math.sin(th[-4])),
                        arrowprops=dict(arrowstyle='-|>', color='magenta', lw=2.2, mutation_scale=14), zorder=6)
            ax.text(x + r * 1.15, y + r * 1.15, f"M={fmt_num(mom, 2)} N·m", color='magenta', fontsize=8.5, weight='bold',
                    bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.75), zorder=7)

    # 5. Carichi distribuiti
    for ql in element_loads:
        eid = ql.get('elem_id')
        asta = next((e for e in elementi if e['id'] == eid), None)
        if not asta: continue
        x1, y1 = nodi_f[int(asta['n1'])]; x2, y2 = nodi_f[int(asta['n2'])]
        dx_a, dy_a = x2 - x1, y2 - y1
        L = math.hypot(dx_a, dy_a)
        if L <= 1e-9: continue
        ux, uy = dx_a / L, dy_a / L
        q_t, q_xx, q_yy = float(ql.get('q_trans', 0)), float(ql.get('q_x', 0)), float(ql.get('q_y', 0))
        wx, wy = q_t * uy + q_xx, -q_t * ux + q_yy
        wmod = math.hypot(wx, wy)
        if wmod <= 1e-9: continue
        ox, oy = wx / wmod * u, wy / wmod * u
        tipo_q = ql.get('tipo', 'costante')
        if tipo_q == 'costante':
            px_q, py_q = [x1, x1 - ox, x2 - ox, x2, x1], [y1, y1 - oy, y2 - oy, y2, y1]
        elif tipo_q == 'triangolare_crescente':
            px_q, py_q = [x1, x2 - ox, x2, x1], [y1, y2 - oy, y2, y1]
        else:
            px_q, py_q = [x1, x1 - ox, x2, x1], [y1, y1 - oy, y2, y1]
        ax.fill(px_q, py_q, color='deepskyblue', alpha=0.16, zorder=4)
        ax.plot(px_q, py_q, color='deepskyblue', lw=1.6, zorder=5)
        for i_q in range(7):
            t_ = i_q / 6.0
            inten = 1.0 if tipo_q == 'costante' else (t_ if tipo_q == 'triangolare_crescente' else 1.0 - t_)
            if inten > 0.05:
                px_, py_ = x1 + t_ * dx_a, y1 + t_ * dy_a
                ax.annotate('', xy=(px_, py_), xytext=(px_ - ox * inten, py_ - oy * inten),
                            arrowprops=dict(facecolor='deepskyblue', edgecolor='deepskyblue', width=1.2, headwidth=5), zorder=5)
        ax.text(0.5 * (x1 + x2) - ox * 1.35, 0.5 * (y1 + y2) - oy * 1.35, f"q = {fmt_num(wmod, 2)} N/m",
                color='darkblue', fontsize=8.5, weight='bold', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.15', fc='lightcyan', ec='deepskyblue', alpha=0.9), zorder=6)

    # 5b. Variazioni termiche sulle aste
    for tl in (thermal_loads or []):
        eid = int(tl.get('elem_id', 0))
        asta = next((e for e in elementi if int(e['id']) == eid), None)
        if not asta:
            continue
        x1, y1 = nodi_f[int(asta['n1'])]; x2, y2 = nodi_f[int(asta['n2'])]
        dx_a, dy_a = x2 - x1, y2 - y1
        L = math.hypot(dx_a, dy_a)
        if L <= 1e-9:
            continue
        nx, ny = -dy_a / L, dx_a / L
        T_p, T_n = float(tl.get('T_pos', 0)), float(tl.get('T_neg', 0))
        xm, ym = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        ax.plot([x1, x2], [y1, y2], color='#fd7e14', lw=5.2, alpha=0.32, zorder=1)
        ax.text(xm + nx * 0.32 * u, ym + ny * 0.32 * u,
                f"ΔT: +{fmt_num(T_p, 1)}/{fmt_num(T_n, 1)} °C",
                fontsize=7.5, weight='bold', color='#b45309', ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.15', fc='#fff3cd', ec='#fd7e14', alpha=0.92), zorder=6)

    # 5c. Cedimenti / rotazioni imposte ai vincoli
    for sup in supports_n:
        nid = int(sup['node_id'])
        if nid not in nodi_f:
            continue
        x, y = nodi_f[nid]
        ced_parti = []
        for k_c, sim_c in (('delta_u', 'Δx'), ('delta_v', 'Δy'), ('delta_n', 'Δn'), ('delta_w', 'Δγ')):
            v_c = float(sup.get(k_c, 0) or 0)
            if abs(v_c) > 1e-12:
                unit_c = 'rad' if k_c == 'delta_w' else 'm'
                val_c = fmt_smart(v_c, fixed_min=1e-4, plain=True)
                ced_parti.append(f"{sim_c}={val_c} {unit_c}")
        if ced_parti:
            d_ced = _dir_terra(nid, 'libero')
            dx_c, dy_c = d_ced
            # Posiziona il badge sul lato libero dietro il vincolo/terra, evitando di coprire aste e nodo
            if abs(dx_c) >= abs(dy_c):
                bx_c = x + (0.42 if dx_c >= 0 else -0.42) * u
                by_c = y - 0.12 * u
                ha_c = 'left' if dx_c >= 0 else 'right'
                va_c = 'center'
            else:
                bx_c = x
                by_c = y + (1.05 if dy_c >= 0 else -1.05) * u
                ha_c = 'center'
                va_c = 'bottom' if dy_c >= 0 else 'top'
            ax.text(bx_c, by_c, ", ".join(ced_parti),
                    fontsize=7.5, weight='bold', color='#6f42c1', ha=ha_c, va=va_c,
                    bbox=dict(boxstyle='round,pad=0.16', fc='#f3e8ff', ec='#6f42c1', alpha=0.94), zorder=7)

    # 6. Nodi
    for nid, (x, y) in nodi_f.items():
        ax.plot(x, y, 'o', color='white', mec='#0d6efd', mew=2, ms=13, zorder=7)
        ax.text(x, y, str(nid), color='black', fontsize=8.5, weight='bold', ha='center', va='center', zorder=8)

    marg_x = 0.30 * span
    marg_y = 0.30 * span
    ax.set_xlim(min(xs) - marg_x, max(xs) + marg_x)
    ax.set_ylim(min(ys) - marg_y, max(ys) + marg_y)
    ax.margins(x=0.12, y=0.22)
    ax.set_aspect('equal')
    ax.axis('off')
    return _salva_figura_png(fig, out_dir, filename, pad_inches=0.16, dpi=140)

def plot_schema_libero_asta(L, n1_id, n2_id, expr_M, expr_T, expr_N, nome_asta, filename, out_dir, dati_json):
    """Genera lo schema di corpo libero completo con vettori (N,T,M) e carichi esterni."""
    plt, patches = _get_mpl()
    fig, ax = plt.subplots(figsize=(8, 2.8), dpi=120)
    
    # 1. Disegno l'asta orizzontale
    ax.plot([0, float(L)], [0, 0], 'k-', lw=4, zorder=2)
    
    # Nodi estremi
    ax.plot(0, 0, 'wo', markeredgecolor='black', markersize=14, mew=2, zorder=5)
    ax.plot(float(L), 0, 'wo', markeredgecolor='black', markersize=14, mew=2, zorder=5)
    ax.text(0, 0, str(n1_id), ha='center', va='center', fontsize=8, weight='bold', zorder=6)
    ax.text(float(L), 0, str(n2_id), ha='center', va='center', fontsize=8, weight='bold', zorder=6)
    
    s_sym = sp.Symbol("s")
    
    # Valutazione sicura
    def valuta_estremo(expr, val_s):
        try:
            expr = sp.sympify(expr)
            for sym in expr.free_symbols:
                if str(sym) != 's': expr = expr.subs(sym, 0.0)
            res = expr.subs(s_sym, val_s).evalf()
            return round(float(res), 2)
        except:
            return 0.0
            
    N0 = valuta_estremo(expr_N, 0); T0 = valuta_estremo(expr_T, 0); M0 = valuta_estremo(expr_M, 0)
    NL = valuta_estremo(expr_N, L); TL = valuta_estremo(expr_T, L); ML = valuta_estremo(expr_M, L)
    
    # 2. DISEGNO CARICHI DISTRIBUITI ESTERNI (q)
    q_val = 0
    for load in dati_json.get('element_loads', []):
        if str(load.get('elem_id')) == str(nome_asta):
            q_val = float(load.get('q_trans', 0))
            
    if q_val != 0:
        num_arrows = max(5, int(float(L) * 3))
        x_vals = np.linspace(0, float(L), num_arrows)
        
        # Convenzione: se q > 0 punta verso il basso (sul grafico)
        start_y = 0.5 if q_val > 0 else -0.5
        end_y = 0.05 if q_val > 0 else -0.05
        
        for x in x_vals:
            ax.annotate('', xy=(x, end_y), xytext=(x, start_y),
                        arrowprops=dict(edgecolor='royalblue', arrowstyle='->', lw=1.5))
        ax.plot([0, float(L)], [start_y, start_y], color='royalblue', lw=1.5)
        ax.text(float(L)/2, start_y + (0.15 if q_val > 0 else -0.2), f"q = {q_val}", 
                color='royalblue', weight='bold', ha='center', fontsize=9)
    q_ax_val = float(dati_json.get('q_ax', 0))
    if q_ax_val != 0:
        ax.text(float(L)/2, 0.8, f"q assiale = {q_ax_val:.4g}", color='royalblue', fontsize=8, ha='center')
                
    # 3. VETTORI SCHEMATICI SOLLECITAZIONI (Convenzione standard)
    # N - Sforzo Normale (blu)
    ax.annotate('', xy=(-0.4, 0), xytext=(0, 0), arrowprops=dict(edgecolor='dodgerblue', arrowstyle='<-', lw=2))
    ax.annotate('', xy=(float(L)+0.4, 0), xytext=(float(L), 0), arrowprops=dict(edgecolor='dodgerblue', arrowstyle='->', lw=2))
    
    # T - Taglio (verde)
    ax.annotate('', xy=(0, 0.4), xytext=(0, 0), arrowprops=dict(edgecolor='forestgreen', arrowstyle='<-', lw=2))
    ax.annotate('', xy=(float(L), -0.4), xytext=(float(L), 0), arrowprops=dict(edgecolor='forestgreen', arrowstyle='->', lw=2))
    
    # M - Momento Flettente (rosso scuro) - Archi stilizzati
    arc_start = patches.Arc((0, 0), 0.6, 0.6, angle=0, theta1=90, theta2=270, color='crimson', lw=2)
    arc_end = patches.Arc((float(L), 0), 0.6, 0.6, angle=0, theta1=-90, theta2=90, color='crimson', lw=2)
    ax.add_patch(arc_start)
    ax.add_patch(arc_end)
    ax.plot(0, 0.3, marker='<', color='crimson', markersize=5) 
    ax.plot(float(L), 0.3, marker='>', color='crimson', markersize=5) 

    # 4. BOX VALORI NUMERICI
    testo_0 = f"s=0\nN: {N0}\nT: {T0}\nM: {M0}"
    ax.text(-0.6, 0, testo_0, ha='right', va='center', 
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='gray', alpha=0.9), 
            fontsize=9, family='monospace')
            
    testo_L = f"s={float(L):.2f}\nN: {NL}\nT: {TL}\nM: {ML}"
    ax.text(float(L)+0.6, 0, testo_L, ha='left', va='center', 
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec='gray', alpha=0.9), 
            fontsize=9, family='monospace')

    # Asse locale s
    ax.annotate('', xy=(float(L)*0.55, -0.6), xytext=(float(L)*0.45, -0.6),
                arrowprops=dict(edgecolor='purple', arrowstyle='->', lw=2))
    ax.text(float(L)*0.5, -0.8, 'asse s', color='purple', weight='bold', fontsize=9, ha='center')

    ax.set_xlim(-1.8, float(L)+1.8)
    ax.set_ylim(-1.2, 1.2)
    ax.margins(x=0.12, y=0.22)
    ax.axis('off')
    
    return _salva_figura_png(fig, out_dir, filename, pad_inches=0.15)
    return filepath.replace("\\", "/").lstrip("/")

def _coeff_poly(expr, s_sym):
    """Coefficienti float (grado decrescente) di un'espressione polinomiale in s, oppure None se non lo e'."""
    try:
        e = sp.sympify(expr)
        e = e.xreplace({n: 0.0 for n in e.atoms(sp.Float) if abs(n) < 1e-10})
        e = e.subs({sym: 0 for sym in e.free_symbols if sym != s_sym})
        if e == 0: return [0.0]
        return [float(c) for c in sp.Poly(sp.expand(e), s_sym).all_coeffs()]
    except Exception:
        return None

def _radici_in(coeffs, L):
    """Radici reali di un polinomio (coefficienti a grado decrescente) cadenti in [0, L]."""
    if coeffs is None: return []
    c = np.trim_zeros(np.array(coeffs, dtype=float), 'f')
    if len(c) < 2: return []
    r = np.roots(c)
    return [float(x.real) for x in r if abs(x.imag) < 1e-9 * max(1.0, abs(x)) and -1e-9 <= x.real <= L + 1e-9]

def candidati_s(M_expr, N_expr, L, pesi=None, n_base=41):
    """Ascisse in cui puo' cadere il massimo di wN*|N| + wM*|M|: estremi, radici di M e N (dove i moduli cambiano
    ramo) e punti stazionari di ogni ramo +-wN*N +-wM*M, piu' una griglia di sicurezza. Cosi' il massimo e' esatto
    (M e N sono polinomi in s) e non dipende dal passo di campionamento.
    pesi = [(wN, wM), ...]: nella verifica (1/A, 1/W); nel dimensionamento la sezione e' incognita e si coprono
    piu' rapporti wN/wM a decadi."""
    s_sym = sp.Symbol('s')
    Lf = float(L)
    cM, cN = _coeff_poly(M_expr, s_sym), _coeff_poly(N_expr, s_sym)
    pts = list(np.linspace(0.0, Lf, n_base))
    pts += _radici_in(cM, Lf) + _radici_in(cN, Lf)
    def der(c):
        if c is None or len(c) < 2: return [0.0]
        return list(np.polyder(np.array(c, dtype=float)))
    dM, dN = der(cM), der(cN)
    if pesi is None:
        pesi = [(w, 1.0) for w in (0.0, 1e-3, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3)]
    def comb(a, b, wa, wb):
        a = np.array(a, dtype=float); b = np.array(b, dtype=float)
        n = max(len(a), len(b))
        a = np.pad(a, (n - len(a), 0)); b = np.pad(b, (n - len(b), 0))
        return list(wa * a + wb * b)
    for wN, wM in pesi:
        for sN in (1.0, -1.0):
            for sM in (1.0, -1.0):
                pts += _radici_in(comb(dN, dM, sN * wN, sM * wM), Lf)
    return np.unique(np.clip(np.array(pts, dtype=float), 0.0, Lf))

def _valuta_poly(expr, s_vals):
    """Valuta un'espressione in s su un array, azzerando rumore e simboli residui."""
    s_sym = sp.Symbol('s')
    e = sp.sympify(expr)
    e = e.xreplace({n: 0.0 for n in e.atoms(sp.Float) if abs(n) < 1e-10})
    e = e.subs({sym: 0 for sym in e.free_symbols if sym != s_sym})
    if not e.free_symbols:
        return np.full_like(s_vals, float(e.evalf()), dtype=float)
    f = sp.lambdify(s_sym, e, 'numpy')
    y = f(s_vals)
    if np.isscalar(y): return np.full_like(s_vals, float(y), dtype=float)
    return np.real(np.asarray(y, dtype=complex)).astype(float)

def _comb_poly(a, b_, wa=1.0, wb=1.0):
    """wa*a + wb*b per polinomi (coefficienti a grado decrescente) di grado diverso."""
    a = np.array(a, dtype=float); b_ = np.array(b_, dtype=float)
    n = max(len(a), len(b_))
    return np.pad(a, (n - len(a), 0)) * wa + np.pad(b_, (n - len(b_), 0)) * wb

def _sigma_max_numerico(barre, X1, A, W, solo_flessione):
    """Massimo di sigma = |N|/A + |M|/W (o |M|/W) su tutta la struttura per una sezione (A, W) e un X1 NUMERICI.
    barre = [(id, L, M0c, M1c, N0c, N1c)] con i coefficienti polinomiali in s. Punti candidati esatti (estremi, radici,
    punti stazionari dei quattro rami +-wN*N +-wM*M) piu' una griglia di sicurezza.
    Ritorna (sigma_max, id asta, s, N, M)."""
    wN = 0.0 if solo_flessione else 1.0 / A
    wM = 1.0 / W
    best = (0.0, None, 0.0, 0.0, 0.0)
    for eid, L, M0c, M1c, N0c, N1c in barre:
        Mc, Nc = _comb_poly(M0c, M1c, 1.0, X1), _comb_poly(N0c, N1c, 1.0, X1)
        pts = list(np.linspace(0.0, L, 41)) + _radici_in(list(Mc), L) + _radici_in(list(Nc), L)
        dM = np.polyder(Mc) if len(Mc) > 1 else np.array([0.0])
        dN = np.polyder(Nc) if len(Nc) > 1 else np.array([0.0])
        for sN in (1.0, -1.0):
            for sM in (1.0, -1.0):
                pts += _radici_in(list(_comb_poly(dN, dM, sN * wN, sM * wM)), L)
        s = np.clip(np.unique(np.array(pts, dtype=float)), 0.0, L)
        Mv, Nv = np.polyval(Mc, s), np.polyval(Nc, s)
        sig = wN * np.abs(Nv) + wM * np.abs(Mv)
        i = int(np.argmax(sig))
        if sig[i] > best[0] * (1 + 1e-9):
            best = (float(sig[i]), eid, float(s[i]), float(Nv[i]), float(Mv[i]))
    return best

def _intervalli_soddisfatti(g, d_lo, d_hi, op, n_grid=600):
    """Insieme delle dimensioni d in [d_lo, d_hi] che soddisfano la richiesta sul margine, dove g(d) = sigma_max(d) - sigma_t.
    MS >= / > t  <=>  g <= 0 ;  MS <= / < t  <=>  g >= 0 ;  MS = t  <=>  g = 0 (radici).
    Ritorna (radici, intervalli): intervalli = [(a, b)] con a/b = None se illimitato dal lato della griglia."""
    ds = np.geomspace(d_lo, d_hi, n_grid)
    gs = np.array([g(x) for x in ds])
    radici = []
    for i in range(len(ds) - 1):
        if gs[i] == 0.0:
            radici.append(float(ds[i]))
        elif gs[i] * gs[i + 1] < 0:
            lo, hi, glo = float(ds[i]), float(ds[i + 1]), gs[i]
            for _ in range(80):                       # bisezione
                mid = 0.5 * (lo + hi)
                gm = g(mid)
                if gm * glo <= 0: hi = mid
                else: lo, glo = mid, gm
            radici.append(0.5 * (lo + hi))
    if gs[-1] == 0.0: radici.append(float(ds[-1]))
    if op == "=":
        return radici, []
    ok = (lambda v: v <= 0) if op in (">=", ">") else (lambda v: v >= 0)
    estremi = [d_lo] + sorted(radici) + [d_hi]
    intervalli = []
    for a_, b_ in zip(estremi[:-1], estremi[1:]):
        if b_ <= a_: continue
        if ok(g(math.sqrt(a_ * b_))):
            if intervalli and intervalli[-1][1] == a_:
                intervalli[-1] = (intervalli[-1][0], b_)
            else:
                intervalli.append((a_, b_))
    out = [(None if a_ == d_lo else a_, None if b_ == d_hi else b_) for a_, b_ in intervalli]
    return radici, out

def _risolvi_struttura_plv_core(dati_json, out_dir, incognita):
    try:
        dati_json = normalizza_dati(dati_json)
    except (KeyError, ValueError, TypeError) as exc_:
        raise ErroreSolutore(_box_errore("Dati di input non validi", f"Campo mancante o non numerico nel JSON: {exc_!r}"))
    _errori_dati = valida_dati(dati_json)
    if _errori_dati:
        raise ErroreSolutore(_box_errore("Dati di input non validi", "<br>".join("&bull; " + m_ for m_ in _errori_dati)))
    if isinstance(incognita, dict):
        nodo_svinc = incognita.get('node')
        tipo_svinc = incognita.get('type')
    else:
        nodo_svinc = incognita[0]
        tipo_svinc = incognita[1]

    report = []
    
    nodi = {int(n['id']): (ex(n['x']), ex(n['y'])) for n in dati_json.get('nodes', [])}
    elementi = [e for e in dati_json.get('elements', [])]
    
    adj = {n_id: [] for n_id in nodi}
    for el in elementi:
        adj[int(el['n1'])].append(int(el['n2']))
        adj[int(el['n2'])].append(int(el['n1']))
    # --- Carichi distribuiti generalizzati (ortogonali all'asta e/o con componenti globali x/y) ---
    def _phi(tipo_q, t, L_t):
        """Profilo di intensita' lungo l'asta (1 = valore massimo)."""
        if tipo_q == 'triangolare_crescente': return t / L_t
        if tipo_q == 'triangolare_decrescente': return 1 - t / L_t
        return sp.Integer(1)

    def _carico(el_load):
        """(asta, L, ux, uy, wx, wy, tipo): (wx, wy) = intensita' MASSIMA per unita' di lunghezza d'asta, componenti
        globali ESATTE. q_trans>0 punta verso la normale destra (giu' per un'asta orizzontale da sinistra a destra);
        q_x, q_y sono componenti globali; riferimento 'proiezione' = intensita' per unita' di lunghezza proiettata."""
        target = next((e for e in elementi if e['id'] == el_load.get('elem_id')), None)
        if target is None: return None
        x1_t, y1_t = nodi[int(target['n1'])]
        x2_t, y2_t = nodi[int(target['n2'])]
        L_t = hyp(x2_t - x1_t, y2_t - y1_t)
        if L_t == 0: return None
        ux, uy = (x2_t - x1_t) / L_t, (y2_t - y1_t) / L_t
        qt, qx, qy = ex(el_load.get('q_trans', 0)), ex(el_load.get('q_x', 0)), ex(el_load.get('q_y', 0))
        wx, wy = qt * uy, -qt * ux
        if qx != 0 or qy != 0:
            f = sp.Integer(1)
            if el_load.get('riferimento', 'lunghezza') == 'proiezione':
                f = sp.Abs(ux * qy - uy * qx) / sp.sqrt(qx ** 2 + qy ** 2)
            wx, wy = wx + f * qx, wy + f * qy
        if wx == 0 and wy == 0: return None
        return target, L_t, ux, uy, wx, wy, el_load.get('tipo', 'costante')

    def _risultante(c):
        """(Rx, Ry, X_baricentro, Y_baricentro) ESATTI del carico distribuito."""
        target, L_t, ux, uy, wx, wy, tipo_q = c
        x1_t, y1_t = nodi[int(target['n1'])]
        x2_t, y2_t = nodi[int(target['n2'])]
        if tipo_q == 'triangolare_crescente': area, k = sp.Rational(1, 2), sp.Rational(2, 3)
        elif tipo_q == 'triangolare_decrescente': area, k = sp.Rational(1, 2), sp.Rational(1, 3)
        else: area, k = sp.Integer(1), sp.Rational(1, 2)
        return wx * L_t * area, wy * L_t * area, x1_t + (x2_t - x1_t) * k, y1_t + (y2_t - y1_t) * k
    
    # ==========================================
    # FASE 1: DATI, SEZIONE E ANALISI CINEMATICA
    # ==========================================
    num_nodi = len(nodi)
    num_aste = len(elementi)
    
    # --- RICONOSCIMENTO DINAMICO (LOGICA MULTI-CORPO) ---
    def _errore_cerniera(titolo, msg):
        report.append(_box_errore(titolo, msg))
        raise ErroreSolutore("\n".join(report))

    def _info_cerniere():
        """Per ogni cerniera interna restituisce {'nodo', 'rami': [(id_asta, set_nodi_ramo), ...], 'corpi_aggiunti'}.

        Una cerniera sul nodo N con d aste incidenti puo' svincolare solo alcune di esse (campo opzionale
        'elems': id delle aste incernierate in N; le altre restano incastrate fra loro e formano un solo corpo).
          - |H| < d : le |H| aste incernierate sono corpi separati  -> |H| corpi aggiunti, |H| equazioni M=0 (una per ramo)
          - |H| = d : cerniera comune a d corpi                     -> d-1 corpi aggiunti, d-1 equazioni M=0
        Senza 'elems' e' ammesso solo il caso d = 2 (cerniera tra due aste): con d > 2 il corpo da svincolare
        sarebbe ambiguo e il calcolo verrebbe sbagliato in silenzio, quindi si blocca con un errore esplicito."""
        info = []
        for cerniera in dati_json.get('internal_hinges', []):
            nc = int(cerniera['node_id'])
            if nc not in nodi:
                _errore_cerniera("Cerniera interna non valida", f"Il nodo {nc} indicato per la cerniera interna non esiste.")
            incid = [el for el in elementi if int(el['n1']) == nc or int(el['n2']) == nc]
            d = len(incid)
            if d < 2:
                _errore_cerniera("Cerniera interna non valida",
                                 f"Il nodo {nc} ha {d} asta/e incidente/i: una cerniera interna richiede almeno 2 aste.")
            elems = cerniera.get('elems')
            if elems:
                ids = {int(x) for x in elems}
                ids_inc = {int(el['id']) for el in incid}
                if not ids <= ids_inc:
                    _errore_cerniera("Cerniera interna non valida",
                                     f"Nel nodo {nc} le aste {sorted(ids - ids_inc)} non sono incidenti (aste del nodo: {sorted(ids_inc)}).")
                hinged = [el for el in incid if int(el['id']) in ids]
            else:
                if d > 2:
                    _errore_cerniera("Cerniera interna ambigua",
                                     f"Il nodo {nc} ha {d} aste incidenti ({sorted(int(e['id']) for e in incid)}): indica nel campo "
                                     f"<b>elems</b> della cerniera le aste che sono realmente incernierate in questo nodo "
                                     f"(es. {{\"node_id\": {nc}, \"elems\": [{int(incid[0]['id'])}]}}). "
                                     "Senza questa informazione le equazioni ausiliarie verrebbero scritte sul corpo sbagliato.")
                hinged = [incid[0]]
            rami_el = hinged[:-1] if len(hinged) == d else hinged
            # 'carichi_su': id dell'asta incernierata il cui corpo possiede coppie concentrate, vincoli (anche rotazionali) e carichi
            # applicati nel nodo-cerniera. Senza il campo restano al corpo collegato rigidamente (o, se tutte le aste sono incernierate,
            # all'ultima asta incidente). Forze Fx, Fy nel nodo non cambiano le equazioni (braccio nullo): cambiano le coppie,
            # i momenti reattivi e le forze trasmesse dalla cerniera.
            attach = None
            if cerniera.get('carichi_su') not in (None, ''):
                att_id = int(cerniera['carichi_su'])
                if att_id not in {int(el['id']) for el in hinged}:
                    _errore_cerniera("Cerniera interna non valida",
                                     f"Nel nodo {nc} 'carichi_su' = {att_id} non e' una delle aste incernierate "
                                     f"({sorted(int(el['id']) for el in hinged)}).")
                if att_id in {int(el['id']) for el in rami_el}:
                    attach = att_id
            rami = []
            for el in rami_el:
                lontano = int(el['n2']) if int(el['n1']) == nc else int(el['n1'])
                visitati, stack, nodi_ramo = {nc}, [lontano], set()
                while stack:
                    cur = stack.pop()
                    if cur not in visitati:
                        visitati.add(cur)
                        nodi_ramo.add(cur)
                        for vic in adj[cur]:
                            if vic not in visitati: stack.append(vic)
                # anello chiuso: dal ramo si rientra nel nodo cerniera attraverso un'altra asta
                for el2 in incid:
                    if el2 is el: continue
                    l2 = int(el2['n2']) if int(el2['n1']) == nc else int(el2['n1'])
                    if l2 in nodi_ramo:
                        _errore_cerniera("Telaio con anello chiuso non supportato",
                                         f"Il ramo dell'asta {el['id']} si ricollega al nodo {nc} tramite l'asta {el2['id']}: "
                                         "la cerniera non separa la struttura in corpi distinti (anello chiuso).")
                rami.append((int(el['id']), nodi_ramo))
            if cerniera.get('carichi_su') not in (None, ''):
                carrier = int(cerniera['carichi_su'])
            elif len(hinged) == d:
                carrier = int(hinged[-1]['id'])          # tutte incernierate: il corpo "base" e' l'ultima asta incidente
            else:
                carrier = None                           # corpo collegato rigidamente al nodo
            info.append({'nodo': nc, 'rami': rami, 'corpi_aggiunti': len(rami_el), 'attach': attach, 'carrier': carrier,
                         'attach_dichiarato': cerniera.get('carichi_su') not in (None, '')})
        return info

    # --- Controllo anelli chiusi: le sezioni di taglio (calcola_sollecitazioni) sono valide solo per telai ad albero.
    # Un anello senza cerniera darebbe risultati sbagliati in silenzio.
    _padre = {n_: n_ for n_ in nodi}
    def _trova(a_):
        while _padre[a_] != a_:
            _padre[a_] = _padre[_padre[a_]]
            a_ = _padre[a_]
        return a_
    for el_ in elementi:
        ra_, rb_ = _trova(int(el_['n1'])), _trova(int(el_['n2']))
        if ra_ == rb_:
            _errore_cerniera("Telaio con anello chiuso non supportato",
                             f"L'asta {el_['id']} chiude un anello nella struttura: il solutore tratta telai ad albero "
                             "(ogni asta separa la struttura in due parti). Spezza l'anello con una cerniera/nodo doppio oppure "
                             "modella la struttura come aste collegate a terra.")
        _padre[ra_] = rb_

    cerniere_info = _info_cerniere()

    def _del_ramo(nodo, nc_, id_ramo, nodi_ramo_, attach_):
        """True se le azioni applicate al nodo appartengono al corpo del ramo: nodi del ramo, oppure il nodo-cerniera stesso
        quando 'carichi_su' lo assegna a quel ramo."""
        return nodo in nodi_ramo_ or (nodo == nc_ and attach_ is not None and attach_ == id_ramo)
    cerniere_interne = dati_json.get('internal_hinges', [])
    for _p in dati_json.get('point_loads', []):
        if _p.get('M', 0) != 0 and any(int(_c['node_id']) == int(_p['node_id']) and _c.get('carichi_su') in (None, '')
                                       for _c in cerniere_interne):
            report.append(_box_avviso("Coppia applicata in una cerniera interna",
                f"Nel Nodo {_p['node_id']} c'&egrave; sia una cerniera interna sia una coppia concentrata: la cerniera non trasmette momento, "
                "quindi la coppia viene assegnata al corpo collegato rigidamente al nodo (non alle aste incernierate). "
                "Se nella prova &egrave; applicata a un altro corpo, indica nella cerniera il campo <b>carichi_su</b> (id dell'asta di quel corpo)."))
    num_cerniere_interne = len(cerniere_interne)
    corpi_aggiunti = sum(ci['corpi_aggiunti'] for ci in cerniere_info)

    # 1. Tronchi: 1 iniziale + 1 per ogni corpo aggiunto dalle cerniere (cerniera tra 2 aste = 1; cerniera comune a k corpi = k-1)
    num_tronchi = 1 + corpi_aggiunti
    gdl_liberi = 3 * num_tronchi

    # 2. Vincoli interni: ogni corpo aggiunto e' collegato da una cerniera (2 traslazioni relative bloccate)
    vincoli_interni = 2 * corpi_aggiunti

    # --- Calcolo Proprietà Geometriche e di Materiale ---
    mat = dati_json.get('material', {})
    sec = dati_json.get('section', {})
    
    E = float(mat.get('E', 2.1e11))  # solo per la stampa; i calcoli usano i valori esatti
    b = float(sec.get('b', 0.2))
    h = float(sec.get('h', 0.4))
    design_mode = str(sec.get('design', 'none'))
    modalita = "dimensionamento" if design_mode in ('h', 'b', 'ratio') else "verifica"
    ratio_k = float(sec.get('ratio', 0) or 0)
    sigma_amm_p1 = float(mat.get('sigma_amm', 235e6))
    
    A = b * h
    I = (b * h**3) / 12
    EI = E * I
    EI_assegnato = 'EI' in mat
    if EI_assegnato:
        EI = float(mat['EI'])
    
    vincoli_reali = 0
    dettagli_vincoli = []
    
    for s in dati_json.get('supports', []):
        nid = int(s.get('node_id'))
        fix_u = s.get('fix_u', False)
        fix_v = s.get('fix_v', False)
        fix_w = s.get('fix_w', False)
        
        v_count = 0
        bloccati = []
        reazioni = []
        nome_vincolo = "Vincolo generico"
        
        if fix_u: 
            v_count += 1
            bloccati.append("trasl. X")
            reazioni.append("R<sub>x</sub>")
        if fix_v: 
            v_count += 1
            bloccati.append("trasl. Y")
            reazioni.append("R<sub>y</sub>")
        if fix_w: 
            v_count += 1
            bloccati.append("rotazione")
            reazioni.append("M")
            
        # Componente lungo la normale del vincolo inclinato (carrello inclinato / pendolo / bipendolo inclinato)
        if s.get('rotated_roller'):
            v_count += 1
            if s.get('pendolo') or s.get('bipendolo'):
                bloccati.append(f"trasl. lungo l'asse delle bielle ({float(s.get('axis_angle', 0)):g}°)")
            else:
                bloccati.append(f"trasl. ortogonale ({float(s.get('angle', 0)):g}°)")
            reazioni.append("R<sub>n</sub>")

        # Classificazione didattica del vincolo
        if s.get('pendolo'):
            nome_vincolo = "Pendolo (biella inclinata)"
        elif s.get('bipendolo'):
            nome_vincolo = "Doppio pendolo / bipendolo inclinato"
        elif s.get('rotated_roller'):
            nome_vincolo = "Carrello Inclinato" + (" con blocco rotazione (bipendolo)" if fix_w else "")
        elif fix_u and fix_v and fix_w:
            nome_vincolo = "Incastro"
        elif fix_u and fix_v and not fix_w:
            nome_vincolo = "Cerniera"
        elif (fix_u and not fix_v and not fix_w) or (not fix_u and fix_v and not fix_w):
            nome_vincolo = "Carrello"
        elif (fix_u or fix_v) and fix_w:
            nome_vincolo = "Bipendolo"

        dettagli_vincoli.append(f"<li><b>Nodo {nid}</b>: {nome_vincolo} &mdash; blocca {', '.join(bloccati)} &rarr; genera ({', '.join(reazioni)}) (<b>s = {v_count}</b>)</li>")
        vincoli_reali += v_count

    # --- STAMPA REPORT E VERIFICA CINEMATICA ---
    vincoli_totali = vincoli_reali + vincoli_interni
    grado_iperstaticita = vincoli_totali - gdl_liberi

    # Verifica rigorosa di non labilita' sull'intera struttura (rango matrice equilibrio con tutti i vincoli attivi)
    struttura_non_labile = True
    if grado_iperstaticita >= 0:
        try:
            _supports_chk = dati_json.get('supports', [])
            _vars_chk, _eqs_chk = [], []
            _polo_chk = list(nodi.keys())[0]
            _xp, _yp = nodi[_polo_chk]
            _eqx_c, _eqy_c, _eqm_c = sp.Integer(0), sp.Integer(0), sp.Integer(0)
            for _is, _sc in enumerate(_supports_chk):
                _nidc = int(_sc.get('node_id'))
                _xc_n, _yc_n = nodi[_nidc]
                if _sc.get('fix_u'):
                    _sy = sp.Symbol(f'_cH_{_is}')
                    _vars_chk.append(_sy); _eqx_c += _sy; _eqm_c -= _sy * (_yc_n - _yp)
                if _sc.get('fix_v'):
                    _sy = sp.Symbol(f'_cV_{_is}')
                    _vars_chk.append(_sy); _eqy_c += _sy; _eqm_c += _sy * (_xc_n - _xp)
                if _sc.get('fix_w'):
                    _sy = sp.Symbol(f'_cM_{_is}')
                    _vars_chk.append(_sy); _eqm_c += _sy
                if _sc.get('rotated_roller'):
                    _sy = sp.Symbol(f'_cRn_{_is}')
                    _ca, _sa = cs_ang(ex(_sc.get('angle', 0)) + 90)
                    _vars_chk.append(_sy); _eqx_c += _sy * _ca; _eqy_c += _sy * _sa
                    _eqm_c += (_sy * _sa) * (_xc_n - _xp) - (_sy * _ca) * (_yc_n - _yp)
            _eqs_chk = [_eqx_c, _eqy_c, _eqm_c]
            for _ci in cerniere_info:
                _nc = _ci['nodo']; _xc, _yc = nodi[_nc]
                for _id_ar, _n_ramo in _ci['rami']:
                    _eqa = sp.Integer(0)
                    for _is, _sc in enumerate(_supports_chk):
                        _nidc = int(_sc.get('node_id'))
                        if _del_ramo(_nidc, _nc, _id_ar, _n_ramo, _ci['attach']):
                            _xc_n, _yc_n = nodi[_nidc]
                            if _sc.get('fix_u'): _eqa -= sp.Symbol(f'_cH_{_is}') * (_yc_n - _yc)
                            if _sc.get('fix_v'): _eqa += sp.Symbol(f'_cV_{_is}') * (_xc_n - _xc)
                            if _sc.get('fix_w'): _eqa += sp.Symbol(f'_cM_{_is}')
                            if _sc.get('rotated_roller'):
                                _ca, _sa = cs_ang(ex(_sc.get('angle', 0)) + 90)
                                _sy = sp.Symbol(f'_cRn_{_is}')
                                _eqa += (_sy * _sa) * (_xc_n - _xc) - (_sy * _ca) * (_yc_n - _yc)
                    _eqs_chk.append(_eqa)
            _A_full, _ = sp.linear_eq_to_matrix(_eqs_chk, _vars_chk)
            struttura_non_labile = (_A_full.rank() == len(_eqs_chk))
        except Exception:
            struttura_non_labile = True

    # ==========================================
    # INTESTAZIONE INIZIALE: BOX RIEPILOGATIVO VERIFICA CINEMATICA (CIMA PAGINA 1)
    # ==========================================
    cern_desc_txt = f"{num_cerniere_interne} {'cerniera interna' if num_cerniere_interne == 1 else 'cerniere interne'}" if num_cerniere_interne else "nessuna cerniera interna"
    report.append("<div style='background-color: #eef4ff; padding: 12px 15px; border: 1px solid #b6d4fe; border-left: 5px solid #0d6efd; border-radius: 4px; margin-bottom: 16px;'>")
    report.append("<h3 style='margin: 0 0 8px 0; color: #0a58ca;'>Verifica Cinematica</h3>")
    report.append("<ul style='margin: 4px 0 8px 0;'>")
    report.append(f"<li><b>Gradi di libert&agrave;:</b> 3 &times; {num_tronchi} {'tronco' if num_tronchi == 1 else 'tronchi'} = <b>{gdl_liberi}</b> ({cern_desc_txt})</li>")
    report.append(f"<li><b>Gradi di vincolo:</b> v<sub>e</sub> = {vincoli_reali} (esterni) + v<sub>i</sub> = {vincoli_interni} (cerniere interne) &rArr; <b>v<sub>tot</sub> = {vincoli_totali}</b></li>")
    if grado_iperstaticita == 0:
        stat_msg = "La struttura &egrave; <b>isostatica</b>."
    elif grado_iperstaticita == 1:
        stat_msg = "La struttura &egrave; <b>1 volta iperstatica</b>."
    elif grado_iperstaticita > 1:
        stat_msg = f"La struttura &egrave; <b>{grado_iperstaticita} volte iperstatica</b>."
    else:
        stat_msg = "La struttura &egrave; <b>labile</b> (h &lt; 0)."
    report.append(f"<li><b>Grado di iperstaticit&agrave;:</b> h = v<sub>tot</sub> &minus; gdl = {vincoli_totali} &minus; {gdl_liberi} = <b>{grado_iperstaticita}</b> &rarr; {stat_msg}</li>")
    if grado_iperstaticita >= 0 and struttura_non_labile:
        report.append("<li><b>Verifica di non labilit&agrave;:</b> sistema statico ben posto / non singolare (vincoli correttamente disposti).</li>")
    report.append("</ul></div>")

    if grado_iperstaticita >= 0 and not struttura_non_labile:
        report.append(_box_errore("Analisi Interrotta (Struttura Labile)",
                                  "I vincoli inseriti non impediscono tutti gli atti di moto rigidi (vincoli mal disposti, ad es. paralleli o concorrenti)."))
        return "\n".join(report), "", ""
    if grado_iperstaticita > 1:
        report.append("<div style='background-color: #fff3cd; padding: 15px; border-left: 5px solid #ffc107; margin-bottom: 20px;'>")
        report.append("<h4 style='color: #856404; margin-top: 0;'>Limite del Solutore Superato</h4>")
        report.append(f"<p style='color: #856404;'>Questa versione del programma supporta un massimo di <b>1 grado di iperstaticit&agrave;</b>. La struttura richiede <b>{grado_iperstaticita}</b> equazioni di congruenza simultanee.</p></div>")
        return "\n".join(report), "", ""
    if grado_iperstaticita < 0:
        report.append(_box_errore("Analisi Interrotta (Struttura Labile)",
                                  "La struttura non &egrave; vincolata in modo sufficiente o efficace. &Egrave; impossibile procedere con l'analisi statica."))
        return "\n".join(report), "", ""

    # Schema generale della struttura reale (generato subito per Sezione 1)
    supports_n = dati_json.get('supports', [])
    percorso_S1 = "plv_S0_S1.png"
    percorso_diag = "plv_diagrammi.png"
    img_schema_gen = plot_schema_generale(
        nodi, elementi, supports_n, cerniere_interne,
        dati_json.get('point_loads', []),
        dati_json.get('element_loads', []),
        dati_json.get('thermal_loads', []),
        out_dir, percorso_S1
    )

    # ==========================================
    # SEZIONE 1: DATI DELLA STRUTTURA E CLASSIFICAZIONE CINEMATICA
    # ==========================================
    report.append("<h2>1. Dati della struttura e classificazione cinematica</h2>")
    cern_pos_txt = (
        "; cerniere interne nei nodi: <b>" + ", ".join(str(ci_['nodo']) for ci_ in cerniere_info) + "</b>"
        if cerniere_info else "; nessuna cerniera interna"
    )
    report.append(f"<p>Telaio piano composto da <b>{num_aste} aste</b> e <b>{num_nodi} nodi</b>{cern_pos_txt}.</p>")
    _w_gen, _h_gen = _scala_dim_png(img_schema_gen, 580, 320)
    _dim_gen = f"width='{_w_gen}' height='{_h_gen}'" if _h_gen else f"width='{_w_gen}'"
    report.append(f"<div style='text-align: center; line-height: 1.0; font-size: 1px; margin: 6px 0 10px 0;'><img src='{_to_file_uri(img_schema_gen)}' {_dim_gen}></div>")

    report.append("<div style='page-break-inside: avoid;'><p><b>Analisi dei vincoli esterni:</b></p><ul>")
    for det in dettagli_vincoli:
        report.append(det)
    report.append("</ul></div>")

    report.append("<div style='background-color: #f8f9fa; padding: 10px 14px; border-left: 4px solid #198754; margin-bottom: 16px;'>")
    report.append("<h4 style='margin: 0 0 6px 0;'>Propriet&agrave; Geometriche e del Materiale</h4><ul style='margin: 4px 0;'>")
    W_f_p1 = (b * (h**2)) / 6
    gpa_note = f" ({fmt_smart(E / 1e9)} GPa)" if E >= 1e9 else ""
    report.append(f"<li>Modulo di Young (E): <b>{fmt_sci_html(E, sig=4)} Pa</b>{gpa_note}</li>")
    report.append(f"<li>Tensione ammissibile (&sigma;<sub>amm</sub>): <b>{fmt_smart(sigma_amm_p1 / 1e6)} MPa</b> ({fmt_sci_html(sigma_amm_p1, sig=4)} Pa)</li>")
    sez_desc = "Dimensioni sezione assegnata (Punto 3)" if modalita == "dimensionamento" else "Dimensioni sezione rettangolare"
    report.append(f"<li>{sez_desc}: b = <b>{fmt_smart(b)} m</b>, h = <b>{fmt_smart(h)} m</b></li>")
    report.append(f"<li>Area (A = b&middot;h): <b>{fmt_smart(A, min_sig=4)} m<sup>2</sup></b></li>")
    report.append(f"<li>Momento d'Inerzia (I = b&middot;h<sup>3</sup>/12): <b>{fmt_sci_html(I, sig=4)} m<sup>4</sup></b></li>")
    report.append(f"<li>Modulo di resistenza flessionale (W = b&middot;h<sup>2</sup>/6): <b>{fmt_sci_html(W_f_p1, sig=4)} m<sup>3</sup></b></li>")
    ei_tag = " (assegnata)" if EI_assegnato else ""
    report.append(f"<li>Rigidezza flessionale (EI){ei_tag}: <b>{fmt_sci_html(EI, sig=4)} N&middot;m<sup>2</sup></b></li>")
    if modalita == "dimensionamento":
        if design_mode == 'ratio':
            report.append(f"<li>Parametro di dimensionamento (Punto 4): rapporto di forma <b>h / b = {fmt_smart(ratio_k)}</b> (b incognita)</li>")
        elif design_mode == 'h':
            report.append(f"<li>Parametro di dimensionamento (Punto 4): base fissata <b>b = {fmt_smart(b)} m</b> (altezza h incognita)</li>")
        elif design_mode == 'b':
            report.append(f"<li>Parametro di dimensionamento (Punto 4): altezza fissata <b>h = {fmt_smart(h)} m</b> (base b incognita)</li>")
    report.append("</ul></div>")

    # ==========================================
    # FASE 2: SCELTA INCOGNITA E SISTEMA PRINCIPALE
    # ==========================================
    reazioni_simboliche = []
    dip_ced = dip_term = False
    supports_n = dati_json.get('supports', [])

    sym_Fn, _cont_Fn = {}, {}
    for i_s, s_ in enumerate(supports_n):
        if s_.get('rotated_roller'):
            nid_ = int(s_['node_id'])
            _cont_Fn[nid_] = _cont_Fn.get(nid_, 0) + 1
            sym_Fn[i_s] = sp.Symbol(f'Rn_{nid_}' if _cont_Fn[nid_] == 1 else f'Rn_{nid_}_{_cont_Fn[nid_]}')

    ang_X1 = sp.Integer(0)
    idx_X1_sup = None
    if grado_iperstaticita > 0:
        if tipo_svinc == 'Fn':
            cand = [i for i, s_ in enumerate(supports_n) if int(s_['node_id']) == nodo_svinc and s_.get('rotated_roller')]
            if not cand:
                _errore_cerniera("Incognita iperstatica non valida",
                                 f"Il Nodo {nodo_svinc} non ha un carrello inclinato / pendolo / bipendolo inclinato: scegli un'altra incognita.")
            idx_X1_sup = cand[0]
            if isinstance(incognita, dict) and incognita.get('sup_idx') is not None and int(incognita['sup_idx']) in cand:
                idx_X1_sup = int(incognita['sup_idx'])
            elif len(cand) > 1 and isinstance(incognita, dict) and incognita.get('angle') is not None:
                for i_c in cand:
                    if abs(float(ex(supports_n[i_c].get('angle', 0))) - float(incognita['angle'])) < 1e-6:
                        idx_X1_sup = i_c
                        break
            ang_X1 = ex(supports_n[idx_X1_sup].get('angle', 0))
        elif tipo_svinc in ('Fx', 'Fy', 'M'):
            chiave_fix = {'Fx': 'fix_u', 'Fy': 'fix_v', 'M': 'fix_w'}[tipo_svinc]
            if not any(int(s_['node_id']) == nodo_svinc and s_.get(chiave_fix) for s_ in supports_n):
                _errore_cerniera("Incognita iperstatica non valida",
                                 f"Nel Nodo {nodo_svinc} non esiste la reazione {tipo_svinc} scelta come incognita.")
    X1_fx, X1_fy, X1_m = sp.Integer(0), sp.Integer(0), sp.Integer(0)
    if tipo_svinc == 'Fx': X1_fx = sp.Integer(1)
    elif tipo_svinc == 'Fy': X1_fy = sp.Integer(1)
    elif tipo_svinc in ('M', 'M_int'): X1_m = sp.Integer(1)
    elif tipo_svinc == 'Fn': X1_fx, X1_fy = cs_ang(ang_X1 + 90)

    def _reazioni_da_vincoli(rimuovi_X1):
        out = []
        for i_s, s in enumerate(supports_n):
            nid = int(s.get('node_id'))
            if s.get('fix_u') and not (rimuovi_X1 and nid == nodo_svinc and tipo_svinc == 'Fx'):
                out.append({'nodo': nid, 'tipo': 'Fx', 'sym': sp.Symbol(f'H_{nid}'), 'sup_idx': i_s})
            if s.get('fix_v') and not (rimuovi_X1 and nid == nodo_svinc and tipo_svinc == 'Fy'):
                out.append({'nodo': nid, 'tipo': 'Fy', 'sym': sp.Symbol(f'V_{nid}'), 'sup_idx': i_s})
            if s.get('fix_w') and not (rimuovi_X1 and nid == nodo_svinc and tipo_svinc == 'M'):
                out.append({'nodo': nid, 'tipo': 'M', 'sym': sp.Symbol(f'M_{nid}'), 'sup_idx': i_s})
            if s.get('rotated_roller') and not (rimuovi_X1 and tipo_svinc == 'Fn' and i_s == idx_X1_sup):
                out.append({'nodo': nid, 'tipo': 'Fn', 'sym': sym_Fn[i_s], 'angle': s.get('angle', 0), 'sup_idx': i_s,
                            'pendolo': bool(s.get('pendolo'))})
        return out

    report.append("<h2>2. Scelta dell'incognita iperstatica e sistema principale</h2>")
    if grado_iperstaticita > 0:
        tipo_azione = {'Fx': "forza orizzontale", 'Fy': "forza verticale", 'Fn': "reazione normale/assiale"}.get(tipo_svinc, "coppia")
        sim_inc_map = {'Fx': f"H<sub>{nodo_svinc}</sub>", 'Fy': f"V<sub>{nodo_svinc}</sub>",
                       'M': f"M<sub>{nodo_svinc}</sub>", 'Fn': sym_reaz_html(sym_Fn[idx_X1_sup]) if idx_X1_sup is not None else f"R<sub>n{nodo_svinc}</sub>"}
        sim_inc_html = sim_inc_map.get(tipo_svinc, "X<sub>1</sub>")
        report.append(f"<p>Si sceglie come incognita iperstatica <b>X<sub>1</sub> = {sim_inc_html}</b> ({tipo_azione} al Nodo {nodo_svinc}), "
                      f"degradando/rimuovendo il vincolo corrispondente per ottenere il sistema isostatico principale <b>S<sub>0</sub></b>.</p>")
        report.append("<p>Formula di sovrapposizione degli effetti: <span class='formula'>S = S<sub>0</sub> + X<sub>1</sub> &middot; S<sub>1</sub></span></p>")
        reazioni_simboliche = _reazioni_da_vincoli(True)
    else:
        report.append("<p>La struttura &egrave; <b>isostatica</b>: non si rimuove alcun vincolo e il sistema principale coincide con quello reale (<b>S = S<sub>0</sub></b>).</p>")
        reazioni_simboliche = _reazioni_da_vincoli(False)

    if cerniere_info:
        report.append(
            "<div style='background-color: #f8f9fa; padding: 8px 12px; border-left: 4px solid #6f42c1; margin: 8px 0 14px 0; font-size: 13px;'>"
            "<b>Legenda reazioni interne nelle cerniere:</b> "
            "<span style='color:#0d6efd; font-weight:bold;'>&#9632; Blu</span> = reazione esercitata dalla cerniera sul tratto precedente (verso +); "
            "<span style='color:#9c36b5; font-weight:bold;'>&#9632; Viola</span> = reazione uguale ed opposta sul tratto successivo (verso &minus;)."
            "</div>"
        )

    Somma_Fx_ext = sum([p.get('Fx', 0) for p in dati_json.get('point_loads', [])])
    Somma_Fy_ext = sum([p.get('Fy', 0) for p in dati_json.get('point_loads', [])])
    polo_id = list(nodi.keys())[0]
    xp, yp = nodi[polo_id]
    Somma_M_ext = sum([p.get('M', 0) + p.get('Fy', 0)*(nodi[int(p['node_id'])][0] - xp) - p.get('Fx', 0)*(nodi[int(p['node_id'])][1] - yp)
                       for p in dati_json.get('point_loads', [])])
    for el in dati_json.get('element_loads', []):
        c = _carico(el)
        if c is None: continue
        Rx, Ry, X_bar, Y_bar = _risultante(c)
        Somma_Fx_ext += Rx; Somma_Fy_ext += Ry
        Somma_M_ext += Ry * (X_bar - xp) - Rx * (Y_bar - yp)

    variabili = [r['sym'] for r in reazioni_simboliche]

    def stop(titolo, msg):
        report.append(_box_errore(titolo, msg))
        raise ErroreSolutore("\n".join(report))

    def risolvi_sicuro(equazioni, vars, nome="S0"):
        return risolvi_sistema_lineare(equazioni, vars, nome, report)

    def _blocco_equazioni_html(righe_eq):
        rows = []
        for etic_, expr_ in righe_eq:
            rows.append(
                f"<tr><td style='border:none; padding:4px 8px 4px 0; text-align:left; vertical-align:middle; white-space:nowrap;'>&bull; {etic_} &rArr;</td>"
                f"<td style='border:none; padding:4px 0; text-align:left; vertical-align:middle;'>{formatta_equazione_html(expr_)}</td></tr>"
            )
        return "<table class='eq-table' style='border:none; width:auto; margin:4px 0; background:transparent;'>" + "".join(rows) + "</table>"

    # ==========================================
    # DECOMPOSIZIONE IN t TRONCHI RIGIDI + REAZIONI INTERNE ALLE CERNIERE
    # ==========================================
    def _decomponi_tronchi():
        el_by_id = {int(e['id']): e for e in elementi}
        hinge_det = {ci['nodo']: {eid for eid, _ in ci['rami']} for ci in cerniere_info}
        el_inc_n = {n_: [] for n_ in nodi}
        for e_ in elementi:
            el_inc_n[int(e_['n1'])].append(int(e_['id']))
            el_inc_n[int(e_['n2'])].append(int(e_['id']))

        par_e = {int(e_['id']): int(e_['id']) for e_ in elementi}
        def _find_e(x_):
            while par_e[x_] != x_:
                par_e[x_] = par_e[par_e[x_]]
                x_ = par_e[x_]
            return x_
        def _union_e(a_, b_):
            ra, rb = _find_e(a_), _find_e(b_)
            if ra != rb:
                par_e[rb] = ra

        for v_ in nodi:
            rig = [eid for eid in el_inc_n[v_] if eid not in hinge_det.get(v_, set())]
            for i_ in range(1, len(rig)):
                _union_e(rig[0], rig[i_])

        gruppi = {}
        for e_ in elementi:
            eid = int(e_['id'])
            gruppi.setdefault(_find_e(eid), []).append(eid)

        tronchi_raw = sorted(gruppi.values(), key=lambda ids: min(ids))
        elem_to_tr = {}
        tronchi = []
        sup_nodi_set = {int(s_['node_id']): s_ for s_ in supports_n}
        cern_nodi_set = {ci['nodo'] for ci in cerniere_info}

        for k_, eids in enumerate(tronchi_raw):
            eids_s = sorted(eids)
            for eid in eids_s:
                elem_to_tr[eid] = k_
            all_nds = set()
            for eid in eids_s:
                el_ = el_by_id[eid]
                all_nds.add(int(el_['n1'])); all_nds.add(int(el_['n2']))
            tronchi.append({
                'idx': k_ + 1,
                'elems': [el_by_id[eid] for eid in eids_s],
                'elem_ids': set(eids_s),
                'all_nodes': all_nds,
                'own_nodes': set(),
            })

        ci_by_node = {ci['nodo']: ci for ci in cerniere_info}
        node_owner = {}
        for v_ in nodi:
            if v_ in ci_by_node:
                ci = ci_by_node[v_]
                if ci['attach'] is not None:
                    tr_k = elem_to_tr[ci['attach']]
                elif ci['carrier'] is not None:
                    tr_k = elem_to_tr[ci['carrier']]
                else:
                    rig_e = next(eid for eid in el_inc_n[v_] if eid not in hinge_det[v_])
                    tr_k = elem_to_tr[rig_e]
            else:
                tr_k = elem_to_tr[el_inc_n[v_][0]]
            node_owner[v_] = tr_k
            tronchi[tr_k]['own_nodes'].add(v_)

        def _dir_verso(nid_s, el_):
            altro = int(el_['n2']) if int(el_['n1']) == nid_s else int(el_['n1'])
            dx = float(nodi[altro][0] - nodi[nid_s][0])
            dy = float(nodi[altro][1] - nodi[nid_s][1])
            L_ = math.hypot(dx, dy)
            return (dx / L_, dy / L_) if L_ > 1e-9 else (1.0, 0.0)

        hinge_pins = []
        for ci in cerniere_info:
            nc = ci['nodo']
            if ci['carrier'] is not None:
                eid_car = ci['carrier']
            else:
                eid_car = next(eid for eid in el_inc_n[nc] if eid not in hinge_det[nc])
            t_carrier = elem_to_tr[eid_car]
            dir_car = _dir_verso(nc, el_by_id[eid_car])
            n_rami = len(ci['rami'])
            for r_i, (id_asta_r, _) in enumerate(ci['rami'], 1):
                t_branch = elem_to_tr[id_asta_r]
                dir_br = _dir_verso(nc, el_by_id[id_asta_r])
                if t_branch <= t_carrier:
                    t_A, t_B = t_branch, t_carrier
                    dir_A, dir_B = dir_br, dir_car
                else:
                    t_A, t_B = t_carrier, t_branch
                    dir_A, dir_B = dir_car, dir_br
                suff = f"{nc}" if n_rami == 1 else f"{nc}_{r_i}"
                hinge_pins.append({
                    'nodo': nc,
                    'sym_H': sp.Symbol(f"Hc_{suff}"),
                    'sym_V': sp.Symbol(f"Vc_{suff}"),
                    't_A': t_A, 't_B': t_B,
                    'dir_A': dir_A, 'dir_B': dir_B,
                    'branch_elem': id_asta_r,
                })

        for tr in tronchi:
            # Sceglie polo naturale per equilibrio del tronco: cerniera interna > incastro/cerniera > primo nodo
            cern_tr = sorted(n_ for n_ in tr['all_nodes'] if n_ in cern_nodi_set)
            if cern_tr:
                tr['polo'] = cern_tr[0]
            else:
                sup_tr = sorted(
                    (n_ for n_ in tr['own_nodes'] if n_ in sup_nodi_set),
                    key=lambda n_: (-(int(bool(sup_nodi_set[n_].get('fix_u'))) + int(bool(sup_nodi_set[n_].get('fix_v'))) + int(bool(sup_nodi_set[n_].get('fix_w')))), n_)
                )
                tr['polo'] = sup_tr[0] if sup_tr else min(tr['all_nodes'])
            aste_s = ", ".join(str(int(e['id'])) for e in tr['elems'])
            nodi_s = ", ".join(str(n_) for n_ in sorted(tr['all_nodes']))
            tr['desc'] = f"aste {aste_s}; nodi {nodi_s}"

        return tronchi, node_owner, elem_to_tr, hinge_pins

    tronchi_sys, node_owner_tr, elem_to_tr, hinge_pins = _decomponi_tronchi()
    vars_cerniere = [s_ for hp in hinge_pins for s_ in (hp['sym_H'], hp['sym_V'])]
    cerniere_plot = []
    for hp in hinge_pins:
        cerniere_plot.append({'nodo': hp['nodo'], 'dir': hp['dir_A'], 'colore': '#0d6efd', 'lato': 'A',
                              'sign': 1, 'sym_H': hp['sym_H'], 'sym_V': hp['sym_V']})
        cerniere_plot.append({'nodo': hp['nodo'], 'dir': hp['dir_B'], 'colore': '#9c36b5', 'lato': 'B',
                              'sign': -1, 'sym_H': hp['sym_H'], 'sym_V': hp['sym_V']})

    def _costruisci_sistema_tronchi(reaz_ext, con_carichi_reali=True, carico_X1=False):
        blocchi_tr = []
        eqs_all = []
        for k_, tr in enumerate(tronchi_sys):
            pk = tr['polo']
            xpk, ypk = nodi[pk]
            eq_x_k = sp.Integer(0)
            eq_y_k = sp.Integer(0)
            eq_m_k = sp.Integer(0)

            for r in reaz_ext:
                if node_owner_tr[r['nodo']] != k_:
                    continue
                xr, yr = nodi[r['nodo']]
                if r['tipo'] == 'Fx':
                    eq_x_k += r['sym']; eq_m_k -= r['sym'] * (yr - ypk)
                elif r['tipo'] == 'Fy':
                    eq_y_k += r['sym']; eq_m_k += r['sym'] * (xr - xpk)
                elif r['tipo'] == 'M':
                    eq_m_k += r['sym']
                elif r['tipo'] == 'Fn':
                    ca, sa = cs_ang(r['angle'] + 90)
                    eq_x_k += r['sym'] * ca; eq_y_k += r['sym'] * sa
                    eq_m_k += (r['sym'] * sa) * (xr - xpk) - (r['sym'] * ca) * (yr - ypk)

            for hp in hinge_pins:
                if hp['t_A'] == k_:
                    sgn = 1
                elif hp['t_B'] == k_:
                    sgn = -1
                else:
                    continue
                xc, yc = nodi[hp['nodo']]
                eq_x_k += sgn * hp['sym_H']
                eq_y_k += sgn * hp['sym_V']
                eq_m_k += sgn * (hp['sym_V'] * (xc - xpk) - hp['sym_H'] * (yc - ypk))

            if con_carichi_reali:
                for p in dati_json.get('point_loads', []):
                    nid_p = int(p['node_id'])
                    if node_owner_tr[nid_p] != k_:
                        continue
                    xp_, yp_ = nodi[nid_p]
                    fx_, fy_, mm_ = p.get('Fx', 0), p.get('Fy', 0), p.get('M', 0)
                    eq_x_k += fx_; eq_y_k += fy_
                    eq_m_k += mm_ + fy_ * (xp_ - xpk) - fx_ * (yp_ - ypk)
                for el_load in dati_json.get('element_loads', []):
                    c = _carico(el_load)
                    if c is None: continue
                    eid_c = int(c[0]['id'])
                    if elem_to_tr[eid_c] != k_:
                        continue
                    Rx, Ry, X_bar, Y_bar = _risultante(c)
                    eq_x_k += Rx; eq_y_k += Ry
                    eq_m_k += Ry * (X_bar - xpk) - Rx * (Y_bar - ypk)

            if carico_X1 and node_owner_tr.get(nodo_svinc) == k_:
                xs_v, ys_v = nodi[nodo_svinc]
                eq_x_k += X1_fx; eq_y_k += X1_fy
                eq_m_k += (X1_m if tipo_svinc == 'M' else 0) + X1_fy * (xs_v - xpk) - X1_fx * (ys_v - ypk)

            eqs_all.extend([eq_x_k, eq_y_k, eq_m_k])
            blocchi_tr.append((tr, eq_x_k, eq_y_k, eq_m_k))
        return eqs_all, blocchi_tr

    # --- Risoluzione S0 ed S1 e verifica labilità sistema principale ---
    eqs_tronchi_S0, blocchi_S0 = _costruisci_sistema_tronchi(reazioni_simboliche, con_carichi_reali=True, carico_X1=False)
    vars_tutte = variabili + vars_cerniere
    try:
        A_chk, _ = sp.linear_eq_to_matrix(eqs_tronchi_S0, vars_tutte)
        if A_chk.rank() < len(vars_tutte):
            report.append(_box_errore(
                "Incognita iperstatica non valida",
                "Rimuovendo questo vincolo il <b>sistema principale risulta labile</b> "
                "(le reazioni restanti non sono determinabili dall'equilibrio). Scegli un'altra incognita iperstatica."
            ))
            return "\n".join(report), "", ""
    except Exception:
        pass
    soluzione_S0 = risolvi_sicuro(eqs_tronchi_S0, vars_tutte, "S0")

    soluzione_S1 = {}
    blocchi_S1 = []
    if grado_iperstaticita > 0:
        eqs_tronchi_S1, blocchi_S1 = _costruisci_sistema_tronchi(reazioni_simboliche, con_carichi_reali=False, carico_X1=True)
        soluzione_S1 = risolvi_sicuro(eqs_tronchi_S1, vars_tutte, "S1")

    # ==========================================
    # CALCOLO SOLLECITAZIONI M, T, N PER S0 E S1
    # ==========================================
    Forze_S0 = {n: {'Fx': 0, 'Fy': 0, 'M': 0} for n in nodi}
    Forze_S1 = {n: {'Fx': 0, 'Fy': 0, 'M': 0} for n in nodi}
    for p in dati_json.get('point_loads', []):
        nid = int(p.get('node_id'))
        Forze_S0[nid]['Fx'] += p.get('Fx', 0); Forze_S0[nid]['Fy'] += p.get('Fy', 0); Forze_S0[nid]['M'] += p.get('M', 0)
    if grado_iperstaticita > 0:
        Forze_S1[nodo_svinc]['Fx'] += X1_fx
        Forze_S1[nodo_svinc]['Fy'] += X1_fy
        Forze_S1[nodo_svinc]['M'] += X1_m
    for r in reazioni_simboliche:
        nid = r['nodo']; tipo_r = r['tipo']
        val_S0 = soluzione_S0.get(r['sym'], 0)
        if tipo_r == 'Fx': Forze_S0[nid]['Fx'] += val_S0
        elif tipo_r == 'Fy': Forze_S0[nid]['Fy'] += val_S0
        elif tipo_r == 'M': Forze_S0[nid]['M'] += val_S0
        elif tipo_r == 'Fn':
            ang_rad = cs_ang(r['angle'] + 90)
            Forze_S0[nid]['Fx'] += val_S0 * ang_rad[0]
            Forze_S0[nid]['Fy'] += val_S0 * ang_rad[1]
        if grado_iperstaticita > 0:
            val_S1 = soluzione_S1.get(r['sym'], 0)
            if tipo_r == 'Fx': Forze_S1[nid]['Fx'] += val_S1
            elif tipo_r == 'Fy': Forze_S1[nid]['Fy'] += val_S1
            elif tipo_r == 'M': Forze_S1[nid]['M'] += val_S1
            elif tipo_r == 'Fn':
                ang_rad = cs_ang(r['angle'] + 90)
                Forze_S1[nid]['Fx'] += val_S1 * ang_rad[0]
                Forze_S1[nid]['Fy'] += val_S1 * ang_rad[1]

    def calcola_sollecitazioni(sistema_forze, el, is_S0=True):
        n1, n2 = int(el['n1']), int(el['n2'])
        x1, y1 = nodi[n1]; x2, y2 = nodi[n2]
        L = hyp(x2 - x1, y2 - y1)
        s = sp.Symbol("s")
        if L == 0:
            return sp.S(0), sp.S(0), sp.S(0)
        cos_alpha = (x2 - x1) / L
        sin_alpha = (y2 - y1) / L
        visited = set(); stack = [n1]
        while stack:
            curr = stack.pop()
            if curr not in visited:
                visited.add(curr)
                for neighbor in adj[curr]:
                    if neighbor != n2 and neighbor not in visited:
                        stack.append(neighbor)
        xc = x1 + cos_alpha * s; yc = y1 + sin_alpha * s
        M_s, Fx_sum, Fy_sum = 0, 0, 0
        for v in visited:
            fx, fy, m = sistema_forze[v]['Fx'], sistema_forze[v]['Fy'], sistema_forze[v]['M']
            Fx_sum += fx; Fy_sum += fy
            M_s += m + fy * (nodi[v][0] - xc) - fx * (nodi[v][1] - yc)
        if is_S0:
            t_ = sp.Symbol('t')
            for el_load in dati_json.get('element_loads', []):
                c = _carico(el_load)
                if c is None: continue
                target_el, L_t, ux_t, uy_t, wx, wy, tipo_q = c
                if target_el['id'] == el['id']:
                    k_ = ux_t * wy - uy_t * wx
                    phi = _phi(tipo_q, t_, L_t)
                    Phi = sp.integrate(phi, (t_, 0, s))
                    M_s += k_ * sp.integrate((t_ - s) * phi, (t_, 0, s))
                    Fx_sum += wx * Phi; Fy_sum += wy * Phi
                elif int(target_el['n1']) in visited and int(target_el['n2']) in visited:
                    Rx, Ry, X_c, Y_c = _risultante(c)
                    Fx_sum += Rx; Fy_sum += Ry
                    M_s += Ry * (X_c - xc) - Rx * (Y_c - yc)
        T_s = -sp.diff(M_s, s)
        N_s = -Fx_sum * cos_alpha - Fy_sum * sin_alpha
        return pulisci_espr(sp.simplify(M_s)), pulisci_espr(sp.simplify(T_s)), pulisci_espr(sp.simplify(N_s))

    M0_dict, M1_dict, T0_dict, N0_dict, T1_dict, N1_dict = {}, {}, {}, {}, {}, {}
    for el in elementi:
        M0, T0, N0 = calcola_sollecitazioni(Forze_S0, el, True)
        if grado_iperstaticita > 0:
            M1, T1, N1 = calcola_sollecitazioni(Forze_S1, el, False)
        else:
            M1, T1, N1 = sp.S(0), sp.S(0), sp.S(0)
        M0_dict[el['id']] = M0; M1_dict[el['id']] = M1
        T0_dict[el['id']] = T0; N0_dict[el['id']] = N0
        T1_dict[el['id']] = T1; N1_dict[el['id']] = N1

    # --- Helper condiviso per Sezione 3 (S0) e Sezione 4 (S1) ---
    nodi_f = {k: (float(v[0]), float(v[1])) for k, v in nodi.items()}
    s_sym_tab = sp.Symbol("s")

    def _tabella_vincoli_esterni(sol_dict, con_unita, include_X1=False):
        righe_sup = {}
        for i_s, s_ in enumerate(supports_n):
            nid_s = int(s_['node_id'])
            note_extra = []
            h_v, v_v, m_v = None, None, None
            # Identifica componenti su questo vincolo
            for r in reazioni_simboliche:
                if r.get('sup_idx') != i_s:
                    continue
                val_r = float(sp.N(sol_dict.get(r['sym'], 0)))
                if r['tipo'] == 'Fx': h_v = (h_v or 0.0) + val_r
                elif r['tipo'] == 'Fy': v_v = (v_v or 0.0) + val_r
                elif r['tipo'] == 'M': m_v = (m_v or 0.0) + val_r
                elif r['tipo'] == 'Fn':
                    ca, sa = cs_ang(r['angle'] + 90)
                    h_v = (h_v or 0.0) + val_r * float(ca)
                    v_v = (v_v or 0.0) + val_r * float(sa)
                    note_extra.append(f"{sym_reaz_html(r['sym'])} = {fmt_num(val_r, 2)}{' N' if con_unita else ''}")
            if include_X1 and grado_iperstaticita > 0:
                if (tipo_svinc == 'Fn' and idx_X1_sup == i_s) or (tipo_svinc in ('Fx', 'Fy', 'M') and nid_s == nodo_svinc and s_.get({'Fx':'fix_u','Fy':'fix_v','M':'fix_w'}[tipo_svinc])):
                    if tipo_svinc == 'Fx': h_v = (h_v or 0.0) + 1.0
                    elif tipo_svinc == 'Fy': v_v = (v_v or 0.0) + 1.0
                    elif tipo_svinc == 'M': m_v = (m_v or 0.0) + 1.0
                    elif tipo_svinc == 'Fn':
                        ca, sa = cs_ang(ang_X1 + 90)
                        h_v = (h_v or 0.0) + float(ca); v_v = (v_v or 0.0) + float(sa)
                    note_extra.append("X<sub>1</sub> = 1")
            elif grado_iperstaticita > 0 and not include_X1:
                if (tipo_svinc == 'Fn' and idx_X1_sup == i_s) or (tipo_svinc in ('Fx', 'Fy', 'M') and nid_s == nodo_svinc and s_.get({'Fx':'fix_u','Fy':'fix_v','M':'fix_w'}[tipo_svinc])):
                    note_extra.append("incognita X<sub>1</sub> rimossa")
            lbl_vin = f"Nodo {nid_s}" + (f" ({', '.join(note_extra)})" if note_extra else "")
            righe_sup[i_s] = (lbl_vin, h_v, v_v, m_v)

        u_f = " [N]" if con_unita else ""
        u_m = " [N&middot;m]" if con_unita else ""
        rows_sup = []
        for i_s in sorted(righe_sup):
            lbl, hv, vv, mv = righe_sup[i_s]
            cH = f"<b>{fmt_num(hv, 2)}</b>" if hv is not None else "&mdash;"
            cV = f"<b>{fmt_num(vv, 2)}</b>" if vv is not None else "&mdash;"
            cM = f"<b>{fmt_num(mv, 2)}</b>" if mv is not None else "&mdash;"
            rows_sup.append(f"<tr><td>{lbl}</td><td>{cH}</td><td>{cV}</td><td>{cM}</td></tr>")
        return _tabella_html(["Vincolo", f"H{u_f}", f"V{u_f}", f"M{u_m}"], rows_sup)

    def _append_sezione_sistema(tag, tit_sec, reaz_plot, sol_dict, sol_plot, blocchi_tr, Md, Td, Nd, con_unita, accent_col):
        report.append(f"<h2 style='page-break-before: always;'>{tit_sec}</h2>")

        # 1. 3 schemi grafici globali (dimensionati per stare tutti insieme sotto il titolo della sezione)
        p_versi = f"{tag.lower()}_versi.png"
        p_calc = f"{tag.lower()}_calcolate.png"
        p_mono = f"{tag.lower()}_monoassiali.png"
        i_versi = plot_sistema_S0(nodi, elementi, reaz_plot, sol_plot, p_versi, out_dir, "versi", con_unita=False, cerniere_plot=cerniere_plot)
        i_calc = plot_sistema_S0(nodi, elementi, reaz_plot, sol_plot, p_calc, out_dir, "calcolate", con_unita=con_unita, cerniere_plot=cerniere_plot)
        i_mono = plot_sistema_S0(nodi, elementi, reaz_plot, sol_plot, p_mono, out_dir, "monoassiali", con_unita=con_unita, cerniere_plot=cerniere_plot)

        sub_tag = "S<sub>0</sub>" if tag == "S0" else "S<sub>1</sub>"
        report.append(_img_block_html(f"<h4 style='margin:0;'>{sub_tag} &mdash; reazioni (scelta dei versi)</h4>", i_versi, 540, 290))
        report.append(_img_block_html(f"<h4 style='margin:0;'>{sub_tag} &mdash; reazioni calcolate</h4>", i_calc, 540, 290))
        report.append(_img_block_html(f"<h4 style='margin:0;'>{sub_tag} &mdash; riferimenti monoassiali (ascissa s)</h4>", i_mono, 540, 290))

        # 2. Tabella riassuntiva dei vincoli esterni
        report.append(f"<h3>Reazioni vincolari esterne ({sub_tag})</h3>")
        report.append(_tabella_vincoli_esterni(sol_dict, con_unita=con_unita, include_X1=(tag == "S1")))

        # 3. Equilibrio esploso tronco per tronco + "Da cui:" su due colonne
        report.append(f"<h3>Equilibrio esploso tronco per tronco ({sub_tag})</h3>")
        report.append(f"<div style='background-color: #f8f9fa; padding: 14px; border-left: 4px solid {accent_col}; margin-bottom: 16px; border-radius: 4px;'>")
        for tr, ex_k, ey_k, em_k in blocchi_tr:
            pk = tr['polo']
            report.append(f"<div style='page-break-inside: avoid; margin-bottom: 6px;'>")
            report.append(f"<p style='margin: 8px 0 4px 0; color: {accent_col};'><b>TRATTO {tr['idx']}</b> &mdash; tronco ({tr['desc']}), polo dei momenti in <b>Nodo {pk}</b>:</p>")
            report.append(_blocco_equazioni_html([
                ("&Sigma;F<sub>x</sub> = 0", ex_k),
                ("&Sigma;F<sub>y</sub> = 0", ey_k),
                (f"&Sigma;M<sub>{pk}</sub> = 0", em_k),
            ]))
            report.append("</div>")

        voci_est = []
        for r in reazioni_simboliche:
            val_r = float(sp.N(sol_dict.get(r['sym'], 0)))
            un_r = (" N&middot;m" if r['tipo'] == 'M' else " N") if con_unita else ""
            voci_est.append(f"<b>{sym_reaz_html(r['sym'])}</b> = {fmt_num(val_r, 2)}{un_r}")
        voci_int = []
        for hp in hinge_pins:
            for sk in (hp['sym_H'], hp['sym_V']):
                val_c = float(sp.N(sol_dict.get(sk, 0)))
                un_c = " N" if con_unita else ""
                voci_int.append(f"<b>{sym_reaz_html(sk)}</b> = {fmt_num(val_c, 2)}{un_c} <small>(su T{hp['t_A']+1}; opposta su T{hp['t_B']+1})</small>")
        if not voci_int:
            mid_e = (len(voci_est) + 1) // 2
            col1, col2 = voci_est[:mid_e], voci_est[mid_e:]
        else:
            col1, col2 = voci_est, voci_int
        report.append("<div style='page-break-inside: avoid;'>")
        report.append("<p style='margin: 10px 0 4px 0;'><b>Da cui:</b></p>")
        report.append(
            "<table class='eq-table' style='width:100%; border:none; background:transparent;'><tr>"
            f"<td style='width:50%; vertical-align:top; text-align:left; border:none; padding:2px 10px 2px 0;'>{'<br>'.join(col1)}</td>"
            f"<td style='width:50%; vertical-align:top; text-align:left; border:none; padding:2px 0 2px 10px;'>{'<br>'.join(col2)}</td>"
            "</tr></table></div></div>"
        )

        # 4. Tabella Azioni di estremita' delle aste
        u_f = " [N]" if con_unita else ""
        u_m = " [N&middot;m]" if con_unita else ""
        idx_s = "0" if tag == "S0" else "1"
        sub_p = f"<sub>{idx_s}</sub>"
        report.append(f"<h3>Azioni di estremit&agrave; delle aste ({sub_tag})</h3>")
        rows_est = []
        for el in elementi:
            eid = el['id']
            L_ = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
            n_i = float(sp.N(Nd[eid].subs(s_sym_tab, 0))); n_j = float(sp.N(Nd[eid].subs(s_sym_tab, L_)))
            t_i = float(sp.N(Td[eid].subs(s_sym_tab, 0))); t_j = float(sp.N(Td[eid].subs(s_sym_tab, L_)))
            m_i = float(sp.N(Md[eid].subs(s_sym_tab, 0))); m_j = float(sp.N(Md[eid].subs(s_sym_tab, L_)))
            rows_est.append(
                f"<tr><td><b>Asta {eid}</b> <small>({el['n1']}&rarr;{el['n2']})</small></td>"
                f"<td>{fmt_num(n_i, 2)}</td><td>{fmt_num(t_i, 2)}</td><td>{fmt_num(m_i, 2)}</td>"
                f"<td>{fmt_num(n_j, 2)}</td><td>{fmt_num(t_j, 2)}</td><td>{fmt_num(m_j, 2)}</td></tr>"
            )
        report.append(_tabella_html([
            "Asta",
            f"N<sub>{idx_s},i</sub>{u_f}", f"T<sub>{idx_s},i</sub>{u_f}", f"M<sub>{idx_s},i</sub>{u_m}",
            f"N<sub>{idx_s},j</sub>{u_f}", f"T<sub>{idx_s},j</sub>{u_f}", f"M<sub>{idx_s},j</sub>{u_m}",
        ], rows_est, keep=True, compact=True))

        # 5. Tabella compatta unica per N, T, M (evita pagine quasi vuote con ultima asta isolata)
        report.append(f"<h3>Caratteristiche della sollecitazione per tratto ({sub_tag})</h3>")
        rows_sol = []
        for el in elementi:
            eid = el['id']
            L_ = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
            Lf = float(L_)
            rows_sol.append(
                f"<tr><td><b>Asta {eid}</b> <small>({el['n1']}&rarr;{el['n2']})</small></td>"
                f"<td>0 &le; s &le; {Lf:.2f} m</td>"
                f"<td><span class='formula'>{pulisci_formula(Nd[eid], decimale=True)}</span></td>"
                f"<td><span class='formula'>{pulisci_formula(Td[eid], decimale=True)}</span></td>"
                f"<td><span class='formula'>{pulisci_formula(Md[eid], decimale=True)}</span></td></tr>"
            )
        report.append(_tabella_html([
            "Asta", "Dominio",
            f"N{sub_p}(s){u_f}", f"T{sub_p}(s){u_f}", f"M{sub_p}(s){u_m}"
        ], rows_sol, keep=True, compact=True))

        # 6. 3 Diagrammi globali del sistema
        sub_u = "₀" if tag == "S0" else "₁"
        tag_u = "S₀" if tag == "S0" else "S₁"
        um_mpl = " [N·m]" if con_unita else ""
        uf_mpl = " [N]" if con_unita else ""
        p_diag_sys = _salva_trio_diagrammi(
            nodi_f, elementi, Md, Td, Nd,
            f"Momento flettente M{sub_u} — {tag_u}{um_mpl} (lato fibre tese)",
            f"Taglio T{sub_u} — {tag_u}{uf_mpl}",
            f"Sforzo normale N{sub_u} — {tag_u}{uf_mpl}",
            f"plv_diag_{tag}.png", out_dir
        )
        # Normalizza etichette titoli nei grafici con notazione pulita
        report.append(_img_block_html(
            f"<h4 style='margin:0;'>Diagrammi globali ({sub_tag}): Momento M, Taglio T, Sforzo Normale N</h4>",
            p_diag_sys, 860, 940, page_break_before=True
        ))

    # ==========================================
    # SEZIONE 3: SISTEMA S0
    # ==========================================
    _append_sezione_sistema(
        "S0", "3. Sistema S<sub>0</sub> &mdash; carichi reali sulla struttura isostatica",
        reazioni_simboliche, soluzione_S0, soluzione_S0, blocchi_S0,
        M0_dict, T0_dict, N0_dict, con_unita=True, accent_col="#0d6efd"
    )

    # ==========================================
    # SEZIONE 4: SISTEMA S1
    # ==========================================
    X1_valore_numerico = 0.0
    X1_exact = sp.Integer(0)
    if grado_iperstaticita > 0:
        X1_sym = sp.Symbol('X_1')
        reazioni_simboliche_S1 = reazioni_simboliche.copy()
        reazioni_simboliche_S1.append({'nodo': nodo_svinc, 'tipo': tipo_svinc, 'sym': X1_sym, 'angle': ang_X1})
        soluzione_plot_S1 = soluzione_S1.copy()
        soluzione_plot_S1[X1_sym] = 1.0
        _append_sezione_sistema(
            "S1", "4. Sistema S<sub>1</sub> &mdash; incognita unitaria (X<sub>1</sub> = 1)",
            reazioni_simboliche_S1, soluzione_S1, soluzione_plot_S1, blocchi_S1,
            M1_dict, T1_dict, N1_dict, con_unita=False, accent_col="#dc3545"
        )
    else:
        report.append("<h2 style='page-break-before: always;'>4. Sistema S<sub>1</sub> &mdash; incognita unitaria</h2><p>Non necessario (la struttura &egrave; isostatica).</p>")

    # ==========================================
    # SEZIONE 5 & 6: COEFFICIENTI PLV (TABELLA COMPATTA) + CONGRUENZA CON CEDIMENTI
    # ==========================================
    if grado_iperstaticita > 0:
        report.append("<h2 style='page-break-before: always;'>5. Coefficienti di deformabilit&agrave; &mdash; Impostazione degli integrali (PLV)</h2>")
        report.append("<div style='background-color: #e2e3e5; padding: 10px 14px; border-left: 4px solid #6c757d; margin-bottom: 14px;'>")
        report.append("<b>Formule generali di Mohr (solo flessione):</b><br>")
        report.append("<span class='formula'>&delta;<sub>10</sub> = (1/EI) &Sigma; &int;<sub>0</sub><sup>L</sup> M<sub>0</sub>(s)&middot;M<sub>1</sub>(s) ds</span>"
                      " &nbsp;&nbsp;e&nbsp;&nbsp; "
                      "<span class='formula'>&delta;<sub>11</sub> = (1/EI) &Sigma; &int;<sub>0</sub><sup>L</sup> M<sub>1</sub>(s)<sup>2</sup> ds</span></div>")

        _mat_d, _sec_d0 = dati_json.get('material', {}), dati_json.get('section', {})
        EI_val = _mat_d['EI'] if 'EI' in _mat_d else (_mat_d.get('E', 2.1e11) * (_sec_d0.get('b', 0.2) * _sec_d0.get('h', 0.4)**3) / 12)
        ei_noto_reale = (modalita == "verifica") or ('EI' in _mat_d)
        delta_10, delta_11 = sp.Integer(0), sp.Integer(0)

        rows_mohr = []
        for el in elementi:
            x_sym = sp.Symbol("s")
            L = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
            M0_expr, M1_expr = M0_dict[el['id']], M1_dict[el['id']]
            d10_t = sp.integrate(M0_expr * M1_expr, (x_sym, 0, L))
            d11_t = sp.integrate(M1_expr * M1_expr, (x_sym, 0, L))
            delta_10 += d10_t; delta_11 += d11_t
            rows_mohr.append(
                f"<tr><td><b>Asta {el['id']}</b></td><td>{float(L):.2f}</td>"
                f"<td><span class='formula'>{pulisci_formula(M0_expr, decimale=True)}</span></td>"
                f"<td><span class='formula'>{pulisci_formula(M1_expr, decimale=True)}</span></td>"
                f"<td>{fmt_smart(float(sp.N(d10_t)), min_sig=4)}</td>"
                f"<td>{fmt_smart(float(sp.N(d11_t)), min_sig=4)}</td></tr>"
            )
        delta_10 = sp.simplify(delta_10); delta_11 = sp.simplify(delta_11)
        d10_num = float(sp.N(delta_10)); d11_num = float(sp.N(delta_11))
        d10_mohr = fmt_smart(d10_num, min_sig=4)
        d11_mohr = fmt_smart(d11_num, min_sig=4)
        rows_mohr.append(
            f"<tr style='background-color:#eef4ff; font-weight:bold;'>"
            f"<td colspan='4' style='text-align:right;'>Totale:</td>"
            f"<td>{d10_mohr}</td>"
            f"<td>{d11_mohr}</td></tr>"
        )
        report.append(_tabella_html([
            "Asta", "L [m]", "M<sub>0</sub>(s) [N&middot;m]", "M<sub>1</sub>(s)",
            "&int;<sub>0</sub><sup>L</sup> M<sub>0</sub>&middot;M<sub>1</sub> ds",
            "&int;<sub>0</sub><sup>L</sup> M<sub>1</sub><sup>2</sup> ds"
        ], rows_mohr, keep=True))

        if delta_11 == 0 or abs(float(delta_11)) < 1e-14 * max(1.0, float(sum(L_**3 for L_ in [hyp(nodi[int(e_['n2'])][0]-nodi[int(e_['n1'])][0], nodi[int(e_['n2'])][1]-nodi[int(e_['n1'])][1]) for e_ in elementi]))):
            stop("Calcolo impossibile: &delta;<sub>11</sub> = 0",
                 "L'incognita iperstatica scelta non produce momento flettente in S<sub>1</sub> (M<sub>1</sub> &equiv; 0). Scegli un'altra incognita.")

        un_d10 = "rad" if tipo_svinc in ('M', 'M_int') else "m"
        un_d11 = "rad/(N&middot;m)" if tipo_svinc in ('M', 'M_int') else "m/N"
        extra_d10 = f" = <b>{fmt_smart(d10_num / float(EI_val), min_sig=4)} {un_d10}</b>" if ei_noto_reale else ""
        extra_d11 = f" = <b>{fmt_smart(d11_num / float(EI_val), min_sig=4)} {un_d11}</b>" if ei_noto_reale else ""
        report.append("<div style='padding: 12px 15px; border: 1px dashed #6c757d; margin-bottom: 18px; background-color: #ffffff;'>")
        report.append(f"<b>&delta;<sub>10</sub></b> = <span class='formula'>{d10_mohr.replace('-', '&minus;')} / EI</span>{extra_d10}<br>")
        report.append(f"<b>&delta;<sub>11</sub></b> = <span class='formula'>{d11_mohr.replace('-', '&minus;')} / EI</span>{extra_d11}</div>")

        # --- DISTORSIONI (cedimenti vincolari e termiche) ---
        alpha_mat = ex(dati_json.get('material', {}).get('alpha', 0))
        h_sez_T = ex(dati_json.get('section', {}).get('h', 0.4))
        s_sym_d = sp.Symbol("s")
        termine_T, termine_N = sp.Integer(0), sp.Integer(0)
        KT_sum = sp.Integer(0)
        righe_term = []
        for tl in dati_json.get('thermal_loads', []):
            eid = int(tl['elem_id'])
            el_t = next((e for e in elementi if int(e['id']) == eid), None)
            if el_t is None: continue
            L_t = hyp(nodi[int(el_t['n2'])][0] - nodi[int(el_t['n1'])][0], nodi[int(el_t['n2'])][1] - nodi[int(el_t['n1'])][1])
            alpha_t = ex(tl.get('alpha', alpha_mat))
            T_pos, T_neg = ex(tl.get('T_pos', 0)), ex(tl.get('T_neg', 0))
            kappa_T = alpha_t * (T_pos - T_neg) / h_sez_T
            eps_T = alpha_t * (T_pos + T_neg) / 2.0
            int_M1 = sp.integrate(M1_dict[el_t['id']], (s_sym_d, 0, L_t))
            int_N1 = sp.integrate(N1_dict[el_t['id']], (s_sym_d, 0, L_t))
            termine_T += kappa_T * int_M1
            KT_sum += alpha_t * (T_pos - T_neg) * int_M1
            termine_N += eps_T * int_N1
            righe_term.append(f"Asta {eid}: &kappa;<sub>T</sub> = {fmt_sci_html(float(kappa_T), sig=4)} m<sup>&minus;1</sup> &rarr; &int;M<sub>1</sub>&kappa;<sub>T</sub>ds = {fmt_smart(float(kappa_T*int_M1), min_sig=4)}"
                              + (f"; &epsilon;<sub>T</sub> = {fmt_sci_html(float(eps_T), sig=4)} &rarr; &int;N<sub>1</sub>&epsilon;<sub>T</sub>ds = {fmt_smart(float(eps_T*int_N1), min_sig=4)}" if abs(eps_T) > 0 else ""))

        sup_map = {}
        for s_ in dati_json.get('supports', []):
            sup_map.setdefault(int(s_['node_id']), {}).update(s_)
        chiavi_delta = {'Fx': 'delta_u', 'Fy': 'delta_v', 'M': 'delta_w', 'Fn': 'delta_n'}
        somma_Rd = sp.Integer(0)
        righe_ced = []
        for r_ in reazioni_simboliche:
            if r_['tipo'] == 'Fn' and r_.get('sup_idx') is not None:
                dval = ex(supports_n[r_['sup_idx']].get('delta_n', 0))
            else:
                dval = ex(sup_map.get(r_['nodo'], {}).get(chiavi_delta[r_['tipo']], 0))
            if dval != 0:
                R1v = ex(soluzione_S1.get(r_['sym'], 0))
                somma_Rd += R1v * dval
                u_d = "rad" if r_['tipo'] == 'M' else "m"
                righe_ced.append(f"Nodo {r_['nodo']} ({sym_reaz_html(r_['sym'])}): R<sub>1</sub> = {fmt_smart(float(R1v))}, &delta;<sub>c</sub> = {fmt_smart(float(dval), fixed_min=1e-4)} {u_d} &rarr; R<sub>1</sub>&middot;&delta;<sub>c</sub> = {fmt_smart(float(R1v*dval), fixed_min=1e-4)} {un_d10}")
        eta_X1 = sp.Integer(0)
        if tipo_svinc == 'Fn' and idx_X1_sup is not None:
            eta_X1 = ex(supports_n[idx_X1_sup].get('delta_n', 0))
        elif tipo_svinc in chiavi_delta:
            eta_X1 = ex(sup_map.get(nodo_svinc, {}).get(chiavi_delta[tipo_svinc], 0))
        if eta_X1 != 0:
            u_x = "rad" if tipo_svinc in ('M', 'M_int') else "m"
            righe_ced.append(f"Vincolo degradato X<sub>1</sub> (Nodo {nodo_svinc}): R<sub>1</sub> = 1, &delta;<sub>c</sub> = {fmt_smart(float(eta_X1), fixed_min=1e-4)} {u_x} &rarr; 1&middot;&delta;<sub>c</sub> = {fmt_smart(float(eta_X1), fixed_min=1e-4)} {un_d10}")

        L_ev_exact = eta_X1 + somma_Rd
        ha_distorsioni = bool(righe_term or righe_ced)
        dip_ced = bool(righe_ced)
        dip_term = bool(righe_term)
        termine_noto_dist = EI_val * (L_ev_exact - termine_T - termine_N)

        report.append("<h2>6. Equazione di congruenza e calcolo di X<sub>1</sub></h2>")
        report.append("<div style='background-color: #f8f9fa; padding: 12px 15px; border-left: 4px solid #dc3545; margin-bottom: 14px;'>")
        if righe_ced:
            report.append("<b>Lavoro esterno virtuale per cedimenti vincolari (L<sub>ev</sub> = &Sigma; R<sub>1,k</sub>&middot;&delta;<sub>c,k</sub>):</b><ul style='margin:4px 0 8px 0;'>")
            for r_txt in righe_ced:
                report.append(f"<li>{r_txt}</li>")
            report.append(f"</ul><b>L<sub>ev</sub> = {fmt_smart(float(L_ev_exact), fixed_min=1e-4)} {un_d10}</b><br><br>")
        if righe_term:
            report.append("<b>Contributo variazioni termiche:</b><ul style='margin:4px 0 8px 0;'>")
            for r_txt in righe_term:
                report.append(f"<li>{r_txt}</li>")
            report.append("</ul>")

        X1_sym_sol = sp.Symbol('X1')
        try:
            soluzione_X1 = sp.solve(sp.Eq(delta_10 + X1_sym_sol * delta_11, termine_noto_dist), X1_sym_sol)
        except Exception as e:
            stop("Errore nell'equazione di congruenza", f"{e}")
        if not soluzione_X1:
            stop("Equazione di congruenza senza soluzione", "&delta;<sub>10</sub> + X<sub>1</sub>&delta;<sub>11</sub> = termine noto non ammette soluzione.")
        X1_exact = pulisci_espr(sp.simplify(soluzione_X1[0]))
        X1_valore_numerico = float(X1_exact)
        unita_X1 = "N&middot;m" if tipo_svinc in ('M', 'M_int') else "N"

        d10_eq = d10_mohr.replace('-', '&minus;')
        d11_eq = d11_mohr.replace('-', '&minus;')
        x1_str = fmt_smart(X1_valore_numerico)
        x1_eq = x1_str.replace('-', '&minus;')

        if righe_term:
            tnoto_num = float(sp.N(termine_noto_dist))
            tnoto_eq = fmt_smart(tnoto_num, min_sig=4).replace('-', '&minus;')
            report.append("<b>Equazione di congruenza:</b><br>"
                          "<span class='formula'>(&delta;<sub>10</sub> + X<sub>1</sub>&middot;&delta;<sub>11</sub>)/EI + &int;M<sub>1</sub>&kappa;<sub>T</sub>ds + &int;N<sub>1</sub>&epsilon;<sub>T</sub>ds = L<sub>ev</sub></span><br>"
                          f"<span class='formula'>X<sub>1</sub> = [EI&middot;(L<sub>ev</sub> &minus; L<sub>th</sub>) &minus; ({d10_eq})] / {d11_eq} = ({tnoto_eq} &minus; ({d10_eq})) / {d11_eq} = {x1_eq} {unita_X1}</span>")
        elif righe_ced:
            lev_ei = float(sp.N(EI_val * L_ev_exact))
            lev_ei_eq = fmt_smart(lev_ei, min_sig=4).replace('-', '&minus;')
            report.append("<b>Equazione di congruenza:</b><br>"
                          "<span class='formula'>(&delta;<sub>10</sub> + X<sub>1</sub>&middot;&delta;<sub>11</sub>)/EI = L<sub>ev</sub> &rArr; &delta;<sub>10</sub> + X<sub>1</sub>&middot;&delta;<sub>11</sub> = EI&middot;L<sub>ev</sub></span><br>"
                          f"<span class='formula'>X<sub>1</sub> = ({lev_ei_eq} &minus; ({d10_eq})) / {d11_eq} = {x1_eq} {unita_X1}</span>")
        else:
            report.append("<b>Equazione di congruenza (vincoli rigidi, nessun cedimento):</b><br>"
                          "<span class='formula'>&delta;<sub>10</sub> + X<sub>1</sub> &middot; &delta;<sub>11</sub> = 0 &rArr; X<sub>1</sub> = &minus;&delta;<sub>10</sub> / &delta;<sub>11</sub></span><br>"
                          f"<span class='formula'>X<sub>1</sub> = &minus;({d10_eq}) / {d11_eq} = {x1_eq} {unita_X1}</span>")
        report.append("</div>")
        report.append(
            "<div style='background-color:#d1e7dd; color:#0f5132; padding:12px 15px; border-left: 5px solid #198754; border-radius: 4px; margin-bottom: 18px;'>"
            f"<span style='font-size:18px;'><b>X<sub>1</sub> = {x1_str} {unita_X1}</b></span></div>"
        )
    else:
        unita_X1 = "N"
        report.append("<h2 style='page-break-before: always;'>5. Coefficienti di deformabilit&agrave; (PLV)</h2><p>Non necessari (la struttura &egrave; isostatica).</p>")
        report.append("<h2>6. Equazione di congruenza</h2><p>Non necessaria (la struttura &egrave; isostatica).</p>")

    # ==========================================
    # SEZIONE 7: SISTEMA REALE — SOVRAPPOSIZIONE S = S0 + X1 * S1
    # ==========================================
    report.append("<h2 style='page-break-before: always;'>7. Sistema reale &mdash; sovrapposizione S = S<sub>0</sub> + X<sub>1</sub> &middot; S<sub>1</sub></h2>")

    nome_sym = {'Fx': 'H', 'Fy': 'V', 'M': 'M', 'Fn': 'R'}
    def _nome_reaz(tipo_, nodo_, sym_=None):
        parti_ = str(sym_).split('_') if sym_ is not None else []
        if tipo_ == 'Fn':
            suff = f"n{nodo_}.{parti_[2]}" if len(parti_) == 3 else f"n{nodo_}"
            return f"R<sub>{suff}</sub>"
        return f"{nome_sym[tipo_]}<sub>{nodo_}</sub>"

    reaz_fin = []
    soluzione_reale_plot = {}
    for r in reazioni_simboliche:
        v = soluzione_S0.get(r['sym'], 0)
        if grado_iperstaticita > 0:
            v = v + X1_exact * soluzione_S1.get(r['sym'], 0)
        v_ex = pulisci_espr(sp.simplify(v))
        soluzione_reale_plot[r['sym']] = float(sp.N(v_ex))
        reaz_fin.append({'nodo': r['nodo'], 'tipo': r['tipo'], 'angle': r.get('angle', 0),
                         'sup_idx': r.get('sup_idx'), 'pendolo': r.get('pendolo', False),
                         'nome': _nome_reaz(r['tipo'], r['nodo'], r['sym'] if r['tipo'] == 'Fn' else None),
                         'val': v_ex})
    reazioni_reale_plot = list(reazioni_simboliche)
    if grado_iperstaticita > 0 and tipo_svinc in nome_sym:
        sym_x_nome = {'Fx': f'H_{nodo_svinc}', 'Fy': f'V_{nodo_svinc}', 'M': f'M_{nodo_svinc}'}.get(tipo_svinc) or str(sym_Fn.get(idx_X1_sup, 'X_1'))
        sym_x_r = sp.Symbol(sym_x_nome)
        soluzione_reale_plot[sym_x_r] = float(X1_exact)
        reazioni_reale_plot.append({'nodo': nodo_svinc, 'tipo': tipo_svinc, 'sym': sym_x_r, 'angle': ang_X1,
                                    'sup_idx': idx_X1_sup,
                                    'pendolo': bool(tipo_svinc == 'Fn' and idx_X1_sup is not None and supports_n[idx_X1_sup].get('pendolo'))})
        reaz_fin.append({'nodo': nodo_svinc, 'tipo': tipo_svinc, 'angle': ex(ang_X1),
                         'sup_idx': idx_X1_sup,
                         'pendolo': bool(tipo_svinc == 'Fn' and supports_n[idx_X1_sup].get('pendolo')),
                         'nome': (_nome_reaz(tipo_svinc, nodo_svinc, sym_Fn.get(idx_X1_sup) if tipo_svinc == 'Fn' else None) + " (= X<sub>1</sub>)"),
                         'val': X1_exact})

    for hp in hinge_pins:
        for sk in (hp['sym_H'], hp['sym_V']):
            v_hp = soluzione_S0.get(sk, 0) + (X1_exact * soluzione_S1.get(sk, 0) if grado_iperstaticita > 0 else 0)
            soluzione_reale_plot[sk] = float(sp.N(v_hp))

    ha_carrello_incl = any(
        rf['tipo'] == 'Fn' and not rf.get('pendolo')
        and not (rf.get('sup_idx') is not None and supports_n[rf['sup_idx']].get('bipendolo'))
        for rf in reaz_fin
    )
    ha_pendoli = any(s_.get('pendolo') or s_.get('bipendolo') for s_ in supports_n)
    conv_parti = ["H positiva verso +x, V positiva verso +y, M positivo antiorario"]
    if ha_carrello_incl:
        conv_parti.append("R<sub>n</sub> positiva lungo la normale del carrello inclinato (angolo + 90&deg;)")
    if ha_pendoli:
        conv_parti.append("per pendoli e doppi pendoli inclinati R<sub>n</sub> &egrave; lungo l'asse delle bielle verso terra (R<sub>n</sub> &gt; 0 = trazione)")
    report.append(f"<p>Convenzione: {'; '.join(conv_parti)}.</p>")

    report.append("<h3>Reazioni vincolari definitive (esterne)</h3>")
    rows_rfin = []
    for rf in reaz_fin:
        unita = "N&middot;m" if rf['tipo'] == 'M' else "N"
        rows_rfin.append(f"<tr><td>{rf['nome']}</td><td>{rf['nodo']}</td><td><b>{fmt_num(float(rf['val']), 2)}</b></td><td>{unita}</td></tr>")
    report.append(_tabella_html(["Reazione", "Nodo", "Valore", "Unit&agrave;"], rows_rfin, keep=True))

    # Verifica equilibrio globale + forze trasmesse nelle cerniere
    def _componenti(rf):
        if rf['tipo'] == 'Fx': return rf['val'], sp.Integer(0), sp.Integer(0)
        if rf['tipo'] == 'Fy': return sp.Integer(0), rf['val'], sp.Integer(0)
        if rf['tipo'] == 'M':  return sp.Integer(0), sp.Integer(0), rf['val']
        c_, s_ = cs_ang(rf['angle'] + 90)
        return rf['val'] * c_, rf['val'] * s_, sp.Integer(0)

    res_x, res_y, res_m = Somma_Fx_ext, Somma_Fy_ext, Somma_M_ext
    for rf in reaz_fin:
        fx_, fy_, m_ = _componenti(rf)
        xn_, yn_ = nodi[rf['nodo']]
        res_x += fx_; res_y += fy_; res_m += m_ + fy_ * (xn_ - xp) - fx_ * (yn_ - yp)
    residui = [float(sp.N(res_x)), float(sp.N(res_y)), float(sp.N(res_m))]
    equilibrio_ok = all(abs(v) < 1e-6 for v in residui)

    righe_cern = []
    cerniere_ok = True
    for ci in cerniere_info:
        nc = ci['nodo']; xc, yc = nodi[nc]
        for id_asta_ramo, ramo_c in ci['rami']:
            Fx_b, Fy_b, M_hinge = sp.Integer(0), sp.Integer(0), sp.Integer(0)
            nodi_corpo = set(ramo_c) | ({nc} if ci['attach'] is not None and ci['attach'] == id_asta_ramo else set())
            for nid_ in nodi_corpo:
                f0, f1 = Forze_S0[nid_], Forze_S1[nid_]
                fx_n = f0['Fx'] + X1_exact * f1['Fx']
                fy_n = f0['Fy'] + X1_exact * f1['Fy']
                m_n = f0['M'] + X1_exact * f1['M']
                Fx_b += fx_n; Fy_b += fy_n
                M_hinge += m_n + fy_n * (nodi[nid_][0] - xc) - fx_n * (nodi[nid_][1] - yc)
            for el_load in dati_json.get('element_loads', []):
                c_ = _carico(el_load)
                if c_ is None: continue
                a_, b_ = int(c_[0]['n1']), int(c_[0]['n2'])
                if (a_ in ramo_c or a_ == nc) and (b_ in ramo_c or b_ == nc) and not (a_ == nc and b_ == nc):
                    fxq, fyq, xq, yq = _risultante(c_)
                    Fx_b += fxq; Fy_b += fyq
                    M_hinge += fyq * (xq - xc) - fxq * (yq - yc)
            Hx, Hy = pulisci_espr(sp.simplify(-Fx_b)), pulisci_espr(sp.simplify(-Fy_b))
            m_res = float(sp.N(M_hinge))
            if abs(m_res) > 1e-6: cerniere_ok = False
            righe_cern.append(
                f"<tr><td>{nc}</td><td>Ramo di Asta {id_asta_ramo} (nodi {', '.join(str(n_) for n_ in sorted(ramo_c))})</td>"
                f"<td><b>{fmt_num(float(Hx), 2)}</b></td><td><b>{fmt_num(float(Hy), 2)}</b></td><td>{fmt_num(m_res, 2)}</td></tr>"
            )
    if righe_cern:
        report.append("<h3>Forze interne definitive nelle cerniere</h3>")
        report.append(_tabella_html(["Nodo cerniera", "Ramo", "H [N]", "V [N]", "M in cerniera [N&middot;m]"], righe_cern, keep=True))

    # Espressioni analitiche definitive N(s), T(s), M(s)
    report.append("<h3>Caratteristiche definitive della sollecitazione per tratto</h3>")
    s_sym_8 = sp.Symbol("s")
    M_fin_all, T_fin_all, N_fin_all = {}, {}, {}
    for el in elementi:
        el_id = el['id']
        L = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
        Lf = float(L)
        M_fin = sp.expand(pulisci_espr(M0_dict[el_id] + X1_exact * M1_dict[el_id]))
        T_fin = sp.expand(pulisci_espr(T0_dict[el_id] + X1_exact * T1_dict[el_id]))
        N_fin = sp.expand(pulisci_espr(N0_dict[el_id] + X1_exact * N1_dict[el_id]))
        M_fin_all[el_id] = M_fin; T_fin_all[el_id] = T_fin; N_fin_all[el_id] = N_fin

        N_0 = float(sp.N(N_fin.subs(s_sym_8, 0))); N_L = float(sp.N(N_fin.subs(s_sym_8, L)))
        T_0 = float(sp.N(T_fin.subs(s_sym_8, 0))); T_L = float(sp.N(T_fin.subs(s_sym_8, L)))
        M_0 = float(sp.N(M_fin.subs(s_sym_8, 0))); M_L = float(sp.N(M_fin.subs(s_sym_8, L)))

        punti_staz_html = []
        cT = _coeff_poly(T_fin, s_sym_8)
        if cT is not None and len(np.trim_zeros(np.array(cT, dtype=float), 'f')) >= 2:
            for s_star in sorted(_radici_in(cT, Lf)):
                if 1e-5 * Lf < s_star < Lf * (1 - 1e-5):
                    M_star = float(sp.N(M_fin.subs(s_sym_8, s_star)))
                    punti_staz_html.append(f"s<sup>*</sup> = {fmt_num(s_star, 2)} m (T = 0) &rArr; <b>M<sub>extr</sub> = {fmt_num(M_star, 2)} N&middot;m</b>")

        report.append(f"<div class='keep-together' style='background-color: #f8f9fa; padding: 10px 12px; border: 1px solid #dee2e6; border-left: 4px solid #0d6efd; margin-bottom: 10px; border-radius: 4px; page-break-inside: avoid;'>")
        report.append(f"<h4 style='color: #0d6efd; margin: 0 0 6px 0;'>Asta {el_id} (Nodo {el['n1']} &rarr; Nodo {el['n2']}, L = {Lf:.2f} m)</h4>")
        report.append("<table class='keep-table' style='width: 100%; margin-bottom: 4px;'><thead><tr>"
                      "<th style='width: 16%;'>Grandezza</th>"
                      "<th style='width: 28%; background-color:#e3f2fd;'>N(s) [N]</th>"
                      "<th style='width: 28%; background-color:#e8f5e9;'>T(s) [N]</th>"
                      "<th style='width: 28%; background-color:#ffebee;'>M(s) [N&middot;m]</th></tr></thead><tbody>")
        report.append(f"<tr><td><b>Espressione</b></td>"
                      f"<td><span class='formula'>{pulisci_formula(N_fin, decimale=True)}</span></td>"
                      f"<td><span class='formula'>{pulisci_formula(T_fin, decimale=True)}</span></td>"
                      f"<td><span class='formula'>{pulisci_formula(M_fin, decimale=True)}</span></td></tr>")
        report.append(f"<tr><td><b>Estremi</b> <small>(s=0 / s={Lf:.2f})</small></td>"
                      f"<td>{fmt_num(N_0, 2)} / {fmt_num(N_L, 2)}</td>"
                      f"<td>{fmt_num(T_0, 2)} / {fmt_num(T_L, 2)}</td>"
                      f"<td>{fmt_num(M_0, 2)} / {fmt_num(M_L, 2)}</td></tr></tbody></table>")
        if punti_staz_html:
            report.append(f"<p style='margin: 4px 0 0 0; font-size: 13px; color: #842029;'><b>Punto stazionario:</b> {'; '.join(punti_staz_html)}</p>")
        report.append("</div>")

    # Schema globale Reale — riferimenti monoassiali + I 3 diagrammi finali non spezzati tra pagine
    img_reale_mono = plot_sistema_S0(
        nodi, elementi, reazioni_reale_plot, soluzione_reale_plot,
        "reale_monoassiali.png", out_dir, "monoassiali", con_unita=True, cerniere_plot=cerniere_plot
    )
    report.append(_img_block_html("<h4 style='margin:0;'>Reale &mdash; riferimenti monoassiali</h4>", img_reale_mono, 560, 300))

    img_fin_diag = _salva_trio_diagrammi(
        nodi_f, elementi, M_fin_all, T_fin_all, N_fin_all,
        "Momento flettente M — reale [N·m] (lato fibre tese)",
        "Taglio T — reale [N]",
        "Sforzo normale N — reale [N]",
        percorso_diag, out_dir
    )
    report.append(_img_block_html(
        "<h3 style='margin:0;'>Diagrammi finali delle sollecitazioni (Sistema Reale)</h3>",
        img_fin_diag, 860, 940, page_break_before=True
    ))

    if equilibrio_ok and cerniere_ok:
        print(f"[VERIFICA EQUILIBRIO GLOBALE] OK: residui Fx={residui[0]:.2e}, Fy={residui[1]:.2e}, M={residui[2]:.2e}")
    else:
        report.append("<div style='background-color:#f8d7da; color:#842029; padding:10px; border-left:4px solid #dc3545; margin-bottom:15px;'>"
                      f"&#10008; Verifica di equilibrio NON soddisfatta: residui &Sigma;F<sub>x</sub> = {fmt_sci_html(residui[0])}, "
                      f"&Sigma;F<sub>y</sub> = {fmt_sci_html(residui[1])}, &Sigma;M = {fmt_sci_html(residui[2])}.</div>")

    ULTIMI_RISULTATI['reazioni'] = [(rf['tipo'], rf['nodo'], float(rf['val'])) for rf in reaz_fin]
    ULTIMI_RISULTATI['X1'] = float(X1_exact)
    ULTIMI_RISULTATI['equilibrio_globale'] = max(abs(v) for v in residui)

    # Controllo coerenza nodale silenzioso (console + UI)
    Lmap = {e_['id']: hyp(nodi[int(e_['n2'])][0] - nodi[int(e_['n1'])][0], nodi[int(e_['n2'])][1] - nodi[int(e_['n1'])][1]) for e_ in elementi}
    hinged_per_nodo = {}
    for c_ in cerniere_interne:
        nc_ = int(c_['node_id'])
        inc_ = [e_ for e_ in elementi if nc_ in (int(e_['n1']), int(e_['n2']))]
        hs_ = {int(x_) for x_ in c_['elems']} if c_.get('elems') else ({int(inc_[0]['id'])} if inc_ else set())
        hinged_per_nodo.setdefault(nc_, set()).update(hs_)
    forze_barra = {n_: [0.0, 0.0] for n_ in nodi}
    mom_rigidi = {n_: 0.0 for n_ in nodi}
    mom_cerniera = []
    for el in elementi:
        eid_ = el['id']; n1_, n2_ = int(el['n1']), int(el['n2']); Lf_ = float(Lmap[eid_])
        if Lf_ == 0: continue
        sv_ = np.array([0.0, Lf_])
        Mv_, Tv_, Nv_ = _valuta_poly(M_fin_all[eid_], sv_), _valuta_poly(T_fin_all[eid_], sv_), _valuta_poly(N_fin_all[eid_], sv_)
        ca_ = float((nodi[n2_][0] - nodi[n1_][0]) / Lmap[eid_]); sa_ = float((nodi[n2_][1] - nodi[n1_][1]) / Lmap[eid_])
        ux_, uy_, nx_, ny_ = ca_, sa_, -sa_, ca_
        forze_barra[n1_][0] += Nv_[0] * ux_ - Tv_[0] * nx_; forze_barra[n1_][1] += Nv_[0] * uy_ - Tv_[0] * ny_
        forze_barra[n2_][0] += -Nv_[1] * ux_ + Tv_[1] * nx_; forze_barra[n2_][1] += -Nv_[1] * uy_ + Tv_[1] * ny_
        for nodo_e, m_e in ((n1_, -Mv_[0]), (n2_, Mv_[1])):
            if int(eid_) in hinged_per_nodo.get(nodo_e, set()):
                mom_cerniera.append((nodo_e, eid_, float(m_e)))
            else:
                mom_rigidi[nodo_e] += float(m_e)
    scala_f = max([1.0] + [abs(float(v_)) for n_ in nodi for v_ in (Forze_S0[n_]['Fx'], Forze_S0[n_]['Fy'])])
    scala_m = max([1.0] + [abs(float(Forze_S0[n_]['M'])) for n_ in nodi] + [scala_f * float(max(Lmap.values()))])
    res_nodi, mt_nodo = [], {}
    carrier_nodo = {ci_['nodo']: ci_['carrier'] for ci_ in cerniere_info}
    for n_ in nodi:
        ft_x = float(sp.N(Forze_S0[n_]['Fx'] + X1_exact * Forze_S1[n_]['Fx']))
        ft_y = float(sp.N(Forze_S0[n_]['Fy'] + X1_exact * Forze_S1[n_]['Fy']))
        mt_ = float(sp.N(Forze_S0[n_]['M'] + X1_exact * Forze_S1[n_]['M']))
        mt_nodo[n_] = mt_
        rx_, ry_ = forze_barra[n_][0] + ft_x, forze_barra[n_][1] + ft_y
        car_ = carrier_nodo.get(n_)
        rm_ = mom_rigidi[n_] if car_ is not None else mom_rigidi[n_] + mt_
        res_nodi.append((n_, rx_, ry_, rm_))
    res_cern = [(n_, e_, m_ + (mt_nodo[n_] if carrier_nodo.get(n_) is not None and int(e_) == int(carrier_nodo[n_]) else 0.0)) for n_, e_, m_ in mom_cerniera]
    max_res_f = max([abs(r_[1]) for r_ in res_nodi] + [abs(r_[2]) for r_ in res_nodi] + [0.0])
    max_res_m = max([abs(r_[3]) for r_ in res_nodi] + [abs(c_[2]) for c_ in res_cern] + [0.0])
    nodi_ko = [r_[0] for r_ in res_nodi if max(abs(r_[1]), abs(r_[2])) > 1e-6 * scala_f or abs(r_[3]) > 1e-6 * scala_m]
    cern_ko = [(c_[0], c_[1]) for c_ in res_cern if abs(c_[2]) > 1e-6 * scala_m]
    ULTIMI_RISULTATI['equilibrio_nodi'] = (max_res_f, max_res_m)
    ULTIMI_RISULTATI['nodi_non_in_equilibrio'] = nodi_ko + [c_[0] for c_ in cern_ko]
    if not nodi_ko and not cern_ko:
        print(f"[VERIFICA EQUILIBRIO NODI] OK: residuo max forze={max_res_f:.2e} N, momenti={max_res_m:.2e} N*m")
    else:
        report.append("<div style='background-color:#f8d7da; color:#842029; padding:10px; border-left:4px solid #dc3545; margin-bottom:15px;'>"
                      f"&#10008; Equilibrio ai nodi NON soddisfatto: residuo max forze = {fmt_sci_html(max_res_f)} N, momenti = {fmt_sci_html(max_res_m)} N&middot;m.</div>")

    pendoli = [rf for rf in reaz_fin if rf['tipo'] == 'Fn' and rf.get('pendolo')]
    ULTIMI_RISULTATI['pendoli'] = []
    if pendoli:
        sigma_amm_p = float(dati_json.get('material', {}).get('sigma_amm', 235e6))
        MS_op_p = dati_json.get('material', {}).get('MS_design_operator', dati_json.get('material', {}).get('MS_operator', '>='))
        MS_t_p = float(dati_json.get('material', {}).get('MS_design_target', dati_json.get('material', {}).get('MS_target', 0.0)))
        op_p, _sig_op_p = _simboli_op(MS_op_p)
        ps = dati_json.get('pendolo_section') or {}
        forma = str(ps.get('shape', 'area')); modo_p = str(ps.get('mode', 'design')); k_rett = float(ps.get('ratio', 1) or 1)
        sigma_t_p = (sigma_amm_p / (1 + MS_t_p)) if (1 + MS_t_p > 1e-9) else float('nan')
        report.append("<h3>Pendoli (bielle): tirante o puntone e sezione</h3>")
        rows_pend, righe_sez = [], []
        for rf in pendoli:
            sup_p = supports_n[rf['sup_idx']] if rf.get('sup_idx') is not None else {}
            ang_ax = float(sup_p.get('axis_angle', float(rf['angle']) + 90))
            Np = float(rf['val'])
            funz, col = ("scarico (N = 0)", "#6c757d") if abs(Np) < 1e-6 else (("TIRANTE (trazione)", "#0d6efd") if Np > 0 else ("PUNTONE (compressione)", "#dc3545"))
            rows_pend.append(f"<tr><td>{rf['nome']}</td><td>{rf['nodo']}</td><td>{ang_ax:.2f}&deg;</td><td><b>{fmt_num(Np, 2)}</b></td><td style='color:{col}; font-weight:bold;'>{funz}</td></tr>")
            info_p = {'nodo': rf['nodo'], 'N': Np, 'funzionamento': funz.split()[0].lower(), 'asse_deg': ang_ax}
            if modo_p == 'verifica':
                A_p = float(ps.get('b',0))*float(ps.get('h',0)) if forma=='rettangolare' else (math.pi*float(ps.get('d',0))**2/4 if forma=='circolare' else float(ps.get('A',0)))
                if A_p > 0:
                    sig_p = abs(Np) / A_p; MS_p = (sigma_amm_p / sig_p - 1) if sig_p > 0 else float('inf')
                    info_p.update({'A': A_p, 'sigma': sig_p, 'MS': MS_p, 'ok': _verifica_op(MS_p, MS_op_p, MS_t_p)})
            elif abs(Np) >= 1e-6 and sigma_t_p == sigma_t_p:
                A_req = abs(Np) / sigma_t_p; info_p['A_req'] = A_req
                if forma == 'rettangolare': info_p.update({'b': math.sqrt(A_req/k_rett), 'h': k_rett*math.sqrt(A_req/k_rett)})
                elif forma == 'circolare': info_p['d'] = math.sqrt(4*A_req/math.pi)
                righe_sez.append(f"<p><b>{rf['nome']} (nodo {rf['nodo']}):</b> {funz} &mdash; A {op_p} {fmt_sci_html(A_req, sig=4)} m&sup2; ({A_req*1e4:.2f} cm&sup2;)</p>")
            ULTIMI_RISULTATI['pendoli'].append(info_p)
        report.append(_tabella_html(["Pendolo", "Nodo", "Asse", "N [N]", "Funzionamento"], rows_pend, keep=True))
        report.extend(righe_sez)

    # Identazione sentinella per il blocco Sezione 8 e 9 sotto
    if True:

        # ==========================================
        # FASE 9: SEZIONE CRITICA E VERIFICA / DIMENSIONAMENTO (Condizionale)
        # ==========================================
        ver_type = dati_json.get('material', {}).get('verification_type', 'combined')
        b_sec = float(dati_json.get('section', {}).get('b', 0.2))
        h_sec = float(dati_json.get('section', {}).get('h', 0.4))
        Area = b_sec * h_sec
        W_f = (b_sec * (h_sec**2)) / 6

        M_fin_dict, N_fin_dict, L_dict = {}, {}, {}
        for el in elementi:
            el_id = el['id']
            L_dict[el_id] = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
            if grado_iperstaticita > 0:
                M_fin_dict[el_id] = sp.simplify(M0_dict[el_id] + X1_exact * M1_dict[el_id])
                N_fin_dict[el_id] = sp.simplify(N0_dict[el_id] + X1_exact * N1_dict[el_id])
            else:
                M_fin_dict[el_id] = M0_dict[el_id]
                N_fin_dict[el_id] = N0_dict[el_id]

        PESI_DECADI = [(w_, 1.0) for w_ in (0.0, 1e-3, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3)]

        def _scansiona(pesi_):
            punti = []
            for el_ in elementi:
                eid_ = el_['id']
                sv = candidati_s(M_fin_dict[eid_], N_fin_dict[eid_], L_dict[eid_], pesi=pesi_, n_base=101)
                punti.append((eid_, sv, _valuta_poly(N_fin_dict[eid_], sv), _valuta_poly(M_fin_dict[eid_], sv)))
            return punti

        pesi_ver = [(0.0, 1.0)] if ver_type == 'bending_only' else [(h_sec / 6.0, 1.0)]
        punti_NM = _scansiona(PESI_DECADI + pesi_ver)

        def _sigma_punti(N_, M_):
            return np.abs(M_) / W_f if ver_type == 'bending_only' else np.abs(N_) / Area + np.abs(M_) / W_f

        max_sigma_calcolata = 0.0
        asta_critica = None
        s_critico = 0.0
        M_in_s_critico = 0.0
        N_in_s_critico = 0.0
        sezioni_equivalenti = []

        for el_id, s_vals, N_vals, M_vals in punti_NM:
            Sigma_vals = _sigma_punti(N_vals, M_vals)
            idx_max = int(np.argmax(Sigma_vals))
            local_max_sigma = float(Sigma_vals[idx_max])
            if local_max_sigma > max_sigma_calcolata * (1 + 1e-9):
                max_sigma_calcolata = local_max_sigma
                asta_critica = el_id
                s_critico = float(s_vals[idx_max])
                M_in_s_critico = float(M_vals[idx_max])
                N_in_s_critico = float(N_vals[idx_max])
        if max_sigma_calcolata > 0:
            for el_id, s_vals, N_vals, M_vals in punti_NM:
                Sigma_vals = _sigma_punti(N_vals, M_vals)
                vicini = np.where(Sigma_vals >= max_sigma_calcolata * (1 - 1e-6))[0]
                if len(vicini):
                    sezioni_equivalenti.append((el_id, float(s_vals[vicini[0]])))
        else:
            asta_critica = None

        sigma_amm = float(dati_json.get('material', {}).get('sigma_amm', 235e6))
        MS = ((sigma_amm / max_sigma_calcolata) - 1) if max_sigma_calcolata > 0 else float('inf')
        MS_target = float(dati_json.get('material', {}).get('MS_target', 0.0))
        MS_op = dati_json.get('material', {}).get('MS_operator', '>=')

        ULTIMI_RISULTATI['MS'] = MS
        ULTIMI_RISULTATI['sigma_max'] = max_sigma_calcolata
        ULTIMI_RISULTATI['sezione_critica'] = (asta_critica, s_critico, N_in_s_critico, M_in_s_critico)

        op_html, sigma_op_html = _simboli_op(MS_op)
        TOL_MS = 1e-9
        if MS_op == ">=": is_verified = (MS >= MS_target - TOL_MS)
        elif MS_op == ">": is_verified = (MS > MS_target + TOL_MS)
        elif MS_op == "<=": is_verified = (MS <= MS_target + TOL_MS)
        elif MS_op == "<": is_verified = (MS < MS_target - TOL_MS)
        else: is_verified = math.isclose(MS, MS_target, abs_tol=1e-2) if MS != float('inf') else False

        _mat_ms = dati_json.get('material', {})
        MS_d_target = float(_mat_ms.get('MS_design_target', MS_target))
        MS_d_op = _mat_ms.get('MS_design_operator', MS_op)
        op_html_d, sigma_op_html_d = _simboli_op(MS_d_op)

        def _raccogli_top_sezioni(punti_, score_fn=None, top_n=3):
            """Restituisce fino a `top_n` sezioni distinte piu' sollecitate (una o piu' aste di picco)."""
            candidati = []
            max_M_g = max([float(np.max(np.abs(p_[3]))) for p_ in punti_] + [1.0])
            max_N_g = max([float(np.max(np.abs(p_[2]))) for p_ in punti_] + [1.0])
            for eid_, sv_, Nv_, Mv_ in punti_:
                # Punti chiave su quest'asta: massimo |M|, massimo |N|, e massimo punteggio combinato
                idxs = {0, len(sv_) - 1, int(np.argmax(np.abs(Mv_))), int(np.argmax(np.abs(Nv_)))}
                if score_fn is not None:
                    sc_arr = score_fn(eid_, sv_, Nv_, Mv_)
                    idxs.add(int(np.argmax(sc_arr)))
                for ii in idxs:
                    nv_a, mv_a = abs(float(Nv_[ii])), abs(float(Mv_[ii]))
                    if nv_a < 1e-6 and mv_a < 1e-6:
                        continue
                    sc = float(score_fn(eid_, np.array([sv_[ii]]), np.array([Nv_[ii]]), np.array([Mv_[ii]]))[0]) if score_fn is not None else (mv_a / max_M_g + 0.25 * nv_a / max_N_g)
                    candidati.append((sc, int(eid_), float(sv_[ii]), nv_a, mv_a))
            candidati.sort(key=lambda x: (-x[0], -x[4], -x[3], x[1]))
            scelti = []
            for item in candidati:
                _, eid_, s_, nv_a, mv_a = item
                # Evita duplicati sulla stessa asta nella stessa zona o valori quasi identici sullo stesso nodo giunto
                if any(c[1] == eid_ and abs(c[2] - s_) < 1e-3 for c in scelti):
                    continue
                scelti.append(item)
                if len(scelti) >= top_n:
                    break
            return scelti

        design = str(dati_json.get('section', {}).get('design', 'none'))
        in_progetto = design in ('h', 'b', 'ratio')

        def eq_table(*elements):
            cells = "".join(
                f"<td valign='middle' style='border:none; padding:0 4px; background:transparent; width:auto;'>{el}</td>"
                for el in elements
            )
            return (
                "<table class='eq-table' cellpadding='0' cellspacing='0' "
                f"style='margin-bottom: 10px; border:none; width:auto; background:transparent; font-family:{_MONO_CSS}; font-size: 15px; color: #000000;'>"
                f"<tr>{cells}</tr></table>"
            )

        report.append("<h2 style='margin-top: 26px;'>8. Sezione pi&ugrave; sollecitata e verifica di resistenza (Punto 3 &mdash; Navier)</h2>")

        # ==========================================
        # PUNTO 3: VERIFICA SEMPRE ESEGUITA SULLA SEZIONE ASSEGNATA
        # ==========================================
        if ver_type == 'bending_only':
            report.append("<p>Criterio: <b>Flessione semplice</b> (&sigma;(s) = |M(s)| / W).</p>")
        else:
            report.append("<p>Criterio: <b>Presso/Tensoflessione combinata</b> (&sigma;(s) = |N(s)| / A + |M(s)| / W).</p>")

        top_ver = _raccogli_top_sezioni(punti_NM, lambda e_, s_, n_, m_: _sigma_punti(n_, m_), top_n=3)
        report.append("<div class='keep-together' style='background-color: #f8f9fa; padding: 14px 16px; border-left: 5px solid #0a58ca; border-radius: 4px; margin-bottom: 16px; page-break-inside: avoid;'>")
        if asta_critica is not None:
            if top_ver:
                rows_tv = []
                for sc_v, e_v, s_v, n_v, m_v in top_ver:
                    is_gov = (e_v == asta_critica and abs(s_v - s_critico) < 1e-4)
                    ruolo = "<b style='color:#0a58ca;'>Sezione critica</b>" if is_gov else "Confronto"
                    rows_tv.append(f"<tr><td>Asta {e_v}</td><td>{fmt_num(s_v, 2)}</td><td>{fmt_num(n_v, 2)}</td><td>{fmt_num(m_v, 2)}</td><td><b>{fmt_smart(sc_v/1e6)}</b></td><td>{ruolo}</td></tr>")
                report.append(_tabella_html(["Asta", "s [m]", "|N| [N]", "|M| [N&middot;m]", "&sigma; [MPa]", "Ruolo"], rows_tv, keep=True, compact=True))
            report.append(f"<p>La tensione normale massima si raggiunge nell'<b>Asta {asta_critica}</b> a <b>s = {fmt_num(s_critico, 2)} m</b> "
                          f"(|N| = {fmt_num(abs(N_in_s_critico), 2)} N, |M| = {fmt_num(abs(M_in_s_critico), 2)} N&middot;m).</p>")
        else:
            report.append("<p>Sollecitazioni nulle su tutta la struttura.</p>")

        _area_h = fmt_smart(Area, min_sig=4)
        report.append(f"<p>Sezione assegnata: b = {fmt_smart(b_sec)} m, h = {fmt_smart(h_sec)} m &rarr; "
                      f"A = {_area_h} m<sup>2</sup>, W = b&middot;h<sup>2</sup>/6 = {fmt_sci_html(W_f, sig=4)} m<sup>3</sup>, "
                      f"&sigma;<sub>amm</sub> = {fmt_smart(sigma_amm/1e6)} MPa.</p>")
        report.append("<div style='margin: 8px 0 12px 0;'>")
        if ver_type == 'bending_only':
            report.append(f"<span class='formula' style='color:#000000;'>&sigma;<sub>max</sub> = |M| / W = {fmt_num(abs(M_in_s_critico), 2)} / ({fmt_sci_html(W_f, sig=4)}) = <b>{fmt_smart(max_sigma_calcolata/1e6)} MPa</b></span><br><br>")
        else:
            report.append(f"<span class='formula' style='color:#000000;'>&sigma;<sub>max</sub> = |N|/A + |M|/W = {fmt_num(abs(N_in_s_critico), 2)}/({_area_h}) + {fmt_num(abs(M_in_s_critico), 2)}/({fmt_sci_html(W_f, sig=4)}) = <b>{fmt_smart(max_sigma_calcolata/1e6)} MPa</b></span><br><br>")
        if MS != float('inf'):
            report.append(f"<span class='formula' style='color:#000000;'>MS = (&sigma;<sub>amm</sub> / &sigma;<sub>max</sub>) &minus; 1 = {fmt_smart(sigma_amm/1e6)} / {fmt_smart(max_sigma_calcolata/1e6)} &minus; 1 = <b>{fmt_num(MS, 2)}</b></span>")
        else:
            report.append("<span class='formula' style='color:#000000;'>MS = &infin;</span>")
        report.append("</div>")

        if is_verified:
            report.append("<div style='text-align:center; padding:12px; background-color:#d1e7dd; color:#0f5132; border: 2px solid #198754; font-size: 16px; border-radius: 6px;'>")
            report.append(f"<b>MS = {fmt_num(MS, 2)} {op_html} {fmt_num(MS_target, 2)}</b> (MS verifica punto 3)<br>VERIFICATA &#10004;</div>")
        else:
            report.append("<div style='text-align:center; padding:12px; background-color:#f8d7da; color:#842029; border: 2px solid #dc3545; font-size: 16px; border-radius: 6px;'>")
            report.append(f"<b>MS = {fmt_num(MS, 2)} (Richiesto {op_html} {fmt_num(MS_target, 2)}, MS verifica punto 3)</b><br>NON VERIFICATA &#10008;</div>")
        report.append("</div>")

        if in_progetto:
            # ==========================================
            # PUNTO 4: DIMENSIONAMENTO (SEZIONE 9)
            # ==========================================
            sigma_amm_val = float(dati_json.get('material', {}).get('sigma_amm', 235e6))
            sigma_target_val = sigma_amm_val / (1 + MS_d_target) if 1 + MS_d_target > 1e-9 else float('nan')
            ratio_k = float(dati_json.get('section', {}).get('ratio', 0) or 0)
            combinato = (ver_type != 'bending_only')
            EI_dir = 'EI' in dati_json.get('material', {})
            sez_dip = (dip_ced and not EI_dir) or (dip_term and (not EI_dir or design in ('h', 'ratio')))
            esito_design = None
            problema_design = None

            if 1 + MS_d_target <= 1e-9:
                problema_design = ('errore', "Dimensionamento impossibile", "Con MS richiesto &le; &minus;1 la tensione obiettivo non &egrave; definita.")
            elif design == 'ratio' and ratio_k <= 0:
                problema_design = ('errore', "Rapporto h/b mancante", "Per il dimensionamento a rapporto fissato indica h/b &gt; 0.")
            elif sez_dip:
                s_d = sp.Symbol('s')
                def _cf(e_):
                    c_ = _coeff_poly(e_, s_d)
                    return c_ if c_ is not None else [0.0]
                barre_d = [(el_['id'], float(L_dict[el_['id']]), _cf(M0_dict[el_['id']]), _cf(M1_dict[el_['id']]),
                            _cf(N0_dict[el_['id']]), _cf(N1_dict[el_['id']])) for el_ in elementi]
                d10_f, d11_f = float(delta_10), float(delta_11)
                KT_f, TN_f, SR_f, eta_f = float(KT_sum), float(termine_N), float(somma_Rd), float(eta_X1)
                E_f = float(dati_json.get('material', {}).get('E', 2.1e11))
                EI_fisso = float(EI_val) if EI_dir else None
                solo_fl = (ver_type == 'bending_only')
                nome_d = 'h' if design == 'h' else 'b'

                def _sez_d(d_):
                    if design == 'b': return d_, h_sec
                    if design == 'h': return b_sec, d_
                    return d_, ratio_k * d_

                def _X1_d(d_):
                    b_, h_ = _sez_d(d_)
                    EI_ = EI_fisso if EI_fisso is not None else E_f * b_ * h_ ** 3 / 12.0
                    return (EI_ * (eta_f - KT_f / h_ - TN_f + SR_f) - d10_f) / d11_f

                def _sig_d(d_):
                    b_, h_ = _sez_d(d_)
                    return _sigma_max_numerico(barre_d, _X1_d(d_), b_ * h_, b_ * h_ ** 2 / 6.0, solo_fl)

                D_LO, D_HI = 1e-3, 5.0
                radici_d, interv_d = _intervalli_soddisfatti(lambda d_: _sig_d(d_)[0] - sigma_target_val, D_LO, D_HI, MS_d_op)

                def _txt_interv(a_, b_):
                    if a_ is None and b_ is None: return f"qualsiasi {nome_d} in [{D_LO*100:g} cm, {D_HI:g} m]"
                    if a_ is None: return f"{nome_d} &le; {b_:.4f} m ({b_*100:.2f} cm)"
                    if b_ is None: return f"{nome_d} &ge; {a_:.4f} m ({a_*100:.2f} cm)"
                    return f"{a_:.4f} m &le; {nome_d} &le; {b_:.4f} m ({a_*100:.2f} cm &le; {nome_d} &le; {b_*100:.2f} cm)"

                if MS_d_op == '=':
                    cand_d = sorted(radici_d)
                    testo_d = (" &nbsp;oppure&nbsp; ".join(f"{nome_d} = {r_:.4f} m ({r_*100:.2f} cm)" for r_ in cand_d)) if cand_d else ""
                else:
                    cand_d = [x_ for iv_ in interv_d for x_ in iv_ if x_ is not None][:1]
                    testo_d = " &nbsp;oppure&nbsp; ".join(_txt_interv(a_, b_) for a_, b_ in interv_d)
                if not cand_d:
                    if interv_d:
                        problema_design = ('avviso', "Dimensionamento: nessun limite",
                                           f"Con le distorsioni assegnate la richiesta &egrave; soddisfatta da {testo_d}.")
                    else:
                        problema_design = ('errore', "Dimensionamento senza soluzione",
                                           f"Nessuna dimensione tra {D_LO*100:g} cm e {D_HI:g} m soddisfa la richiesta.")
                else:
                    d_ref = float(cand_d[0])
                    b_ref, h_ref = _sez_d(d_ref)
                    sig_t_, asta_t_, s_t_, N_t_, M_t_ = _sig_d(d_ref)
                    MS_c = (sigma_amm / sig_t_ - 1) if sig_t_ > 0 else float('inf')
                    X1_ref = _X1_d(d_ref)
                    esito_design = {'b': b_ref, 'h': h_ref, 'd': d_ref, 'MS_c': MS_c, 'asta': int(asta_t_), 's': s_t_,
                                    'N': abs(N_t_), 'M': abs(M_t_), 'sig_t': sigma_target_val, 'distorsioni': True,
                                    'X1': X1_ref, 'testo': testo_d, 'radici': list(radici_d), 'intervalli': list(interv_d)}
            else:
                sig_t = sigma_target_val

                def _dimensioni(punti_):
                    n_tot = sum(len(p_[1]) for p_ in punti_)
                    N_all = np.concatenate([np.abs(p_[2]) for p_ in punti_]) if combinato else np.zeros(n_tot)
                    M_all = np.concatenate([np.abs(p_[3]) for p_ in punti_])
                    s_all = np.concatenate([p_[1] for p_ in punti_])
                    id_all = np.concatenate([np.full(len(p_[1]), p_[0]) for p_ in punti_])
                    if design == 'b':
                        d_all = (N_all * h_sec + 6 * M_all) / (h_sec ** 2 * sig_t)
                    elif design == 'h':
                        d_all = (N_all + np.sqrt(N_all ** 2 + 24 * sig_t * b_sec * M_all)) / (2 * sig_t * b_sec)
                    else:
                        d_all = np.zeros(len(M_all))
                        for i_, (n_i, m_i) in enumerate(zip(N_all, M_all)):
                            if n_i == 0 and m_i == 0: continue
                            rad = np.roots([sig_t, 0.0, -n_i / ratio_k, -6.0 * m_i / ratio_k ** 2])
                            rad = [r_.real for r_ in rad if abs(r_.imag) < 1e-9 * max(1.0, abs(r_)) and r_.real > 0]
                            d_all[i_] = max(rad) if rad else 0.0
                    return N_all, M_all, s_all, id_all, d_all

                def _sezione(d_):
                    if design == 'b': return d_, h_sec
                    if design == 'h': return b_sec, d_
                    return d_, ratio_k * d_

                N_all, M_all, s_all, id_all, d_all = _dimensioni(punti_NM)
                d_ref = float(np.max(d_all))
                for _giro in range(3):
                    if not combinato or d_ref <= 0: break
                    b_r, h_r = _sezione(d_ref)
                    punti2 = _scansiona(PESI_DECADI + [(h_r / 6.0, 1.0)])
                    r2 = _dimensioni(punti2)
                    if float(np.max(r2[4])) <= d_ref * (1 + 1e-12): break
                    N_all, M_all, s_all, id_all, d_all = r2
                    d_ref = float(np.max(d_all))
                i_c = int(np.argmax(d_all))
                b_ref, h_ref = _sezione(d_ref)
                A_c, W_c = b_ref * h_ref, b_ref * h_ref ** 2 / 6
                sig_c = float(np.max(N_all / A_c + M_all / W_c)) if d_ref > 0 else 0.0
                MS_c = (sigma_amm / sig_c - 1) if sig_c > 0 else float('inf')
                esito_design = {'b': b_ref, 'h': h_ref, 'd': d_ref, 'MS_c': MS_c, 'asta': int(id_all[i_c]), 's': float(s_all[i_c]),
                                'N': float(N_all[i_c]), 'M': float(M_all[i_c]), 'sig_t': sig_t}

            if esito_design is not None:
                N_E_abs, M_E_abs = esito_design['N'], esito_design['M']
                asta_gov, s_gov = esito_design['asta'], esito_design['s']
            else:
                N_E_abs, M_E_abs = abs(N_in_s_critico), abs(M_in_s_critico)
                asta_gov, s_gov = asta_critica, s_critico
            sei_Me = 6 * M_E_abs

            # ==========================================
            # SEZIONE 9: DIMENSIONAMENTO DELLA SEZIONE (PUNTO 4)
            # ==========================================
            report.append("<h2 style='margin-top: 24px;'>9. Dimensionamento della sezione (Punto 4)</h2>")
            report.append("<div class='keep-together' style='background-color: #f8f9fa; padding: 14px 16px; border-left: 5px solid #0a58ca; border-radius: 4px; margin-bottom: 18px; page-break-inside: avoid;'>")
            sigma_target_html = fmt_sci_html(sigma_target_val, sig=4) if sigma_target_val == sigma_target_val else "n.d."
            bl_style = "border-left: 3px solid #0d6efd; background-color: #ffffff; padding: 10px 15px; margin: 10px 0 14px 0;"

            if asta_gov is not None:
                stessa_crit = (asta_critica is not None and asta_gov == asta_critica and abs(s_gov - s_critico) < 1e-3)
                if stessa_crit:
                    report.append(
                        f"<p style='margin:0 0 8px 0;'>Sezione critica governante (invariata rispetto al Punto 3): "
                        f"<b>Asta {asta_gov}</b> a <b>s = {fmt_num(s_gov, 2)} m</b> "
                        f"(|N<sub>E</sub>| = {fmt_num(N_E_abs, 2)} N, |M<sub>E</sub>| = {fmt_num(M_E_abs, 2)} N&middot;m).</p>"
                    )
                else:
                    report.append(
                        f"<p style='margin:0 0 8px 0; color:#664d03; background-color:#fff3cd; padding:6px 10px; border-left:3px solid #ffc107;'>"
                        f"<b>Cambio della sezione critica rispetto al Punto 3:</b> al variare delle dimensioni della sezione "
                        f"il dimensionamento &egrave; governato dall'<b>Asta {asta_gov}</b> a <b>s = {fmt_num(s_gov, 2)} m</b> "
                        f"(|N<sub>E</sub>| = {fmt_num(N_E_abs, 2)} N, |M<sub>E</sub>| = {fmt_num(M_E_abs, 2)} N&middot;m) "
                        f"invece che dall'Asta {asta_critica} (s = {fmt_num(s_critico, 2)} m).</p>"
                    )

            report.append(f"<p>Imponendo il requisito sul margine di sicurezza richiesto (Punto 4: <b>MS {op_html_d} {fmt_num(MS_d_target, 2)}</b>) e la formula di Navier con <i>A = b&middot;h</i> e <i>W = b&middot;h<sup>2</sup>/6</i>:</p>")
            report.append(f"<div style=\"{bl_style}\">")
            report.append(eq_table(f"MS {op_html_d} {fmt_num(MS_d_target, 2)} &rArr;",
                                   _frac_html("&sigma;<sub>amm</sub>", "&sigma;<sub>max</sub>", "#000000"),
                                   f"&minus; 1 {op_html_d} {fmt_num(MS_d_target, 2)} &rArr; &sigma;<sub>max</sub> {sigma_op_html_d}",
                                   _frac_html("&sigma;<sub>amm</sub>", f"1 + {fmt_num(MS_d_target, 2)}", "#000000"),
                                   f"= {fmt_smart(sigma_target_val/1e6)} MPa ({sigma_target_html} Pa)"))
            if ver_type == 'bending_only':
                report.append(eq_table("&sigma;<sub>max</sub> =", _frac_html("|M<sub>E</sub>|", "W", "#000000"), "=",
                                       _frac_html("6&middot;|M<sub>E</sub>|", "b&middot;h<sup>2</sup>", "#000000"), "=",
                                       _frac_html(f"{fmt_num(sei_Me, 2)}", "b&middot;h<sup>2</sup>", "#000000"),
                                       sigma_op_html_d, f"{sigma_target_html} Pa"))
            else:
                report.append(eq_table("&sigma;<sub>max</sub> =", _frac_html("|N<sub>E</sub>|", "A", "#000000"), "+",
                                       _frac_html("|M<sub>E</sub>|", "W", "#000000"), "=",
                                       _frac_html(f"h&middot;{fmt_num(N_E_abs, 2)} + {fmt_num(sei_Me, 2)}", "b&middot;h<sup>2</sup>", "#000000"),
                                       sigma_op_html_d, f"{sigma_target_html} Pa"))
            report.append("</div>")

            if not sez_dip and esito_design is not None:
                report.append(f"<div style=\"{bl_style}\">")
                if design == 'ratio':
                    c3 = sigma_target_val
                    c1 = N_E_abs / ratio_k if combinato else 0.0
                    c0 = 6.0 * M_E_abs / (ratio_k ** 2)
                    ineq_pol = "&ge;" if MS_d_op in (">", ">=") else ("&le;" if MS_d_op in ("<", "<=") else "=")
                    report.append(f"<p style='margin:0 0 6px 0;'>Sostituendo <b>h = {fmt_num(ratio_k, 2)}&middot;b</b> si ottiene l'equazione cubica in <i>b</i>:</p>")
                    if combinato and c1 > 1e-9:
                        report.append(eq_table("&sigma;<sub>t</sub>&middot;b<sup>3</sup> &minus;",
                                               _frac_html("|N<sub>E</sub>|", "k", "#000000"),
                                               "&middot;b &minus;",
                                               _frac_html("6&middot;|M<sub>E</sub>|", "k<sup>2</sup>", "#000000"),
                                               f"{ineq_pol} 0 &rArr; {sigma_target_html}&middot;b<sup>3</sup> &minus; {fmt_num(c1, 2)}&middot;b &minus; {fmt_num(c0, 2)} {ineq_pol} 0"))
                    else:
                        report.append(eq_table("&sigma;<sub>t</sub>&middot;b<sup>3</sup> &minus;",
                                               _frac_html("6&middot;|M<sub>E</sub>|", "k<sup>2</sup>", "#000000"),
                                               f"{ineq_pol} 0 &rArr; {sigma_target_html}&middot;b<sup>3</sup> &minus; {fmt_num(c0, 2)} {ineq_pol} 0"))
                elif design == 'h':
                    ineq_pol = "&ge;" if MS_d_op in (">", ">=") else ("&le;" if MS_d_op in ("<", "<=") else "=")
                    report.append(f"<p style='margin:0 0 6px 0;'>Con <b>b = {fmt_num(b_sec, 3)} m</b> fissata si ottiene l'equazione in <i>h</i>:</p>")
                    if combinato and N_E_abs > 1e-9:
                        report.append(eq_table(f"(&sigma;<sub>t</sub>&middot;b)&middot;h<sup>2</sup> &minus; |N<sub>E</sub>|&middot;h &minus; 6|M<sub>E</sub>| {ineq_pol} 0 &rArr; {fmt_sci_html(sigma_target_val * b_sec, sig=4)}&middot;h<sup>2</sup> &minus; {fmt_num(N_E_abs, 2)}&middot;h &minus; {fmt_num(sei_Me, 2)} {ineq_pol} 0"))
                    else:
                        report.append(eq_table(f"(&sigma;<sub>t</sub>&middot;b)&middot;h<sup>2</sup> &minus; 6|M<sub>E</sub>| {ineq_pol} 0 &rArr; {fmt_sci_html(sigma_target_val * b_sec, sig=4)}&middot;h<sup>2</sup> &minus; {fmt_num(sei_Me, 2)} {ineq_pol} 0"))
                else:
                    report.append(f"<p style='margin:0 0 6px 0;'>Con <b>h = {fmt_num(h_sec, 3)} m</b> fissata si ricava direttamente <i>b</i>:</p>")
                    num_b = f"{fmt_num(h_sec, 3)}&middot;{fmt_num(N_E_abs, 2)} + {fmt_num(sei_Me, 2)}" if combinato else f"{fmt_num(sei_Me, 2)}"
                    report.append(eq_table(f"b {op_html_d}", _frac_html(num_b, f"{fmt_sci_html(sigma_target_val * h_sec**2, sig=4)}", "#000000")))
                report.append("</div>")

            lett_attiva = (dip_ced and not EI_dir) or dip_term
            if lett_attiva and 1 + MS_d_target > 1e-9:
                try:
                    s_l = sp.Symbol('s')
                    def _cfl(e_):
                        c_ = _coeff_poly(e_, s_l)
                        return c_ if c_ is not None else [0.0]
                    barre_l = [(el_['id'], float(L_dict[el_['id']]), _cfl(M0_dict[el_['id']]), _cfl(M1_dict[el_['id']]),
                                _cfl(N0_dict[el_['id']]), _cfl(N1_dict[el_['id']])) for el_ in elementi]
                    d10_l, d11_l = float(delta_10), float(delta_11)
                    KT_l, TN_l, SR_l, eta_l = float(KT_sum), float(termine_N), float(somma_Rd), float(eta_X1)
                    E_l = float(dati_json.get('material', {}).get('E', 2.1e11))
                    b_r, h_r = (esito_design['b'], esito_design['h']) if (esito_design is not None and esito_design.get('distorsioni')) else (b_sec, h_sec)
                    EI_l = float(EI_val) if EI_dir else None
                    def _X1_num(b_, h_):
                        EI_ = EI_l if EI_l is not None else E_l * b_ * h_ ** 3 / 12.0
                        return (EI_ * (eta_l - KT_l / h_ - TN_l + SR_l) - d10_l) / d11_l
                    _, asta_r, s_r, _Nr, _Mr = _sigma_max_numerico(barre_l, _X1_num(b_r, h_r), b_r * h_r, b_r * h_r ** 2 / 6.0, ver_type == 'bending_only')
                    bar_r = next(x_ for x_ in barre_l if x_[0] == asta_r)
                    n0_, n1_ = float(np.polyval(bar_r[4], s_r)), float(np.polyval(bar_r[5], s_r))
                    m0_, m1_ = float(np.polyval(bar_r[2], s_r)), float(np.polyval(bar_r[3], s_r))
                    b_s, h_s = sp.symbols('b h', positive=True)
                    EI_s = sp.Float(EI_l) if EI_l is not None else sp.Float(E_l) * b_s * h_s ** 3 / 12
                    X1_s = (EI_s * (eta_l - KT_l / h_s - TN_l + SR_l) - d10_l) / d11_l
                    N_s, M_s = n0_ + X1_s * n1_, m0_ + X1_s * m1_
                    sN_, sM_ = (1.0 if float(_Nr) >= 0 else -1.0), (1.0 if float(_Mr) >= 0 else -1.0)
                    lhs_s = sp.expand((6 * sM_ * M_s) if ver_type == 'bending_only' else (sN_ * N_s * h_s + 6 * sM_ * M_s))
                    def _fmt_bh(e_):
                        e_ = sp.expand(e_)
                        termini_ = [(float(t_.as_coeff_Mul()[0]), t_.as_coeff_Mul()[1]) for t_ in sp.Add.make_args(e_)]
                        out_ = ""
                        for c_, m_ in sorted(termini_, key=lambda x_: -sp.degree(x_[1], h_s) - sp.degree(x_[1], b_s)):
                            if abs(c_) < 1e-9: continue
                            vp_ = []
                            for f_ in m_.as_ordered_factors():
                                if isinstance(f_, sp.Pow): vp_.append(f"{f_.base}<sup>{f_.exp}</sup>")
                                elif f_ != 1: vp_.append(str(f_))
                            mon_ = "".join(vp_)
                            c_str = fmt_sci_html(abs(c_), sig=4)
                            sep = "&middot;" if mon_ else ""
                            out_ += f" {'+' if c_ >= 0 else '&minus;'} {c_str}{sep}{mon_}"
                        out_ = out_.strip()
                        return out_[1:].strip() if out_.startswith('+') else out_
                    report.append("<p>Con cedimenti o variazioni termiche l'incognita iperstatica dipende dalla sezione: "
                                  f"nella sezione critica (Asta {asta_r}, s = {fmt_num(s_r, 2)} m) la relazione esplicita tra <i>b</i> e <i>h</i> diventa:</p>")
                    report.append(f"<div style=\"{bl_style}\">")
                    report.append(eq_table("X<sub>1</sub>(b,h) =", _fmt_bh(X1_s)))
                    report.append(eq_table(f"{_fmt_bh(lhs_s)} {sigma_op_html_d} {sigma_target_html}&middot;b&middot;h<sup>2</sup>"))
                    report.append("</div>")
                    ULTIMI_RISULTATI['relazione_letterale'] = {'lhs': str(lhs_s), 'sigma_t': sigma_target_val, 'asta': int(asta_r), 's': float(s_r),
                                                               'X1': str(sp.expand(X1_s)), 'op': MS_d_op}
                except Exception as _e_l:
                    report.append(_box_avviso("Relazione letterale non generata", f"Non e' stato possibile scrivere la relazione in b e h: {_e_l}"))

            if problema_design is not None:
                report.append(_box_avviso(problema_design[1], problema_design[2]) if problema_design[0] == 'avviso'
                              else _box_errore(problema_design[1], problema_design[2]))
            elif esito_design is not None:
                b_ref, h_ref, MS_c = esito_design['b'], esito_design['h'], esito_design['MS_c']
                if esito_design.get('distorsioni'):
                    report.append(f"<p>Per la sezione progettata con distorsioni imposte: <b>X<sub>1</sub> = {fmt_smart(esito_design['X1'])} {unita_X1}</b> "
                                  f"(sezione critica: Asta {esito_design['asta']}, s = {fmt_num(esito_design['s'], 2)} m).</p>")
                    esito = esito_design['testo']
                    if design == 'ratio':
                        esito += f"<br><small>(con h = {fmt_smart(ratio_k)}&middot;b; per b<sub>min</sub> = {b_ref:.4f} m: h<sub>min</sub> = {h_ref:.4f} m)</small>"
                    elif design == 'b':
                        esito += f"<br><small>(con h = {fmt_smart(h_sec)} m assegnata)</small>"
                    else:
                        esito += f"<br><small>(con b = {fmt_smart(b_sec)} m assegnata)</small>"
                elif design == 'b':
                    esito = f"b<sub>min</sub> = {b_ref:.4f} m ({b_ref*100:.2f} cm) &rArr; b {op_html_d} {b_ref:.4f} m &nbsp;(h = {fmt_smart(h_sec)} m)"
                elif design == 'h':
                    esito = f"h<sub>min</sub> = {h_ref:.4f} m ({h_ref*100:.2f} cm) &rArr; h {op_html_d} {h_ref:.4f} m &nbsp;(b = {fmt_smart(b_sec)} m)"
                else:
                    esito = f"b<sub>min</sub> = {b_ref:.4f} m ({b_ref*100:.2f} cm), &nbsp;h<sub>min</sub> = {fmt_smart(ratio_k)}&middot;b<sub>min</sub> = {h_ref:.4f} m ({h_ref*100:.2f} cm)"
                report.append("<div style='text-align:center; padding:12px 15px; background-color:#d1e7dd; color:#0f5132; border: 2px solid #198754; font-size: 16px; margin-top: 10px; border-radius: 6px;'>")
                report.append(f"<b>{esito}</b><br><small>Controllo finale: MS = {fmt_num(MS_c, 2)} (MS<sub>req</sub> {op_html_d} {fmt_num(MS_d_target, 2)})</small></div>")
                ULTIMI_RISULTATI['dimensionamento'] = dict(esito_design)
            report.append("</div>")

    return "\n".join(report), percorso_S1, percorso_diag

def risolvi_struttura_plv(dati_json, out_dir, incognita):
    """Entry point (stessa firma di prima). Gli errori bloccanti tornano come report HTML con box rosso,
    senza immagini (img_s1 = img_diag = "")."""
    ULTIMI_RISULTATI.clear()
    try:
        return _risolvi_struttura_plv_core(dati_json, out_dir, incognita)
    except ErroreSolutore as e:
        return str(e), "", ""

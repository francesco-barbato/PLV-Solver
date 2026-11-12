import os
import json
import math

import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

import sympy as sp
import numpy as np

# =====================================================================
# FASE 1 - ARITMETICA ESATTA, PULIZIA RUMORE, ERRORI BLOCCANTI
# =====================================================================
TOL_RUMORE = 1e-10
_COSTANTI = [sp.sqrt(2), sp.sqrt(3), sp.sqrt(5)]

ULTIMI_RISULTATI = {}   # riempito a ogni calcolo: reazioni finali, X1, MS, pendoli, residui dei controlli

class ErroreSolutore(Exception):
    """Errore bloccante: il messaggio contiene il report HTML parziale + box d'errore."""
    pass

def ex(v, default=0):
    """Converte un numero (int/float/str/SymPy) in valore ESATTO SymPy.
    4.0 -> 4 ; 0.5 -> 1/2 ; 2.8284271247461903 -> 2*sqrt(2) ; 1.7320508 -> sqrt(3)."""
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

def pulisci_formula(expr):
    """Pulisce le espressioni SymPy per renderle leggibili in HTML nativo."""
    testo = str(expr)
    testo = testo.replace('**', '^')
    testo = testo.replace('*', '·')
    testo = testo.replace('sqrt', '√')
    # Sostituzioni per abbellire le variabili
    testo = testo.replace('X1', 'X<sub>1</sub>')
    return testo

def plot_sistema_S0(nodi, elementi, reazioni, valori_calcolati, filename, out_dir, tipo_grafico="versi"):
    nodi = {k: (float(v[0]), float(v[1])) for k, v in nodi.items()}   # solo per il disegno
    plt.figure(figsize=(10, 5), dpi=150) # Immagine più grande e definita
    
    # 1. Disegna le travi e gli assi monoassiali (ascissa s)
    for el in elementi:
        n1, n2 = nodi[int(el['n1'])], nodi[int(el['n2'])]
        plt.plot([n1[0], n2[0]], [n1[1], n2[1]], 'k-', linewidth=3.5, zorder=1)
        
        if tipo_grafico == "monoassiali":
            dx, dy = n2[0] - n1[0], n2[1] - n1[1]
            L = math.hypot(dx, dy)
            if L > 0:
                ux, uy = dx/L, dy/L
                plt.annotate('', xy=(n1[0] + ux*L*0.4, n1[1] + uy*L*0.4), xytext=(n1[0], n1[1]),
                             arrowprops=dict(edgecolor='purple', arrowstyle='->', lw=2.5), zorder=4)
                plt.text(n1[0] + ux*L*0.2, n1[1] + uy*L*0.2 + 0.15, 's', color='purple', weight='bold', fontsize=12)

    # 2. Disegna i nodi
    for nid, coords in nodi.items():
        plt.scatter(coords[0], coords[1], c='white', s=150, zorder=2, edgecolors='black', linewidth=2.5)
        plt.text(coords[0], coords[1], f"{nid}", fontsize=11, weight='bold', ha='center', va='center', color='black')

    # 3. Disegna le reazioni con riquadro bianco per massima visibilità
    scale = max([math.hypot(nodi[int(e['n2'])][0]-nodi[int(e['n1'])][0], nodi[int(e['n2'])][1]-nodi[int(e['n1'])][1]) for e in elementi]) * 0.18
    if scale == 0: scale = 0.5 
    
    for r in reazioni:
        xn, yn = nodi[r['nodo']]
        if tipo_grafico == "versi":
            label = f"{r['sym']}"
            colore_freccia = 'darkorange'
        else:
            valore = valori_calcolati.get(r['sym'], 0)
            try:
                num_val = float(valore)
                valore_arrotondato = round(num_val, 2)
                colore_freccia = 'red' if num_val < 0 else 'forestgreen'
            except Exception:
                valore_arrotondato = str(valore) 
                colore_freccia = 'forestgreen' 
            label = f"{r['sym']} = {valore_arrotondato}"

        # Stile box bianco dietro il testo per evitare sovrapposizioni illeggibili
        box_style = dict(boxstyle='round,pad=0.3', fc='white', ec=colore_freccia, alpha=0.9, lw=1.5)

        if r['tipo'] == 'Fx':
            plt.annotate(label, xy=(xn, yn), xytext=(xn - scale*1.2, yn), 
                         arrowprops=dict(facecolor=colore_freccia, width=2.5, headwidth=9), 
                         color=colore_freccia, weight='bold', fontsize=10, bbox=box_style, ha='right')
        elif r['tipo'] == 'Fy':
            plt.annotate(label, xy=(xn, yn), xytext=(xn, yn - scale*1.2), 
                         arrowprops=dict(facecolor=colore_freccia, width=2.5, headwidth=9), 
                         color=colore_freccia, weight='bold', fontsize=10, bbox=box_style, ha='center')
        elif r['tipo'] == 'M':
            plt.plot(xn, yn, marker='o', markersize=24, color=colore_freccia, zorder=3, fillstyle='none', markeredgewidth=2.5)
            plt.text(xn + 0.25, yn + 0.25, label, color=colore_freccia, weight='bold', fontsize=10, bbox=box_style)
        elif r['tipo'] == 'Fn':
            _ang = math.radians(float(r.get('angle', 0)) + 90)
            _dx, _dy = math.cos(_ang), math.sin(_ang)
            plt.annotate(label, xy=(xn, yn), xytext=(xn - scale*1.2*_dx, yn - scale*1.2*_dy),
                         arrowprops=dict(facecolor=colore_freccia, width=2.5, headwidth=9),
                         color=colore_freccia, weight='bold', fontsize=10, bbox=box_style, ha='center')
            if r.get('pendolo'):
                # biella dal nodo verso il vincolo a terra (direzione dell'asse), cerniera a terra con tratteggio
                xa, ya = xn + scale * 1.5 * _dx, yn + scale * 1.5 * _dy
                plt.plot([xn, xa], [yn, ya], color='dimgray', lw=3, zorder=1, solid_capstyle='round')
                plt.plot([xa], [ya], marker='o', markersize=9, markerfacecolor='white', markeredgecolor='dimgray', zorder=3)
                nx_, ny_ = -_dy, _dx
                for k_ in (-2, -1, 0, 1, 2):
                    bx_, by_ = xa + nx_ * 0.12 * scale * k_, ya + ny_ * 0.12 * scale * k_
                    plt.plot([bx_, bx_ + _dx * 0.18 * scale + nx_ * 0.05 * scale], [by_, by_ + _dy * 0.18 * scale + ny_ * 0.05 * scale],
                             color='dimgray', lw=1, zorder=1)

    plt.axis('equal'); plt.axis('off')
    filepath = os.path.join(out_dir, filename)
    plt.savefig(filepath, bbox_inches='tight', dpi=150)
    plt.close()
    return filepath.replace("\\", "/").lstrip("/")

def plot_diagrammi_singola_asta(L, expr_M, expr_T, expr_N, nome_asta, filename, out_dir, nome_sistema="S₀"):
    """Genera i 3 grafici (M, T, N) affiancati per una singola asta (Versione Ant-Crash Termico)."""
    fig, axes = plt.subplots(1, 3, figsize=(10, 2.5), dpi=120)
    fig.suptitle(f"Diagrammi Sistema {nome_sistema} - Asta {nome_asta}", fontsize=11, weight='bold')
    
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
    
    # Momento (rosso): M positivo = fibre tese a sinistra dell'asta (n1->n2). Stessa convenzione del diagramma globale:
    # il diagramma si disegna dal lato teso, quindi per un'asta orizzontale M<0 (fibre inferiori tese) va verso il basso.
    axes[0].plot(s_vals, M_vals, color='crimson', lw=2)
    axes[0].fill_between(s_vals, 0, M_vals, color='crimson', alpha=0.2)
    axes[0].axhline(0, color='black', lw=1)
    axes[0].set_title("Momento M₀(s)", fontsize=10)
    
    # Taglio (verde)
    axes[1].plot(s_vals, T_vals, color='forestgreen', lw=2)
    axes[1].fill_between(s_vals, 0, T_vals, color='forestgreen', alpha=0.2)
    axes[1].axhline(0, color='black', lw=1)
    axes[1].set_title("Taglio T₀(s)", fontsize=10)
    
    # Sforzo Normale (blu)
    axes[2].plot(s_vals, N_vals, color='dodgerblue', lw=2)
    axes[2].fill_between(s_vals, 0, N_vals, color='dodgerblue', alpha=0.2)
    axes[2].axhline(0, color='black', lw=1)
    axes[2].set_title("Sforzo Normale N₀(s)", fontsize=10)
    
    for ax in axes:
        ax.set_xlabel("ascissa s [m]", fontsize=8)
        ax.grid(True, linestyle='--', alpha=0.6)
        
    try:
        plt.tight_layout()
        filepath = os.path.join(out_dir, filename)
        plt.savefig(filepath, bbox_inches='tight')
    finally:
        plt.close(fig) # Garantisce la chiusura in memoria anche in caso di errore
        
    return filepath.replace("\\", "/").lstrip("/")

import matplotlib.patches as patches

def plot_schema_libero_asta(L, n1_id, n2_id, expr_M, expr_T, expr_N, nome_asta, filename, out_dir, dati_json):
    """Genera lo schema di corpo libero completo con vettori (N,T,M) e carichi esterni."""
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
    ax.axis('off')
    
    filepath = os.path.join(out_dir, filename)
    plt.savefig(filepath, bbox_inches='tight')
    plt.close()
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

    report.append("<h2>1. Dati della Struttura e Verifica Cinematica</h2>")
    report.append(f"<p>Telaio piano composto da <b>{num_aste} aste</b> e <b>{num_nodi} nodi</b>.</p>")
    
    # Stampa del box Proprietà Sezione
    report.append("<div style='background-color: #f8f9fa; padding: 10px; border-left: 4px solid #198754; margin-bottom: 15px;'>")
    report.append("<h4>Proprietà del Materiale e della Sezione</h4>")
    report.append("<ul>")
    report.append(f"<li>Modulo di Young (E): <b>{E:.2e} Pa</b></li>")
    report.append(f"<li>Dimensioni sezione rettangolare: b = <b>{b} m</b>, h = <b>{h} m</b></li>")
    report.append(f"<li>Area della sezione (A): b &middot; h = <b>{A:.4e} m<sup>2</sup></b></li>")
    report.append(f"<li>Momento d'Inerzia (I): (b &middot; h<sup>3</sup>) / 12 = <b>{I:.4e} m<sup>4</sup></b></li>")
    if EI_assegnato:
        report.append(f"<li>Rigidezza Flessionale (EI) <b>assegnata</b>: <b>{EI:.4e} N&middot;m<sup>2</sup></b></li>")
    else:
        report.append(f"<li>Rigidezza Flessionale (EI): E &middot; I = <b>{EI:.4e} N&middot;m<sup>2</sup></b></li>")
    report.append("</ul></div>")
    
    report.append("<p><b>Analisi dei vincoli esterni:</b></p>")
    report.append("<ul>")
    for det in dettagli_vincoli:
        report.append(det)
    report.append("</ul>")

    report.append("<div style='background-color: #e9ecef; padding: 15px; border-left: 4px solid #6c757d; margin-top: 10px; margin-bottom: 20px;'>")
    report.append(f"<p><b>Verifica Cinematica (Modello Multi-Corpo):</b></p>")
    report.append("<ul>")
    report.append(f"<li>Tronchi rigidi individuati (1 base + {corpi_aggiunti} aggiunti da {num_cerniere_interne} cerniere interne): <b>{num_tronchi}</b></li>")
    report.append(f"<li>Gradi di libertà liberi (3 &middot; tronchi): 3 &middot; {num_tronchi} = <b>{gdl_liberi}</b></li>")
    report.append(f"<li>Gradi di vincolo esterni a terra (v<sub>e</sub>): <b>{vincoli_reali}</b></li>")
    report.append(f"<li>Gradi di vincolo interni (v<sub>i</sub>): <b>{vincoli_interni}</b> (la cerniera blocca 2 traslazioni relative)</li>")
    report.append(f"<li>Grado di iperstaticità (h = v<sub>tot</sub> - gdl): {vincoli_totali} - {gdl_liberi} = <b>{grado_iperstaticita}</b></li>")
    report.append("</ul>")
    
    if grado_iperstaticita == 0:
        report.append("<p style='color:green;'>&rarr; La struttura è <b>isostatica</b>.</p>")
        report.append("</div>")
        
    elif grado_iperstaticita == 1:
        report.append(f"<p style='color:#0d6efd;'>&rarr; La struttura è <b>1 volta iperstatica</b>.</p>")
        report.append("<p><b>Conseguenza operativa:</b><br>")
        report.append("Le sole equazioni cardinali della statica non bastano. Si adotta il <i>Metodo delle Forze</i> con il Principio dei Lavori Virtuali, considerando il solo contributo flessionale.</p>")
        report.append("</div>")
        
    elif grado_iperstaticita > 1:
        report.append(f"<p style='color:darkorange;'>&rarr; La struttura è <b>{grado_iperstaticita} volte iperstatica</b>.</p>")
        report.append("</div>") # Chiude il box grigio principale
        
        # Genera un box giallo di allerta e ferma il solutore
        report.append("<div style='background-color: #fff3cd; padding: 15px; border-left: 5px solid #ffc107; margin-bottom: 20px;'>")
        report.append("<h4 style='color: #856404; margin-top: 0;'>⚠️ Limite del Solutore Superato</h4>")
        report.append(f"<p style='color: #856404;'>Questa versione del programma supporta un massimo di <b>1 grado di iperstaticità</b>.<br>")
        report.append(f"La struttura inserita richiede <b>{grado_iperstaticita} equazioni di congruenza simultanee</b>, funzionalità attualmente non implementata.</p>")
        report.append("</div>")
        return "\n".join(report), "", ""
        
    else:
        report.append(f"<p style='color:red;'>&rarr; La struttura è <b>labile</b> (h &lt; 0).</p>")
        report.append("</div>") # Chiude il box grigio principale
        
        # Genera un box rosso di errore e ferma il solutore
        report.append("<div style='background-color: #f8d7da; padding: 15px; border-left: 5px solid #dc3545; margin-bottom: 20px;'>")
        report.append("<h4 style='color: #721c24; margin-top: 0;'>🛑 Analisi Interrotta</h4>")
        report.append("<p style='color: #721c24;'>La struttura non è vincolata in modo sufficiente o efficace. È impossibile procedere con l'analisi statica.</p>")
        report.append("</div>")
        return "\n".join(report), "", ""

    # ==========================================
    # FASE 2: SCELTA INCOGNITA E SISTEMA PRINCIPALE
    # ==========================================
    reazioni_simboliche = []
    dip_ced = dip_term = False   # distorsioni che dipendono dalla sezione (usato dal dimensionamento)
    supports_n = dati_json.get('supports', [])

    # Simboli univoci per le reazioni dei vincoli inclinati (carrello inclinato, pendolo, bipendolo inclinato):
    # piu' pendoli sullo stesso nodo (doppio pendolo) devono avere incognite distinte.
    sym_Fn, _cont_Fn = {}, {}
    for i_s, s_ in enumerate(supports_n):
        if s_.get('rotated_roller'):
            nid_ = int(s_['node_id'])
            _cont_Fn[nid_] = _cont_Fn.get(nid_, 0) + 1
            sym_Fn[i_s] = sp.Symbol(f'Rn_{nid_}' if _cont_Fn[nid_] == 1 else f'Rn_{nid_}_{_cont_Fn[nid_]}')

    # Direzione dell'incognita iperstatica X1 (forza unitaria / coppia unitaria nel sistema S1)
    ang_X1 = sp.Integer(0)
    idx_X1_sup = None          # indice (in 'supports') del vincolo inclinato che porta X1 = Rn
    if grado_iperstaticita > 0:
        if tipo_svinc == 'Fn':
            cand = [i for i, s_ in enumerate(supports_n) if int(s_['node_id']) == nodo_svinc and s_.get('rotated_roller')]
            if not cand:
                _errore_cerniera("Incognita iperstatica non valida",
                                 f"Il Nodo {nodo_svinc} non ha un carrello inclinato / pendolo / bipendolo inclinato: "
                                 "scegli un'altra incognita.")
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
                                 f"Nel Nodo {nodo_svinc} non esiste la reazione {tipo_svinc} scelta come incognita "
                                 "(il vincolo e' stato modificato o cancellato?). Riseleziona l'incognita nel tab Tutor PLV.")
    X1_fx, X1_fy, X1_m = sp.Integer(0), sp.Integer(0), sp.Integer(0)
    if tipo_svinc == 'Fx': X1_fx = sp.Integer(1)
    elif tipo_svinc == 'Fy': X1_fy = sp.Integer(1)
    elif tipo_svinc in ('M', 'M_int'): X1_m = sp.Integer(1)
    elif tipo_svinc == 'Fn': X1_fx, X1_fy = cs_ang(ang_X1 + 90)

    def _reazioni_da_vincoli(rimuovi_X1):
        """Reazioni incognite dei vincoli a terra; con rimuovi_X1 si toglie quella scelta come incognita iperstatica."""
        out = []
        for i_s, s in enumerate(supports_n):
            nid = int(s.get('node_id'))
            if s.get('fix_u') and not (rimuovi_X1 and nid == nodo_svinc and tipo_svinc == 'Fx'):
                out.append({'nodo': nid, 'tipo': 'Fx', 'sym': sp.Symbol(f'H_{nid}')})
            if s.get('fix_v') and not (rimuovi_X1 and nid == nodo_svinc and tipo_svinc == 'Fy'):
                out.append({'nodo': nid, 'tipo': 'Fy', 'sym': sp.Symbol(f'V_{nid}')})
            if s.get('fix_w') and not (rimuovi_X1 and nid == nodo_svinc and tipo_svinc == 'M'):
                out.append({'nodo': nid, 'tipo': 'M', 'sym': sp.Symbol(f'M_{nid}')})
            if s.get('rotated_roller') and not (rimuovi_X1 and tipo_svinc == 'Fn' and i_s == idx_X1_sup):
                out.append({'nodo': nid, 'tipo': 'Fn', 'sym': sym_Fn[i_s], 'angle': s.get('angle', 0), 'sup_idx': i_s,
                            'pendolo': bool(s.get('pendolo'))})
        return out

    if grado_iperstaticita > 0:
        tipo_azione = {'Fx': "Forza Orizzontale", 'Fy': "Forza Verticale", 'Fn': "Forza Normale"}.get(tipo_svinc, "Momento")
        report.append("<h2>2. Scelta dell'incognita iperstatica e sistema principale</h2>")
        report.append(f"<p>Si assume come incognita iperstatica <b>X<sub>1</sub></b> = {tipo_azione.lower()} al Nodo {nodo_svinc} (X<sub>1</sub> = {tipo_svinc}<sub>{nodo_svinc}</sub>). Il sistema isostatico principale si ottiene rimuovendo il vincolo corrispondente.</p>")
        reazioni_simboliche = _reazioni_da_vincoli(True)
    else:
        report.append("<h2>2. Sistema Isostatico</h2>")
        report.append("<p>La struttura è <b>isostatica</b>. Non viene rimosso alcun vincolo e si procede al calcolo diretto delle reazioni e delle sollecitazioni.</p>")
        reazioni_simboliche = _reazioni_da_vincoli(False)   # manteniamo tutti i vincoli

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
        Somma_Fx_ext += Rx
        Somma_Fy_ext += Ry
        Somma_M_ext += Ry * (X_bar - xp) - Rx * (Y_bar - yp)

    variabili = [r['sym'] for r in reazioni_simboliche]
    
    def stop(titolo, msg):
        """Errore bloccante: aggiunge il box rosso al report ed interrompe il solutore."""
        report.append(_box_errore(titolo, msg))
        raise ErroreSolutore("\n".join(report))

    def risolvi_sicuro(equazioni, vars, nome="S0"):
        return risolvi_sistema_lineare(equazioni, vars, nome, report)

    # ==========================================
    # FASE 3: SISTEMA S0 (CARICHI REALI) E GRAFICI REAZIONI
    # ==========================================
    report.append("<h2>3. Sistema S<sub>0</sub> &mdash; carichi reali</h2>")
    report.append(f"<p>Si impone l'equilibrio globale (&Sigma;F<sub>x</sub>=0, &Sigma;F<sub>y</sub>=0, &Sigma;M=0).</p>")

    eq_x = sum([r['sym'] for r in reazioni_simboliche if r['tipo'] == 'Fx']) + Somma_Fx_ext
    eq_y = sum([r['sym'] for r in reazioni_simboliche if r['tipo'] == 'Fy']) + Somma_Fy_ext
    for r in reazioni_simboliche:
        if r['tipo'] == 'Fn':
            ang_rad = cs_ang(r['angle'] + 90)
            eq_x += r['sym'] * ang_rad[0]
            eq_y += r['sym'] * ang_rad[1]

    eq_m = sum([r['sym'] for r in reazioni_simboliche if r['tipo'] == 'M'])
    for r in reazioni_simboliche:
        xn, yn = nodi[r['nodo']]
        if r['tipo'] == 'Fx': eq_m -= r['sym'] * (yn - yp)
        if r['tipo'] == 'Fy': eq_m += r['sym'] * (xn - xp)
        if r['tipo'] == 'Fn':
            ang_rad = cs_ang(r['angle'] + 90)
            eq_m -= (r['sym'] * ang_rad[0]) * (yn - yp)
            eq_m += (r['sym'] * ang_rad[1]) * (xn - xp)
    eq_m += Somma_M_ext
    
    equazioni_S0 = [eq_x, eq_y, eq_m]
    
    for ci in cerniere_info:
      nc = ci['nodo']
      xc, yc = nodi[nc]
      for _id_asta_ramo, nodi_ramo in ci['rami']:
        eq_aux = 0
        for r in reazioni_simboliche:
            if _del_ramo(r['nodo'], nc, _id_asta_ramo, nodi_ramo, ci['attach']):
                xr, yr = nodi[r['nodo']]
                if r['tipo'] == 'Fx': eq_aux -= r['sym'] * (yr - yc)
                elif r['tipo'] == 'Fy': eq_aux += r['sym'] * (xr - xc)
                elif r['tipo'] == 'M': eq_aux += r['sym']
                elif r['tipo'] == 'Fn':
                    ang_rad = cs_ang(r['angle'] + 90)
                    eq_aux -= (r['sym'] * ang_rad[0]) * (yr - yc)
                    eq_aux += (r['sym'] * ang_rad[1]) * (xr - xc)

        for p in dati_json.get('point_loads', []):
            if _del_ramo(int(p['node_id']), nc, _id_asta_ramo, nodi_ramo, ci['attach']):
                xp_, yp_ = nodi[int(p['node_id'])]
                eq_aux += p.get('M', 0) + p.get('Fy', 0)*(xp_ - xc) - p.get('Fx', 0)*(yp_ - yc)

        for el_load in dati_json.get('element_loads', []):
            c = _carico(el_load)
            if c is None: continue
            n1_t, n2_t = int(c[0]['n1']), int(c[0]['n2'])
            if (n1_t in nodi_ramo or n1_t == nc) and (n2_t in nodi_ramo or n2_t == nc) and not (n1_t == nc and n2_t == nc):
                Rx, Ry, X_bar, Y_bar = _risultante(c)
                eq_aux += Ry * (X_bar - xc) - Rx * (Y_bar - yc)
        equazioni_S0.append(eq_aux)
        
    soluzione_S0 = risolvi_sicuro(equazioni_S0, variabili, "S0")

    img_s0_versi = plot_sistema_S0(nodi, elementi, reazioni_simboliche, soluzione_S0, "s0_versi.png", out_dir, "versi")
    img_s0_calc = plot_sistema_S0(nodi, elementi, reazioni_simboliche, soluzione_S0, "s0_calcolate.png", out_dir, "calcolate")
    img_s0_mono = plot_sistema_S0(nodi, elementi, reazioni_simboliche, soluzione_S0, "s0_monoassiali.png", out_dir, "monoassiali")

    report.append("<div style='text-align: center;'>")
    report.append(f"<h4>S<sub>0</sub> &mdash; reazioni (scelta dei versi)</h4><img src='file:///{img_s0_versi}' width='600'><br>")
    report.append(f"<h4>S<sub>0</sub> &mdash; reazioni calcolate</h4><img src='file:///{img_s0_calc}' width='600'><br>")
    report.append(f"<h4>S<sub>0</sub> &mdash; riferimenti monoassiali (ascissa s)</h4><img src='file:///{img_s0_mono}' width='600'>")
    report.append("</div>")

    report.append("<div style='background-color: #f8f9fa; padding: 15px; border-left: 4px solid #0d6efd; margin-bottom: 15px;'>")
    report.append(f"<h4>Equilibrio Globale S<sub>0</sub> (Polo Nodo {polo_id})</h4>")
    report.append("<ul>")
    report.append(f"<li>&Sigma;F<sub>x</sub> = 0 &rArr; <span class='formula'>{pulisci_formula(eq_x)} = 0</span></li>")
    report.append(f"<li>&Sigma;F<sub>y</sub> = 0 &rArr; <span class='formula'>{pulisci_formula(eq_y)} = 0</span></li>")
    report.append(f"<li>&Sigma;M<sub>{polo_id}</sub> = 0 &rArr; <span class='formula'>{pulisci_formula(eq_m)} = 0</span></li>")
    if len(equazioni_S0) > 3:
        report.append(f"<li><i>(+ {len(equazioni_S0)-3} equazioni ausiliarie per le cerniere interne)</i></li>")
    report.append("</ul></div>")

    # ==========================================
    # FASE 4: SISTEMA S1 (INCOGNITA UNITARIA)
    # ==========================================
    X1_valore_numerico = 0.0  # solo per la stampa
    X1_exact = sp.Integer(0)   # valore ESATTO usato in tutte le espressioni (isostatiche: 0)
    soluzione_S1 = {}
    
    if grado_iperstaticita > 0:
        if tipo_svinc == 'M_int':
            azione_S1_txt = "coppia di momenti interni unitari (flessione relativa)"
        else:
            azione_S1_txt = "coppia unitaria" if tipo_svinc == 'M' else f"forza {tipo_azione.lower().split()[1]} unitaria"
            
        report.append(f"<h2>4. Sistema S<sub>1</sub> &mdash; {azione_S1_txt} (X<sub>1</sub> = 1)</h2>")
        
        xs, ys = nodi[nodo_svinc]
        Somma_Fx_S1 = X1_fx
        Somma_Fy_S1 = X1_fy
        if tipo_svinc == 'M_int':
            Somma_M_S1 = 0
        else:
            Somma_M_S1 = X1_m + X1_fy * (xs - xp) - X1_fx * (ys - yp)

        eq_x_S1 = sum([r['sym'] for r in reazioni_simboliche if r['tipo'] == 'Fx']) + Somma_Fx_S1
        eq_y_S1 = sum([r['sym'] for r in reazioni_simboliche if r['tipo'] == 'Fy']) + Somma_Fy_S1
        for r in reazioni_simboliche:
            if r['tipo'] == 'Fn':
                ang_rad = cs_ang(r['angle'] + 90)
                eq_x_S1 += r['sym'] * ang_rad[0]
                eq_y_S1 += r['sym'] * ang_rad[1]

        eq_m_S1 = sum([r['sym'] for r in reazioni_simboliche if r['tipo'] == 'M'])
        for r in reazioni_simboliche:
            xn, yn = nodi[r['nodo']]
            if r['tipo'] == 'Fx': eq_m_S1 -= r['sym'] * (yn - yp)
            if r['tipo'] == 'Fy': eq_m_S1 += r['sym'] * (xn - xp)
            if r['tipo'] == 'Fn':
                ang_rad = cs_ang(r['angle'] + 90)
                eq_m_S1 -= (r['sym'] * ang_rad[0]) * (yn - yp)
                eq_m_S1 += (r['sym'] * ang_rad[1]) * (xn - xp)
        eq_m_S1 += Somma_M_S1

        equazioni_S1 = [eq_x_S1, eq_y_S1, eq_m_S1]
        
        for ci in cerniere_info:
          nc = ci['nodo']
          xc, yc = nodi[nc]
          for _id_asta_ramo, nodi_ramo in ci['rami']:
            eq_aux = 0
            for r in reazioni_simboliche:
                if _del_ramo(r['nodo'], nc, _id_asta_ramo, nodi_ramo, ci['attach']):
                    xr, yr = nodi[r['nodo']]
                    if r['tipo'] == 'Fx': eq_aux -= r['sym'] * (yr - yc)
                    elif r['tipo'] == 'Fy': eq_aux += r['sym'] * (xr - xc)
                    elif r['tipo'] == 'M': eq_aux += r['sym']
                    elif r['tipo'] == 'Fn':
                        ang_rad = cs_ang(r['angle'] + 90)
                        eq_aux -= (r['sym'] * ang_rad[0]) * (yr - yc)
                        eq_aux += (r['sym'] * ang_rad[1]) * (xr - xc)
                    
            if _del_ramo(nodo_svinc, nc, _id_asta_ramo, nodi_ramo, ci['attach']):
                xs_v, ys_v = nodi[nodo_svinc]
                eq_aux += (X1_m if tipo_svinc == 'M' else 0) + X1_fy * (xs_v - xc) - X1_fx * (ys_v - yc)
                
            if tipo_svinc == 'M_int' and nodo_svinc == nc:
                eq_aux += 1
                
            equazioni_S1.append(eq_aux)

        # --- Controllo: il sistema isostatico principale deve essere staticamente determinato ---
        try:
            A_chk, _b_chk = sp.linear_eq_to_matrix(equazioni_S1, variabili)
            if A_chk.rank() < len(variabili):
                report.append("<div style='background-color: #f8d7da; padding: 15px; border-left: 5px solid #dc3545; margin-bottom: 20px;'>")
                report.append("<h4 style='color: #721c24; margin-top: 0;'>🛑 Incognita iperstatica non valida</h4>")
                report.append("<p style='color: #721c24;'>Rimuovendo questo vincolo il <b>sistema principale risulta labile</b> (le reazioni restanti non sono determinabili dall'equilibrio, es. tre vincoli concorrenti in un punto o paralleli). Scegli un'altra incognita iperstatica.</p></div>")
                return "\n".join(report), "", ""
        except Exception:
            pass
        soluzione_S1 = risolvi_sicuro(equazioni_S1, variabili, "S1")

        X1_sym = sp.Symbol('X_1')
        reazioni_simboliche_S1 = reazioni_simboliche.copy()
        reazioni_simboliche_S1.append({'nodo': nodo_svinc, 'tipo': tipo_svinc, 'sym': X1_sym, 'angle': ang_X1})
        
        soluzione_plot_S1 = soluzione_S1.copy()
        soluzione_plot_S1[X1_sym] = 1.0 

        img_s1_versi = plot_sistema_S0(nodi, elementi, reazioni_simboliche_S1, soluzione_plot_S1, "s1_versi.png", out_dir, "versi")
        img_s1_calc = plot_sistema_S0(nodi, elementi, reazioni_simboliche_S1, soluzione_plot_S1, "s1_calcolate.png", out_dir, "calcolate")
        img_s1_mono = plot_sistema_S0(nodi, elementi, reazioni_simboliche_S1, soluzione_plot_S1, "s1_monoassiali.png", out_dir, "monoassiali")

        report.append("<div style='text-align: center; margin-bottom: 25px;'>")
        report.append(f"<h4>S<sub>1</sub> &mdash; forze applicate e reazioni incognite</h4><img src='file:///{img_s1_versi}' width='600'><br>")
        report.append(f"<h4>S<sub>1</sub> &mdash; reazioni calcolate (equilibrio con X<sub>1</sub> = 1)</h4><img src='file:///{img_s1_calc}' width='600'><br>")
        report.append(f"<h4>S<sub>1</sub> &mdash; riferimenti monoassiali (ascissa s)</h4><img src='file:///{img_s1_mono}' width='600'>")
        report.append("</div>")

        report.append("<div style='background-color: #f8f9fa; padding: 15px; border-left: 4px solid #dc3545; margin-bottom: 20px; border-radius: 4px;'>")
        report.append(f"<h4>Equilibrio Globale S<sub>1</sub> (Polo Nodo {polo_id})</h4>")
        report.append("<ul>")
        report.append(f"<li>&Sigma;F<sub>x</sub> = 0 &rArr; <span class='formula'>{pulisci_formula(eq_x_S1)} = 0</span></li>")
        report.append(f"<li>&Sigma;F<sub>y</sub> = 0 &rArr; <span class='formula'>{pulisci_formula(eq_y_S1)} = 0</span></li>")
        report.append(f"<li>&Sigma;M<sub>{polo_id}</sub> = 0 &rArr; <span class='formula'>{pulisci_formula(eq_m_S1)} = 0</span></li>")
        report.append("</ul></div>")

    # ==========================================
    # FASE 5: SOLLECITAZIONI E RIFERIMENTI MONOASSIALI
    # ==========================================
    report.append("<h2>5. Caratteristiche della sollecitazione (Studio Analitico M, T, N)</h2>")
    
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
        nid = r['nodo']
        tipo_r = r['tipo']
        val_S0 = soluzione_S0.get(r['sym'], 0)
        
        if tipo_r == 'Fx': 
            Forze_S0[nid]['Fx'] += val_S0
        elif tipo_r == 'Fy': 
            Forze_S0[nid]['Fy'] += val_S0
        elif tipo_r == 'M': 
            Forze_S0[nid]['M'] += val_S0
        elif tipo_r == 'Fn':
            ang_rad = cs_ang(r['angle'] + 90)
            Forze_S0[nid]['Fx'] += val_S0 * ang_rad[0]
            Forze_S0[nid]['Fy'] += val_S0 * ang_rad[1]
            
        if grado_iperstaticita > 0:
            val_S1 = soluzione_S1.get(r['sym'], 0)
            if tipo_r == 'Fx': 
                Forze_S1[nid]['Fx'] += val_S1
            elif tipo_r == 'Fy': 
                Forze_S1[nid]['Fy'] += val_S1
            elif tipo_r == 'M': 
                Forze_S1[nid]['M'] += val_S1
            elif tipo_r == 'Fn':
                ang_rad = cs_ang(r['angle'] + 90)
                Forze_S1[nid]['Fx'] += val_S1 * ang_rad[0]
                Forze_S1[nid]['Fy'] += val_S1 * ang_rad[1]

    M0_dict, M1_dict, T0_dict, N0_dict, T1_dict, N1_dict = {}, {}, {}, {}, {}, {}
    
    def calcola_sollecitazioni(sistema_forze, el, is_S0=True):
        n1, n2 = int(el['n1']), int(el['n2'])
        x1, y1 = nodi[n1]
        x2, y2 = nodi[n2]
        L = hyp(x2 - x1, y2 - y1)
        s = sp.Symbol("s")
        if L == 0: return sp.S(0), sp.S(0), sp.S(0)
            
        cos_alpha = (x2 - x1) / L
        sin_alpha = (y2 - y1) / L
        
        visited = set()
        stack = [n1]
        while stack:
            curr = stack.pop()
            if curr not in visited:
                visited.add(curr)
                for neighbor in adj[curr]:
                    if neighbor != n2 and neighbor not in visited: stack.append(neighbor)
                        
        xc = x1 + cos_alpha * s
        yc = y1 + sin_alpha * s
        
        M_s, Fx_sum, Fy_sum = 0, 0, 0
        for v in visited:
            fx, fy, m = sistema_forze[v]['Fx'], sistema_forze[v]['Fy'], sistema_forze[v]['M']
            Fx_sum += fx
            Fy_sum += fy
            M_s += m + fy * (nodi[v][0] - xc) - fx * (nodi[v][1] - yc)
            
        if is_S0:
            t_ = sp.Symbol('t')
            for el_load in dati_json.get('element_loads', []):
                c = _carico(el_load)
                if c is None: continue
                target_el, L_t, ux_t, uy_t, wx, wy, tipo_q = c
                if target_el['id'] == el['id']:
                    # carico sulla stessa asta: porzione tra 0 e s (contributo a M, T e N)
                    k_ = ux_t * wy - uy_t * wx
                    phi = _phi(tipo_q, t_, L_t)
                    Phi = sp.integrate(phi, (t_, 0, s))
                    M_s += k_ * sp.integrate((t_ - s) * phi, (t_, 0, s))
                    Fx_sum += wx * Phi
                    Fy_sum += wy * Phi
                elif int(target_el['n1']) in visited and int(target_el['n2']) in visited:
                    Rx, Ry, X_c, Y_c = _risultante(c)
                    Fx_sum += Rx
                    Fy_sum += Ry
                    M_s += Ry * (X_c - xc) - Rx * (Y_c - yc)

        T_s = -sp.diff(M_s, s)
        N_s = -Fx_sum * cos_alpha - Fy_sum * sin_alpha
        return pulisci_espr(sp.simplify(M_s)), pulisci_espr(sp.simplify(T_s)), pulisci_espr(sp.simplify(N_s))

    report.append("<h3 style='color: #198754; margin-top: 30px; border-bottom: 2px solid #198754;'>5.1 Sistema Principale (S₀) - Carichi Reali</h3>")
    for el in elementi:
        L = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
        
        M0, T0, N0 = calcola_sollecitazioni(Forze_S0, el, True)
        if grado_iperstaticita > 0:
            M1, T1, N1 = calcola_sollecitazioni(Forze_S1, el, False)
        else:
            M1, T1, N1 = sp.S(0), sp.S(0), sp.S(0)
            
        M0_dict[el['id']] = M0; M1_dict[el['id']] = M1
        T0_dict[el['id']] = T0; N0_dict[el['id']] = N0
        T1_dict[el['id']] = T1; N1_dict[el['id']] = N1
        
        file_libero = f"s0_libero_{el['id']}.png"
        q_perp_p, q_ax_p = 0, 0
        for _el_load in dati_json.get('element_loads', []):
            _c = _carico(_el_load)
            if _c is not None and _c[0]['id'] == el['id']:
                q_perp_p += _c[4] * _c[3] - _c[5] * _c[2]     # componente trasversale (positiva verso la normale destra)
                q_ax_p += _c[4] * _c[2] + _c[5] * _c[3]       # componente assiale
        dati_plot = {'element_loads': [{'elem_id': el['id'], 'q_trans': q_perp_p}], 'q_ax': q_ax_p}
        img_libero = plot_schema_libero_asta(L, el['n1'], el['n2'], M0, T0, N0, el['id'], file_libero, out_dir, dati_plot)
        file_asta = f"s0_asta_{el['id']}.png"
        img_asta = plot_diagrammi_singola_asta(L, M0, T0, N0, el['id'], file_asta, out_dir, "S₀")
        
        report.append(f"<div style='background-color: #ffffff; padding: 15px; border: 1px solid #ced4da; margin-bottom: 25px; border-radius: 8px; box-shadow: 2px 2px 8px rgba(0,0,0,0.05);'>")
        report.append(f"<h4 style='color: #198754; margin-top:0; border-bottom: 1px solid #dee2e6; padding-bottom: 5px;'>Asta {el['id']} (Nodo {el['n1']} &rarr; Nodo {el['n2']})</h4>")
        report.append(f"<div style='text-align: center; margin: 15px 0;'><img src='file:///{img_libero}' style='max-width: 90%;'></div>")
        report.append(f"<p style='color: #6c757d; font-size: 13px;'><i>Dominio asse locale:</i> 0 &le; s &le; {float(L):.2f} m</p>")
        report.append("<table style='width: 100%; margin-bottom: 15px;'><tr>")
        report.append("<th style='width: 33%; background-color:#e3f2fd; border-bottom: 2px solid #0d6efd;'>Sforzo Normale N<sub>0</sub>(s)</th>")
        report.append("<th style='width: 33%; background-color:#e8f5e9; border-bottom: 2px solid #198754;'>Taglio T<sub>0</sub>(s)</th>")
        report.append("<th style='width: 33%; background-color:#ffebee; border-bottom: 2px solid #dc3545;'>Momento Flettente M<sub>0</sub>(s)</th></tr><tr>")
        report.append(f"<td><span class='formula'>{pulisci_formula(N0)}</span></td>")
        report.append(f"<td><span class='formula'>{pulisci_formula(T0)}</span></td>")
        report.append(f"<td><span class='formula'>{pulisci_formula(M0)}</span></td></tr></table>")
        report.append(f"<div style='text-align: center; margin-top: 15px;'><img src='file:///{img_asta}' style='max-width: 100%; border: 1px solid #dee2e6; border-radius: 4px;'></div></div>")

    if grado_iperstaticita > 0:
        report.append("<h3 style='color: #dc3545; margin-top: 30px; border-bottom: 2px solid #dc3545;'>5.2 Sistema Unitario (S₁) - Incognita Virtuale Unitaria</h3>")
        for el in elementi:
            L = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
            M1, T1, N1 = M1_dict[el['id']], T1_dict[el['id']], N1_dict[el['id']]
            
            file_libero_s1 = f"s1_libero_{el['id']}.png"
            img_libero_s1 = plot_schema_libero_asta(L, el['n1'], el['n2'], M1, T1, N1, el['id'], file_libero_s1, out_dir, {})
            file_asta_s1 = f"s1_asta_{el['id']}.png"
            img_asta_s1 = plot_diagrammi_singola_asta(L, M1, T1, N1, el['id'], file_asta_s1, out_dir, "S₁")
            
            report.append(f"<div style='background-color: #ffffff; padding: 15px; border: 1px solid #ced4da; margin-bottom: 25px; border-radius: 8px; box-shadow: 2px 2px 8px rgba(0,0,0,0.05);'>")
            report.append(f"<h4 style='color: #dc3545; margin-top:0; border-bottom: 1px solid #dee2e6; padding-bottom: 5px;'>Asta {el['id']} (Nodo {el['n1']} &rarr; Nodo {el['n2']})</h4>")
            report.append(f"<div style='text-align: center; margin: 15px 0;'><img src='file:///{img_libero_s1}' style='max-width: 90%;'></div>")
            report.append("<table style='width: 100%; margin-bottom: 15px;'><tr>")
            report.append("<th style='width: 33%; background-color:#e3f2fd; border-bottom: 2px solid #0d6efd;'>Sforzo Normale N<sub>1</sub>(s)</th>")
            report.append("<th style='width: 33%; background-color:#e8f5e9; border-bottom: 2px solid #198754;'>Taglio T<sub>1</sub>(s)</th>")
            report.append("<th style='width: 33%; background-color:#ffebee; border-bottom: 2px solid #dc3545;'>Momento Flettente M<sub>1</sub>(s)</th></tr><tr>")
            report.append(f"<td><span class='formula'>{pulisci_formula(N1)}</span></td>")
            report.append(f"<td><span class='formula'>{pulisci_formula(T1)}</span></td>")
            report.append(f"<td><span class='formula'>{pulisci_formula(M1)}</span></td></tr></table>")
            report.append(f"<div style='text-align: center; margin-top: 15px;'><img src='file:///{img_asta_s1}' style='max-width: 100%; border: 1px solid #dee2e6; border-radius: 4px;'></div></div>")

    # ==========================================
    # FASE 6 E 7: PLV E CONGRUENZA
    # ==========================================
    if grado_iperstaticita > 0:
        report.append("<h2>6. Coefficienti di deformabilità, impostazione degli integrali (PLV)</h2>")
        report.append("<div style='background-color: #e2e3e5; padding: 10px; border-left: 4px solid #6c757d; margin-bottom: 15px;'>")
        report.append("<b>Formula generale Integrali di Mohr (solo flessione):</b><br>")
        report.append("<span class='formula'>&delta;<sub>ij</sub> = &int; [M<sub>i</sub>(s) &middot; M<sub>j</sub>(s) / EI] ds</span>")
        report.append("</div>")
        
        _mat_d, _sec_d0 = dati_json.get('material', {}), dati_json.get('section', {})
        if 'EI' in _mat_d:
            EI_val = _mat_d['EI']
        else:
            EI_val = (_mat_d.get('E', 2.1e11) * (_sec_d0.get('b', 0.2) * _sec_d0.get('h', 0.4)**3) / 12)
        delta_10, delta_11 = 0, 0
        
        report.append("<div style='background-color: #f8f9fa; padding: 15px; border-left: 4px solid #6c757d; margin-bottom: 20px; border-radius: 4px;'>")
        report.append("<h4>Impostazione degli Integrali di Mohr tratto per tratto:</h4><ul>")
        for el in elementi:
            x_sym = sp.Symbol("s")
            L = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
            M0_expr, M1_expr = M0_dict[el['id']], M1_dict[el['id']]
            
            d10_tratto = sp.integrate(M0_expr * M1_expr, (x_sym, 0, L))
            d11_tratto = sp.integrate(M1_expr * M1_expr, (x_sym, 0, L))
            
            report.append(f"<li style='margin-bottom: 15px;'><b>Asta {el['id']} (L = {float(L):.2f} m):</b><br>")
            report.append(f"&delta;<sub>10</sub> = (1/EI) &int;<sub>0</sub><sup>{float(L):.2f}</sup> [ {pulisci_formula(M0_expr)} ] &middot; [ {pulisci_formula(M1_expr)} ] ds = <span class='formula'>{pulisci_formula(d10_tratto)} / EI</span><br>")
            report.append(f"&delta;<sub>11</sub> = (1/EI) &int;<sub>0</sub><sup>{float(L):.2f}</sup> [ {pulisci_formula(M1_expr)} ]<sup>2</sup> ds = <span class='formula'>{pulisci_formula(d11_tratto)} / EI</span></li>")
            
            delta_10 += d10_tratto
            delta_11 += d11_tratto

        report.append("</ul></div>")
        delta_10 = sp.simplify(delta_10); delta_11 = sp.simplify(delta_11)
        if delta_11 == 0 or abs(float(delta_11)) < 1e-14 * max(1.0, float(sum(L_ ** 3 for L_ in [hyp(nodi[int(e_['n2'])][0] - nodi[int(e_['n1'])][0], nodi[int(e_['n2'])][1] - nodi[int(e_['n1'])][1]) for e_ in elementi]))):
            stop("Calcolo impossibile: &delta;<sub>11</sub> = 0",
                 "L'incognita iperstatica scelta non produce momento flettente nel sistema S<sub>1</sub> (M<sub>1</sub> &equiv; 0): "
                 "l'equazione di congruenza non può determinare X<sub>1</sub>. Scegli un'altra incognita iperstatica. "
                 "Nessun risultato numerico viene mostrato per evitare valori falsi.")
        report.append(f"<p>Sommando i contributi flessionali di tutti i tratti si ottiene:</p>")
        report.append(f"<div style='padding: 15px; border: 1px dashed #6c757d; display: inline-block; margin-bottom: 20px; background-color: #ffffff;'>")
        report.append(f"<b>&delta;<sub>10</sub></b> = <span class='formula'>{pulisci_formula(sp.simplify(delta_10))} / EI</span><br><br>")
        report.append(f"<b>&delta;<sub>11</sub></b> = <span class='formula'>{pulisci_formula(sp.simplify(delta_11))} / EI</span>")
        report.append("</div>")

        # ------------------------------------------------------------------
        # DISTORSIONI: variazioni termiche, cedimenti vincolari, rotazioni imposte
        #   eta = (d10 + X1*d11)/EI + Dtermico + Dassiale - sum(R1_i * delta_i)
        # ------------------------------------------------------------------
        alpha_mat = ex(dati_json.get('material', {}).get('alpha', 0))
        h_sez_T = ex(dati_json.get('section', {}).get('h', 0.4))
        s_sym_d = sp.Symbol("s")
        termine_T, termine_N = sp.Integer(0), sp.Integer(0)
        KT_sum = sp.Integer(0)      # sum(alpha*(T+ - T-)*int M1 ds): il termine termico flessionale vale KT_sum/h (serve al dimensionamento)
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
            righe_term.append(f"Asta {eid}: &kappa;<sub>T</sub> = &alpha;(T<sub>+</sub>&minus;T<sub>&minus;</sub>)/h = {float(kappa_T):.4e} 1/m &rarr; &int;M<sub>1</sub>&kappa;<sub>T</sub>ds = {float(kappa_T*int_M1):.4e}"
                              + (f"; &epsilon;<sub>T</sub> = {float(eps_T):.3e} &rarr; &int;N<sub>1</sub>&epsilon;<sub>T</sub>ds = {float(eps_T*int_N1):.4e}" if abs(eps_T) > 0 else ""))

        sup_map = {}
        for s_ in dati_json.get('supports', []):
            sup_map.setdefault(int(s_['node_id']), {}).update(s_)
        chiavi_delta = {'Fx': 'delta_u', 'Fy': 'delta_v', 'M': 'delta_w', 'Fn': 'delta_n'}
        somma_Rd = sp.Integer(0)
        righe_ced = []
        for r_ in reazioni_simboliche:
            if r_['tipo'] == 'Fn' and r_.get('sup_idx') is not None:
                # piu' bielle/vincoli inclinati sullo stesso nodo: il cedimento appartiene a QUELLO specifico vincolo
                dval = ex(supports_n[r_['sup_idx']].get('delta_n', 0))
            else:
                dval = ex(sup_map.get(r_['nodo'], {}).get(chiavi_delta[r_['tipo']], 0))
            if dval != 0:
                R1v = ex(soluzione_S1.get(r_['sym'], 0))
                somma_Rd += R1v * dval
                righe_ced.append(f"Nodo {r_['nodo']} ({r_['tipo']}): R<sub>1</sub> = {float(R1v):.4f}, &delta; imposto = {float(dval):.4e} &rarr; R<sub>1</sub>&middot;&delta; = {float(R1v*dval):.4e}")
        eta_X1 = sp.Integer(0)
        if tipo_svinc == 'Fn' and idx_X1_sup is not None:
            eta_X1 = ex(supports_n[idx_X1_sup].get('delta_n', 0))
        elif tipo_svinc in chiavi_delta:
            eta_X1 = ex(sup_map.get(nodo_svinc, {}).get(chiavi_delta[tipo_svinc], 0))
        ha_distorsioni = bool(righe_term or righe_ced or eta_X1 != 0)
        dip_ced = bool(righe_ced or eta_X1 != 0)
        dip_term = bool(righe_term)
        termine_noto_dist = EI_val * (eta_X1 - termine_T - termine_N + somma_Rd)
        if ha_distorsioni:
            report.append("<h2>6-bis. Distorsioni: variazioni termiche, cedimenti vincolari, rotazioni imposte</h2>")
            report.append("<div style='background-color:#fff3cd; padding:12px; border-left:4px solid #ffc107; margin-bottom:15px;'>")
            report.append("<b>PLV con distorsioni (sistema S<sub>1</sub> sulle deformazioni reali):</b><br>")
            report.append("<span class='formula'>&eta;<sub>1</sub> = (&delta;<sub>10</sub> + X<sub>1</sub>&delta;<sub>11</sub>)/EI + &int;M<sub>1</sub>&kappa;<sub>T</sub>ds + &int;N<sub>1</sub>&epsilon;<sub>T</sub>ds &minus; &Sigma;R<sub>1,i</sub>&middot;&delta;<sub>i</sub></span><br>")
            report.append("dove &eta;<sub>1</sub> &egrave; lo spostamento imposto nel punto/direzione di X<sub>1</sub> (cedimento del vincolo rimosso, altrimenti 0), "
                          "R<sub>1,i</sub> sono le reazioni di S<sub>1</sub> sui vincoli che cedono, &delta;<sub>i</sub> i cedimenti/rotazioni imposti.<ul>")
            for r_txt in righe_term: report.append(f"<li>{r_txt}</li>")
            for r_txt in righe_ced: report.append(f"<li>{r_txt}</li>")
            if eta_X1 != 0: report.append(f"<li>Vincolo rimosso (X<sub>1</sub>) con cedimento imposto &eta;<sub>1</sub> = {float(eta_X1):.4e}</li>")
            report.append("</ul></div>")

        report.append("<h2>7. Equazione di Congruenza e Calcolo Incognita Iperstatica</h2>")
        report.append("<div style='background-color: #e2e3e5; padding: 10px; border-left: 4px solid #dc3545; margin-bottom: 15px;'>")
        report.append("<b>Equazione di congruenza e formula risolutiva:</b><br>")
        if ha_distorsioni:
            report.append("<span class='formula'>&delta;<sub>10</sub> + X<sub>1</sub>&delta;<sub>11</sub> = EI&middot;(&eta;<sub>1</sub> &minus; &int;M<sub>1</sub>&kappa;<sub>T</sub>ds &minus; &int;N<sub>1</sub>&epsilon;<sub>T</sub>ds + &Sigma;R<sub>1,i</sub>&delta;<sub>i</sub>)</span><br>")
            report.append(f"<span class='formula'>EI&middot;(&hellip;) = {float(termine_noto_dist):.4f} N&middot;m&sup3;-equiv. &rArr; X<sub>1</sub> = [EI&middot;(&hellip;) &minus; &delta;<sub>10</sub>] / &delta;<sub>11</sub></span></div>")
        else:
            report.append("<span class='formula'>&delta;<sub>10</sub> + X<sub>1</sub> &middot; &delta;<sub>11</sub> = 0  &rArr;  X<sub>1</sub> = - &delta;<sub>10</sub> / &delta;<sub>11</sub></span></div>")
        
        X1_sym = sp.Symbol('X1')
        try:
            soluzione_X1 = sp.solve(sp.Eq(delta_10 + X1_sym * delta_11, termine_noto_dist), X1_sym)
        except Exception as e:
            stop("Errore nell'equazione di congruenza", f"{e}")
        if not soluzione_X1:
            stop("Equazione di congruenza senza soluzione", "&delta;<sub>10</sub> + X<sub>1</sub>&delta;<sub>11</sub> = termine noto non ammette soluzione.")
        X1_exact = pulisci_espr(sp.simplify(soluzione_X1[0]))
        X1_valore_numerico = float(X1_exact)
        report.append(f"<div style='background-color:#d1e7dd; color:#0f5132; padding:15px; border-left: 5px solid #198754; border-radius: 4px; margin-top: 15px;'>")
        report.append(f"<b>Soluzione dell'equazione (Valore dell'incognita iperstatica):</b><br>")
        report.append(f"<span style='font-size:20px;'><b>X<sub>1</sub> = {X1_valore_numerico:.3f}</b></span>"
                      + (f" &nbsp;(esatto: <span class='formula'>{pulisci_formula(X1_exact)}</span>)" if not X1_exact.is_Integer else "") + "</div>")
    else:
        report.append("<h2>6. Integrazione PLV non necessaria</h2>")
        report.append("<p>La struttura è <b>isostatica</b>. Non ci sono incognite iperstatiche da risolvere con il Principio dei Lavori Virtuali.</p>")
        report.append("<h2>7. Congruenza non necessaria</h2>")
        report.append("<p>Le sollecitazioni calcolate nel sistema reale S<sub>0</sub> sono definitive e verranno usate direttamente per il dimensionamento e la verifica.</p>")

# ==========================================
# FASE 8: SOLLECITAZIONI FINALI TRATTO PER TRATTO (Sovrapposizione)
# ==========================================
    if X1_valore_numerico is not None:
        report.append("<h2>8. Caratteristiche finali della sollecitazione per tratto (M, T, N reali)</h2>")
        report.append(f"<p>Avendo ricavato <b>X<sub>1</sub> = {X1_valore_numerico:.3f}</b>, calcoliamo le sollecitazioni reali definitive sovrapponendo gli effetti:<br>")
        report.append("M(s) = M<sub>0</sub>(s) + X<sub>1</sub> &middot; M<sub>1</sub>(s)<br>")
        report.append("T(s) = T<sub>0</sub>(s) + X<sub>1</sub> &middot; T<sub>1</sub>(s)<br>")
        report.append("N(s) = N<sub>0</sub>(s) + X<sub>1</sub> &middot; N<sub>1</sub>(s)</p>")
        
        for el in elementi:
            el_id = el['id']
            L = hyp(nodi[int(el['n2'])][0] - nodi[int(el['n1'])][0], nodi[int(el['n2'])][1] - nodi[int(el['n1'])][1])
            
            # Calcolo delle funzioni reali finali 
            M_fin = sp.simplify(M0_dict[el_id] + X1_exact * M1_dict[el_id])
            T_fin = sp.simplify(T0_dict[el_id] + X1_exact * T1_dict[el_id])
            N_fin = sp.simplify(N0_dict[el_id] + X1_exact * N1_dict[el_id])
            
            report.append(f"<div style='background-color: #f8f9fa; padding: 15px; border: 1px solid #dee2e6; border-left: 4px solid #0d6efd; margin-bottom: 12px; border-radius: 4px;'>")
            report.append(f"<h4 style='color: #0d6efd; margin: 0 0 10px 0;'>Asta {el_id} (Nodo {el['n1']} &rarr; Nodo {el['n2']})</h4>")
            report.append("<table style='width: 100%;'><tr>")
            report.append("<th style='width: 33%; background-color:#e3f2fd; border-bottom: 2px solid #0d6efd;'>N(s)</th>")
            report.append("<th style='width: 33%; background-color:#e8f5e9; border-bottom: 2px solid #198754;'>T(s)</th>")
            report.append("<th style='width: 33%; background-color:#ffebee; border-bottom: 2px solid #dc3545;'>M(s)</th></tr><tr>")
            report.append(f"<td><span class='formula'>{pulisci_formula(N_fin)}</span></td>")
            report.append(f"<td><span class='formula'>{pulisci_formula(T_fin)}</span></td>")
            report.append(f"<td><span class='formula'>{pulisci_formula(M_fin)}</span></td></tr></table></div>")

        # ==========================================
        # FASE 8-bis: REAZIONI FINALI, FORZE NELLE CERNIERE, VERIFICA DI EQUILIBRIO
        # ==========================================
        report.append("<h2>8-bis. Reazioni vincolari finali e forze nelle cerniere interne</h2>")


        # --- Reazioni finali: R = R0 + X1 * R1 (piu' la stessa X1, che e' la reazione del vincolo rimosso) ---
        nome_sym = {'Fx': 'H', 'Fy': 'V', 'M': 'M', 'Fn': 'Rn'}
        reaz_fin = []
        def _nome_reaz(tipo_, nodo_, sym_=None):
            parti_ = str(sym_).split('_') if sym_ is not None else []
            suff = f"{nodo_}.{parti_[2]}" if len(parti_) == 3 else f"{nodo_}"     # 2o/3o pendolo sullo stesso nodo
            return f"{nome_sym[tipo_]}<sub>{suff}</sub>"
        for r in reazioni_simboliche:
            v = soluzione_S0.get(r['sym'], 0)
            if grado_iperstaticita > 0:
                v = v + X1_exact * soluzione_S1.get(r['sym'], 0)
            reaz_fin.append({'nodo': r['nodo'], 'tipo': r['tipo'], 'angle': r.get('angle', 0),
                             'sup_idx': r.get('sup_idx'), 'pendolo': r.get('pendolo', False),
                             'nome': _nome_reaz(r['tipo'], r['nodo'], r['sym'] if r['tipo'] == 'Fn' else None),
                             'val': pulisci_espr(sp.simplify(v))})
        if grado_iperstaticita > 0 and tipo_svinc in nome_sym:
            reaz_fin.append({'nodo': nodo_svinc, 'tipo': tipo_svinc, 'angle': ex(ang_X1),
                             'sup_idx': idx_X1_sup,
                             'pendolo': bool(tipo_svinc == 'Fn' and supports_n[idx_X1_sup].get('pendolo')),
                             'nome': (_nome_reaz(tipo_svinc, nodo_svinc, sym_Fn.get(idx_X1_sup) if tipo_svinc == 'Fn' else None)
                                      + " (= X<sub>1</sub>)"),
                             'val': X1_exact})

        report.append("<p>Convenzione: H positiva verso +x, V positiva verso +y, M positivo antiorario; "
                      "R<sub>n</sub> positiva lungo la normale del carrello (angolo + 90&deg;); per pendoli e doppi pendoli (bipendoli) inclinati "
                      "R<sub>n</sub> &egrave; lungo l'asse delle bielle, orientato verso il vincolo a terra (R<sub>n</sub> &gt; 0 = trazione). "
                      "Valori ottenuti per sovrapposizione R = R<sub>0</sub> + X<sub>1</sub>&middot;R<sub>1</sub>.</p>")
        report.append("<table><tr><th>Reazione</th><th>Nodo</th><th>Valore</th><th>Unit&agrave;</th></tr>")
        for rf in reaz_fin:
            unita = "N&middot;m" if rf['tipo'] == 'M' else "N"
            report.append(f"<tr><td>{rf['nome']}</td><td>{rf['nodo']}</td><td><b>{float(rf['val']):.2f}</b></td><td>{unita}</td></tr>")
        report.append("</table>")

        # --- Verifica di equilibrio globale con le reazioni finali ---
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
            res_x += fx_
            res_y += fy_
            res_m += m_ + fy_ * (xn_ - xp) - fx_ * (yn_ - yp)
        residui = [float(sp.N(res_x)), float(sp.N(res_y)), float(sp.N(res_m))]
        equilibrio_ok = all(abs(v) < 1e-6 for v in residui)

        # --- Forze trasmesse dalle cerniere interne ---
        righe_cern = []
        cerniere_ok = True
        for ci in cerniere_info:
          nc = ci['nodo']
          xc, yc = nodi[nc]
          for id_asta_ramo, ramo_c in ci['rami']:
            Fx_b, Fy_b, M_hinge = sp.Integer(0), sp.Integer(0), sp.Integer(0)
            nodi_corpo = set(ramo_c) | ({nc} if ci['attach'] is not None and ci['attach'] == id_asta_ramo else set())
            for nid_ in nodi_corpo:
                f0, f1 = Forze_S0[nid_], Forze_S1[nid_]
                fx_n = f0['Fx'] + X1_exact * f1['Fx']
                fy_n = f0['Fy'] + X1_exact * f1['Fy']
                m_n = f0['M'] + X1_exact * f1['M']
                Fx_b += fx_n
                Fy_b += fy_n
                M_hinge += m_n + fy_n * (nodi[nid_][0] - xc) - fx_n * (nodi[nid_][1] - yc)
            for el_load in dati_json.get('element_loads', []):
                c_ = _carico(el_load)
                if c_ is None: continue
                a_, b_ = int(c_[0]['n1']), int(c_[0]['n2'])
                if (a_ in ramo_c or a_ == nc) and (b_ in ramo_c or b_ == nc) and not (a_ == nc and b_ == nc):
                    fxq, fyq, xq, yq = _risultante(c_)
                    Fx_b += fxq
                    Fy_b += fyq
                    M_hinge += fyq * (xq - xc) - fxq * (yq - yc)
            Hx, Hy = pulisci_espr(sp.simplify(-Fx_b)), pulisci_espr(sp.simplify(-Fy_b))
            m_res = float(sp.N(M_hinge))
            if abs(m_res) > 1e-6: cerniere_ok = False
            righe_cern.append(
                f"<tr><td>{nc}</td><td>Ramo di Asta {id_asta_ramo} (nodi {', '.join(str(n_) for n_ in sorted(ramo_c))})</td>"
                f"<td><b>{float(Hx):.2f}</b></td><td><b>{float(Hy):.2f}</b></td><td>{m_res:.2e}</td></tr>")

        if righe_cern:
            report.append("<h3>Forze trasmesse dalle cerniere interne</h3>")
            report.append("<p>Forza (H, V) che il resto della struttura esercita, attraverso la cerniera, sul ramo indicato "
                          "(sull'altro lato agisce la forza opposta). L'ultima colonna &egrave; il momento in cerniera "
                          "calcolato dall'equilibrio del ramo: deve risultare nullo.</p>")
            report.append("<table><tr><th>Nodo cerniera</th><th>Ramo</th><th>H [N]</th><th>V [N]</th><th>M in cerniera [N&middot;m]</th></tr>")
            report.extend(righe_cern)
            report.append("</table>")

        if equilibrio_ok and cerniere_ok:
            report.append("<div style='background-color:#d1e7dd; color:#0f5132; padding:10px; border-left:4px solid #198754; margin-bottom:15px;'>"
                          f"&#10004; Verifica di equilibrio globale con le reazioni finali: residui &Sigma;F<sub>x</sub> = {residui[0]:.1e}, "
                          f"&Sigma;F<sub>y</sub> = {residui[1]:.1e}, &Sigma;M = {residui[2]:.1e}"
                          + (" &mdash; momento nullo in tutte le cerniere." if righe_cern else ".") + "</div>")
        else:
            report.append("<div style='background-color:#f8d7da; color:#842029; padding:10px; border-left:4px solid #dc3545; margin-bottom:15px;'>"
                          f"&#10008; Verifica di equilibrio NON soddisfatta: residui &Sigma;F<sub>x</sub> = {residui[0]:.3e}, "
                          f"&Sigma;F<sub>y</sub> = {residui[1]:.3e}, &Sigma;M = {residui[2]:.3e}"
                          + ("; momento non nullo in almeno una cerniera." if not cerniere_ok else ".")
                          + " I risultati vanno considerati inaffidabili.</div>")

        ULTIMI_RISULTATI['reazioni'] = [(rf['tipo'], rf['nodo'], float(rf['val'])) for rf in reaz_fin]
        ULTIMI_RISULTATI['X1'] = float(X1_exact)
        ULTIMI_RISULTATI['equilibrio_globale'] = max(abs(v) for v in residui)

        # ==========================================
        # FASE 8-ter: PENDOLI (BIELLE INCLINATE) - TIRANTE O PUNTONE E DIMENSIONAMENTO DELLA SEZIONE
        # ==========================================
        pendoli = [rf for rf in reaz_fin if rf['tipo'] == 'Fn' and rf.get('pendolo')]
        ULTIMI_RISULTATI['pendoli'] = []
        if pendoli:
            sigma_amm_p = float(dati_json.get('material', {}).get('sigma_amm', 235e6))
            MS_op_p = dati_json.get('material', {}).get('MS_design_operator', dati_json.get('material', {}).get('MS_operator', '>='))
            MS_t_p = float(dati_json.get('material', {}).get('MS_design_target', dati_json.get('material', {}).get('MS_target', 0.0)))
            op_p, _sig_op_p = _simboli_op(MS_op_p)
            ps = dati_json.get('pendolo_section') or {}
            forma = str(ps.get('shape', 'area'))
            modo_p = str(ps.get('mode', 'design'))
            k_rett = float(ps.get('ratio', 1) or 1)
            report.append("<h2>8-ter. Pendoli (bielle): tirante o puntone e sezione</h2>")
            report.append("<p>Il pendolo &egrave; una biella cerniera-cerniera: trasmette solo una forza assiale N lungo il suo asse "
                          "(deformabilit&agrave; assiale trascurata, quindi &egrave; un vincolo rigido). "
                          "N = R<sub>n</sub> &egrave; positiva se la biella <b>tira</b> la struttura verso il vincolo a terra "
                          "(<b>tirante</b>, trazione); N negativa = la biella <b>spinge</b> (<b>puntone</b>, compressione).</p>")
            if sigma_t_ok := (1 + MS_t_p > 1e-9):
                sigma_t_p = sigma_amm_p / (1 + MS_t_p)
            else:
                sigma_t_p = float('nan')
            report.append("<table><tr><th>Pendolo</th><th>Nodo</th><th>Asse (da +x, verso terra)</th><th>N [N]</th><th>Funzionamento</th></tr>")
            righe_sez = []
            for rf in pendoli:
                sup_p = supports_n[rf['sup_idx']] if rf.get('sup_idx') is not None else {}
                ang_ax = float(sup_p.get('axis_angle', float(rf['angle']) + 90))
                Np = float(rf['val'])
                if abs(Np) < 1e-6:
                    funz, col = "scarico (N = 0)", "#6c757d"
                elif Np > 0:
                    funz, col = "TIRANTE (trazione)", "#0d6efd"
                else:
                    funz, col = "PUNTONE (compressione)", "#dc3545"
                report.append(f"<tr><td>{rf['nome']}</td><td>{rf['nodo']}</td><td>{ang_ax:.2f}&deg;</td>"
                              f"<td><b>{Np:.2f}</b></td><td style='color:{col}; font-weight:bold;'>{funz}</td></tr>")
                info_p = {'nodo': rf['nodo'], 'N': Np, 'funzionamento': funz.split()[0].lower(), 'asse_deg': ang_ax}

                # --- sezione del pendolo: sigma = |N|/A, MS = sigma_amm/sigma - 1 ---
                blocco_sez = [f"<h4>{rf['nome']} (nodo {rf['nodo']}): {funz}</h4>"]
                L_pend = None
                if sup_p.get('pend_dx') is not None and sup_p.get('pend_dy') is not None:
                    L_pend = math.hypot(float(sup_p['pend_dx']), float(sup_p['pend_dy']))
                    blocco_sez.append(f"<p>Lunghezza del pendolo: L = {L_pend:.4f} m.</p>")
                if modo_p == 'verifica':
                    if forma == 'rettangolare':
                        A_p = float(ps.get('b', 0)) * float(ps.get('h', 0)); desc = f"b = {float(ps.get('b',0))} m, h = {float(ps.get('h',0))} m"
                    elif forma == 'circolare':
                        A_p = math.pi * float(ps.get('d', 0)) ** 2 / 4; desc = f"d = {float(ps.get('d',0))} m"
                    else:
                        A_p = float(ps.get('A', 0)); desc = ""
                    if A_p <= 0:
                        blocco_sez.append(_box_errore("Sezione del pendolo mancante", "Per la verifica indica le dimensioni (A, oppure b e h, oppure d) della sezione del pendolo."))
                    else:
                        sig_p = abs(Np) / A_p
                        MS_p = (sigma_amm_p / sig_p - 1) if sig_p > 0 else float('inf')
                        ok_p = _verifica_op(MS_p, MS_op_p, MS_t_p)
                        blocco_sez.append(f"<p>Sezione assegnata {desc}: A = {A_p:.4e} m&sup2;. "
                                          f"&sigma; = |N|/A = {abs(Np):.2f}/{A_p:.4e} = <b>{sig_p/1e6:.3f} MPa</b>; "
                                          f"MS = &sigma;<sub>lim</sub>/&sigma; &minus; 1 = <b>{MS_p:.3f}</b> (richiesto {op_p} {MS_t_p:g}) &rarr; "
                                          + ("<b style='color:#198754;'>verifica soddisfatta &#10004;</b>" if ok_p else "<b style='color:#dc3545;'>verifica NON soddisfatta &#10008;</b>") + "</p>")
                        info_p.update({'A': A_p, 'sigma': sig_p, 'MS': MS_p, 'ok': ok_p})
                else:
                    if not sigma_t_ok:
                        blocco_sez.append(_box_errore("Dimensionamento impossibile", "Con MS richiesto &le; &minus;1 la tensione obiettivo non &egrave; definita."))
                    elif abs(Np) < 1e-6:
                        blocco_sez.append("<p>Forza nulla: &sigma; = 0 per qualsiasi sezione, nessun limite sull'area.</p>")
                    else:
                        A_req = abs(Np) / sigma_t_p
                        blocco_sez.append(f"<p>MS {op_p} {MS_t_p:g} &rArr; &sigma; = |N|/A {_sig_op_p} &sigma;<sub>amm</sub>/(1+MS) = {sigma_amm_p/1e6:.2f}/{1 + MS_t_p:g} = {sigma_t_p/1e6:.3f} MPa "
                                          f"&rArr; <b>A {op_p} |N|&middot;(1+MS)/&sigma;<sub>amm</sub> = {A_req:.4e} m&sup2; = {A_req*1e4:.3f} cm&sup2;</b></p>")
                        dim_txt = ""
                        if forma == 'rettangolare':
                            b_p = math.sqrt(A_req / k_rett)
                            dim_txt = f"Sezione rettangolare con h/b = {k_rett:g}: <b>b {op_p} {b_p*100:.3f} cm, h = {k_rett:g}&middot;b {op_p} {k_rett*b_p*100:.3f} cm</b>"
                            info_p.update({'b': b_p, 'h': k_rett * b_p})
                        elif forma == 'circolare':
                            d_p = math.sqrt(4 * A_req / math.pi)
                            dim_txt = f"Sezione circolare piena: <b>d {op_p} {d_p*100:.3f} cm</b>"
                            info_p['d'] = d_p
                        if dim_txt: blocco_sez.append(f"<p>{dim_txt}</p>")
                        info_p['A_req'] = A_req
                # --- stabilita' (carico di punta di Eulero) per il pendolo che lavora da puntone: P_cr = pi^2 E I_min / L^2 (cerniera-cerniera) ---
                if Np < -1e-6:
                    E_p = float(dati_json.get('material', {}).get('E', 2.1e11))
                    I_p = None
                    if modo_p == 'verifica':
                        if forma == 'rettangolare':
                            bb_, hh_ = float(ps.get('b', 0)), float(ps.get('h', 0))
                            I_p = max(bb_, hh_) * min(bb_, hh_) ** 3 / 12.0
                        elif forma == 'circolare':
                            I_p = math.pi * float(ps.get('d', 0)) ** 4 / 64.0
                    else:
                        if info_p.get('b') is not None:
                            I_p = max(info_p['b'], info_p['h']) * min(info_p['b'], info_p['h']) ** 3 / 12.0
                        elif info_p.get('d') is not None:
                            I_p = math.pi * info_p['d'] ** 4 / 64.0
                    if L_pend is None or L_pend <= 0:
                        blocco_sez.append("<p><small>Stabilit&agrave; del puntone: lunghezza del pendolo non nota (mancano pend_dx, pend_dy), carico critico non calcolato.</small></p>")
                    elif I_p is None or I_p <= 0:
                        I_req = abs(Np) * L_pend ** 2 / (math.pi ** 2 * E_p)
                        blocco_sez.append(f"<p>Stabilit&agrave; del puntone (Eulero, cerniera-cerniera): per P<sub>cr</sub> = &pi;&sup2;EI/L&sup2; &ge; |N| serve "
                                          f"<b>I<sub>min</sub> &ge; |N|L&sup2;/(&pi;&sup2;E) = {I_req:.4e} m<sup>4</sup></b> "
                                          "(forma della sezione non assegnata: non si pu&ograve; calcolare P<sub>cr</sub>).</p>")
                        info_p['I_min_richiesto'] = I_req
                    else:
                        P_cr = math.pi ** 2 * E_p * I_p / L_pend ** 2
                        lam = P_cr / abs(Np)
                        ok_e = P_cr >= abs(Np)
                        blocco_sez.append(f"<p>Stabilit&agrave; del puntone (Eulero, cerniera-cerniera, L = {L_pend:.4f} m): I<sub>min</sub> = {I_p:.4e} m<sup>4</sup>, "
                                          f"P<sub>cr</sub> = &pi;&sup2;EI<sub>min</sub>/L&sup2; = <b>{P_cr:.2f} N</b>; P<sub>cr</sub>/|N| = <b>{lam:.3f}</b> &rarr; "
                                          + ("<b style='color:#198754;'>stabile &#10004;</b>" if ok_e else "<b style='color:#dc3545;'>INSTABILE (P<sub>cr</sub> &lt; |N|) &#10008;</b>") + "</p>")
                        info_p.update({'P_cr': P_cr, 'lambda_euler': lam, 'stabile': ok_e, 'I_min': I_p, 'L': L_pend})
                righe_sez.append("".join(blocco_sez))
                ULTIMI_RISULTATI['pendoli'].append(info_p)
            report.append("</table>")
            report.extend(righe_sez)
            if any(i_['funzionamento'] == 'puntone' for i_ in ULTIMI_RISULTATI['pendoli']):
                report.append("<p><small>Per il puntone &egrave; riportata anche la verifica di stabilit&agrave; (carico critico euleriano, cerniera-cerniera): "
                              "&egrave; indicativa, la prova potrebbe richiedere solo la verifica a compressione semplice.</small></p>")

        # ==========================================
        # FASE 8-quater: COERENZA M, T, N FINALI - EQUILIBRIO DI OGNI NODO
        # Le forze e i momenti trasmessi dalle aste ai nodi (calcolati dalle espressioni finali N, T, M) devono equilibrare
        # carichi concentrati e reazioni applicati al nodo. Un errore di segno in N, T o M farebbe fallire questo controllo.
        # ==========================================
        s_nodo = sp.Symbol("s")
        Lmap = {e_['id']: hyp(nodi[int(e_['n2'])][0] - nodi[int(e_['n1'])][0], nodi[int(e_['n2'])][1] - nodi[int(e_['n1'])][1]) for e_ in elementi}
        hinged_per_nodo = {}
        for c_ in cerniere_interne:
            nc_ = int(c_['node_id'])
            inc_ = [e_ for e_ in elementi if nc_ in (int(e_['n1']), int(e_['n2']))]
            hs_ = {int(x_) for x_ in c_['elems']} if c_.get('elems') else ({int(inc_[0]['id'])} if inc_ else set())
            hinged_per_nodo.setdefault(nc_, set()).update(hs_)
        forze_barra = {n_: [0.0, 0.0] for n_ in nodi}      # somma delle forze delle aste sul nodo
        mom_rigidi = {n_: 0.0 for n_ in nodi}              # momenti delle aste collegate rigidamente al nodo
        mom_cerniera = []                                  # (nodo, asta, momento) per le aste incernierate
        for el in elementi:
            eid_ = el['id']
            n1_, n2_ = int(el['n1']), int(el['n2'])
            Lf_ = float(Lmap[eid_])
            if Lf_ == 0: continue
            if grado_iperstaticita > 0:
                Mf_ = M0_dict[eid_] + X1_exact * M1_dict[eid_]; Tf_ = T0_dict[eid_] + X1_exact * T1_dict[eid_]; Nf_ = N0_dict[eid_] + X1_exact * N1_dict[eid_]
            else:
                Mf_, Tf_, Nf_ = M0_dict[eid_], T0_dict[eid_], N0_dict[eid_]
            sv_ = np.array([0.0, Lf_])
            Mv_, Tv_, Nv_ = _valuta_poly(Mf_, sv_), _valuta_poly(Tf_, sv_), _valuta_poly(Nf_, sv_)
            ca_ = float((nodi[n2_][0] - nodi[n1_][0]) / Lmap[eid_]); sa_ = float((nodi[n2_][1] - nodi[n1_][1]) / Lmap[eid_])
            ux_, uy_, nx_, ny_ = ca_, sa_, -sa_, ca_
            # estremo iniziale (nodo n1): forza dell'asta sul nodo = N u - T n ; momento = -M
            forze_barra[n1_][0] += Nv_[0] * ux_ - Tv_[0] * nx_
            forze_barra[n1_][1] += Nv_[0] * uy_ - Tv_[0] * ny_
            # estremo finale (nodo n2): forza = -N u + T n ; momento = +M
            forze_barra[n2_][0] += -Nv_[1] * ux_ + Tv_[1] * nx_
            forze_barra[n2_][1] += -Nv_[1] * uy_ + Tv_[1] * ny_
            for nodo_e, m_e in ((n1_, -Mv_[0]), (n2_, Mv_[1])):
                if int(eid_) in hinged_per_nodo.get(nodo_e, set()):
                    mom_cerniera.append((nodo_e, eid_, float(m_e)))
                else:
                    mom_rigidi[nodo_e] += float(m_e)
        scala_f = max([1.0] + [abs(float(v_)) for n_ in nodi for v_ in (Forze_S0[n_]['Fx'], Forze_S0[n_]['Fy'])])
        scala_m = max([1.0] + [abs(float(Forze_S0[n_]['M'])) for n_ in nodi] + [scala_f * float(max(Lmap.values()))])
        res_nodi = []
        mt_nodo = {}
        carrier_nodo = {ci_['nodo']: ci_['carrier'] for ci_ in cerniere_info}
        for n_ in nodi:
            ft_x = float(sp.N(Forze_S0[n_]['Fx'] + X1_exact * Forze_S1[n_]['Fx']))
            ft_y = float(sp.N(Forze_S0[n_]['Fy'] + X1_exact * Forze_S1[n_]['Fy']))
            mt_ = float(sp.N(Forze_S0[n_]['M'] + X1_exact * Forze_S1[n_]['M']))
            mt_nodo[n_] = mt_
            rx_, ry_ = forze_barra[n_][0] + ft_x, forze_barra[n_][1] + ft_y
            car_ = carrier_nodo.get(n_)
            if car_ is not None:
                # coppie e momenti reattivi del nodo appartengono a un'asta incernierata: il suo momento d'estremo deve equilibrarli,
                # mentre le aste collegate rigidamente al nodo (se esistono) si equilibrano da sole
                rm_ = mom_rigidi[n_]
            else:
                rm_ = mom_rigidi[n_] + mt_
            res_nodi.append((n_, rx_, ry_, rm_))
        res_cern = [(n_, e_, m_ + (mt_nodo[n_] if carrier_nodo.get(n_) is not None and int(e_) == int(carrier_nodo[n_]) else 0.0))
                    for n_, e_, m_ in mom_cerniera]
        max_res_f = max([abs(r_[1]) for r_ in res_nodi] + [abs(r_[2]) for r_ in res_nodi] + [0.0])
        max_res_m = max([abs(r_[3]) for r_ in res_nodi] + [abs(c_[2]) for c_ in res_cern] + [0.0])
        nodi_ko = [r_[0] for r_ in res_nodi if max(abs(r_[1]), abs(r_[2])) > 1e-6 * scala_f or abs(r_[3]) > 1e-6 * scala_m]
        cern_ko = [(c_[0], c_[1]) for c_ in res_cern if abs(c_[2]) > 1e-6 * scala_m]
        ULTIMI_RISULTATI['equilibrio_nodi'] = (max_res_f, max_res_m)
        ULTIMI_RISULTATI['nodi_non_in_equilibrio'] = nodi_ko + [c_[0] for c_ in cern_ko]
        report.append("<h3>Controllo di coerenza: equilibrio di ogni nodo (N, T, M finali)</h3>")
        if not nodi_ko and not cern_ko:
            report.append("<div style='background-color:#d1e7dd; color:#0f5132; padding:10px; border-left:4px solid #198754; margin-bottom:15px;'>"
                          f"&#10004; Le forze (N, T) e i momenti M che le aste trasmettono ai nodi equilibrano carichi e reazioni in tutti i nodi: "
                          f"residuo massimo forze = {max_res_f:.1e} N, momenti = {max_res_m:.1e} N&middot;m; momento nullo agli estremi incernierati.</div>")
        else:
            elenco_ko = ", ".join(f"N{n_}" for n_ in nodi_ko) + (", cerniere: " + ", ".join(f"asta {e_} al nodo N{n_}" for n_, e_ in cern_ko) if cern_ko else "")
            report.append("<div style='background-color:#f8d7da; color:#842029; padding:10px; border-left:4px solid #dc3545; margin-bottom:15px;'>"
                          f"&#10008; Equilibrio ai nodi NON soddisfatto ({elenco_ko}): residuo massimo forze = {max_res_f:.3e} N, momenti = {max_res_m:.3e} N&middot;m. "
                          "Le sollecitazioni calcolate sono incoerenti con carichi e reazioni: non fidarti dei risultati.</div>")

        # ==========================================
        # FASE 9: SEZIONE CRITICA E VERIFICA DI RESISTENZA (Navier)
        # ==========================================
        report.append("<h2 style='margin-top: 30px;'>9. Sezione Critica e Verifica di Resistenza (Navier)</h2>")
        
        # Acquisizione tipo verifica dal JSON
        ver_type = dati_json.get('material', {}).get('verification_type', 'combined')
        
        if ver_type == 'bending_only':
            report.append("<p>Si individua la sezione critica cercando il punto in cui il <b>Momento Flettente |M(s)|</b> raggiunge il suo massimo valore assoluto (Criterio: Flessione Semplice).</p>")
        else:
            report.append("<p>Si individua la sezione critica cercando il punto in cui la <b>tensione normale combinata</b> &sigma;(s) = |N(s)|/A + |M(s)|/W raggiunge il suo massimo valore assoluto (Criterio: Presso/Tensoflessione).</p>")
        
        # 2. Geometria delle Masse
        b_sec = float(dati_json.get('section', {}).get('b', 0.2))
        h_sec = float(dati_json.get('section', {}).get('h', 0.4))
        Area = b_sec * h_sec
        W_f = (b_sec * (h_sec**2)) / 6
        
        max_sigma_calcolata = 0.0
        asta_critica = None
        s_critico = 0.0
        M_in_s_critico = 0.0
        N_in_s_critico = 0.0
        sezioni_equivalenti = []

        s_sym = sp.Symbol("s")

        # 1. Espressioni finali M(s), N(s) di ogni asta (M0 + X1*M1, N0 + X1*N1)
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
            """Punti (id asta, s, N, M) in cui puo' cadere il massimo di sigma: estremi, radici e punti stazionari ESATTI
            di M(s), N(s) (sono polinomi), piu' una griglia di sicurezza. Non dipende dal passo di campionamento."""
            punti = []
            for el_ in elementi:
                eid_ = el_['id']
                sv = candidati_s(M_fin_dict[eid_], N_fin_dict[eid_], L_dict[eid_], pesi=pesi_, n_base=101)
                punti.append((eid_, sv, _valuta_poly(N_fin_dict[eid_], sv), _valuta_poly(M_fin_dict[eid_], sv)))
            return punti

        # nella verifica la sezione e' nota: wN/wM = (1/A)/(1/W) = W/A = h/6. Nel dimensionamento si coprono piu' rapporti.
        pesi_ver = [(0.0, 1.0)] if ver_type == 'bending_only' else [(h_sec / 6.0, 1.0)]
        punti_NM = _scansiona(PESI_DECADI + pesi_ver)   # servono anche al dimensionamento

        def _sigma_punti(N_, M_):
            return np.abs(M_) / W_f if ver_type == 'bending_only' else np.abs(N_) / Area + np.abs(M_) / W_f

        for el_id, s_vals, N_vals, M_vals in punti_NM:
            Sigma_vals = _sigma_punti(N_vals, M_vals)
            idx_max = int(np.argmax(Sigma_vals))
            local_max_sigma = float(Sigma_vals[idx_max])
            if local_max_sigma > max_sigma_calcolata * (1 + 1e-9):   # parita': resta la prima asta (id piu' basso)
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

        # 3. Calcolo Margine di Sicurezza (MS) effettivo
        sigma_amm = float(dati_json.get('material', {}).get('sigma_amm', 235e6))
        
        if max_sigma_calcolata > 0:
            MS = (sigma_amm / max_sigma_calcolata) - 1
        else:
            MS = float('inf')
            
        MS_target = float(dati_json.get('material', {}).get('MS_target', 0.0))
        MS_op = dati_json.get('material', {}).get('MS_operator', '>=')
        
        ULTIMI_RISULTATI['MS'] = MS
        ULTIMI_RISULTATI['sigma_max'] = max_sigma_calcolata
        ULTIMI_RISULTATI['sezione_critica'] = (asta_critica, s_critico, N_in_s_critico, M_in_s_critico)

        # Gestione segni
        TOL_MS = 1e-9   # evita che l'arrotondamento faccia fallire un MS pari al valore richiesto
        if MS_op == ">=":
            op_html = "&ge;"
            sigma_op_html = "&le;"
            is_verified = (MS >= MS_target - TOL_MS)
        elif MS_op == ">":
            op_html = "&gt;"
            sigma_op_html = "&lt;"
            is_verified = (MS > MS_target + TOL_MS)
        elif MS_op == "<=":
            op_html = "&le;"
            sigma_op_html = "&ge;"
            is_verified = (MS <= MS_target + TOL_MS)
        elif MS_op == "<":
            op_html = "&lt;"
            sigma_op_html = "&gt;"
            is_verified = (MS < MS_target - TOL_MS)
        else:
            op_html = "="
            sigma_op_html = "="
            is_verified = math.isclose(MS, MS_target, abs_tol=1e-2) if MS != float('inf') else False

        # MS della RICHIESTA DI PROGETTO (punto 4 dell'esame: relazione fra i parametri di sezione). E' indipendente dall'MS di
        # verifica del punto 3; se i campi dedicati mancano (vecchi JSON) si ricade sui valori di verifica.
        _mat_ms = dati_json.get('material', {})
        MS_d_target = float(_mat_ms.get('MS_design_target', MS_target))
        MS_d_op = _mat_ms.get('MS_design_operator', MS_op)
        op_html_d, sigma_op_html_d = _simboli_op(MS_d_op)

        # Output HTML
        report.append("<div style='background-color: #f8f9fa; padding: 20px; border-left: 5px solid #0a58ca; border-radius: 4px; margin-bottom: 20px;'>")
        report.append("<h4 style='color: #0a58ca; margin-top:0;'>1. Individuazione della Sezione Critica</h4>")
        if asta_critica is not None:
            report.append(f"<ul style='margin-bottom: 15px;'><li>La tensione normale massima (secondo il criterio scelto) si trova nell'<b>Asta {asta_critica}</b> all'ascissa locale <b>s = {s_critico:.3f} m</b>.</li>")
            report.append(f"<li>Sollecitazioni associate: |N| = {abs(N_in_s_critico):.2f} N, |M| = {abs(M_in_s_critico):.2f} N&middot;m</li>")
            if len(sezioni_equivalenti) > 1:
                elenco_eq = "; ".join(f"asta {a_} (s = {s_:.3f} m)" for a_, s_ in sezioni_equivalenti[:6])
                report.append(f"<li>Sezioni con la stessa tensione massima (equivalenti): {elenco_eq}.</li>")
            report.append("</ul>")
        else:
            report.append(f"<ul style='margin-bottom: 15px;'><li>Sollecitazioni nulle su tutta la struttura.</li></ul>")
        
        report.append("<h4 style='color: #0a58ca;'>2. Caratteristiche della Sezione e del Materiale</h4>")
        report.append(f"<ul style='margin-bottom: 15px;'><li>Base (b) = {b_sec} m, Altezza (h) = {h_sec} m</li>")
        report.append(f"<li>Area della sezione A = b&middot;h = <b>{Area:.4e} m&sup2;</b></li>")
        report.append(f"<li>Modulo di resistenza flessionale W<sub>f</sub> = (b&middot;h&sup2;) / 6 = <b>{W_f:.4e} m&sup3;</b></li>")
        report.append(f"<li>Tensione limite/ammissibile del materiale: &sigma;<sub>lim</sub> = <b>{(sigma_amm/1e6):.2f} MPa</b></li></ul>")
        
        design_mode = str(dati_json.get('section', {}).get('design', 'none'))
        in_progetto = design_mode in ('h', 'b', 'ratio')
        if in_progetto:
            report.append("<h4 style='color: #0a58ca;'>3. Margine di Sicurezza con la sezione di prova (b, h assegnati)</h4>")
            report.append("<p><i>La sezione e' da dimensionare: qui b e h sono solo valori di prova per individuare la sezione critica. "
                          "La verifica che conta &egrave; il dimensionamento al punto 10.</i></p>")
        else:
            report.append("<h4 style='color: #0a58ca;'>3. Verifica di Resistenza e Margine di Sicurezza (MS)</h4>")
        report.append("<p style='margin-bottom: 10px;'>Applicando la formula di Navier per la sollecitazione:</p>")
        report.append(f"<div style='margin-bottom: 15px;'>")
        
        if ver_type == 'bending_only':
            report.append(f"<span class='formula' style='color: #000000;'>&sigma;<sub>max</sub> = |M| / W<sub>f</sub> = {abs(M_in_s_critico):.2f} / {W_f:.4e} = {(max_sigma_calcolata/1e6):.2f} MPa</span><br><br>")
        else:
            report.append(f"<span class='formula' style='color: #000000;'>&sigma;<sub>max</sub> = |N| / A + |M| / W<sub>f</sub> = {abs(N_in_s_critico):.2f} / {Area:.4e} + {abs(M_in_s_critico):.2f} / {W_f:.4e} = {(max_sigma_calcolata/1e6):.2f} MPa</span><br><br>")
        
        if MS != float('inf'):
            report.append(f"<span class='formula' style='color: #000000;'>MS = (&sigma;<sub>lim</sub> / &sigma;<sub>max</sub>) - 1 = {(sigma_amm/1e6):.2f} / {(max_sigma_calcolata/1e6):.2f} - 1 = <b>{MS:.2f}</b></span>")
        else:
            report.append(f"<span class='formula' style='color: #000000;'>MS = &infin;</span>")
        report.append(f"</div>")
        
        if in_progetto:
            report.append("<div style='text-align:center; padding:12px; background-color:#e9ecef; color:#343a40; border: 1px solid #adb5bd; font-size: 15px; margin-top: 10px; border-radius: 6px;'>")
            report.append(f"MS (sezione di prova) = {MS:.2f}" if MS != float('inf') else "MS = &infin;")
            report.append("<br><small>Nessun esito di verifica: vedi dimensionamento (punto 10).</small></div>")
        elif is_verified:
            report.append("<div style='text-align:center; padding:15px; background-color:#d1e7dd; color:#0f5132; border: 2px solid #198754; font-size: 16px; margin-top: 10px; border-radius: 6px;'>")
            report.append(f"<b>MS = {MS:.2f} {op_html} {MS_target:.2f}</b><br>")
            report.append("LA VERIFICA È SODDISFATTA &#10004;")
            report.append("</div>")
        else:
            report.append("<div style='text-align:center; padding:15px; background-color:#f8d7da; color:#842029; border: 2px solid #dc3545; font-size: 16px; margin-top: 10px; border-radius: 6px;'>")
            report.append(f"<b>MS = {MS:.2f} (Richiesto {op_html} {MS_target:.2f})</b><br>")
            report.append("LA VERIFICA NON È SODDISFATTA &#10008;")
            report.append("</div>")
        report.append("</div>")

        # ==========================================
        # FASE 10: DIMENSIONAMENTO DELLA SEZIONE (Dinamico)
        # ==========================================
        def frac(num, den):
            return f"<table cellpadding='0' cellspacing='0' style='margin: 0; display: inline-table; vertical-align: middle;'><tr><td style='border-bottom: 1px solid #000000; text-align: center; padding: 0 4px;'>{num}</td></tr><tr><td style='text-align: center; padding: 0 4px;'>{den}</td></tr></table>"
            
        def eq_table(*elements):
            html = "<table cellpadding='0' cellspacing='0' style='margin-bottom: 10px; font-family: \"Courier New\", Courier, monospace; font-size: 15px; color: #000000;'><tr>"
            for el in elements:
                html += f"<td valign='middle' style='padding: 0 4px;'>{el}</td>"
            html += "</tr></table>"
            return html

        report.append("<h2>10. Dimensionamento della sezione</h2>")

        sigma_amm_val = float(dati_json.get('material', {}).get('sigma_amm', 235e6))
        sigma_target_val = sigma_amm_val / (1 + MS_d_target) if 1 + MS_d_target > 1e-9 else float('nan')

        # ---- calcolo numerico del dimensionamento (b, h): fatto PRIMA di scrivere la disequazione letterale, cosi' la
        # disequazione usa le sollecitazioni della sezione che governa davvero il progetto (non quelle della sezione di prova)
        design = str(dati_json.get('section', {}).get('design', 'none'))
        ratio_k = float(dati_json.get('section', {}).get('ratio', 0) or 0)
        combinato = (ver_type != 'bending_only')
        EI_dir = 'EI' in dati_json.get('material', {})
        sez_dip = (dip_ced and not EI_dir) or (dip_term and (not EI_dir or design in ('h', 'ratio')))
        esito_design = None            # None = non eseguito; altrimenti dict con i risultati
        problema_design = None         # (tipo, titolo, messaggio) se il dimensionamento non si puo' fare
        if design in ('h', 'b', 'ratio'):
            if 1 + MS_d_target <= 1e-9:
                problema_design = ('errore', "Dimensionamento impossibile", "Con MS richiesto &le; &minus;1 la tensione obiettivo non &egrave; definita.")
            elif design == 'ratio' and ratio_k <= 0:
                problema_design = ('errore', "Rapporto h/b mancante", "Per il dimensionamento a rapporto fissato indica h/b &gt; 0.")
            elif sez_dip:
                # Cedimenti / distorsioni termiche con sezione incognita: X1 dipende da EI (e da h per il gradiente termico), quindi
                # le sollecitazioni cambiano con la dimensione d. Si risolve sigma_max(d) = sigma_t numericamente:
                #   X1(d) = [EI(d)*(eta1 - KT/h(d) - int N1 eps_T + sum R1*delta) - d10] / d11
                # sigma_max(d) puo' NON essere monotona (le distorsioni crescono con la rigidezza), quindi si cercano tutte le radici.
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
                    if a_ is None: return f"{nome_d} &le; {b_:.5f} m ({b_*100:.2f} cm)"
                    if b_ is None: return f"{nome_d} &ge; {a_:.5f} m ({a_*100:.2f} cm)"
                    return f"{a_:.5f} m &le; {nome_d} &le; {b_:.5f} m ({a_*100:.2f} cm &le; {nome_d} &le; {b_*100:.2f} cm)"

                if MS_d_op == '=':
                    cand_d = sorted(radici_d)
                    testo_d = (" &nbsp;oppure&nbsp; ".join(f"{nome_d} = {r_:.5f} m ({r_*100:.2f} cm)" for r_ in cand_d)) if cand_d else ""
                else:
                    cand_d = [x_ for iv_ in interv_d for x_ in iv_ if x_ is not None][:1]
                    testo_d = " &nbsp;oppure&nbsp; ".join(_txt_interv(a_, b_) for a_, b_ in interv_d)
                if not cand_d:
                    if interv_d:
                        problema_design = ('avviso', "Dimensionamento: nessun limite",
                                           f"Con le distorsioni assegnate la richiesta &egrave; soddisfatta da {testo_d} (nessun limite nell'intervallo esplorato).")
                    else:
                        problema_design = ('errore', "Dimensionamento senza soluzione",
                                           f"Nessuna dimensione tra {D_LO*100:g} cm e {D_HI:g} m soddisfa la richiesta sul margine di sicurezza "
                                           "(con cedimenti/distorsioni termiche &sigma;<sub>max</sub> non decresce sempre aumentando la sezione).")
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
                    if design == 'b':       # h fissata -> b
                        d_all = (N_all * h_sec + 6 * M_all) / (h_sec ** 2 * sig_t)
                    elif design == 'h':     # b fissata -> h
                        d_all = (N_all + np.sqrt(N_all ** 2 + 24 * sig_t * b_sec * M_all)) / (2 * sig_t * b_sec)
                    else:                   # h = ratio_k * b -> b : sig_t b^3 - (N/k) b - 6M/k^2 = 0
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
                # raffinamento: con la sezione trovata il rapporto W/A e' noto, quindi si ricalcolano i punti critici esatti
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
        else:
            N_E_abs, M_E_abs = abs(N_in_s_critico), abs(M_in_s_critico)
        sei_Me = 6 * M_E_abs

        if sigma_target_val == sigma_target_val:     # non NaN
            mantissa, exp = f"{sigma_target_val:.3e}".split('e')
            sigma_target_str = f"{float(mantissa):.3f}&middot;10^{int(exp)}"
        else:
            sigma_target_str = "n.d."

        bl_style = "border-left: 3px solid #0d6efd; background-color: #f8f9fa; padding-left: 15px; padding-top: 10px; padding-bottom: 10px; margin-bottom: 25px; margin-top: 10px;"

        report.append(f"<p>Si mantiene la sezione in forma <b>letterale</b> (base b e altezza h come incognite) e si impone la richiesta di progetto sul margine di sicurezza (MS {op_html_d} {MS_d_target}):</p>")
        report.append(f"<div style=\"{bl_style}\">")
        report.append(eq_table(f"MS {op_html_d} {MS_d_target} &rarr;", frac("&sigma;<sub>amm</sub>", "&sigma;<sub>eq,max</sub>"), f"- 1 {op_html_d} {MS_d_target} &rarr; &sigma;<sub>eq,max</sub> {sigma_op_html_d}", frac("&sigma;<sub>amm</sub>", f"1 + {MS_d_target}")))
        report.append("</div>")

        report.append("<p>Esprimendo la tensione massima con A = b&middot;h e W = b&middot;h&sup2;/6 (b e h letterali):</p>")

        report.append(f"<div style=\"{bl_style}\">")
        if ver_type == 'bending_only':
            report.append(eq_table("&sigma;<sub>eq,max</sub> =", frac("|M<sub>E</sub>|", "W"), "=", frac("6&middot;|M<sub>E</sub>|", "b&middot;h&sup2;")))
        else:
            report.append(eq_table("&sigma;<sub>eq,max</sub> =", frac("|N<sub>E</sub>|", "A"), "+", frac("|M<sub>E</sub>|", "W"), "=", frac("|N<sub>E</sub>|", "b&middot;h"), "+", frac("6&middot;|M<sub>E</sub>|", "b&middot;h&sup2;"), "=", frac("h&middot;|N<sub>E</sub>| + 6&middot;|M<sub>E</sub>|", "b&middot;h&sup2;")))
        report.append("</div>")

        if esito_design is not None:
            report.append(f"<p>Sezione che governa il progetto: <b>asta {esito_design['asta']}, s = {esito_design['s']:.3f} m</b> "
                          f"(|N<sub>E</sub>| = {N_E_abs:.2f} N, |M<sub>E</sub>| = {M_E_abs:.2f} N&middot;m).</p>")
        report.append(f"<p>Sostituendo i valori numerici (in modulo) e riducendo il termine destro <span class='formula'>&sigma;<sub>amm</sub> / (1 + {MS_d_target})</span> :</p>")

        report.append(f"<div style=\"{bl_style}\">")
        if ver_type == 'bending_only':
            report.append(eq_table(frac(f"{sei_Me:.2f}", "b&middot;h&sup2;"), sigma_op_html_d, f"{sigma_target_str} N/m&sup2;"))
        else:
            report.append(eq_table(frac(f"h&middot;{N_E_abs:.2f} + {sei_Me:.2f}", "b&middot;h&sup2;"), sigma_op_html_d, f"{sigma_target_str} N/m&sup2;"))
        report.append("</div>")

        # ---- 10-ter: relazione LETTERALE in b e h quando le distorsioni rendono X1 dipendente dalla sezione ----
        # X1(b,h) = [EI*(eta1 - KT/h - int N1 epsT + sum R1 delta) - d10]/d11 con EI = E*b*h^3/12  =>  N(b,h) e M(b,h) della sezione
        # critica non sono numeri: la disequazione sigma <= sigma_t va scritta come polinomio in b e h.
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
                if esito_design is not None and esito_design.get('distorsioni'):
                    b_r, h_r = esito_design['b'], esito_design['h']
                else:
                    b_r, h_r = b_sec, h_sec
                if EI_dir:
                    EI_l = float(EI_val)
                else:
                    EI_l = None
                def _X1_num(b_, h_):
                    EI_ = EI_l if EI_l is not None else E_l * b_ * h_ ** 3 / 12.0
                    return (EI_ * (eta_l - KT_l / h_ - TN_l + SR_l) - d10_l) / d11_l
                sig_r, asta_r, s_r, _Nr, _Mr = _sigma_max_numerico(barre_l, _X1_num(b_r, h_r), b_r * h_r, b_r * h_r ** 2 / 6.0,
                                                                 ver_type == 'bending_only')
                bar_r = next(x_ for x_ in barre_l if x_[0] == asta_r)
                n0_, n1_ = float(np.polyval(bar_r[4], s_r)), float(np.polyval(bar_r[5], s_r))
                m0_, m1_ = float(np.polyval(bar_r[2], s_r)), float(np.polyval(bar_r[3], s_r))
                b_s, h_s = sp.symbols('b h', positive=True)
                EI_s = sp.Float(EI_l) if EI_l is not None else sp.Float(E_l) * b_s * h_s ** 3 / 12
                X1_s = (EI_s * (eta_l - KT_l / h_s - TN_l + SR_l) - d10_l) / d11_l
                N_s, M_s = n0_ + X1_s * n1_, m0_ + X1_s * m1_
                sN_, sM_ = (1.0 if float(_Nr) >= 0 else -1.0), (1.0 if float(_Mr) >= 0 else -1.0)
                if ver_type == 'bending_only':
                    lhs_s = 6 * sM_ * M_s
                else:
                    lhs_s = sN_ * N_s * h_s + 6 * sM_ * M_s
                lhs_s = sp.expand(lhs_s)
                def _fmt_s(e_):
                    e_ = sp.expand(e_)
                    termini_ = []
                    for t_ in sp.Add.make_args(e_):
                        c_, m_ = t_.as_coeff_Mul()
                        termini_.append((float(c_), m_))
                    out_ = ""
                    for c_, m_ in sorted(termini_, key=lambda x_: -sp.degree(x_[1], h_s) - sp.degree(x_[1], b_s)):
                        if abs(c_) < 1e-9: continue
                        mon_ = str(m_).replace('**', '^').replace('*', '&middot;')
                        mon_ = "" if mon_ == "1" else "&middot;" + mon_
                        out_ += f" {'+' if c_ >= 0 else '&minus;'} {abs(c_):.4g}{mon_}"
                    return out_.strip().lstrip('+ ').strip() if out_.strip().startswith('+') else out_.strip()
                report.append("<h3>10-ter. Relazione letterale fra b e h (con cedimenti / distorsioni termiche)</h3>")
                report.append("<p>Con cedimenti o &Delta;T l'incognita iperstatica dipende dalla sezione: "
                              "X<sub>1</sub>(b,h) = [EI&middot;(&eta;<sub>1</sub> &minus; K<sub>T</sub>/h &minus; &int;N<sub>1</sub>&epsilon;<sub>T</sub>ds + &Sigma;R<sub>1</sub>&delta;) &minus; &delta;<sub>10</sub>] / &delta;<sub>11</sub>"
                              + ("" if EI_l is not None else " con EI = E&middot;b&middot;h&sup3;/12") + ". "
                              f"Nella sezione critica (asta {asta_r}, s = {s_r:.3f} m) N = N<sub>0</sub> + X<sub>1</sub>N<sub>1</sub> e M = M<sub>0</sub> + X<sub>1</sub>M<sub>1</sub> "
                              "diventano quindi funzioni di b e h.</p>")
                report.append(f"<div style=\"{bl_style}\">")
                report.append(eq_table("X<sub>1</sub>(b,h) =", _fmt_s(X1_s)))
                report.append("</div>")
                if ver_type == 'bending_only':
                    rel_txt = f"{_fmt_s(lhs_s)} {sigma_op_html_d} {sigma_target_str}&middot;b&middot;h&sup2;"
                else:
                    rel_txt = f"{_fmt_s(lhs_s)} {sigma_op_html_d} {sigma_target_str}&middot;b&middot;h&sup2;"
                report.append("<p>Relazione fra i parametri geometrici (|N| e |M| esplicitati con i segni che hanno nella sezione critica di riferimento: "
                              f"N {'&ge;' if sN_ > 0 else '&le;'} 0, M {'&ge;' if sM_ > 0 else '&le;'} 0; "
                              "&sigma;<sub>amm</sub>/(1+MS) = " + sigma_target_str + " Pa):</p>")
                report.append(f"<div style=\"{bl_style}\">")
                report.append(eq_table(rel_txt))
                report.append("</div>")
                report.append("<p><small>La relazione vale finch&eacute; la sezione critica resta quella indicata e i segni di N ed M non cambiano al variare di b e h; "
                              "il dimensionamento numerico della sezione 10-bis controlla tutte le sezioni.</small></p>")
                ULTIMI_RISULTATI['relazione_letterale'] = {'lhs': str(lhs_s), 'sigma_t': sigma_target_val, 'asta': int(asta_r), 's': float(s_r),
                                                           'X1': str(sp.expand(X1_s)), 'op': MS_d_op}
            except Exception as _e_l:
                report.append(_box_avviso("Relazione letterale non generata", f"Non e' stato possibile scrivere la relazione in b e h: {_e_l}"))

        # ---- 10-bis: esito del dimensionamento (b, h) ----
        if problema_design is not None:
            report.append("<h3>10-bis. Risoluzione del dimensionamento</h3>")
            report.append(_box_avviso(problema_design[1], problema_design[2]) if problema_design[0] == 'avviso'
                          else _box_errore(problema_design[1], problema_design[2]))
        elif esito_design is not None:
            report.append("<h3>10-bis. Risoluzione del dimensionamento</h3>")
            b_ref, h_ref, MS_c = esito_design['b'], esito_design['h'], esito_design['MS_c']
            if esito_design.get('distorsioni'):
                report.append(f"<p>Tensione obiettivo: &sigma;<sub>t</sub> = &sigma;<sub>amm</sub>/(1+MS) = {sigma_amm/1e6:.2f}/{1 + MS_d_target:g} = <b>{esito_design['sig_t']/1e6:.3f} MPa</b>. "
                              "Con cedimenti / distorsioni termiche l'incognita iperstatica dipende dalla sezione: "
                              "X<sub>1</sub>(d) = [EI(d)&middot;(&eta;<sub>1</sub> &minus; K<sub>T</sub>/h &minus; &int;N<sub>1</sub>&epsilon;<sub>T</sub>ds + &Sigma;R<sub>1</sub>&delta;) &minus; &delta;<sub>10</sub>] / &delta;<sub>11</sub>. "
                              "Si risolve numericamente &sigma;<sub>max</sub>(d) = &sigma;<sub>t</sub> (con punti critici esatti per ogni d); "
                              "&sigma;<sub>max</sub>(d) pu&ograve; non essere monotona, quindi si cercano tutte le soluzioni in 1 cm &divide; 5 m.</p>")
                report.append(f"<p>Le sollecitazioni e X<sub>1</sub> mostrati nelle sezioni 5-9 sono relativi alla sezione di prova; per la sezione progettata "
                              f"X<sub>1</sub> = <b>{esito_design['X1']:.4f}</b> (sezione critica: asta {esito_design['asta']}, s = {esito_design['s']:.3f} m, "
                              f"|N| = {N_E_abs:.2f} N, |M| = {M_E_abs:.2f} N&middot;m).</p>")
            else:
                report.append(f"<p>Tensione obiettivo: &sigma;<sub>t</sub> = &sigma;<sub>amm</sub>/(1+MS) = {sigma_amm/1e6:.2f}/{1 + MS_d_target:g} = <b>{esito_design['sig_t']/1e6:.3f} MPa</b>. "
                              f"Si cerca, su tutte le sezioni candidate, la dimensione che porta <i>quella</i> sezione a &sigma;<sub>t</sub>; "
                              f"vale quella pi&ugrave; gravosa (asta {esito_design['asta']}, s = {esito_design['s']:.3f} m, |N| = {N_E_abs:.2f} N, |M| = {M_E_abs:.2f} N&middot;m).</p>")
            if esito_design.get('distorsioni'):
                esito = esito_design['testo']
                if design == 'ratio':
                    esito += f"<br><small>(con h = {ratio_k:g}&middot;b; per b = {b_ref:.5f} m: h = {h_ref:.5f} m)</small>"
                elif design == 'b':
                    esito += f"<br><small>(con h = {h_sec} m assegnata)</small>"
                else:
                    esito += f"<br><small>(con b = {b_sec} m assegnata)</small>"
            elif design == 'b':
                esito = f"b {op_html_d} {b_ref:.5f} m ({b_ref*100:.2f} cm) &nbsp;con h = {h_sec} m assegnata"
            elif design == 'h':
                esito = f"h {op_html_d} {h_ref:.5f} m ({h_ref*100:.2f} cm) &nbsp;con b = {b_sec} m assegnata"
            else:
                esito = (f"b {op_html_d} {b_ref:.5f} m ({b_ref*100:.2f} cm), &nbsp;h = {ratio_k:g}&middot;b {op_html_d} {h_ref:.5f} m ({h_ref*100:.2f} cm)<br>"
                         f"<small>(da &sigma;<sub>t</sub>&middot;b&sup3; &minus; (N/k)&middot;b &minus; 6M/k&sup2; = 0 con k = h/b = {ratio_k:g})</small>")
            report.append("<div style='text-align:center; padding:15px; background-color:#d1e7dd; color:#0f5132; border: 2px solid #198754; font-size: 16px; margin-top: 10px; border-radius: 6px;'>")
            report.append(f"<b>{esito}</b><br><small>Controllo con la sezione trovata: MS = {MS_c:.4f}</small></div>")
            report.append("<p><small>I punti critici sono calcolati in forma esatta (estremi, radici e punti stazionari di M(s) e N(s)), non per campionamento.</small></p>")
            ULTIMI_RISULTATI['dimensionamento'] = dict(esito_design)

    # ==========================================
    # FASE 11: GRAFICI MATPLOTLIB
    # ==========================================
    percorso_S1 = "plv_S0_S1.png"
    percorso_diag = "plv_diagrammi.png"
    nodi = {k: (float(v[0]), float(v[1])) for k, v in nodi.items()}   # da qui solo disegno
    
    # Plot S1 (Mantenuto identico)
    fig_s1 = plt.figure(figsize=(8, 4))
    try:
        for el in elementi: plt.plot([nodi[int(el['n1'])][0], nodi[int(el['n2'])][0]], [nodi[int(el['n1'])][1], nodi[int(el['n2'])][1]], 'k-', linewidth=3, zorder=1)
        for nid, coords in nodi.items():
            plt.scatter(coords[0], coords[1], c='white', s=100, zorder=2, edgecolors='black', linewidth=2)
            plt.text(coords[0], coords[1], f"{nid}", fontsize=10, weight='bold', horizontalalignment='center', verticalalignment='center', color='black')
        max_x = max(coords[0] for coords in nodi.values())
        min_x = min(coords[0] for coords in nodi.values())
        L_scale = max((max_x - min_x) * 0.15, 0.5) 
        if grado_iperstaticita > 0:   # xs, ys esistono solo se c'e' un'incognita iperstatica
            if tipo_svinc == 'Fx': plt.annotate('X1', xy=(xs, ys), xytext=(xs + L_scale, ys), arrowprops=dict(facecolor='red', width=3, headwidth=10), color='red', weight='bold')
            elif tipo_svinc == 'Fy': plt.annotate('X1', xy=(xs, ys), xytext=(xs, ys + L_scale), arrowprops=dict(facecolor='red', width=3, headwidth=10), color='red', weight='bold')
            elif tipo_svinc == 'M': 
                plt.plot(xs, ys, marker='o', markersize=20, color='red', zorder=3, fillstyle='none', markeredgewidth=3)
                plt.text(xs + 0.2, ys + 0.2, "X1", color='red', weight='bold')
            elif tipo_svinc == 'Fn':
                _a = math.radians(float(ang_X1) + 90)
                plt.annotate('X1', xy=(xs, ys), xytext=(xs - L_scale*math.cos(_a), ys - L_scale*math.sin(_a)), arrowprops=dict(facecolor='red', width=3, headwidth=10), color='red', weight='bold')
        plt.axis('equal'); plt.axis('off')
        plt.savefig(os.path.join(out_dir, percorso_S1), bbox_inches='tight', dpi=120)
    finally:
        plt.close(fig_s1)

    # Plot Diagrammi
    fig_diag, axes = plt.subplots(3, 1, figsize=(9, 14))
    titoli = ["Momento Flettente M(x)", "Taglio T(x)", "Sforzo Normale N(x)"]
    colori = ['crimson', 'forestgreen', 'dodgerblue']
    
    try:
        if X1_valore_numerico is not None:
            etichette_fatte = set()
            for idx, (ax, titolo, col) in enumerate(zip(axes, titoli, colori)):
                ax.set_title(titolo, fontsize=13, weight='bold', pad=10)
                
                for el in elementi:
                    n1, n2 = int(el['n1']), int(el['n2'])
                    x_c, y_c = [nodi[n1][0], nodi[n2][0]], [nodi[n1][1], nodi[n2][1]]
                    ax.plot(x_c, y_c, color='#adb5bd', linestyle='--', linewidth=1.5, zorder=1)
                
                for el in elementi:
                    n1, n2 = int(el['n1']), int(el['n2'])
                    x_c, y_c = [nodi[n1][0], nodi[n2][0]], [nodi[n1][1], nodi[n2][1]]
                    
                    if idx == 0: 
                        expr = M0_dict[el['id']] + X1_exact * M1_dict[el['id']] if grado_iperstaticita > 0 else M0_dict[el['id']]
                    elif idx == 1: 
                        expr = T0_dict[el['id']] + X1_exact * T1_dict[el['id']] if grado_iperstaticita > 0 else T0_dict[el['id']]
                    else: 
                        expr = N0_dict[el['id']] + X1_exact * N1_dict[el['id']] if grado_iperstaticita > 0 else N0_dict[el['id']]
                    
                    x_sym = sp.Symbol("s")
                    L = math.hypot(nodi[n2][0] - nodi[n1][0], nodi[n2][1] - nodi[n1][1])
                    if L == 0: continue
                    x_vals = np.linspace(0, L, 121)
                    
                    # FILTRO TERMICO SUI DIAGRAMMI FINALI
                    expr = sp.sympify(expr)
                    expr = expr.xreplace({n: 0.0 for n in expr.atoms(sp.Float) if abs(n) < 1e-10})
                    
                    try:
                        func = sp.lambdify(x_sym, expr, 'numpy')
                        y_vals = func(x_vals)
                        if np.isscalar(y_vals): y_vals = np.full(len(x_vals), float(y_vals))
                        y_vals = np.real(np.asarray(y_vals, dtype=complex)).astype(float)
                    except:
                        y_vals = np.array([float(sp.re(expr.subs(x_sym, v).evalf())) for v in x_vals])
                    
                    dx, dy = x_c[1] - x_c[0], y_c[1] - y_c[0]
                    ux, uy = -dy/L, dx/L
                    bx = x_c[0] + (dx/L) * x_vals
                    by = y_c[0] + (dy/L) * x_vals
                    
                    max_val = max([np.max(np.abs(y_vals)), 1e-5])
                    max_span = max([max(n[0] for n in nodi.values()) - min(n[0] for n in nodi.values()), 1])
                    scale_f = (max_span * 0.15) / max_val
                    
                    px = bx + ux * y_vals * scale_f
                    py = by + uy * y_vals * scale_f
                    
                    ax.fill(np.append(bx, px[::-1]), np.append(by, py[::-1]), col, alpha=0.25, zorder=2)
                    ax.plot(px, py, col, linewidth=2, zorder=3)
                    
                    # Punti chiave ESATTI (non campionati): estremi dell'asta, massimi/minimi interni (derivata nulla) e
                    # passaggi per lo zero. Il valore agli estremi e' quello che serve per tracciare il diagramma a mano.
                    coef_k = _coeff_poly(expr, x_sym)
                    chiavi = [(0.0, 'estremo'), (L, 'estremo')]
                    if coef_k is not None:
                        c_k = np.trim_zeros(np.array(coef_k, dtype=float), 'f')
                        if len(c_k) > 1:
                            for r_ in _radici_in(list(np.polyder(c_k)), L):
                                if 1e-9 * L < r_ < L * (1 - 1e-9): chiavi.append((float(r_), 'estremo interno'))
                            for r_ in _radici_in(list(c_k), L):
                                if 1e-9 * L < r_ < L * (1 - 1e-9): chiavi.append((float(r_), 'zero'))
                        val_k = lambda t_: float(np.polyval(c_k, t_)) if len(c_k) else 0.0
                    else:
                        chiavi.append((float(x_vals[int(np.argmax(np.abs(y_vals)))]), 'estremo interno'))
                        val_k = lambda t_: float(np.interp(t_, x_vals, y_vals))
                    soglia_lab = max(0.1, 1e-3 * max_val)
                    for t_k, tipo_k in chiavi:
                        v_k = val_k(t_k)
                        bxk = x_c[0] + (dx / L) * t_k
                        byk = y_c[0] + (dy / L) * t_k
                        pxk, pyk = bxk + ux * v_k * scale_f, byk + uy * v_k * scale_f
                        if tipo_k == 'zero':
                            ax.plot([bxk], [byk], marker='o', markersize=5, markerfacecolor='white', markeredgecolor=col, zorder=5)
                            continue
                        if abs(v_k) <= soglia_lab: continue
                        chiave_lab = (idx, round(pxk, 3), round(pyk, 3), round(v_k, 1))
                        if chiave_lab in etichette_fatte: continue
                        etichette_fatte.add(chiave_lab)
                        ax.plot([pxk], [pyk], marker='o', markersize=3.5, color=col, zorder=5)
                        sg_ = 1.0 if v_k >= 0 else -1.0
                        ax.annotate(f"{v_k:.1f}", (pxk, pyk), xytext=(ux * sg_ * 9, uy * sg_ * 9), textcoords='offset points',
                                    color=col, weight='bold', fontsize=8, ha='center', va='center',
                                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec=col, alpha=0.85), zorder=6)

                ax.set_aspect('equal')
                ax.axis('off')
                
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, percorso_diag), bbox_inches='tight', dpi=130)
    finally:
        plt.close(fig_diag)
    
    return "\n".join(report), "plv_S0_S1.png", "plv_diagrammi.png"

def risolvi_struttura_plv(dati_json, out_dir, incognita):
    """Entry point (stessa firma di prima). Gli errori bloccanti tornano come report HTML con box rosso,
    senza immagini (img_s1 = img_diag = "")."""
    ULTIMI_RISULTATI.clear()
    try:
        return _risolvi_struttura_plv_core(dati_json, out_dir, incognita)
    except ErroreSolutore as e:
        return str(e), "", ""

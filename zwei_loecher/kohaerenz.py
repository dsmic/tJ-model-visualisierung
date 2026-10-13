"""Kohaerenz und Frustration: Vorzeichen geschlossener Prozesse, Paarbewegungen, J-assistierte Schritte.

Reproduziert die Tabellen in README.md, Abschnitt 7 (Ausgabe als Markdown).
Start: python kohaerenz.py   (braucht tkinter wie tj_gui.py und numpy; Laufzeit etwa eine Minute)

Begriffe:
- Geschlossener Prozess: Endzustand = Anfangszustand. Sein Vorzeichen haengt nicht von der Fermion-Ordnung ab.
- Effektive Huepfmatrix M: M[Q, P] = Summe der Vorzeichen aller betrachteten Prozesse von Lage P nach Lage Q.
  M ist unfrustriert, wenn jede geschlossene Schleife das Produkt +1 hat; dann laesst sich M durch Umeichen
  (Vorzeichenwahl der Basiszustaende) vorzeichenfrei machen und hat denselben groessten Eigenwert wie |M|.
  Das Verhaeltnis lambda_max(M) / lambda_max(|M|) misst die Frustration (1 = keine).
"""
import random
from collections import Counter, defaultdict

import numpy as np

import zwei_loch as Z
from zwei_loch import T

S = T.site
BONDS = [tuple(sorted(b)) for b in T.BONDS]
BIDX = {b: i for i, b in enumerate(BONDS)}


def perm_sign(p):
    sg, seen = 1, set()
    for i in range(len(p)):
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            j = p[j]
            length += 1
        if length and length % 2 == 0:
            sg = -sg
    return sg


def run(conf, steps):
    """Fuehrt Vertauschungen aus: (Vorzeichen, Endkonfiguration)."""
    c, sg = list(conf), 1
    for a, b in steps:
        sg *= T.swap_sign(c, a, b)
        c[a], c[b] = c[b], c[a]
    return sg, c


def moves_to_hops(holes, seq):
    """holes: {Name: Platz}, seq: Liste (Name, Zug R/L/U/D) -> Liste der Lochspruenge und Endplaetze."""
    pos, hops = dict(holes), []
    for h, m in seq:
        x, y = divmod(pos[h], T.L)
        n = S(x + Z.DELTA[m][0], y + Z.DELTA[m][1])
        hops.append((pos[h], n))
        pos[h] = n
    return hops, pos


# --- 1. Regel fuer geschlossene Prozesse ---

def random_closed_process(conf, k, rnd):
    """k zufaellige Vertauschungen, danach gierig zurueck zum Anfang (nur benachbarte, ungleiche Inhalte)."""
    c, steps = list(conf), []
    for _ in range(k):
        a, b = rnd.choice(T.BONDS)
        if c[a] != c[b]:
            c[a], c[b] = c[b], c[a]
            steps.append((a, b))
    for _ in range(200000):
        if c == conf:
            return steps
        a, b = rnd.choice(T.BONDS)
        if c[a] != c[b] and (c[a] == conf[b] or c[b] == conf[a]) and (c[a] != conf[a] or c[b] != conf[b]):
            c[a], c[b] = c[b], c[a]
            steps.append((a, b))
    return None


def section_rule(n=300, seed=1):
    rnd = random.Random(seed)
    ok = bad = 0
    for _ in range(n):
        conf = [rnd.choice([T.UP, T.DN, T.HOLE]) for _ in range(T.N)]
        steps = random_closed_process(conf, rnd.randint(2, 40), rnd)
        if steps is None:
            continue
        ident, c, s = list(range(T.N)), list(conf), 0
        for a, b in steps:
            s += T.HOLE not in (c[a], c[b])
            c[a], c[b] = c[b], c[a]
            ident[a], ident[b] = ident[b], ident[a]
        holes = [i for i in range(T.N) if conf[i] == T.HOLE]
        pred = (-1) ** s * perm_sign([holes.index(ident[h]) for h in holes])
        ok += T.apply_steps(conf, steps)[1] == pred
        bad += T.apply_steps(conf, steps)[1] != pred
    print("## 7.1 Regel für geschlossene Prozesse\n")
    print(f"Vorzeichen = (−1)^s · sgn(π_Loch) in {ok} von {ok + bad} geschlossenen Zufallsprozessen "
          f"(beliebige Startkonfigurationen, Seed {seed}); falsch: {bad}\n")


# --- 2. geschlossene Paar- und Lochprozesse ---

def caterpillar_loop(a, M):
    """A laeuft die geschlossene Schleife M, B steht auf dem letzten Schleifenplatz und folgt direkt."""
    pa = Z.path_sites(a, M)
    b = pa[-2]
    seq = []
    mb = M[-1] + M[:-1]
    for k in range(len(M)):
        seq += [("A", M[k]), ("B", mb[k])]
    return {"A": a, "B": b}, seq


def closed_row(name, holes, seq):
    hops, pos = moves_to_hops(holes, seq)
    start = Z.neel_with_holes(list(holes.values()))
    assert Z.legal(start, hops), name
    sg, end = run(start, hops)
    swapped = [h for h in holes if pos[h] != holes[h]]
    found, _, cut = T.search_signs(start, hops, cap=300000)
    reach = " und ".join(f"{s:+d}" for s in sorted(found)) + (" (Suche abgebrochen)" if cut else "")
    return (f"| {name} | {'ja' if end == start else 'nein'} | "
            f"{'ja' if swapped else 'nein'} | {f'{sg:+d}' if end == start else '–'} | {reach} |")


def section_closed():
    print("## 7.2 Geschlossene Prozesse von Paaren und Löchern\n")
    print("Volle Néel-Konfiguration mit nur den beteiligten Löchern. „ohne Reparatur“ = nur die Lochsprünge "
          "(Vorzeichen nur angegeben, wenn das schon ein geschlossener Prozess ist); „erreichbar“ = alle Vorzeichen, "
          "die mit beliebig eingestreuten Reparaturen bis zum Anfangszustand erreichbar sind.\n")
    print("| Prozess | Endzustand ohne Reparatur = Anfang | Löcher vertauscht | Vorzeichen ohne Reparatur "
          "| erreichbar |\n|---|---|---|---|---|")
    rows = []
    for M in ["RDLU", "RRDLLU", "RRDDLLUU", "RRDDLULU", "RRRDDDLLLUUU"]:
        for a in (S(3, 2), S(3, 3)):
            holes, seq = caterpillar_loop(a, M)
            rows.append(closed_row(f"Raupe, Umlauf {M}, Start Untergitter {Z.sublattice(a)}", holes, seq))
    # Wende: Kopf und Schwanz tauschen die Plaetze
    rows.append(closed_row("Raupe mit Wende (A: URDLL, B: RURDL)", {"A": S(3, 4), "B": S(3, 3)},
                           [("A", "U"), ("B", "R"), ("A", "R"), ("B", "U"), ("A", "D"), ("B", "R"),
                            ("A", "L"), ("B", "D"), ("A", "L"), ("B", "L")]))
    ring = "RRDDLLUU"
    p = Z.path_sites(S(2, 2), ring)
    holes = {"A1": p[0], "B1": p[7], "A2": p[4], "B2": p[3]}
    seq = []
    for k in range(4):
        seq += [("A1", ring[k]), ("A2", ring[(4 + k) % 8]), ("B1", ring[(7 + k) % 8]), ("B2", ring[(3 + k) % 8])]
    rows.append(closed_row("zwei Raupenpaare tauschen die Plätze", holes, seq))
    r1, r2 = "RRDDLLUU", "LLUURRDD"
    q1, q2 = Z.path_sites(S(1, 1), r1), Z.path_sites(S(4, 4), r2)
    holes = {"A1": q1[0], "B1": q1[7], "A2": q2[0], "B2": q2[7]}
    seq = []
    for k in range(8):
        seq += [("A1", r1[k]), ("B1", r1[(7 + k) % 8])]
    for k in range(8):
        seq += [("A2", r2[k]), ("B2", r2[(7 + k) % 8])]
    rows.append(closed_row("zwei Raupenpaare, Umläufe kreuzen sich", holes, seq))
    seq = []
    for k in range(4):
        seq += [("A", ring[k]), ("B", ring[(4 + k) % 8])]
    rows.append(closed_row("zwei einzelne Löcher tauschen die Plätze", {"A": p[0], "B": p[4]}, seq))
    rows.append(closed_row("zwei einzelne Löcher, Umläufe kreuzen sich", {"A": q1[0], "B": q2[0]},
                           [("A", m) for m in r1] + [("B", m) for m in r2]))
    rows.append(closed_row("ein Loch, Umlauf RDLU", {"A": S(3, 3)}, [("A", m) for m in "RDLU"]))
    print("\n".join(rows) + "\n")


# --- 3. effektive Paarbewegung ---

def pair_matrices():
    """Raupe (s=0), Pivot und Domino (s=1, Reparatur danach) auf allen Bindungen des 8x8-Bretts."""
    nb = len(BONDS)
    Mcat, Mj = np.zeros((nb, nb)), np.zeros((nb, nb))
    for (u, v) in BONDS:
        conf = Z.neel_with_holes([u, v])
        for lead, tail in ((u, v), (v, u)):
            for n in T.neighbours(lead):
                if n == tail:
                    continue
                sg, c = run(conf, [(lead, n), (tail, lead)])           # Raupe
                assert not T.neel_violations(c)
                Mcat[BIDX[tuple(sorted((lead, n)))], BIDX[(u, v)]] += sg
                for m in T.neighbours(tail):                           # Pivot um lead: tail -> n ueber die Ecke m
                    if m == lead or m not in T.neighbours(n):
                        continue
                    sg, c = run(conf, [(tail, m), (m, n), (tail, m)])
                    assert not T.neel_violations(c)
                    Mj[BIDX[tuple(sorted((lead, n)))], BIDX[(u, v)]] += sg
        (xu, yu), (xv, yv) = divmod(u, T.L), divmod(v, T.L)
        for e in ([(1, 0), (-1, 0)] if xu == xv else [(0, 1), (0, -1)]):    # Domino: beide seitlich
            a2, b2 = (xu + e[0], yu + e[1]), (xv + e[0], yv + e[1])
            if not all(0 <= q < T.L for q in a2 + b2):
                continue
            ua, vb = S(*a2), S(*b2)
            sgs = {run(conf, order + [(u, v)])[0] for order in ([(u, ua), (v, vb)], [(v, vb), (u, ua)])}
            assert len(sgs) == 1
            Mj[BIDX[tuple(sorted((ua, vb)))], BIDX[(u, v)]] += sgs.pop()
    return Mcat, Mj


def delta_signs():
    """Vorzeichen von Delta(b)|Neel>, Delta(i,j) = c_i,up c_j,dn - c_i,dn c_j,up (Singulett auf Bindung b)."""
    neel = Z.neel_with_holes([])
    out = []
    for i, j in BONDS:
        for si, sj, pref in ((T.UP, T.DN, 1), (T.DN, T.UP, -1)):
            if neel[i] == si and neel[j] == sj:
                occ = T.to_occ(neel)
                out.append(pref * T.apply_c(occ, T.orbital(j, sj), False) * T.apply_c(occ, T.orbital(i, si), False))
    return np.array(out)


def lam_max(M):
    return np.linalg.eigvalsh(M)[-1]


def section_pairs():
    Mcat, Mj = pair_matrices()
    sd = delta_signs()
    horiz = np.array([1.0 if u // T.L == v // T.L else -1.0 for u, v in BONDS])
    w, V = np.linalg.eigh(Mj)
    D = np.sign(V[:, -1])
    nb = len(BONDS)
    print("## 7.3 Effektive Paarbewegung (Paar aus Nachbarn, Néel-Hintergrund, 8×8)\n")
    print("| Schrittarten | λ_max(M) | λ_max(\\|M\\|) | Verhältnis | Entartung |\n|---|---|---|---|---|")
    for name, M in (("nur Raupe (s = 0)", Mcat), ("nur Pivot + Domino (s = 1)", Mj)):
        wm = np.linalg.eigvalsh(M)
        deg = int(np.sum(np.abs(wm - wm[-1]) < 1e-8))
        print(f"| {name} | {wm[-1]:.3f} | {lam_max(np.abs(M)):.3f} | {wm[-1] / lam_max(np.abs(M)):.3f} | {deg} |")
    print()
    G = D[:, None] * Mcat * D[None, :]
    straight = np.array([[horiz[i] == horiz[j] for j in range(nb)] for i in range(nb)])
    nz = Mcat != 0
    shared = nz & (Mj != 0)
    print(f"- Pivot und Domino werden vorzeichenfrei mit der Eichung D = sgn(Perron-Vektor); "
          f"D · Δ-Vorzeichen · (+1 waagerecht, −1 senkrecht) ist konstant: "
          f"{'ja' if len(set(D * sd * horiz)) == 1 else 'nein'} (d-Wellen-Muster).")
    print(f"- Übergänge, die Raupe mit Wende und Pivot gemeinsam haben: {int(shared.sum())}, "
          f"davon gleiches Vorzeichen: {int((np.sign(Mcat[shared]) == np.sign(Mj[shared])).sum())}.")
    print(f"- In dieser Eichung: Raupe mit Wende positiv {int((G[nz & ~straight] > 0).sum())}, "
          f"negativ {int((G[nz & ~straight] < 0).sum())}; gerade Raupenschritte positiv "
          f"{int((G[nz & straight] > 0).sum())}, negativ {int((G[nz & straight] < 0).sum())}.\n")
    print("Gemischt: M = Raupe + λ · (Pivot + Domino). Überlapp des obersten Zustands mit dem s- bzw. "
          "d-Wellen-Paar Σ_b f(b) Δ(b)|Néel⟩ (f = 1 bzw. f = ±1 für waagerecht/senkrecht), normiert:\n")
    print("| λ | λ_max(M) / λ_max(\\|M\\|) | Entartung | Überlapp s | Überlapp d |\n|---|---|---|---|---|")
    for lam in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
        M = Mcat + lam * Mj
        wm, Vm = np.linalg.eigh(M)
        deg = int(np.sum(np.abs(wm - wm[-1]) < 1e-8))
        ov = np.mean([((sd * Vm[:, -1 - k]).sum() ** 2 / nb, (horiz * sd * Vm[:, -1 - k]).sum() ** 2 / nb)
                      for k in range(deg)], axis=0)
        print(f"| {lam} | {wm[-1] / lam_max(np.abs(M)):.3f} | {deg} | {ov[0]:.3f} | {ov[1]:.3f} |")
    print()


# --- 4. einzelnes Loch: J-assistierte Schritte ---

def single_processes(h):
    """Alle Prozesse aus 2 Lochspruengen und 1 Spin-Vertauschung (vor, zwischen oder nach den Spruengen),
    die wieder bei Neel enden: {Zielplatz: [(Vorzeichen, Art)]}."""
    conf = Z.neel_with_holes([h])
    out = defaultdict(list)
    kinds = {0: "vorher", 1: "dazwischen", 2: "nachher"}
    for a in T.neighbours(h):
        for b in T.neighbours(a):
            if b == h:
                continue
            hops = [(h, a), (a, b)]
            for pos in range(3):
                c0 = run(conf, hops[:pos])[1]
                for p, q in T.BONDS:
                    if T.HOLE in (c0[p], c0[q]) or c0[p] == c0[q]:
                        continue
                    sg, end = run(conf, hops[:pos] + [(p, q)] + hops[pos:])
                    if end[b] == T.HOLE and not T.neel_violations(end):
                        out[b].append((sg, kinds[pos]))
    return out


def section_single():
    h = S(3, 3)
    print("## 7.4 Ein Loch: Schritt aus 2 Sprüngen und 1 Vertauschung\n")
    print("| Ziel relativ zum Loch | Prozesse | vorher (berührt künftigen Weg) | dazwischen | nachher (Reparatur) "
          "| Vorzeichen |\n|---|---|---|---|---|---|")
    for b, lst in sorted(single_processes(h).items()):
        k = Counter(kind for _, kind in lst)
        d = (b // T.L - 3, b % T.L - 3)
        print(f"| {d} | {len(lst)} | {k['vorher']} | {k['dazwischen']} | {k['nachher']} "
              f"| {', '.join(f'{s:+d}' for s in sorted(set(s for s, _ in lst)))} |")
    sites = [s for s in range(T.N) if Z.sublattice(s) == Z.sublattice(h)]
    ix = {s: i for i, s in enumerate(sites)}
    M = np.zeros((len(sites), len(sites)))
    for s in sites:
        for b, lst in single_processes(s).items():
            M[ix[b], ix[s]] += sum(sg for sg, _ in lst)
    loops = {"Dreieck (diagonal, diagonal, gerade zurück)": [S(3, 3), S(4, 4), S(5, 3), S(3, 3)],
             "Raute aus 4 Diagonalschritten": [S(3, 3), S(4, 4), S(5, 3), S(4, 2), S(3, 3)]}
    print()
    for name, loop in loops.items():
        p = np.prod([np.sign(M[ix[loop[i + 1]], ix[loop[i]]]) for i in range(len(loop) - 1)])
        print(f"- {name}: Produkt der Vorzeichen {int(p):+d}")
    w = np.linalg.eigvalsh(M)
    print(f"- Hüpfmatrix auf dem Untergitter: λ_max(M) = {w[-1]:.3f}, λ_max(|M|) = {lam_max(np.abs(M)):.3f}, "
          f"λ_min(M) = {w[0]:.3f} (Spektrum gespiegelt: M ist bis auf Umeichen −|M|)\n")


if __name__ == "__main__":
    section_rule()
    section_closed()
    section_pairs()
    section_single()

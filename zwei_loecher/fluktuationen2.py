"""Paar mit bis zu ZWEI J-Fluktuationen: vollstaendige Prozessaufzaehlung wie in fluktuationen.py.

Reproduziert README.md, Abschnitt 7.8 (Ausgabe als Markdown).
Start: python fluktuationen2.py   (braucht tkinter wie tj_gui.py, numpy und scipy; Laufzeit je nach Rechner
       etwa 10-30 Minuten, nutzt alle Prozessorkerne)

Zustaende: Paar (zwei benachbarte Loecher) im 6x6-Fenster, Hintergrund Neel mit 0, 1 oder 2 Fluktuationen
(umgeklappten Nachbarpaaren), jede hoechstens DMAX vom Paar entfernt. Prozesse: alle Folgen aus 1 bis 3
Elementarschritten zwischen zwei Zustaenden ohne Zwischenstopp in einem Zustand.
"""
import itertools
import multiprocessing as mp
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

import fluktuationen as F

N, WIN, NB, DMAX = F.N, F.WIN, F.NB, F.DMAX
MAXFLUCT = 2


def matched_pairs(w, up):
    """Zerlegung der falschen Spins in benachbarte, entgegengesetzte Paare (oder None)."""
    if not w:
        return []
    a = w[0]
    for b in w[1:]:
        if F.dist(a, b) == 1 and ((up >> a & 1) != (up >> b & 1)):
            rest = matched_pairs([x for x in w if x not in (a, b)], up)
            if rest is not None:
                return [(a, b)] + rest
    return None


def kind(up, dn):
    """None (kein Zustand) oder Zahl der Fluktuationen (0 bis MAXFLUCT)."""
    h = F.holes(up, dn)
    if len(h) != 2 or F.dist(*h) != 1 or not set(h) <= WIN:
        return None
    w = F.wrong(up, dn)
    if len(w) > 2 * MAXFLUCT or len(w) % 2 or not set(w) <= WIN:
        return None
    pairs = matched_pairs(w, up)
    if pairs is None:
        return None
    if any(min(F.dist(x, y) for x in p for y in h) > DMAX for p in pairs):
        return None
    return len(pairs)


def build_states():
    bonds = [(a, b) for a in sorted(WIN) for b in NB[a] if b in WIN and a < b]
    states = set()
    for a, b in bonds:
        up, dn = F.NEEL_UP & ~(1 << a) & ~(1 << b), F.NEEL_DN & ~(1 << a) & ~(1 << b)
        states.add((up, dn))
        near = [(p, q) for p, q in bonds if not {p, q} & {a, b}
                and min(F.dist(x, y) for x in (p, q) for y in (a, b)) <= DMAX]
        for k in range(1, MAXFLUCT + 1):
            for combo in itertools.combinations(near, k):
                if len({s for bd in combo for s in bd}) != 2 * k:
                    continue
                u, d = up, dn
                for p, q in combo:
                    _, u, d = F.swap(u, d, p, q)
                if kind(u, d) == k:
                    states.add((u, d))
    return sorted(states)


IDX = {}


def _init(states):
    global IDX
    IDX = {s: i for i, s in enumerate(states)}


def _process_one(item):
    i, (up, dn) = item
    out = []
    stack = [(up, dn, 1, 0, 0, 0)]
    while stack:
        u, d, sg, m, s, depth = stack.pop()
        if depth == 3:
            continue
        for typ, sg1, u2, d2 in F.steps(u, d):
            m2, s2 = m + (typ == "h"), s + (typ == "s")
            j = IDX.get((u2, d2))
            if j is not None:
                out.append((j, i, sg * sg1, m2, s2))
            elif len(F.wrong(u2, d2)) <= 2 * MAXFLUCT + 2 * (2 - depth):
                stack.append((u2, d2, sg * sg1, m2, s2, depth + 1))
    return out


def enumerate_processes(states, workers=None):
    with mp.Pool(workers or mp.cpu_count(), initializer=_init, initargs=(states,)) as pool:
        chunks = pool.map(_process_one, list(enumerate(states)), chunksize=50)
    return [p for c in chunks for p in c]


def build(procs, n, keep, ws, sub, randomize=False):
    """Prozessmatrix; randomize: Übergänge mit Lochsprüngen bekommen ein zufälliges Vorzeichen (Gegenprobe)."""
    rng = np.random.default_rng(1)
    flip = {}
    rows, cols, vals = [], [], []
    for j, i, sg, m, s in procs:
        if sub[i] and sub[j] and keep(m, s):
            if randomize and m > 0:
                key = (min(i, j), max(i, j))
                if key not in flip:
                    flip[key] = rng.choice([-1, 1])
                sg = flip[key] * abs(sg)
            rows.append(j)
            cols.append(i)
            vals.append(sg * ws ** s)
    M = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    M.sum_duplicates()
    used = np.where((abs(M).sum(axis=0).A1 > 0) & sub)[0]
    return M[used][:, used], used


def top(M):
    w, v = sla.eigsh(M, k=1, which="LA", tol=1e-10)
    return w[0], v[:, 0]


def main():
    t0 = time.time()
    states = build_states()
    kinds = np.array([kind(*s) for s in states])
    n = len(states)
    print("## 7.8 Zwei Fluktuationen\n")
    print(f"{n} Zustände: " + ", ".join(f"{int(np.sum(kinds == k))} mit {k} Fluktuation(en)" for k in range(MAXFLUCT + 1))
          + f" ({time.time() - t0:.0f} s)")
    procs = enumerate_processes(states)
    print(f"{len(procs)} Prozesse aufgezählt ({time.time() - t0:.0f} s).")
    C = sp.csr_matrix(([1] * len(procs), ([p[0] for p in procs], [p[1] for p in procs])), shape=(n, n))
    S_ = sp.csr_matrix(([p[2] for p in procs], ([p[0] for p in procs], [p[1] for p in procs])), shape=(n, n))
    asym = int(((abs(C - C.T) + abs(S_ - S_.T)) > 0).sum())
    print(f"Symmetrieprüfung: {asym} Einträge ohne passenden Rückweg.\n")
    variants = [("Raupe, sauberes Néel", lambda m, s: s == 0, 0, False),
                ("Raupe + freie J-Fluktuationen, höchstens 1", lambda m, s: s == 0 or m == 0, 1, False),
                ("Raupe + freie J-Fluktuationen, höchstens 2", lambda m, s: s == 0 or m == 0, 2, False),
                ("Gegenprobe: höchstens 2, zufällige Vorzeichen der Sprungprozesse", lambda m, s: s == 0 or m == 0, 2, True),
                ("alle Prozesse, höchstens 1 Fluktuation", lambda m, s: True, 1, False),
                ("alle Prozesse, höchstens 2 Fluktuationen", lambda m, s: True, 2, False)]
    print("| Prozesse | Gewicht Vertauschung | Verhältnis λ_max(M)/λ_max(\\|M\\|) | Anteil 0 / 1 / 2 Fluktuationen |")
    print("|---|---|---|---|")
    for ws in (1.0, 0.3, 0.1):
        for name, keep, kmax, randomize in variants:
            M, used = build(procs, n, keep, ws, kinds <= kmax, randomize)
            a, v = top(M)
            b, _ = top(abs(M))
            share = " / ".join(f"{np.sum(v[kinds[used] == k] ** 2):.2f}" for k in range(MAXFLUCT + 1))
            print(f"| {name} | {ws} | {a / b:.3f} | {share} |")
    print(f"\nGesamtlaufzeit {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

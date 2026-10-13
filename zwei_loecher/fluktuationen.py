"""Paarbewegung mit gelockertem Neel: Vertauschungen, die keine Reparaturen sind (J-Fluktuationen).

Reproduziert README.md, Abschnitt 7.6 (Ausgabe als Markdown).
Start: python fluktuationen.py   (braucht tkinter wie tj_gui.py und numpy; Laufzeit etwa 5 Minuten)

Zustaende: zwei Loecher auf Nachbarplaetzen (Paar) im 6x6-Fenster des 8x8-Bretts, Hintergrund Neel oder Neel
mit genau einer Fluktuation (benachbartes, umgeklapptes Spinpaar) hoechstens DMAX vom Paar entfernt.
Prozesse: alle Folgen aus 1 bis 3 Elementarschritten (Lochsprung, Vertauschung benachbarter ungleicher Spins),
die zwischen zwei Zustaenden laufen, ohne unterwegs einen Zustand zu treffen. Vertauschungen nur in der Naehe
von Loechern oder falschen Spins. Vorzeichen wie tj_gui.swap_sign, ueber Bitmasken berechnet (zu Beginn geprueft).
"""
import random
import time
from collections import Counter, defaultdict, deque

import numpy as np

import kohaerenz as K
import zwei_loch as Z
from zwei_loch import T

L, N = T.L, T.N
WIN = {T.site(x, y) for x in range(1, 7) for y in range(1, 7)}
NB = {s: T.neighbours(s) for s in range(N)}
SUB = [(s // L + s % L) % 2 for s in range(N)]
NEEL_UP = sum(1 << s for s in range(N) if SUB[s] == 0)
NEEL_DN = sum(1 << s for s in range(N) if SUB[s] == 1)
DMAX = 2


def dist(a, b):
    return abs(a // L - b // L) + abs(a % L - b % L)


def csign(o, k):
    return -1 if bin(o & ((1 << k) - 1)).count("1") & 1 else 1


def hop(up, dn, src, dst):
    """Spin von src auf das Loch dst (c^dag_dst c_src): (Vorzeichen, up, dn)."""
    is_up = up >> src & 1
    k1, k2 = (N + src, N + dst) if is_up else (src, dst)
    o = dn | (up << N)
    s1 = csign(o, k1)
    o ^= 1 << k1
    s2 = csign(o, k2)
    if is_up:
        return s1 * s2, up ^ (1 << src) ^ (1 << dst), dn
    return s1 * s2, up, dn ^ (1 << src) ^ (1 << dst)


def swap(up, dn, a, b):
    """Vertauschung ungleicher Spins auf a und b, Operatorfolge wie tj_gui.swap_sign."""
    va, vb = (1 if up >> a & 1 else -1), (1 if up >> b & 1 else -1)
    o, sg = dn | (up << N), 1
    for k, create in (((N if vb == 1 else 0) + b, False), ((N if vb == 1 else 0) + a, True),
                      ((N if va == 1 else 0) + a, False), ((N if va == 1 else 0) + b, True)):
        sg *= csign(o, k)
        o ^= 1 << k
    return sg, up ^ (1 << a) ^ (1 << b), dn ^ (1 << a) ^ (1 << b)


def holes(up, dn):
    return [s for s in range(N) if not ((up | dn) >> s & 1)]


def wrong(up, dn):
    return [s for s in range(N) if (up >> s & 1 and SUB[s] == 1) or (dn >> s & 1 and SUB[s] == 0)]


def to_masks(conf):
    return (sum(1 << s for s in range(N) if conf[s] == T.UP), sum(1 << s for s in range(N) if conf[s] == T.DN))


def check_signs(n=3000, seed=3):
    rnd, ok, bad = random.Random(seed), 0, 0
    for _ in range(n):
        conf = [rnd.choice([T.UP, T.DN, T.HOLE]) for _ in range(N)]
        a, b = rnd.choice(T.BONDS)
        if conf[a] == conf[b]:
            continue
        up, dn = to_masks(conf)
        if T.HOLE in (conf[a], conf[b]):
            src, dst = (a, b) if conf[b] == T.HOLE else (b, a)
            got = hop(up, dn, src, dst)[0]
        else:
            got = swap(up, dn, a, b)[0]
        ok += got == T.swap_sign(conf, a, b)
        bad += got != T.swap_sign(conf, a, b)
    return ok, bad


def kind(up, dn):
    """None (kein Zustand), 0 (sauberes Neel) oder 1 (eine Fluktuation nahe am Paar)."""
    h = holes(up, dn)
    if len(h) != 2 or dist(*h) != 1 or not set(h) <= WIN:
        return None
    w = wrong(up, dn)
    if not w:
        return 0
    if len(w) == 2 and dist(*w) == 1 and set(w) <= WIN and min(dist(a, b) for a in w for b in h) <= DMAX:
        return 1
    return None


def build_states():
    bonds = [(a, b) for a in sorted(WIN) for b in NB[a] if b in WIN and a < b]
    states = []
    for a, b in bonds:
        up, dn = NEEL_UP & ~(1 << a) & ~(1 << b), NEEL_DN & ~(1 << a) & ~(1 << b)
        states.append((up, dn))
        for p, q in bonds:
            if not {p, q} & {a, b}:
                _, u2, d2 = swap(up, dn, p, q)
                if kind(u2, d2) == 1:
                    states.append((u2, d2))
    return states


def steps(up, dn):
    """Elementarschritte: alle Lochsprünge im Fenster und Vertauschungen nahe an den Löchern. Die Auswahl hängt nur
    von den Lochplätzen ab; die ändern sich durch eine Vertauschung nicht, also werden Hin- und Rückweg gleich
    behandelt und die Prozessmatrix ist symmetrisch."""
    focus = set(holes(up, dn))
    for x in holes(up, dn):
        for n in NB[x]:
            if n in WIN and (up | dn) >> n & 1:
                yield "h", *hop(up, dn, n, x)
    for p in WIN:
        if min(dist(p, f) for f in focus) > DMAX + 1:
            continue
        for q in NB[p]:
            if q in WIN and p < q and ((up >> p & 1 and dn >> q & 1) or (dn >> p & 1 and up >> q & 1)):
                yield "s", *swap(up, dn, p, q)


def enumerate_processes(states):
    idx = {s: i for i, s in enumerate(states)}
    procs = defaultdict(list)            # (Ziel, Start) -> [(Vorzeichen, m, s)]
    for i, (up, dn) in enumerate(states):
        stack = [(up, dn, 1, 0, 0, 0)]
        while stack:
            u, d, sg, m, s, depth = stack.pop()
            if depth == 3:
                continue
            for typ, sg1, u2, d2 in steps(u, d):
                m2, s2 = m + (typ == "h"), s + (typ == "s")
                if (u2, d2) in idx:
                    procs[(idx[(u2, d2)], i)].append((sg * sg1, m2, s2))
                elif len(wrong(u2, d2)) <= 2 * (2 - depth) + 2:
                    stack.append((u2, d2, sg * sg1, m2, s2, depth + 1))
    return procs


def ratio(M):
    used = np.abs(M).sum(axis=0) > 0
    M = M[np.ix_(used, used)]
    return np.linalg.eigvalsh(M)[-1] / np.linalg.eigvalsh(np.abs(M))[-1], M.shape[0]


def describe(states, k):
    up, dn = states[k]
    f = lambda s: f"({s // L},{s % L})"
    w = wrong(up, dn)
    return f"Paar {' '.join(f(s) for s in holes(up, dn))}, " + (f"Fluktuation {' '.join(f(s) for s in w)}" if w else "sauber")


def heal_signs(conf, sign, region, depth=5):
    """Alle Vorzeichen, mit denen reine Spin-Vertauschungen (in region) conf bis Neel heilen
    (Breitensuche; 'Widerspruch', falls derselbe Zustand mit zwei Vorzeichen erreicht wird)."""
    seen, out = {tuple(conf): sign}, set()
    queue = deque([(list(conf), sign, 0)])
    while queue:
        c, sg, d = queue.popleft()
        if not T.neel_violations(c):
            out.add(sg)
            continue
        if d == depth:
            continue
        for a, b in T.BONDS:
            if a in region and b in region and T.HOLE not in (c[a], c[b]) and c[a] != c[b]:
                n2 = list(c)
                s2 = sg * T.swap_sign(c, a, b)
                n2[a], n2[b] = n2[b], n2[a]
                key = tuple(n2)
                if key in seen:
                    if seen[key] != s2:
                        out.add("Widerspruch")
                    continue
                seen[key] = s2
                queue.append((n2, s2, d + 1))
    return out


def section_path_fluctuation():
    """Weg UU eines Lochs bzw. eines Raupenpaars, mit und ohne vorher umgeklapptes Nachbarpaar."""
    S = T.site
    h, b, p1, p2 = S(5, 3), S(6, 3), S(4, 3), S(3, 3)
    region = {S(x, y) for x in range(2, 8) for y in range(1, 6)}
    single, pair = [(h, p1), (p1, p2)], [(h, p1), (b, h), (p1, p2), (h, p1)]

    def signs(holes_, hops, flip):
        conf = Z.neel_with_holes(holes_)
        c, sg = list(conf), 1
        for a, x in ([flip] if flip else []) + hops:
            sg *= T.swap_sign(c, a, x)
            c[a], c[x] = c[x], c[a]
        return heal_signs(c, sg, region)
    fmt = lambda st: " und ".join(f"{s:+d}" if isinstance(s, int) else s for s in sorted(st, key=str))
    print("Weg UU ab (5,3); vorher wird ein richtiges Nachbarpaar umgeklappt (keine Reparatur). Danach wird mit "
          "beliebigen reinen Spin-Vertauschungen bis Néel geheilt; angegeben sind alle dabei erreichbaren Vorzeichen.\n")
    print("| umgeklapptes Paar | ein Loch | Raupenpaar (B von (6,3) folgt) |\n|---|---|---|")
    for name, flip in (("keines", None), ("U mit UR (berührt den Weg)", (p1, S(4, 4))),
                       ("U mit UL (berührt den Weg)", (p1, S(4, 2))), ("UU mit UUR (Zielplatz)", (p2, S(3, 4))),
                       ("R mit UR (neben dem Weg)", (S(5, 4), S(4, 4)))):
        print(f"| {name} | {fmt(signs([h], single, flip))} | {fmt(signs([h, b], pair, flip))} |")
    print()


def main():
    print("## 7.6 Gelockertes Néel: J-Fluktuationen\n")
    section_path_fluctuation()
    ok, bad = check_signs()
    print(f"Bitmasken-Vorzeichen gegen tj_gui.swap_sign: {ok} gleich, {bad} verschieden.\n")
    rnd, ok, bad = random.Random(2), 0, 0
    for _ in range(300):
        conf = [rnd.choice([T.UP, T.DN, T.HOLE]) for _ in range(N)]
        st = K.random_closed_process(conf, rnd.randint(2, 40), rnd)
        if st is None:
            continue
        ident, c, nup = list(range(N)), list(conf), 0
        for a, b in st:
            nup += T.HOLE in (c[a], c[b]) and T.UP in (c[a], c[b])
            c[a], c[b] = c[b], c[a]
            ident[a], ident[b] = ident[b], ident[a]
        hl = [i for i in range(N) if conf[i] == T.HOLE]
        pred = (-1) ** nup * K.perm_sign([hl.index(ident[h]) for h in hl])
        ok += T.apply_steps(conf, st)[1] == pred
        bad += T.apply_steps(conf, st)[1] != pred
    print(f"Gleichwertige Form der Regel 7.1: Vorzeichen = (−1)^(n↑) · sgn(π_Loch), n↑ = Zahl der Lochsprünge "
          f"mit einem ↑-Spin: {ok} von {ok + bad} geschlossenen Zufallsprozessen.\n")
    t0 = time.time()
    states = build_states()
    procs = enumerate_processes(states)
    kinds = np.array([kind(*s) for s in states])
    print(f"{len(states)} Zustände ({int(np.sum(kinds == 0))} ohne Fluktuation), {len(procs)} Übergänge, "
          f"Aufzählung {time.time() - t0:.0f} s.")
    asym = sum(1 for (j, i), lst in procs.items() if len(procs.get((i, j), [])) != len(lst)
               or sum(p[0] for p in procs.get((i, j), [])) != sum(p[0] for p in lst))
    print(f"Symmetrieprüfung: {asym} Übergänge ohne passenden Rückweg.\n")
    n = len(states)
    clean = np.where(kinds == 0)[0]

    def matrix(keep, ws):
        M = np.zeros((n, n))
        for (j, i), lst in procs.items():
            for sg, m, s in lst:
                if keep(i, j, m, s):
                    M[j, i] += sg * ws ** s
        return M
    rows = [("sauberes Néel, alle Prozesse", lambda i, j, m, s: True, clean),
            ("sauberes Néel, nur Prozesse mit Vertauschung", lambda i, j, m, s: s >= 1, clean),
            ("mit Fluktuation, nur Prozesse mit Vertauschung", lambda i, j, m, s: s >= 1, None),
            ("mit Fluktuation, alle Prozesse", lambda i, j, m, s: True, None)]
    print("| Zustandsraum, Prozesse | Verhältnis (Gewicht Vertauschung 1) | Verhältnis (Gewicht 0,3) | Zustände |")
    print("|---|---|---|---|")
    for name, keep, sub in rows:
        out = []
        for ws in (1.0, 0.3):
            M = matrix(keep, ws)
            if sub is not None:
                M = M[np.ix_(sub, sub)]
            out.append(ratio(M))
        print(f"| {name} | {out[0][0]:.3f} | {out[1][0]:.3f} | {out[0][1]} |")
    print()
    print("Raupe allein und Raupe mit freien J-Fluktuationen (Prozesse ohne Vertauschung oder ohne Lochsprung):\n")
    print("Gegenprobe: dieselben Prozesse, aber jeder Übergang mit Lochsprüngen bekommt ein zufälliges Vorzeichen "
          "(gleich für Hin- und Rückweg). Ist das Verhältnis dann ähnlich, kommt es nur von der Verdünnung durch "
          "reine Spinprozesse.\n")
    print("| Prozesse | Gewicht Vertauschung | Verhältnis | Entartung | Anteil des Grundzustands mit Fluktuation |")
    print("|---|---|---|---|---|")
    rnd = np.random.default_rng(1)
    flip = {}

    def random_matrix(keep, ws):
        M = np.zeros((n, n))
        for (j, i), lst in procs.items():
            for sg, m, s in lst:
                if keep(i, j, m, s):
                    if m > 0:
                        key = (min(i, j), max(i, j))
                        flip.setdefault(key, rnd.choice([-1, 1]))
                        sg = flip[key] * abs(sg)
                    M[j, i] += sg * ws ** s
        return M
    for ws in (1.0, 0.3, 0.1):
        for name, keep, sub in (("Raupe, sauberes Néel", lambda i, j, m, s: s == 0, clean),
                                ("Raupe + freie J-Fluktuationen", lambda i, j, m, s: s == 0 or m == 0, None),
                                ("Gegenprobe: zufällige Vorzeichen der Sprungprozesse",
                                 lambda i, j, m, s: s == 0 or m == 0, "zufall")):
            if isinstance(sub, str):
                M, sub = random_matrix(keep, ws), None
            else:
                M = matrix(keep, ws)
            idx = np.arange(n) if sub is None else sub
            M = M[np.ix_(idx, idx)]
            used = np.abs(M).sum(axis=0) > 0
            M, idx = M[np.ix_(used, used)], idx[used]
            w, V = np.linalg.eigh(M)
            deg = int(np.sum(np.abs(w - w[-1]) < 1e-8))
            frac = float(np.sum(V[:, -1] ** 2 * (kinds[idx] == 1))) if deg == 1 else float("nan")
            print(f"| {name} | {ws} | {w[-1] / np.linalg.eigvalsh(np.abs(M))[-1]:.3f} | {deg} | "
                  f"{'–' if np.isnan(frac) or sub is not None else f'{frac:.2f}'} |")
    print()
    amp = {}
    for (j, i), lst in procs.items():
        a = sum(sg for sg, m, s in lst if s >= 1)
        if a:
            amp[(i, j)] = int(np.sign(a))
    nbr = defaultdict(dict)
    for (i, j), sg in amp.items():
        nbr[i][j] = sg
    for i in clean:
        for j, s1 in nbr[i].items():
            for k, s2 in nbr[j].items():
                if k != i and i in nbr[k] and kinds[j] == 1 and kinds[k] == 1 and holes(*states[j]) == holes(*states[k]) \
                        and s1 * s2 * nbr[k][i] < 0:
                    print("Kleinste frustrierte Schleife (Beispiel), Produkt −1:\n")
                    for a, b in ((i, j), (j, k), (k, i)):
                        c = Counter((m, s) for sg, m, s in procs[(b, a)] if s >= 1)
                        typ = ", ".join(f"{cnt}× ({m} Sprünge, {s} Vertauschung{'en' if s != 1 else ''})"
                                        for (m, s), cnt in sorted(c.items()))
                        print(f"- {describe(states, a)} → {describe(states, b)}: {typ}")
                    print()
                    return


if __name__ == "__main__":
    main()

"""Ist das Raupenpaar (mit Spinwolke) beweglich? Bandstruktur auf dem Torus.

Reproduziert README.md, Abschnitt 7.9 (Ausgabe als Markdown).
Teil 1: Pivot und Domino gegen Raupenwege zum selben Ziel, mit einer Fluktuation an verschiedenen Stellen (8×8-Torus).
Teil 2: Bandstruktur des Paares.
Start: python beweglichkeit.py            (L = 12, bis zu zwei Fluktuationen; braucht numpy und scipy, etwa 4 Minuten)
       python beweglichkeit.py --pruefen  (kleiner Torus: Bloch-Blöcke gegen die volle Prozessmatrix; braucht tkinter)

Zustände wie in fluktuationen2.py (Paar aus zwei benachbarten Löchern, Néel, bis zu KMAX Fluktuationen höchstens DMAX
vom Paar entfernt), aber auf einem L×L-Torus ohne Fenster. Prozesse: alle Folgen aus 1 bis 3 Elementarschritten
zwischen zwei Zuständen ohne Zwischenstopp. Verschiebungen um (dx, dy) mit dx + dy gerade erhalten das Néel-Muster;
mit ihnen zerfällt die Prozessmatrix in Blöcke M(K), einen für jeden Impuls K des Paares. Das oberste Band
λ(K) = größter Eigenwert von M(K) zeigt, ob das Paar sich über große Entfernungen kohärent bewegen kann: flach heißt
lokalisiert (unendlich schwer), gekrümmt heißt beweglich.
"""
import multiprocessing as mp
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

DMAX, KMAX = 2, 2
L = N = FULL = ROW = 0
NB, SUB, A_MASK, B_MASK, NEAR = [], [], 0, 0, []


def setup(size):
    global L, N, FULL, ROW, NB, SUB, A_MASK, B_MASK, NEAR
    L, N = size, size * size
    FULL, ROW = (1 << N) - 1, (1 << L) - 1
    NB = [[site(x + dx, y + dy) for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))] for x, y in map(coords, range(N))]
    SUB = [sum(coords(s)) % 2 for s in range(N)]
    A_MASK = sum(1 << s for s in range(N) if SUB[s] == 0)        # Néel: ↑ auf Untergitter A
    B_MASK = FULL ^ A_MASK
    NEAR = [[q for q in range(N) if dist(p, q) <= DMAX + 1] for p in range(N)]


def site(x, y):
    return (x % L) * L + y % L


def coords(s):
    return divmod(s, L)


def dist(a, b):
    (ax, ay), (bx, by) = coords(a), coords(b)
    dx, dy = abs(ax - bx), abs(ay - by)
    return min(dx, L - dx) + min(dy, L - dy)


def bits(m):
    out = []
    while m:
        low = m & -m
        out.append(low.bit_length() - 1)
        m ^= low
    return out


# --- Vorzeichen wie fluktuationen.hop / swap (NetKet-Ordnung: Orbital = Spin·N + Platz, ↓ = 0, ↑ = 1) ---

def csign(o, k):
    return -1 if (o & ((1 << k) - 1)).bit_count() & 1 else 1


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


# --- Zustände ---

def wrong_mask(up, dn):
    return (up & B_MASK) | (dn & A_MASK)


def decompositions(w):
    """Alle Zerlegungen der falschen Plätze in benachbarte Paare (eines auf A, eines auf B)."""
    if not w:
        return [[]]
    a, out = w[0], []
    for b in w[1:]:
        if dist(a, b) == 1:
            out += [[(a, b)] + r for r in decompositions([x for x in w if x not in (a, b)])]
    return out


def kind(up, dn):
    """None (kein Zustand) oder Zahl der Fluktuationen (0 bis KMAX)."""
    h = bits(FULL ^ (up | dn))
    if len(h) != 2 or dist(*h) != 1:
        return None
    w = bits(wrong_mask(up, dn))
    if len(w) > 2 * KMAX or len(w) % 2:
        return None
    if any(all(min(dist(x, y) for x in p for y in h) <= DMAX for p in dec) for dec in decompositions(w)):
        return len(w) // 2
    return None


def states_on(bond_list):
    all_bonds = [(s, NB[s][0]) for s in range(N)] + [(s, NB[s][2]) for s in range(N)]
    out = []
    for a, b in bond_list:
        up, dn = A_MASK & ~(1 << a) & ~(1 << b), B_MASK & ~(1 << a) & ~(1 << b)
        st = {(up, dn)}
        near = [(p, q) for p, q in all_bonds if not {p, q} & {a, b}
                and min(dist(x, y) for x in (p, q) for y in (a, b)) <= DMAX]
        frontier = [((up, dn), ())]
        for k in range(1, KMAX + 1):
            new = []
            for (u, d), used in frontier:
                for p, q in near:
                    if (p, q) <= (used[-1] if used else (-1, -1)) or {p, q} & {x for bd in used for x in bd}:
                        continue
                    _, u2, d2 = swap(u, d, p, q)
                    new.append(((u2, d2), used + ((p, q),)))
                    if kind(u2, d2) == k:
                        st.add((u2, d2))
            frontier = new
        out += sorted(st)
    return out


# --- Verschiebungen (Vorzeichen aus der Umordnung der Erzeuger) ---

def shift_x(m):
    last = m >> ((L - 1) * L)
    nl = last.bit_count()
    return ((m << L) & FULL) | last, nl * (m.bit_count() - nl)


def shift_y(m):
    out, inv = 0, 0
    for r in range(L):
        row = (m >> (r * L)) & ROW
        if row >> (L - 1) & 1:
            inv += row.bit_count() - 1
            row = ((row << 1) & ROW) | 1
        else:
            row <<= 1
        out |= row << (r * L)
    return out, inv


def translate(up, dn, sx, sy):
    inv = 0
    for _ in range(sx % L):
        up, i1 = shift_x(up)
        dn, i2 = shift_x(dn)
        inv += i1 + i2
    for _ in range(sy % L):
        up, i1 = shift_y(up)
        dn, i2 = shift_y(dn)
        inv += i1 + i2
    return up, dn, -1 if inv & 1 else 1


def rep_bonds():
    return [(site(0, 0), site(1, 0)), (site(1, 0), site(2, 0)), (site(0, 0), site(0, 1)), (site(1, 0), site(1, 1))]


def canon(up, dn):
    """Verschiebung S (dx + dy gerade), die den Zustand auf einen Vertreter bringt: (Index, Vorzeichen, S)."""
    a, b = bits(FULL ^ (up | dn))
    (ax, ay), (bx, by) = coords(a), coords(b)
    if ay == by:
        anchor = (ax, ay) if (bx - ax) % L == 1 else (bx, by)
    else:
        anchor = (ax, ay) if (by - ay) % L == 1 else (bx, by)
    c = sum(anchor) % 2
    sx, sy = (c - anchor[0]) % L, (-anchor[1]) % L
    u2, d2, eta = translate(up, dn, sx, sy)
    return IDX[(u2, d2)], eta, (sx, sy)


# --- Prozesse ---

STATES, IDX, MAPPER, BONDCACHE = [], {}, None, {}


def swap_bonds(hmask):
    if hmask not in BONDCACHE:
        sites = {q for h in bits(hmask) for q in NEAR[h]}
        BONDCACHE[hmask] = sorted({(min(p, q), max(p, q)) for p in sites for q in NB[p]})
    return BONDCACHE[hmask]


def steps(up, dn):
    """Lochsprünge und Vertauschungen nahe an den Löchern; die Auswahl hängt nur von den Lochplätzen ab."""
    hmask = FULL ^ (up | dn)
    for x in bits(hmask):
        for n in NB[x]:
            if not hmask >> n & 1:
                yield "h", *hop(up, dn, n, x)
    for p, q in swap_bonds(hmask):
        if (up >> p & 1 and dn >> q & 1) or (dn >> p & 1 and up >> q & 1):
            yield "s", *swap(up, dn, p, q)


def _process_one(i):
    up, dn = STATES[i]
    out, cache = [], {}
    stack = [(up, dn, 1, 0, 0, 0)]
    while stack:
        u, d, sg, m, s, depth = stack.pop()
        if depth == 3:
            continue
        for typ, sg1, u2, d2 in steps(u, d):
            m2, s2 = m + (typ == "h"), s + (typ == "s")
            key = (u2, d2)
            if key not in cache:
                cache[key] = MAPPER(u2, d2) if kind(u2, d2) is not None else None
            hit = cache[key]
            if hit is not None:
                j, eta, S = hit
                out.append((j, i, sg * sg1 * eta, S[0], S[1], m2, s2))
            elif wrong_mask(u2, d2).bit_count() <= 2 * KMAX + 2 * (2 - depth):
                stack.append((u2, d2, sg * sg1, m2, s2, depth + 1))
    return out


def enumerate_processes(states, mapper, workers=None):
    global STATES, IDX, MAPPER
    STATES, MAPPER = states, mapper
    IDX = {s: i for i, s in enumerate(states)}
    with mp.get_context("fork").Pool(workers or mp.cpu_count()) as pool:
        chunks = pool.map(_process_one, range(len(states)), chunksize=4)
    return np.array([p for c in chunks for p in c], dtype=np.int64).reshape(-1, 7)


# --- Bloch-Blöcke ---

def prepare(P, keep, ws, absolute=False):
    """Amplituden der Übergänge; absolute: Beträge der Gesamtamplituden (ohne Vorzeichen)."""
    J, I, SG, SX, SY, Mh, Sw = P[keep].T
    vals = SG * ws ** Sw.astype(float)
    if absolute:
        key = np.stack([J, I, SX, SY], axis=1)
        uniq, inv = np.unique(key, axis=0, return_inverse=True)
        vals = np.abs(np.bincount(inv.ravel(), weights=vals))
        J, I, SX, SY = uniq.T
    return J, I, vals, SX, SY


def bloch(prep, n, K):
    J, I, vals, SX, SY = prep
    ph = np.exp(1j * 2 * np.pi / L * (K[0] * SX + K[1] * SY))
    H = sp.csr_matrix((vals * ph, (J, I)), shape=(n, n))
    H.sum_duplicates()
    return H


def top(H, k=1):
    if H.shape[0] <= 300:
        w, v = np.linalg.eigh(H.toarray())
        return w[::-1][:k], v[:, ::-1][:, :k]
    w, v = sla.eigsh(H, k=k, which="LA", tol=1e-10)
    o = np.argsort(w)[::-1]
    return w[o], v[:, o]


def kpoints():
    """Impulse (Vielfache von 2π/L); K und K + (π, π) sind gleichwertig, daher nur die Hälfte."""
    seen, out = set(), []
    for i in range(L):
        for j in range(L):
            if (i, j) not in seen:
                seen |= {(i, j), ((i + L // 2) % L, (j + L // 2) % L)}
                out.append((i, j))
    return out


def analyse(P, kinds, keep_ms, kmax, ws):
    """keep_ms(m, s, kj, ki): Auswahl nach Zahl der Sprünge m, Vertauschungen s und Fluktuationen in Ziel und Start."""
    n = len(kinds)
    sub = kinds <= kmax
    keep = keep_ms(P[:, 5], P[:, 6], kinds[P[:, 0]], kinds[P[:, 1]]) & sub[P[:, 0]] & sub[P[:, 1]]
    ks = kpoints()
    prep, prep_abs = prepare(P, keep, ws), prepare(P, keep, ws, absolute=True)
    band, vecs = {}, {}
    for K in ks:
        w, v = top(bloch(prep, n, K), k=2)
        band[K], vecs[K] = w, v[:, 0]
    absband = {K: top(bloch(prep_abs, n, K))[0][0] for K in ks}
    lam = np.array([band[K][0] for K in ks])
    kbest = ks[int(np.argmax(lam))]
    lam_abs = absband[(0, 0)]
    abs_vals = np.array([absband[K] for K in ks])

    def curvature(f, K0):
        nbrs = [((K0[0] + dx) % L, (K0[1] + dy) % L) for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))]
        get = lambda K: f[K] if K in f else f[((K[0] + L // 2) % L, (K[1] + L // 2) % L)]
        return (get(K0) - np.mean([get(K) for K in nbrs])) / (2 * np.pi / L) ** 2
    lamd = {K: band[K][0] for K in ks}
    curv, curv_abs = curvature(lamd, kbest), curvature(absband, (0, 0))
    nflat = int(np.sum(lam > lam.max() - 1e-6 * abs(lam.max())))
    v = vecs[kbest]
    share = [float(np.sum(np.abs(v[kinds == k]) ** 2)) for k in range(kmax + 1)]
    return dict(ratio=lam.max() / lam_abs, kbest=kbest, width=lam.max() - lam.min(),
                width_abs=abs_vals.max() - abs_vals.min(), curv=curv, curv_abs=curv_abs, nflat=nflat, nk=len(ks),
                gap=band[kbest][0] - band[kbest][1], share=share)


def kname(K):
    f = lambda i: "0" if i == 0 else ("π" if 2 * i == L else f"{2 * i}π/{L}" if 2 * i < L else f"−{2 * (L - i)}π/{L}")
    return f"({f(K[0])}, {f(K[1])})"


def pruefen():
    """Kleiner Torus: Spektrum der vollen Prozessmatrix = Vereinigung der Spektren aller Bloch-Blöcke."""
    global KMAX
    sys.path.insert(0, ".")
    import random
    import fluktuationen as F
    rnd, ok, bad = random.Random(5), 0, 0
    setup(8)
    for _ in range(3000):
        up = dn = 0
        for s in range(64):
            r = rnd.random()
            up |= (r < 0.4) << s
            dn |= (0.4 <= r < 0.8) << s
        a = rnd.randrange(64)
        b = rnd.choice(F.NB[a])
        filled = lambda s: (up | dn) >> s & 1
        if filled(a) and not filled(b):
            ok += hop(up, dn, a, b) == F.hop(up, dn, a, b)
            bad += hop(up, dn, a, b) != F.hop(up, dn, a, b)
        elif filled(a) and filled(b) and (up >> a & 1) != (up >> b & 1):
            ok += swap(up, dn, a, b) == F.swap(up, dn, a, b)
            bad += swap(up, dn, a, b) != F.swap(up, dn, a, b)
    print(f"Vorzeichen gegen fluktuationen.py: {ok} gleich, {bad} verschieden.")
    for size, km in ((6, 0), (8, 0), (8, 1)):
        KMAX = km
        setup(size)
        allb = [(s, NB[s][0]) for s in range(N)] + [(s, NB[s][2]) for s in range(N)]
        full = states_on(allb)
        Pf = enumerate_processes(full, lambda u, d: (IDX[(u, d)], 1, (0, 0)))
        nf = len(full)
        Mf = sp.csr_matrix((Pf[:, 2] * 0.3 ** Pf[:, 6], (Pf[:, 0], Pf[:, 1])), shape=(nf, nf)).toarray()
        reps = states_on(rep_bonds())
        Pr = enumerate_processes(reps, canon)
        prep = prepare(Pr, np.ones(len(Pr), bool), 0.3)
        ev = []
        for i in range(L):
            for j in range(L):
                ev += list(np.linalg.eigvalsh(bloch(prep, len(reps), (i, j)).toarray()))
        ev = np.sort(ev)[::2]                         # K und K + (π, π) haben denselben Block
        ef = np.linalg.eigvalsh(Mf)
        print(f"L = {size}, höchstens {km} Fluktuation(en): {nf} Zustände, symmetrisch: "
              f"{np.abs(Mf - Mf.T).max() == 0}, größte Abweichung der Spektren {np.abs(np.sort(ev) - ef).max():.1e}")


def heal(up, dn, sign, target, region, depth=5):
    """Alle Vorzeichen, mit denen reine Vertauschungen in region den Zustand in target überführen
    ('Widerspruch', falls derselbe Zustand mit zwei Vorzeichen erreicht wird)."""
    bonds = [(p, q) for p in region for q in NB[p] if q in region and p < q]
    seen, frontier, out = {(up, dn): sign}, [(up, dn)], set()
    for level in range(depth + 1):
        new = []
        for u, d in frontier:
            sg = seen[(u, d)]
            if (u, d) == target:
                out.add(sg)
                continue
            if level == depth:
                continue
            for p, q in bonds:
                if (u >> p & 1 and d >> q & 1) or (d >> p & 1 and u >> q & 1):
                    s1, u2, d2 = swap(u, d, p, q)
                    if (u2, d2) in seen:
                        if seen[(u2, d2)] != sg * s1:
                            out.add("Widerspruch")
                        continue
                    seen[(u2, d2)] = sg * s1
                    new.append((u2, d2))
        frontier = new
    return out


def section_turn_rule():
    """Pivot und Domino gegen Raupenwege zum selben Ziel, mit einer Fluktuation an verschiedenen Stellen."""
    setup(8)
    S = site
    a, b, c, a2 = S(3, 3), S(4, 3), S(3, 4), S(4, 4)          # Paar a–b; Pivot und Domino bewegen die Bindung c–a2
    region = {s for s in range(N) if min(dist(s, x) for x in (a, b, c, a2)) <= 2}
    routes = {"Pivot": ([(c, a), (a2, c)], [(a2, b), (b, a)], (a2, b)),
              "Domino": ([(c, a), (a2, b)], [(a2, b), (b, a), (c, a2), (a2, b)], (c, a2))}

    def run(flip, hops, end):
        up, dn = A_MASK & ~(1 << a) & ~(1 << b), B_MASK & ~(1 << a) & ~(1 << b)
        sg = 1
        if flip:
            sg, up, dn = swap(up, dn, *flip)
        for src, dst in hops:
            s1, up, dn = hop(up, dn, src, dst)
            sg *= s1
        hm = (1 << end[0]) | (1 << end[1])
        return heal(up, dn, sg, (A_MASK & ~hm, B_MASK & ~hm), region)
    cases = [("keine", None), ("genau die bewegte Bindung c–a′ umgeklappt", (c, a2)),
             ("c mit seinem äußeren Nachbarn (bewegte Bindung gebrochen)", (c, S(2, 4))),
             ("a′ mit seinem oberen Nachbarn (bewegte Bindung gebrochen)", (a2, S(4, 5))),
             ("Bindung auf der anderen Seite des Paares (nicht berührt)", (S(3, 2), S(4, 2))),
             ("Bindung zwei Plätze entfernt (nicht berührt)", (S(3, 1), S(4, 1)))]
    print("Paar a = (3,3), b = (4,3). Pivot: a läuft über c = (3,4) nach a′ = (4,4); Domino: beide Löcher nach "
          "oben auf c, a′. Beide bewegen die Spins der Nachbarbindung c–a′ je einmal. Verglichen wird mit Raupenwegen "
          "zum selben Ziel (Raupe mit Wende bzw. zwei Raupenschritte), die jeden Spin zweimal bewegen. Danach wird mit "
          "reinen Vertauschungen bis Néel geheilt. Angegeben: relatives Vorzeichen (J-Wende) · (Raupenweg).\n")
    print("| Fluktuation vorher | Pivot gegen Raupe mit Wende | Domino gegen zwei Raupenschritte |")
    print("|---|---|---|")
    for name, flip in cases:
        cells = []
        for key in ("Pivot", "Domino"):
            turn, cat, end = routes[key]
            x, y = run(flip, turn, end), run(flip, cat, end)
            cells.append(" und ".join(sorted({f"{p * q:+d}".replace("-", "−") for p in x for q in y if isinstance(p, int)
                                              and isinstance(q, int)})) or "–" if "Widerspruch" not in x | y
                         else "Widerspruch")
        print(f"| {name} | {cells[0]} | {cells[1]} |")
    print()


def main():
    global KMAX
    t0 = time.time()
    print("## 7.9 Fällt die Paarkohärenz mit der Entfernung ab?\n")
    print("### Wann behält eine J-Wende ihr Vorzeichen?\n")
    section_turn_rule()
    print("### Bandstruktur auf dem Torus\n")
    setup(12)
    KMAX = 0
    reps0 = states_on(rep_bonds())
    P0 = enumerate_processes(reps0, canon)
    kinds0 = np.zeros(len(reps0), int)
    KMAX = 2
    reps = states_on(rep_bonds())
    kinds = np.array([kind(*s) for s in reps])
    print(f"Torus {L}×{L}, {len(kpoints())} Impulse. Sauberes Néel: {len(reps0)} Zustände je Impuls, {len(P0)} Prozesse. "
          f"Mit Fluktuationen: {len(reps)} Zustände je Impuls (" +
          ", ".join(f"{int(np.sum(kinds == k))} mit {k}" for k in range(KMAX + 1)) + " Fluktuationen), ", end="")
    P = enumerate_processes(reps, canon)
    print(f"{len(P)} Prozesse ({time.time() - t0:.0f} s).")
    for K in ((0, 0), (1, 2)):
        H = bloch(prepare(P, np.ones(len(P), bool), 0.3), len(reps), K)
        print(f"Hermitesch bei K = {kname(K)}: größte Abweichung {abs(H - H.getH()).max():.1e}")
    print()
    pd = lambda m, s: (m == 2) & (s == 1)
    variants = [("sauber: nur Raupe", P0, kinds0, lambda m, s, kj, ki: s == 0, 0),
                ("sauber: nur Pivot + Domino", P0, kinds0, lambda m, s, kj, ki: pd(m, s), 0),
                ("sauber: alle Prozesse (Raupe, Pivot, Domino)", P0, kinds0, lambda m, s, kj, ki: m >= 0, 0),
                ("Raupe + Fluktuationen, Sprünge ändern die Fluktuationszahl nicht, höchstens 2", P, kinds,
                 lambda m, s, kj, ki: ((s == 0) & (kj == ki)) | (m == 0), 2),
                ("nur Sprünge oder nur Vertauschungen (7.7/7.8), höchstens 1", P, kinds,
                 lambda m, s, kj, ki: (s == 0) | (m == 0), 1),
                ("nur Sprünge oder nur Vertauschungen (7.7/7.8), höchstens 2", P, kinds,
                 lambda m, s, kj, ki: (s == 0) | (m == 0), 2),
                ("alle Prozesse, höchstens 1", P, kinds, lambda m, s, kj, ki: m >= 0, 1),
                ("alle Prozesse, höchstens 2", P, kinds, lambda m, s, kj, ki: m >= 0, 2)]
    print("| Prozesse | Gewicht Vertauschung | Verhältnis | Maximum bei K | Impulse mit Maximum | "
          "Breite oberes Band (ohne Vorzeichen) | Krümmung am Maximum (ohne Vorzeichen) | Lücke zum 2. Band | "
          "Anteil 0 / 1 / 2 Fluktuationen |")
    print("|---|---|---|---|---|---|---|---|---|")
    for ws in (1.0, 0.3, 0.1):
        for name, PP, kk, keep_ms, kmax in variants:
            r = analyse(PP, kk, keep_ms, kmax, ws)
            print(f"| {name} | {ws} | {r['ratio']:.3f} | {kname(r['kbest'])} | {r['nflat']} von {r['nk']} | "
                  f"{r['width']:.3f} ({r['width_abs']:.3f}) | {r['curv']:.3f} ({r['curv_abs']:.3f}) | "
                  f"{r['gap']:.3f} | {' / '.join(f'{x:.2f}' for x in r['share'])} |")
    print(f"\nGesamtlaufzeit {time.time() - t0:.0f} s")


if __name__ == "__main__":
    if "--pruefen" in sys.argv:
        pruefen()
    else:
        main()

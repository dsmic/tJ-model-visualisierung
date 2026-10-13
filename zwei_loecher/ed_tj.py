"""Exakte Diagonalisierung des t-J-Modells (Standard: 4x4, periodisch) mit 0, 1 und 2 Loechern.

Verglichen werden
  tJ   : das t-J-Modell  H = -t sum P(c^dag_is c_js + h.c.)P + J sum (S_i.S_j - n_i n_j / 4)
  frei : seine vorzeichenfreie Version |H| (alle Nebendiagonal-Elemente -|H_fi|).
|H| enthaelt dieselben Prozesse mit denselben Betraegen, aber ohne Vorzeichen. Die Vorzeichen, die im
Repo vom Reparatur-Zeitpunkt abhaengen, sind genau das, was |H| weglaesst: das Gewicht eines
Prozesses in exp(-tau H) ist t^m (J/2)^s * (positive Faktoren) * Repo-Vorzeichen (siehe README.md).

Fermion-Ordnung wie NetKet / tj_gui.py: Orbital = spin_idx * N + site, down = 0, up = 1, site = x * L + y.

Start:  python ed_tj.py                      (4x4 periodisch, J/t = 0.2 0.4 0.6 1.0; etwa 20 Minuten)
        python ed_tj.py --L 3 --offen --J 0.4
        python ed_tj.py --pruefen            (Matrixelemente gegen tj_gui.swap_sign, braucht tkinter)
        python ed_tj.py --netket             (Matrixelemente und Energien gegen NetKet, kleines Gitter)
        python ed_tj.py --reihenfolge --J 0.4  (Entropie fuer verschiedene Platzreihenfolgen)
Braucht numpy (>= 2.0) und scipy.
"""
import argparse
import itertools
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = os.path.dirname(os.path.abspath(__file__))


def popcount(x):
    return np.bitwise_count(x).astype(np.int64)


def lattice_bonds(L, pbc):
    bonds = set()
    for x in range(L):
        for y in range(L):
            for dx, dy in ((0, 1), (1, 0)):
                nx, ny = x + dx, y + dy
                if pbc:
                    nx, ny = nx % L, ny % L
                elif nx >= L or ny >= L:
                    continue
                a, b = x * L + y, nx * L + ny
                if a != b:
                    bonds.add((min(a, b), max(a, b)))
    return sorted(bonds)


def wrap(L, pbc, c):
    """Komponente einer Verschiebung; periodisch in den Bereich (-L/2, L/2] gelegt."""
    if not pbc:
        return c
    h = (L - 1) // 2
    return (c + h) % L - h


def displacement(L, pbc, a, b):
    """Verschiebung b - a als (dx, dy); periodisch im Bereich (-L/2, L/2]."""
    (xa, ya), (xb, yb) = divmod(a, L), divmod(b, L)
    return wrap(L, pbc, xb - xa), wrap(L, pbc, yb - ya)


def canon(L, pbc, vecs):
    """Kanonische Form einer Menge von Verschiebungen unter den 8 Gittersymmetrien (D4)."""
    best = None
    for sx, sy, sw in itertools.product((1, -1), (1, -1), (False, True)):
        out = []
        for dx, dy in vecs:
            v = (sx * dx, sy * dy)
            v = (v[1], v[0]) if sw else v
            out.append((wrap(L, pbc, v[0]), wrap(L, pbc, v[1])))
        key = tuple(sorted(out))
        best = key if best is None or key < best else best
    return best


def manhattan(v):
    return abs(v[0]) + abs(v[1])


class Sector:
    """Basis und Operatoren fuer festes (n_up, n_dn). Zustand = Bitmuster der 2N Orbitale (int64)."""

    def __init__(self, L, pbc, n_up, n_dn):
        self.L, self.pbc, self.N = L, pbc, L * L
        self.n_up, self.n_dn = n_up, n_dn
        N = self.N
        self.bonds = lattice_bonds(L, pbc)
        states = []
        for up in itertools.combinations(range(N), n_up):
            ub = sum(1 << s for s in up)
            rest = [s for s in range(N) if not ub >> s & 1]
            for dn in itertools.combinations(rest, n_dn):
                states.append((ub << N) | sum(1 << s for s in dn))
        self.states = np.array(sorted(states), dtype=np.int64)
        self.dim = len(self.states)
        mask = (1 << N) - 1
        self.up = (self.states >> N) & mask
        self.dn = self.states & mask
        self.holes = ~(self.up | self.dn) & mask
        self.hop = self._hopping()
        self.flips = [self._flip(i, j) for i, j in self.bonds]
        self.zz = [0.25 * self._sz(i) * self._sz(j) for i, j in self.bonds]
        self.nn = [self._n(i) * self._n(j) for i, j in self.bonds]

    def _bit(self, arr, s):
        return (arr >> s) & 1

    def _sz(self, s):
        return (self._bit(self.up, s) - self._bit(self.dn, s)).astype(float)

    def _n(self, s):
        return (self._bit(self.up, s) + self._bit(self.dn, s)).astype(float)

    def _apply(self, ops):
        """ops: Liste (Orbital, erzeugen) in Anwendungsreihenfolge -> (Spalten, neue Zustaende, Vorzeichen)."""
        occ = self.states.copy()
        sign = np.ones(self.dim, dtype=np.int64)
        ok = np.ones(self.dim, dtype=bool)
        for k, create in ops:
            bit = np.int64(1) << k
            has = (occ & bit) != 0
            ok &= ~has if create else has
            sign *= 1 - 2 * (popcount(occ & (bit - 1)) % 2)
            occ ^= bit
        cols = np.nonzero(ok)[0]
        return cols, occ[ok], sign[ok]

    def _matrix(self, terms):
        rows, cols, vals = [], [], []
        for c, new, v in terms:
            r = np.searchsorted(self.states, new)
            assert np.all(self.states[r] == new)
            rows.append(r)
            cols.append(c)
            vals.append(v)
        rows, cols, vals = (np.concatenate(x) if x else np.zeros(0) for x in (rows, cols, vals))
        return sp.csr_matrix((vals.astype(float), (rows, cols)), shape=(self.dim, self.dim))

    def _hopping(self):
        """Matrix von sum_<ij>,s (c^dag_is c_js + h.c.) im Raum ohne Doppelbesetzung (Eintraege +-1)."""
        N, terms = self.N, []
        for i, j in self.bonds:
            for off in (0, N):
                for src, dst in ((i, j), (j, i)):
                    c, new, sg = self._apply([(off + src, False), (off + dst, True)])
                    other = (N - off) + dst  # Orbital der anderen Spinsorte am Zielplatz
                    free = ((new >> other) & 1) == 0
                    terms.append((c[free], new[free], sg[free]))
        return self._matrix(terms)

    def _flip(self, i, j):
        """1/2 (S+_i S-_j + S-_i S+_j) mit Fermion-Vorzeichen; S+_a S-_b = c^dag_a,up c_a,dn c^dag_b,dn c_b,up."""
        N, terms = self.N, []
        for a, b in ((i, j), (j, i)):
            c, new, sg = self._apply([(N + b, False), (b, True), (a, False), (N + a, True)])
            terms.append((c, new, 0.5 * sg))
        return self._matrix(terms)

    def hamiltonian(self, t, J, free=False):
        hop, flips = self.hop, self.flips
        if free:
            hop, flips = abs(hop), [abs(f) for f in flips]
            H = -t * hop - J * sum(flips)
        else:
            H = -t * hop + J * sum(flips)
        diag = J * sum(zz - 0.25 * nn for zz, nn in zip(self.zz, self.nn))
        return (H + sp.diags(diag)).tocsr()

    def hole_sites(self):
        """Lochplaetze je Zustand als Array (dim, Anzahl Loecher), aufsteigend."""
        nh = self.N - self.n_up - self.n_dn
        out = np.zeros((self.dim, nh), dtype=np.int64)
        h = self.holes.copy()
        for k in range(nh):
            low = h & -h
            out[:, k] = np.log2(low).astype(np.int64)
            h ^= low
        return out

    def site_major_sign(self, order=None):
        """Vorzeichen fuer den Wechsel von der NetKet-Ordnung zur platzweisen Ordnung (Orbitale nach Platz
        sortiert, Plaetze in der Reihenfolge order; Standard: site = x * L + y).
        Zaehlt die Paare besetzter Orbitale, deren Reihenfolge sich dabei umkehrt."""
        order = list(range(self.N)) if order is None else list(order)
        pos = np.empty(self.N, dtype=np.int64)
        pos[order] = np.arange(self.N)
        up = [self._bit(self.up, s) for s in range(self.N)]
        dn = [self._bit(self.dn, s) for s in range(self.N)]
        cnt = np.zeros(self.dim, dtype=np.int64)
        for i in range(self.N):
            for j in range(self.N):
                if i < j and pos[i] > pos[j]:      # gleiche Spinsorte, i vor j in NetKet
                    cnt += up[i] * up[j] + dn[i] * dn[j]
                if pos[i] < pos[j]:                # dn(j) vor up(i) in NetKet
                    cnt += up[i] * dn[j]
        return 1 - 2 * (cnt % 2)

    def squeezed(self, order=None):
        """(Index der Lochkonfiguration, Index der Spinfolge) je Zustand. Spinfolge = Spins der besetzten
        Plaetze in der Reihenfolge order (Ogata-Shiba-'squeezed space')."""
        order = range(self.N) if order is None else order
        pos = np.zeros(self.dim, dtype=np.int64)
        seq = np.zeros(self.dim, dtype=np.int64)
        for s in order:
            occ = self._bit(self.up | self.dn, s)
            seq |= (self._bit(self.up, s) * occ) << pos
            pos += occ
        _, xi = np.unique(self.holes, return_inverse=True)
        _, si = np.unique(seq, return_inverse=True)
        return xi, si


def lowest(H, k=16, tol=1e-7):
    """Niedrigste Eigenwerte; gibt (E0, Eigenvektoren aller entarteten Grundzustaende) zurueck.
    Lanczos kann einzelne entartete Partner verpassen, daher wird mit groesserem k nachgerechnet,
    bis die Entartung stabil ist."""
    n = H.shape[0]
    if n <= 3000:
        w, v = np.linalg.eigh(H.toarray())
        deg = int(np.sum(w < w[0] + tol * max(1.0, abs(w[0]))))
        return w[0], v[:, :deg]
    last = None
    while True:
        w, v = sla.eigsh(H, k=min(k, n - 2), which="SA", tol=1e-12, ncv=max(4 * k, 80))
        o = np.argsort(w)
        w, v = w[o], v[:, o]
        deg = int(np.sum(w < w[0] + tol * max(1.0, abs(w[0]))))
        if deg <= k - 4 and deg == last:
            return w[0], v[:, :deg]
        last, k = deg, k + 8


class Result:
    pass


def solve(L, pbc, n_holes, t, J, free, all_sz):
    """Grundzustand mit n_holes Loechern. tJ: Sektor mit minimalem Sz (SU(2) genuegt).
    frei: |H| ist nicht SU(2)-symmetrisch, daher alle Sektoren Sz >= 0 (Spiegelung up<->dn)."""
    N = L * L
    ne = N - n_holes
    ups = range((ne + 1) // 2, ne + 1) if (free and all_sz) else [(ne + 1) // 2]
    best = None
    for n_up in ups:
        sec = SECTORS.get((L, pbc, n_up, ne - n_up))
        if sec is None:
            sec = SECTORS[(L, pbc, n_up, ne - n_up)] = Sector(L, pbc, n_up, ne - n_up)
        E, vecs = lowest(sec.hamiltonian(t, J, free))
        if best is None or E < best.E - 1e-9:
            best = Result()
            best.E, best.vecs, best.sec, best.sz = E, vecs, sec, (2 * n_up - ne) / 2
    best.deg = best.vecs.shape[1]
    return best


SECTORS = {}


def expect(vecs, M):
    """Mittelwert ueber die entarteten Grundzustaende."""
    return float(np.mean([v @ (M @ v) for v in vecs.T]))


class Classes:
    """Ordnet Schluessel fortlaufende Nummern zu (fuer np.bincount)."""

    def __init__(self):
        self.keys, self.ids = [], {}

    def id(self, key):
        if key is None:
            return -1
        if key not in self.ids:
            self.ids[key] = len(self.keys)
            self.keys.append(key)
        return self.ids[key]


def pair_distribution(res):
    """Wahrscheinlichkeit je Lochpaar-Verschiebung (bis auf Symmetrie) und Zufallserwartung."""
    sec, L, pbc = res.sec, res.sec.L, res.sec.pbc
    hs = sec.hole_sites()
    cl = Classes()
    table = np.array([[cl.id(canon(L, pbc, [displacement(L, pbc, a, b)])[0]) if a != b else -1
                       for b in range(sec.N)] for a in range(sec.N)])
    prob = np.mean(res.vecs ** 2, axis=1)
    p = np.bincount(table[hs[:, 0], hs[:, 1]], weights=prob, minlength=len(cl.keys))
    ref = np.bincount(table[0, 1:], minlength=len(cl.keys))  # Zahl der Partnerplaetze je Klasse
    return {k: p[i] for i, k in enumerate(cl.keys)}, {k: ref[i] for i, k in enumerate(cl.keys)}


def bond_correlations(res, free=False):
    """<S_p.S_q> (xy- und zz-Teil) je Bindungsklasse: relativ zum Loch (1 Loch) bzw. nach Lochpaar-
    Verschiebung und Abstand der Bindung zum naechsten Loch (2 Loecher). Bindungen mit Loch zaehlen nicht.
    Rueckgabe {Klasse: (xy, zz, Gewicht)}; Gewicht = Wahrscheinlichkeit x Zahl der Bindungen.
    frei: xy-Teil mit dem Bindungsoperator von |H| (-|1/2 (S+S- + S-S+)|), also das Marshall-gedrehte S.S."""
    sec, L, pbc = res.sec, res.sec.L, res.sec.pbc
    hs = sec.hole_sites()
    nh = hs.shape[1]
    cl = Classes()
    nb = len(sec.bonds)
    if nh == 1:
        table = np.array([[cl.id(None if h in (p, q) else
                                 canon(L, pbc, [displacement(L, pbc, h, p), displacement(L, pbc, h, q)]))
                           for p, q in sec.bonds] for h in range(sec.N)])
        kid = table[hs[:, 0]]                                   # (dim, nb)
    else:
        table = np.full((sec.N, sec.N, nb), -1)
        for h1, h2 in itertools.combinations(range(sec.N), 2):
            dpair = canon(L, pbc, [displacement(L, pbc, h1, h2)])[0]
            for b, (p, q) in enumerate(sec.bonds):
                if {p, q} & {h1, h2}:
                    continue
                d = min(manhattan(displacement(L, pbc, hh, e)) for hh in (h1, h2) for e in (p, q))
                table[h1, h2, b] = cl.id((dpair, d))
        kid = table[hs[:, 0], hs[:, 1]]
    n = len(cl.keys)
    sums = np.zeros((3, n))
    for b in range(nb):
        k = kid[:, b]
        m = k >= 0
        for v in res.vecs.T:
            xy = v * ((-abs(sec.flips[b]) if free else sec.flips[b]) @ v)
            sums[0] += np.bincount(k[m], weights=xy[m], minlength=n) / res.deg
            sums[1] += np.bincount(k[m], weights=(v * v * sec.zz[b])[m], minlength=n) / res.deg
            sums[2] += np.bincount(k[m], weights=(v * v)[m], minlength=n) / res.deg
    return {key: (sums[0, i] / sums[2, i], sums[1, i] / sums[2, i], sums[2, i]) for i, key in enumerate(cl.keys)}


def charge_spin_entropy(res, free, order=None):
    """Verschraenkungsentropie zwischen Lochpositionen und Spinfolge (squeezed space), je Grundzustand.
    tJ: Wellenfunktion in platzweiser Fermion-Ordnung; frei: |psi| (Grundzustand des bosonischen |H|).
    order: Reihenfolge der Plaetze (Standard zeilenweise)."""
    sec = res.sec
    xi, si = sec.squeezed(order)
    g = sec.site_major_sign(order)
    out = []
    for v in res.vecs.T:
        amp = np.abs(v) if free else g * v
        M = np.zeros((xi.max() + 1, si.max() + 1))
        M[xi, si] = amp
        lam = np.linalg.svd(M, compute_uv=False) ** 2
        lam = lam[lam > 1e-14]
        out.append(float(-np.sum(lam * np.log(lam))))
    return out


def fmt_vec(v):
    return f"({v[0]},{v[1]})"


def run(args):
    L, pbc, t = args.L, not args.offen, 1.0
    lines = [f"# ED-Ergebnisse t-J {L}x{L} {'periodisch' if pbc else 'offen'}, t = 1\n"]
    lines.append("Varianten: **tJ** = t-J-Modell, **frei** = |H| (gleiche Prozesse, alle Vorzeichen entfernt). "
                 "Entartete Grundzustände werden gemittelt.\n")
    summary = ["| J/t | Variante | E(0) | E(1) | E(2) | Sz(1) | Sz(2) | Bindungsenergie E(2)+E(0)−2E(1) "
               "| ⟨Hüpfen⟩ 2 Löcher | Entartung 1/2 Loch |", "|---|---|---|---|---|---|---|---|---|---|"]
    details = []
    for J in args.J:
        for free in (False, True):
            name = "frei" if free else "tJ"
            t0 = time.time()
            r = [solve(L, pbc, nh, t, J, free, not args.nur_min_sz) for nh in (0, 1, 2)]
            eb = r[2].E + r[0].E - 2 * r[1].E
            kin = -t * expect(r[2].vecs, abs(r[2].sec.hop) if free else r[2].sec.hop)
            summary.append(f"| {J} | {name} | {r[0].E:.5f} | {r[1].E:.5f} | {r[2].E:.5f} | {r[1].sz} | {r[2].sz} "
                           f"| {eb:+.5f} | {kin:.4f} | {r[1].deg}/{r[2].deg} |")
            print(f"J={J} {name}: E={[round(float(x.E), 6) for x in r]} Eb={eb:+.5f} ({time.time() - t0:.0f}s)",
                  file=sys.stderr)
            d = [f"\n## J/t = {J}, {name}\n"]
            cls, ref = pair_distribution(r[2])
            tot = sum(ref.values())
            d.append("Lochabstand im 2-Loch-Grundzustand (Klasse = Verschiebung (dx,dy) bis auf Symmetrie):\n")
            d.append("| Verschiebung | Manhattan | P | P zufällig | P / zufällig |\n|---|---|---|---|---|")
            for key in sorted(ref, key=lambda k: (manhattan(k), k)):
                p = cls.get(key, 0.0)
                d.append(f"| {fmt_vec(key)} | {manhattan(key)} | {p:.4f} | {ref[key] / tot:.4f} "
                         f"| {p / (ref[key] / tot):.3f} |")
            if r[2].sec.dim <= 500000:
                S2 = charge_spin_entropy(r[2], free)
                S1 = charge_spin_entropy(r[1], free)
                d.append(f"\nVerschränkungsentropie Lochpositionen ↔ Spinfolge (squeezed space): "
                         f"2 Löcher {', '.join(f'{s:.4f}' for s in S2)}; "
                         f"1 Loch {', '.join(f'{s:.4f}' for s in S1)} "
                         f"(je Grundzustand; Maximum ln(Zahl der Lochkonfigurationen) = "
                         f"{math.log(math.comb(L * L, 2)):.3f} bzw. {math.log(L * L):.3f})\n")
            bc1 = bond_correlations(r[1], free)
            d.append("\n1 Loch: ⟨S_p·S_q⟩ je Bindung (p,q), Lage relativ zum Loch (Verschiebungen der Endpunkte)"
                     + (", xy-Teil mit dem Bindungsoperator von |H|" if free else "") + ":\n")
            d.append("| Bindung (Endpunkte rel. zum Loch) | xy-Teil | zz-Teil | S·S |\n|---|---|---|---|")
            for key in sorted(bc1, key=lambda k: (manhattan(k[0]) + manhattan(k[1]), k)):
                xy, zz, _ = bc1[key]
                d.append(f"| {fmt_vec(key[0])}–{fmt_vec(key[1])} | {xy:+.4f} | {zz:+.4f} | {xy + zz:+.4f} |")
            bc2 = bond_correlations(r[2], free)
            d.append("\n2 Löcher: ⟨S_p·S_q⟩ nach Lochpaar-Verschiebung und Abstand der Bindung zum nächsten Loch "
                     "(Gewicht = Wahrscheinlichkeit × Zahl der Bindungen):\n")
            d.append("| Lochpaar | Abstand Bindung–Loch | xy-Teil | zz-Teil | S·S | Gewicht |\n|---|---|---|---|---|---|")
            for key in sorted(bc2, key=lambda k: (manhattan(k[0]), k)):
                xy, zz, w = bc2[key]
                d.append(f"| {fmt_vec(key[0])} | {key[1]} | {xy:+.4f} | {zz:+.4f} | {xy + zz:+.4f} | {w:.3f} |")
            details += d
    lines += summary + details
    text = "\n".join(lines) + "\n"
    if args.out:
        with open(args.out, "w") as f:
            f.write(text)
    print(text)


def order_check(L, pbc, J):
    """Verschraenkungsentropie Loecher <-> Spins (2 Loecher) fuer verschiedene Platzreihenfolgen."""
    orders = {"zeilenweise": list(range(L * L)),
              "spaltenweise": [y * L + x for x in range(L) for y in range(L)],
              "Schlangenlinie": [x * L + (y if x % 2 == 0 else L - 1 - y) for x in range(L) for y in range(L)]}
    for free in (False, True):
        r = solve(L, pbc, 2, 1.0, J, free, False)
        for name, order in orders.items():
            S = charge_spin_entropy(r, free, order)
            print(f"J/t = {J}, {'frei' if free else 'tJ'}, {name}: {', '.join(f'{s:.4f}' for s in S)}")


def check_tj_gui(n=400, seed=1):
    """Vergleicht alle Nebendiagonal-Elemente (4x4 periodisch, 1 und 2 Loecher) mit tj_gui.swap_sign
    (4x4 eingebettet ins 8x8-Brett, restliche Plaetze leer; die Platzreihenfolge bleibt erhalten)."""
    sys.path.insert(0, os.path.join(HERE, ".."))
    import tj_gui as T
    rng = np.random.default_rng(seed)
    L, t, J, ok, bad = 4, 1.0, 0.4, 0, 0
    for n_up, n_dn in ((8, 7), (7, 7)):
        sec = Sector(L, True, n_up, n_dn)
        H = sec.hamiltonian(t, J).tocsc()
        for i in rng.choice(sec.dim, size=n, replace=False):
            conf = [T.HOLE] * T.N
            for s in range(sec.N):
                x, y = divmod(s, L)
                if sec.up[i] >> s & 1:
                    conf[T.site(x, y)] = T.UP
                elif sec.dn[i] >> s & 1:
                    conf[T.site(x, y)] = T.DN
            col = H[:, i]
            for f, val in zip(col.indices, col.data):
                if f == i:
                    continue
                diff = [s for s in range(sec.N) if ((sec.up[f] ^ sec.up[i]) | (sec.dn[f] ^ sec.dn[i])) >> s & 1]
                a, b = (T.site(*divmod(s, L)) for s in diff)
                sg = T.swap_sign(conf, a, b)
                hop = T.HOLE in (conf[a], conf[b])
                expected = -t * sg if hop else -(J / 2) * sg
                if abs(val - expected) < 1e-12:
                    ok += 1
                else:
                    bad += 1
    print(f"Matrixelemente gegen tj_gui.swap_sign: {ok} gleich, {bad} verschieden")
    return bad == 0


def check_netket(J=0.4):
    """Vergleicht die komplette Hamiltonmatrix (3x3 offen und periodisch, 1 und 2 Loecher) mit einem
    NetKet-Operator aus Fermion-Operatoren: projiziertes Huepfen c^dag_is (1 - n_i,-s) c_js (1 - n_j,-s) und
    J (S_i.S_j - n_i n_j / 4). Der Raum ohne Doppelbesetzung ist invariant; die NetKet-Matrix wird
    ueber get_conn auf ihm aufgebaut. Verglichen werden alle Spalten und die Grundzustandsenergie."""
    import netket as nk
    from netket.operator.fermion import create as cdag, destroy as c, number as nop
    t, good = 1.0, True
    for L, pbc in ((3, False), (3, True)):
        N = L * L
        for nh in (1, 2):
            ne = N - nh
            n_up, n_dn = (ne + 1) // 2, ne // 2
            sec = Sector(L, pbc, n_up, n_dn)
            hi = nk.hilbert.SpinOrbitalFermions(N, s=1 / 2, n_fermions_per_spin=(n_dn, n_up))
            H = 0.0
            for i, j in sec.bonds:
                for s in (-1, 1):
                    for a, b in ((i, j), (j, i)):
                        H += -t * (cdag(hi, a, s) @ (1 - nop(hi, a, -s)) @ c(hi, b, s) @ (1 - nop(hi, b, -s)))
                spm = cdag(hi, i, 1) @ c(hi, i, -1) @ cdag(hi, j, -1) @ c(hi, j, 1)
                smp = cdag(hi, j, 1) @ c(hi, j, -1) @ cdag(hi, i, -1) @ c(hi, i, 1)
                szi = 0.5 * (nop(hi, i, 1) - nop(hi, i, -1))
                szj = 0.5 * (nop(hi, j, 1) - nop(hi, j, -1))
                ni = nop(hi, i, 1) + nop(hi, i, -1)
                nj = nop(hi, j, 1) + nop(hi, j, -1)
                H += J * (szi @ szj + 0.5 * (spm + smp) - 0.25 * (ni @ nj))
            weights = np.int64(1) << np.arange(2 * N, dtype=np.int64)
            x = ((sec.states[:, None] >> np.arange(2 * N)) & 1).astype(np.int8)
            rows, cols, vals = [], [], []
            for k in range(sec.dim):
                xp, mels = H.get_conn(x[k])
                key = (np.asarray(xp, dtype=np.int64) * weights).sum(axis=1)
                key, inv = np.unique(key, return_inverse=True)       # gleiche Zielzustaende addieren
                mel = np.bincount(inv.ravel(), weights=np.real(mels), minlength=len(key))
                keep = np.abs(mel) > 1e-12
                key, mel = key[keep], mel[keep]
                r = np.searchsorted(sec.states, key)
                if np.any(r >= sec.dim) or np.any(sec.states[np.minimum(r, sec.dim - 1)] != key):
                    good = False  # Verbindung aus dem Raum ohne Doppelbesetzung heraus
                    continue
                rows.append(np.full(len(r), k))  # get_conn liefert Zeile k: H[k, x']
                cols.append(r)
                vals.append(mel)
            ref = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                shape=(sec.dim, sec.dim))
            mine = sec.hamiltonian(t, J)
            diff = abs(ref - mine).max()
            e_nk, e_me = lowest(ref)[0], lowest(mine)[0]
            print(f"NetKet {L}x{L} {'periodisch' if pbc else 'offen'}, {nh} {'Loch' if nh == 1 else 'Löcher'} (n_up={n_up}, "
                  f"n_dn={n_dn}, dim={sec.dim}): max |H_NetKet - H_hier| = {diff:.1e}, "
                  f"E0 NetKet {e_nk:.10f}, hier {e_me:.10f}")
            good &= diff < 1e-10
    return good


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--L", type=int, default=4)
    ap.add_argument("--offen", action="store_true", help="offene statt periodische Raender")
    ap.add_argument("--J", type=float, nargs="+", default=[0.2, 0.4, 0.6, 1.0])
    ap.add_argument("--nur-min-sz", action="store_true", help="auch fuer |H| nur den Sektor mit minimalem Sz")
    ap.add_argument("--out", help="Ergebnisse zusaetzlich als Markdown-Datei schreiben")
    ap.add_argument("--pruefen", action="store_true", help="Matrixelemente gegen tj_gui.swap_sign pruefen")
    ap.add_argument("--netket", action="store_true", help="Matrixelemente und Energien gegen NetKet pruefen")
    ap.add_argument("--reihenfolge", action="store_true",
                    help="Entropie Loecher <-> Spins fuer verschiedene Platzreihenfolgen (erstes J)")
    a = ap.parse_args()
    if a.reihenfolge:
        order_check(a.L, not a.offen, a.J[0])
        sys.exit(0)
    if a.pruefen or a.netket:
        good = True
        if a.pruefen:
            good &= check_tj_gui()
        if a.netket:
            good &= check_netket()
        sys.exit(0 if good else 1)
    run(a)

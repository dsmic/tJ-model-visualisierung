"""Zwei Loecher im t-J-Modell: Wege, Sprung-Reihenfolgen, Reparatur-Zeitpunkte, Vorzeichen.

Baut auf ../tj_gui.py auf (gleiche Vorzeichenkonvention, NetKet-Fermion-Ordnung, 8x8, offene Raender).
Ein Prozess ist eine Liste von Vertauschungen (a, b): Lochspruenge und Reparaturen
(Vertauschen benachbarter falscher up/dn-Paare).
"""
import itertools
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import tj_gui as T  # noqa: E402

L, N, HOLE = T.L, T.N, T.HOLE
DELTA = {"R": (0, 1), "L": (0, -1), "U": (-1, 0), "D": (1, 0)}
WHEN = {"ende": "Reparatur am Ende", "nach_a": "Reparatur nach A (bevor B weiterläuft)",
        "sofort": "Reparatur sofort"}


def sublattice(s):
    return (s // L + s % L) % 2


def parse_moves(text):
    return "".join(c for c in text.upper() if c in "RLUD")


def parse_offset(text):
    """'-1,1' -> (-1, 1): Startplatz von B relativ zu A (Zeile, Spalte)."""
    dx, dy = (int(v) for v in text.replace(";", ",").split(","))
    return dx, dy


def neel_with_holes(holes):
    """Volle Neel-Konfiguration, Loecher nur auf den angegebenen Plaetzen."""
    conf = [T.UP if sublattice(s) == 0 else T.DN for s in range(N)]
    for h in holes:
        conf[h] = HOLE
    return conf


def random_start_with(holes):
    """Zufaellige Neel-Konfiguration wie T.neel_random_start (28 up, 28 dn, 8 Loecher), in der die
    angegebenen Plaetze Loecher sind. Fehlende Loecher werden mit zufaelligen anderen Loechern desselben
    Untergitters vertauscht; der Spin passt dann auf beiden Plaetzen zum Neel-Muster."""
    conf = T.neel_random_start()
    for h in holes:
        if conf[h] != HOLE:
            cand = [s for s in range(N) if conf[s] == HOLE and s not in holes
                    and sublattice(s) == sublattice(h)]
            s = random.choice(cand)
            conf[h], conf[s] = conf[s], conf[h]
    return conf


def path_sites(h, moves):
    """Plaetze des Weges ab h (inklusive Start); None, wenn der Weg das Brett verlaesst."""
    x, y = divmod(h, L)
    out = [h]
    for m in moves:
        x, y = x + DELTA[m][0], y + DELTA[m][1]
        if not (0 <= x < L and 0 <= y < L):
            return None
        out.append(T.site(x, y))
    return out


def hops_of(sites):
    return [(sites[i], sites[i + 1]) for i in range(len(sites) - 1)]


def default_order(na, nb):
    return "A" * na + "B" * nb


def merge(ha, hb, order=None):
    """Mischt die Spruenge von A und B nach order (Folge aus 'A'/'B'; None = erst A, dann B)."""
    order = (order or default_order(len(ha), len(hb))).upper()
    if sorted(order) != sorted(default_order(len(ha), len(hb))):
        raise ValueError(f"Reihenfolge braucht {len(ha)}x A und {len(hb)}x B")
    ia, ib = iter(ha), iter(hb)
    return [next(ia) if o == "A" else next(ib) for o in order]


def all_orders(na, nb):
    for pos in itertools.combinations(range(na + nb), na):
        yield "".join("A" if k in pos else "B" for k in range(na + nb))


def legal(conf, hops):
    """Jeder Sprung geht von einem Loch auf einen Spin (kein Loch im Weg)."""
    c = list(conf)
    for a, b in hops:
        if c[a] != HOLE or c[b] == HOLE:
            return False
        c[a], c[b] = c[b], c[a]
    return True


def classify(pa, pb):
    """Geometrie der beiden Wege (Plaetze inkl. Start)."""
    out = []
    if set(pa) & set(pb):
        out.append("Kreuzung A/B")
    if len(pa) != len(set(pa)) or len(pb) != len(set(pb)):
        out.append("Schleife")
    return ", ".join(out) or "kein Platz doppelt"


def step_counts(start, steps):
    """(m, s, n-): Lochspruenge, Reparaturen, Spruenge eines falschen Spins."""
    c, m, nm = list(start), 0, 0
    for a, b in steps:
        if HOLE in (c[a], c[b]):
            m += 1
            src = b if c[a] == HOLE else a
            nm += src in T.neel_violations(c)
        c[a], c[b] = c[b], c[a]
    return m, len(steps) - m, nm


def repair_available(conf):
    """Repariert so viele benachbarte falsche up/dn-Paare wie moeglich (zufaellige Wahl);
    aendert conf, gibt die Vertauschungen zurueck. Einzelne falsche Spins koennen uebrig bleiben."""
    swaps = []
    while True:
        pairs = T.wrong_pairs(conf)
        if not pairs:
            return swaps
        a, b = random.choice(pairs)
        conf[a], conf[b] = conf[b], conf[a]
        swaps.append((a, b))


def strategy_steps(start, hops, order, when):
    """Prozess mit fester Reparatur-Strategie, None wenn am Ende nicht per Paartausch reparierbar.
    when: 'ende' (nur nach dem letzten Sprung), 'sofort' (vor jedem Sprung alles Reparierbare),
    'nach_a' (nach dem letzten Sprung von A alles Reparierbare, dann B weiter, Rest am Ende)."""
    c, steps = list(start), []
    last_a = max((k for k, o in enumerate(order) if o == "A"), default=-1)
    for k, (a, b) in enumerate(hops):
        if when == "sofort":
            steps += repair_available(c)
        c[a], c[b] = c[b], c[a]
        steps.append((a, b))
        if when == "nach_a" and k == last_a:
            steps += repair_available(c)
    rest = T.repair_neel(c)
    return None if rest is None else steps + rest


def search_orders(start, ha, hb, order=None, cap=20000):
    """Erschoepfende Suche nach Reparaturen (beliebig zwischen den Spruengen) fuer Endvorzeichen +1/-1.
    order: feste Reihenfolge, None (erst A, dann B) oder 'ALLE' (Vereinigung ueber alle legalen
    Reihenfolgen). Gibt ({Vorzeichen: Schrittliste}, Zahl legaler Reihenfolgen, abgebrochen)."""
    if order and order.upper() == "ALLE":
        orders = list(all_orders(len(ha), len(hb)))
    else:
        orders = [order or default_order(len(ha), len(hb))]
    found, n_legal, cut = {}, 0, False
    for o in orders:
        hops = merge(ha, hb, o)
        if not legal(start, hops):
            continue
        n_legal += 1
        f, _, c = T.search_signs(start, hops, cap=cap)
        cut |= c
        for sg, steps in f.items():
            found.setdefault(sg, steps)
        if len(found) == 2:
            break
    return found, n_legal, cut


def random_trial(ma, mb, offset, order=None, tries=20000):
    """Zufaelliger gueltiger Start fuer Weg A (ab A) und Weg B (ab A + offset).
    Gibt (start, a, b, Spruenge A, Spruenge B) oder None."""
    na, nb = len(ma), len(mb)
    for _ in range(tries):
        a = random.randrange(N)
        x, y = divmod(a, L)
        bx, by = x + offset[0], y + offset[1]
        if not (0 <= bx < L and 0 <= by < L) or (bx, by) == (x, y):
            continue
        b = T.site(bx, by)
        pa, pb = path_sites(a, ma), path_sites(b, mb)
        if pa is None or pb is None:
            continue
        start = random_start_with([a, b])
        ha, hb = hops_of(pa), hops_of(pb)
        if order and order.upper() == "ALLE":
            if not any(legal(start, merge(ha, hb, o)) for o in all_orders(na, nb)):
                continue
        elif not legal(start, merge(ha, hb, order)):
            continue
        return start, a, b, ha, hb
    return None


def j_term_check(n=3000, seed=3):
    """Prueft: T.swap_sign = -<f| S+_i S-_j |i> (Spin-Flip-Term in Hamilton-Ordnung).
    Dann ist das Repo-Vorzeichen eines Prozesses genau das Vorzeichen seines Gewichts in exp(-tau H):
    Lochsprung: -H = +t c^dagger c, Paartausch: -H = -(J/2)(S+S- + S-S+) = (J/2) * swap_sign."""
    rnd = random.Random(seed)
    ok = 0
    for _ in range(n):
        conf = T.neel_random_start()
        for _ in range(rnd.randint(0, 30)):
            a, b = rnd.choice(T.BONDS)
            conf[a], conf[b] = conf[b], conf[a]
        pairs = [(a, b) for a, b in T.BONDS if HOLE not in (conf[a], conf[b]) and conf[a] != conf[b]]
        a, b = rnd.choice(pairs)
        i, j = (a, b) if conf[a] == T.DN else (b, a)
        occ, sg = T.to_occ(conf), 1
        # S+_i S-_j = c^dagger_{i up} c_{i dn} c^dagger_{j dn} c_{j up}, von rechts angewendet
        for k, create in [(T.orbital(j, T.UP), False), (T.orbital(j, T.DN), True),
                          (T.orbital(i, T.DN), False), (T.orbital(i, T.UP), True)]:
            sg *= T.apply_c(occ, k, create)
        ok += sg == -T.swap_sign(conf, a, b)
    return ok, n

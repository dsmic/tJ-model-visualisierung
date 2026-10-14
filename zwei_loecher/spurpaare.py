"""Paare auf gemeinsamer Spur: Wie beweglich ist ein Paar, wenn nur Prozesse mit leerer ungerader Menge überleben?

Reproduziert README.md, Abschnitt 7.18 (Ausgabe als Markdown).
Start: python spurpaare.py   (braucht numpy und scipy wie beweglichkeit.py, unter einer Minute)

Sauberes Néel auf dem 12×12-Torus, zwei benachbarte Löcher. Erlaubt sind nur Lochsprünge; gezählt werden alle
Folgen von höchstens TIEFE Sprüngen von einem sauberen Paarzustand zum nächsten, ohne Zwischenstopp. Ohne
Vertauschungen endet ein Prozess genau dann wieder im sauberen Néel, wenn jedes Elektron gerade oft bewegt wurde
(O leer): Das zweite Loch hebt auf, was das erste verschoben hat, auch mit Abstand dazwischen. TIEFE 2 ist die reine
Raupe; größere Tiefe erlaubt, dass sich die Löcher entlang der Spur voneinander entfernen. Jeder Prozess bekommt das
Gewicht w^(Zahl der Sprünge). Ausgewertet wird das oberste Band wie in beweglichkeit.py.
"""
import multiprocessing as mp
from collections import Counter

import numpy as np

import beweglichkeit as B

TIEFE = 2


def _one(i):
    up, dn = B.STATES[i]
    out, stack = [], [(up, dn, 1, 0)]
    while stack:
        u, d, sg, m = stack.pop()
        if m == TIEFE:
            continue
        hm = B.FULL ^ (u | d)
        for x in B.bits(hm):
            for n in B.NB[x]:
                if hm >> n & 1:
                    continue
                s1, u2, d2 = B.hop(u, d, n, x)
                if B.kind(u2, d2) == 0:
                    j, eta, S = B.canon(u2, d2)
                    out.append((j, i, sg * s1 * eta, S[0], S[1], m + 1, 0))
                elif B.wrong_mask(u2, d2).bit_count() <= TIEFE - m - 1:
                    stack.append((u2, d2, sg * s1, m + 1))
    return out


def band(P, n, w):
    vals = P[:, 2] * w ** P[:, 5].astype(float)
    key = np.stack([P[:, 0], P[:, 1], P[:, 3], P[:, 4]], axis=1)
    uniq, inv = np.unique(key, axis=0, return_inverse=True)
    prep = (P[:, 0], P[:, 1], vals, P[:, 3], P[:, 4])
    prep_abs = (uniq[:, 0], uniq[:, 1], np.abs(np.bincount(inv.ravel(), weights=vals)), uniq[:, 2], uniq[:, 3])
    ks = B.kpoints()
    lam = np.array([np.linalg.eigvalsh(B.bloch(prep, n, K).toarray())[-1] for K in ks])
    lab = np.array([np.linalg.eigvalsh(B.bloch(prep_abs, n, K).toarray())[-1] for K in ks])

    def curv(arr, K0):
        f = dict(zip(ks, arr))
        get = lambda K: f[K] if K in f else f[((K[0] + B.L // 2) % B.L, (K[1] + B.L // 2) % B.L)]
        nb = [((K0[0] + a) % B.L, (K0[1] + b) % B.L) for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1))]
        return get(K0) - np.mean([get(K) for K in nb])
    kb = ks[int(np.argmax(lam))]
    nflat = int(np.sum(lam > lam.max() - 1e-9 * max(1.0, abs(lam.max()))))
    return lam.max() / lab.max(), nflat, len(ks), curv(lam, kb) / curv(lab, (0, 0))


def main():
    global TIEFE
    B.setup(12)
    B.KMAX = 0
    reps = B.states_on(B.rep_bonds())
    B.STATES, B.IDX, B.MAPPER = reps, {s: i for i, s in enumerate(reps)}, B.canon
    print("## 7.18 Paare auf gemeinsamer Spur\n")
    print("| höchstens Sprünge | Prozesse (nach Sprungzahl) | Gewicht w | Verhältnis | Impulse am Maximum | "
          "Beweglichkeit |")
    print("|---|---|---|---|---|---|")
    for TIEFE in (2, 4, 6, 8):
        with mp.get_context("fork").Pool(mp.cpu_count()) as pool:
            P = np.array([p for c in pool.map(_one, range(len(reps))) for p in c], dtype=np.int64).reshape(-1, 7)
        counts = ", ".join(f"{k}: {v}" for k, v in sorted(Counter(P[:, 5].tolist()).items()))
        for w in (0.3, 0.5, 0.7):
            ratio, nflat, nk, mob = band(P, len(reps), w)
            print(f"| {TIEFE} | {len(P)} ({counts}) | {w} | {ratio:.3f} | {nflat} von {nk} | {mob:.3f} |")


if __name__ == "__main__":
    main()

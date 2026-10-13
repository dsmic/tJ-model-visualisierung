"""Getrennte Partner und viele Löcher: Löcher als zwei Sorten Fermionen.

Reproduziert README.md, Abschnitte 7.12 und 7.13 (Ausgabe als Markdown).
Start: python partner.py   (braucht numpy und scipy wie beweglichkeit.py, wenige Sekunden)

Néel auf einem 10×10-Torus mit zwei Löchern: eine ↑-Lücke auf Untergitter A und eine ↓-Lücke auf B, am Anfang
benachbart oder 4 bis 5 Plätze entfernt. Jede Lücke läuft eine zufällige geschlossene Schleife aus J-Schritten
(zwei Sprünge diagonal oder gerade, danach sofort die Reparatur). Die Schritte beider Lücken werden zufällig
verschränkt; die Wege dürfen sich kreuzen. Verglichen wird das Vorzeichen des Gesamtprozesses mit dem Produkt der
Vorzeichen, die jede Schleife allein hat.
"""
import random

import beweglichkeit as B

MOVES = [(1, 1), (1, -1), (-1, 1), (-1, -1), (2, 0), (-2, 0), (0, 2), (0, -2)]


def neel(holes):
    hm = sum(1 << h for h in holes)
    return B.A_MASK & ~hm, B.B_MASK & ~hm


def jstep(up, dn, h, d, rnd, caterpillar=False):
    """Loch h macht einen J-Schritt um d (zwei Sprünge, danach Reparatur). Ist der Zwischenplatz ein Loch der anderen
    Sorte und caterpillar gesetzt, wird es ein Raupenschritt (h springt über dieses Loch, ohne J).
    None, falls nicht möglich."""
    x, y = B.coords(h)
    p2 = B.site(x + d[0], y + d[1])
    mids = [B.site(x + d[0], y), B.site(x, y + d[1])] if d[0] and d[1] else [B.site(x + d[0] // 2, y + d[1] // 2)]
    p1 = rnd.choice(mids)
    if caterpillar and not (up | dn) >> p1 & 1 and (up | dn) >> p2 & 1:
        s1, up, dn = B.hop(up, dn, p2, p1)              # Kopf p1 springt nach p2
        s2, up, dn = B.hop(up, dn, p1, h)               # Schwanz h folgt nach p1
        return None if B.wrong_mask(up, dn) else (s1 * s2, up, dn, p2, "R")
    if not ((up | dn) >> p1 & 1 and (up | dn) >> p2 & 1):
        return None
    s1, up, dn = B.hop(up, dn, p1, h)
    s2, up, dn = B.hop(up, dn, p2, p1)
    sg = s1 * s2
    w = B.wrong_mask(up, dn)
    if w:
        if not (w >> h & 1 and w >> p1 & 1):
            return None
        s3, up, dn = B.swap(up, dn, h, p1)
        sg *= s3
    if B.wrong_mask(up, dn):
        return None
    return sg, up, dn, p2, "J"


def loop(rnd, n):
    """Zufällige geschlossene Schleife aus n J-Schritten."""
    while True:
        out = [rnd.choice(MOVES) for _ in range(n)]
        if sum(a for a, _ in out) == 0 and sum(b for _, b in out) == 0:
            return out


def run(holes, plans, order, rnd):
    up, dn = neel(holes)
    pos, sg, idx = list(holes), 1, [0] * len(holes)
    for who in order:
        r = jstep(up, dn, pos[who], plans[who][idx[who]], rnd)
        if r is None:
            return None
        s1, up, dn, pos[who], _ = r
        sg *= s1
        idx[who] += 1
    assert (up, dn) == neel(holes) and pos == list(holes)
    return sg


def perm_sign(p):
    p, sg = list(p), 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            sg = -sg
    return sg


def many_trial(rnd, n_each=3, k=12):
    """Zufälliger geschlossener Prozess mit n_each ↑-Lücken und n_each ↓-Lücken. (ok, Platztausch?, Schrittzahlen,
    Treffer falscher Regeln) oder None."""
    A = [s for s in range(B.N) if B.SUB[s] == 0]
    Bs = [s for s in range(B.N) if B.SUB[s] == 1]
    ups, dns = [], []
    while len(ups) < n_each:                                   # teils Nachbarpaare, teils einzeln
        a = rnd.choice(A)
        b = rnd.choice(B.NB[a]) if rnd.random() < 0.6 else rnd.choice(Bs)
        if a not in ups and b not in dns:
            ups.append(a)
            dns.append(b)
    start = ups + dns
    targets = (start[:n_each], start[n_each:])
    pos, (up, dn) = list(start), neel(start)
    sg, kinds = 1, {"J": 0, "R": 0}

    def move(i, d):
        nonlocal up, dn, sg
        r = jstep(up, dn, pos[i], d, rnd, caterpillar=True)
        if r is None:
            return False
        s1, up, dn, new, kind = r
        kinds[kind] += 1
        sg *= s1
        pos[i] = new
        return True
    for _ in range(k):
        move(rnd.randrange(len(pos)), rnd.choice(MOVES))
    for _ in range(400):                                       # zurück, Zuordnung innerhalb der Sorte beliebig
        todo = [i for i in range(len(pos)) if pos[i] not in targets[i >= n_each]]
        if not todo:
            break
        i = rnd.choice(todo)
        tg = min((t for t in targets[i >= n_each] if t not in pos), key=lambda t: B.dist(pos[i], t))
        x, y = B.coords(pos[i])
        for d in sorted(MOVES, key=lambda d: (B.dist(B.site(x + d[0], y + d[1]), tg), rnd.random()))[:3]:
            if (B.dist(B.site(x + d[0], y + d[1]), tg) < B.dist(pos[i], tg) or rnd.random() < 0.2) and move(i, d):
                break
    else:
        return None
    if (up, dn) != neel(start):
        return None
    ps = perm_sign([targets[0].index(p) for p in pos[:n_each]]) * perm_sign([targets[1].index(p) for p in pos[n_each:]])
    steps = kinds["J"] + kinds["R"]
    wrong = {"ohne (−1)^Schritte": ps == sg, "ohne Permutation": (-1) ** steps == sg,
             "nur J-Schritte gezählt": (-1) ** kinds["J"] * ps == sg}
    return sg == (-1) ** steps * ps, ps == -1, kinds, wrong


def main(trials=3000, seed=4):
    B.setup(10)
    rnd = random.Random(seed)
    a = B.site(4, 4)                                                    # ↑-Lücke (Untergitter A)
    partners = [B.site(5, 4), B.site(4, 5), B.site(1, 6), B.site(7, 3)]  # ↓-Lücke (Untergitter B)
    prod = same = 0
    singles = {1: 0, -1: 0}
    n_all = n_same = 0
    for _ in range(trials):
        b = rnd.choice(partners)
        la = loop(rnd, rnd.randint(3, 6))
        lb = la if rnd.random() < 0.5 else loop(rnd, rnd.randint(3, 6))
        order = [0] * len(la) + [1] * len(lb)
        rnd.shuffle(order)
        both = run([a, b], [la, lb], order, rnd)
        sa = run([a, b], [la, []], [0] * len(la), rnd)
        sb = run([a, b], [[], lb], [1] * len(lb), rnd)
        if None in (both, sa, sb):
            continue
        n_all += 1
        prod += both == sa * sb
        singles[sa] += 1
        if la is lb:
            n_same += 1
            same += both == 1
    print("## 7.12 Getrennte Partner\n")
    print(f"Vorzeichen des Gesamtprozesses = Produkt der Einzelvorzeichen: {prod} von {n_all}.")
    print(f"Deckungsgleiche Schleifen (↓-Lücke läuft dieselbe Form wie die ↑-Lücke): Gesamtvorzeichen +1 in "
          f"{same} von {n_same}.")
    print(f"Einzelschleifen der ↑-Lücke: {singles[1]}-mal +1, {singles[-1]}-mal −1 (frustriert).\n")
    print("## 7.13 Viele Löcher\n")
    rnd = random.Random(7)
    ok = n = swapped = 0
    steps = {"J": 0, "R": 0}
    wrong = {}
    for _ in range(1500):
        r = many_trial(rnd)
        if r is None:
            continue
        good, sw, kinds, wr = r
        n += 1
        ok += good
        swapped += sw
        for key in steps:
            steps[key] += kinds[key]
        for key, val in wr.items():
            wrong[key] = wrong.get(key, 0) + val
    print(f"Vorzeichen = (−1)^(Schritte) · sgn(π der ↑-Lücken) · sgn(π der ↓-Lücken): {ok} von {n} "
          f"({swapped} mit Platztausch innerhalb einer Sorte; {steps['J']} J-Schritte, {steps['R']} Raupenschritte).")
    print("Falsche Regeln zum Vergleich, Treffer: " + ", ".join(f"{k} {v}" for k, v in wrong.items()) + ".")


if __name__ == "__main__":
    main()

"""Getrennte Partner: Das Vorzeichen zweier unabhängig laufender Löcher ist ein Produkt.

Reproduziert README.md, Abschnitt 7.12 (Ausgabe als Markdown).
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


def jstep(up, dn, h, d, rnd):
    """Loch h macht einen J-Schritt um d (zwei Sprünge, danach Reparatur). None, falls nicht möglich."""
    x, y = B.coords(h)
    p2 = B.site(x + d[0], y + d[1])
    mids = [B.site(x + d[0], y), B.site(x, y + d[1])] if d[0] and d[1] else [B.site(x + d[0] // 2, y + d[1] // 2)]
    p1 = rnd.choice(mids)
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
    return sg, up, dn, p2


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
        s1, up, dn, pos[who] = r
        sg *= s1
        idx[who] += 1
    assert (up, dn) == neel(holes) and pos == list(holes)
    return sg


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
    print(f"Einzelschleifen der ↑-Lücke: {singles[1]}-mal +1, {singles[-1]}-mal −1 (frustriert).")


if __name__ == "__main__":
    main()

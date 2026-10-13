"""Reproduziert die Zwei-Loch-Ergebnisse aus README.md (Ausgabe als Markdown-Tabellen).

Start: python experimente.py   (braucht nur tkinter wie tj_gui.py; Laufzeit einige Minuten)
"""
import random
from collections import Counter

import zwei_loch as Z
from zwei_loch import T

S = T.site


def fmt_sites(*ss):
    return ", ".join(f"({s // Z.L},{s % Z.L})" for s in ss)


def kinds(found):
    return "beide" if len(found) == 2 else ("keines" if not found else f"nur {next(iter(found)):+d}")


def section_jcheck():
    ok, n = Z.j_term_check()
    print("## J-Term-Check\n")
    print(f"swap_sign = −⟨f|S⁺S⁻|i⟩: {ok}/{n} Übereinstimmungen\n")


def section_crossing():
    print("## Kreuzungsbeispiel (volle Néel-Konfiguration, nur die Löcher A und B)\n")
    a, b = S(4, 2), S(3, 3)
    start = Z.neel_with_holes([a, b])
    ha, hb = Z.hops_of(Z.path_sites(a, "RR")), Z.hops_of(Z.path_sites(b, "DD"))
    order = Z.default_order(2, 2)
    hops = Z.merge(ha, hb, order)
    print(f"A ab {fmt_sites(a)}: RR, B ab {fmt_sites(b)}: DD, Kreuzung bei {fmt_sites(S(4, 3))}\n")
    print("| Ablauf | m | s | n₋ | Vorzeichen |\n|---|---|---|---|---|")
    for when, label in [("nach_a", "Spur von A repariert, bevor B kreuzt"),
                        ("ende", "Spur von A repariert, nachdem B gekreuzt hat")]:
        steps = Z.strategy_steps(start, hops, order, when)
        m, s, nm = Z.step_counts(start, steps)
        print(f"| {label} | {m} | {s} | {nm} | {T.apply_steps(start, steps)[1]:+d} |")
    print()


CASES = [
    ("getrennt", S(4, 2), "RR", S(1, 5), "DD"),
    ("getrennt", S(5, 1), "RRRR", S(0, 6), "DDDD"),
    ("Kreuzung", S(4, 2), "RR", S(3, 3), "DD"),
    ("Kreuzung", S(4, 1), "RRRR", S(2, 3), "DDDD"),
    ("B folgt A, gerade", S(4, 3), "RR", S(4, 2), "RR"),
    ("B folgt A, gerade", S(4, 3), "RRR", S(4, 2), "RRR"),
    ("B folgt A, Ecke", S(4, 3), "RU", S(4, 2), "RR"),
    ("B folgt A, Ecke", S(4, 3), "RRUU", S(4, 2), "RRRU"),
    ("A umrundet einen Spin", S(3, 3), "RRDDLLUU", S(0, 7), ""),
    ("A umrundet Loch B", S(3, 3), "RRDDLLUU", S(4, 4), ""),
    ("Austausch über Diagonale", S(3, 3), "RD", S(4, 4), "LU"),
    ("kein Austausch (hin und zurück)", S(3, 3), "RL", S(4, 4), "LR"),
]


def section_controlled():
    print("## Kontrollierte Fälle (alle Sprung-Reihenfolgen, Reparaturen beliebig dazwischen)\n")
    print("| Fall | A | B | Geometrie | Ergebnis je Reihenfolge | erreichbar insgesamt |")
    print("|---|---|---|---|---|---|")
    for name, a, ma, b, mb in CASES:
        start = Z.neel_with_holes([a, b])
        pa, pb = Z.path_sites(a, ma), Z.path_sites(b, mb)
        ha, hb = Z.hops_of(pa), Z.hops_of(pb)
        per, union = Counter(), set()
        for o in Z.all_orders(len(ha), len(hb)):
            hops = Z.merge(ha, hb, o)
            if not Z.legal(start, hops):
                continue
            found, _, _ = T.search_signs(start, hops, cap=200000)
            per[kinds(found)] += 1
            union |= set(found)
        res = ", ".join(f"{k}: {v}" for k, v in sorted(per.items()))
        print(f"| {name} | {fmt_sites(a)} {ma} | {fmt_sites(b)} {mb or '–'} | {Z.classify(pa, pb)} | {res} "
              f"| {' und '.join(f'{s:+d}' for s in sorted(union)) or 'keines'} |")
    print()


def section_encircle():
    print("## Umrundung mit fester Strategie (Reparatur sofort)\n")
    print("| Fall | Vorzeichen | m | s | n₋ |\n|---|---|---|---|---|")
    for name, holes in [("A umrundet einen Spin", [S(3, 3), S(0, 7)]),
                        ("A umrundet Loch B", [S(3, 3), S(4, 4)])]:
        start = Z.neel_with_holes(holes)
        hops = Z.hops_of(Z.path_sites(holes[0], "RRDDLLUU"))
        steps = Z.strategy_steps(start, hops, "A" * len(hops), "sofort")
        m, s, nm = Z.step_counts(start, steps)
        print(f"| {name} | {T.apply_steps(start, steps)[1]:+d} | {m} | {s} | {nm} |")
    print()


def section_random(n_cases=600, seed=7):
    print(f"## Zufallstest ({n_cases} Fälle, Seed {seed})\n")
    print("Zufällige Wege (Länge A 1–5, B 0–5, Gesamtlänge gerade, kein direktes Zurückgehen), Starts im "
          "Zentrum einer vollen Néel-Konfiguration mit nur A und B, zufällige Sprung-Reihenfolge, "
          "Reparaturen beliebig dazwischen.\n")
    rnd = random.Random(seed)
    stat, tested = {}, 0
    while tested < n_cases:
        a = S(rnd.randint(2, 5), rnd.randint(2, 5))
        b = S(rnd.randint(2, 5), rnd.randint(2, 5))
        la, lb = rnd.randint(1, 5), rnd.randint(0, 5)
        if a == b or (la + lb) % 2:
            continue
        ma = "".join(rnd.choice("RLUD") for _ in range(la))
        mb = "".join(rnd.choice("RLUD") for _ in range(lb))
        if any(x + y in ("RL", "LR", "UD", "DU") for m in (ma, mb) for x, y in zip(m, m[1:])):
            continue
        pa, pb = Z.path_sites(a, ma), Z.path_sites(b, mb)
        if pa is None or pb is None:
            continue
        order = list("A" * la + "B" * lb)
        rnd.shuffle(order)
        start = Z.neel_with_holes([a, b])
        hops = Z.merge(Z.hops_of(pa), Z.hops_of(pb), "".join(order))
        if not Z.legal(start, hops):
            continue
        found, _, cut = T.search_signs(start, hops, cap=100000)
        if cut:
            continue
        tested += 1
        k = "beide" if len(found) == 2 else ("keines" if not found else "eines")
        stat.setdefault(Z.classify(pa, pb), Counter())[k] += 1
    print("| Geometrie | beide Vorzeichen | eines | nicht reparierbar |\n|---|---|---|---|")
    for cat, c in sorted(stat.items()):
        print(f"| {cat} | {c['beide']} | {c['eines']} | {c['keines']} |")
    print()


def section_strategies(n=100, seed=11):
    print(f"## Strategievergleich auf zufälligen Starts ({n} Starts je Zeile, Seed {seed})\n")
    print("Zufällige Néel-Starts (28 up, 28 dn, 8 Löcher) wie in tj_gui.py. In jedem Start werden "
          "dieselben Sprünge mit drei Reparatur-Strategien ausgeführt; gezählt wird, wie oft sich das "
          "Endvorzeichen zwischen zwei Strategien im selben Start unterscheidet.\n")
    random.seed(seed)
    rows = [("Kreuzung", "RR", "DD", (-1, 1)), ("Kreuzung", "RRRR", "DD", (-1, 2)),
            ("Kreuzung, A mit Ecke", "RRUU", "LL", (-1, 3)),
            ("B folgt A", "RR", "RR", (0, -1)), ("B folgt A, Ecke", "RRUU", "RRRU", (0, -1)),
            ("getrennt", "RR", "DD", (-3, 3))]
    print("| Fall | A | B | B-Start rel. zu A | Ende ≠ nach A | Ende ≠ sofort | n₋ (Ende / nach A / sofort) "
          "| verworfen |")
    print("|---|---|---|---|---|---|---|---|")
    for name, ma, mb, off in rows:
        diff_na = diff_so = done = skipped = 0
        nms = {w: Counter() for w in ("ende", "nach_a", "sofort")}
        while done < n and skipped < 20 * n:
            t = Z.random_trial(ma, mb, off)
            if t is None:
                break
            start, a, b, ha, hb = t
            order = Z.default_order(len(ha), len(hb))
            hops = Z.merge(ha, hb, order)
            res = {w: Z.strategy_steps(start, hops, order, w) for w in nms}
            if any(r is None for r in res.values()):
                skipped += 1
                continue
            done += 1
            sg = {w: T.apply_steps(start, r)[1] for w, r in res.items()}
            for w, r in res.items():
                nms[w][Z.step_counts(start, r)[2]] += 1
            diff_na += sg["ende"] != sg["nach_a"]
            diff_so += sg["ende"] != sg["sofort"]
        nm_txt = " / ".join(",".join(str(k) for k in sorted(nms[w])) for w in ("ende", "nach_a", "sofort"))
        print(f"| {name} | {ma} | {mb} | {off} | {diff_na}/{done} | {diff_so}/{done} | {nm_txt} | {skipped} |")
    print()


if __name__ == "__main__":
    section_jcheck()
    section_crossing()
    section_controlled()
    section_encircle()
    section_strategies()
    section_random()

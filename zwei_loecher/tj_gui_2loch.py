"""Erweiterung von ../tj_gui.py fuer zwei Loecher A und B (alles aus tj_gui.py bleibt verfuegbar).

Start: python tj_gui_2loch.py   (nur tkinter noetig)

Zusaetzliche Leisten unten:
- Weg A, Weg B (Folgen aus R/L/U/D), Startplatz von B relativ zu A (Zeile,Spalte) und Reihenfolge der
  Spruenge (Folge aus A/B; leer = erst A, dann B; 'alle' = alle Reihenfolgen, nur fuer die Suche).
- "2 Löcher: Strategien vergleichen": 100 zufaellige Starts; in jedem Start dieselben Spruenge mit
  Reparatur am Ende / nach A (bevor B weiterlaeuft) / sofort. Ausgabe: Vorzeichen, n- und wie oft sich
  die Strategien im selben Start unterscheiden.
- "Variante laden": laedt nacheinander Start und Strategie (Starts mit verschiedenen Vorzeichen zuerst);
  mit "Vor" wird der Prozess abgespielt.
- "2 Löcher: Suche +1/-1": erschoepfende Suche ueber Reparaturen zwischen den Spruengen;
  danach laden "Beispiel +1" / "Beispiel -1" (zweite Leiste) die gefundenen Prozesse.
- "Beispiel Kreuzung": A RR, B DD, B startet bei (-1,1): Reparatur vor oder nach dem Kreuzen.
Die Loecher A und B sind beschriftet; die Beschriftung wandert mit dem Loch.
"""
import tkinter as tk

import zwei_loch as Z
from zwei_loch import T


class App2(T.App):
    def __init__(self, root):
        self.labels = {}
        super().__init__(root)
        root.title("t-J-Modell 8x8 – zwei Löcher")
        main = self.canvas.master
        bar4 = tk.Frame(main)
        bar4.grid(row=4, column=0, sticky="ew")
        self.path_a = tk.StringVar(value="RR")
        self.path_b = tk.StringVar(value="DD")
        self.offset = tk.StringVar(value="-1,1")
        self.order = tk.StringVar(value="")
        for text, var, width in (("Weg A:", self.path_a, 10), ("Weg B:", self.path_b, 10),
                                 ("B-Start relativ zu A (Zeile,Spalte):", self.offset, 6),
                                 ("Reihenfolge (A/B, leer, alle):", self.order, 10)):
            tk.Label(bar4, text=text).pack(side="left")
            tk.Entry(bar4, textvariable=var, width=width).pack(side="left", padx=(0, 6))
        bar5 = tk.Frame(main)
        bar5.grid(row=5, column=0, sticky="ew")
        tk.Button(bar5, text="2 Löcher: Strategien vergleichen", command=self.run_two_test).pack(side="left")
        tk.Button(bar5, text="Variante laden", command=self.load_variant).pack(side="left")
        tk.Button(bar5, text="2 Löcher: Suche +1/-1", command=self.run_two_search).pack(side="left")
        tk.Button(bar5, text="Beispiel Kreuzung", command=self.crossing_example).pack(side="left")
        self.two_results = []
        self.variant_index = 0

    # --- Beschriftung der Loecher (wandert mit, ist Teil von Zurueck/Vor) ---

    def new(self):
        self.labels = {}
        super().new()

    def snapshot(self):
        return super().snapshot() + (dict(self.labels),)

    def restore(self, state):
        self.labels = dict(state[3]) if len(state) > 3 else {}
        super().restore(state)

    def load_saved(self, i):
        self.labels = {}
        super().load_saved(i)

    def click_swap(self, a, b):
        self.push_history()
        if T.HOLE in (self.conf[a], self.conf[b]):
            self.trail.append((a, b))
        self.sign *= T.swap_sign(self.conf, a, b)
        self.conf[a], self.conf[b] = self.conf[b], self.conf[a]
        la, lb = self.labels.pop(a, None), self.labels.pop(b, None)
        if la:
            self.labels[b] = la
        if lb:
            self.labels[a] = lb
        self.draw()

    def draw(self):
        super().draw()
        scale, ox, oy = T.geometry(self.canvas)
        for s, name in self.labels.items():
            px, py = T.site_pos(s, scale, ox, oy)
            self.canvas.create_text(px, py, text=name, fill="#000",
                                    font=("TkDefaultFont", max(7, int(16 * scale)), "bold"))

    def load_result(self, r, label):
        self.labels = dict(r.get("labels", {}))
        super().load_result(r, label)

    def load_found(self, sign):
        rs = [r for r in self.found_results if sign in r["found"]]
        if not rs:
            self.test_info.config(text=self.test_head + f"\nKein Beispiel mit Vorzeichen {sign:+d}")
            return
        i = self.found_index[sign] % len(rs)
        self.found_index[sign] += 1
        r, steps = rs[i], rs[i]["found"][sign]
        self.load_result({"start": r["start"], "steps": list(steps), "sign": sign, "labels": r.get("labels", {}),
                          "n_swaps": len(steps) - r["hops"]}, f"Beispiel {i + 1}/{len(rs)} fuer {sign:+d}")

    # --- Zwei-Loch-Tests ---

    def two_inputs(self, allow_all=False):
        ma, mb = Z.parse_moves(self.path_a.get()), Z.parse_moves(self.path_b.get())
        order = self.order.get().strip().upper() or None
        try:
            off = Z.parse_offset(self.offset.get())
        except ValueError:
            self.test_info.config(text="B-Start relativ zu A als 'Zeile,Spalte' eingeben, z.B. -1,1")
            return None
        if not ma:
            self.test_info.config(text="Weg A leer: Folge aus R, L, U, D eingeben (z.B. RR)")
            return None
        if order == "ALLE" and not allow_all:
            self.test_info.config(text="Reihenfolge 'alle' gibt es nur für die Suche; hier eine feste Folge "
                                       "aus A/B oder leer (erst A, dann B) eingeben")
            return None
        if order and order != "ALLE" and sorted(order) != sorted("A" * len(ma) + "B" * len(mb)):
            self.test_info.config(text=f"Reihenfolge braucht genau {len(ma)}x A und {len(mb)}x B")
            return None
        return ma, mb, off, order

    def head_line(self, ma, mb, off, order, geometry):
        return (f"A {ma}, B {mb or '–'}, B-Start {off} relativ zu A, Reihenfolge "
                f"{order or 'erst A, dann B'}: {geometry}")

    def run_two_test(self, n_trials=100):
        inp = self.two_inputs()
        if inp is None:
            return
        ma, mb, off, order = inp
        results, skipped, geometry = [], 0, ""
        while len(results) < n_trials and skipped < 50 * n_trials:
            t = Z.random_trial(ma, mb, off, order, tries=2000)
            if t is None:
                break
            start, a, b, ha, hb = t
            o = order or Z.default_order(len(ha), len(hb))
            hops = Z.merge(ha, hb, o)
            runs = {w: Z.strategy_steps(start, hops, o, w) for w in Z.WHEN}
            if any(v is None for v in runs.values()):
                skipped += 1
                continue
            geometry = Z.classify(Z.path_sites(a, ma), Z.path_sites(b, mb))
            labels = {a: "A", b: "B"}
            variants = {}
            for w, steps in runs.items():
                _, s, nm = Z.step_counts(start, steps)
                variants[w] = {"start": start, "steps": steps, "sign": T.apply_steps(start, steps)[1],
                               "n_swaps": s, "n_minus": nm, "labels": labels}
            results.append(variants)
        if not results:
            self.test_info.config(text="Kein gültiger, reparierbarer Lauf gefunden (Weg verlässt das Brett, "
                                       "trifft ein anderes Loch oder lässt sich nicht per Paartausch reparieren)")
            return
        n = len(results)
        lines = [self.head_line(ma, mb, off, order, geometry) + f"; {n} Starts"
                 + (f" ({skipped} nicht reparierbare verworfen)" if skipped else "")]
        for w, name in Z.WHEN.items():
            plus = sum(r[w]["sign"] > 0 for r in results)
            nms = sorted({r[w]["n_minus"] for r in results})
            lines.append(f"{name}: +1: {plus}x, -1: {n - plus}x, n₋ ∈ {{{', '.join(map(str, nms))}}}")
        differs = lambda r, w: r["ende"]["sign"] != r[w]["sign"]  # noqa: E731
        lines.append(f"Vorzeichen im selben Start verschieden: Ende/nach A {sum(differs(r, 'nach_a') for r in results)}x, "
                     f"Ende/sofort {sum(differs(r, 'sofort') for r in results)}x   (\"Variante laden\" zeigt sie)")
        results.sort(key=lambda r: not (differs(r, "nach_a") or differs(r, "sofort")))
        self.two_results, self.variant_index = results, 0
        self.test_head = "\n".join(lines)
        self.test_info.config(text=self.test_head)

    def load_variant(self):
        if not self.two_results:
            self.test_info.config(text="Erst \"2 Löcher: Strategien vergleichen\" oder \"Beispiel Kreuzung\" ausführen")
            return
        k = self.variant_index
        self.variant_index += 1
        whens = list(Z.WHEN)
        i, w = (k // len(whens)) % len(self.two_results), whens[k % len(whens)]
        v = self.two_results[i][w]
        self.load_result(v, f"Start {i + 1}/{len(self.two_results)}, {Z.WHEN[w]}, n₋ = {v['n_minus']}")

    def run_two_search(self, n_trials=50):
        inp = self.two_inputs(allow_all=True)
        if inp is None:
            return
        ma, mb, off, order = inp
        results, capped, geometry = [], 0, ""
        while len(results) < n_trials:
            t = Z.random_trial(ma, mb, off, order, tries=2000)
            if t is None:
                break
            start, a, b, ha, hb = t
            geometry = Z.classify(Z.path_sites(a, ma), Z.path_sites(b, mb))
            found, _, cut = Z.search_orders(start, ha, hb, order)
            capped += cut
            results.append({"start": start, "hops": len(ha) + len(hb), "found": found, "labels": {a: "A", b: "B"}})
        if not results:
            self.test_info.config(text="Kein gültiger Lauf gefunden")
            return
        n = len(results)
        cnt = lambda f: sum(f(set(r["found"])) for r in results)  # noqa: E731
        plus, minus = cnt(lambda f: 1 in f), cnt(lambda f: -1 in f)
        both, none = cnt(lambda f: len(f) == 2), cnt(lambda f: not f)
        self.found_results, self.found_index = results, {1: 0, -1: 0}
        self.test_head = (self.head_line(ma, mb, off, order, geometry)
                          + f"; Suche nach Reparaturen (zwischen den Sprüngen erlaubt), {n} Starts:\n"
                          f"+1 erreichbar: {plus}x, -1 erreichbar: {minus}x, beide: {both}x, "
                          f"nur +1: {plus - both}x, nur -1: {minus - both}x, keines: {none}x"
                          + (f"   ({capped}x Suche abgebrochen)" if capped else "")
                          + "\n\"Beispiel +1\" / \"Beispiel -1\" laden gefundene Prozesse")
        self.test_info.config(text=self.test_head)

    def crossing_example(self):
        self.path_a.set("RR")
        self.path_b.set("DD")
        self.offset.set("-1,1")
        self.order.set("")
        self.run_two_test(n_trials=20)
        self.load_variant()


if __name__ == "__main__":
    root = tk.Tk()
    App2(root)
    root.mainloop()

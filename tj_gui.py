"""Interaktives 8x8-Gitter fuers t-J-Modell (keine Doppelbesetzung).

Rot = up, blau = down, leerer Kreis = Loch.
Hintergrund: gruen = Vorzeichen +1, lila = Vorzeichen -1.
Klick zwischen zwei Nachbar-Sites vertauscht deren Inhalt; das Vorzeichen
folgt der NetKet-Konvention (SpinOrbitalFermions): Orbital-Index
k = spin_idx * N + site  mit spin_idx(down) = 0, spin_idx(up) = 1,
site = x * L + y, Jordan-Wigner-Vorzeichen (-1)^(Anzahl besetzter Orbitale < k).
Start-Konfiguration hat Vorzeichen +1.
Graue Linien zeigen die Wege der Loecher seit dem letzten "Neu".
Rechts: "Sichern" (Taste s) legt ein Bild des Boards in einem 4x4-Raster ab;
Klick auf ein Bild macht es samt Lochspur zur neuen Start-Konfiguration.
Schalter "Néel im Hintergrund" zeigt die Néel-Farben hinter den Kreisen.
"Check Néel" markiert Spins auf dem falschen Untergitter orange.
"Lochspur löschen" entfernt die grauen Linien (mit "Zurück" rueckgaengig).
"""
import random
import tkinter as tk

L = 8
N = L * L
UP, DN, HOLE = 1, -1, 0

COL_UP, COL_DN = "#d62728", "#1f5fd6"
BG_PLUS, BG_MINUS = "#c8f0c8", "#dcc8f0"
CELL = 70
MARGIN = 40
BOARD_PX = 2 * MARGIN + (L - 1) * CELL
THUMB_GRID = 4
COL_TRAIL = "#9a9a9a"
COL_WRONG = "#ff8c00"
COL_UP_LIGHT, COL_DN_LIGHT = "#f2a0a0", "#9fbcf2"


def site(x, y):
    return x * L + y


def orbital(s, spin):
    return (1 if spin == UP else 0) * N + s


def to_occ(conf):
    """conf: Liste (N) mit +1/-1/0  ->  NetKet-Besetzungsvektor (2N)."""
    occ = [0] * (2 * N)
    for s, v in enumerate(conf):
        if v != HOLE:
            occ[orbital(s, v)] = 1
    return occ


def apply_c(occ, k, create):
    """Wendet c bzw. c^dagger auf Orbital k an. Gibt Vorzeichen zurueck (0 = verschwindet)."""
    if occ[k] == (1 if create else 0):
        return 0
    sign = -1 if sum(occ[:k]) % 2 else 1
    occ[k] = 1 if create else 0
    return sign


def swap_sign(conf, a, b):
    """Vorzeichen der Vertauschung der Inhalte von Site a und b (NetKet-Konvention)."""
    va, vb = conf[a], conf[b]
    if va == vb:
        return 1
    occ = to_occ(conf)
    if va == HOLE or vb == HOLE:
        src, dst = (a, b) if vb == HOLE else (b, a)
        spin = conf[src]
        ops = [(orbital(src, spin), False), (orbital(dst, spin), True)]
    else:
        # va auf a, vb auf b -> vb auf a, va auf b
        ops = [(orbital(b, vb), False), (orbital(a, vb), True),
               (orbital(a, va), False), (orbital(b, va), True)]
    sign = 1
    for k, create in ops:
        sign *= apply_c(occ, k, create)
    return sign


def neel_random_start():
    """Neel-Konfiguration, auf jeder Untergitter-Haelfte 4 zufaellige Loecher -> 28 up, 28 dn, 8 Loecher."""
    conf = [UP if (s // L + s % L) % 2 == 0 else DN for s in range(N)]
    for sub in (UP, DN):
        sites = [s for s in range(N) if conf[s] == sub]
        for s in random.sample(sites, L // 2):
            conf[s] = HOLE
    return conf


def neighbours(s):
    x, y = divmod(s, L)
    return [site(x + dx, y + dy) for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))
            if 0 <= x + dx < L and 0 <= y + dy < L]


def walk_steps(conf, h, moves):
    """Lochspruenge (a, b) fuer Zugfolge aus R/L/U/D ab Loch h; None, wenn Rand oder anderes Loch im Weg."""
    delta = {"R": (0, 1), "L": (0, -1), "U": (-1, 0), "D": (1, 0)}
    conf = list(conf)
    steps, pos = [], h
    for m in moves:
        x, y = divmod(pos, L)
        nx, ny = x + delta[m][0], y + delta[m][1]
        if not (0 <= nx < L and 0 <= ny < L) or conf[site(nx, ny)] == HOLE:
            return None
        nxt = site(nx, ny)
        steps.append((pos, nxt))
        conf[pos], conf[nxt] = conf[nxt], conf[pos]
        pos = nxt
    return steps


def apply_steps(conf, steps):
    """Fuehrt Vertauschungen aus; gibt (neue Konfiguration, Vorzeichen) zurueck."""
    conf, sign = list(conf), 1
    for a, b in steps:
        sign *= swap_sign(conf, a, b)
        conf[a], conf[b] = conf[b], conf[a]
    return conf, sign


def wrong_pairs(conf):
    wrong = neel_violations(conf)
    return [(a, b) for a in wrong for b in neighbours(a) if b in wrong and a < b
            and conf[a] != conf[b]]


def interleaved_steps(conf, hops, p=0.5):
    """Wie hops + repair_neel, aber vor jedem Sprung werden mit Wahrscheinlichkeit p zufaellig
    falsche Nachbarpaare repariert. Gibt (Schrittliste, Anzahl Reparaturen) oder None."""
    conf, steps, n_rep = list(conf), [], 0
    for a, b in hops:
        while random.random() < p:
            pairs = wrong_pairs(conf)
            if not pairs:
                break
            x, y = random.choice(pairs)
            conf[x], conf[y] = conf[y], conf[x]
            steps.append((x, y))
            n_rep += 1
        conf[a], conf[b] = conf[b], conf[a]
        steps.append((a, b))
    rest = repair_neel(conf)
    if rest is None:
        return None
    return steps + rest, n_rep + len(rest)


def search_signs(conf, hops, cap=20000):
    """Sucht Reparatur-Vertauschungen (benachbarte falsche up/dn-Paare), beliebig zwischen die Lochspruenge
    gemischt, so dass am Ende Neel erreicht ist. Gibt ({Vorzeichen: Schrittliste}, Zustaende, abgebrochen)."""
    found, seen = {}, set()
    stack = [(0, tuple(conf), 1, ())]
    while stack and len(found) < 2:
        k, c, sign, steps = stack.pop()
        if (k, c, sign) in seen:
            continue
        if len(seen) >= cap:
            return found, len(seen), True
        seen.add((k, c, sign))
        cl = list(c)
        pairs = wrong_pairs(cl)
        if k == len(hops) and not neel_violations(cl):
            found.setdefault(sign, steps)
            continue
        moves = list(pairs) + hops[k:k + 1]
        for a, b in moves:
            new = list(c)
            sg = sign * swap_sign(cl, a, b)
            new[a], new[b] = new[b], new[a]
            stack.append((k + (1 if (a, b) not in pairs else 0), tuple(new), sg, steps + ((a, b),)))
    return found, len(seen), False


def repair_neel(conf):
    """Repariert falsche Spins durch Vertauschen benachbarter falscher up/dn-Paare (zufaellige Wahl).
    Gibt Liste der Vertauschungen zurueck oder None, wenn nicht so reparierbar."""
    conf = list(conf)
    swaps = []
    while True:
        wrong = neel_violations(conf)
        if not wrong:
            return swaps
        pairs = [(a, b) for a in wrong for b in neighbours(a) if b in wrong and a < b
                 and conf[a] != conf[b]]
        if not pairs:
            return None
        a, b = random.choice(pairs)
        swaps.append((a, b))
        conf[a], conf[b] = conf[b], conf[a]


def neel_violations(conf):
    """Sites mit Spin auf dem falschen Untergitter (Loecher zaehlen nicht)."""
    return {s for s, v in enumerate(conf)
            if v != HOLE and v != (UP if (s // L + s % L) % 2 == 0 else DN)}


def bond_list():
    out = []
    for s in range(N):
        x, y = divmod(s, L)
        if y + 1 < L:
            out.append((s, site(x, y + 1)))
        if x + 1 < L:
            out.append((s, site(x + 1, y)))
    return out


BONDS = bond_list()


def geometry(canvas):
    """(scale, ox, oy): Board wird quadratisch in den Canvas eingepasst und zentriert."""
    w, h = max(canvas.winfo_width(), 1), max(canvas.winfo_height(), 1)
    scale = min(w, h) / BOARD_PX
    return scale, (w - BOARD_PX * scale) / 2, (h - BOARD_PX * scale) / 2


def site_pos(s, scale, ox, oy):
    x, y = divmod(s, L)
    return ox + (MARGIN + y * CELL) * scale, oy + (MARGIN + x * CELL) * scale


def render(canvas, conf, sign, trail, hover=None, neel=False, wrong=()):
    """Zeichnet das Board skaliert auf den Canvas (Hauptbrett und Vorschaubilder)."""
    scale, ox, oy = geometry(canvas)
    cell = CELL * scale

    def pos(s):
        return site_pos(s, scale, ox, oy)

    canvas.delete("all")
    bg = BG_PLUS if sign > 0 else BG_MINUS
    canvas.configure(bg=bg)
    for a, b in BONDS:
        hot = hover == (a, b)
        canvas.create_line(*pos(a), *pos(b), fill="#222" if hot else "#999",
                           width=max(2, 5 * scale) if hot else 1)
    if neel:
        rn = cell * 0.38
        for s in range(N):
            px, py = pos(s)
            x, y = divmod(s, L)
            canvas.create_oval(px - rn, py - rn, px + rn, py + rn, outline="",
                               fill=COL_UP_LIGHT if (x + y) % 2 == 0 else COL_DN_LIGHT)
    for a, b in trail:
        canvas.create_line(*pos(a), *pos(b), fill=COL_TRAIL, width=max(2, 6 * scale),
                           capstyle="round")
    r = cell * 0.3
    for s in wrong:
        px, py = pos(s)
        rw = cell * 0.43
        canvas.create_oval(px - rw, py - rw, px + rw, py + rw, outline=COL_WRONG,
                           width=max(2, 5 * scale), fill="")
    for s, v in enumerate(conf):
        px, py = pos(s)
        if v == HOLE:
            canvas.create_oval(px - r, py - r, px + r, py + r, outline="#333",
                               width=max(1, 2 * scale), fill=bg)
        else:
            canvas.create_oval(px - r, py - r, px + r, py + r, outline="",
                               fill=COL_UP if v == UP else COL_DN)


class App:
    def __init__(self, root):
        self.root = root
        root.title("t-J-Modell 8x8")
        side_px = min(2 * BOARD_PX, root.winfo_screenheight() - 160,
                      (root.winfo_screenwidth() - 30) // 2)
        root.geometry(f"{2 * side_px + 30}x{side_px + 40}")
        root.columnconfigure(0, weight=1)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(0, weight=1)

        main = tk.Frame(root)
        main.grid(row=0, column=0, sticky="nsew")
        main.rowconfigure(0, weight=1)
        main.columnconfigure(0, weight=1)
        self.canvas = tk.Canvas(main, width=10, height=10, highlightthickness=0)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        bar = tk.Frame(main)
        bar.grid(row=1, column=0, sticky="ew")
        tk.Button(bar, text="Neu (zufällig)", command=self.new).pack(side="left")
        tk.Button(bar, text="Zurück", command=self.undo).pack(side="left")
        tk.Button(bar, text="Vor", command=self.redo).pack(side="left")
        tk.Button(bar, text="Alles zurück", command=self.undo_all).pack(side="left")
        tk.Button(bar, text="Check Néel", command=self.check_neel).pack(side="left")
        tk.Button(bar, text="Lochspur löschen", command=self.clear_trail).pack(side="left")
        self.neel_var = tk.BooleanVar(value=False)
        tk.Checkbutton(bar, text="Néel im Hintergrund", variable=self.neel_var,
                       command=self.toggle_neel).pack(side="left", padx=6)
        self.info = tk.Label(bar, anchor="e")
        self.info.pack(side="right", padx=8)
        bar2 = tk.Frame(main)
        bar2.grid(row=2, column=0, sticky="ew")
        self.path_var = tk.StringVar(value="RRULLD")
        tk.Entry(bar2, textvariable=self.path_var, width=14).pack(side="left", padx=(0, 4))
        tk.Button(bar2, text="Test, Reparatur am Ende", command=self.run_loop_test).pack(side="left")
        tk.Button(bar2, text="Test, Reparatur zwischendurch",
                  command=lambda: self.run_loop_test(interleave=0.5)).pack(side="left")
        tk.Button(bar2, text="Test, Reparatur sofort",
                  command=lambda: self.run_loop_test(interleave=1.0)).pack(side="left")
        tk.Button(bar2, text="Suche +1/-1 Reparatur", command=self.run_sign_search).pack(side="left")
        tk.Button(bar2, text="Beispiel +1", command=lambda: self.load_found(1)).pack(side="left")
        tk.Button(bar2, text="Beispiel -1", command=lambda: self.load_found(-1)).pack(side="left")
        tk.Button(bar2, text="Abweichung laden", command=self.load_deviation).pack(side="left")
        tk.Button(bar2, text="Zufallsbeispiel", command=self.load_random_example).pack(side="left")
        bar3 = tk.Frame(main)
        bar3.grid(row=3, column=0, sticky="ew")
        tk.Label(bar3, fg="#555", anchor="w",
                 text="Wegeingabe: R/L/U/D = rechts/links/oben/unten ab zufälligem Loch, z.B. RRULLD; "
                      "danach Néel-Reparatur durch Vertauschen.").pack(anchor="w")
        self.test_info = tk.Label(bar3, anchor="w", justify="left", text="")
        self.test_info.pack(anchor="w")
        self.test_results = []
        self.found_results = []
        self.test_head = ""
        self.deviants = []
        self.dev_index = 0

        side = tk.Frame(root)
        side.grid(row=0, column=1, sticky="nsew", padx=(6, 0))
        side.rowconfigure(1, weight=1)
        side.columnconfigure(0, weight=1)
        btns = tk.Frame(side)
        btns.grid(row=0, column=0, sticky="ew")
        tk.Button(btns, text="Sichern (s)", command=self.save).pack(side="left")
        tk.Button(btns, text="Bilder löschen", command=self.clear_saved).pack(side="left")
        grid = tk.Frame(side)
        grid.grid(row=1, column=0, sticky="nsew")
        for k in range(THUMB_GRID):
            grid.rowconfigure(k, weight=1, uniform="r")
            grid.columnconfigure(k, weight=1, uniform="c")
        self.thumbs = []
        for i in range(THUMB_GRID * THUMB_GRID):
            t = tk.Canvas(grid, width=10, height=10, bg="#eee",
                          highlightthickness=1, highlightbackground="#bbb")
            t.grid(row=i // THUMB_GRID, column=i % THUMB_GRID, sticky="nsew")
            t.bind("<Configure>", lambda e, i=i: self.draw_thumb(i))
            t.bind("<Button-1>", lambda e, i=i: self.load_saved(i))
            self.thumbs.append(t)
        self.saved = []

        self.canvas.bind("<Configure>", lambda e: self.draw())
        self.canvas.bind("<Button-1>", self.click)
        self.canvas.bind("<Motion>", self.hover)
        self.canvas.bind("<Leave>", lambda e: self.set_hover(None))
        root.bind("s", lambda e: self.save())
        self.hover_bond = None
        self.new()

    def new(self):
        self.conf = neel_random_start()
        self.sign = 1
        self.trail = []
        self.history = []
        self.redo_stack = []
        self.wrong = None
        self.draw()

    def check_neel(self):
        self.wrong = neel_violations(self.conf)
        self.draw()

    def snapshot(self):
        return (list(self.conf), self.sign, list(self.trail))

    def push_history(self):
        self.history.append(self.snapshot())
        self.redo_stack = []
        self.wrong = None

    def restore(self, state):
        self.conf, self.sign, self.trail = list(state[0]), state[1], list(state[2])
        self.wrong = None
        self.draw()

    def toggle_neel(self):
        self.draw()
        self.draw_saved()

    def load_saved(self, i):
        if i < len(self.saved):
            conf, sign, trail = self.saved[i]
            self.conf, self.sign, self.trail = list(conf), sign, list(trail)
            self.history = []
            self.redo_stack = []
            self.wrong = None
            self.draw()

    def clear_trail(self):
        if self.trail:
            self.push_history()
            self.trail = []
            self.draw()

    def undo_all(self):
        while self.history:
            self.undo()

    def undo(self):
        if self.history:
            self.redo_stack.append(self.snapshot())
            self.restore(self.history.pop())

    def redo(self):
        if self.redo_stack:
            self.history.append(self.snapshot())
            self.restore(self.redo_stack.pop())

    def save(self):
        if len(self.saved) == len(self.thumbs):
            self.saved.pop(0)
        self.saved.append((list(self.conf), self.sign, list(self.trail)))
        self.draw_saved()

    def clear_saved(self):
        self.saved = []
        self.draw_saved()

    def draw_thumb(self, i):
        t = self.thumbs[i]
        if i < len(self.saved):
            render(t, *self.saved[i], neel=self.neel_var.get())
        else:
            t.delete("all")
            t.configure(bg="#eee")

    def draw_saved(self):
        for i in range(len(self.thumbs)):
            self.draw_thumb(i)

    def bond_at(self, px, py):
        """Naechstgelegene Bindung, wenn der Klick nahe ihrer Mitte liegt."""
        geo = geometry(self.canvas)
        best, best_d = None, 0.3 * CELL * geo[0]
        for s, t in BONDS:
            (x1, y1), (x2, y2) = site_pos(s, *geo), site_pos(t, *geo)
            d = ((px - (x1 + x2) / 2) ** 2 + (py - (y1 + y2) / 2) ** 2) ** 0.5
            if d < best_d:
                best, best_d = (s, t), d
        return best

    def hover(self, e):
        self.set_hover(self.bond_at(e.x, e.y))

    def set_hover(self, bond):
        if bond != self.hover_bond:
            self.hover_bond = bond
            self.draw()

    def run_loop_test(self, n_trials=100, interleave=None):
        moves = "".join(c for c in self.path_var.get().upper() if c in "RLUD")
        if not moves:
            self.test_info.config(text="Weg leer: Folge aus R, L, U, D eingeben (z.B. RRULLD)")
            return
        results, fails, attempts = [], 0, 0
        while len(results) < n_trials and attempts < 100000:
            attempts += 1
            start = neel_random_start()
            h = random.choice([s for s in range(N) if start[s] == HOLE])
            hops = walk_steps(start, h, moves)
            if hops is None:
                continue
            if interleave is not None:
                res = interleaved_steps(start, hops, interleave)
            else:
                swaps = repair_neel(apply_steps(start, hops)[0])
                res = None if swaps is None else (hops + swaps, len(swaps))
            if res is None:
                fails += 1
                continue
            steps, n_swaps = res
            _, sign = apply_steps(start, steps)
            results.append({"start": start, "steps": steps, "sign": sign, "n_swaps": n_swaps})
        if not results:
            self.test_info.config(text="Kein gueltiger/reparierbarer Lauf gefunden")
            return
        n = len(results)
        n_plus = sum(r["sign"] > 0 for r in results)
        major = 1 if 2 * n_plus >= n else -1
        self.deviants = [r for r in results if r["sign"] != major]
        lines = [f"Weg {moves}{' (Reparatur am Ende)' if interleave is None else ' (Reparatur sofort)' if interleave >= 1 else ' (Reparatur zwischendurch)'}: {n} Laeufe, Vorzeichen +1: {n_plus}x, -1: {n - n_plus}x"
                 + (f"   ({fails} nicht reparierbare Laeufe verworfen)" if fails else "")]
        lines.append("immer das gleiche Vorzeichen" if not self.deviants
                     else f"{len(self.deviants)} Abweichler (Vorzeichen {-major:+d}, Mehrheit {major:+d})")
        self.test_results, self.dev_index = results, 0
        self.test_head = "\n".join(lines)
        self.test_info.config(text=self.test_head)

    def run_sign_search(self, n_trials=100):
        moves = "".join(c for c in self.path_var.get().upper() if c in "RLUD")
        if not moves:
            self.test_info.config(text="Weg leer: Folge aus R, L, U, D eingeben (z.B. RRULLD)")
            return
        results, attempts, capped = [], 0, 0
        while len(results) < n_trials and attempts < 100000:
            attempts += 1
            start = neel_random_start()
            h = random.choice([s for s in range(N) if start[s] == HOLE])
            hops = walk_steps(start, h, moves)
            if hops is None:
                continue
            found, _, cut = search_signs(start, hops)
            capped += cut
            results.append({"start": start, "hops": len(hops), "found": found})
        if not results:
            self.test_info.config(text="Kein gueltiger Lauf gefunden")
            return
        n = len(results)
        cnt = lambda f: sum(f(set(r["found"])) for r in results)
        plus, minus = cnt(lambda f: 1 in f), cnt(lambda f: -1 in f)
        both, none = cnt(lambda f: len(f) == 2), cnt(lambda f: not f)
        self.found_results = results
        self.found_index = {1: 0, -1: 0}
        self.test_head = (f"Weg {moves}, Suche nach Reparatur-Vertauschungen (zwischen den Spruengen erlaubt), "
                          f"{n} Starts:\n+1 erreichbar: {plus}x, -1 erreichbar: {minus}x, "
                          f"beide: {both}x, nur +1: {plus - both}x, nur -1: {minus - both}x, keines: {none}x"
                          + (f"   ({capped}x Suche abgebrochen)" if capped else ""))
        self.test_info.config(text=self.test_head)

    def load_found(self, sign):
        rs = [r for r in self.found_results if sign in r["found"]]
        if not rs:
            self.test_info.config(text=self.test_head + f"\nKein Beispiel mit Vorzeichen {sign:+d}")
            return
        i = self.found_index[sign] % len(rs)
        self.found_index[sign] += 1
        r, steps = rs[i], rs[i]["found"][sign]
        self.load_result({"start": r["start"], "steps": list(steps), "sign": sign,
                          "n_swaps": len(steps) - r["hops"]}, f"Beispiel {i + 1}/{len(rs)} fuer {sign:+d}")

    def load_result(self, r, label):
        self.conf, self.sign, self.trail = list(r["start"]), 1, []
        self.history, self.redo_stack, self.wrong = [], [], None
        for a, b in r["steps"]:
            self.click_swap(a, b)
        self.undo_all()
        self.test_info.config(
            text=self.test_head + f"\n{label}: Endvorzeichen {r['sign']:+d}, "
            f"{r['n_swaps']} Reparatur-Vertauschungen (\"Vor\" spielt Weg + Reparatur ab)")

    def load_deviation(self):
        if self.deviants:
            self.dev_index += 1
            n = (self.dev_index - 1) % len(self.deviants) + 1
            self.load_result(self.deviants[n - 1], f"Abweichler {n}/{len(self.deviants)}")

    def load_random_example(self):
        if self.test_results:
            self.load_result(random.choice(self.test_results), "Zufallsbeispiel")

    def click_swap(self, a, b):
        self.push_history()
        if HOLE in (self.conf[a], self.conf[b]):
            self.trail.append((a, b))
        self.sign *= swap_sign(self.conf, a, b)
        self.conf[a], self.conf[b] = self.conf[b], self.conf[a]
        self.draw()

    def click(self, e):
        bond = self.bond_at(e.x, e.y)
        if bond is None:
            return
        a, b = bond
        if self.conf[a] == self.conf[b]:
            return
        self.click_swap(a, b)

    def draw(self):
        render(self.canvas, self.conf, self.sign, self.trail, hover=self.hover_bond,
               neel=self.neel_var.get(), wrong=self.wrong or ())
        n_up = self.conf.count(UP)
        n_dn = self.conf.count(DN)
        n_wrong = len(neel_violations(self.conf))
        check = f"Néel falsch: {n_wrong}    "
        self.info.config(text=f"{check}Vorzeichen {'+' if self.sign > 0 else '-'}1    "
                              f"up {n_up}  dn {n_dn}  Löcher {N - n_up - n_dn}")


if __name__ == "__main__":
    root = tk.Tk()
    App(root)
    root.mainloop()

# Zwei Löcher: Reparatur-Zeitpunkt, Verschränkung, Paarbildung

Dieses Verzeichnis erweitert die Vorzeichen-Untersuchung aus dem Hauptverzeichnis
([../README.md](../README.md), [../tj_gui.py](../tj_gui.py)) auf **zwei Löcher A und B**. Die Frage:
Kann die Abhängigkeit des Vorzeichens vom *Zeitpunkt* der Néel-Reparatur zu einer Verschränkung von
Löchern oder Elektronenpaaren führen?

Kurzfassung:

- Das Repo-Vorzeichen eines Prozesses ist genau das Vorzeichen seines Gewichts in exp(−τH)
  (Abschnitt 1). Prozesse, die sich nur im Reparatur-Zeitpunkt unterscheiden, interferieren also.
- **Satz:** Besuchen die Wege aller Löcher zusammen keinen Platz zweimal, ist das Vorzeichen eindeutig
  (Abschnitt 2). Ein Zwei-Loch-Vorzeichen kann deshalb nur dann vom Reparatur-Zeitpunkt abhängen,
  wenn sich die Wege kreuzen oder ein Weg eine Schleife hat.
- **Kreuzung:** Repariert J die Spur von A, *bevor* B sie kreuzt, ist das Vorzeichen entgegengesetzt
  zu „Reparatur danach“ (in 100 von 100 Starts, Abschnitt 4). Folgt B dagegen direkt hinter A
  („Spur teilen“, Raupe), ist das Vorzeichen eindeutig.
- **Kohärenz und ODLRO (Abschnitt 7):** Für geschlossene Prozesse gilt Vorzeichen = (−1)^s · sgn(π_Loch)
  (s = Zahl der Spin-Vertauschungen). Die reine Raupe ist dadurch frustriert (jede Wende ist ein
  Lochaustausch, −1). Paarschritte mit genau einer Reparatur (Pivot, Domino) sind dagegen ganz unfrustriert,
  und zwar im d-Wellen-Muster: Die Reparatur gleicht das Austausch-Vorzeichen aus. Ein einzelnes Loch wird
  durch solche Schritte kohärent, bleibt aber frustriert (−1 pro Schritt). Lockert man Néel und erlaubt
  J-Fluktuationen (Vertauschungen, die keine Reparaturen sind), kommt auch für Paare Frustration zurück
  (Verhältnis 0,92 statt 1), und zwar dort, wo das Paar eine Fluktuation in seiner Nähe umordnet.
- **Exakte Diagonalisierung (4×4, periodisch):** Mit den Vorzeichen sind zwei Löcher gebunden
  (E_b = −0,35 t bei J/t = 0,4), in der vorzeichenfreien Version |H| nicht (+0,19 t). Die Vorzeichen
  erzeugen außerdem eine starke Verschränkung zwischen Lochpositionen und Spins (Entropie ≈ 3 statt
  < 0,9) und schwächen die Singulett-Bindungen neben dem Loch (Abschnitt 6).

## Inhalt

| Datei | Zweck | Start |
|---|---|---|
| `zwei_loch.py` | Kernfunktionen: Wege, Sprung-Reihenfolgen, Reparatur-Strategien, Suche (baut auf `../tj_gui.py` auf) | – |
| `tj_gui_2loch.py` | GUI: alles aus `tj_gui.py` plus zwei Leisten für zwei Löcher | `python tj_gui_2loch.py` |
| `tj_gui_web.html` | dieselben Funktionen als Webseite (JavaScript, für Touch bedienbar, z.B. auf dem iPad); Vorzeichen gegen `tj_gui.py` geprüft (400/400 Prozesse) | im Browser öffnen |
| `experimente.py` | reproduziert die Tabellen aus Abschnitt 4 (Markdown-Ausgabe) | `python experimente.py` |
| `kohaerenz.py` | Regel für geschlossene Prozesse, Paar- und Einzelloch-Bewegung, Frustration (Abschnitt 7, etwa 20 s) | `python kohaerenz.py` |
| `fluktuationen.py` | Paarbewegung mit einer J-Fluktuation, vollständige Prozessaufzählung (Abschnitt 7.6, etwa 3 min) | `python fluktuationen.py` |
| `ed_tj.py` | exakte Diagonalisierung, t-J-Modell gegen vorzeichenfreie Version (etwa 17 min) | `python ed_tj.py --out ed_ergebnisse.md` |
| `ed_ergebnisse.md` | vollständige Ausgabe von `ed_tj.py` (4×4 periodisch, J/t = 0,2 / 0,4 / 0,6 / 1,0) | – |

`tj_gui_2loch.py`, `zwei_loch.py` und `experimente.py` brauchen nur tkinter (wie `tj_gui.py`),
`ed_tj.py` braucht numpy (≥ 2.0) und scipy, für `--netket` zusätzlich NetKet.

## Bedienung der GUI

`tj_gui_2loch.py` enthält das komplette Hauptprogramm; darunter kommen zwei Leisten hinzu:

- **Weg A**, **Weg B**: Folgen aus `R`, `L`, `U`, `D` wie im Hauptprogramm.
- **B-Start relativ zu A (Zeile,Spalte)**: z.B. `-1,1` = eine Zeile höher, eine Spalte weiter rechts.
  Der Start von A wird zufällig gewählt; A und B werden zu Löchern gemacht (Rest wie
  `neel_random_start`: 28 up, 28 dn, 8 Löcher).
- **Reihenfolge**: Folge aus `A`/`B`, z.B. `ABAB`; leer = erst alle Sprünge von A, dann B;
  `alle` = alle Reihenfolgen (nur für die Suche).
- **2 Löcher: Strategien vergleichen**: 100 zufällige Starts. In jedem Start werden dieselben Sprünge mit
  drei Strategien ausgeführt: *Reparatur am Ende*, *Reparatur nach A* (alles Reparierbare nach dem
  letzten Sprung von A, bevor B weiterläuft) und *Reparatur sofort* (vor jedem Sprung). Ausgabe:
  Vorzeichen, Werte von n₋ und wie oft sich die Strategien *im selben Start* unterscheiden.
  Starts, die sich mit einer Strategie nicht per Paartausch reparieren lassen, werden verworfen.
- **Variante laden**: lädt nacheinander Start und Strategie (Starts mit verschiedenen Vorzeichen zuerst);
  **Vor** spielt Weg und Reparatur ab. Die Löcher sind mit A und B beschriftet.
- **2 Löcher: Suche +1/-1**: für 50 Starts erschöpfende Suche nach Reparaturen (beliebig zwischen die
  Sprünge gemischt) mit Endvorzeichen +1 bzw. −1; **Beispiel +1** / **Beispiel -1** (zweite Leiste)
  laden die gefundenen Prozesse.
- **Beispiel Kreuzung**: setzt A = `RR`, B = `DD`, B-Start `-1,1` und vergleicht die Strategien.

Koordinaten in diesem README: (Zeile, Spalte), `U` = Zeile − 1, `R` = Spalte + 1, wie in `tj_gui.py`.

## 1. Repo-Vorzeichen = Vorzeichen des Quantengewichts

In `tj_gui.py` ist ein Lochsprung der Operator c†_dst c_src. Für den Paartausch gilt (geprüft in
`zwei_loch.j_term_check`, 3000/3000):

    swap_sign = −⟨f| S⁺_i S⁻_j + S⁻_i S⁺_j |i⟩ .

Mit H = −t Σ c†c + J Σ (S_i·S_j − n_i n_j/4) hat ein Prozess aus m Sprüngen und s Paartauschen in der
Entwicklung von exp(−τH) das Gewicht

    t^m · (J/2)^s · (positive Diagonalfaktoren) · Repo-Vorzeichen .

Zwei Prozesse mit gleichem Anfangs- und Endzustand haben ein relatives Vorzeichen, das nicht von der
Fermion-Ordnung abhängt (eine andere Ordnung multipliziert beide mit demselben Faktor). „Reparatur zu
verschiedenen Zeiten“ heißt also: verschiedene Weltlinien, die in der Quantendynamik alle zugleich
beitragen und miteinander interferieren. Das vorzeichenfreie Modell |H| (Abschnitt 6) enthält dieselben
Prozesse mit denselben Beträgen, aber ohne diese Interferenz.

## 2. Satz für zwei (und mehr) Löcher

Die Schritte 1–3 der Beweis-Skizze im Haupt-README gelten unverändert für mehrere Löcher, wenn m die
Gesamtzahl der Lochsprünge ist und F vom Start und den Endplätzen der (benannten) Löcher abhängt:

    Vorzeichen = (−1)^(m+s) · F = (−1)^(3m/2 + n₋) · F ,   s = m/2 − n₋ ,

n₋ = Zahl der Sprünge, bei denen ein Loch einen *falschen* (noch nicht reparierten) Spin bewegt.
Der Reparatur-Zeitpunkt wirkt nur über n₋.

**Satz.** Besuchen die Wege aller Löcher zusammen keinen Platz zweimal (Startplätze mitgezählt),
ist n₋ = 0 für jede Reparatur-Strategie und jede Sprung-Reihenfolge; das Vorzeichen ist eindeutig.

*Beweis* wie Schritt 4 im Haupt-README: Invariante „alle falschen Spins stehen auf bereits besuchten
Plätzen“ (ein Sprung legt den Spin auf den eben verlassenen Platz, Reparaturen tauschen nur falsche Spins).
Das gilt für beliebig viele Löcher. Der Zielplatz jedes Sprungs ist noch von keinem Loch besucht, enthält
also seinen ursprünglichen, richtigen Spin. □

Numerisch (Zufallstest, Abschnitt 4.5): 302 reparierbare Fälle ohne doppelten Platz, immer genau ein
Vorzeichen.

**Die Umkehrung gilt für zwei Löcher nicht.** Ein doppelt besuchter Platz (Kreuzung von A und B oder
Schleife eines Weges) ist notwendig, aber nicht hinreichend: von 35 reparierbaren Kreuzungsfällen ohne
Schleife waren nur 9 zweideutig (Abschnitt 4.5; warum, erklärt Abschnitt 3).

## 3. Paritätsregel (Erklärung, kein Beweis)

Das Vorzeichen hängt nur von der Parität von n₋ ab. Wird ein falsches Paar (u, v) repariert, *bevor*
Löcher wieder auf u oder v springen, statt danach, springen diese Löcher dort richtige statt falsche
Spins. n₋ ändert sich dabei um die Zahl dieser späteren Sprünge auf u, v.

- **Kreuzung** (A `RR`, B `DD` quer durch die Spur): B springt auf genau einen Platz des Paares →
  n₋ ändert sich um 1 → Vorzeichen kippt.
- **B folgt A** (B direkt hinter A): B springt auf beide Plätze des Paares → n₋ ändert sich um 2 →
  kein Vorzeichenwechsel. B „repariert“ die Spur von A selbst.

Die Strategie-Tabelle (4.2) bestätigt das: Kreuzung n₋ = 1 / 0 / 0 (verschiedene Vorzeichen),
B folgt A n₋ = 2 / 0 / 0 bzw. 4 / {0, 2} / 0 (gleiche Vorzeichen).

## 4. Numerische Ergebnisse (`experimente.py`)

### 4.1 Kreuzungsbeispiel

Volle Néel-Konfiguration, nur die Löcher A (4,2) und B (3,3). A läuft `RR`, B läuft `DD`; die Wege
kreuzen sich bei (4,3). Gleicher Start, gleiche Sprünge, gleicher Endzustand:

| Ablauf | m | s | n₋ | Vorzeichen |
|---|---|---|---|---|
| Spur von A repariert, **bevor** B kreuzt | 4 | 2 | 0 | +1 |
| Spur von A repariert, **nachdem** B gekreuzt hat | 4 | 1 | 1 | −1 |

### 4.2 Strategievergleich auf zufälligen Starts

100 zufällige Néel-Starts je Zeile (8 Löcher), Sprünge erst A, dann B. Gezählt wird, wie oft sich das
Endvorzeichen zwischen zwei Strategien *im selben Start* unterscheidet. „verworfen“ = Starts, in denen
mindestens eine Strategie nicht per Paartausch reparierbar war.

| Fall | A | B | B-Start rel. zu A | Ende ≠ nach A | Ende ≠ sofort | n₋ (Ende / nach A / sofort) | verworfen |
|---|---|---|---|---|---|---|---|
| Kreuzung | RR | DD | (−1, 1) | 100/100 | 100/100 | 1 / 0 / 0 | 0 |
| Kreuzung | RRRR | DD | (−1, 2) | 100/100 | 100/100 | 1 / 0 / 0 | 143 |
| Kreuzung, A mit Ecke | RRUU | LL | (−1, 3) | 100/100 | 100/100 | 1 / 0 / 0 | 98 |
| B folgt A | RR | RR | (0, −1) | 0/100 | 0/100 | 2 / 0 / 0 | 0 |
| B folgt A, Ecke | RRUU | RRRU | (0, −1) | 0/100 | 0/100 | 4 / 0,2 / 0 | 27 |
| getrennt | RR | DD | (−3, 3) | 0/100 | 0/100 | 0 / 0 / 0 | 0 |

Die absoluten Vorzeichen sind bei nicht geschlossenen Wegen je Start verschieden (etwa halb +1, halb −1,
wie `RUULDDRU` im Haupt-README); aussagekräftig ist der Vergleich im selben Start.

### 4.3 Kontrollierte Fälle

Volle Néel-Konfiguration mit nur A und B, *alle* Sprung-Reihenfolgen, Reparaturen beliebig dazwischen
(erschöpfende Suche):

| Fall | A | B | Geometrie | Ergebnis je Reihenfolge | erreichbar insgesamt |
|---|---|---|---|---|---|
| getrennt | (4,2) RR | (1,5) DD | kein Platz doppelt | nur −1: 6 | −1 |
| getrennt | (5,1) RRRR | (0,6) DDDD | kein Platz doppelt | nur +1: 70 | +1 |
| Kreuzung | (4,2) RR | (3,3) DD | Kreuzung | beide: 2 | +1 und −1 |
| Kreuzung | (4,1) RRRR | (2,3) DDDD | Kreuzung | nur −1: 10, keines: 24 | −1 |
| B folgt A, gerade | (4,3) RR | (4,2) RR | Kreuzung | nur +1: 2 | +1 |
| B folgt A, gerade | (4,3) RRR | (4,2) RRR | Kreuzung | nur +1: 5 | +1 |
| B folgt A, Ecke | (4,3) RU | (4,2) RR | Kreuzung | nur +1: 2 | +1 |
| B folgt A, Ecke | (4,3) RRUU | (4,2) RRRU | Kreuzung | nur +1: 14 | +1 |
| A umrundet einen Spin | (3,3) RRDDLLUU | (0,7) – | Schleife | beide: 1 | +1 und −1 |
| A umrundet Loch B | (3,3) RRDDLLUU | (4,4) – | Schleife | beide: 1 | +1 und −1 |
| Austausch über Diagonale | (3,3) RD | (4,4) LU | Kreuzung | nur −1: 4 | −1 |
| kein Austausch (hin und zurück) | (3,3) RL | (4,4) LR | Schleife | nur +1: 6 | +1 |

„keines“ = in dieser Reihenfolge nicht per Paartausch bis Néel reparierbar.

### 4.4 Umrundung: Löcher sehen einander nicht direkt

A umrundet den 3×3-Block (`RRDDLLUU`, geschlossen, Reparatur sofort). In der Mitte steht einmal ein
Spin, einmal Loch B:

| Fall | Vorzeichen | m | s | n₋ |
|---|---|---|---|---|
| A umrundet einen Spin | +1 | 8 | 4 | 0 |
| A umrundet Loch B | +1 | 8 | 4 | 0 |

Das Vorzeichen tragen die unreparierten (falschen) Spins, nicht die Löcher selbst (nur dieses Beispiel
gerechnet).

### 4.5 Zufallstest

600 Fälle (Seed 7): zufällige Wege (A 1–5, B 0–5 Sprünge, Gesamtlänge gerade, kein direktes
Zurückgehen), Starts im Zentrum einer vollen Néel-Konfiguration mit nur A und B, zufällige
Sprung-Reihenfolge, Reparaturen beliebig dazwischen:

| Geometrie | beide Vorzeichen | eines | nicht reparierbar |
|---|---|---|---|
| kein Platz doppelt | 0 | 302 | 237 |
| Kreuzung A/B | 9 | 26 | 10 |
| Schleife (nur ein Loch besucht einen Platz zweimal) | 9 | 0 | 6 |
| Kreuzung und Schleife | 0 | 0 | 1 |

## 5. Mechanismen für Verschränkung

**M1 – Kreuzende Spuren (gegenseitige Statistik).** Loch A hinterlässt eine Spur falscher Spins. Kreuzt B
sie, bevor J sie repariert, kippt das Vorzeichen (4.1, 4.2). In der Zwei-Loch-Amplitude ψ(x_A, x_B)
tragen kreuzende Weltlinien deshalb ein Vorzeichen, das vom Zeitpunkt des Kreuzens *relativ* zum
Zeitpunkt der Reparatur abhängt, also von beiden Lochgeschichten zugleich. Eine solche Phase wirkt wie ein
Controlled-Z-Gatter: ψ zerfällt dann im Allgemeinen nicht mehr in φ(x_A)·χ(x_B). Für B ist die
unreparierte Spur von A wie eine π-Flusslinie, die J mit der Zeit abbaut; die Löcher selbst sehen
einander nicht (4.4).

**M2 – Spur teilen → Lochpaar.** Folgt B direkt hinter A, repariert B die Spur von A selbst; die
Paarbewegung kostet keine Spur-Energie und hat keinen Vorzeichenkonflikt, solange das Paar nicht wendet
(4.2, 4.3). Einschränkung aus Abschnitt 7: Bei einer Wende tauschen Kopf und Schwanz die Plätze, das kostet
−1; reine Raupenbewegung ist deshalb frustriert. Erst Paarschritte mit einer Reparatur heben das auf. Unabhängig laufende
Löcher mit kreuzenden Spuren haben dagegen Weltlinien beider Vorzeichen, die sich teilweise auslöschen.
Das bevorzugt den Paarkanal; ein gebundenes Paar ψ(R, r) = Φ(R)·φ(r) ist in den Einzelpositionen
verschränkt.

**M3 – Elektronenpaar speichert den Zeitpunkt (Singulett/Triplett).** Die Reparatur ist der J-Term, also
der Austausch eines ↑↓-Elektronenpaars. Die Bindungskorrelation ⟨S⁺_p S⁻_q + h.c.⟩ im Endzustand wird
durch das relative Vorzeichen zweier Weltlinien bestimmt: „Paar (p,q) am Ende getauscht“ gegen „früher
getauscht, bevor ein Loch vorbeikam“. Gleiches Vorzeichen ergibt einen singulettartigen Beitrag
(J-Energiegewinn), entgegengesetztes einen triplettartigen. Prägen zwei Lochgeschichten mit gleichem
Endplatz verschiedene Paarkorrelationen auf, sind Lochbewegung und Spins verschränkt. Wird der Zeitpunkt
anderswo festgehalten (dritte Teilchen, Phononen, Messung), verschwindet die Interferenz.

**M4 – Austauschstatistik (spekulativ).** Der Platztausch zweier Löcher über die Diagonale einer
Plakette liefert eindeutig −1, „kein Tausch“ +1 (4.3). Laufen die Löcher dabei über unreparierte Spuren,
käme nach M1 ein zeitpunktabhängiges Vorzeichen hinzu; die Symmetrie eines Lochpaares könnte dann von
J/t abhängen. Nur dieses eine Beispiel ist gerechnet.

## 6. Exakte Diagonalisierung (`ed_tj.py`)

Exakte Diagonalisierung des t-J-Modells auf dem **4×4-Gitter mit periodischen Rändern**, t = 1,
J/t = 0,2 / 0,4 / 0,6 / 1,0, mit 0, 1 und 2 Löchern (bis 411 840 Zustände, Lanczos aus scipy).
Verglichen werden

- **tJ**: das t-J-Modell H = −t Σ P(c†_iσ c_jσ + h.c.)P + J Σ (S_i·S_j − n_i n_j/4),
  in derselben Fermion-Ordnung wie NetKet und `tj_gui.py`;
- **frei**: |H|, dieselben Matrixelemente ohne Vorzeichen (alle Nebendiagonal-Elemente −|H_fi|).
  |H| enthält dieselben Prozesse mit denselben Beträgen wie H (Abschnitt 1), aber keine Interferenz.

Prüfungen: `python ed_tj.py --netket` vergleicht die komplette Hamiltonmatrix (3×3, offen und
periodisch, 1 und 2 Löcher) mit einem aus NetKet-Fermion-Operatoren gebauten t-J-Operator: maximale
Abweichung 2·10⁻¹⁶, gleiche Grundzustandsenergien. `python ed_tj.py --pruefen` vergleicht 15 831
Nebendiagonal-Elemente auf 4×4 mit `tj_gui.swap_sign`: alle gleich. Ohne Loch ergibt sich
E(0) = −11,2285 J − 8 J, also −0,7018 J pro Platz für den Heisenberg-Anteil (bekannter 4×4-Wert).

Technisches: tJ wird im Sektor Sz = 0 bzw. 1/2 gerechnet (SU(2)-Symmetrie). |H| ist nicht
SU(2)-symmetrisch, daher wurden alle Sektoren Sz ≥ 0 gerechnet; der Grundzustand lag immer bei
Sz = 0 bzw. 1/2. Entartete Grundzustände (tJ: 1 Loch 6-fach, 2 Löcher 3-fach) werden gemittelt.
Auf dem 4×4-Gitter sind die Lochabstände (1,1) und (2,0) durch eine Zusatzsymmetrie gleichwertig.
Laufzeit des Standardlaufs: etwa 17 Minuten. Vollständige Ausgabe: [`ed_ergebnisse.md`](ed_ergebnisse.md).

### 6.1 Paarbindung und Lochabstand

Bindungsenergie E_b = E(2) + E(0) − 2E(1) (negativ = gebunden), Hüpfenergie ⟨−t Σ c†c⟩ im
2-Loch-Grundzustand, Wahrscheinlichkeit des Lochabstands relativ zu zufällig verteilten Löchern
(Nachbarn = Abstand (1,0), größter Abstand = (2,2)):

| J/t | E_b tJ | E_b frei | Hüpfen tJ | Hüpfen frei | Nachbarn tJ | Nachbarn frei | (2,2) tJ | (2,2) frei |
|---|---|---|---|---|---|---|---|---|
| 0,2 | −0,169 | +0,274 | −5,88 | −7,50 | 0,85 | 0,66 | 0,87 | 1,35 |
| 0,4 | −0,348 | +0,186 | −5,55 | −7,44 | 1,08 | 0,68 | 0,51 | 1,32 |
| 0,6 | −0,527 | +0,104 | −5,28 | −7,40 | 1,29 | 0,71 | 0,34 | 1,29 |
| 1,0 | −0,887 | −0,052 | −4,80 | −7,31 | 1,67 | 0,75 | 0,16 | 1,24 |

- Mit Vorzeichen sind zwei Löcher bei allen J/t gebunden, der Abstand ist kurz. Ohne Vorzeichen sind sie
  bis J/t = 0,6 ungebunden und meiden sich (Nachbarn 0,66–0,75 mal so häufig wie zufällig); bei
  J/t = 1,0 ist E_b nur schwach negativ, die Löcher meiden sich weiterhin.
- Die Vorzeichen kosten Hüpfenergie (−5,5 statt −7,4 bei J/t = 0,4): Weltlinien verschiedenen Vorzeichens
  löschen sich teilweise aus.
- Das passt zu M2: Ein Paar verringert diesen Verlust, ohne Vorzeichen gibt es nichts zu gewinnen.
  Die Rechnung trennt aber nicht die vom Reparatur-Zeitpunkt abhängigen Vorzeichen von den übrigen
  (z. B. beim Austausch der Löcher).

### 6.2 Verschränkung und Bindungen am Loch

S = Verschränkungsentropie zwischen den Lochpositionen und der Spinfolge im „squeezed space“
(Ogata–Shiba: Spins der besetzten Plätze in Platzreihenfolge; tJ in platzweiser Fermion-Ordnung,
frei mit |ψ| als Grundzustand des bosonischen |H|). Bei entarteten Grundzuständen ist der Bereich über
die berechneten Eigenvektoren angegeben. Maximal möglich: ln 120 = 4,79 (2 Löcher), ln 16 = 2,77 (1 Loch).
⟨S_p·S_q⟩ im 1-Loch-Grundzustand für Bindungen mit einem Endpunkt neben dem Loch („am Loch“) und für die
am weitesten entfernten Bindungen („fern“); ohne Loch: −0,351. Für |H| ist der xy-Teil mit dessen
eigenem Bindungsoperator gerechnet (Marshall-gedrehtes S·S).

| J/t | S 2 Löcher tJ | S 2 Löcher frei | S 1 Loch tJ | S 1 Loch frei | S·S am Loch tJ | am Loch frei | fern tJ | fern frei |
|---|---|---|---|---|---|---|---|---|
| 0,2 | 2,93–3,31 | 0,39 | 2,12–2,21 | 0,42 | −0,242 | −0,334 | −0,331 | −0,341 |
| 0,4 | 3,05–3,30 | 0,59 | 2,10–2,23 | 0,55 | −0,270 | −0,341 | −0,339 | −0,344 |
| 0,6 | 3,14–3,23 | 0,72 | 2,12–2,19 | 0,63 | −0,286 | −0,344 | −0,342 | −0,345 |
| 1,0 | 2,86–2,94 | 0,88 | 1,82–2,21 | 0,72 | −0,307 | −0,348 | −0,344 | −0,346 |

- **Loch–Spin-Verschränkung:** Mit Vorzeichen sind Lochpositionen und Spins stark verschränkt
  (S ≈ 3 von maximal 4,8), ohne Vorzeichen kaum (S < 0,9). Die Vorzeichen, also die Interferenz der
  Weltlinien, erzeugen diese Verschränkung (M1, M3). Die Zahl hängt etwas von der Platzreihenfolge im
  squeezed space ab (Tabelle: zeilenweise). `python ed_tj.py --reihenfolge --J 0.4` vergleicht
  zeilenweise, spaltenweise und Schlangenlinie: tJ 2,75–3,30, frei 0,59–0,60. Der Unterschied
  zwischen tJ und frei bleibt also.
- **Bindungen am Loch:** Mit Vorzeichen ist die Singulett-Korrelation der Bindungen neben dem Loch um
  10–30 % geschwächt, ohne Vorzeichen kaum (≤ 3 %). Der xy-Teil bleibt in jeder Bindungsklasse negativ
  (singulettartig). Im Bild von M3: Beiträge verschiedener Reparatur-Zeitpunkte heben sich am Loch
  teilweise auf; eine Bindung, die im Mittel triplettartig wird, gibt es auf 4×4 nicht.

## 7. Kohärenz, Frustration und ODLRO (`kohaerenz.py`)

Frage: Gibt es Bewegungen, bei denen die Reparaturen die Phasen *von selbst* so setzen, dass sich alle
Weltlinien addieren (Voraussetzung für ODLRO, also langreichweitige Paarkohärenz)?

### Was „frustriert“ und „unfrustriert“ heißt

Für eine Bewegungsart (z. B. „Paar macht Raupenschritte“) wird eine effektive Hüpfmatrix M aufgestellt:
M[Q, P] = Summe der Vorzeichen aller Prozesse, die von Lage P nach Lage Q führen. Das Vorzeichen eines
einzelnen Eintrags ist Konvention (es hängt davon ab, welches Vorzeichen man einem Basiszustand gibt).
Konventionsfrei ist nur das Produkt der Vorzeichen entlang einer **geschlossenen Schleife** P → Q → … → P.

- **Unfrustriert:** Jede geschlossene Schleife hat das Produkt +1. Dann kann man die Basiszustände so
  umeichen, dass alle Einträge positiv sind. Alle Wege zwischen zwei Lagen haben dasselbe Vorzeichen und
  addieren sich; der Grundzustand ist knotenfrei (Perron–Frobenius) und gewinnt die volle Hüpfenergie:
  λ_max(M) = λ_max(|M|).
- **Frustriert:** Mindestens eine Schleife hat −1. Keine Umeichung macht alle Einträge positiv, so wie man
  auf einem Dreieck nicht alle drei Bindungen antiferromagnetisch befriedigen kann (daher der Name).
  Wege löschen sich teilweise aus, λ_max(M) < λ_max(|M|); im Extremfall gibt es ein flaches Band, das
  Teilchen ist durch Interferenz lokalisiert.
- Maß: λ_max(M) / λ_max(|M|) (1 = unfrustriert). Für ODLRO in einem Kanal braucht man, dass die
  Bewegung in diesem Kanal unfrustriert ist (Feynman-Argument für Bosonen mit positiven Gewichten).

### 7.1 Regel für geschlossene Prozesse

    Vorzeichen eines geschlossenen Prozesses = (−1)^s · sgn(π_Loch)

s = Zahl der Spin-Vertauschungen (J-Prozesse), π_Loch = Permutation der Löcher. **Lochsprünge tragen kein
Vorzeichen.** Herleitung aus den Schritten 1–2 der Beweis-Skizze im Haupt-README: Das Vorzeichen ist
sgn(π↑)·sgn(π↓); alle m + s Schritte sind Transpositionen, also sgn(π↑)·sgn(π↓)·sgn(π_Loch) = (−1)^(m+s);
bei einem geschlossenen Prozess auf dem bipartiten Gitter ist m gerade. Numerisch: 130 von 130 geschlossenen
Zufallsprozessen auf beliebigen Konfigurationen. Folgerung: Kohärent ist eine Bewegung genau dann, wenn auf
jeder geschlossenen Weltlinie die Zahl der Reparaturen und die Zahl der Lochvertauschungen dieselbe Parität
haben.

### 7.2 Geschlossene Prozesse von Paaren und Löchern

| Prozess | Löcher vertauscht | erreichbare Vorzeichen (Reparaturen beliebig eingestreut) |
|---|---|---|
| Raupe (B folgt A), Umläufe `RDLU`, `RRDLLU`, `RRDDLLUU`, `RRDDLULU`, `RRRDDDLLLUUU`, beide Untergitter | nein | nur +1 (10/10), Néel ohne J wiederhergestellt |
| Raupe mit Wende (Kopf und Schwanz tauschen) | ja | nur **−1** |
| zwei Raupenpaare tauschen die Plätze | ja (zwei Vertauschungen) | nur +1 |
| zwei Raupenpaare, Umläufe kreuzen sich | nein | nur +1 |
| zwei einzelne Löcher tauschen die Plätze | ja | nur −1 |
| zwei einzelne Löcher, Umläufe kreuzen sich | nein | +1 und −1 |
| ein Loch, Umlauf `RDLU` | – | +1 und −1 |

Raupenpaare sehen einander nicht (keine Spur) und sind untereinander bosonisch. Die Wende eines einzelnen
Paares ist aber ein Lochaustausch ohne Reparatur und kostet −1.

### 7.3 Effektive Paarbewegung

Paar aus Nachbarn im Néel-Hintergrund (8×8, alle 112 Bindungen), drei Schrittarten:
**Raupe** (Kopf springt, Schwanz folgt, s = 0), **Pivot** (ein Loch läuft über die Plakettenecke zum anderen
Nachbarn des Partners, danach eine Reparatur, s = 1), **Domino** (beide Löcher seitlich, danach eine
Reparatur, s = 1). Alle Prozesse mit Gewicht 1.

| Schrittarten | λ_max(M) | λ_max(\|M\|) | Verhältnis | Entartung |
|---|---|---|---|---|
| nur Raupe | 2,000 | 5,696 | 0,351 | 49 (flaches Band) |
| nur Pivot + Domino | 5,719 | 5,719 | **1,000** | 1 |

- Die reine Raupe ist stark frustriert, das Paar ist durch Interferenz lokalisiert (Wenden kosten −1).
- Pivot und Domino sind **unfrustriert**. Die Umeichung, die alle Einträge positiv macht, ist genau das
  d-Wellen-Muster: Vorzeichen von Δ(b)|Néel⟩ × (+1 waagerecht, −1 senkrecht), mit
  Δ(i,j) = c_i↑ c_j↓ − c_i↓ c_j↑. Ein Pivot hat gegenüber der Raupe mit Wende zum selben Ziel eine Reparatur
  *und* eine Lochvertauschung mehr; die beiden Minuszeichen heben sich auf (392 von 392 gemeinsamen
  Übergängen mit gleichem Vorzeichen).
- In dieser Eichung sind die Raupenschritte mit Wende alle positiv (392/392), die geraden Raupenschritte
  alle negativ (192/192). Die Frustration kommt also aus den geraden Raupenschritten.
- Varianten mit der Vertauschung *vor* den Sprüngen (auf dem künftigen Weg) haben dasselbe Vorzeichen wie
  die Reparatur danach (7.4); sie ändern nur Gewichte, nicht diese Ergebnisse.

Gemischt, M = Raupe + λ·(Pivot + Domino), oberster Zustand im Überlapp mit dem s- bzw. d-Wellen-Paar
Σ_b f(b) Δ(b)|Néel⟩:

| λ | 0 | 0,25 | 0,5 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|
| λ_max(M) / λ_max(\|M\|) | 0,35 | 0,48 | 0,57 | 0,68 | 0,78 | 0,87 |
| Überlapp d-Welle | 0,02 | 0,77 | 0,78 | 0,78 | 0,78 | 0,78 |
| Überlapp s-Welle | 0 | 0 | 0 | 0 | 0 | 0 |

### 7.4 Ein Loch: Schritt aus 2 Sprüngen und 1 Vertauschung

Alle Prozesse aus zwei Lochsprüngen und genau einer Spin-Vertauschung, die wieder bei Néel enden:

| Ziel | Prozesse | Vertauschung vorher (berührt künftigen Weg) | dazwischen | nachher (Reparatur) | Vorzeichen je Ziel |
|---|---|---|---|---|---|
| diagonal (±1, ±1) | 4 | 2 | 0 | 2 | einheitlich |
| gerade (±2, 0), (0, ±2) | 2 | 1 | 0 | 1 | einheitlich |

- Die Vertauschung *vor* den Sprüngen dreht ein richtiges Paar auf dem künftigen Weg um; das Loch schluckt
  beide falschen Spins und hinterlässt keine Spur. Sie hat immer dasselbe Vorzeichen wie die Reparatur
  danach (beide s = 1): Der Zeitpunkt spielt keine Rolle, die Beiträge addieren sich.
- Es bleibt −1 pro Schritt: Dreieck (diagonal, diagonal, gerade zurück) −1, Raute aus vier Diagonalschritten
  +1. Die Hüpfmatrix auf dem Untergitter hat λ_max(M) = 7,87 gegen λ_max(|M|) = 20,49; das Spektrum ist
  gespiegelt, M ist bis auf Umeichen −|M|.
- Bei gleichen Gewichten folgt daraus das Band E(k) ∝ (cos kₓ + cos k_y)² mit Minimum auf der ganzen Linie
  |kₓ| + |k_y| = π (darauf (π/2, π/2) und (π, 0)), die bekannte Form der Einzelloch-Dispersion im
  t-J-Modell. Das Loch ist kohärent, aber frustriert, und sein tiefster Zustand ist eine entartete Linie.

### 7.5 Deutung

| | Einzelloch | Paar |
|---|---|---|
| keine Spur; Reparatur vorher und nachher konstruktiv | ja | ja |
| Vorzeichen geschlossener Wege | (−1)^Schritte → frustriert | +1 in der d-Wellen-Eichung → unfrustriert (Pivot, Domino) |
| bosonisch untereinander | nein (Löcher tauschen mit −1) | ja (Paartausch +1) |

Die Reparatur im Paarschritt gleicht das Vorzeichen des Lochaustauschs von selbst aus. Dadurch wird aus
zwei fermionischen Löchern ein Objekt, dessen Bewegung keine Vorzeichen hat – in der d-Wellen-Basis. Das ist
der qualitative Mechanismus, der hier für ODLRO im d-Wellen-Paarkanal spricht. Dagegen wirken die geraden
Raupenschritte (reine Lochreparatur ohne J).

Das gilt im strengen Néel-Hintergrund. Lässt man J-Fluktuationen zu, wird auch die Paarbewegung
frustriert (7.6).

Offen: ob der Tausch zweier Paare auch über Pivot- und Domino-Schritte +1 ergibt (für das Vielpaar-Argument
nötig), und wie groß λ physikalisch ist (Energienenner der Zwischenzustände).

### 7.6 Gelockertes Néel: Vertauschungen, die keine Reparaturen sind (`fluktuationen.py`)

Bisher war der Hintergrund zwischen zwei Paarschritten immer exakt Néel. Jetzt darf J auch ein richtiges
Nachbarpaar umklappen (Fluktuation: zwei benachbarte falsche Spins), höchstens eine Fluktuation, höchstens
zwei Plätze vom Paar entfernt. Alle Prozesse aus 1 bis 3 Elementarschritten zwischen solchen Zuständen werden
vollständig aufgezählt (6×6-Fenster, 1456 Zustände, 72 192 Übergänge); keine Schrittart ist vorgegeben.

| Zustandsraum, Prozesse | Verhältnis (Vertauschung : Sprung = 1 : 1) | Verhältnis (0,3 : 1) |
|---|---|---|
| sauberes Néel, alle Prozesse | 0,902 | 0,749 |
| sauberes Néel, nur Prozesse mit Vertauschung | **1,000** | **1,000** |
| mit Fluktuation, nur Prozesse mit Vertauschung | 0,923 | 0,881 |
| mit Fluktuation, alle Prozesse | 0,917 | 0,867 |

- Die vollständige Aufzählung bestätigt 7.3: Im sauberen Néel sind alle Prozesse mit Vertauschung
  unfrustriert. Mit einer Fluktuation sind sie frustriert.
- **Kleinste frustrierte Schleife** (das Paar bleibt stehen): J klappt ein Nachbarpaar neben dem Paar um;
  das Paar dreht diese Fluktuation mit 2 Sprüngen und 1 Vertauschung um 90°; J klappt sie zurück.
  3 Vertauschungen, kein Lochaustausch → −1.
- **Phase-String-Form der Regel:** Bei einem geschlossenen Prozess ändert sich die Zahl der ↑-Spins auf
  einem Untergitter nicht. Jede Vertauschung und jeder Sprung eines ↑-Spins ändert sie um 1, also ist
  (−1)^s = (−1)^(n↑) und Vorzeichen = (−1)^(n↑) · sgn(π_Loch), n↑ = Zahl der Lochsprünge mit einem ↑-Spin
  (129 von 129 geschlossenen Zufallsprozessen). Im sauberen Néel bewegen Pivot und Domino ↑- und ↓-Spins
  so, dass n↑ und der Lochaustausch dieselbe Parität haben. Ordnet das Paar eine Fluktuation um, kann es auf
  einer geschlossenen Weltlinie eine ungerade Zahl von ↑-Spins bewegen, ohne dass Löcher tauschen.
- **Die Frustration ist lokal.** Eine Fluktuation, die die Wege der Löcher nicht berührt, ändert nichts:
  Die Schleife zerfällt dann in eine reine Spinschleife (gerade Zahl von Vertauschungen, Marshall-Regel) und
  eine saubere Paarschleife.
- **Deutung (Hypothese):** Das Paar bleibt kohärent, solange Spinfluktuationen in seiner Nähe selten sind
  oder als feste Singulett-Paare mitlaufen; viele bewegliche Fluktuationen um das Paar herum frustrieren
  es. Das passt zum Bild der Phase-String-Theorie, nach dem Supraleitung kurzreichweitige, gepaarte
  Spinkorrelationen braucht (aus dem Gedächtnis). Geprüft ist hier nur der Fall einer einzelnen Fluktuation.


## 8. Grenzen

- Die Abschnitte 2–4 behandeln Vorzeichen *einzelner* Prozesse. Die Schlüsse auf Verschränkung in
  Abschnitt 5 sind Argumente über die Struktur der Amplitude; quantitativ geprüft ist nur, was in
  Abschnitt 6 steht.
- Die Reparatur tauscht nur benachbarte *falsche* up/dn-Paare. Der J-Term tauscht auch richtige Paare
  (Zusatz im Haupt-README); viele „nicht reparierbar“-Fälle sind ein Artefakt dieser Einschränkung.
- Die ED läuft auf 4×4 mit periodischen Rändern. Dieses Gitter hat zusätzliche Symmetrien und starke
  Endlichkeitseffekte; die Bindungsenergien sind dort nur qualitativ aussagekräftig.
- |H| entfernt *alle* Vorzeichen, auch die beim Austausch zweier Löcher, nicht nur die vom
  Reparatur-Zeitpunkt abhängigen. Der Vergleich zeigt die Wirkung der Vorzeichen insgesamt.
- Abschnitt 7 nutzt effektive Modelle: fester Néel-Hintergrund (in 7.6 mit höchstens einer Fluktuation),
  Paare nur aus Nachbarn, Prozesse bis 3 Elementarschritte, alle Prozesse mit einfachen Gewichten, keine
  Energienenner. Die Vorzeichen sind exakt, die Gewichte nicht.

## 9. Stichworte zur Literatur

Aus dem Gedächtnis, vor dem Zitieren bitte prüfen:

- Phase-String-Effekt (Z.-Y. Weng und Mitarbeiter, ab etwa 1996): das Vorzeichen eines Lochweges zählt
  Austausche mit Spins in der Marshall-Basis; n₋ hier ist eng damit verwandt.
- Mutual-Chern-Simons-Theorie (Kou, Qi, Weng, 2005): Löcher und Spinonen sehen einander als π-Fluss
  (passt zu 4.4: Löcher sehen nur unreparierte Spins).
- σ·t-J-Modell (Weng-Gruppe): t-J ohne Phase-String; dort wurde nach meiner Erinnerung eine stark
  geschwächte Paarbindung gefunden.
- Vorzeichenstruktur des t-J-Modells (K. Wu, Z.-Y. Weng, J. Zaanen, 2008): Vorzeichen einer Weltlinie aus
  Loch-Spin-Austauschen und Loch-Loch-Austauschen; die Regel in 7.1 ist dieselbe Struktur in einer anderen
  Eichung.
- Einzelloch-Dispersion ∝ (cos kₓ + cos k_y)² mit Minimum bei (π/2, π/2) (u. a. Kane, Lee, Read;
  Martinez, Horsch; Dagotto-Übersichtsartikel).
- Trugman-Schleifen (S. A. Trugman, 1988); String-Bild und Paarbildung über geteilte Strings
  (u. a. Grusdt, Bohrdt und Mitarbeiter).

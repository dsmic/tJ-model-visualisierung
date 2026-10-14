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
  (Verhältnis 0,92 statt 1), und zwar dort, wo das Paar eine Fluktuation in seiner Nähe umordnet. Die Raupe
  ist dagegen immun gegen Fluktuationen auf ihrem Weg (der Schwanz bewegt jeden Spin ein zweites Mal), ein
  einzelnes Loch nicht.
- **Abfall mit der Entfernung (7.9):** Die Raupe allein ist unbeweglich (exakt flaches Band). Beweglich wird
  das Paar nur durch J-Wenden (Pivot, Domino), und deren Vorzeichen hängt davon ab, ob die Nachbarbindung,
  deren Spins sie bewegen, antiparallel ist. Bei *eingefrorenen* zufälligen Fluktuationen fällt die
  Paarkohärenz deshalb mit der Entfernung ab (mit Swap-Dynamik nicht, siehe 7.11). Kein Abfall, wenn diese Bindungen in jeder Komponente
  antiparallel sind (z. B. Singuletts, kurzreichweitiges RVB): Dann ist jedes Vorzeichen geometrisch
  festgelegt, Paare sind Bosonen mit endlicher Masse, und ODLRO folgt aus der Bose-Kondensation (Argument,
  nicht gerechnet). Die frühere Deutung „Raupe + gebundene Fluktuation“ (7.7, 7.8) ist dadurch überholt.
- **Idee für ODLRO (7.10):** Ein Singulett ist eine Bindung, die nie repariert werden muss. Tauscht ein Paar mit
  ihm den Platz (zwei Sprünge, kein J), liefert das Singulett von selbst das Vorzeichen, das sonst die
  Reparatur liefert. Im Singulett-Hintergrund bewegt sich das Paar daher ohne Vorzeichen, ohne Reparatur und
  ohne Spur; ein einzelnes Loch dagegen nicht. Das ist das Gegenstück zu Nagaoka (ein Loch ist im Ferromagneten
  vorzeichenfrei) und im Kern Andersons RVB-Bild, hier aus der Vorzeichenregel hergeleitet.
- **Swaps als Filter (7.11):** Vertauschungen tragen kein Vorzeichen, sie erzeugen alternative Wege. Im
  Quanten-Grundzustand wird für jeden Paarschritt über alle Varianten summiert; übrig bleibt der konstruktive
  Überschuss mit festem, überall gleichem Vorzeichen. Die Auslöschung macht das Paar schwerer, aber nicht
  inkohärent über große Entfernungen. Für die Raupe sind alle Varianten konstruktiv, für ein einzelnes Loch
  nicht: Das begünstigt Paare, und Paare sind Bosonen.
- **Getrennte Partner (7.12):** Im Néel ist ein Loch eine ↑-Lücke, das andere eine ↓-Lücke. Bei sofort
  reparierten J-Schritten vertauscht jede nur Elektronen ihrer Sorte; das Vorzeichen zweier unabhängig laufender
  Löcher ist daher das Produkt ihrer Einzelvorzeichen (1563 von 1563), deckungsgleiche Wege geben immer +1. Die
  Amplitude für Trennen und Wiedertreffen ist G(R)² ≥ 0: Die Löcher müssen nicht zusammenbleiben (Cooper-Paarung
  im Ortsraum).
- **Viele Löcher (7.13):** Mit J-Schritten und Raupenschritten gilt für beliebig viele Löcher
  Vorzeichen = (−1)^Schritte · sgn(π der ↑-Lücken) · sgn(π der ↓-Lücken) (1500 von 1500). Die Löcher sind dann zwei
  Sorten gewöhnlicher Fermionen ohne Vorzeichen zwischen den Sorten und ohne Phase-String; der Raupenschritt ist
  korreliertes Hüpfen einer Sorte über die andere. Daraus folgt das BCS-Bild mit überlappenden Paaren, und die
  Lage der Fermi-Fläche lässt nur die d-Welle zu (Argument).
- **Exakte Diagonalisierung (4×4, periodisch):** Mit den Vorzeichen sind zwei Löcher gebunden
  (E_b = −0,35 t bei J/t = 0,4), in der vorzeichenfreien Version |H| nicht (+0,19 t). Die Vorzeichen
  erzeugen außerdem eine starke Verschränkung zwischen Lochpositionen und Spins (Entropie ≈ 3 statt
  < 0,9) und schwächen die Singulett-Bindungen neben dem Loch (Abschnitt 6).

## Inhalt

| Datei | Zweck | Start |
|---|---|---|
| `zwei_loch.py` | Kernfunktionen: Wege, Sprung-Reihenfolgen, Reparatur-Strategien, Suche (baut auf `../tj_gui.py` auf) | – |
| `tj_gui_2loch.py` | GUI: alles aus `tj_gui.py` plus zwei Leisten für zwei Löcher | `python tj_gui_2loch.py` |
| `tj_gui_web.html` | dieselben Funktionen als Webseite (JavaScript, für Touch bedienbar, z.B. auf dem iPad); Vorzeichen gegen `tj_gui.py` geprüft (400/400 Prozesse). Zusätzlich: zwei Bewegungen A und B merken, abspielen und vergleichen (relatives Vorzeichen, Regelprüfung); auf der veröffentlichten Seite speichert „Für Claude speichern“ den Vergleich im Datenspeicher des Artifacts (Sammlung `vergleiche`) | im Browser öffnen |
| `experimente.py` | reproduziert die Tabellen aus Abschnitt 4 (Markdown-Ausgabe) | `python experimente.py` |
| `kohaerenz.py` | Regel für geschlossene Prozesse, Paar- und Einzelloch-Bewegung, Frustration (Abschnitt 7, etwa 20 s) | `python kohaerenz.py` |
| `fluktuationen2.py` | wie `fluktuationen.py`, aber bis zu zwei Fluktuationen, mit Gegenprobe (Abschnitt 7.8, etwa 6 min auf 4 Kernen; numpy, scipy) | `python fluktuationen2.py` |
| `beweglichkeit.py` | J-Wende über intakte und gebrochene Bindungen; Bandstruktur des Paares auf dem 12×12-Torus (Bloch-Blöcke mit Verschiebungsvorzeichen), sauber und mit bis zu zwei Fluktuationen (Abschnitt 7.9, etwa 4 min auf 4 Kernen; numpy, scipy; `--pruefen` vergleicht auf kleinen Tori mit der vollen Prozessmatrix) | `python beweglichkeit.py` |
| `partner.py` | zwei getrennt laufende Löcher: Vorzeichen = Produkt der Einzelvorzeichen (7.12); viele Löcher als zwei Sorten Fermionen (7.13); wenige Sekunden | `python partner.py` |
| `spurpaare.py` | Paare auf gemeinsamer Spur (nur Prozesse mit leerer ungerader Menge, Löcher dürfen sich entlang der Spur trennen): Bandstruktur auf dem Torus (Abschnitt 7.18, unter einer Minute) | `python spurpaare.py` |
| `fluktuationen.py` | Fluktuationen auf dem Weg (ein Loch, Raupe) und Paarbewegung mit einer J-Fluktuation, vollständige Prozessaufzählung (Abschnitt 7.6, etwa 5 min) | `python fluktuationen.py` |
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
  Genauer: Kehrt die Raupe zur selben Bindung zurück, sind die Löcher genau dann vertauscht, wenn die Zahl
  der Schritte ungerade ist (jeder Schritt wechselt das Untergitter beider Löcher). Die Raupe ist damit ein
  Teilchen mit umgekehrtem Hüpfvorzeichen auf dem Liniengraphen des Quadratgitters (Bindungen als Knoten,
  verbunden bei gemeinsamem Platz): gleiches Spektrum wie −A(Liniengraph), geprüft. Der Liniengraph besteht
  aus Dreiecken (drei Bindungen an einem Platz); daher das flache Band mit Vielfachheit
  |Bindungen| − |Plätze| + 1 = 112 − 64 + 1 = 49, wie beim Kagome-Gitter.
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
vollständig aufgezählt (6×6-Fenster, 1456 Zustände, 70 600 Übergänge); keine Schrittart ist vorgegeben.
Vertauschungen werden nur nahe an den Löchern betrachtet; diese Auswahl hängt nur von den Lochplätzen ab, so
dass Hin- und Rückweg gleich behandelt werden (die Prozessmatrix ist symmetrisch, das Skript prüft es).

| Zustandsraum, Prozesse | Verhältnis (Vertauschung : Sprung = 1 : 1) | Verhältnis (0,3 : 1) |
|---|---|---|
| sauberes Néel, alle Prozesse | 0,902 | 0,749 |
| sauberes Néel, nur Prozesse mit Vertauschung | **1,000** | **1,000** |
| mit Fluktuation, nur Prozesse mit Vertauschung | 0,923 | 0,881 |
| mit Fluktuation, alle Prozesse | 0,917 | 0,866 |

- **Korrektur (7.9):** Die Zeile „sauberes Néel, nur Prozesse mit Vertauschung“ ist trivial. Ein Pivot oder
  Domino landet nach seinen zwei Sprüngen in einem Zustand mit einer Fluktuation neben dem Paar und wird in
  dieser Aufzählung dort zerlegt (2 Sprünge, dann 1 Vertauschung). Zwischen sauberen Zuständen bleibt nur
  „umklappen und zurück“ (Diagonalbeiträge), daher 1,000. Dass Pivot und Domino im sauberen Néel unfrustriert
  sind, zeigt 7.3 (auf dem Torus bestätigt in 7.9). Mit einer Fluktuation sind die Prozesse mit Vertauschung
  frustriert.
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

**Fluktuation auf dem Weg: einzelnes Loch gegen Raupe.** Weg UU ab (5,3). Vorher wird ein richtiges
Nachbarpaar umgeklappt (keine Reparatur), danach wird mit beliebigen reinen Spin-Vertauschungen bis Néel
geheilt (Breitensuche); angegeben sind alle dabei erreichbaren Vorzeichen:

| umgeklapptes Paar | ein Loch | Raupenpaar (B folgt von (6,3)) |
|---|---|---|
| keines | −1 | +1 |
| U mit UR (berührt den Weg) | **+1** | +1 |
| U mit UL (berührt den Weg) | **+1** | +1 |
| UU mit UUR (Zielplatz) | **+1** | +1 |
| R mit UR (neben dem Weg) | −1 | +1 |

- **Ein Loch:** Eine Fluktuation, die den Weg berührt, dreht das Vorzeichen um; eine daneben nicht. In der
  Form (−1)^(n↑): Das Loch bewegt jeden Spin auf seinem Weg einmal. Die Fluktuation tauscht aus, welcher
  Spin (↑ oder ↓) dort steht, und ändert damit die Parität von n↑. Reine Spin-Vertauschungen zählen nicht.
- **Raupe:** immer +1. Der Schwanz bewegt jeden Spin, den der Kopf verschoben hat, ein zweites Mal, also ist
  n↑ immer gerade – in jedem Hintergrund. Die Raupe trägt keinen Phase-String; als Vorzeichen bleibt nur der
  Lochaustausch (Wende). Voraussetzung: Zwischen Kopf- und Schwanzsprung bewegt kein J-Prozess den
  verschobenen Spin.
- **Pivot und Domino:** Die beiden Partner bewegen verschiedene Spins je einmal. Der Lochaustausch ist
  ausgeglichen (7.3), aber das Vorzeichen hängt vom Hintergrund ab; daher die Frustration in der Tabelle oben.

| Paarschritt | geschützt gegen Wende (Lochaustausch) | geschützt gegen Fluktuationen im Hintergrund |
|---|---|---|
| Raupe | nein | ja |
| Pivot, Domino | ja (d-Wellen-Eichung) | nein |
- **Deutung (Hypothese, siehe aber 7.7):** Das Paar bleibt kohärent, solange Spinfluktuationen in seiner Nähe selten sind
  oder als feste Singulett-Paare mitlaufen, oder solange es sich überwiegend als Raupe bewegt; viele
  bewegliche Fluktuationen um ein Paar, das sich über Pivot und Domino bewegt, frustrieren es. Das passt zum Bild der Phase-String-Theorie, nach dem Supraleitung kurzreichweitige, gepaarte
  Spinkorrelationen braucht (aus dem Gedächtnis). Geprüft ist hier nur der Fall einer einzelnen Fluktuation.


### 7.7 Mechanismus bei aufgeweichtem Néel: Raupenpaar mit gebundener Fluktuation

Ausgangspunkt: Schon bei einem einzelnen Loch ändern Vertauschungen, die den Weg berühren, das Vorzeichen
(7.6, und in der Web-GUI gespeicherte Vergleiche: Weg UU bzw. DD mit vorherigem Umklappen und Zurückklappen,
jeweils entgegengesetztes Vorzeichen). Bei aufgeweichtem Néel sind für ein einzelnes Loch also praktisch
immer beide Vorzeichen erreichbar. Gesucht ist ein Mechanismus, der trotzdem Paarkohärenz (ODLRO) erlaubt.

1. **Die Raupe ist immun gegen freie Fluktuationen.** Der Schwanz bewegt jeden Spin, den der Kopf verschoben
   hat, ein zweites Mal; n↑ ist immer gerade, die Phase-Strings der beiden Löcher heben sich auf. Fluktuationen
   auf dem Weg ändern das Vorzeichen der Raupe nicht (Tabelle in 7.6). Ein Raupenpaar hat also, anders als ein
   einzelnes Loch, in jedem Hintergrund ein festes Vorzeichen.
2. **Ihr einziges Vorzeichen ist der Lochaustausch**, und der macht die reine Raupe frustriert (Liniengraph,
   flaches Band, 7.3).
3. **Freie J-Fluktuationen in der Nähe des Paares verringern die Frustration deutlich.** Vollständige
   Aufzählung wie in 7.6 (6×6-Fenster, höchstens eine Fluktuation), Prozesse ohne Vertauschung (Raupe, auch
   durch eine Fluktuation hindurch oder sie schluckend) und Prozesse ohne Lochsprung (Fluktuation entsteht,
   wandert, vergeht). **Gegenprobe:** dieselben Prozesse, aber jeder Übergang mit Lochsprüngen bekommt ein
   zufälliges Vorzeichen. Wird das Verhältnis dadurch kaum schlechter, kommt die Verbesserung nur von der
   Verdünnung durch die (vorzeichenfreien) reinen Spinprozesse.

| Gewicht Vertauschung | Raupe, sauberes Néel | Raupe + freie J-Fluktuationen | Gegenprobe: zufällige Vorzeichen | Anteil des Grundzustands mit Fluktuation |
|---|---|---|---|---|
| 1 | 0,708 | 0,992 | 0,976 | 71 % |
| 0,3 | 0,708 | **0,898** | 0,761 | 66 % |
| 0,1 | 0,708 | **0,881** | 0,698 | 54 % |

   Bei Gewicht 1 ist der hohe Wert also überwiegend Verdünnung; bei 0,3 und 0,1 (physikalisch liegt das Gewicht
   eines J-Prozesses relativ zu einem Sprung eher bei J/2t ≈ 0,2) sind die echten Vorzeichen deutlich besser als
   zufällige. In diesem Fenster zählen auch Hin-und-zurück-Sprünge als Diagonalbeiträge; deshalb ist die
   saubere Raupe hier 0,708 statt 0,351 wie in 7.3. Mit allen Prozessen, also auch Pivot und Domino, ist es
   etwas frustrierter (0,917 bzw. 0,866, Tabelle in 7.6).
   (Eine frühere Fassung nannte hier 0,993 / 0,901 / 0,882 aus einer Aufzählung, deren Prozessmatrix nicht ganz
   symmetrisch war; die Abweichung beträgt höchstens 0,003. Die Gegenprobe ist neu.)

**Korrektur (7.9):** Die Auswahl „Prozesse ohne Vertauschung oder ohne Lochsprung“ enthält Pivot und Domino,
zerlegt in zwei Hälften über den Zustand mit einer Fluktuation neben dem Paar. „Raupe + freie J-Fluktuationen“
heißt also eigentlich „Raupe + J-Wenden + Fluktuationen“, und die „gebundene Fluktuation“ ist zum großen Teil der
Zwischenzustand einer J-Wende. Nimmt man die Wenden heraus, ist das Paar fast unbeweglich. Außerdem misst das
Verhältnis nur die lokale Frustration, nicht den Abfall mit der Entfernung. Die folgende Hypothese ist deshalb
durch 7.9 ersetzt.

**Mechanismus (Hypothese, überholt durch 7.9):** Das kohärente Objekt ist ein Raupenpaar, das eine
Spinfluktuation an sich bindet. Freie Fluktuationen auf seinem Weg ändern sein Vorzeichen nicht (1.), und die
gebundene Fluktuation liefert beim Wenden das Vorzeichen, das den Lochaustausch ausgleicht (3.). Was beim
einzelnen Loch die Vorzeichen zufällig macht, wird beim Paar zur Ressource. Physikalisch ist das ein
Lochpaar mit Spinwolke (einem umgeklappten Paar, also einer Singulett- bzw. Triplett-artigen Anregung), passend
zum Bild des Cooper-Paares im RVB-Bild.

Offen bzw. nicht geprüft: welche Paarsymmetrie dieser Zustand hat (Überlapp mit der s-Welle 0, mit einer
gleichförmigen d-Welle nur 0,05 bis 0,37), wie es mit realistischen Energien aussieht, und ob der Tausch zweier
solcher Paare +1 ergibt. Zwei Fluktuationen: 7.8.

### 7.8 Zwei Fluktuationen (`fluktuationen2.py`)

Wie 7.7, aber bis zu zwei Fluktuationen in der Nähe des Paares: 15 218 Zustände (60 ohne, 1396 mit einer,
13 762 mit zwei Fluktuationen), 2,58 Millionen Prozesse, Prozessmatrix symmetrisch (geprüft). Innerhalb dieses
Laufs ist „höchstens 1“ etwas anders abgeschnitten als in 7.7 (Prozesse über Zustände mit zwei Fluktuationen
fehlen dort), daher die leicht anderen Werte.

| Gewicht Vertauschung | Raupe, sauber | Raupe + Fluktuationen, höchstens 1 | höchstens 2 | Gegenprobe, höchstens 2, zufällige Vorzeichen | alle Prozesse, höchstens 1 | alle Prozesse, höchstens 2 |
|---|---|---|---|---|---|---|
| 1 | 0,708 | 0,991 | 0,997 | 0,993 | 0,889 | 0,926 |
| 0,3 | 0,708 | 0,880 | **0,932** | 0,721 | 0,813 | 0,887 |
| 0,1 | 0,708 | 0,878 | **0,903** | 0,725 | 0,827 | 0,881 |

Anteile des Grundzustands (Raupe + Fluktuationen, höchstens 2) mit 0 / 1 / 2 Fluktuationen: 7 / 34 / 59 %
(Gewicht 1), 9 / 36 / 54 % (0,3), 25 / 52 / 23 % (0,1).

- **Mehr Aufweichung, weniger Frustration:** Mit einer zweiten erlaubten Fluktuation steigt das Verhältnis in
  allen Fällen weiter (Raupe + Fluktuationen bei 0,3: 0,880 → 0,932; alle Prozesse: 0,813 → 0,887).
- **Gegenprobe:** Mit zufälligen Vorzeichen der Sprungprozesse bleibt es bei 0,3 und 0,1 auf dem Niveau der
  sauberen Raupe (0,72). Die Verbesserung kommt also aus der Vorzeichenstruktur, nicht aus der Verdünnung.
- **Der Grundzustand trägt die Fluktuationen:** Bei 0,3 hat er zu 54 % zwei Fluktuationen in der Nähe des
  Paares. Das Paar ist von einer Spinwolke umgeben, und diese Wolke macht seine Bewegung kohärenter.

Diese Zahlen gelten für Fluktuationen in der Nähe des Paares, die mit ihm wandern. Sie sagen nichts darüber, ob
die Kohärenz mit der Entfernung abfällt, und die Auswahl enthält J-Wenden (Korrektur in 7.7). Die frühere
Schlussfolgerung, ein aufgeweichtes Néel mache die Paarbewegung insgesamt kohärenter, ist durch 7.9 ersetzt: Für
freie Fluktuationen im Hintergrund ist ein Abfall mit der Entfernung zu erwarten. Grenzen: kleines Fenster,
Fluktuationen höchstens zwei Plätze vom Paar entfernt, einfache Gewichte statt Energien.

### 7.9 Fällt die Paarkohärenz mit der Entfernung ab? (`beweglichkeit.py`)

ODLRO heißt: G(R) = ⟨Δ†(R) Δ(0)⟩ bleibt für große Abstände R endlich. G(R) summiert über alle Prozesse, die
ein Paar von 0 nach R bringen. Das Verhältnis λ_max(M)/λ_max(|M|) aus 7.3–7.8 ist lokal und sagt darüber nichts.
Beispiel unten: Raupe mit Fluktuationen, aber ohne J-Wenden, hat das Verhältnis 0,92 und ist trotzdem fast
unbeweglich.

**Exakte Grundlage.** Aus 7.1 folgt für zwei beliebige Prozesse mit gleichem Anfang und gleichem Ende:

    relatives Vorzeichen = (−1)^(n↑ des einen + n↑ des anderen) · sgn(relative Lochpermutation)

n↑ = Zahl der Lochsprünge, die einen ↑-Spin bewegen; Spin-Vertauschungen zählen in dieser Form nicht. Ein Abfall
mit R entsteht, wenn die Wege nach R sich in der Parität von n↑ unterscheiden, je nachdem, was unterwegs liegt.

1. **Einzelnes Loch:** n↑ = Zahl der ↑-Spins, die es verschoben hat. Jede Fluktuation am Weg dreht das
   Vorzeichen (7.6). Abfall mit der Entfernung (Phase-String).
2. **Raupe:** n↑ ist immer gerade; aus dem Hintergrund kommt kein Abfall. Aber die Raupe allein ist unbeweglich:
   Auf dem Torus ist ihr oberstes Band exakt flach (alle 72 Impulse am Maximum). Es gibt dann räumlich begrenzte
   Eigenzustände; das Paar kommt durch Interferenz nicht vom Fleck, auch so fällt G(R) ab.
   Fluktuationen ohne J-Wenden helfen kaum (Beweglichkeit 6–9 % des vorzeichenfreien Werts, Tabelle unten).
3. **Beweglich wird das Paar durch J-Wenden** (Pivot, Domino). Im sauberen Néel ist das Band dann gekrümmt, mit
   eindeutigem Maximum bei K = 0. Die Krümmung ist praktisch dieselbe wie mit Pivot und Domino allein (2,887
   gegen 2,881 bei Gewicht 1); die Raupe trägt zur Beweglichkeit nichts bei.
4. **Vorzeichen einer J-Wende.** Pivot und Domino bewegen die beiden Spins einer Nachbarbindung je einmal. Ihr
   Beitrag zu n↑ ist die Zahl der ↑-Spins auf dieser Bindung: bei antiparallelen Spins fest, bei parallelen
   umgedreht. Paar a = (3,3), b = (4,3); Pivot: a läuft über c = (3,4) nach a′ = (4,4); Domino: beide Löcher auf
   c, a′. Vergleich mit Raupenwegen zum selben Ziel (jeder Spin zweimal bewegt), danach mit reinen
   Vertauschungen bis Néel geheilt (alle Heilwege geben dasselbe Vorzeichen):

| Fluktuation vorher | Pivot gegen Raupe mit Wende | Domino gegen zwei Raupenschritte |
|---|---|---|
| keine | +1 | +1 |
| genau die bewegte Bindung c–a′ umgeklappt | +1 | +1 |
| c mit seinem äußeren Nachbarn (bewegte Bindung gebrochen) | **−1** | **−1** |
| a′ mit seinem oberen Nachbarn (bewegte Bindung gebrochen) | **−1** | **−1** |
| Bindung auf der anderen Seite des Paares (nicht berührt) | +1 | +1 |
| Bindung zwei Plätze entfernt (nicht berührt) | +1 | +1 |

   Eine umgeklappte Bindung ist weiterhin antiparallel und stört nicht. Es kommt nur darauf an, ob die bewegte
   Bindung gebrochen ist, also parallele Spins trägt.

**Folgerung.**

- **Aufgeweichtes Néel mit zufälligen, eingefrorenen Fluktuationen** (dieser Punkt setzt voraus, dass jeder
  Weg einen festen Hintergrund sieht; mit Swap-Dynamik gilt er nicht, siehe 7.11): Ein Teil der J-Wenden trifft gebrochene
  Bindungen. Welche das sind, hängt vom Weg und vom Zeitpunkt ab. Je weiter das Paar läuft, desto mehr Wenden
  braucht es, desto stärker heben sich die Beiträge auf. G(R) fällt etwa exponentiell ab, mit einer Länge von der
  Größenordnung des Abstands zwischen gebrochenen Bindungen am Weg (abgeschätzt, nicht gerechnet). Die
  Formulierung aus 7.7/7.8 führt also zu einem Abfall; das Argument „bei aufgeweichtem Néel sind beide Vorzeichen
  erreichbar“ gilt auch für die Wenden des Paares.
- **Kein Abfall**, wenn jede Bindung, über die das Paar wendet, in jeder Komponente der Wellenfunktion
  antiparallel ist (S^z = 0 auf der Bindung, z. B. ein Singulett (↑↓ − ↓↑)/√2). Dann gibt die Wende in jeder
  Komponente dasselbe Vorzeichen. Das ist das Bild eines kurzreichweitigen RVB-Hintergrunds: Das Paar wendet,
  indem es mit einem Singulett die Plätze tauscht. Nach einer Wende liegt die bewegte Bindung wieder neben dem
  Paar (Pivot dreht Paar und Bindung im 2×2-Quadrat, Domino tauscht sie). Für eine Verschiebung braucht das Paar
  aber immer wieder Bindungen aus dem Hintergrund.

**Begründungskette für ODLRO** (unter dieser Bedingung):

1. Das Vorzeichen jedes Prozesses hängt nur von der Geometrie ab, nicht vom Hintergrund und nicht von R:
   Raupenschritte haben gerades n↑, Wenden über antiparallele Bindungen einen festen Beitrag, Vertauschungen
   zählen nicht (Regel oben, Tabelle).
2. Paare sind untereinander bosonisch: Platztausch zweier Paare = zwei Lochvertauschungen = +1 (7.2).
3. Ein Paar ist beweglich: eindeutiges Bandmaximum bei K = 0, endliche Masse (Tabelle unten). Die verbleibende
   Frustration zwischen geraden Raupenschritten und Wenden ist überall dieselbe. Sie macht das Paar schwerer,
   erzeugt aber keinen Abfall mit R.
4. Ein dünnes Gas von Bosonen mit endlicher Masse und ohne Vorzeichen untereinander kondensiert bei T = 0 (in
   2D bei T = 0 echte langreichweitige Ordnung), also G(R) → konstant ≠ 0. Paarsymmetrie aus der Eichung:
   d-Welle (7.3). Dieser Schritt ist das Standardargument für Bosonen und hier nicht gerechnet.

**Bandstruktur auf dem Torus.** 12×12, periodisch. Verschiebungen um (dx, dy) mit dx + dy gerade erhalten das
Néel-Muster. Mit ihnen (samt Vorzeichen aus der Umordnung der fermionischen Erzeuger) zerfällt die
Prozessmatrix in Blöcke M(K), 72 verschiedene Impulse. Prozesse wie in 7.8 (1–3 Elementarschritte zwischen
Zuständen), aber ohne Fenster. Geprüft auf 6×6 und 8×8: Die Spektren aller Blöcke zusammen stimmen mit dem
Spektrum der vollen Prozessmatrix überein (Abweichung ≤ 6·10⁻¹³, `--pruefen`). Die Zahlen für das saubere Néel
stammen aus einer eigenen Aufzählung nur mit sauberen Zuständen, damit Pivot und Domino ganz enthalten sind.
**Beweglichkeit** = Krümmung des obersten Bandes am Maximum geteilt durch dieselbe Größe ohne Vorzeichen
(= Masse ohne Vorzeichen / Masse mit Vorzeichen; 0 = flaches Band).

| Prozesse | Impulse am Maximum (von 72) | Beweglichkeit bei Gewicht 1 / 0,3 / 0,1 | Verhältnis bei 0,3 |
|---|---|---|---|
| sauber, nur Raupe | 72 (flach) | 0 / 0 / 0 | 0,667 |
| sauber, nur Pivot + Domino | 1 | 1,00 / 1,00 / 1,00 | 1,000 |
| sauber, alle (Raupe, Pivot, Domino) | 1 | 0,74 / 0,46 / 0,23 | 0,828 |
| Raupe + bis zu 2 Fluktuationen, ohne J-Wenden | 1 | 0,06 / 0,09 / 0,06 | 0,924 |
| Auswahl aus 7.7/7.8 (nur Sprünge oder nur Vertauschungen), bis zu 2 | 1 | 0,38 / 0,46 / 0,39 | 0,932 |
| alle Prozesse, bis zu 1 Fluktuation | 1 | 0,20 / 0,25 / 0,27 | 0,785 |
| alle Prozesse, bis zu 2 Fluktuationen | 1 | 0,45 / 0,49 / 0,46 | 0,867 |

„Ohne J-Wenden“ heißt: Lochsprünge nur, wenn sich die Zahl der Fluktuationen nicht ändert; damit fallen die
Hälften von Pivot und Domino weg. Die Fluktuationen sind hier an das Paar gebunden (höchstens zwei Plätze
entfernt); freie Fluktuationen im übrigen Hintergrund enthält die Rechnung nicht. Sie zeigt also, dass ein
Paar mit Spinwolke in einem sonst sauberen Néel beweglich ist (endliche Masse, eindeutiges Maximum bei K = 0),
nicht aber, wie schnell die Kohärenz in einem fluktuierenden Hintergrund abfällt.

**Stand der Hypothese:** Die Raupe sorgt für ein Vorzeichen, das nicht vom Hintergrund abhängt, die J-Wenden
für Beweglichkeit. Die Wenden sind nur dann vom Hintergrund unabhängig, wenn die Bindungen, über die sie laufen,
antiparallel sind. ODLRO verlangt in diesem Bild also kurzreichweitig gepaarte Spins (Singuletts) dort, wo das
Paar wendet. Ein Néel mit zufälligen Fluktuationen reicht nicht. Offen: (a) den Abfall in einem fluktuierenden
Hintergrund direkt rechnen (Paar plus eigene Bindung in einem Hintergrund mit zufälligen Fluktuationen);
(b) ob die Dynamik selbst dafür sorgt, dass neben einem Paar bevorzugt Singuletts liegen (Energie);
(c) Quantenhintergrund statt einzelner S^z-Konfigurationen.

### 7.10 Idee: Warum ODLRO herauskommen sollte

Simulation allein liefert kein ODLRO; man braucht einen Grund. Der folgende Gedankengang verbindet die Ergebnisse
aus 7.1–7.9.

**1. Ein Singulett ist eine dauerhaft reparierte Bindung.** Eine J-Wende (Domino, Pivot) bewegt die zwei Spins
einer Nachbarbindung je einmal. Im Néel landen sie mit falscher Orientierung und brauchen eine Reparatur; deren
Vertauschung gleicht den Lochaustausch aus (7.3). In der Form (−1)^(n↑) · sgn(π_Loch) übernimmt diesen Ausgleich
der eine ↑-Spin auf der Bindung. Die Tabelle in 7.9 zeigt: Die Bindung in Néel-Orientierung und dieselbe Bindung
umgeklappt geben **dasselbe** Vorzeichen; nur muss im einen Fall danach repariert werden, im anderen nicht. Das
Umklappen vorher ist genau die „Vertauschung, die den künftigen Weg berührt“ aus 7.4, die dasselbe Vorzeichen
hat wie die Reparatur danach. Ein Singulett (↑↓ − ↓↑)/√2 ist die Überlagerung beider Fälle mit festem
Vorzeichen. Es hat keine falsche Orientierung, wird beim Wenden einfach mitgenommen (zwei Sprünge, kein J; die
Sprünge sind spinrotationsinvariant, das Singulett bleibt ein Singulett) und liefert in jeder Komponente genau
einen ↑-Spin, also den Ausgleich. **Die Reparatur passiert von selbst.**

**2. Im Singulett-Hintergrund bewegt sich das Paar vorzeichenfrei.** Wenn das Paar nur mit Singuletts den Platz
tauscht, haben alle seine Prozesse die Vorzeichenstruktur von Pivot und Domino im sauberen Néel: unfrustriert in
der d-Wellen-Eichung (7.3; Torus 7.9: Verhältnis 1,000, eindeutiges Maximum bei K = 0). Es hinterlässt keine
Spur: Das Singulett liegt danach dort, wo das Paar war, der Hintergrund ist wieder eine Singulett-Bedeckung. Die
Raupe dagegen verschiebt einen Spin um zwei Plätze und reißt dabei ein Singulett auseinander; sie gehört im
Singulett-Hintergrund nicht zur niederenergetischen Bewegung.

**3. Ein einzelnes Loch kann das nicht.** Es bewegt von einem Singulett nur einen Spin. In den beiden Komponenten
ist dieser einmal ↑, einmal ↓, das Vorzeichen also zufällig (Phase-String), und das Singulett bleibt zerrissen
zurück. Löcher gewinnen im Singulett-Hintergrund kohärente Bewegungsenergie nur als Paar: ein kinetischer Grund
für Paarbildung.

**4. Gegenstück zu Nagaoka.** Aus derselben Regel: Ein einzelnes Loch ist genau dann vorzeichenfrei, wenn alle
Spins, die es bewegt, gleich orientiert sind, also im Ferromagneten (Satz von Nagaoka). Ein Paar ist
vorzeichenfrei, wenn die Spins, die es einmal bewegt, in antiparallelen Bindungen stecken, also im
Singulett-Hintergrund. So wie ein Loch sich ein ferromagnetisches Umfeld schafft, sollten Lochpaare sich ein
Singulett-Umfeld schaffen. Ein Singulett-Hintergrund kostet gegenüber dem Antiferromagneten wenig
Austauschenergie (resonierende Singulett-Zustände liegen energetisch nahe am Néel-Grundzustand; aus dem
Gedächtnis), und die Paare gewinnen dort Bewegungsenergie in der Ordnung t²/ΔE statt t²J/ΔE² mit Reparatur.

**5. Von da zu ODLRO.** Im Singulett-Hintergrund sind die Paare (a) vorzeichenfrei beweglich (2.), (b)
untereinander Bosonen (Platztausch +1, 7.2) und (c) unabhängig vom Hintergrund, weil die Bindungen, die sie
bewegen, in jeder Komponente antiparallel sind. Das effektive Modell sind harte Bosonen mit positivem Hüpfen
(d-Wellen-Eichung). Dann sind nach Perron–Frobenius alle Amplituden des Grundzustands positiv, und G(R) ist eine
Summe positiver Terme: Auslöschung ist unmöglich. G(R) kann nur noch durch Lokalisierung abfallen (Kristall,
Valenzbindungskristall); in einer resonierenden Singulett-Flüssigkeit kondensieren die Bosonen, G(R) → konstant.
Die Eichung macht daraus ein d-Wellen-Kondensat.

**Grenzfall, keine Bedingung (nach 7.11):** In einem Quanten-Antiferromagneten trägt eine Nachbarbindung zu
gut einem Viertel parallele Spins. Das macht die Wenden nicht inkohärent, sondern nur schwächer: Nach 7.11 bleibt
der Überschuss der antiparallelen Komponenten, und der ist positiv. Das Singulett ist der Grenzfall mit
größtem Überschuss (1); ein Singulett-Hintergrund ist günstig, aber nicht nötig.

**Was davon geprüft ist und was nicht.** Geprüft: Vorzeichenregel (7.1), Pivot/Domino unfrustriert in der
d-Wellen-Eichung (7.3, 7.9), beide Orientierungen der bewegten Bindung geben dasselbe Vorzeichen (7.9), Paartausch
+1 (7.2), Raupe flach (7.9). Argument, nicht gerechnet: dass der Hintergrund sich tatsächlich in Singuletts
ordnet (Energiebilanz), dass das effektive Bosonenmodell in einer nicht-orthogonalen Singulett-Basis wirklich
vorzeichenfrei bleibt, und der Schritt vom vorzeichenfreien Bosonenmodell zur Kondensation. Neu ist das Bild
nicht: Es ist im Kern Andersons RVB-Idee (Singulett-Flüssigkeit, Dotierung → Supraleitung) bzw. das dotierte
Dimermodell von Rokhsar und Kivelson (aus dem Gedächtnis). Neu ist hier nur, dass es aus der Vorzeichenregel
und dem Reparatur-Bild folgt.

### 7.11 Idee: Die Swaps summieren über Hintergründe, nur der konstruktive Überschuss bleibt

Vorschlag aus der Diskussion. Er ersetzt die Folgerung „Abfall bei aufgeweichtem Néel“ aus 7.9 und macht 7.10
zum Grenzfall. Schwachpunkt: Punkt 5 setzt ein fest gebundenes Paar als Quasiteilchen voraus. 7.12 kommt ohne
das aus.

1. **Swaps machen Wege, keine Vorzeichen.** In der Form (−1)^(n↑) · sgn(π_Loch) tragen Vertauschungen kein
   Vorzeichen. Sie erzeugen Alternativen: Sie ändern, welcher Spin dort sitzt, wo ein Loch vorbeikommt, und damit
   n↑. Zu jedem Lochweg gibt es viele Varianten mit Vertauschungen vorher, dazwischen, nachher, daneben.
2. **Eingefroren gegen dynamisch.** 7.9 hat Fluktuationen als eingefroren behandelt: Jeder Weg sieht einen festen
   Hintergrund, sein Vorzeichen ist zufällig, Wege verschiedener Geometrie heben sich auf, und G(R) fällt ab wie
   bei einem Teilchen im Zufallsfluss. Im Quanten-Grundzustand ist nichts eingefroren. Für jeden Elementarschritt
   des Paares werden alle Varianten kohärent summiert. Die destruktiven ziehen von den konstruktiven ab; übrig
   bleibt eine Nettoamplitude mit festem Vorzeichen. Weil der Hintergrund überall derselbe ist, gilt sie für jeden
   gleichartigen Schritt. Über große Entfernungen bleibt ein Hüpfproblem mit festen Vorzeichen (das Paar mit
   seiner Spinwolke als Quasiteilchen, 7.9: eindeutiges Bandmaximum bei K = 0). Die Auslöschung macht das Paar
   schwerer, erzeugt aber keinen Abfall mit R.
3. **Größe des Überschusses.**
   - Raupe: alle Varianten konstruktiv (7.6), Faktor 1.
   - J-Wende: Anteil der antiparallelen minus Anteil der parallelen Komponenten der bewegten Bindung,
     −4⟨S^z S^z⟩ = −(4/3)⟨S·S⟩. Das ergibt 0,47 im 4×4-Grundzustand ohne Loch (⟨S·S⟩ = −0,351) und 0,36 an
     einem Loch (−0,270 bei J/t = 0,4, beides 6.2). Der Überschuss ist positiv, solange die Bindungen
     antiferromagnetisch korreliert sind.
   - Singulett: 1, der Grenzfall aus 7.10.
4. **Warum Paare.** Ein einzelnes Loch bewegt jeden Spin einmal. Varianten, die seinen Weg berühren, sind
   destruktiv (Beispiel UU mit vorherigem Umklappen, 7.6). Übrig bleiben J-Schritte (7.4) mit einem Überschuss
   wie bei der Wende, also Bewegung nur mit J. Das Paar behält mit der Raupe Bewegung ohne J, bei der alle
   Varianten konstruktiv sind. Das ist ein kinetischer Grund für die Bindung (passt zu 6.1: mit Vorzeichen
   gebunden, ohne nicht).
5. **ODLRO.** Nettoamplituden mit festem, überall gleichem Vorzeichen machen das Paar zu einem Quasiteilchen
   mit Band. Paare sind Bosonen (Platztausch +1, 7.2). Ein dünnes Bose-Gas kondensiert, also G(R) → konstant.

**Wann trotzdem Abfall:** eingefrorene Unordnung (Störstellen, eingefrorene Spins); verschwindender
Überschuss (Bindungen ohne antiferromagnetische Korrelation, z. B. bei hoher Temperatur); Lokalisierung
(flaches Band, Kristall); eine Spinwolke, die nicht gebunden bleibt (Quasiteilchengewicht null).

**Offen:** ob die Spinwolke des Paares endlich bleibt (beim einzelnen Loch ist das als Spin-Polaron bekannt, aus
dem Gedächtnis); die Wechselwirkung zwischen Paaren; der Vielteilchenschritt zur Kondensation.

### 7.12 Idee: Getrennte Partner, das Vorzeichen als Produkt (`partner.py`)

Einwand gegen 7.11: Zwei Löcher bleiben nicht wirklich zusammen. Der Mechanismus muss erlauben, dass sich Paare
trennen und in Phase wieder verbinden.

1. **Die beiden Löcher eines Paares sind verschiedene Sorten.** Im Néel sitzt ein Loch auf A (dort fehlt ein ↑),
   das andere auf B (dort fehlt ein ↓). Ein J-Schritt (zwei Sprünge, sofort repariert) tauscht die ↑-Lücke mit
   genau einem ↑-Elektron. Das ↓-Elektron, das sie unterwegs verschoben hat, bringt die Reparatur an seinen Platz
   zurück. Eine ↑-Lücke vertauscht also nur ↑-Elektronen, eine ↓-Lücke nur ↓-Elektronen.
2. **Darum zerfällt das Vorzeichen in ein Produkt.** Das Vorzeichen eines geschlossenen Prozesses ist
   sgn(π↑) · sgn(π↓) (Herleitung in 7.1). Jeder Faktor gehört zu einer Lücke und hängt nur von deren eigenem Weg
   ab, nicht davon, wo die andere ist und wann sie läuft. Geprüft an 1563 geschlossenen Zufallsprozessen (beide
   Lücken laufen zufällige Schleifen aus 3–6 J-Schritten, Schritte zufällig verschränkt, Anfangsabstand 1 bis 5,
   Wege dürfen sich kreuzen): immer Vorzeichen = s(Weg der ↑-Lücke) · s(Weg der ↓-Lücke). Die Einzelschleifen
   sind frustriert (764-mal +1, 799-mal −1).
3. **Gleiche Wege heben sich auf.** Die ↓-Lücke hat dieselbe Vorzeichenfunktion wie die ↑-Lücke (Symmetrie:
   alle Spins umklappen und um einen Platz verschieben). Laufen beide Lücken deckungsgleiche Wege, ist das
   Produkt +1, wie weit sie auch voneinander entfernt sind und in welcher Reihenfolge sie laufen (729 von 729).
4. **Auftrennen und in Phase verbinden.** Ein Paar trennt sich, beide Löcher laufen unabhängig, und sie treffen
   sich um R verschoben wieder. Die Amplitude dafür ist das Produkt der Einlochamplituden,
   G_↑(R) · G_↓(R) = G(R)², bis auf einen Faktor, der nur von Anfangs- und Endbindung abhängt (die Form des
   Paares). Ein Quadrat ist nie negativ, so stark G(R) selbst durch die Frustration im Vorzeichen schwankt.
   Gemischte Beiträge (einer +, einer −) gibt es, sie können das Quadrat aber nicht umdrehen; übrig bleibt der
   konstruktive Teil. Die Löcher müssen dafür nicht zusammenbleiben.
5. **Das ist Cooper-Paarung im Ortsraum.** Gepaart werden Partner mit derselben Dynamik (↑- und ↓-Lücke sind
   Zeitumkehr-Partner); ihre Phasen heben sich auf, wie im Anderson-Theorem. Tauschen zwei Paare ihre Partner, gibt
   das das fermionische −1 (7.2: zwei einzelne Löcher tauschen mit −1). Das ist die Determinanten-Struktur eines
   BCS-Zustands und verträgt sich mit ODLRO. Kommt eine Anziehung zwischen ↑- und ↓-Lücke hinzu (6.1: mit
   Vorzeichen gebunden, ohne nicht), folgt nach Cooper die Paarbildung. Weil der Paarkanal nie destruktiv ist,
   reicht jede Anziehung. Das ist das Standardargument und hier nicht gerechnet.
6. **Wo der Reparatur-Zeitpunkt eingreift.** Das Produkt ist für sofort reparierte Schritte geprüft. Wird eine
   Reparatur aufgeschoben, bis die andere Lücke die Spur kreuzt, ändert sich das Vorzeichen (4.1). Dann bewegt
   eine Lücke Elektronen der anderen Sorte, und die Faktoren sind nicht mehr unabhängig. An dieser Stelle
   verschränkt der Reparatur-Zeitpunkt die beiden Löcher (Abschnitt 5).
7. **Was Paare bricht.** Alles, was den beiden Partnern verschiedene Vorzeichen gibt, z. B. eine eingefrorene
   Fluktuation, die nur den Weg des einen Partners berührt. Das wäre das Gegenstück zu magnetischen Störstellen
   bei Abrikosov–Gor'kov. Nicht geprüft.

### 7.13 Viele Löcher: zwei Sorten Fermionen (`partner.py`)

7.12 behandelt zwei Löcher; ODLRO ist eine Vielteilchenfrage.

1. **Test.** Drei ↑-Lücken (auf A) und drei ↓-Lücken (auf B) auf dem 10×10-Torus, teils als Nachbarpaare, teils
   einzeln. Zufällige Folgen aus J-Schritten und Raupenschritten, danach zurück auf dieselben Plätze; Lücken
   derselben Sorte dürfen dabei die Plätze tauschen (Partnertausch). In dieser Sprache ist ein **Raupenschritt**
   eine Lücke, die über eine benachbarte Lücke der anderen Sorte hinwegspringt: zwei Sprünge, kein J. Ergebnis in
   1500 von 1500 Prozessen:

       Vorzeichen = (−1)^(Zahl der Schritte) · sgn(π der ↑-Lücken) · sgn(π der ↓-Lücken)

   642 davon mit Platztausch innerhalb einer Sorte; 25 416 J-Schritte und 3163 Raupenschritte. Falsche Regeln
   treffen nur in gut der Hälfte der Fälle: ohne (−1)^Schritte 789, ohne Permutation 858, nur J-Schritte gezählt 833.
2. **Bedeutung.** Mit diesen Schritten sind die Löcher zwei Sorten gewöhnlicher Fermionen:
   - innerhalb einer Sorte fermionisch (Platztausch −1);
   - zwischen den Sorten kein Vorzeichen;
   - kein Phase-String vom Hintergrund;
   - jeder Schritt mit festem Hüpfvorzeichen.

   Néel und Reparaturen verschwinden aus der Vorzeichenbuchhaltung. Übrig bleibt ein Zweisortenproblem:
   - Jede Sorte hüpft auf ihrem Untergitter (J-Schritte, diagonal und gerade um zwei Plätze). Das ergibt das Band
     ∝ (cos kₓ + cos k_y)² mit Minimum auf der Linie |kₓ| + |k_y| = π (7.4).
   - Der Raupenschritt ist korreliertes Hüpfen über die andere Sorte hinweg, ohne J. Er wirkt nur, wenn eine ↑- und
     eine ↓-Lücke benachbart sind. Das ist Bewegungsenergie, die nur Nachbarpaare haben, also eine Anziehung. Dazu
     kommt die Austausch-Anziehung: Benachbarte Löcher zerstören zusammen eine Bindung weniger.

   Die Zwei-Loch-Regel aus 7.12 ist der Fall einer Lücke je Sorte.
3. **ODLRO bei endlicher Dotierung (Argument).** Jede Sorte füllt einen Fermi-See im selben Band. Gepaart werden
   eine ↑-Lücke bei k und eine ↓-Lücke bei −k, die dieselbe Energie haben. Der Paarkanal ist nie destruktiv (7.12:
   G² ≥ 0; im Impulsraum divergiert die Paarsuszeptibilität logarithmisch). Darum führt jede Anziehung zu
   Cooper-Paaren und zu einem BCS-Zustand mit ODLRO. Die Paare überlappen, die Partner wechseln ständig; der
   Platztausch innerhalb einer Sorte (−1) ist die Antisymmetrie des BCS-Zustands. Kein Loch muss bei einem
   bestimmten Partner bleiben.
4. **Paarsymmetrie (Argument).** Die Anziehung wirkt auf Nachbarbindungen mit einer A- und einer B-Lücke. Für sie
   gibt es den erweiterten s-Kanal (cos kₓ + cos k_y) und den d-Kanal (cos kₓ − cos k_y). Die Fermi-Fläche liegt
   nahe der Linie |kₓ| + |k_y| = π, und dort ist cos kₓ + cos k_y = 0. Der s-Kanal hat an der Fermi-Fläche also
   kein Gewicht; es bleibt die d-Welle. Das passt zu 7.3 (Überlapp mit der s-Welle 0).
5. **Wo das Bild bricht.** Bei Schritten, die die Sorten mischen:
   - Sprünge ohne sofortige Reparatur, deren Spur die andere Sorte kreuzt (4.1);
   - Fluktuationen im Hintergrund (kein exaktes Néel).

   Dort entstehen Vorzeichen zwischen den Sorten (gegenseitige Statistik, Phase-String). Wie häufig das bei
   endlicher Dotierung ist und ob es die Paarung zerstört, ist offen.

### 7.14 Der Hintergrund als Filter: die ungerade Menge

Gesucht ist ein Mechanismus, bei dem der Hintergrund selbst durch Interferenz entscheidet, welche
Transferoperatoren langreichweitig kohärent bleiben, und nicht die Statistik der Löcher.

**Exakte Regel (Umformung von 7.6).** n↑ zählt Lochsprünge, die ein ↑-Elektron bewegen. Nach Elektronen sortiert:

    Vorzeichen = (−1)^(Zahl der ↑-Elektronen in O) · sgn(π_Loch)

Dabei ist O die **ungerade Menge**: die Elektronen, die von Lochsprüngen (aller Löcher zusammen) ungerade oft
bewegt wurden. Elektronen, die zweimal bewegt wurden, und alle Vertauschungen spielen keine Rolle. Für zwei
Prozesse mit gleichem Anfang und Ende ist das relative Vorzeichen entsprechend (−1)^(↑ in O₁ + ↑ in O₂) mal die
relative Lochpermutation.

Bekannte Fälle als Probe:
- Raupe: Jedes Elektron wird zweimal bewegt, O ist leer, also immun (7.6).
- Einzelnes Loch: O ist seine Spur. Eine Fluktuation, die ein Spur-Elektron mit einem von außerhalb tauscht,
  ändert die Parität (7.6).
- J-Wende: O ist die bewegte Bindung c–a′. Umklappen genau dieser Bindung lässt die Zahl der ↑ gleich, Brechen
  ändert sie (Tabelle in 7.9).
- J-Schritt einer Lücke: O enthält ein ↑ und ein ↓, also −1 pro Schritt (7.4, 7.13).

**Wie der Hintergrund filtert (Idee).** Der Hintergrund bestimmt das Vorzeichen nur über die Spins in O:
- Kommen die Spins in O in Paaren vor, die in jeder Komponente des Hintergrunds antiparallel sind (z. B. ein
  Singulett mit beiden Enden in O), ist die Parität fest, und alle Varianten interferieren konstruktiv.
- Jeder Spin in O, dessen Partner nicht in O liegt, trägt in den Komponenten des Hintergrunds verschiedene
  Paritäten. Varianten mit entgegengesetztem Vorzeichen löschen sich dann teilweise aus.

**Folgen:**
- *Ein Loch* kann die Partner der Spins, die es bewegt, nicht mitbewegen (außer entlang der eigenen Spur). Seine
  Varianten löschen sich aus.
- *Zwei oder mehr Löcher* können sich gegenseitig neutralisieren: Ein zweites Loch bewegt die Partner der Spins,
  die das erste bewegt hat. Wann und wo das geschieht, ist gleichgültig. Die Phase ist so lange in der Bindung des
  Hintergrunds gespeichert, bis der Partner kommt. Die Löcher können sich also trennen, und welches Loch welches
  neutralisiert, entscheiden die Bindungen des Hintergrunds, nicht die Löcher. Auch Ketten sind möglich: Ein Loch
  übernimmt, wo ein anderes aufgehört hat.
- Nach der Summe über die Varianten des Hintergrunds bleiben Prozesse übrig, deren ungerade Menge neutral ist. Ihre
  Nettowirkung auf die Ladung ist ein Transfer von Lochpaaren. Ein Einloch-Transfer über große Entfernungen bleibt
  nicht übrig.

**Was für ODLRO noch fehlt.** Neutral heißt nur: Jeder Prozess hat ein festes Vorzeichen, unabhängig von der
Komponente des Hintergrunds. Verschiedene neutrale Prozesse können aber verschiedene feste Vorzeichen haben. Jedes
Singulett in O trägt genau ein ↑ bei, also gilt

    Vorzeichen eines neutralen Prozesses = (−1)^(Zahl der Singuletts in O) · sgn(π_Loch)

Zwei neutrale Paartransfers zum selben Ziel addieren sich nur, wenn diese Größe für beide gleich ist (bis auf
einen Faktor, der nur vom Ziel abhängt). Genau hier sitzt die bekannte Frustration: Die Raupe hat O = ∅, eine
J-Wende ein Singulett in O. Gerade Raupenschritte und J-Wenden passen nicht zusammen (7.3, 7.9). Die eigentliche
Frage lautet also: Gibt es einen Grund, weshalb (−1)^(Singuletts in O) · sgn(π_Loch) für alle neutralen
Paartransfers nur vom Ziel abhängt? Wie weit sich Partner seitlich trennen können, hängt zwar von der Länge der
Bindungen ab, ist für ODLRO aber nicht entscheidend: Entlang der Spur können Partner beliebig weit
auseinanderliegen (Folgen auf derselben Spur gibt O = ∅ ohne jede Bedingung).

### 7.15 Welche Schritte stören, und warum ein Singulett-Hintergrund genau sie aussortiert

Die Frage aus 7.14: Gibt es einen Grund, weshalb (−1)^(Singuletts in O) · sgn(π_Loch) für alle verbleibenden
Paartransfers nur vom Ziel abhängt?

1. **Nur ein Schritttyp stört.** In der d-Wellen-Eichung aus 7.3 sind Pivot, Domino und alle Raupenschritte mit
   Wende positiv (392/392). Negativ sind genau die geraden Raupenschritte (192/192). Die gesamte Frustration der
   Paarbewegung kommt also aus einem einzigen Schritttyp.
2. **Was ein Raupenschritt mit einem Singulett macht.** Ein Raupenschritt verschiebt genau ein Elektron um zwei
   Plätze, vom Zielplatz des Kopfes auf den alten Platz des Schwanzes.
   - Beim geraden Schritt geht es zwei Plätze geradeaus. Sein Singulett-Partner (ein Nachbar des Ausgangsplatzes)
     ist danach drei oder √5 Plätze entfernt, nie ein Nachbar. Das Singulett wird gestreckt.
   - Beim Schritt mit Wende geht es diagonal. Das Singulett bleibt genau dann ein Nachbarpaar, wenn der Partner an
     der vierten Ecke der Plakette sitzt. Dann ist der Schritt derselbe Übergang wie ein Pivot (Paar und Singulett
     drehen sich auf der Plakette).
   - Pivot und Domino nehmen das Singulett als Ganzes mit; es bleibt ein Nachbarpaar.
3. **Folge.** In einem Hintergrund aus Nachbar-Singuletts bleiben genau die Plakettenbewegungen im Singulett-Raum:
   Paar und paralleles Singulett auf einer Plakette drehen sich (Pivot bzw. Raupenschritt mit Wende) oder
   tauschen (Domino bzw. zwei Raupenschritte mit Wende). Alle sind in der d-Wellen-Eichung positiv; Übergänge zum
   selben Ziel addieren sich (+1 in der Tabelle von 7.9). Die störenden geraden Raupenschritte führen aus dem
   Singulett-Raum heraus. **Der Singulett-Hintergrund sortiert also genau die frustrierenden Schritte aus.**
   Im Néel ist es umgekehrt: Dort ist der gerade Raupenschritt der natürlichste Schritt, weil er Néel ohne J
   wiederherstellt. Deshalb ist die Paarbewegung im Néel frustriert (flaches Band).
4. **Grenze des Arguments, und wo die Interferenz hineinkommt.** Das Aussortieren in 3. beruht nicht auf
   Vorzeichen, sondern darauf, dass ein gestrecktes Singulett den Singulett-Raum verlässt (es kostet J-Energie und
   überlappt weniger). Verboten ist der gerade Raupenschritt nicht: J kann das gestreckte Singulett wieder in
   Nachbar-Singuletts umordnen. Ob dieser Umweg mit dem d-Wellen-Vorzeichen zurückkommt (harmlos) oder mit dem
   falschen (Restfrustration der Ordnung J/ΔE), ist eine reine Interferenzfrage. Ebenso offen ist, ob die
   Resonanzen der Singuletts untereinander mit der d-Wellen-Eichung der Paarbewegung verträglich sind. Beides lässt
   sich auf einem kleinen Gitter exakt prüfen.

### 7.16 Idee: Vorzeichen als dynamisches Feld des Hintergrunds

Ausgangspunkt (Beobachtung aus der Diskussion, 7.6): Vertauschungen des Hintergrunds machen schon für ein
einzelnes Loch beide Vorzeichen erreichbar, mit demselben Anfang und demselben Ende.

1. **Kein Vorzeichen ist fest.** Nach 7.14 ist das Vorzeichen (−1)^(↑ in O) · sgn(π_Loch). Jede Vertauschung,
   die ein Elektron aus O gegen eines von außerhalb tauscht, dreht es um. Das gilt für jeden Prozess mit
   nichtleerem O (Einzelloch, J-Wende). Selbst die Raupe ist nur in niedrigster Ordnung starr: Tauscht eine
   Vertauschung zwischen Kopf- und Schwanzsprung das verschobene Elektron gegen einen Nachbarn aus, werden zwei
   entgegengesetzte Elektronen je einmal bewegt, und das Vorzeichen kehrt sich um. Bei genügend Bewegung im
   Hintergrund sind die Vorzeichen also keine Eigenschaft der Lochwege mehr, sondern hängen davon ab, wie der
   Hintergrund mit den Löchern korreliert ist. Der Hintergrund wirkt wie ein dynamisches Z2-Eichfeld:
   Welcher Fluss (−1) auf einer geschlossenen Bahn liegt, ist nicht vorgegeben.
2. **Was der Grundzustand wählt.** Für Bosonen gilt die diamagnetische Ungleichung: Jeder Fluss erhöht die
   Energie, am tiefsten liegt der flussfreie Fall, in dem sich alle Wege addieren. Paare sind untereinander Bosonen
   (Platztausch +1, 7.2). Können die Vorzeichen vom Hintergrund eingestellt werden, gewinnen die Paare
   Bewegungsenergie, wenn der Hintergrund sich so korreliert, dass alle Paarwege konstruktiv sind. Dann bewegen
   sich die Paare ohne Vorzeichen; nach Perron–Frobenius sind alle Amplituden positiv, und der Paartransfer
   ⟨Δ†(R)Δ(0)⟩ ist eine Summe positiver Terme. Das ist ODLRO, sobald die Paare beweglich sind.
3. **Warum der Kanal das Paar ist** (korrigiert, siehe 7.17, Schritt 8). Das Tauschvorzeichen der Löcher ist
   nicht fest; der Hintergrund kann es ausgleichen. Dass ODLRO im Paarkanal erscheint, folgt aus dem
   Pauli-Prinzip für die Elektronen, nicht aus der Statistik der Löcher.
4. **Hinweis in den bisherigen Zahlen.** Je mehr Freiheit der Hintergrund hat, desto weniger frustriert ist die
   Paarbewegung: 0,708 (sauber) → 0,880 (höchstens eine Fluktuation) → 0,932 (höchstens zwei), Gewicht 0,3, 7.8.
   Mit zufälligen statt echten Vorzeichen bleibt es bei 0,721. Die Verbesserung kommt also aus der echten
   Vorzeichenfreiheit, nicht aus Verdünnung. Das ist der Trend, den 1.–2. vorhersagen.
5. **Was offen ist.** (a) Kosten: Die nötigen Korrelationen des Hintergrunds kosten J-Energie; sie müssen
   billiger sein als der Gewinn an Paarbewegung. (b) Läuft der Trend aus 4. gegen 1, wenn der Hintergrund noch mehr
   Freiheit bekommt (drei und mehr Fluktuationen)? (c) Die Vergleichsrechnung für viele einzelne Löcher.

### 7.17 Argumentation: Fluktuationen existieren ⇒ ODLRO

Markiert ist, was exakt ist, was ein Satz ist und was Annahme oder Argument ist.

**Schritt 1 (exakt): Woher Vorzeichen kommen.** Für jeden Prozess gilt
Vorzeichen = (−1)^(↑-Elektronen in O) · sgn(π_Loch). O sind die Elektronen, die von Lochsprüngen ungerade oft
bewegt wurden. Vertauschungen tragen kein Vorzeichen; reine Spinprozesse sind vorzeichenfrei (Marshall). Jeder
Lochsprung trägt also den Faktor −1 genau dann, wenn er ein ↑ bewegt: Die Spins des Hintergrunds sind ein
Z2-Feld, an das die Löcher koppeln (7.1, 7.6, 7.14).

**Schritt 2 (exakt): Ohne Fluktuationen ist dieses Feld starr, und die Paarbewegung ist frustriert.** Im sauberen
Néel liegt das Feld fest. Die Raupe hat O = ∅, ihr Vorzeichen ist allein der Lochaustausch. Nach dem Kopfsprung
sitzt das verschobene Elektron zwischen Nachbarn mit gleichem Spin, keine Vertauschung kann es erreichen, und
das Vorzeichen ist nicht veränderbar. Ergebnis: flaches Band, die gerade Raupe widerspricht den J-Wenden
(7.3, 7.9). Kein Grund für ODLRO.

**Schritt 3 (exakt): Mit Fluktuationen wird jedes Vorzeichen verstellbar.** Tauscht eine Fluktuation ein
Elektron aus O gegen eines von außerhalb, ändert sich die Parität. Anfang und Ende bleiben gleich, wenn sie
zurückgetauscht wird (7.6). Für die Raupe gilt dasselbe, sobald neben ihrer Spur eine Fluktuation existiert: Eine
Vertauschung zwischen Kopf- und Schwanzsprung macht O zu zwei entgegengesetzten Elektronen und dreht das
Vorzeichen. Mit Fluktuationen ist das Z2-Feld also dynamisch. Das Vorzeichen eines Paartransfers ist keine
Eigenschaft der Geometrie mehr, sondern hängt davon ab, wie der Hintergrund mit den Löchern korreliert ist.

**Schritt 4 (exakt): Paare sind Bosonen.** Tauschen zwei Paare die Plätze, ist O leer und die Lochpermutation
gerade: +1 (7.2).

**Schritt 5 (Satz, dann Annahme): Bosonen verdrängen Frustration.** Diamagnetische Ungleichung: Für Bosonen
mit positiven Hüpfamplituden erhöht jeder Fluss (jedes −1 auf einer geschlossenen Bahn) die Grundzustandsenergie.
Frustration bringt den Paaren also nur Nachteile. Weil das Feld nach Schritt 3 verstellbar ist, kann der
Grundzustand den Hintergrund so mit den Paaren korrelieren, dass ihre Bahnen konstruktiv werden; der Preis ist
Austauschenergie des Hintergrunds. *Annahme:* Der Gewinn an Paarbewegung übersteigt diesen Preis. Das ist das
Gegenstück zur Verdrängung eines Magnetfelds durch ein Kondensat.

**Schritt 6 (Hinweis aus den Rechnungen): Der Grundzustand tut das tatsächlich.** Ein Paar bindet im
Grundzustand Fluktuationen (54 % mit zwei bei Gewicht 0,3, 7.8). Mit wachsender Freiheit des Hintergrunds sinkt
die Frustration der Paarbewegung: 0,708 → 0,880 → 0,932, mit zufälligen Vorzeichen dagegen nur 0,721. Auf dem
Torus hat das Paar mit Fluktuationen ein Band mit eindeutigem Maximum bei K = 0 (7.9). Die Spinwolke ist die
Korrelation aus Schritt 5.

**Schritt 7 (Standardargument): ODLRO.** Ohne Frustration ist die Paardynamik vorzeichenfrei. Nach
Perron–Frobenius sind dann alle Amplituden des Grundzustands positiv, und ⟨Δ†(R)Δ(0)⟩ ist eine Summe positiver
Terme. Auslöschung ist auf keiner Entfernung möglich. Sind die Paare beweglich (kein Kristall), kondensieren sie
und G(R) → konstant ≠ 0.

**Schritt 8: Warum der Kanal das Paar ist (korrigiert).** Mit „Tausch zweier Löcher“ ist ein geschlossener
Prozess gemeint, nach dem zwei Löcher auf den Startplätzen des jeweils anderen stehen: sgn(π_Loch) = −1. Dieses
−1 tritt aber nie allein auf, sondern immer zusammen mit (−1)^(↑ in O); physikalisch ist nur das Produkt.
- Tauschen die Löcher auf demselben Weg, ist O leer, und es bleibt −1 (7.2).
- Tauschen sie auf verschiedenen Wegen, hängt das Vorzeichen von O ab. Schon im sauberen Néel gibt es Lochtausch
  mit Gesamtvorzeichen +1: Pivot hin und Raupe mit Wende zurück vertauschen die beiden Löcher, und das ↑ der
  bewegten Bindung gleicht das aus (7.3).
- Mit Fluktuationen ist es nach Schritt 3 verstellbar.

Die frühere Fassung behauptete, das Tauschvorzeichen einzelner Löcher könne kein Hintergrund beseitigen. Das ist
falsch: Die Statistik der Löcher ist durch den Hintergrund bekleidet. Das passt zum Phase-String-Bild, in dem
Löcher als Bosonen beschrieben werden (aus dem Gedächtnis). Schritt 4 und 5 zeichnen Paare daher nicht vor
einzelnen Löchern aus.

Der Kanal, in dem ODLRO erscheinen kann, ist trotzdem das Paar, und zwar aus einem anderen Grund. Ein Loch
entsteht durch c_iσ, einen fermionischen Operator. Die Einteilchen-Dichtematrix ⟨c†(R)c(0)⟩ hat nach dem
Pauli-Prinzip Eigenwerte ≤ 1. Sie kann also keine makroskopische Besetzung und kein ODLRO haben (C. N. Yang, 1962,
aus dem Gedächtnis). Das erste mögliche ODLRO liegt in der Zweiteilchen-Dichtematrix, dem Paartransfer
⟨Δ†(R)Δ(0)⟩. Schritte 1–7 dürfen auch die Bewegung einzelner Löcher konstruktiv machen; langreichweitig sichtbar
wird das nur im Paarkanal. Die Löcher eines Paares müssen dafür nicht beisammenbleiben. Offen ist, ob neben der
Verdrängung der Frustration noch eine Paarung im Spinhintergrund nötig ist, damit der Paarkanal tatsächlich
kondensiert.

**Kurz:** Fluktuationen machen die Vorzeichen zu einem dynamischen Feld (exakt). Paare sind Bosonen (exakt); das
Tauschvorzeichen auch einzelner Löcher ist durch den Hintergrund ausgleichbar (Schritt 8). Bosonen senken ihre
Energie, indem sie Frustration verdrängen (Satz für feste Felder). Ein frustrationsfreies Bosonensystem hat
positive Amplituden und kondensiert (Standard). Dass das ODLRO im Paarkanal liegt, folgt aus dem Pauli-Prinzip
für die Elektronen. Die einzige echte Annahme steckt in Schritt 5:
Die Verdrängung muss sich energetisch lohnen. Schritt 6 zeigt, dass der Grundzustand im kleinen System genau in
diese Richtung geht.

### 7.18 Annahme „jede Lochbewegung gibt es mit beiden Vorzeichen“: was daraus folgt (`spurpaare.py`)

**Annahme A.** Fluktuationen des Hintergrunds machen jede Lochbewegung mit beiden Vorzeichen möglich.

**Präzisierung über die Regel.** Das Vorzeichen ist (−1)^(↑ in O) · sgn(π_Loch). „Beide Vorzeichen“ heißt: Für ein
Elektron in O ist offen, ob es ↑ oder ↓ war. Zwei Ausnahmen bleiben fest:
- Elektronen, die gar nicht in O liegen (zweimal bewegt);
- zwei Elektronen in O, die im selben Prozess miteinander vertauscht wurden. Eine Vertauschung gibt es nur zwischen
  entgegengesetzten Spins, sie *beglaubigt* also: genau ein ↑.

Summiert man über die Varianten, bleibt ein fester Beitrag nur von Prozessen, deren O leer ist oder nur aus
beglaubigten Paaren besteht. Alle anderen heben sich weg, umso stärker, je mehr unbeglaubigte Elektronen O enthält.

**Schritt 1: Einzelne Löcher.** Ein Loch, das ohne Reparatur über R Plätze läuft, hat R unbeglaubigte Elektronen in O.
Sein Beitrag fällt mit R weg. Übrig bleibt der J-Schritt: Seine Reparatur vertauscht genau die beiden bewegten
Elektronen und beglaubigt sie. Einzelne Löcher bewegen sich unter A also nur langsam (über J) und frustriert (7.4).

**Schritt 2: Zwei Löcher auf gemeinsamer Spur.** Bewegt das zweite Loch dieselben Elektronen wie das erste, ist O
leer, gleich mit welchem Abstand und welcher Verzögerung. Solche Prozesse kann der Hintergrund nicht verändern.
Schneller Ladungstransport ohne Auslöschung bleibt unter A nur für zwei Löcher auf gemeinsamer Spur. Das wählt
Paare aus, ohne dass sie Nachbarn bleiben müssen.

**Schritt 3: Das allein reicht nicht (gerechnet).** Bei leerem O ist das Vorzeichen allein sgn(π_Loch). Der Tausch der
beiden Löcher innerhalb des Paares gibt −1. Das Paar mit nur solchen Prozessen ist stark frustriert:

| höchstens Sprünge pro Prozess (Abstand entlang der Spur) | Impulse am Maximum (von 72) | Beweglichkeit (w = 0,3 / 0,5 / 0,7) | Verhältnis (w = 0,5) |
|---|---|---|---|
| 2 (reine Raupe) | 72, flach | 0 / 0 / 0 | 0,667 |
| 4 (Abstand bis 1) | 72, flach | 0 / 0 / 0 | 0,583 |
| 6 (Abstand bis 2) | 1 | 0,012 / 0,016 / 0,017 | 0,525 |
| 8 (Abstand bis 3) | 1 | 0,024 / 0,028 / 0,029 | 0,494 |

Beweglichkeit wie in 7.9 (Krümmung des obersten Bandes, geteilt durch dieselbe Größe ohne Vorzeichen); jeder Prozess
hat das Gewicht w^(Sprünge). Das Band bleibt bis Abstand 1 exakt flach. Danach ist es kaum gekrümmt (2–3 % des
vorzeichenfreien Werts), und die Frustration nimmt mit dem Abstand sogar zu. Auftrennen entlang der Spur macht
das Paar also nicht beweglich.

**Schritt 4: Was das Paar beweglich macht: beglaubigte Vertauschungen.** Eine J-Wende mit ihrer Reparatur
vertauscht die beiden bewegten Elektronen. Die Vertauschung beglaubigt genau ein ↑ in O und liefert das −1, das den
Lochtausch ausgleicht (7.3). Solche Prozesse überleben A. Ist die Bindung gerade parallel, gibt es diesen Prozess
einfach nicht; er tritt seltener auf, bekommt aber kein falsches Vorzeichen. Mit ihnen ist das Paar beweglich,
mit eindeutigem Bandmaximum bei K = 0 (7.9: Beweglichkeit 0,74 / 0,46 / 0,23).

**Schritt 5: Der letzte Widerspruch, und was A daran ändert.** Gegen die J-Wenden steht nur der gerade
Raupenschritt (O leer, falsches Vorzeichen in der d-Wellen-Eichung, 7.3). Unter A hat auch er eine beglaubigte
Variante: Eine Vertauschung zwischen Kopf- und Schwanzsprung tauscht das verschobene Elektron gegen einen
entgegengesetzten Nachbarn. Das beglaubigt ein ↑ in O und dreht das Vorzeichen in die passende Richtung. Im
sauberen Néel gibt es diese Variante nicht; dort steht das verschobene Elektron zwischen gleichen Spins. Mit
Fluktuationen gibt es sie. Fluktuationen schwächen also genau den einen störenden Schritt. Das passt zum Trend in
7.8 (0,708 → 0,880 → 0,932 mit wachsender Freiheit des Hintergrunds, mit zufälligen Vorzeichen 0,721).

**Schritt 6: ODLRO.** Unter A bleiben als schnelle, kohärente Ladungsträger Paare: Ihr Tausch untereinander gibt +1,
und sie haben ein Band mit eindeutigem Maximum bei K = 0. Ein dünnes Gas solcher Bosonen kondensiert dort, also
⟨Δ†(R)Δ(0)⟩ → konstant. Einzelne Löcher bleiben langsam und frustriert, und nach Pauli kann ihr Einteilchenkanal
ohnehin kein ODLRO tragen (7.17, Schritt 8).

**Zusammengefasst.** A wirkt doppelt:
- Fluktuationen löschen alles Unbeglaubigte aus und wählen damit Paare aus (Schritte 1–2).
- Vertauschungen, die Teil des Prozesses sind, beglaubigen Spinpaare und liefern das Vorzeichen, das den Lochtausch
  ausgleicht (Schritte 4–5).

Ohne den zweiten Teil wäre das ausgewählte Paar lokalisiert (Schritt 3). Offen ist, ob die beglaubigte Variante des
geraden Raupenschritts den Widerspruch ganz aufhebt oder nur abschwächt; das hängt von den Gewichten ab.

### 7.19 Die Löcher wandern einzeln zur neuen Position

Einwand gegen 7.18: Die beiden Löcher müssen nicht als Paar wandern. Es reicht, dass jedes einzeln an die neue
Position kommt.

1. **Exakt: Ein Schritt mit eigener Reparatur hat ein festes Vorzeichen.** Vorzeichen = (−1)^(n↑) · sgn(π_Loch),
   n↑ = Zahl der Sprünge, die ein ↑ bewegen. Ein J-Schritt besteht aus zwei Sprüngen; danach vertauscht die
   Reparatur genau die beiden bewegten Elektronen. Eine Vertauschung gibt es nur zwischen entgegengesetzten Spins.
   Also bewegt genau einer der beiden Sprünge ein ↑, und jeder solche Schritt trägt −1, in jedem Hintergrund,
   auch mit Fluktuationen. Eine Bewegung aus solchen Schritten hat das Vorzeichen
   (−1)^(Zahl der Schritte) · sgn(π_Loch). Der Hintergrund entscheidet nur, welche Schritte möglich sind (Gewichte),
   nicht ihr Vorzeichen.
2. **Unter A** (7.18) löschen sich alle anderen Beiträge aus. Übrig bleiben für jedes Loch genau diese Schritte.
3. **Zwei Löcher wandern unabhängig.** Ein Loch sitzt auf A (↑-Lücke), das andere auf B (↓-Lücke). Wege und Zeiten
   sind beliebig, die Löcher können weit auseinander sein. Jedes bleibt auf seinem Untergitter (Schrittweite 2),
   tauschen können die beiden also nicht. Vorzeichen = (−1)^(n₁) · (−1)^(n₂). Die Summe über alle Wege zerfällt in
   ein Produkt: Transferamplitude = G_A(R) · G_B(R). Beide Lücken haben dieselbe Regel (−1 pro Schritt) und in einem
   gleichförmigen Hintergrund dieselben Gewichte, also G_A = G_B und

       Paartransfer = G(R)² ≥ 0.

   Wie stark G(R) selbst im Vorzeichen schwankt (das einzelne Loch ist frustriert, −1 pro Dreieck, 7.4), spielt keine
   Rolle: Der Paartransfer ist auf jeder Entfernung nicht negativ. Numerisch gestützt durch 7.12 (Produktregel 1563
   von 1563, deckungsgleiche Wege immer +1).
4. **Rolle der Fluktuationen.** Sie löschen genau die Beiträge aus, deren Vorzeichen vom Hintergrund abhängt. Das
   sind die Beiträge, die den beiden Partnern verschiedene Vorzeichen geben und das Quadrat zerstören könnten
   (Paarbrecher). Übrig bleibt ein Vorzeichen, das für beide Partner gleich ist; ihr Produkt ist ein Quadrat.
5. **Von „nie negativ“ zu ODLRO.** Positivität auf jeder Entfernung ist nötig, reicht aber allein nicht: Für zwei
   freie Löcher ist G(R)² zwar positiv, fällt aber mit R ab, und eine makroskopische Besetzung des Paarkanals entsteht
   nicht. Dazu braucht es eine Anziehung zwischen ↑- und ↓-Lücke. Weil der Paarkanal eine Summe von Quadraten ist und
   nie destruktiv interferiert, reicht jede Anziehung (Cooper-Argument in dieser Sprache). Hinweise auf die
   Anziehung: 6.1 (mit Vorzeichen gebunden, ohne nicht) und der Raupenschritt, den nur benachbarte Löcher haben
   (7.13).

Damit ist 7.18, Schritt 2 („schneller Transport nur auf gemeinsamer Spur“) für die Frage nach ODLRO nicht
entscheidend. Auch einzeln, langsam wandernde Löcher geben einen Paartransfer ohne Auslöschung.

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
- Abschnitt 7 nutzt effektive Modelle: fester Néel-Hintergrund (7.6–7.9 mit höchstens einer bzw. zwei Fluktuationen nahe am Paar),
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

## 10. Was gefunden wurde

Stand am Ende der Diskussion. „Gerechnet“ heißt: mit den Skripten in diesem Verzeichnis geprüft. „Argument“
heißt: begründet, aber nicht gerechnet.

**Gerechnet**

1. Das Repo-Vorzeichen eines Prozesses ist das Vorzeichen seines Gewichts in exp(−τH) (Abschnitt 1).
   Prozesse, die sich nur im Reparatur-Zeitpunkt unterscheiden, interferieren.
2. Besuchen die Wege aller Löcher keinen Platz zweimal, hängt das Vorzeichen nicht vom Reparatur-Zeitpunkt ab
   (Abschnitt 2). Eine Abhängigkeit gibt es nur bei Kreuzungen oder Schleifen (Abschnitt 4).
3. Regel für geschlossene Prozesse: Vorzeichen = (−1)^s · sgn(π_Loch) = (−1)^(n↑) · sgn(π_Loch) (7.1, 7.6).
   In der zweiten Form tragen Vertauschungen kein Vorzeichen.
4. Ein einzelnes Loch wechselt das Vorzeichen, sobald eine Fluktuation seinen Weg berührt (Phase-String). Das
   Raupenpaar nicht: Der Schwanz bewegt jeden Spin ein zweites Mal (7.6).
5. Die Raupe allein ist durch Interferenz lokalisiert: flaches Band, bei dem sich der Lochaustausch beim Wenden
   (−1) nicht ausgleicht (7.3; auf dem Torus liegen alle 72 Impulse am Maximum, 7.9).
6. Pivot und Domino (Paarschritte mit J) sind unfrustriert in der d-Wellen-Eichung (7.3; Torus: Verhältnis
   1,000). Erst sie machen das Paar beweglich; die Raupe trägt zur Beweglichkeit nichts bei (7.9).
7. Eine J-Wende hat ein festes Vorzeichen, wenn die Nachbarbindung, deren Spins sie bewegt, antiparallel ist
   (in beiden Orientierungen dasselbe). Trägt die Bindung parallele Spins, kehrt es sich um (7.9).
8. Tauschen zwei Raupenpaare die Plätze, ergibt das +1: Paare sind untereinander bosonisch (7.2).
9. Exakte Diagonalisierung 4×4: Mit Vorzeichen sind zwei Löcher gebunden (E_b = −0,35 t bei J/t = 0,4), ohne
   Vorzeichen nicht (+0,19 t) (6.1).
10. Das Paar mit Spinwolke (bis zu zwei Fluktuationen) hat auf dem Torus ein gekrümmtes Band mit eindeutigem
    Maximum bei K = 0 (7.9).
11. Zwei getrennt laufende Löcher (↑-Lücke und ↓-Lücke, sofort reparierte J-Schritte): Das Vorzeichen ist das
    Produkt der Einzelvorzeichen, bei deckungsgleichen Wegen immer +1, auch bei großem Abstand und gekreuzten
    Wegen (7.12).
12. Viele Löcher mit J-Schritten und Raupenschritten: Vorzeichen = (−1)^Schritte · sgn(π der ↑-Lücken) ·
    sgn(π der ↓-Lücken), also zwei Sorten Fermionen ohne Vorzeichen zwischen den Sorten (7.13).

**Korrigiert**

- 7.6: Die Zeile „sauberes Néel, nur Prozesse mit Vertauschung: 1,000“ war trivial, weil die Aufzählung Pivot
  und Domino zerlegt hatte.
- 7.7/7.8: Die Auswahl „Raupe + freie J-Fluktuationen“ enthielt zerlegte J-Wenden. Die Deutung als „gebundene
  Fluktuation“ ist überholt.
- Das Verhältnis λ_max(M)/λ_max(|M|) misst nur lokale Frustration und ist kein Maß für ODLRO. Beispiel: Raupe
  mit Fluktuationen, aber ohne J-Wenden, hat 0,92 und ist trotzdem fast unbeweglich.
- 7.9: Der dort vorhergesagte Abfall mit der Entfernung gilt nur für eingefrorene Fluktuationen (7.11).

**Idee für ODLRO (Argument, 7.12 und 7.13)**

Mit sofort reparierten J-Schritten und Raupenschritten sind die Löcher zwei Sorten gewöhnlicher Fermionen
(↑-Lücken auf A, ↓-Lücken auf B): fermionisch innerhalb einer Sorte, kein Vorzeichen zwischen den Sorten, kein
Phase-String. Der Raupenschritt ist korreliertes Hüpfen benachbarter ↑- und ↓-Lücken, also eine Anziehung. Bei
endlicher Dotierung folgt das BCS-Bild: Cooper-Paare aus ↑-Lücke bei k und ↓-Lücke bei −k, überlappend, mit
ständig wechselnden Partnern, und wegen der Fermi-Fläche nahe |kₓ| + |k_y| = π in der d-Welle. Für zwei Löcher
heißt das: Ein Paar muss nicht zusammenbleiben. Seine beiden Löcher sind verschiedene Sorten (↑-Lücke, ↓-Lücke). Jede
vertauscht nur Elektronen ihrer Sorte, und das Vorzeichen ist das Produkt zweier Faktoren, je einer pro Loch.
Trennen sich die Löcher und treffen sich wieder, ist die Amplitude G(R)² und damit nie negativ: Die Phasen der
Partner heben sich auf, wie bei Cooper-Paaren aus Zeitumkehr-Partnern. Partnertausch gibt das fermionische −1
eines BCS-Zustands. Mit einer Anziehung (6.1) folgt Cooper-Paarung und daraus ODLRO; dieser letzte Schritt ist das
Standardargument. Der Reparatur-Zeitpunkt greift genau dort ein, wo eine Lücke die noch unreparierte Spur der
anderen kreuzt: Dann sind die Faktoren nicht mehr unabhängig (4.1, 7.12).

Vorstufe (7.11): Vertauschungen tragen kein Vorzeichen, sie erzeugen alternative Wege. Im Quanten-Grundzustand wird für jeden
Paarschritt über alle Varianten summiert. Übrig bleibt der konstruktive Überschuss, mit festem und überall
gleichem Vorzeichen: bei der Raupe 1, bei einer J-Wende −(4/3)⟨S·S⟩ der bewegten Bindung, also etwa 0,4 bis
0,5. Der Schritt von dort zu einem fest gebundenen Paar als Quasiteilchen ist schwach; 7.12 braucht ihn nicht. Ein
Singulett-Hintergrund (7.10) ist der Grenzfall mit größtem Überschuss.

**Der Hintergrund als Filter (7.14)**

Exakt: Vorzeichen = (−1)^(↑-Elektronen in O) · sgn(π_Loch), mit O = Elektronen, die ungerade oft von
Lochsprüngen bewegt wurden. Idee: Der Hintergrund löscht Varianten aus, in denen ein Spin in O keinen
antiparallelen Partner in O hat, und verstärkt die anderen. Löcher können sich gegenseitig neutralisieren, auch
getrennt und zeitversetzt, weil die Bindung des Hintergrunds die Phase speichert. Übrig bleiben Paartransfers.
Offen ist, ob alle neutralen Paartransfers zum selben Ziel auch dasselbe Vorzeichen haben: Das Vorzeichen ist
(−1)^(Singuletts in O) · sgn(π_Loch), und gerade Raupenschritte passen damit nicht zu J-Wenden.

**Warum ein Singulett-Hintergrund die störenden Schritte aussortiert (7.15)**

In der d-Wellen-Eichung sind alle Paarschritte positiv außer den geraden Raupenschritten. Ein gerader
Raupenschritt verschiebt ein Elektron um zwei Plätze geradeaus und streckt dabei sein Singulett. In einem
Hintergrund aus Nachbar-Singuletts bleiben nur die Plakettenbewegungen (Paar und Singulett drehen oder tauschen),
und die sind alle positiv. Offen: ob der Umweg über ein gestrecktes Singulett und J zurück das richtige Vorzeichen
trägt, und ob die Resonanzen der Singuletts mit der d-Wellen-Eichung verträglich sind.

**Vorzeichen als dynamisches Feld (7.16)**

Weil Vertauschungen des Hintergrunds für jeden Prozess beide Vorzeichen erreichbar machen (in höherer Ordnung
selbst für die Raupe), sind die Vorzeichen keine feste Eigenschaft der Lochwege; der Hintergrund wirkt wie ein
dynamisches Z2-Eichfeld. Für Bosonen ist der flussfreie Fall energetisch am günstigsten (diamagnetische
Ungleichung). Paare sind Bosonen, also sollte sich der Hintergrund so korrelieren, dass alle Paarwege konstruktiv
sind; der Paartransfer ist dann eine Summe positiver Terme (ODLRO). Das Tauschvorzeichen auch einzelner Löcher ist
durch den Hintergrund ausgleichbar (korrigiert in 7.17); dass ODLRO im Paarkanal liegt, folgt aus dem Pauli-Prinzip
(Einteilchen-Dichtematrix der Elektronen hat Eigenwerte ≤ 1). Hinweis: Mit mehr Freiheit im Hintergrund sinkt die Frustration der Paarbewegung (7.8), und
zwar wegen der echten Vorzeichen. Offen: die Energiekosten, der Grenzwert bei mehr Freiheit, der Vergleich mit
einzelnen Löchern.

**Argumentationskette (7.17)**

Fluktuationen machen die Vorzeichen zu einem dynamischen Z2-Feld; ohne sie ist die Raupe starr und die
Paarbewegung frustriert (exakt). Paare sind Bosonen (exakt). Bosonen senken ihre Energie, indem sie Frustration
verdrängen (diamagnetische Ungleichung, für feste Felder ein Satz). Ohne Frustration sind alle Amplituden positiv,
und der Paartransfer zeigt ODLRO (Standard). Das Tauschvorzeichen einzelner Löcher ist ebenfalls ausgleichbar;
der Kanal ist das Paar, weil die Einteilchen-Dichtematrix der Elektronen nach Pauli kein ODLRO haben kann. Die einzige Annahme: Die Verdrängung lohnt sich energetisch. Die
gebundene Spinwolke und die sinkende Frustration in 7.8 zeigen, dass der Grundzustand diesen Weg geht.

**Annahme „beide Vorzeichen“ (7.18)**

Unter der Annahme, dass Fluktuationen jede Lochbewegung mit beiden Vorzeichen möglich machen, bleiben nur Prozesse
mit leerer ungerader Menge oder mit *beglaubigten* Paaren in O (zwei Elektronen, die im selben Prozess vertauscht
wurden, also sicher eines ↑). Das wählt Paare aus; einzelne Löcher bleiben langsam. Ein Paar mit nur leerer
ungerader Menge (gemeinsame Spur, auch mit Abstand) ist aber praktisch lokalisiert: Das Band ist bis Abstand 1 exakt
flach, danach liegt die Beweglichkeit bei 2–3 % (gerechnet). Beweglich wird es erst durch beglaubigte
Vertauschungen (J-Wende mit Reparatur), die den Lochtausch ausgleichen. Fluktuationen geben auch dem geraden
Raupenschritt, dem einzigen störenden Schritt, eine beglaubigte Variante mit passendem Vorzeichen.

**Einzeln wandernde Löcher (7.19)**

Ein J-Schritt mit eigener Reparatur hat immer genau einen Sprung, der ein ↑ bewegt, also das feste Vorzeichen −1,
in jedem Hintergrund (exakt). Wandern die beiden Löcher einzeln zur neuen Position, zerfällt die Summe über alle
Wege in ein Produkt G_A(R) · G_B(R) = G(R)² ≥ 0: Der Paartransfer ist auf jeder Entfernung nicht negativ, wie
stark G(R) selbst schwankt. Fluktuationen löschen genau die Beiträge aus, die den Partnern verschiedene
Vorzeichen geben würden. Für ODLRO braucht es zusätzlich eine Anziehung; wegen der Positivität reicht jede.

**Offen**

- wie häufig Schritte sind, die die Sorten mischen (unreparierte Spuren, Fluktuationen), und ob sie die
  Paarung zerstören (7.12, Punkt 7; 7.13, Punkt 5);
- wie stark die Anziehung zwischen ↑- und ↓-Lücke ist und welche Paarsymmetrie daraus folgt;
- der Vielteilchenschritt (Cooper-Paarung bzw. Kondensation);
- ob sich destruktive Wege sogar paarweise exakt aufheben (stärkere Form der Idee);
- echte Energien statt einfacher Gewichte, größere Systeme.

**Werkzeuge:** `tj_gui_web.html` (Web-GUI; zwei Bewegungen vergleichen und für Claude speichern),
`kohaerenz.py`, `fluktuationen.py`, `fluktuationen2.py`, `beweglichkeit.py`, `partner.py`, `ed_tj.py` (siehe Inhalt).

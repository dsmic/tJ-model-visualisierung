# tj_gui – t-J-Modell auf dem 8x8-Gitter

Interaktive tkinter-Oberfläche, um das Fermionen-Vorzeichen von Konfigurationen des t-J-Modells
(hardcore, keine Doppelbesetzung) nachzuvollziehen. Das Vorzeichen folgt der
Fermion-Ordnung von NetKet (`SpinOrbitalFermions`, NetKet-Dokumentation; die Konvention wurde in einem eigenen Notebook nachvollzogen):
Orbital-Index = `spin_idx*N + site`, `site = x*L + y`, down = Spin-Index 0, up = Spin-Index 1.

Start: `python tj_gui.py` (nur tkinter nötig; die Vorzeichenregel ist in Python nachgebaut. Für NetKet-Vergleiche eine Conda-Umgebung mit NetKet verwenden).

## Was das Programm zeigt

- 8x8-Brett mit offenen Rändern. Rot = up, blau = down, leerer Kreis = Loch.
- Hintergrund: hellgrün = Vorzeichen +1, hellviolett = Vorzeichen −1.
- Start: zufällige Néel-artige Konfiguration mit 7/8-Füllung (28 up, 28 dn, 8 Löcher;
  je 4 Löcher pro Untergitter).
- Klick auf eine Bindung zwischen zwei Nachbarplätzen vertauscht deren Inhalt. Das Vorzeichen
  wird dabei mit dem Jordan-Wigner-Vorzeichen des Prozesses multipliziert (Lochsprung =
  `c†_dst c_src`, up/dn-Austausch = vier Operatoren). Das Vorzeichen ist ein Produkt der
  Schrittvorzeichen ab +1, also pfadabhängig, keine absolute Funktion der Konfiguration.
  Die Implementierung wurde gegen die NetKet-Operatoren geprüft (1738/1738 Übereinstimmungen).

## Bedienelemente

Untere Leiste:

- **Neu (zufällig)**, **Zurück**, **Vor**, **Alles zurück** (bis zum Laden/Neu-Start).
- **Check Néel**: markiert Spins auf dem falschen Untergitter orange.
- **Lochspur löschen**: die graue Linie zeigt die Lochbewegungen seit dem Start.
- **Néel im Hintergrund**: hellrote/hellblaue Kreise, etwas größer als die Spins, zeigen, welche
  Farbe die Néel-Konfiguration an dem Platz hätte.
- Rechts unten: Anzahl falscher Néel-Plätze, Vorzeichen, Teilchenzahlen.

Rechts: **Sichern (s)** legt das aktuelle Bild (Konfiguration, Vorzeichen, Spur) in einem 4x4-Raster
ab; Klick auf ein Bild lädt es als neue Startkonfiguration. **Bilder löschen** leert das Raster.
Das Fenster skaliert beim Ändern der Größe mit.

Zweite Leiste (Vorzeichentest):

- Eingabefeld: Weg des Lochs als Folge aus `R`, `L`, `U`, `D` (rechts/links/oben/unten),
  z.B. `RRULLD`.
- **Test, Reparatur am Ende**: 100 zufällige Startkonfigurationen, der Weg wird ab einem zufälligen Loch
  gelaufen (Läufe über den Rand oder auf ein anderes Loch werden neu gewürfelt). Danach wird die
  Néel-Konfiguration durch Vertauschen benachbarter falscher up/dn-Paare (zufällige Paarwahl)
  repariert. Läufe, die sich so nicht reparieren lassen, werden verworfen und gezählt.
  Ausgewertet wird das Endvorzeichen.
- **Test, Reparatur zwischendurch**: wie oben, aber vor jedem Lochsprung werden mit
  Wahrscheinlichkeit 0,5 zufällig falsche Nachbarpaare repariert; am Ende wird der Rest repariert.
- **Test, Reparatur sofort**: vor jedem Lochsprung werden alle gerade vorhandenen benachbarten
  falschen up/dn-Paare repariert (zufällige Paarwahl); am Ende wird der Rest repariert.
- **Suche +1/-1 Reparatur**: für jeden der 100 Starts wird die Reparatur gezielt gesucht, einmal
  für Endvorzeichen +1 und einmal für −1. Die Reparatur-Vertauschungen (benachbarte falsche
  up/dn-Paare) dürfen beliebig zwischen die Lochsprünge gemischt werden. Die Suche ist erschöpfend
  (Tiefensuche über alle Reihenfolgen, Abbruch bei 20000 Zuständen, dann Hinweis). Ausgabe: wie oft
  +1, −1, beide, nur eines oder keines erreichbar war. **Beispiel +1** / **Beispiel -1** laden
  nacheinander Starts samt gefundener Schrittfolge (mit „Vor“ abspielbar).
- **Abweichung laden**: lädt eine Startkonfiguration, deren Vorzeichen von der Mehrheit abweicht.
  **Zufallsbeispiel**: lädt einen beliebigen Lauf. Mit **Vor** wird Weg plus Reparatur Schritt für
  Schritt abgespielt (Lochspur grau).

## Ergebnisse der 100er-Tests

Je Zelle 100 gültige Läufe mit zufälligen Startkonfigurationen (ein Durchlauf, ein Zufallsstand);
in Klammern die Anzahl der Läufe, die sich per Paar-Austausch nicht reparieren ließen und
verworfen wurden.

| Weg | Reparatur am Ende | Reparatur zwischendurch (50 %) | Reparatur sofort |
|---|---|---|---|
| `RULD` (Plakette) | 0× +1 / 100× −1 (0) | 60× +1 / 40× −1 (0) | 100× +1 / 0× −1 (0) |
| `RDLU` | 0× +1 / 100× −1 (0) | 61× +1 / 39× −1 (0) | 100× +1 / 0× −1 (0) |
| `RRULLD` | 100× +1 / 0× −1 (0) | 35× +1 / 65× −1 (14) | 0× +1 / 100× −1 (0) |
| `RRDDLLUU` | 0× +1 / 100× −1 (114) | 74× +1 / 26× −1 (47) | 100× +1 / 0× −1 (0) |
| `RURDDLUL` | 100× +1 / 0× −1 (71) | 75× +1 / 25× −1 (10) | 100× +1 / 0× −1 (0) |
| `RRULLDRRULLD` | 0× +1 / 100× −1 (0) | 38× +1 / 62× −1 (16) | 100× +1 / 0× −1 (0) |
| `RUULDDRU` (nicht geschlossen) | 51× +1 / 49× −1 (0) | 59× +1 / 41× −1 (17) | 49× +1 / 51× −1 (0) |

Kurzfassung:

- Geschlossene Wege (Loch kehrt zum Startplatz zurück), Reparatur nur am Ende oder immer sofort:
  je Weg immer dasselbe Vorzeichen, unabhängig von Startkonfiguration, Lochposition und Wahl der
  Reparatur-Vertauschungen. Das Vorzeichen bei „sofort“ ist dabei bei jedem getesteten Weg das
  entgegengesetzte von „am Ende“; einzige Ausnahme ist `RURDDLUL` (beide +1).
- Reparatur zwischendurch (zufälliger Zeitpunkt): das Vorzeichen ist nicht mehr einheitlich.
  Der Zeitpunkt der Reparatur ändert also das Ergebnis.
- Nicht geschlossener Weg `RUULDDRU`: bei allen Varianten beide Vorzeichen etwa gleich oft.
- Bei längeren Wegen ließen sich Läufe mit reinem Paar-Austausch am Ende teils nicht reparieren;
  diese Läufe sind verworfen und in Klammern gezählt.

### Suche nach +1/−1-Reparaturen

Je 100 Starts, Reparaturen beliebig zwischen den Lochsprüngen:

| Weg | +1 erreichbar | −1 erreichbar | beide | keines |
|---|---|---|---|---|
| `UUUR` | 48 | 52 | 0 | 0 |
| `RULD` | 100 | 100 | 100 | 0 |
| `RRULLD` | 100 | 100 | 100 | 0 |
| `RRDDLLUU` | 100 | 100 | 100 | 0 |

Bei den geschlossenen Wegen lässt sich das Vorzeichen durch die Wahl der Reparaturen also
frei bestimmen (obwohl Reparatur am Ende bzw. sofort jeweils nur ein Vorzeichen liefert).
Bei `UUUR` (Loch endet nicht am Startplatz) legt dagegen schon der Startzustand das Vorzeichen
fest: in jedem Start ist genau ein Vorzeichen erreichbar, und zwar etwa je zur Hälfte +1 bzw. −1.
Es hängt nicht allein von der Lochposition ab.

Systematischer Test (Suche nach +1/−1 mit je ca. 5–40 Starts pro Weg; alle Wege der Längen 2, 4, 6
ohne direktes Zurückgehen, für Länge 8 eine Zufallsauswahl von 138 Wegen, insgesamt 1182 Wege):
**Beide Vorzeichen sind genau dann erreichbar, wenn der Weg einen Platz zweimal besucht**
(also eine Schleife enthält; nicht nur geschlossene Wege, z.B. `RRULDD`). Bei selbstvermeidenden
Wegen, wie `UUUR`, ist in jedem Start nur ein Vorzeichen erreichbar. In keinem der 1182 Wege gab es
eine Ausnahme von dieser Regel, und es gab keine Mischfälle innerhalb eines Weges. Das ist ein
empirischer Befund, kein Beweis. Wege ungerader Länge wurden nicht getestet.

## Beweis-Skizze: wann ist das Vorzeichen eindeutig?

Setting: Weg aus m Lochsprüngen ab Loch p₀ (m gerade, Loch endet in p_m, andere Löcher bleiben),
dazwischen s Reparatur-Vertauschungen benachbarter *falscher* up/dn-Paare, am Ende Néel.

1. **Vorzeichen = Permutationsparität.** In der NetKet-Ordnung stehen alle down-Orbitale vor allen
   up-Orbitalen. Ein Sprung eines Spins σ hat daher das Vorzeichen (−1)^(Zahl der σ-Teilchen zwischen
   Start- und Zielplatz). Die Vertauschung up(a)/dn(b) ist ein up-Sprung a→b plus ein dn-Sprung b→a.
   Verfolgt man die Teilchen mit Namen, ist das Gesamtvorzeichen sgn(π↑)·sgn(π↓), wobei π_σ die
   Permutation der σ-Teilchen (in Platzreihenfolge) zwischen Anfang und Ende ist. Es hängt nur von
   der Endzuordnung der Teilchen ab, nicht vom Weg dorthin. (Numerisch geprüft: 3000 zufällige Prozesse.)
2. **Zerlegung.** Alle m+s Schritte sind Transpositionen von Platzinhalten, also hat die
   Gesamtpermutation aller 64 Inhalte (Löcher mitgezählt) das Vorzeichen (−1)^(m+s). Sie zerfällt in
   eine "stabile" Typ-Umordnung (hängt nur vom Typmuster vorher/nachher ab) mal sgn(π↑)·sgn(π↓)·sgn(π_Loch).
   Das Typmuster vorher/nachher ist Néel mit Loch in p₀ bzw. p_m. Damit gilt
   **Vorzeichen = (−1)^(m+s) · F(Start, p_m)**, F unabhängig von Weg und Reparatur.
   (Numerisch geprüft: in 3347 Gruppen mit gleichem Start und Endplatz nie ein Widerspruch.)
3. **Zähle falsche Spins.** Sei W die Zahl der Spins auf dem falschen Untergitter. Jeder Lochsprung
   ändert W um ε=±1 (+1, wenn ein richtiger Spin springt, −1, wenn ein falscher springt), jede
   Reparatur ändert W um −2, W beginnt und endet bei 0. Also s = (Σε)/2 = m/2 − n₋, mit n₋ = Zahl der
   Sprünge eines falschen Spins. Damit **Vorzeichen = (−1)^(3m/2 + n₋) · F**.
4. **Satz:** Besucht der Weg keinen Platz zweimal, ist n₋ = 0, das Vorzeichen also eindeutig.
   Beweis: Invariant "alle falschen Spins stehen auf bereits besuchten Wegplätzen" (Induktion: ein
   Sprung legt den Spin auf einen besuchten Platz, Reparaturen vertauschen nur falsche Spins).
   Der nächste Wegplatz p_{k+1} ist unbesucht, also unberührt, enthält den ursprünglichen, richtigen
   Spin, und dieser springt als richtiger: ε=+1. Also n₋=0 und s=m/2 in jeder Reparaturstrategie. □
   **Zusatz (Extra-Vertauschungen):** Ein up/dn-Paar benachbarter Plätze ist immer entweder beide
   richtig oder beide falsch, jede Vertauschung ändert W also um ±2. Mit s₊ Vertauschungen (+2) und
   s₋ (−2) folgt s = 2·s₊ + m/2 − n₋, die Parität von s hängt also nur von n₋ ab, egal wie viele
   Extra-Vertauschungen eingestreut werden. Berühren die Extra-Vertauschungen keinen Wegplatz, bleibt
   der Beweis von Schritt 4 gültig (n₋ = 0): das Vorzeichen ist weiterhin eindeutig (numerisch
   bestätigt, `UUUR`). Vertauschungen, die einen noch nicht besuchten Wegplatz berühren, können das
   Vorzeichen dagegen ändern.
5. **Umkehrung (nicht bewiesen).** Mit einem Wiederbesuch kann ein Sprung einen falschen Spin
   bewegen (n₋ ≥ 1); je nachdem ob man den Spin vorher repariert oder nicht, ändert sich die
   Parität von n₋ und damit das Vorzeichen (Beispiel `RULD`: Reparatur am Ende n₋=1, s=1;
   sofort n₋=0, s=2; Vorzeichen entgegengesetzt). Dass beide Paritäten für jeden Weg mit Wiederbesuch
   und jeden Start auch tatsächlich bis zum Néel-Zustand vervollständigbar sind, ist nur empirisch
   belegt (1182 Wege), nicht bewiesen.

## Vorzeichen eines einfachen Weges gezielt ändern

Für einen einfachen Weg (kein Platz doppelt besucht, z.B. `UUUR`) ist das Vorzeichen bei reiner
Reparatur (nur falsche up/dn-Paare tauschen) eindeutig: jeder Start liefert nur ein Vorzeichen.
Es lässt sich aber durch **zusätzliche Vertauschungen auf noch nicht besuchten Wegplätzen**
ändern:

- Vor dem Lochsprung wird ein up/dn-Paar vertauscht, das einen späteren Wegplatz enthält. Dort
  stehen dann zwei falsche Spins; wenn das Loch ankommt, springt ein falscher Spin (n₋ ≥ 1).
  Das ändert die Parität von n₋ und damit das Vorzeichen (Schritt 3 der Beweis-Skizze).
- Danach wird wie üblich bis zur Néel-Konfiguration repariert. Start, Weg und Endkonfiguration sind
  dieselben, nur die Zuordnung der Teilchen am Ende ist eine andere. Es ist also ein anderer
  Prozess (anderer Spin-Austausch, im t-J-Modell ein J-Term), kein Vorzeichenwechsel desselben Prozesses.
- Extra-Vertauschungen, die keinen Wegplatz berühren (auch wenn sie nicht zurückgenommen, sondern
  anders repariert werden), ändern das Vorzeichen nicht (Zusatz in der Beweis-Skizze: s ändert sich nur um
  gerade Zahlen). Dasselbe gilt für Vertauschungen, die einfach zurückgetauscht werden.

Numerischer Befund (`UUUR`, 6 gültige Starts, erschöpfende Suche mit bis zu 3 Extra-Vertauschungen
nahe dem Weg): mit Extra-Vertauschungen auf Wegplätzen sind in jedem Start beide Vorzeichen
erreichbar, ohne sie (bis 5 Extra-Vertauschungen) nur eines. Die Stichprobe ist klein. Das Programm
bietet dafür noch keinen eigenen Button; „Suche +1/-1 Reparatur“ erlaubt nur Reparatur-Paare.

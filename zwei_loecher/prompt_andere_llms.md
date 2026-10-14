Bitte um kritische Prüfung: Kann ein fluktuierender Spinhintergrund im t-J-Modell Paar-ODLRO erzwingen?

Ich möchte deine ehrliche, kritische Einschätzung zu einer Argumentation. Bitte stimme nicht einfach zu. Benenne
konkret, welche Schritte stimmen, welche falsch, lückenhaft oder zirkulär sind, und wo es bekannte Ergebnisse gibt,
die dafür oder dagegen sprechen.

== 1. Modell und Sichtweise ==

2D-t-J-Modell auf dem Quadratgitter, Löcher in einem antiferromagnetischen Hintergrund (Néel als Ausgangsbild).
Ein Lochsprung verschiebt ein Elektron (mit seinem Spin) um einen Platz. Ein Loch, das läuft, hinterlässt eine
Spur verschobener Spins; der J-Term (Vertauschung benachbarter entgegengesetzter Spins, S⁺S⁻) kann diese Spur
„reparieren“, aber auch neue Fluktuationen erzeugen (Vertauschungen, die keine Reparatur sind). Jeder Prozess
(Folge aus Lochsprüngen und Vertauschungen) trägt zu exp(−τH) mit einem Vorzeichen bei. Diese Vorzeichen werden
exakt aus der fermionischen Operatoralgebra berechnet (feste Orbitalordnung); geprüft ist, dass das so berechnete
Vorzeichen das Vorzeichen des Gewichts in exp(−τH) ist.

Bezeichnungen:
- O = Menge der Elektronen, die von Lochsprüngen (aller Löcher zusammen) ungerade oft bewegt wurden.
- n↑ = Zahl der Lochsprünge, die ein ↑-Elektron bewegen.
- π_Loch = Permutation der Löcher in einem geschlossenen Prozess.
- J-Schritt = zwei Lochsprünge, danach die Vertauschung genau der beiden bewegten Elektronen (Reparatur).
- Raupe = zwei benachbarte Löcher, das zweite folgt genau der Spur des ersten.
- Pivot/Domino = Paarbewegungen aus zwei Sprüngen und einer Reparatur (Paar dreht sich bzw. verschiebt sich
  seitlich über eine Nachbarbindung).

== 2. Was exakt bzw. numerisch geprüft ist ==

(R1) Für jeden geschlossenen Prozess gilt
     Vorzeichen = (−1)^s · sgn(π_Loch) = (−1)^(n↑) · sgn(π_Loch) = (−1)^(Zahl der ↑ in O) · sgn(π_Loch),
     s = Zahl der Vertauschungen. In der zweiten und dritten Form tragen Vertauschungen kein Vorzeichen
     (Marshall-Regel). Geprüft an zufälligen geschlossenen Prozessen auf beliebigen Konfigurationen.
(R2) Einzelnes Loch: Wird vor seinem Weg ein Spin auf dem Weg mit einem Nachbarn neben dem Weg vertauscht und
     danach zurückgetauscht, kehrt sich das Vorzeichen bei gleichem Anfang und Ende um. In einem fluktuierenden
     Hintergrund gibt es Lochbewegungen also mit beiden Vorzeichen.
(R3) Raupe: Jedes Elektron wird zweimal bewegt, O ist leer, ihr Vorzeichen hängt nicht vom Hintergrund ab. Allein
     ist sie aber lokalisiert: Ihre effektive Hüpfmatrix hat dasselbe Spektrum wie −A(Liniengraph des
     Quadratgitters), oberstes Band exakt flach (Grund: Wenden vertauschen die beiden Löcher, −1).
(R4) Pivot und Domino sind im sauberen Néel unfrustriert, und zwar in einer d-Wellen-Eichung. In dieser Eichung
     sind alle Paarbewegungen positiv außer den geraden Raupenschritten (alle negativ). Das Vorzeichen von
     Pivot/Domino kehrt sich um, wenn die Nachbarbindung, deren zwei Spins sie bewegen, parallel ist.
(R5) Tauschen zwei Paare die Plätze, ergibt das +1 (Paare sind untereinander bosonisch).
(R6) Ein J-Schritt mit eigener Reparatur hat immer genau einen Sprung, der ein ↑ bewegt (vertauscht werden
     kann nur ein antiparalleles Paar). Er trägt also in jedem Hintergrund das Vorzeichen −1.
(R7) Eine ↑-Lücke (Loch auf Untergitter A) und eine ↓-Lücke (auf B), die unabhängig mit solchen J-Schritten
     laufen, auch weit getrennt und mit sich kreuzenden Wegen: Das Vorzeichen ist das Produkt der
     Einzelvorzeichen (1563 von 1563 Zufallsprozessen), deckungsgleiche Wege geben immer +1. Viele Löcher mit
     J-Schritten und Raupenschritten: Vorzeichen = (−1)^(Schritte) · sgn(π_↑) · sgn(π_↓) (1500 von 1500), d. h.
     zwei Sorten Fermionen ohne Vorzeichen zwischen den Sorten.
(R8) Kleine Systeme, effektive Prozessmatrizen (bis 3 Elementarschritte, einfache Gewichte, keine Energienenner):
     Je mehr Fluktuationen nahe einem Paar erlaubt sind, desto weniger frustriert ist die Paarbewegung
     (λ_max(M)/λ_max(|M|) = 0,708 → 0,880 → 0,932 für 0, ≤1, ≤2 Fluktuationen; mit zufälligen Vorzeichen nur
     0,721). Der Grundzustand bindet Fluktuationen an das Paar. Auf einem 12×12-Torus hat das Paar mit
     Fluktuationen ein Band mit eindeutigem Maximum bei K = 0. Paare, deren Löcher nur auf gemeinsamer Spur
     laufen (O leer, Abstand entlang der Spur bis 3), bleiben praktisch lokalisiert (Beweglichkeit 2–3 % des
     vorzeichenfreien Werts).
(R9) Exakte Diagonalisierung 4×4 periodisch, J/t = 0,4: Zwei Löcher sind mit Vorzeichen gebunden
     (E_b = −0,35 t), im vorzeichenfreien Modell |H| nicht (+0,19 t). Die Vorzeichen kosten Hüpfenergie
     (−5,5 statt −7,4).

== 3. Die Argumentation (bitte prüfen) ==

Annahme A: Fluktuationen des Hintergrunds machen jede Bewegung eines einzelnen Lochs mit beiden Vorzeichen
möglich (gestützt durch R2).

Schritt 1 (Auswahl durch Interferenz): Nach R1 hängt das Vorzeichen nur von den ↑ in O ab. Unter A ist für ein
Elektron in O offen, ob es ↑ war; summiert man über die Varianten, löschen sich solche Beiträge aus. Fest bleiben
nur Beiträge, deren O leer ist oder nur aus „beglaubigten“ Paaren besteht (zwei Elektronen, die im selben Prozess
miteinander vertauscht wurden, also sicher eines ↑, eines ↓).

Schritt 2 (Löcher müssen nicht zusammenbleiben): Zwei Löcher, die einzeln (mit beglaubigten J-Schritten, R6) zur
neuen Position laufen: Die Summe über alle Wege zerfällt in G_A(R) · G_B(R). Beide Sorten haben dieselbe
Vorzeichenregel, in einem gleichförmigen Hintergrund also G_A = G_B und Paartransfer = G(R)² ≥ 0, auf jeder
Entfernung nicht negativ, egal wie stark G(R) selbst im Vorzeichen schwankt. Die Fluktuationen löschen genau die
Beiträge aus, die den beiden Partnern verschiedene Vorzeichen geben würden (Paarbrecher).

Schritt 3 (Interferenz erzeugt die Anziehung): Interferenz vernichtet kein Gewicht, sie verteilt es um (die Spur
der Zweiteilchen-Dichtematrix ist N↑ · N↓). Im Grundzustand fließt Gewicht dorthin, wo die Wege konstruktiv
interferieren. Unabhängige Bewegung verliert durch Auslöschung (Schritt 1, Frustration der Einzellochbewegung);
Bewegung, die es nur für korrelierte Löcher gibt (Raupe, gemeinsame Spur; O leer, immun), ist vollständig
konstruktiv. Also fließt Gewicht in gepaarte Konfigurationen; diese Umverteilung ist die effektive Anziehung
(passt zu R9).

Schritt 4 (ODLRO): Das Gewicht sammelt sich im konstruktivsten Paarzustand (Bandmaximum bei K = 0, R8). Paare
tauschen untereinander mit +1 (R5), viele können ihn gemeinsam besetzen: makroskopische Besetzung (nach Yangs
Schranke ∝ N), also ⟨Δ†(R)Δ(0)⟩ → konstant ≠ 0. ODLRO im Einteilchenkanal ist für Elektronen nach Pauli/Yang
ausgeschlossen, deshalb erscheint es im Paarkanal.

== 4. Bekannte Schwachstellen ==

- Innerhalb einer Sorte sind die Löcher Fermionen (R7); bei endlicher Dichte gilt Pauli. Dass die Umverteilung
  aus Schritt 3 tatsächlich einen festen Anteil des Gewichts in einen Paarzustand bündelt, ist nicht gezeigt.
- Schritt 1 ist heuristisch (Gewichte der Varianten, Quanten-Hintergrund statt fester Konfigurationen).
- Rest-Frustration zwischen geraden Raupenschritten und J-Wenden; eine beglaubigte Variante des geraden
  Raupenschritts mit passendem Vorzeichen gibt es nur mit Fluktuationen, ihr Gewicht ist offen.
- Die Numerik ist klein (4×4 ED, 6×6-Fenster, ein Paar auf 12×12, höchstens zwei Fluktuationen, einfache
  Gewichte statt Energien).

== 5. Fragen ==

1. Welche Schritte sind korrekt, welche falsch, lückenhaft oder zirkulär? Bitte konkret.
2. Ist Schritt 3 („Interferenz verteilt um, das ist die Anziehung“) tragfähig? Wie verhält er sich zu bekannten
   Mechanismen kinetischer Paarbildung?
3. Wie verhält sich das zur Literatur (Phase-String-Theorie und gegenseitige Chern-Simons-Statistik nach Weng,
   Vorzeichenstruktur nach Wu–Weng–Zaanen, String-Bild nach Grusdt/Bohrdt, RVB nach Anderson, Trugman-Schleifen)?
   Ist irgendetwas daran neu?
4. Gibt es ein Gegenbeispiel: Systeme mit starken Spinfluktuationen, in denen Einzellochbewegung beide Vorzeichen
   hat, aber keine Supraleitung entsteht? Woran scheitert die Argumentation dort?
5. Welche konkrete Rechnung (Größe, Methode, Observable) würde die Argumentation entscheiden?

Details und Skripte (falls zugänglich): https://github.com/dsmic/tJ-model-visualisierung/tree/nachdenken/zwei_loecher
(README, Abschnitt 7).

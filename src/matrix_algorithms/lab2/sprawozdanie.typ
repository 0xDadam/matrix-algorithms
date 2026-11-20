#set document(
  title: "Rekurencyjne algorytmy macierzowe",
  author: "Zespół projektowy",
  date: datetime.today(),
)
#set page(paper: "a4", margin: 2cm)
#set text(font: "Libertinus Serif", lang: "pl", size: 10.5pt)
#set heading(numbering: "1.1")
#set par(justify: true, leading: 0.65em)

#align(center)[
  #text(size: 18pt, weight: "bold")[Rekurencyjne algorytmy macierzowe]
  #v(0.7cm)
  #text(size: 13pt)[Odwracanie macierzy, faktoryzacja LU i mnożenie]
  #v(1.2cm)
  Laboratorium z algorytmów macierzowych
  #v(0.3cm)
  #datetime.today().display("[day].[month].[year]")
]

= Cel i zakres ćwiczenia

Celem zadania było zaimplementowanie rekurencyjnych algorytmów macierzowych:

- odwracania macierzy metodą Bineta,
- odwracania macierzy metodą Strassena,
- faktoryzacji $L U$ z użyciem mnożenia Bineta i Strassena,
- odwracania macierzy metodą Ladermana,
- oraz, opcjonalnie, rekurencyjnego mnożenia metodą Ladermana.

Macierze wejściowe są kwadratowe i mogą mieć dowolny rozmiar. Dla rozmiarów
niebędących potęgą podstawy algorytmu stosowane jest dopełnienie macierzy:
dla Bineta i Strassena do najbliższej potęgi $2$, a dla Ladermana do
najbliższej potęgi $3$. Elementy dopełnienia są dobierane tak, aby nie
zmienić lewego górnego bloku wyniku.

Wszystkie operacje arytmetyczne wykonywane są na obiekcie `CounterFloat`.
Dzięki temu osobno zliczane są dodawania, odejmowania, mnożenia i dzielenia,
a liczba FLOP jest ich sumą.

#pagebreak()
= Opis algorytmów

== Rekurencyjne odwracanie przez dopełnienie Schura

Macierz dzielimy na cztery bloki:

$ A = mat( A_11, A_12; A_21, A_22 ). $

Najpierw rekurencyjnie obliczamy $P = A_11^(-1)$, a następnie

$ S = A_22 - A_21 P A_12. $

Po obliczeniu $S^(-1)$ wynik ma postać

$ A^(-1) =
  mat(
    P + P A_12 S^(-1) A_21 P, & -P A_12 S^(-1);
    -S^(-1) A_21 P, & S^(-1)
  ). $

W implementacji ta sama procedura jest używana dla Bineta i Strassena;
zmienia się wyłącznie funkcja mnożąca bloki. Przypadkiem bazowym jest
macierz $1 times 1$.

Pseudokod:

```text
inverse(A, multiply):
    B := dopełnij A macierzą jednostkową do rozmiaru 2^k
    inverse_recursive(B)
    zwróć lewy górny blok rozmiaru A

inverse_recursive(A):
    jeśli rozmiar(A) = 1:
        zwróć 1 / A[1,1]
    podziel A na A11, A12, A21, A22
    P := inverse_recursive(A11)
    p_b := multiply(P, A12)
    c_p := multiply(A21, P)
    S := A22 - multiply(c_p, A12)
    Q := inverse_recursive(S)
    zwróć blokową macierz:
        P + multiply(multiply(p_b,Q),c_p),
        -multiply(p_b,Q)
        -multiply(Q,c_p), Q
```

== Rekurencyjna faktoryzacja LU

Dla tego samego podziału wyznaczane są kolejno:

$ L_11 U_11 = A_11, quad U_12 = L_11^(-1) A_12, quad
L_21 = A_21 U_11^(-1), $

$ S = A_22 - L_21 U_12, quad L_22 U_22 = S. $

Wynikiem są macierze blokowe
$ L = mat(L_11, 0; L_21, L_22) $ oraz
$ U = mat(U_11, U_12; 0, U_22). $
Rozwiązywanie układów trójkątnych odbywa się iteracyjnie, natomiast
mnożenie aktualizacji Schura wykorzystuje wybrany algorytm rekurencyjny.

#pagebreak()
== Algorytm Ladermana

Laderman dzieli macierz na dziewięć bloków i wykonuje 23 rekurencyjne
mnożenia zamiast 27 mnożeń klasycznych. Dla macierzy $3 times 3$ wzór
rekurencyjny ma postać

$ T(n) = 23 T(n/3) + O(n^2), $

czyli $T(n) = O(n^(log_3 23)) approx O(n^"2.855")$.
Odwracanie metodą Ladermana korzysta z tego mnożenia blokowego i z
dwublokowego wzoru Schura dla pozostałej części macierzy.

== Złożoność z uwzględnieniem paddingu

Niech

$ p_2(n) = 2^ceil(log_2 n) quad "oraz" quad
  p_3(n) = 3^ceil(log_3 n) $

oznaczają rozmiary macierzy po dopełnieniu. To właśnie $p_2(n)$ lub
$p_3(n)$, a nie bezpośrednio $n$, jest rozmiarem przekazywanym do
rekurencji. Dopełnienie wymaga $Theta(p_b(n)^2)$ pamięci i pracy
organizacyjnej, natomiast zasadnicze operacje arytmetyczne mają koszt:

#table(
  columns: (1.3fr, 2fr, 2.4fr),
  inset: 5pt,
  stroke: .5pt,
  table.header([Metoda], [Rozmiar po paddingu], [Koszt mnożenia]),
  [Binet], [$p_2(n)$], [$Theta(p_2(n)^3)$],
  [Strassen], [$p_2(n)$], [$Theta(p_2(n)^(log_2 7))$],
  [Laderman], [$p_3(n)$], [$Theta(p_3(n)^(log_3 23))$],
)

Dla odwracania Bineta i Strassena rekurencja zawiera dwa odwracania
bloków oraz stałą liczbę mnożeń blokowych. Jej koszt jest więc taki jak
koszt dominującego mnożenia: odpowiednio
$Theta(p_2(n)^3)$ i $Theta(p_2(n)^(log_2 7))$. Odwracanie Ladermana ma
koszt $Theta(p_3(n)^(log_3 23))$.

W faktoryzacji LU występują dodatkowo iteracyjne rozwiązania układów
trójkątnych o koszcie $Theta(p_2(n)^3)$. Dlatego dla LU Bineta koszt
wynosi $Theta(p_2(n)^3)$, a dla LU Strassena również
$Theta(p_2(n)^3)$, mimo że samo mnożenie Strassena ma lepszy wykładnik.
Ponieważ $n <= p_b(n) < b n$, padding nie zmienia klasy asymptotycznej
względem $n$, ale powoduje schodkowy przebieg czasu i liczby FLOP na
wykresach.

#pagebreak()
= Implementacja

Najważniejszy wspólny fragment odwracania:

```python
def _inverse_recursive(matrix, multiply):
    size = matrix.shape[0]
    if size == 1:
        if float(matrix[0, 0]) == 0.0:
            raise ValueError("Macierz jest osobliwa: zerowy pivot.")
        return np.array([[CounterFloat(1.0) / matrix[0, 0]]],
                        dtype=object)

    half = size // 2
    a11 = matrix[:half, :half]
    a12 = matrix[:half, half:]
    a21 = matrix[half:, :half]
    a22 = matrix[half:, half:]

    p = _inverse_recursive(a11, multiply)
    p_b = multiply(p, a12)
    c_p = multiply(a21, p)
    schur = a22 - multiply(c_p, a12)
    schur_inverse = _inverse_recursive(schur, multiply)

    top_right = -1.0 * multiply(p_b, schur_inverse)
    bottom_left = -1.0 * multiply(schur_inverse, c_p)
    correction = multiply(multiply(p_b, schur_inverse), c_p)
    top_left = p + correction
    return np.block([[top_left, top_right],
                     [bottom_left, schur_inverse]])
```

Zliczanie operacji jest realizowane bez zmiany kodu algorytmów:

```python
def operation_counts(cls):
    return {
        "additions": cls.additions,
        "subtractions": cls.subtractions,
        "multiplications": cls.multiplications,
        "divisions": cls.divisions,
        "total": cls.additions + cls.subtractions
                 + cls.multiplications + cls.divisions,
    }
```

W przypadku wejścia o rozmiarze niebędącym potęgą podstawy wynik jest
przycinany do oryginalnego rozmiaru. Dopełnienie macierzą jednostkową dla
odwracania zachowuje blok $A^(-1)$, ponieważ
$ mat(A, 0; 0, I)^(-1) = mat(A^(-1), 0; 0, I) $.

#pagebreak()
= Środowisko eksperymentalne

Kod źródłowy zapisano w Pythonie 3.12. Python jest językiem
interpretowanym, dlatego nie spełnia literalnie wymagania zadania o języku
nieinterpretowanym. Do uruchomienia benchmarku przygotowano jednak
wykonywalny program `benchmark.exe` zbudowany za pomocą Nuitka; pomiary
wykonano z tego programu. Jest to istotne ograniczenie metodologiczne:
algorytmy pozostają zaimplementowane w Pythonie, a kompilacja zmniejsza
narzut interpretacji, lecz nie zmienia języka źródłowego.

Pomiary wykonano na:

#table(
  columns: (1fr, 2fr),
  inset: 5pt,
  stroke: .5pt,
  table.header([Parametr], [Wartość]),
  [Procesor], [Intel(R) Core(TM) Ultra 5 250KF Plus],
  [Rdzenie / wątki], [18 / 18],
  [Maksymalne taktowanie], [4.2 GHz],
  [Pamięć RAM], [około 16 GiB],
  [System], [Windows],
)

Macierze losowe generowano ze stałym ziarnem, z wartościami praktycznie
z przedziału $(0.00000001, 1.0)$ (`random.uniform` teoretycznie może
zwrócić także końce przedziału). Dla każdego rozmiaru mierzono czas
`time.perf_counter()`, a licznik `CounterFloat` zerowano przed pomiarem.

= Wyniki

Poniższe wykresy dotyczą rekurencyjnego mnożenia z pierwszego zestawu,
które jest podstawową operacją wykorzystywaną we wszystkich algorytmach
z tego sprawozdania. Zmierzono rozmiary będące potęgami $2$ lub $3$:
$2, 3, 4, 8, 9, 16, 27, 32, 64, 81, 128, 243, 256, 512$.
Pojedynczy pomiar miał limit 600 s. Binet i Strassen zakończyły pomiary
na $n=512$, natomiast Laderman przekroczył limit przy $n=256$.
Punkt timeout jest oznaczony na wykresie osobno i nie jest traktowany jako
poprawny czas wykonania. Linie ciągłe łączą wyłącznie rzeczywiste punkty
pomiarowe.

#figure(
  image("czas.svg", width: 100%),
  caption: [Czas rekurencyjnego mnożenia z punktami rzeczywistych pomiarów. Skala jest liniowa.],
) <czas>

#figure(
  image("flop.svg", width: 100%),
  caption: [Liczba operacji zmiennoprzecinkowych dla udanych pomiarów. Skala jest liniowa.],
) <flop>

#figure(
  image("czas_log.svg", width: 100%),
  caption: [Czas rekurencyjnego mnożenia w skali logarytmicznej dla rzeczywistych pomiarów.],
) <czas-log>

#figure(
  image("flop_log.svg", width: 100%),
  caption: [Liczba operacji zmiennoprzecinkowych w skali logarytmicznej dla udanych pomiarów.],
) <flop-log>

Wybrane wartości z pomiaru:

#table(
  columns: (1fr, 1fr, 1fr, 1fr, 1fr),
  inset: 4pt,
  stroke: .5pt,
  table.header([Algorytm], [$n$], [Czas [s]], [FLOP], [Mnożenia]),
  [Binet], [4], [0.000291], [112], [64],
  [Strassen], [4], [0.000411], [247], [49],
  [Laderman], [9], [0.003624], [3\,665], [529],
  [Binet], [64], [0.697908], [520\,192], [262\,144],
  [Strassen], [64], [0.740104], [798\,967], [117\,649],
  [Laderman], [64], [2.010205], [2\,192\,801], [279\,841],
  [Binet], [128], [5.905380], [4\,177\,920], [2\,097\,152],
  [Strassen], [128], [4.504200], [5\,666\,497], [823\,543],
  [Laderman], [128], [38.867784], [51\,077\,401], [6\,436\,343],
  [Binet], [243], [37.927962], [33\,488\,896], [16\,777\,216],
  [Strassen], [243], [32.087173], [39\,960\,391], [5\,764\,801],
  [Laderman], [243], [38.997924], [51\,077\,401], [6\,436\,343],
  [Binet], [512], [400.942317], [268\,173\,312], [134\,217\,728],
  [Strassen], [512], [307.714916], [280\,902\,385], [40\,353\,607],
)

Teoria zgadza się z pomiarem liczby mnożeń. Dla Bineta jest to
$p_2(n)^3$, dla Strassena $p_2(n)^(log_2 7)$, a dla Ladermana
$p_3(n)^(log_3 23)$ dla rozmiarów będących odpowiednią potęgą. Dla pozostałych rozmiarów koszt obejmuje również dopełnienie do
najbliższej potęgi, co wyjaśnia odcinki o prawie stałym koszcie.
W nowych pomiarach widać to bezpośrednio: wartości dla $n=243$ i $n=256$
mają ten sam rozmiar po paddingu dla Bineta i Strassena, a wartości dla
$n=128$ i $n=243$ mają analogiczny wspólny koszt po paddingu dla
Ladermana. Dla Ladermana przy $n=256$ odnotowano timeout po 600 s.

= Sprawdzenie poprawności

== Dowód indukcyjny

Teza: procedura `inverse_recursive` zwraca odwrotność macierzy
nieosobliwej, jeśli wszystkie rekurencyjne pivoty są niezerowe.

Dla rozmiaru $1$ algorytm zwraca $[1/a]$, więc teza zachodzi. Załóżmy,
że teza zachodzi dla bloków mniejszych niż $n$. Algorytm poprawnie wyznacza
$P=A_11^(-1)$ oraz $Q=S^(-1)$ na mocy założenia indukcyjnego. Podstawienie
tych wartości do wzoru na odwrotność blokową daje macierz, której iloczyn
z $A$ jest macierzą jednostkową. Zatem teza zachodzi dla $n$.
Dopełnienie macierzą jednostkową nie zmienia lewego górnego bloku
odwrotności, więc dowód obejmuje także dowolny rozmiar wejścia.

Dla LU podstawą indukcji jest rozkład $[a]=[1][a]$. Przy założeniu
poprawności rozkładów $A_11$ i $S$, wzory na $U_12$, $L_21$ i aktualizację
Schura zapewniają równość $L U=A$. To kończy dowód indukcyjny.

== Przykład ręczny

Dla
$ A = mat(2, 1; 1, 1) $
otrzymujemy $A_11=[2]$, $P=[1/2]$ oraz

$ S = [1] - [1][1/2][1] = [1/2], quad S^(-1)=[2]. $

Wzór blokowy daje

$ A^(-1) = mat(1, -1; -1, 2). $

Rzeczywiście:
$ mat(2,1;1,1) mat(1,-1;-1,2) = mat(1,0;0,1). $

Kod zwraca tę samą macierz, z dokładnością arytmetyki zmiennoprzecinkowej.
W teście losowej macierzy $3 times 3$ dla wszystkich trzech metod
sprawdzenie `A * A_inv approx I` zakończyło się wynikiem `True`.
Dodatkowo dla faktoryzacji LU rozmiaru $3$ otrzymano
$ max_(i,j) |(L U-A)_(i,j)| = 0.0 $ po konwersji do liczb
zmiennoprzecinkowych.

= Wnioski

Zaimplementowano wszystkie wymagane rekurencyjne procedury. Wspólny wzór
Schura pozwolił użyć tych samych testów dla Bineta i Strassena, a różnice
w liczbie operacji wynikają wyłącznie z użytego mnożenia bloków. Strassen
wykonuje mniej mnożeń, ale dodatkowe dodawania i odejmowania powodują, że
dla małych rozmiarów nie musi być szybszy od Bineta. Strassen ma najlepszy
wykładnik asymptotyczny spośród trzech metod, natomiast Laderman nadal
poprawia wykładnik klasycznego mnożenia. Jego 23 mnożenia i duża liczba
operacji pomocniczych powodują jednak znaczny koszt stały.

Wyniki potwierdzają poprawność implementacji dla testów jednostkowych,
przykładu ręcznego i losowej macierzy $3 times 3$. Aby w pełni spełnić
wymóg osobnych wykresów dla czasu odwracania i faktoryzacji LU, należałoby
wykonać dodatkową serię pomiarów tych dwóch procedur; dostarczony plik
CSV zawiera pomiary rekurencyjnego mnożenia, dlatego wykresy w tym raporcie
opisują właśnie ten zakres.

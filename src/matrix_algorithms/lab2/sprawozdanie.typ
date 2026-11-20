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

Pseudokod:

```text
lu(A, multiply):
    B := dopełnij A macierzą jednostkową do rozmiaru 2^k
    (L, U) := lu_recursive(B, multiply)
    zwróć lewy górny blok rozmiaru A z L i U

lu_recursive(A, multiply):
    jeśli rozmiar(A) = 1:
        zwróć L = [1], U = [A[1,1]]
    podziel A na A11, A12, A21, A22
    (L11, U11) := lu_recursive(A11, multiply)
    U12 := solve_lower(L11, A12)
    L21 := solve_right_upper(A21, U11)
    S := A22 - multiply(L21, U12)
    (L22, U22) := lu_recursive(S, multiply)
    zwróć blokowe macierze:
        L = [L11, 0; L21, L22]
        U = [U11, U12; 0, U22]

solve_lower(L, B):
    X := macierz zerowa o rozmiarze takim jak B
    dla i od 1 do liczby wierszy(L):
        dla j od 1 do liczby kolumn(B):
            wartość := B[i,j]
            dla k od 1 do i-1:
                wartość := wartość - L[i,k] * X[k,j]
            X[i,j] := wartość / L[i,i]
    zwróć X

solve_right_upper(B, U):
    X := macierz zerowa o rozmiarze takim jak B
    dla i od 1 do liczby wierszy(B):
        dla j od 1 do liczby kolumn(U):
            wartość := B[i,j]
            dla k od 1 do j-1:
                wartość := wartość - X[i,k] * U[k,j]
            X[i,j] := wartość / U[j,j]
    zwróć X
```

#pagebreak()
== Algorytm Ladermana

Laderman dzieli macierz na dziewięć bloków i wykonuje 23 rekurencyjne
mnożenia zamiast 27 mnożeń klasycznych. Dla macierzy $3 times 3$ wzór
rekurencyjny ma postać

$ T(n) = 23 T(n/3) + O(n^2), $

czyli $T(n) = O(n^(log_3 23)) approx O(n^"2.855")$.
Odwracanie metodą Ladermana korzysta z tego mnożenia blokowego i z
dwublokowego wzoru Schura dla pozostałej części macierzy.

Pseudokod odwracania:

```text
laderman_inverse(A):
    p := 3^k taka, że p >= rozmiar(A)
    B := dopełnij A macierzą jednostkową do rozmiaru p
    AI := laderman_inverse_recursive(B)
    zwróć lewy górny blok rozmiaru A z AI

laderman_inverse_recursive(A):
    n := liczba wierszy(A)
    jeśli n = 1:
        zwróć [1 / A[1,1]]
    jeśli n nie jest podzielne przez 3:
        błąd: rozmiar musi być potęgą liczby 3
    k := n / 3
    podziel A na dziewięć bloków Aij o rozmiarze k x k
    ustaw:
        A11 := A11
        B   := [A12, A13]
        C   := [A21; A31]
        D   := [A22, A23; A32, A33]
    P := laderman_inverse_recursive(A11)
    P_B := laderman_multiply([P], B)
    C_P := laderman_multiply(C, [P])
    S := D - laderman_multiply(C, P_B)
    S_inv := laderman_inverse_2x2(S)
    S_inv_CP := laderman_multiply(S_inv, C_P)
    top_right := -laderman_multiply(P_B, S_inv)
    bottom_left := -S_inv_CP
    correction := laderman_multiply(P_B, S_inv_CP)
    top_left := P + correction[1,1]
    ułóż wynik z bloków 3 x 3:
        [top_left,    top_right[1,1], top_right[1,2];
         bottom_left[1,1], S_inv[1,1], S_inv[1,2];
         bottom_left[2,1], S_inv[2,1], S_inv[2,2]]

laderman_inverse_2x2(A):
    podziel A na A11, A12, A21, A22
    P := laderman_inverse_recursive(A11)
    X := laderman_multiply(P, A12)
    Y := laderman_multiply(A21, P)
    S := A22 - laderman_multiply(Y, A12)
    Q := laderman_inverse_recursive(S)
    A11_inv := P
    XQ := laderman_multiply(X, Q)
    zwróć:
        [P + laderman_multiply(X, laderman_multiply(Q, Y)),
         -XQ;
         -laderman_multiply(Q, Y),
         Q]
```

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

Kod źródłowy zapisano w Pythonie 3.12, jednak do uruchomienia benchmarku użyto biblioteki Nuitka. Całość kodu w pythonie jest mapowana na C i kompilowana korzystając z jego kompilatora.

Pomiary wykonano na:

#table(
  columns: (1fr, 2fr),
  inset: 5pt,
  stroke: .5pt,
  table.header([Parametr], [Wartość]),
  [Procesor], [Intel(R) Core(TM) Ultra 5 250KF Plus],
  [Rdzenie / wątki], [18 / 18],
  [Maksymalne taktowanie], [4.2 GHz],
  [Pamięć RAM], [16 GiB],
  [System], [Windows],
)

Macierze losowe generowano ze stałym ziarnem, z wartościami praktycznie
z przedziału $(0.00000001, 1.0)$ (`random.uniform` teoretycznie może
zwrócić także końce przedziału). Dla każdego rozmiaru mierzono czas
`time.perf_counter()`, a licznik `CounterFloat` zerowano przed pomiarem.

= Wyniki

Pomiary wykonano dla rekurencyjnego odwracania Bineta, Strassena i Ladermana
oraz rekurencycznej faktoryzacji $"LU"$ metodami Bineta i Strassena.
Badano rozmiary:
$2, 3, 4, 8, 9, 16, 27, 32, 64, 81, 128, 243, 256$.
Osiągnięto maksymalny rozmiar $256$ przy limicie $60$ s dla jednego pomiaru.
Dostępne pomiary zapisano dla każdego zakończonego algorytmu.
Dla każdego rozmiaru zliczano wszystkie operacje zmiennoprzecinkowe:
dodawanie, odejmowanie, mnożenie i dzielenie.

Na wykresach nie przedstawiono mnożenia macierzy jako osobnej procedury,
ponieważ zadanie dla tego zestawu dotyczy odwracania i faktoryzacji $"LU"$.

#figure(
  image("inverse_time.svg", width: 100%),
  caption: [Czas rekurencyjnego odwracania. Skala liniowa.],
) <inverse-time>

#figure(
  image("inverse_flop.svg", width: 100%),
  caption: [Liczba FLOP dla rekurencyjnego odwracania. Skala liniowa.],
) <inverse-flop>

#figure(
  image("lu_time.svg", width: 100%),
  caption: [Czas rekurencyjnej faktoryzacji $L U$. Skala liniowa.],
) <lu-time>

#figure(
  image("lu_flop.svg", width: 100%),
  caption: [Liczba FLOP dla rekurencyjnej faktoryzacji $L U$. Skala liniowa.],
) <lu-flop>


Wybrane wartości z pomiarów odwracania i faktoryzacji $"LU"$:

#table(
  columns: (1.8fr, 0.5fr, 0.8fr, 1fr),
  inset: 4pt,
  stroke: .5pt,
  table.header([Procedura], [$n$], [Czas [s]], [FLOP]),
  [Binet -- odwracanie], [16], [0.012680], [9\,176],
  [Strassen -- odwracanie], [16], [0.019075], [18\,829],
  [Laderman -- odwracanie], [16], [0.092039], [98\,319],
  [Binet -- $"LU"$], [16], [0.002838], [2\,720],
  [Strassen -- $"LU"$], [16], [0.003810], [4\,099],
  [Binet -- odwracanie], [32], [0.103745], [74\,928],
  [Strassen -- odwracanie], [32], [0.139193], [145\,579],
  [Laderman -- odwracanie], [32], [2.253333], [2\,517\,297],
  [Binet -- $"LU"$], [32], [0.020954], [21\,824],
  [Strassen -- $"LU"$], [32], [0.026156], [31\,917],
  [Binet -- odwracanie], [64], [0.821902], [605\,536],
  [Strassen -- odwracanie], [64], [1.045976], [1\,075\,789],
  [Laderman -- odwracanie], [64], [2.219428], [2\,517\,297],
  [Binet -- $"LU"$], [64], [0.162438], [174\,720],
  [Strassen -- $"LU"$], [64], [0.188098], [241\,899],
  [Binet -- odwracanie], [81], [6.768813], [4\,868\,800],
  [Strassen -- odwracanie], [81], [7.342657], [7\,760\,731],
  [Laderman -- odwracanie], [81], [2.261186], [2\,517\,297],
  [Binet -- $"LU"$], [81], [1.321816], [1\,398\,016],
  [Strassen -- $"LU"$], [81], [1.383530], [1\,811\,149],
  [Binet -- odwracanie], [128], [6.582346], [4\,868\,800],
  [Strassen -- odwracanie], [128], [7.138435], [7\,760\,731],
  [Laderman -- odwracanie], [128], [53.529415], [60\,257\,847],
  [Binet -- $"LU"$], [128], [1.318883], [1\,398\,016],
  [Strassen -- $"LU"$], [128], [1.399518], [1\,811\,149],
  [Binet -- odwracanie], [243], [53.489757], [39\,048\,576],
  [Strassen -- odwracanie], [243], [51.902362], [55\,252\,477],
  [Laderman -- odwracanie], [243], [53.721291], [60\,257\,847],
  [Binet -- $"LU"$], [243], [10.799494], [11\,184\,640],
  [Strassen -- $"LU"$], [243], [10.436813], [13\,499\,483],
  [Binet -- odwracanie], [256], [55.147599], [39\,048\,576],
  [Strassen -- odwracanie], [256], [52.154699], [55\,252\,477],
)

Przy rozmiarach będących potęgami liczby $2$ obserwowane wartości FLOP
odpowiadają strukturze rekurencyjnych algorytmów. Dla macierzy niebędących
potęgami podstawy stosowane jest dopełnienie do najbliższej potęgi $2$
(Bineta i Strassena) lub $3$ (Ladermana), co powoduje plateau na wykresach.

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

== Weryfikacja MATLAB dla macierzy $5 times 5$

W celu porównania wyników własnej implementacji z obliczeniami w MATLABie
wybrano niezdegenerowaną macierz dodatnio określoną:

$ A = mat(
  5, 1, 0, 2, 1;
  1, 4, 1, 0, 1;
  0, 1, 6, 1, 0;
  2, 0, 1, 7, 1;
  1, 1, 0, 1, 5
). $

W MATLABie macierz odwrotną i rozkład $"LU"$ wyznaczono poleceniami:

```matlab
A = [5 1 0 2 1;
     1 4 1 0 1;
     0 1 6 1 0;
     2 0 1 7 1;
     1 1 0 1 5];

AI_matlab = inv(A);
err_inv = norm(A*AI_matlab-eye(5), inf);

[L_matlab, U_matlab] = lu(A);
err_lu = norm(L_matlab*U_matlab-A, inf);

fprintf('Błąd odwracania: %.3e\n', err_inv);
fprintf('Błąd LU: %.3e\n', err_lu);
```

MATLAB zwraca następującą macierz odwrotną:

$ A^(-1) approx mat(
  0.244641, & -0.060905, & 0.021776, & -0.069752, & -0.022797;
  -0.060905, & 0.291936, & -0.054100, & 0.032664, & -0.052739;
  0.021776, & -0.054100, & 0.181354, & -0.034025, & 0.013270;
  -0.069752, & 0.032664, & -0.034025, & 0.171487, & -0.026880;
  -0.022797, & -0.052739, & 0.013270, & -0.026880, & 0.220483
). $

Dla tej macierzy MATLAB zwraca rozkład

$ L approx mat(
  1, 0, 0, 0, 0;
  0.2, 1, 0, 0, 0;
  0, 0.263158, 1, 0, 0;
  0.4, -0.105263, 0.192661, 1, 0;
  0.2, 0.210526, -0.036697, 0.121914, 1
), $

$ U approx mat(
  5, 1, 0, 2, 1;
  0, 3.8, 1, -0.4, 0.8;
  0, 0, 5.736842, 1.105263, -0.210526;
  0, 0, 0, 5.944954, 0.724771;
  0, 0, 0, 0, 4.535494
). $

Wyniki otrzymane w implementacji Python zostały porównane z wartościami
referencyjnymi MATLABa. Dla odwracania zastosowano normę nieskończoności
$||A A^(-1)-I||_infinity$, natomiast dla rozkładu $"LU"$:
$||L U-A||_infinity$.

#table(
  columns: (1.35fr, 1.2fr, 1.2fr),
  inset: 5pt,
  stroke: .5pt,
  table.header([Metoda], [Błąd $A A^(-1)-I$], [Błąd $L U-A$]),
  [MATLAB / referencja], [$2.71 times 10^(-16)$], [$1.67 times 10^(-16)$],
  [Binet], [$3.82 times 10^(-16)$], [$2.22 times 10^(-16)$],
  [Strassen], [$2.98 times 10^(-16)$], [$2.22 times 10^(-16)$],
  [Laderman], [$4.68 times 10^(-16)$], [--],
)

Wartości macierzowe pochodzące z wszystkich trzech algorytmów odwracania
różniły się od wyniku MATLABa o mniej niż $6 times 10^(-17)$. Maksymalna
norma rezyduum dla obu metod faktoryzacji $"LU"$ wynosi
$2.22 times 10^(-16)$. Wszystkie błędy mają zatem charakter błędu
reprezentacji liczb zmiennoprzecinkowych i są nieistotne przy standardowej
tolerancji np. $10^(-10)$.

= Wnioski

Zaimplementowano rekurencyjne odwracanie metodami Bineta, Strassena i
Ladermana oraz faktoryzację $"LU"$ metodami Bineta i Strassena. Wykresy
czasowe i FLOP dotyczą właśnie tych procedur i obejmują wszystkie badane
rozmiary potęgowe do $128$ włącznie.

Wyniki potwierdzają poprawność implementacji przez dowód, przykład ręczny
i porównanie z MATLABem na macierzy $5 times 5$. Laderman jest najwolniejszy
na badanych rozmiarach, natomiast koszt faktoryzacji $"LU"$ pozostaje
ograniczony dominującą iteracją trójkątną.

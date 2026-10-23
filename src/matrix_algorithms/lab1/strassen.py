import argparse
import time
import numpy as np

from src.matrix_algorithms.core import Matrix, pad_to_power_of_two, random_matrix
from src.matrix_algorithms.counterfloat import CounterFloat


def strassen_multiplication(A: Matrix, B: Matrix) -> Matrix:
    orig_n, _ = A.shape
    _, orig_p = B.shape

    padded_A, padded_B, _ = pad_to_power_of_two(A, B)
    res_padded = _strassen_recursive(padded_A, padded_B)

    return res_padded[:orig_n, :orig_p]


def _strassen_recursive(A: Matrix, B: Matrix) -> Matrix:
    n = A.shape[0]

    if n == 1:
        return A * B

    mid = n // 2

    A11, A12 = A[:mid, :mid], A[:mid, mid:]
    A21, A22 = A[mid:, :mid], A[mid:, mid:]

    B11, B12 = B[:mid, :mid], B[:mid, mid:]
    B21, B22 = B[mid:, :mid], B[mid:, mid:]

    M1 = _strassen_recursive(A11 + A22, B11 + B22)
    M2 = _strassen_recursive(A21 + A22, B11)
    M3 = _strassen_recursive(A11, B12 - B22)
    M4 = _strassen_recursive(A22, B21 - B11)
    M5 = _strassen_recursive(A11 + A12, B22)
    M6 = _strassen_recursive(A21 - A11, B11 + B12)
    M7 = _strassen_recursive(A12 - A22, B21 + B22)

    C11 = M1 + M4 - M5 + M7
    C12 = M3 + M5
    C21 = M2 + M4
    C22 = M1 - M2 + M3 + M6

    return np.block([[C11, C12], [C21, C22]])


def main() -> None:
    parser = argparse.ArgumentParser(description="Mnożenie macierzy metodą Strassena")
    parser.add_argument("size", type=int, help="Rozmiar kwadratowych macierzy N x N")
    parser.add_argument("seed", type=int, nargs="?", default=None, help="Seed dla generatora")
    args = parser.parse_args()

    A = random_matrix(args.size, args.seed)
    seed_B = args.seed + 1 if args.seed is not None else None
    B = random_matrix(args.size, seed_B)

    CounterFloat.reset_counter()
    start_time = time.perf_counter()
    _ = strassen_multiplication(A, B)
    elapsed = time.perf_counter() - start_time

    print(f"=== Strassen ({args.size}x{args.size}) ===")
    print(f"Czas: {elapsed:.6f} s")
    print(f"FLOP: {CounterFloat.operation_counts()}")


if __name__ == "__main__":
    main()
import argparse
import time
import numpy as np

from src.matrix_algorithms.core import Matrix, pad_to_power_of_two, random_matrix
from src.matrix_algorithms.counterfloat import CounterFloat


def binet_multiplication(A: Matrix, B: Matrix) -> Matrix:
    orig_n, _ = A.shape
    _, orig_p = B.shape

    padded_A, padded_B, _ = pad_to_power_of_two(A, B)
    res_padded = _binet_recursive(padded_A, padded_B)

    return res_padded[:orig_n, :orig_p]


def _binet_recursive(A: Matrix, B: Matrix) -> Matrix:
    n = A.shape[0]

    if n == 1:
        return A * B

    mid = n // 2

    A11, A12 = A[:mid, :mid], A[:mid, mid:]
    A21, A22 = A[mid:, :mid], A[mid:, mid:]

    B11, B12 = B[:mid, :mid], B[:mid, mid:]
    B21, B22 = B[mid:, :mid], B[mid:, mid:]

    P1 = _binet_recursive(A11, B11)
    P2 = _binet_recursive(A12, B21)
    P3 = _binet_recursive(A11, B12)
    P4 = _binet_recursive(A12, B22)
    P5 = _binet_recursive(A21, B11)
    P6 = _binet_recursive(A22, B21)
    P7 = _binet_recursive(A21, B12)
    P8 = _binet_recursive(A22, B22)

    C11 = P1 + P2
    C12 = P3 + P4
    C21 = P5 + P6
    C22 = P7 + P8

    return np.block([[C11, C12], [C21, C22]])


def main() -> None:
    parser = argparse.ArgumentParser(description="Mnożenie macierzy metodą Binéta")
    parser.add_argument("size", type=int, help="Rozmiar kwadratowych macierzy N x N")
    parser.add_argument("seed", type=int, nargs="?", default=None, help="Seed dla generatora")
    args = parser.parse_args()

    A = random_matrix(args.size, args.seed)
    seed_B = args.seed + 1 if args.seed is not None else None
    B = random_matrix(args.size, seed_B)

    CounterFloat.reset_counter()
    start_time = time.perf_counter()
    _ = binet_multiplication(A, B)
    elapsed = time.perf_counter() - start_time

    print(f"=== Binét ({args.size}x{args.size}) ===")
    print(f"Czas: {elapsed:.6f} s")
    print(f"FLOP: {CounterFloat.operation_counts()}")


if __name__ == "__main__":
    main()
from __future__ import annotations

import argparse
import time
from typing import TypeAlias

import numpy as np

from src.matrix_algorithms.core import Matrix, random_matrix, zero_matrix
from src.matrix_algorithms.counterfloat import CounterFloat
from src.matrix_algorithms.lab1.laderman import laderman_multiply


BlockMatrix: TypeAlias = list[list[Matrix]]


def _block_add(left: BlockMatrix, right: BlockMatrix) -> BlockMatrix:
    return [
        [left[i][j] + right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def _block_sub(left: BlockMatrix, right: BlockMatrix) -> BlockMatrix:
    return [
        [left[i][j] - right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def _block_neg(matrix: BlockMatrix) -> BlockMatrix:
    return [
        [-1.0 * matrix[i][j] for j in range(len(matrix[0]))]
        for i in range(len(matrix))
    ]


def _block_matmul(left: BlockMatrix, right: BlockMatrix) -> BlockMatrix:
    rows = len(left)
    inner = len(right)
    columns = len(right[0])

    if len(left[0]) != inner:
        raise ValueError("niezgodne wymiary macierzy blokowych")

    result: BlockMatrix = []

    for i in range(rows):
        result_row: list[Matrix] = []

        for j in range(columns):
            A = left[i][0]
            B = right[0][j]

            if A.shape == (1, 1):
                value = np.array([[A[0, 0] * B[0, 0]]], dtype=object)
            else:
                value = laderman_multiply(A, B)

            for k in range(1, inner):
                value = value + laderman_multiply(left[i][k], right[k][j])

            result_row.append(value)

        result.append(result_row)

    return result


def _invert_2x2_blocks(matrix: BlockMatrix) -> BlockMatrix:
    A11, A12 = matrix[0]
    A21, A22 = matrix[1]

    A11_inv = _invert_power_of_three(A11)

    A11_inv_A12 = laderman_multiply(A11_inv, A12)
    A21_A11_inv_A12 = laderman_multiply(A21, A11_inv_A12)
    schur = A22 - A21_A11_inv_A12

    schur_inv = _invert_power_of_three(schur)

    X = A11_inv_A12
    Y = laderman_multiply(schur_inv, A21)

    top_left_correction = laderman_multiply(
        laderman_multiply(X, Y),
        A11_inv,
    )
    result_11 = A11_inv + top_left_correction

    result_12 = -1.0 * laderman_multiply(X, schur_inv)
    result_21 = -1.0 * laderman_multiply(Y, A11_inv)
    result_22 = schur_inv

    return [
        [result_11, result_12],
        [result_21, result_22],
    ]


def _invert_power_of_three(matrix: Matrix) -> Matrix:
     n = matrix.shape[0]

     if n == 1:
          value = float(matrix[0, 0])
          if value == 0.0:
               raise ValueError(
                    "zerowy pivot"
               )
          return np.array(
               [[CounterFloat(1.0) / matrix[0, 0]]],
               dtype=object,
          )

     if n % 3 != 0:
          raise ValueError("rozmiar rekurencyjnej macierzy musi byc podzielny przez 3")

     k = n // 3

     blocks = [
          [
               matrix[i * k:(i + 1) * k, j * k:(j + 1) * k]
               for j in range(3)
          ]
          for i in range(3)
     ]

     A11 = blocks[0][0]

     B: BlockMatrix = [[blocks[0][1], blocks[0][2]]]
     C: BlockMatrix = [
          [blocks[1][0]],
          [blocks[2][0]],
     ]
     D: BlockMatrix = [
          [blocks[1][1], blocks[1][2]],
          [blocks[2][1], blocks[2][2]],
          ]

     P = _invert_power_of_three(A11)

     P_B = _block_matmul([[P]], B)
     C_P = _block_matmul(C, [[P]])

     
     C_P_B = _block_matmul(C, P_B)
     schur = _block_sub(D, C_P_B)

     schur_inv = _invert_2x2_blocks(schur)

     top_right = _block_neg(_block_matmul(P_B, schur_inv))

     bottom_left = _block_neg(_block_matmul(schur_inv, C_P))

     correction = _block_matmul(P_B, _block_matmul(schur_inv, C_P), )

     top_left = P + correction[0][0]

     return np.block(
          [
               [top_left, top_right[0][0], top_right[0][1]],
               [bottom_left[0][0], schur_inv[0][0], schur_inv[0][1]],
               [bottom_left[1][0], schur_inv[1][0], schur_inv[1][1]],
          ]
     )


def laderman_inverse(matrix: Matrix) -> Matrix:
    if not isinstance(matrix, np.ndarray):
        matrix = np.asarray(matrix, dtype=object)

    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("macierz musi byc kwadratowa")

    size = matrix.shape[0]
    if size == 0:
        raise ValueError("macierz nie moze byc pusta")

    matrix_cf = np.array(
        [
            [
                value if isinstance(value, CounterFloat)
                else CounterFloat(float(value))
                for value in row
            ]
            for row in matrix
        ],
        dtype=object,
    )

    padded_size = 1
    while padded_size < size:
        padded_size *= 3

    if padded_size == size:
        padded_matrix = matrix_cf
    else:
        padded_matrix = zero_matrix(padded_size)
        padded_matrix[:size, :size] = matrix_cf

        for i in range(size, padded_size):
            padded_matrix[i, i] = CounterFloat(1.0)

    try:
        inverse_padded = _invert_power_of_three(padded_matrix)
    except ZeroDivisionError as error:
        raise ValueError(
            "dzielenie przez zero"
        ) from error

    return inverse_padded[:size, :size]


def verify_inverse(
    matrix: Matrix,
    inverse: Matrix,
    *,
    rtol: float = 1e-7,
    atol: float = 1e-7,
) -> bool:

    product = laderman_multiply(matrix, inverse)

    product_float = np.array(
        [[float(value) for value in row] for row in product],
        dtype=float,
    )

    return bool(
        np.allclose(
            product_float,
            np.eye(len(matrix)),
            rtol=rtol,
            atol=atol,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rekurencyjne odwracanie macierzy podziałem 3x3 "
                    "z mnożeniem Ladermana"
    )
    parser.add_argument("size", type=int, help="Rozmiar macierzy N x N")
    parser.add_argument("seed", type=int, nargs="?", default=None)
    args = parser.parse_args()

    if args.size <= 0:
        parser.error("Rozmiar macierzy musi być dodatni.")

    A = random_matrix(args.size, args.seed)

    for i in range(args.size):
        A[i, i] = A[i, i] + CounterFloat(float(args.size))

    CounterFloat.reset_counter()
    start_time = time.perf_counter()

    A_inv = laderman_inverse(A)

    elapsed = time.perf_counter() - start_time
    counts = CounterFloat.operation_counts()

    is_correct = verify_inverse(A, A_inv)

    print(f"=== Laderman Inverse 3x3 ({args.size}x{args.size}) ===")
    print(f"Czas: {elapsed:.6f} s")
    print(f"FLOP: {counts}")
    print(f"Poprawność (A * A_inv ≈ I): {is_correct}")


if __name__ == "__main__":
    main()

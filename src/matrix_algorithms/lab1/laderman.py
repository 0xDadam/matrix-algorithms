from __future__ import annotations

import argparse

import numpy as np

from src.matrix_algorithms.core import Matrix, random_matrix, zero_matrix
from src.matrix_algorithms.counterfloat import CounterFloat


def _split(matrix: Matrix) -> tuple[Matrix, ...]:
    size = len(matrix) // 3
    return tuple(
        matrix[
            block_row * size : (block_row + 1) * size,
            block_column * size : (block_column + 1) * size,
        ]
        for block_row in range(3)
        for block_column in range(3)
    )


def _join(blocks: tuple[Matrix, ...]) -> Matrix:
    return np.block(
        [
            [blocks[block_row * 3 + block_column] for block_column in range(3)]
            for block_row in range(3)
        ]
    )


def _laderman_power_of_three(left: Matrix, right: Matrix) -> Matrix:
    if len(left) == 1:
        return np.array([[left[0, 0] * right[0, 0]]], dtype=object)

    a11, a12, a13, a21, a22, a23, a31, a32, a33 = _split(left)
    b11, b12, b13, b21, b22, b23, b31, b32, b33 = _split(right)
    multiply = _laderman_power_of_three

    m1 = multiply(a11 + a12 + a13 - a21 - a22 - a32 - a33, b22)
    m2 = multiply(a11 - a21, b22 - b12)
    m3 = multiply(a22, b12 + b21 + b33 - b11 - b22 - b23 - b31)
    m4 = multiply(a21 - a11 + a22, b11 - b12 + b22)
    m5 = multiply(a21 + a22, b12 - b11)
    m6 = multiply(a11, b11)
    m7 = multiply(a31 - a11 + a32, b11 - b13 + b23)
    m8 = multiply(a31 - a11, b13 - b23)
    m9 = multiply(a31 + a32, b13 - b11)
    m10 = multiply(a11 + a12 + a13 - a22 - a23 - a31 - a32, b23)
    m11 = multiply(a32, b13 + b21 + b32 - b11 - b22 - b23 - b31)
    m12 = multiply(a32 - a13 + a33, b22 + b31 - b32)
    m13 = multiply(a13 - a33, b22 - b32)
    m14 = multiply(a13, b31)
    m15 = multiply(a32 + a33, b32 - b31)
    m16 = multiply(a22 - a13 + a23, b23 + b31 - b33)
    m17 = multiply(a13 - a23, b23 - b33)
    m18 = multiply(a22 + a23, b33 - b31)
    m19 = multiply(a12, b21)
    m20 = multiply(a23, b32)
    m21 = multiply(a21, b13)
    m22 = multiply(a31, b12)
    m23 = multiply(a33, b33)

    return _join(
        (
            m6 + m14 + m19,
            m1 + m4 + m5 + m6 + m12 + m14 + m15,
            m6 + m7 + m9 + m10 + m14 + m16 + m18,
            m2 + m3 + m4 + m6 + m14 + m16 + m17,
            m2 + m4 + m5 + m6 + m20,
            m14 + m16 + m17 + m18 + m21,
            m6 + m7 + m8 + m11 + m12 + m13 + m14,
            m12 + m13 + m14 + m15 + m22,
            m6 + m7 + m8 + m9 + m23,
        )
    )


def laderman_multiply(left: Matrix, right: Matrix) -> Matrix:
    """Mnoży dwie kwadratowe macierze dowolnego rozmiaru metodą Ladermana."""
    if left.ndim != 2 or right.ndim != 2 or left.shape != right.shape:
        raise ValueError("Macierze muszą być niepuste i mieć ten sam rozmiar.")
    size = len(left)
    if size == 0 or left.shape[1] != size:
        raise ValueError("Macierze muszą być kwadratowe.")

    padded_size = 1
    while padded_size < size:
        padded_size *= 3
    padded_left = zero_matrix(padded_size)
    padded_right = zero_matrix(padded_size)
    padded_left[:size, :size] = left
    padded_right[:size, :size] = right

    return _laderman_power_of_three(padded_left, padded_right)[:size, :size]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("size", type=int, help="rozmiar kwadratowych macierzy")
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    CounterFloat.reset_counter()
    laderman_multiply(
        random_matrix(args.size, args.seed),
        random_matrix(
            args.size,
            args.seed + 1 if args.seed is not None else None,
        ),
    )
    print(CounterFloat.operation_counts())

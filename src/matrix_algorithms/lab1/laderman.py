from __future__ import annotations

import random
from typing import TypeAlias

from src.matrix_algorithms.counterfloat import CounterFloat


Matrix: TypeAlias = list[list[CounterFloat]]


def _zero_matrix(size: int) -> Matrix:
    return [[CounterFloat(0.0) for _ in range(size)] for _ in range(size)]


def _add(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    return [[left[row][column] + right[row][column] for column in range(size)] for row in range(size)]


def _subtract(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    return [[left[row][column] - right[row][column] for column in range(size)] for row in range(size)]


def _split(matrix: Matrix) -> tuple[Matrix, ...]:
    size = len(matrix) // 3
    return tuple(
        [
            [
                matrix[block_row * size + row][block_column * size : (block_column + 1) * size]
                for row in range(size)
            ]
            for block_row in range(3)
            for block_column in range(3)
        ]
    )


def _join(blocks: tuple[Matrix, ...]) -> Matrix:
    size = len(blocks[0])
    return [
        [
            blocks[block_row * 3 + block_column][row][column]
            for block_column in range(3)
            for column in range(size)
        ]
        for block_row in range(3)
        for row in range(size)
    ]


def _laderman_power_of_three(left: Matrix, right: Matrix) -> Matrix:
    if len(left) == 1:
        return [[left[0][0] * right[0][0]]]

    a11, a12, a13, a21, a22, a23, a31, a32, a33 = _split(left)
    b11, b12, b13, b21, b22, b23, b31, b32, b33 = _split(right)
    multiply = _laderman_power_of_three

    m1 = multiply(_subtract(_subtract(_subtract(_subtract(_add(_add(a11, a12), a13), a21), a22), a32), a33), b22)
    m2 = multiply(_subtract(a11, a21), _subtract(b22, b12))
    m3 = multiply(a22, _subtract(_subtract(_subtract(_subtract(_add(_add(b12, b21), b33), b11), b22), b23), b31))
    m4 = multiply(_add(_subtract(a21, a11), a22), _add(_subtract(b11, b12), b22))
    m5 = multiply(_add(a21, a22), _subtract(b12, b11))
    m6 = multiply(a11, b11)
    m7 = multiply(_add(_subtract(a31, a11), a32), _add(_subtract(b11, b13), b23))
    m8 = multiply(_subtract(a31, a11), _subtract(b13, b23))
    m9 = multiply(_add(a31, a32), _subtract(b13, b11))
    m10 = multiply(_subtract(_subtract(_subtract(_subtract(_add(_add(a11, a12), a13), a22), a23), a31), a32), b23)
    m11 = multiply(a32, _subtract(_subtract(_subtract(_subtract(_add(_add(b13, b21), b32), b11), b22), b23), b31))
    m12 = multiply(_add(_subtract(a32, a13), a33), _subtract(_add(b22, b31), b32))
    m13 = multiply(_subtract(a13, a33), _subtract(b22, b32))
    m14 = multiply(a13, b31)
    m15 = multiply(_add(a32, a33), _subtract(b32, b31))
    m16 = multiply(_add(_subtract(a22, a13), a23), _subtract(_add(b23, b31), b33))
    m17 = multiply(_subtract(a13, a23), _subtract(b23, b33))
    m18 = multiply(_add(a22, a23), _subtract(b33, b31))
    m19 = multiply(a12, b21)
    m20 = multiply(a23, b32)
    m21 = multiply(a21, b13)
    m22 = multiply(a31, b12)
    m23 = multiply(a33, b33)

    return _join(
        (
            _add(_add(m6, m14), m19),
            _add(_add(_add(_add(_add(_add(m1, m4), m5), m6), m12), m14), m15),
            _add(_add(_add(_add(_add(_add(m6, m7), m9), m10), m14), m16), m18),
            _add(_add(_add(_add(_add(_add(m2, m3), m4), m6), m14), m16), m17),
            _add(_add(_add(_add(m2, m4), m5), m6), m20),
            _add(_add(_add(_add(m14, m16), m17), m18), m21),
            _add(_add(_add(_add(_add(_add(m6, m7), m8), m11), m12), m13), m14),
            _add(_add(_add(_add(m12, m13), m14), m15), m22),
            _add(_add(_add(_add(m6, m7), m8), m9), m23),
        )
    )


def laderman_multiply(left: Matrix, right: Matrix) -> Matrix:
    """Mnoży dwie kwadratowe macierze dowolnego rozmiaru metodą Ladermana."""
    if not left or len(left) != len(right):
        raise ValueError("Macierze muszą być niepuste i mieć ten sam rozmiar.")
    size = len(left)
    if any(len(row) != size for row in left + right):
        raise ValueError("Macierze muszą być kwadratowe.")

    padded_size = 1
    while padded_size < size:
        padded_size *= 3
    padded_left = _zero_matrix(padded_size)
    padded_right = _zero_matrix(padded_size)
    for row in range(size):
        padded_left[row][:size] = left[row]
        padded_right[row][:size] = right[row]

    result = _laderman_power_of_three(padded_left, padded_right)
    return [row[:size] for row in result[:size]]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("size", type=int, help="rozmiar kwadratowych macierzy")
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    CounterFloat.reset_counter()
    laderman_multiply(random_matrix(args.size, args.seed), random_matrix(args.size, args.seed + 1 if args.seed is not None else None))
    print(CounterFloat.operation_counts())
from __future__ import annotations

from collections.abc import Callable

import numpy as np

from ..core import Matrix, pad_to_power_of_two, zero_matrix
from ..counterfloat import CounterFloat


Multiply = Callable[[Matrix, Matrix], Matrix]


def _validate_matrix(matrix: Matrix) -> Matrix:
    array = np.asarray(matrix)
    if array.ndim != 2 or array.shape[0] != array.shape[1] or array.shape[0] == 0:
        raise ValueError("Macierz musi być niepusta i kwadratowa.")
    return array.astype(object, copy=False)


def _pad_with_identity(matrix: Matrix) -> Matrix:
    size = matrix.shape[0]
    padded, _, target = pad_to_power_of_two(matrix, matrix)
    if target == size:
        return matrix

    for index in range(size, target):
        padded[index, index] = CounterFloat(1.0)
    return padded


def _inverse_recursive(matrix: Matrix, multiply: Multiply) -> Matrix:
    size = matrix.shape[0]
    if size == 1:
        if float(matrix[0, 0]) == 0.0:
            raise ValueError("Macierz jest osobliwa: zerowy pivot.")
        return np.array([[CounterFloat(1.0) / matrix[0, 0]]], dtype=object)

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

    return np.block(
        [[top_left, top_right], [bottom_left, schur_inverse]]
    )


def recursive_inverse(matrix: Matrix, multiply: Multiply) -> Matrix:
    original = _validate_matrix(matrix)
    padded = _pad_with_identity(original)
    return _inverse_recursive(padded, multiply)[: original.shape[0], : original.shape[1]]


def _solve_lower(left: Matrix, right: Matrix) -> Matrix:
    size = left.shape[0]
    result = zero_matrix(size)
    for row in range(size):
        for column in range(size):
            value = right[row, column]
            for index in range(row):
                value = value - left[row, index] * result[index, column]
            result[row, column] = value / left[row, row]
    return result


def _solve_right_upper(left: Matrix, right: Matrix) -> Matrix:
    size = right.shape[0]
    result = zero_matrix(size)
    for row in range(size):
        for column in range(size):
            value = left[row, column]
            for index in range(column):
                value = value - result[row, index] * right[index, column]
            result[row, column] = value / right[column, column]
    return result


def _lu_recursive(matrix: Matrix, multiply: Multiply) -> tuple[Matrix, Matrix]:
    size = matrix.shape[0]
    if size == 1:
        if float(matrix[0, 0]) == 0.0:
            raise ValueError("Macierz jest osobliwa: zerowy pivot.")
        return (
            np.array([[CounterFloat(1.0)]], dtype=object),
            np.array([[matrix[0, 0]]], dtype=object),
        )

    half = size // 2
    a11 = matrix[:half, :half]
    a12 = matrix[:half, half:]
    a21 = matrix[half:, :half]
    a22 = matrix[half:, half:]

    l11, u11 = _lu_recursive(a11, multiply)
    u12 = _solve_lower(l11, a12)
    l21 = _solve_right_upper(a21, u11)
    schur = a22 - multiply(l21, u12)
    l22, u22 = _lu_recursive(schur, multiply)

    zero = zero_matrix(half)
    lower = np.block([[l11, zero], [l21, l22]])
    upper = np.block([[u11, u12], [zero, u22]])
    return lower, upper


def recursive_lu(matrix: Matrix, multiply: Multiply) -> tuple[Matrix, Matrix]:
    original = _validate_matrix(matrix)
    padded = _pad_with_identity(original)
    lower, upper = _lu_recursive(padded, multiply)
    size = original.shape[0]
    return lower[:size, :size], upper[:size, :size]

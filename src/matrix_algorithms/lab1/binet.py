import numpy as np
from src.matrix_algorithms.core import pad_to_power_of_two


def binet_multiplication(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, int]:
    orig_n, _ = A.shape
    _, orig_p = B.shape

    padded_A, padded_B, _ = pad_to_power_of_two(A, B)
    res_padded, ops = _binet_recursive(padded_A, padded_B)

    return res_padded[:orig_n, :orig_p], ops


def _binet_recursive(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, int]:
    n = A.shape[0]

    if n == 1:
        return A * B, 1

    mid = n // 2

    A11, A12 = A[:mid, :mid], A[:mid, mid:]
    A21, A22 = A[mid:, :mid], A[mid:, mid:]

    B11, B12 = B[:mid, :mid], B[:mid, mid:]
    B21, B22 = B[mid:, :mid], B[mid:, mid:]

    P1, ops1 = _binet_recursive(A11, B11)
    P2, ops2 = _binet_recursive(A12, B21)
    P3, ops3 = _binet_recursive(A11, B12)
    P4, ops4 = _binet_recursive(A12, B22)
    P5, ops5 = _binet_recursive(A21, B11)
    P6, ops6 = _binet_recursive(A22, B21)
    P7, ops7 = _binet_recursive(A21, B12)
    P8, ops8 = _binet_recursive(A22, B22)

    C11 = P1 + P2
    C12 = P3 + P4
    C21 = P5 + P6
    C22 = P7 + P8

    add_ops = 4 * (mid * mid)
    total_ops = ops1 + ops2 + ops3 + ops4 + ops5 + ops6 + ops7 + ops8 + add_ops

    C = np.block([[C11, C12], [C21, C22]])
    return C, total_ops
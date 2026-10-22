import numpy as np
from src.matrix_algorithms.core import pad_to_power_of_two


def strassen_multiplication(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, int]:
    orig_n, _ = A.shape
    _, orig_p = B.shape

    padded_A, padded_B, _ = pad_to_power_of_two(A, B)
    res_padded, ops = _strassen_recursive(padded_A, padded_B)

    return res_padded[:orig_n, :orig_p], ops


def _strassen_recursive(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, int]:
    n = A.shape[0]

    if n == 1:
        return A * B, 1

    mid = n // 2

    A11, A12 = A[:mid, :mid], A[:mid, mid:]
    A21, A22 = A[mid:, :mid], A[mid:, mid:]

    B11, B12 = B[:mid, :mid], B[:mid, mid:]
    B21, B22 = B[mid:, :mid], B[mid:, mid:]

    sub_size = mid * mid

    M1, ops1 = _strassen_recursive(A11 + A22, B11 + B22)
    M2, ops2 = _strassen_recursive(A21 + A22, B11)
    M3, ops3 = _strassen_recursive(A11, B12 - B22)
    M4, ops4 = _strassen_recursive(A22, B21 - B11)
    M5, ops5 = _strassen_recursive(A11 + A12, B22)
    M6, ops6 = _strassen_recursive(A21 - A11, B11 + B12)
    M7, ops7 = _strassen_recursive(A12 - A22, B21 + B22)

    prep_ops = 10 * sub_size

    C11 = M1 + M4 - M5 + M7
    C12 = M3 + M5
    C21 = M2 + M4
    C22 = M1 - M2 + M3 + M6

    combine_ops = 8 * sub_size

    total_ops = (
        ops1 + ops2 + ops3 + ops4 + ops5 + ops6 + ops7
        + prep_ops + combine_ops
    )

    C = np.block([[C11, C12], [C21, C22]])
    return C, total_ops
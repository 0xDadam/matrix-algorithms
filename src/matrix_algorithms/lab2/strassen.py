from __future__ import annotations

from ..core import Matrix
from ..lab1.strassen import strassen_multiplication
from ._recursive_factorization import recursive_inverse, recursive_lu


def strassen_inverse(matrix: Matrix) -> Matrix:
    return recursive_inverse(matrix, strassen_multiplication)


def strassen_lu(matrix: Matrix) -> tuple[Matrix, Matrix]:
    return recursive_lu(matrix, strassen_multiplication)


inverse = strassen_inverse
lu_factorization = strassen_lu

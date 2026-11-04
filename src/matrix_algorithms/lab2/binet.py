from __future__ import annotations

from ..core import Matrix
from ..lab1.binet import binet_multiplication
from ._recursive_factorization import recursive_inverse, recursive_lu


def binet_inverse(matrix: Matrix) -> Matrix:
    return recursive_inverse(matrix, binet_multiplication)


def binet_lu(matrix: Matrix) -> tuple[Matrix, Matrix]:
    return recursive_lu(matrix, binet_multiplication)


inverse = binet_inverse
lu_factorization = binet_lu

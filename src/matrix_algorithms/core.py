from typing import TypeAlias
import random

import numpy as np

from src.matrix_algorithms.counterfloat import CounterFloat


Matrix: TypeAlias = np.ndarray


def random_matrix(size: int, seed: int | None = None) -> Matrix:
    """Tworzy kwadratową macierz NumPy z elementami ``CounterFloat``."""
    if size <= 0:
        raise ValueError("Rozmiar macierzy musi być dodatni.")

    generator = random.Random(seed)
    low, high = 0.00000001, 1.0
    values = [
        [
            CounterFloat(generator.uniform(low, high))
            for _ in range(size)
        ]
        for _ in range(size)
    ]
    return np.array(values, dtype=object)


def pad_to_power_of_two(
    A: np.ndarray,
    B: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int]:
    n, m = A.shape
    m_b, p = B.shape
    max_dim = max(n, m, m_b, p)
    target_size = 1 if max_dim == 0 else 1 << (max_dim - 1).bit_length()

    padded_A = np.zeros((target_size, target_size), dtype=np.float64)
    padded_B = np.zeros((target_size, target_size), dtype=np.float64)
    padded_A[:n, :m] = A
    padded_B[:m_b, :p] = B

    return padded_A, padded_B, target_size


def zero_matrix(size: int) -> Matrix:
    """Tworzy macierz ``dtype=object`` wyzerowaną przez ``CounterFloat``."""
    result = np.empty((size, size), dtype=object)
    result.fill(CounterFloat(0.0))
    return result

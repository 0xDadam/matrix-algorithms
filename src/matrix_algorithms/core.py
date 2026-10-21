import numpy as np


def generate_random_matrix(rows: int, cols: int, seed: int | None = None) -> np.ndarray:
    if seed is not None:
        np.random.seed(seed)
    
    low = 0.00000001
    high = 1.0
    matrix = np.random.uniform(low=low, high=high, size=(rows, cols))
    
    return np.clip(matrix, low + 1e-12, high - 1e-12)


def pad_to_power_of_two(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    n, m = A.shape
    m_b, p = B.shape
    max_dim = max(n, m, m_b, p)
    
    target_size = 1 if max_dim == 0 else 1 << (max_dim - 1).bit_length()

    padded_A = np.zeros((target_size, target_size), dtype=np.float64)
    padded_B = np.zeros((target_size, target_size), dtype=np.float64)

    padded_A[:n, :m] = A
    padded_B[:m_b, :p] = B

    return padded_A, padded_B, target_size
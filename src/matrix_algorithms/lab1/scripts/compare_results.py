from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ...core import random_matrix
from ..binet import binet_multiplication
from ..laderman import laderman_multiply
from ..strassen import strassen_multiplication


def as_float(matrix: np.ndarray) -> np.ndarray:
    return np.asarray(matrix, dtype=float)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Porównuje wyniki trzech algorytmów mnożenia macierzy."
    )
    parser.add_argument("size", type=int, help="Rozmiar macierzy N x N")
    parser.add_argument("--seed", type=int, default=42, help="Seed generatora")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/matlab_comparison"),
        help="Katalog na macierze i wyniki",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Wyświetla macierze wejściowe i wyniki",
    )
    args = parser.parse_args()

    left = random_matrix(args.size, args.seed)
    right = random_matrix(args.size, args.seed + 1)
    results = {
        "Binet": binet_multiplication(left, right),
        "Strassen": strassen_multiplication(left, right),
        "Laderman": laderman_multiply(left, right),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output_dir / "A.csv", as_float(left), delimiter=",", fmt="%.17g")
    np.savetxt(args.output_dir / "B.csv", as_float(right), delimiter=",", fmt="%.17g")
    for name, result in results.items():
        np.savetxt(
            args.output_dir / f"{name.lower()}.csv",
            as_float(result),
            delimiter=",",
            fmt="%.17g",
        )

    print(f"Macierze: {args.size}x{args.size}, seed={args.seed}")
    print(f"Zapisano A, B i wyniki do: {args.output_dir}")

    if args.show:
        np.set_printoptions(precision=6, suppress=True)
        print("\nA:\n", as_float(left))
        print("\nB:\n", as_float(right))
        for name, result in results.items():
            print(f"\n{name}:\n", as_float(result))


if __name__ == "__main__":
    main()

import argparse
import csv
import os
import time

from src.matrix_algorithms.counterfloat import CounterFloat
from src.matrix_algorithms.core import random_matrix
from src.matrix_algorithms.lab2.binet import binet_inverse, binet_lu
from src.matrix_algorithms.lab2.laderman_inverse import laderman_inverse
from src.matrix_algorithms.lab2.strassen import strassen_inverse, strassen_lu


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Benchmark rekurencyjnego odwracania i faktoryzacji LU"
    )
    parser.add_argument("--max-size", type=int, default=64)
    parser.add_argument("--step", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument(
        "--output",
        type=str,
        default="results/lab2_benchmark.csv",
    )
    args = parser.parse_args()

    algorithms = {
        "Binet - odwracanie": ("odwracanie", binet_inverse),
        "Strassen - odwracanie": ("odwracanie", strassen_inverse),
        "Laderman - odwracanie": ("odwracanie", laderman_inverse),
        "Binet - LU": ("LU", binet_lu),
        "Strassen - LU": ("LU", strassen_lu),
    }

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    active = {name: True for name in algorithms}
    fieldnames = [
        "operation",
        "algorithm",
        "size",
        "time_seconds",
        "flops",
        "additions",
        "subtractions",
        "multiplications",
        "divisions",
    ]

    with open(args.output, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()

        for size in range(1, args.max_size + 1, args.step):
            if not any(active.values()):
                break

            print(f"\n--- Rozmiar N = {size} ---")
            for name, (operation, function) in algorithms.items():
                if not active[name]:
                    continue

                matrix = random_matrix(size, seed=args.seed + size)
                for index in range(size):
                    matrix[index, index] = (
                        matrix[index, index] + CounterFloat(float(size))
                    )

                CounterFloat.reset_counter()
                start = time.perf_counter()
                try:
                    function(matrix)
                    elapsed = time.perf_counter() - start
                    counts = CounterFloat.operation_counts()
                    writer.writerow({
                        "operation": operation,
                        "algorithm": name,
                        "size": size,
                        "time_seconds": elapsed,
                        "flops": counts["total"],
                        "additions": counts["additions"],
                        "subtractions": counts["subtractions"],
                        "multiplications": counts["multiplications"],
                        "divisions": counts["divisions"],
                    })
                    output.flush()
                    print(
                        f"{name:24s} czas={elapsed:.6f}s "
                        f"FLOP={counts['total']}"
                    )

                    if elapsed > args.timeout:
                        active[name] = False
                        print(f"[!] {name}: przekroczono limit czasu.")
                except Exception as error:
                    active[name] = False
                    print(f"[!] {name}: {error}")


if __name__ == "__main__":
    main()

import argparse
import csv
import os
import time

from src.matrix_algorithms.core import random_matrix
from src.matrix_algorithms.counterfloat import CounterFloat
from src.matrix_algorithms.lab1.binet import binet_multiplication
from src.matrix_algorithms.lab1.laderman import laderman_multiply
from src.matrix_algorithms.lab1.strassen import strassen_multiplication


def main() -> None:
    parser = argparse.ArgumentParser(description="Strategic benchmark on pure powers")
    parser.add_argument("--output", type=str, default="results/benchmark_powers.csv")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    powers_of_two = [2**k for k in range(1, 1)]
    powers_of_three = [3**k for k in range(5, 7)]


    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["algorithm", "size", "time_seconds", "flops"])
        writer.writeheader()

        # 1. Test Binét & Strassen
        for N in powers_of_two:
            for name, func in [("Binet", binet_multiplication), ("Strassen", strassen_multiplication)]:
                print(f"Testing {name} for N = {N}...", end="", flush=True)
                
                try:
                    A = random_matrix(N, seed=42)
                    B = random_matrix(N, seed=43)
                    
                    CounterFloat.reset_counter()
                    start = time.perf_counter()
                    _ = func(A, B)
                    elapsed = time.perf_counter() - start
                    counts = CounterFloat.operation_counts()
                    
                    writer.writerow({
                        "algorithm": name, "size": N, "time_seconds": elapsed, "flops": counts["total"]
                    })
                    f.flush()
                    print(f" [OK] Czas: {elapsed:.4f}s")
                    
                except Exception as e:
                    print(f" [FAILED]: {e}")
                    break

        # 2. Test Laderman
        for N in powers_of_three:
            print(f"Testing Laderman for N = {N}...", end="", flush=True)
            try:
                A = random_matrix(N, seed=42)
                B = random_matrix(N, seed=43)

                CounterFloat.reset_counter()
                start = time.perf_counter()
                _ = laderman_multiply(A, B)
                elapsed = time.perf_counter() - start
                counts = CounterFloat.operation_counts()
                
                writer.writerow({
                    "algorithm": "Laderman", "size": N, "time_seconds": elapsed, "flops": counts["total"]
                })
                f.flush()
                print(f" [OK] Czas: {elapsed:.4f}s")
                
            except Exception as e:
                print(f" [FAILED]: {e}")
                break

    print(f"\nZakończono. Wyniki zapisane w {args.output}")


if __name__ == "__main__":
    main()
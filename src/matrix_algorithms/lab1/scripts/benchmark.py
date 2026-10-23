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
    parser = argparse.ArgumentParser(description="Pętla pomiarowa algorytmów macierzowych")
    parser.add_argument("--max-size", type=int, default=128, help="Maksymalny rozmiar macierzy (np. 128, 500, 1000)")
    parser.add_argument("--step", type=int, default=1, help="Krok przyrostu rozmiaru macierzy")
    parser.add_argument("--seed", type=int, default=42, help="Seed generatora losowego")
    parser.add_argument("--output", type=str, default="results/benchmark_results.csv", help="Plik wyjściowy CSV")
    parser.add_argument("--timeout", type=float, default=10.0, help="Max czas w sekundach na pojedyncze mnożenie")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    algorithms = {
        "Binet": binet_multiplication,
        "Strassen": strassen_multiplication,
        "Laderman": laderman_multiply,
    }

    active_algos = {name: True for name in algorithms}
    fieldnames = ["algorithm", "size", "time_seconds", "flops", "additions", "subtractions", "multiplications", "divisions"]

    print(f"Rozpoczynanie benchmarku od 1 do {args.max_size} (krok: {args.step})...")

    with open(args.output, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for size in range(1, args.max_size + 1, args.step):
            if not any(active_algos.values()):
                print("Wszystkie algorytmy przekroczyły limit czasowy. Kończenie testu.")
                break

            print(f"\n--- Test dla rozmiaru N = {size} ---")

            for name, algo_func in algorithms.items():
                if not active_algos[name]:
                    continue

                A = random_matrix(size, seed=args.seed)
                B = random_matrix(size, seed=args.seed + 1)

                CounterFloat.reset_counter()
                start_time = time.perf_counter()

                try:
                    _ = algo_func(A, B)
                    elapsed = time.perf_counter() - start_time
                    counts = CounterFloat.operation_counts()

                    writer.writerow({
                        "algorithm": name,
                        "size": size,
                        "time_seconds": elapsed,
                        "flops": counts["total"],
                        "additions": counts["additions"],
                        "subtractions": counts["subtractions"],
                        "multiplications": counts["multiplications"],
                        "divisions": counts["divisions"],
                    })
                    f.flush()

                    print(f"  {name:10s} | Czas: {elapsed:8.4f} s | FLOP: {counts['total']:,} \
                          | ADD: {counts['additions']} | SUB: {counts['subtractions']} | MULT: {counts['multiplications']}")

                    if elapsed > args.timeout:
                        print(f"  [!] {name} przekroczył limit czasu {args.timeout}s. Zatracanie wyższych rozmiarów dla tego algorytmu.")
                        active_algos[name] = False

                except Exception as e:
                    print(f"  [!] Błąd w {name} dla rozmiaru {size}: {e}")
                    active_algos[name] = False


if __name__ == "__main__":
    main()
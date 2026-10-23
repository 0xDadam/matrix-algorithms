import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Generowanie wykresów do raportu")
    parser.add_argument("--input", type=str, default="results/benchmark_results.csv", help="Plik CSV z wynikami")
    parser.add_argument("--output-dir", type=str, default="plots", help="Katalog wyjściowy na wykresy")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Błąd: Plik '{args.input}' nie istnieje. Uruchom najpierw benchmark.")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    df = pd.read_csv(args.input)

    plt.figure(figsize=(9, 6))
    for algo in df["algorithm"].unique():
        sub_df = df[df["algorithm"] == algo]
        plt.plot(sub_df["size"], sub_df["time_seconds"], label=algo, marker="o", markersize=2)

    plt.xlabel("Rozmiar macierzy (N)")
    plt.ylabel("Czas wykonania [s]")
    plt.title("Czas działania algorytmów rekurencyjnych")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    time_plot_path = os.path.join(args.output_dir, "time_vs_size.png")
    plt.savefig(time_plot_path, dpi=300)
    plt.close()

    plt.figure(figsize=(9, 6))
    for algo in df["algorithm"].unique():
        sub_df = df[df["algorithm"] == algo]
        plt.plot(sub_df["size"], sub_df["flops"], label=algo, marker="s", markersize=2)

    plt.xlabel("Rozmiar macierzy (N)")
    plt.ylabel("Liczba operacji zmiennoprzecinkowych (FLOP)")
    plt.title("Liczba operacji zmiennoprzecinkowych")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    flops_plot_path = os.path.join(args.output_dir, "flops_vs_size.png")
    plt.savefig(flops_plot_path, dpi=300)
    plt.close()

    print(f"Pomyślnie wygenerowano wykresy w katalogu '{args.output_dir}/':")
    print(f" - {time_plot_path}")
    print(f" - {flops_plot_path}")


if __name__ == "__main__":
    main()
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd


def make_plots(input_path: str, output_dir: str) -> None:
    frame = pd.read_csv(input_path)
    os.makedirs(output_dir, exist_ok=True)
    is_powers_series = "powers" in os.path.basename(input_path)
    name_prefix = "powers_" if is_powers_series else ""
    sizes = sorted(frame["size"].unique())

    for operation, title_prefix, time_filename, flop_filename in (
        ("odwracanie", "odwracania macierzy", f"{name_prefix}inverse_time.svg", f"{name_prefix}inverse_flop.svg"),
        ("LU", "faktoryzacji LU", f"{name_prefix}lu_time.svg", f"{name_prefix}lu_flop.svg"),
    ):
        subset = frame[frame["operation"] == operation]

        figure, axis = plt.subplots(figsize=(9, 6))
        for name in subset["algorithm"].unique():
            values = subset[subset["algorithm"] == name]
            axis.plot(values["size"], values["time_seconds"], label=name)
        axis.set_xlabel("Rozmiar macierzy (N)")
        axis.set_ylabel("Czas [s]")
        axis.set_title(f"Czas {title_prefix}")
        if is_powers_series:
            axis.set_xticks(sizes)
            axis.tick_params(axis="x", labelrotation=90, labelsize=8)
        axis.grid(True, linestyle="--", alpha=0.6)
        axis.legend()
        figure.tight_layout()
        figure.savefig(os.path.join(output_dir, time_filename), format="svg")
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(9, 6))
        for name in subset["algorithm"].unique():
            values = subset[subset["algorithm"] == name]
            axis.plot(values["size"], values["flops"], label=name)
        axis.set_xlabel("Rozmiar macierzy (N)")
        axis.set_ylabel("Liczba operacji zmiennoprzecinkowych (FLOP)")
        axis.set_title(f"Liczba FLOP przy {title_prefix}")
        if is_powers_series:
            axis.set_xticks(sizes)
            axis.tick_params(axis="x", labelrotation=90, labelsize=8)
        axis.grid(True, linestyle="--", alpha=0.6)
        axis.legend()
        figure.tight_layout()
        figure.savefig(os.path.join(output_dir, flop_filename), format="svg")
        plt.close(figure)

    if is_powers_series:
        return

    figure, axis = plt.subplots(figsize=(9, 6))
    for name in frame["algorithm"].unique():
        values = frame[frame["algorithm"] == name]
        axis.plot(values["size"], values["time_seconds"], label=name)
    axis.set_xlabel("Rozmiar macierzy (N)")
    axis.set_ylabel("Czas [s]")
    axis.set_title("Czas rekurencyjnego odwracania i faktoryzacji LU")
    axis.grid(True, linestyle="--", alpha=0.6)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(os.path.join(output_dir, "czas.svg"), format="svg")
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(9, 6))
    for name in frame["algorithm"].unique():
        values = frame[frame["algorithm"] == name]
        axis.plot(values["size"], values["flops"], label=name)
    axis.set_xlabel("Rozmiar macierzy (N)")
    axis.set_ylabel("Liczba operacji zmiennoprzecinkowych (FLOP)")
    axis.set_title("FLOP rekurencyjnego odwracania i faktoryzacji LU")
    axis.grid(True, linestyle="--", alpha=0.6)
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(os.path.join(output_dir, "flop.svg"), format="svg")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generowanie wykresów labu 2")
    parser.add_argument("--input", default="results/lab2_benchmark.csv")
    parser.add_argument("--output-dir", default="src/matrix_algorithms/lab2")
    args = parser.parse_args()
    make_plots(args.input, args.output_dir)


if __name__ == "__main__":
    main()

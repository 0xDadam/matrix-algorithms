import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()

    matlab = np.atleast_2d(np.loadtxt(args.directory / "matlab.csv", delimiter=","))
    for name in ("binet", "strassen", "laderman"):
        result = np.loadtxt(args.directory / f"{name}.csv", delimiter=",")
        error = np.max(np.abs(result - matlab))
        print(f"{name:8s}: {'OK' if np.allclose(result, matlab) else 'BŁĄD'} ({error:.3e})")


if __name__ == "__main__":
    main()

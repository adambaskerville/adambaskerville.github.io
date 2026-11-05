from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

cache_dir = Path(tempfile.gettempdir()) / "matplotlib"
os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HISTORICAL_MONTHS = 2.0
SECONDS_PER_MONTH = 30.0 * 24.0 * 3600.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot convergence and timing data from original_pekeris.py benchmark output."
    )
    parser.add_argument(
        "json_file",
        type=Path,
        help="JSON file created by original_pekeris.py --json",
    )
    parser.add_argument(
        "--energy-plot",
        type=Path,
        default=Path("pekeris_energy_convergence.png"),
        help="Output path for the energy convergence plot.",
    )
    parser.add_argument(
        "--time-plot",
        type=Path,
        default=Path("pekeris_timing.png"),
        help="Output path for the timing plot.",
    )
    parser.add_argument(
        "--historical-months",
        type=float,
        default=HISTORICAL_MONTHS,
        help="Historical WEIZAC runtime in months used for the headline speedup estimate.",
    )
    return parser.parse_args()


def load_results(path: Path) -> list[dict]:
    data = json.loads(path.read_text())
    if not isinstance(data, list) or not data:
        raise ValueError("Benchmark JSON must contain a non-empty list of result dictionaries.")
    return data


def save_energy_plot(results: list[dict], path: Path) -> None:
    omegas = [item["omega"] for item in results]
    energies = [item["energy_hartree"] for item in results]

    plt.figure(figsize=(7, 4.5))
    plt.plot(omegas, energies, marker="o", linewidth=1.8)
    plt.xlabel(r"Matrix order $\omega$")
    plt.ylabel("Energy / Hartree")
    plt.title("Pekeris Helium Convergence")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def save_timing_plot(results: list[dict], path: Path) -> None:
    sizes = [item["basis_size"] for item in results]
    total_seconds = [item["total_seconds"] for item in results]

    plt.figure(figsize=(7, 4.5))
    plt.plot(sizes, total_seconds, marker="o", linewidth=1.8)
    plt.xlabel("Basis size")
    plt.ylabel("Runtime / s")
    plt.title("Pekeris Helium Runtime")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def print_summary(results: list[dict], historical_months: float) -> None:
    final = results[-1]
    historical_seconds = historical_months * SECONDS_PER_MONTH
    speedup = historical_seconds / final["total_seconds"]

    print("Final benchmark point")
    print(
        f"omega={final['omega']}, size={final['basis_size']}, "
        f"epsilon={final['epsilon']:.9f}, energy={final['energy_hartree']:.12f} Ha"
    )
    print(f"runtime={final['total_seconds']:.4f} s")
    print(
        f"if WEIZAC took about {historical_months:.1f} months, "
        f"the speedup is roughly {speedup:,.0f}x"
    )


def main() -> None:
    args = parse_args()
    results = load_results(args.json_file)
    save_energy_plot(results, args.energy_plot)
    save_timing_plot(results, args.time_plot)
    print_summary(results, args.historical_months)
    print(f"Saved {args.energy_plot}")
    print(f"Saved {args.time_plot}")


if __name__ == "__main__":
    main()

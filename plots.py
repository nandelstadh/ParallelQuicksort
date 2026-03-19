#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt

DISTRIBUTION_LABELS = {
    "u": "Uniform",
    "n": "Normal",
    "e": "Exponential",
}


def load_results(csv_path: Path):
    data = {}
    with csv_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            n = int(row["n"].strip())
            distribution = row["d"].strip()
            threads = int(row["t"].strip())
            runtime = float(row["result"].strip())
            data.setdefault((n, distribution), {})[threads] = runtime
    return data


def compute_speedups(runtimes_by_threads):
    if 1 not in runtimes_by_threads:
        raise ValueError("Missing single-thread baseline (t=1).")
    baseline = runtimes_by_threads[1]
    return {
        threads: baseline / runtime
        for threads, runtime in sorted(runtimes_by_threads.items())
    }


def plot_speedup_for_n(n, all_results, output_dir: Path, show: bool):
    fig, ax = plt.subplots(figsize=(8, 5))

    plotted_lines = 0
    for distribution in ("u", "n", "e"):
        runtimes = all_results.get((n, distribution))
        if not runtimes:
            continue
        speedups = compute_speedups(runtimes)
        threads = list(speedups.keys())
        values = list(speedups.values())
        ax.plot(
            threads,
            values,
            marker="o",
            linewidth=2,
            label=DISTRIBUTION_LABELS[distribution],
        )
        plotted_lines += 1

    if plotted_lines == 0:
        plt.close(fig)
        raise ValueError(f"No data found for n={n}.")

    ax.set_title(f"Parallel quicksort speedup for N={n}")
    ax.set_xlabel("Number of threads")
    ax.set_ylabel("Speedup (times)")
    ax.set_xticks(sorted({1, 2, 4, 8}))
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend()
    fig.tight_layout()

    output_path = output_dir / f"speedup_n_{n}.png"
    fig.savefig(output_path, dpi=150)

    if show:
        plt.show()
    plt.close(fig)

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot speedup vs thread count from results.csv. "
            "Each plot overlays uniform, normal, and exponential distributions."
        )
    )
    parser.add_argument(
        "--csv", default="results.csv", help="Path to results CSV file."
    )
    parser.add_argument(
        "--n", type=int, help="Fixed n to plot. If omitted, generates one plot per n."
    )
    parser.add_argument(
        "--output-dir", default=".", help="Directory where plot image(s) are written."
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plot window(s) in addition to saving.",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = load_results(csv_path)
    available_ns = sorted({n for (n, _) in all_results})
    if not available_ns:
        raise ValueError("No rows were parsed from the CSV file.")

    ns_to_plot = [args.n] if args.n is not None else available_ns
    for n in ns_to_plot:
        output_path = plot_speedup_for_n(n, all_results, output_dir, args.show)
        print(f"Saved {output_path}")


if __name__ == "__main__":
    main()

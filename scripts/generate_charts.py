#!/usr/bin/env python3
"""Generate SVG benchmark charts from per-compiler JSON result files.

Usage: generate_charts.py <results-dir> <output-dir>

<results-dir> layout (from actions/download-artifact with merge-multiple=false):
  results-dir/
    bench-results-gcc-14/
      bench_horner_cxx17.json
      bench_horner_cxx20.json
      bench1D_cxx17.json
      ...
    bench-results-llvm-21/
      ...

Generates:
  horner_performance.svg
  fitting_performance.svg
  cross_compiler_overview.svg
  average_improvement.svg
"""
import json
import math
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("matplotlib/numpy not installed — skipping chart generation", file=sys.stderr)
    sys.exit(0)

COMPILER_COLORS = {
    "gcc-14": "#4CAF50",
    "gcc-15": "#2E7D32",
    "llvm-20": "#2196F3",
    "llvm-21": "#0D47A1",
}
DEFAULT_COLORS = ["#4CAF50", "#2E7D32", "#2196F3", "#0D47A1", "#FF9800", "#9C27B0"]


def load_results(results_dir: Path) -> dict[str, dict[str, list[dict]]]:
    """Return {compiler: {bench_name: rows}}."""
    data: dict[str, dict[str, list[dict]]] = {}
    for compiler_dir in sorted(results_dir.iterdir()):
        if not compiler_dir.is_dir():
            continue
        # Strip "bench-results-" prefix if present
        compiler = compiler_dir.name.removeprefix("bench-results-")
        data[compiler] = {}
        for json_file in sorted(compiler_dir.glob("*.json")):
            payload = json.loads(json_file.read_text())
            data[compiler][payload["title"]] = payload["rows"]
    return data


def _extract_numeric_label(name: str) -> float:
    """Try to extract a leading number from a benchmark name (e.g. degree)."""
    import re
    m = re.search(r"\d+", name)
    return float(m.group()) if m else float("nan")


def plot_horner_performance(data: dict, out_dir: Path) -> None:
    """ns/op by benchmark variant for horner benches across cxx17/cxx20."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
    fig.suptitle("Horner Evaluation Performance", fontsize=14, fontweight="bold")

    for ax, std in zip(axes, ["cxx17", "cxx20"]):
        bench_key = f"bench_horner_{std}"
        ax.set_title(f"C++{std[3:]}")
        ax.set_xlabel("Benchmark")
        ax.set_ylabel("ns / op")

        compilers = [c for c in data if bench_key in data[c]]
        if not compilers:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue

        # Gather all benchmark names from first compiler
        first = compilers[0]
        names = [r["name"] for r in data[first][bench_key]]
        x = np.arange(len(names))
        width = 0.8 / max(len(compilers), 1)

        for i, compiler in enumerate(compilers):
            rows = {r["name"]: r for r in data[compiler].get(bench_key, [])}
            vals = [rows.get(n, {}).get("ns_per_op", float("nan")) for n in names]
            color = COMPILER_COLORS.get(compiler, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
            ax.bar(x + i * width, vals, width, label=compiler, color=color, alpha=0.85)

        ax.set_xticks(x + width * (len(compilers) - 1) / 2)
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out_path = out_dir / "horner_performance.svg"
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_fitting_performance(data: dict, out_dir: Path) -> None:
    """Fitting time across compilers for 1D benches."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
    fig.suptitle("1D Fitting Performance", fontsize=14, fontweight="bold")

    for ax, std in zip(axes, ["cxx17", "cxx20"]):
        bench_key = f"bench1D_{std}"
        ax.set_title(f"C++{std[3:]}")
        ax.set_xlabel("Benchmark")
        ax.set_ylabel("ns / op")

        compilers = [c for c in data if bench_key in data[c]]
        if not compilers:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue

        first = compilers[0]
        names = [r["name"] for r in data[first][bench_key]]
        x = np.arange(len(names))
        width = 0.8 / max(len(compilers), 1)

        for i, compiler in enumerate(compilers):
            rows = {r["name"]: r for r in data[compiler].get(bench_key, [])}
            vals = [rows.get(n, {}).get("ns_per_op", float("nan")) for n in names]
            color = COMPILER_COLORS.get(compiler, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
            ax.bar(x + i * width, vals, width, label=compiler, color=color, alpha=0.85)

        ax.set_xticks(x + width * (len(compilers) - 1) / 2)
        ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out_path = out_dir / "fitting_performance.svg"
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_cross_compiler_overview(data: dict, out_dir: Path) -> None:
    """Relative performance vs GCC baseline across all bench suites."""
    # Pick a baseline compiler (prefer gcc-14, else first available)
    baseline = next((c for c in ["gcc-14", "gcc-15"] if c in data), None)
    if baseline is None or len(data) < 2:
        print("Skipping cross-compiler overview: no baseline or < 2 compilers")
        return

    bench_suites = [
        "bench_horner_cxx17", "bench_horner_cxx20",
        "bench1D_cxx17", "bench1D_cxx20",
        "benchMany_cxx17", "benchMany_cxx20",
        "bench_ND_cxx17", "bench_ND_cxx20",
    ]

    compilers = [c for c in data if c != baseline]
    fig, ax = plt.subplots(figsize=(12, 6))
    fig.suptitle(f"Relative Performance vs {baseline} (lower = faster)", fontsize=13, fontweight="bold")

    suite_labels = []
    compiler_ratios: dict[str, list[float]] = {c: [] for c in compilers}

    for suite in bench_suites:
        base_rows = data[baseline].get(suite)
        if not base_rows:
            continue
        base_map = {r["name"]: r["ns_per_op"] for r in base_rows}
        suite_labels.append(suite.replace("bench_", "").replace("_cxx", " C++"))

        for compiler in compilers:
            comp_rows = data[compiler].get(suite, [])
            comp_map = {r["name"]: r["ns_per_op"] for r in comp_rows}
            # Geometric mean of ratios over shared benchmark names
            ratios = []
            for name, base_val in base_map.items():
                comp_val = comp_map.get(name)
                if comp_val and base_val and base_val > 0:
                    ratios.append(comp_val / base_val)
            geomean = math.exp(sum(math.log(r) for r in ratios) / len(ratios)) if ratios else float("nan")
            compiler_ratios[compiler].append(geomean)

    x = np.arange(len(suite_labels))
    width = 0.8 / max(len(compilers), 1)

    for i, compiler in enumerate(compilers):
        color = COMPILER_COLORS.get(compiler, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        ax.bar(x + i * width, compiler_ratios[compiler], width, label=compiler, color=color, alpha=0.85)

    ax.axhline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.5, label=f"{baseline} (baseline)")
    ax.set_xticks(x + width * (len(compilers) - 1) / 2)
    ax.set_xticklabels(suite_labels, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Relative ns/op (ratio to baseline)")
    ax.legend(fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out_path = out_dir / "cross_compiler_overview.svg"
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_average_improvement(data: dict, out_dir: Path) -> None:
    """Geometric mean speedup bar chart per compiler vs gcc-14 baseline."""
    baseline = next((c for c in ["gcc-14", "gcc-15"] if c in data), None)
    if baseline is None or len(data) < 2:
        print("Skipping average improvement chart: no baseline or < 2 compilers")
        return

    compilers = [c for c in data if c != baseline]
    geomeans: list[float] = []

    for compiler in compilers:
        ratios: list[float] = []
        for suite, rows in data[baseline].items():
            base_map = {r["name"]: r["ns_per_op"] for r in rows}
            comp_rows = data[compiler].get(suite, [])
            comp_map = {r["name"]: r["ns_per_op"] for r in comp_rows}
            for name, base_val in base_map.items():
                comp_val = comp_map.get(name)
                if comp_val and base_val and base_val > 0:
                    # speedup = base/comp (>1 means comp is faster)
                    ratios.append(base_val / comp_val)
        geomean = math.exp(sum(math.log(r) for r in ratios) / len(ratios)) if ratios else float("nan")
        geomeans.append(geomean)

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle(f"Average Speedup vs {baseline}", fontsize=13, fontweight="bold")

    colors = [COMPILER_COLORS.get(c, DEFAULT_COLORS[i % len(DEFAULT_COLORS)]) for i, c in enumerate(compilers)]
    bars = ax.bar(compilers, geomeans, color=colors, alpha=0.85)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1, alpha=0.5)

    for bar, val in zip(bars, geomeans):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{val:.3f}x", ha="center", va="bottom", fontsize=10)

    ax.set_ylabel("Geometric mean speedup (>1 = faster than baseline)")
    ax.set_xlabel("Compiler")
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    out_path = out_dir / "average_improvement.svg"
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def main() -> None:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <results-dir> <output-dir>", file=sys.stderr)
        sys.exit(1)

    results_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    data = load_results(results_dir)
    if not data:
        print("No benchmark data found.", file=sys.stderr)
        sys.exit(0)

    print(f"Loaded data for compilers: {list(data.keys())}")

    plot_horner_performance(data, out_dir)
    plot_fitting_performance(data, out_dir)
    plot_cross_compiler_overview(data, out_dir)
    plot_average_improvement(data, out_dir)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Cross-compiler benchmark analysis.

Usage: analyze_bench.py <results-dir> <output-dir>

Generates:
  summary.md  — Markdown table: bench name × compiler → ns/op
  summary.csv — CSV equivalent
"""
import csv
import json
import sys
from pathlib import Path


def load_results(results_dir: Path) -> dict[str, dict[str, list[dict]]]:
    """Return {compiler: {bench_name: rows}}."""
    data: dict[str, dict[str, list[dict]]] = {}
    for compiler_dir in sorted(results_dir.iterdir()):
        if not compiler_dir.is_dir():
            continue
        compiler = compiler_dir.name.removeprefix("bench-results-")
        data[compiler] = {}
        for json_file in sorted(compiler_dir.glob("*.json")):
            payload = json.loads(json_file.read_text())
            data[compiler][payload["title"]] = payload["rows"]
    return data


def build_table(data: dict) -> tuple[list[str], list[str], dict[tuple[str, str], float]]:
    """Return (compilers, row_keys, {(bench_suite::row_name, compiler): ns_per_op})."""
    compilers = sorted(data.keys())

    # Collect all (suite, row_name) pairs in stable order
    row_keys: list[str] = []
    seen: set[str] = set()
    for compiler in compilers:
        for suite, rows in sorted(data[compiler].items()):
            for row in rows:
                key = f"{suite}::{row['name']}"
                if key not in seen:
                    seen.add(key)
                    row_keys.append(key)

    values: dict[tuple[str, str], float] = {}
    for compiler in compilers:
        for suite, rows in data[compiler].items():
            for row in rows:
                key = f"{suite}::{row['name']}"
                values[(key, compiler)] = row["ns_per_op"]

    return compilers, row_keys, values


def write_markdown(
    compilers: list[str],
    row_keys: list[str],
    values: dict[tuple[str, str], float],
    out_path: Path,
) -> None:
    lines: list[str] = [
        "# Benchmark Summary\n",
        "All values are **ns/op** (lower is better).\n",
        "",
    ]

    header = "| Benchmark | " + " | ".join(compilers) + " |"
    sep = "| --- |" + " ---: |" * len(compilers)
    lines += [header, sep]

    for key in row_keys:
        parts = [f"`{key}`"]
        for compiler in compilers:
            val = values.get((key, compiler))
            parts.append(f"{val:.2f}" if val is not None else "—")
        lines.append("| " + " | ".join(parts) + " |")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Saved {out_path}")


def write_csv(
    compilers: list[str],
    row_keys: list[str],
    values: dict[tuple[str, str], float],
    out_path: Path,
) -> None:
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["benchmark"] + compilers)
        for key in row_keys:
            row = [key]
            for compiler in compilers:
                val = values.get((key, compiler))
                row.append(f"{val:.4f}" if val is not None else "")
            writer.writerow(row)
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

    compilers, row_keys, values = build_table(data)
    write_markdown(compilers, row_keys, values, out_dir / "summary.md")
    write_csv(compilers, row_keys, values, out_dir / "summary.csv")


if __name__ == "__main__":
    main()

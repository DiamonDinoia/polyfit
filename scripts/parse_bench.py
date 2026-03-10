#!/usr/bin/env python3
"""Parse nanobench stdout into JSON.

Usage: parse_bench.py <input.txt> <output.json>

Nanobench output format (text table):
|               ns/op |                op/s |    err% |          ins/op |          bra/op |   miss% |     total | benchmark
|--------------------:|--------------------:|--------:|----------------:|----------------:|--------:|----------:|:----------
|               12.34 |       81,037,277.15 |    0.1% |           42.00 |            8.00 |    0.0% |      0.00 | `my bench name`
"""
import json
import re
import sys
from pathlib import Path


def parse_nanobench_text(text: str) -> list[dict]:
    """Parse nanobench pipe-delimited table output into a list of row dicts."""
    results = []
    header_seen = False

    for line in text.splitlines():
        line = line.strip()
        # Skip separator lines (all dashes and pipes)
        if re.match(r"^\|[-: |]+\|$", line):
            header_seen = True
            continue
        if not line.startswith("|"):
            continue
        # Split on pipe
        parts = [p.strip() for p in line.split("|")]
        # parts[0] is empty (before first |), parts[-1] is empty (after last |)
        parts = [p for p in parts if p != ""]

        # Header row: contains "ns/op"
        if "ns/op" in parts[0]:
            header_seen = True
            continue
        if not header_seen:
            continue

        # Data row: expect at least 9 fields
        # [ns/op, op/s, err%, ins/op, bra/op, miss%, total, benchmark]
        # (some builds omit bra/op and miss% — handle both 7-col and 9-col)
        if len(parts) < 4:
            continue

        try:
            ns_per_op_str = parts[0].replace(",", "")
            ops_per_sec_str = parts[1].replace(",", "")
            err_pct_str = parts[2].rstrip("%")
            # benchmark name is last column, strip backticks
            name = parts[-1].strip("`").strip()

            row = {
                "name": name,
                "ns_per_op": float(ns_per_op_str),
                "ops_per_sec": float(ops_per_sec_str),
                "err_pct": float(err_pct_str),
            }
            # Optionally parse ins/op if present
            if len(parts) >= 5:
                try:
                    row["ins_per_op"] = float(parts[3].replace(",", ""))
                except ValueError:
                    pass
            results.append(row)
        except (ValueError, IndexError):
            continue

    return results


def main() -> None:
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.txt> <output.json>", file=sys.stderr)
        sys.exit(1)

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    text = input_path.read_text(encoding="utf-8", errors="replace")
    rows = parse_nanobench_text(text)

    # Use the input filename (minus extension) as the title
    title = input_path.stem

    result = {"title": title, "rows": rows}
    output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"Parsed {len(rows)} rows from '{input_path}' -> '{output_path}'")


if __name__ == "__main__":
    main()

"""
Normalize / deduplicate mammalian ortholog targets per (Query ID, gene prefix) in SEA results.

- If *_HUMAN exists for a prefix, drop *_RAT, *_MOUSE, *_BOVIN, *_GORGO for that prefix.
- If no *_HUMAN but mammal ortholog rows exist, rewrite Target ID to {prefix}_HUMAN and keep a
  single best row (lowest P-Value, then highest Z-Score).
- Other species (e.g. ARATH, YEAST) on the same prefix are kept as-is.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

MAMMAL_ORTHOLOG_SUFFIXES = frozenset({"RAT", "MOUSE", "BOVIN", "GORGO"})


def split_target_id(target_id: str) -> tuple[str, str | None]:
    if "_" not in target_id:
        return target_id, None
    prefix, species = target_id.rsplit("_", 1)
    return prefix, species


def _p_value_float(v: str) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("inf")


def _z_score_float(v: str) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("-inf")


def pick_representative_mammal_row(rows: list[dict], p_col: str, z_col: str) -> dict:
    best = min(rows, key=lambda r: (_p_value_float(r[p_col]), -_z_score_float(r[z_col])))
    return dict(best)


def process_prefix_block(rows: list[dict], target_col: str, p_col: str, z_col: str) -> list[dict]:
    human_rows: list[dict] = []
    mammal_rows: list[dict] = []
    other_rows: list[dict] = []

    for row in rows:
        _, sp = split_target_id(row[target_col])
        if sp == "HUMAN":
            human_rows.append(row)
        elif sp in MAMMAL_ORTHOLOG_SUFFIXES:
            mammal_rows.append(row)
        else:
            other_rows.append(row)

    if human_rows:
        return human_rows + other_rows

    if mammal_rows:
        rep = pick_representative_mammal_row(mammal_rows, p_col, z_col)
        prefix, _ = split_target_id(rep[target_col])
        rep[target_col] = f"{prefix}_HUMAN"
        return [rep] + other_rows

    return list(rows)


def filter_sea_results(
    rows: list[dict],
    fieldnames: list[str],
    query_col: str = "Query ID",
    target_col: str = "Target ID",
    p_col: str = "P-Value",
    z_col: str = "Z-Score",
) -> list[dict]:
    # Group by (Query ID, gene prefix)
    groups: dict[tuple[str, str], list[dict]] = {}
    for row in rows:
        prefix, _ = split_target_id(row[target_col])
        key = (row[query_col], prefix)
        groups.setdefault(key, []).append(row)

    out: list[dict] = []
    # Preserve global order: iterate rows in original order, emit each group once when first seen
    seen: set[tuple[str, str]] = set()
    for row in rows:
        prefix, _ = split_target_id(row[target_col])
        key = (row[query_col], prefix)
        if key in seen:
            continue
        seen.add(key)
        out.extend(process_prefix_block(groups[key], target_col, p_col, z_col))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter mammal ortholog duplicates in SEA CSV.")
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent / "sea_results.csv",
        help="Input CSV path",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "sea_results_filtered.csv",
        help="Output CSV path",
    )
    args = parser.parse_args()

    with open(args.input, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        if fieldnames is None:
            raise SystemExit("CSV has no header row")
        rows = list(reader)

    filtered = filter_sea_results(rows, fieldnames)
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(filtered)

    print(f"Read {len(rows)} rows from {args.input}")
    print(f"Wrote {len(filtered)} rows to {args.output} ({len(rows) - len(filtered)} rows removed)")


if __name__ == "__main__":
    main()

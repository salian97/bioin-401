#!/usr/bin/env python3
"""
Build SEA enrichment tables: attach DrugBank IDs from structure_links.csv and fetch
target protein sequences + UniProt accessions via the UniProt REST API (same approach
as get_protein_sequences.fetch_uniprot_data).

- Default: process all of sea_results_filtered.csv -> sea_data.csv (sorted by DrugBank ID).
- ``--sample``: seven-drug subset -> sea_sample.csv (legacy).
"""
from __future__ import annotations

import argparse
import csv
import io
import re
import time
from pathlib import Path

import requests

SCRIPT_DIR = Path(__file__).resolve().parent
SEA_PATH = SCRIPT_DIR / "sea_results_filtered.csv"
STRUCTURE_LINKS = (
    SCRIPT_DIR.parent / "yaml-generation" / "structure_links.csv"
)
OUT_SAMPLE = SCRIPT_DIR / "sea_sample.csv"
OUT_FULL = SCRIPT_DIR / "sea_data.csv"

DRUG_NAMES = [
    "Goserelin",
    "Desmopressin",
    "Pravastatin",
    "Ramipril",
    "Baclofen",
    "Pentagastrin",
    "Nicotine",
]

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
CHUNK_SIZE = 50
SLEEP_S = 0.2


def drugbank_id_sort_key(drugbank_id: str) -> tuple[int, str]:
    drugbank_id = (drugbank_id or "").strip()
    m = re.match(r"DB(\d+)", drugbank_id, re.I)
    if m:
        return (int(m.group(1)), drugbank_id.upper())
    return (10**9, drugbank_id)


def load_drugbank_by_name(structure_links_path: Path) -> dict[str, str]:
    """Map exact DrugBank `Name` -> `DrugBank ID`."""
    out: dict[str, str] = {}
    with structure_links_path.open(newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = (row.get("Name") or "").strip()
            dbid = (row.get("DrugBank ID") or "").strip()
            if name and dbid:
                out[name] = dbid
    return out


def fetch_uniprot_batch(entry_names: list[str]) -> dict[str, tuple[str, str]]:
    """
    For each UniProt entry name (e.g. GABR1_HUMAN), return (accession, sequence).
    Some names map to multiple accessions (e.g. reviewed + unreviewed); prefer a
    row with a non-empty sequence, then the longest sequence.
    """
    result: dict[str, tuple[str, str]] = {}
    if not entry_names:
        return result
    query = " OR ".join(f"id:{name}" for name in entry_names)
    # Request extra rows so duplicates (same Entry Name) are not truncated.
    size = min(500, max(len(entry_names) * 3, len(entry_names)))
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,id,sequence",
        "size": size,
    }
    r = requests.get(UNIPROT_SEARCH, params=params, timeout=60)
    r.raise_for_status()
    reader = csv.DictReader(io.StringIO(r.text), delimiter="\t")
    for row in reader:
        ename = row.get("Entry Name", "").strip()
        acc = row.get("Entry", "").strip()
        seq = row.get("Sequence", "").strip()
        if not ename:
            continue
        prev = result.get(ename)
        if prev is None:
            result[ename] = (acc, seq)
            continue
        _, prev_seq = prev
        if len(seq) > len(prev_seq):
            result[ename] = (acc, seq)
    return result


OUT_FIELDS_DB_FIRST = [
    "DrugBank ID",
    "Query ID",
    "Target ID",
    "UniProt ID",
    "Query Smiles",
    "Target Sequence",
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Build SEA table with DrugBank IDs and UniProt sequences.")
    parser.add_argument(
        "--sample",
        action="store_true",
        help=f"Only {len(DRUG_NAMES)} named drugs; write {OUT_SAMPLE.name} (legacy).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output CSV path (defaults: sea_data.csv or sea_sample.csv with --sample).",
    )
    args = parser.parse_args()

    sample_mode = args.sample
    out_path = args.output
    if out_path is None:
        out_path = OUT_SAMPLE if sample_mode else OUT_FULL

    drug_set = set(DRUG_NAMES) if sample_mode else None
    db_by_name = load_drugbank_by_name(STRUCTURE_LINKS)
    if sample_mode:
        missing_db = drug_set - set(db_by_name.keys())
        if missing_db:
            raise SystemExit(f"DrugBank ID not found in structure_links for: {sorted(missing_db)}")

    rows_out: list[dict[str, str]] = []
    with SEA_PATH.open(newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise SystemExit("SEA file has no header")
        for row in reader:
            qid = (row.get("Query ID") or "").strip()
            if sample_mode and qid not in drug_set:
                continue
            rows_out.append(row)

    if not rows_out:
        raise SystemExit("No matching rows in sea_results_filtered.csv")

    # Primary organization: DrugBank ID, then drug name, then target (tie-break only).
    rows_out.sort(
        key=lambda r: (
            drugbank_id_sort_key(db_by_name.get((r.get("Query ID") or "").strip(), "")),
            (r.get("Query ID") or "").strip().lower(),
            (r.get("Target ID") or "").strip(),
        )
    )

    unique_targets: list[str] = []
    seen: set[str] = set()
    for r in rows_out:
        tid = (r.get("Target ID") or "").strip()
        if tid and tid not in seen:
            seen.add(tid)
            unique_targets.append(tid)

    seq_map: dict[str, tuple[str, str]] = {}
    for i in range(0, len(unique_targets), CHUNK_SIZE):
        chunk = unique_targets[i : i + CHUNK_SIZE]
        try:
            batch = fetch_uniprot_batch(chunk)
            seq_map.update(batch)
        except Exception as e:
            print(f"Warning: UniProt chunk {i}-{i+len(chunk)} failed: {e}")
        time.sleep(SLEEP_S)

    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=OUT_FIELDS_DB_FIRST)
        w.writeheader()
        for r in rows_out:
            qid = (r.get("Query ID") or "").strip()
            tid = (r.get("Target ID") or "").strip()
            smiles = r.get("Query Smiles") or ""
            acc, seq = seq_map.get(tid, ("", ""))
            w.writerow(
                {
                    "DrugBank ID": db_by_name.get(qid, ""),
                    "Query ID": qid,
                    "Target ID": tid,
                    "UniProt ID": acc,
                    "Query Smiles": smiles,
                    "Target Sequence": seq,
                }
            )

    n_seq = sum(
        1
        for r in rows_out
        if seq_map.get((r.get("Target ID") or "").strip(), ("", ""))[1]
    )
    print(f"Wrote {len(rows_out)} rows to {out_path}")
    print(f"Unique targets: {len(unique_targets)}; rows with non-empty sequence: {n_seq}")


if __name__ == "__main__":
    main()

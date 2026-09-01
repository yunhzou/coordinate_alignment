#!/usr/bin/env python3
"""Recover flushed shard records and mark unfinished catalog rows."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import zlib
from collections import Counter
from pathlib import Path

from rdkit import Chem

from rxn_core.fragment_matching.serialization import FRAGMENT_DETECTION_SCHEMA


def _recovered_records(path):
    decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
    text = decoder.decompress(path.read_bytes()).decode("utf-8")
    for line in text.splitlines():
        if line.strip():
            yield json.loads(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parts", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--id-column", default="Inventory ID")
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--branch-limit", type=int, required=True)
    parser.add_argument("--seed-limit", type=int, required=True)
    parser.add_argument("--elapsed-seconds", type=float, required=True)
    args = parser.parse_args()

    recovered = {}
    for path in sorted(Path(args.parts).glob("part_*.jsonl.gz")):
        for record in _recovered_records(path):
            recovered[record["source_id"]] = record

    with gzip.open(args.catalog, "rt", encoding="utf-8") as stream:
        catalog_rows = list(csv.DictReader(stream))
    records = []
    unresolved = []
    for row_index, row in enumerate(catalog_rows):
        source_id = row[args.id_column]
        record = recovered.get(source_id)
        if record is None:
            unresolved.append(source_id)
            record = {
                "schema": FRAGMENT_DETECTION_SCHEMA,
                "row_index": row_index,
                "source_id": source_id,
                "representation": row["SMILES"],
                "status": "timeout",
                "complete": False,
                "branch_limit": args.branch_limit,
                "maximum_branch_count": 0,
                "capped_seed_count": 0,
                "best_fragment_size": 0,
                "candidates": [],
            }
        elif (Chem.AddHs(Chem.MolFromSmiles(row["SMILES"])).GetNumAtoms()
              > args.seed_limit and record["status"] != "capped"):
            record["status"] = "seed_limited"
            record["complete"] = False
        records.append(record)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(record, separators=(",", ":")) + "\n")

    counts = Counter()
    counts["rows"] = len(records)
    counts["searched"] = len(recovered)
    counts["unresolved"] = len(unresolved)
    counts["matched_precursors"] = sum(
        bool(record["candidates"]) for record in records)
    counts["fragment_candidates"] = sum(
        len(record["candidates"]) for record in records)
    counts["capped"] = sum(
        record["status"] == "capped" for record in records)
    counts["seed_limited"] = sum(
        record["status"] == "seed_limited" for record in records)
    summary = {
        "schema": "rxn_core.retro_catalog_summary/v3",
        "target_smiles": args.target_smiles,
        "catalog": args.catalog,
        "output": str(output),
        "shard_index": 0,
        "shard_count": 1,
        "workers": 96,
        "config": {
            "minimum_fragment_size": 1,
            "iso_tolerance": 0.5,
            "branch_limit": args.branch_limit,
            "candidate_limit": 100,
            "seed_limit": args.seed_limit,
            "maximum_boundary_bonds": None,
            "maximum_leftover_fragments": None,
        },
        "minimum_target_coverage_fraction": None,
        "derived_minimum_fragment_size": 1,
        "explicit_hydrogens": True,
        "saved_all_results": True,
        "elapsed_seconds": args.elapsed_seconds,
        "counts": dict(counts),
        "unresolved_source_ids": unresolved,
    }
    output.with_suffix(output.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

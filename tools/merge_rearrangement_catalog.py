#!/usr/bin/env python3
"""Merge and rank no-cut one-precursor rearrangement shards."""
from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parts", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--expected-id")
    args = parser.parse_args()

    paths = sorted(Path(args.parts).glob("part_*.jsonl.gz"))
    if not paths:
        raise RuntimeError("no rearrangement shards found")
    records = []
    scan_counts = Counter()
    for path in paths:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            records.extend(json.loads(line) for line in stream)
        summary_path = path.with_suffix(path.suffix + ".summary.json")
        if summary_path.exists():
            scan_counts.update(json.loads(summary_path.read_text())["counts"])
    records.sort(key=lambda item: (
        not item["complete"],
        tuple(item["score"]),
        len(item["excess_elements"]),
        len(item["broken"]),
        item["precursor_id"],
    ))
    known_rank = next(
        (index for index, item in enumerate(records, 1)
         if item["precursor_id"] == args.expected_id), None)
    report = {
        "schema": "rxn_core.rearrangement_catalog/v1",
        "scan_counts": dict(scan_counts),
        "result_count": len(records),
        "expected_id": args.expected_id,
        "expected_rank": known_rank,
        "expected_recovered": known_rank is not None,
        "results": records,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "output": str(output),
        "scan_counts": dict(scan_counts),
        "result_count": len(records),
        "expected_rank": known_rank,
    }, indent=2))


if __name__ == "__main__":
    main()

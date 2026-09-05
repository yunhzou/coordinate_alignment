#!/usr/bin/env python3
"""Audit completed writes and unfinished sources from scanner progress logs."""
import argparse
import csv
import gzip
import json
from pathlib import Path


def audit(run, catalog, id_column, continuations=()):
    completed = {}
    paths = [path for directory in (run, *continuations)
             for path in sorted((directory / "logs").glob("*.out"))]
    for path in paths:
        for line in path.read_text().splitlines():
            if not line.startswith('{"'):
                continue
            record = json.loads(line)
            if all(key in record for key in ('source_id', 'counts', 'output')):
                record = {**record, 'precursor_id': record['source_id'],
                          'pair_elapsed_seconds': record['elapsed_seconds']}
            if "precursor_id" in record:
                source = record["precursor_id"]
                if source in completed:
                    raise ValueError(f"source completed more than once: {source}")
                completed[source] = record
    with gzip.open(catalog, "rt") as stream:
        inventory = list(csv.DictReader(stream))
    missing = [{"source_id": row[id_column], "smiles": row["SMILES"], "row_index": i}
               for i, row in enumerate(inventory) if row[id_column] not in completed]
    summaries = [json.loads(path.read_text()) for path in (run / "parts").glob("*.summary.json")]
    all_sources_saved = set(completed) == {row[id_column] for row in inventory}
    initial_complete = bool(summaries) and len(summaries) == summaries[0]["shard_count"]
    return {
        "expected_rows": len(inventory), "saved_rows": len(completed),
        "completed_shards": len(summaries),
        "scan_complete": all_sources_saved and (initial_complete or bool(continuations)),
        "initial_scan_complete": all_sources_saved and initial_complete and not continuations,
        "continuation_directories": [str(path) for path in continuations],
        "unfinished_sources": missing,
        "slowest_saved_sources": sorted(completed.values(), key=lambda r: -r["pair_elapsed_seconds"])[:20],
        "timing_scope": "progress emitted only after complete compressed records are flushed; pair times include persistence",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--id-column", default="Inventory ID")
    parser.add_argument("--continuation", type=Path, action='append', default=[])
    parser.add_argument("--output", type=Path,
                        help="separate combined audit without changing the frozen resume input")
    args = parser.parse_args()
    report = audit(args.run, args.catalog, args.id_column, args.continuation)
    (args.output or args.run / "progress_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k not in ("unfinished_sources", "slowest_saved_sources")}))

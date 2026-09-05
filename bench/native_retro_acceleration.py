#!/usr/bin/env python3
"""Bounded same-host detection benchmark with full saved evidence."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

TARGET = "CC(C)C1=CC=CC(C(C)C)=C1/N=C2/C(C3=C4C(C=CC=C42)=CC=C3)=N/C5=C(C(C)C)C=CC=C5C(C)C"


def run_case(source_id, output):
    from rxn_core.fragment_matching import (
        detect_fragments, prepare_fragment_target, FragmentDetectionConfig)
    from rxn_core.fragment_matching.serialization import fragment_detection_to_record
    from rxn_core.smiles import smiles_to_weighted_graph
    with gzip.open("data/inventory/processed/inventory_structure_bank.csv.gz", "rt") as stream:
        row = next(r for r in csv.DictReader(stream) if r["Inventory ID"] == source_id)
    config = FragmentDetectionConfig(branch_limit=100, seed_mode="orbit_representatives")
    target = prepare_fragment_target(smiles_to_weighted_graph(TARGET, expand_hydrogens=True), config=config)
    source = smiles_to_weighted_graph(row["SMILES"], expand_hydrogens=True)
    started = time.perf_counter()
    result = detect_fragments(source, target, source_id=source_id, config=config)
    elapsed = time.perf_counter() - started
    record = fragment_detection_to_record(result, row_index=0, representation=row["SMILES"])
    output.mkdir(parents=True, exist_ok=True)
    with gzip.open(output / (source_id + ".json.gz"), "wt") as stream:
        json.dump(record, stream, separators=(",", ":"))
    # Archive references may change; compare the complete expanded evidence.
    normalized = {k: v for k, v in record.items()
                  if k not in ("schema", "hierarchy_fragments", "generators", "search_graphs", "candidates")}
    normalized["search_graphs"] = [g.to_record() for g in result.search_graphs]
    normalized["candidates"] = [dict(c, aam_hierarchy=typed.aam_hierarchy.to_record())
        for c, typed in zip(record["candidates"], result.candidates, strict=True)]
    digest = hashlib.sha256(json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    metrics = dict(source_id=source_id, detection_seconds=elapsed, candidates=len(result.candidates),
                   status=result.status, maximum_branch_count=result.maximum_branch_count,
                   full_evidence_sha256=digest)
    (output / (source_id + ".metrics.json")).write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-id", action="append")
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    ids = args.source_id or ["INVENTORY-001161", "INVENTORY-001301", "INVENTORY-000435"]
    if args.worker:
        run_case(ids[0], args.output)
    else:
        for source_id in ids:
            try:
                subprocess.run([sys.executable, __file__, "--worker", "--source-id", source_id,
                                "--output", str(args.output)], check=True, timeout=args.timeout)
            except subprocess.TimeoutExpired:
                print(json.dumps(dict(source_id=source_id, timeout_seconds=args.timeout)), flush=True)

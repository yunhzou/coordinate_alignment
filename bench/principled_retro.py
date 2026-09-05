#!/usr/bin/env python3
"""Saved end-to-end smoke experiment: methanol from a two-entry test bank.

This is a geometric regression, not a claim of chemical feasibility. Existing
detections are replayed; use a new output directory for a fresh timing run.
"""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys
import time

from rxn_core.fragment_matching import detect_fragments, FragmentDetectionConfig
from rxn_core.fragment_matching.serialization import fragment_detection_to_record, FRAGMENT_DETECTION_SCHEMA
from rxn_core.smiles import smiles_to_weighted_graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    directory = Path(args.output_dir).resolve()
    parts = directory / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    target_smiles = "CO"
    target = smiles_to_weighted_graph(target_smiles, expand_hydrogens=True)
    bank = (("bromomethane", "CBr"), ("water", "O"))
    timings = {}
    for index, (source_id, smiles) in enumerate(bank):
        path = parts / f"part_{index:03d}.jsonl.gz"
        summary_path = path.with_suffix(path.suffix + ".summary.json")
        if path.exists():
            with gzip.open(path, "rt") as stream:
                record = json.loads(stream.readline())
            if record["schema"] != FRAGMENT_DETECTION_SCHEMA:
                raise ValueError("choose a fresh directory for the current detection schema")
            timings[source_id] = json.loads(summary_path.read_text())["seconds"]
            continue
        started = time.perf_counter()
        result = detect_fragments(smiles_to_weighted_graph(smiles, expand_hydrogens=True),
            target, source_id=source_id, config=FragmentDetectionConfig(branch_limit=100))
        elapsed = time.perf_counter() - started
        record = fragment_detection_to_record(result, row_index=index, representation=smiles)
        with gzip.open(path, "wt") as stream:
            stream.write(json.dumps(record) + "\n")
        summary = {"schema": "rxn_core.retro_catalog_summary/v3", "target_smiles": target_smiles,
                   "seconds": elapsed, "counts": {"rows": 1, "searched": 1,
                       "matched_precursors": bool(result.candidates),
                       "fragment_candidates": len(result.candidates), "capped": not result.complete}}
        summary_path.write_text(json.dumps(summary) + "\n")
        timings[source_id] = elapsed
        print(json.dumps({"source": source_id, "seconds": elapsed,
                          "candidates": len(result.candidates), "status": result.status}), flush=True)
    tools = Path(__file__).resolve().parents[1] / "tools"
    output = directory / "methanol.json"
    subprocess.run([sys.executable, str(tools / "merge_retro_catalog.py"),
        "--parts", str(parts), "--target-smiles", target_smiles, "--output", str(output),
        "--pattern-limit", "2", "--recommendations-per-pattern", "1",
        "--expected-id", "bromomethane", "--expected-id", "water"], check=True)
    subprocess.run([sys.executable, str(tools / "build_retro_db_viewer.py"),
        "--results", str(output), "--output", str(directory / "methanol.html"),
        "--title", "Methanol: exact geometric assembly smoke test"], check=True)
    (directory / "timings.json").write_text(json.dumps(timings, indent=2) + "\n")
    print(json.dumps({"html": str(directory / "methanol.html"), "detection_seconds": timings}))


if __name__ == "__main__":
    main()

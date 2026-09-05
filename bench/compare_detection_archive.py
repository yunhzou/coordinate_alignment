#!/usr/bin/env python3
"""Time one saved source and compare full detection evidence, not just coverage."""
import argparse
import gzip
import json
from pathlib import Path
import time

from rxn_core.fragment_matching import detect_fragments, FragmentDetectionConfig
from rxn_core.fragment_matching.serialization import fragment_detection_to_record
from rxn_core.smiles import smiles_to_weighted_graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-part", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    with gzip.open(args.baseline_part, "rt") as stream:
        for line in stream:
            prior = json.loads(line)
            if prior["source_id"] == args.source_id:
                break
        else:
            raise ValueError("source not found in baseline")
    source = smiles_to_weighted_graph(prior["representation"], expand_hydrogens=True)
    target = smiles_to_weighted_graph(args.target_smiles, expand_hydrogens=True)
    start = time.perf_counter()
    result = detect_fragments(source, target, source_id=args.source_id,
        config=FragmentDetectionConfig(branch_limit=100, seed_mode="orbit_representatives"))
    detection_seconds = time.perf_counter() - start
    start = time.perf_counter()
    record = fragment_detection_to_record(result, row_index=prior["row_index"],
                                          representation=prior["representation"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / (args.source_id + ".json.gz")
    encoded = json.dumps(record, separators=(",", ":"))
    with gzip.open(output, "wt") as stream:
        stream.write(encoded)
    persistence_seconds = time.perf_counter() - start
    # Expand only for this differential oracle. Production keeps the archive shared.
    normalized = dict(record)
    normalized.pop("schema")
    normalized.pop("hierarchy_fragments")
    normalized.pop("generators")
    normalized["search_graphs"] = [g.to_record() for g in result.search_graphs]
    normalized["candidates"] = [dict(c, aam_hierarchy=typed.aam_hierarchy.to_record())
        for c, typed in zip(record["candidates"], result.candidates, strict=True)]
    prior.pop("schema")
    baseline_seconds = prior.pop("timing")["detection_seconds"]
    normalized = json.loads(json.dumps(normalized))
    changed = [key for key in prior if prior[key] != normalized[key]]
    metrics = {"source_id": args.source_id, "baseline_detection_seconds": baseline_seconds,
        "detection_seconds": detection_seconds, "persistence_seconds": persistence_seconds,
        "candidate_count": len(result.candidates), "maximum_branch_count": result.maximum_branch_count,
        "identical_full_evidence": not changed, "different_fields": changed,
        "baseline_json_bytes": len(line), "shared_json_bytes": len(encoded),
        "compressed_bytes": output.stat().st_size,
        "timing_note": "baseline cluster worker and current host may differ; correctness comparison is exact"}
    (args.output_dir / (args.source_id + ".metrics.json")).write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2), flush=True)
    assert not changed, changed


if __name__ == "__main__":
    main()

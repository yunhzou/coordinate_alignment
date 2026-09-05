#!/usr/bin/env python3
"""Compare saved full-bank shard timings; never repeat matching."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


def summarize(directory):
    summaries = [json.loads(p.read_text()) for p in sorted(
        (directory / "parts").glob("*.summary.json"))]
    if not summaries:
        return {"completed_shards": 0, "complete": False}
    counts = Counter()
    for summary in summaries:
        counts.update(summary["counts"])
    return {
        "completed_shards": len(summaries),
        "expected_shards": summaries[0]["shard_count"],
        "complete": len(summaries) == summaries[0]["shard_count"],
        "counts": dict(counts),
        "slowest_shard_seconds": max(s["elapsed_seconds"] for s in summaries),
        "sum_shard_seconds": sum(s["elapsed_seconds"] for s in summaries),
        "workers_per_shard": sorted({s["workers"] for s in summaries}),
        "target_smiles": summaries[0]["target_smiles"],
        "config": summaries[0]["config"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--job-id", required=True)
    args = parser.parse_args()
    baseline, current = summarize(args.baseline), summarize(args.run)
    report = {
        "revision": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "job_id": args.job_id, "catalog": str(args.catalog.resolve()),
        "catalog_sha256": hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
        "baseline_directory": str(args.baseline.resolve()),
        "run_directory": str(args.run.resolve()),
        "baseline": baseline, "current": current,
        "timing_scope": "scan including serialization; slowest shard excludes scheduler queue/startup and assembly",
        "policy_difference": "current detection has no candidate cap; baseline candidate cap was 100",
    }
    if current["complete"]:
        report["scan_slowdown"] = current["slowest_shard_seconds"] / baseline["slowest_shard_seconds"]
    slowest = []
    for summary_path in sorted((args.run / "parts").glob("*.summary.json")):
        path = Path(str(summary_path).removesuffix(".summary.json"))
        with gzip.open(path, "rt") as stream:
            for line in stream:
                record = json.loads(line)
                slowest.append({
                    "row_index": record["row_index"], "source_id": record["source_id"],
                    "seconds": record["timing"]["detection_seconds"],
                    "candidates": len(record["candidates"]),
                    "initial_families": record["initial_family_count"],
                    "status": record["status"],
                })
    report["slowest_completed_precursors"] = sorted(slowest, key=lambda r: -r["seconds"])[:20]
    output = args.run / "benchmark.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

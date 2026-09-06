#!/usr/bin/env python3
"""Summarize active bank writes without loading molecular result archives."""
import argparse
from collections import Counter
import json
from pathlib import Path
import statistics
import time


def summarize(run, expected_rows):
    saved, detected, elapsed = {}, {}, {}
    for path in (run / "logs").glob("*.out"):
        for line in path.read_text().splitlines():
            # Slurm streams may end in a not-yet-complete JSON log line.
            if not (line.startswith('{"') and line.endswith('}')):
                continue
            record = json.loads(line)
            if "precursor_id" in record:
                saved[record["precursor_id"]] = record
                elapsed[record["shard"]] = record["elapsed_seconds"]
            if record.get("stage") == "detection_checkpoint":
                detected[record["source_id"]] = record
    summaries = [json.loads(path.read_text()) for path in (run / "parts").glob("*.summary.json")
                 if path.exists()]
    counts = Counter()
    for summary in summaries:
        counts.update(summary["counts"])
    pair_times = [v["pair_elapsed_seconds"] for v in saved.values()]
    detection_times = [v["detection_seconds"] for v in detected.values()]
    report = {"measured_unix": time.time(), "expected_rows": expected_rows,
        "saved_rows": len(saved), "typed_checkpoint_events": len(detected),
        "completed_shards": len(summaries), "completed_shard_counts": dict(counts),
        "slowest_saved": sorted(saved.values(), key=lambda v: -v["pair_elapsed_seconds"])[:10],
        "pair_seconds_median": statistics.median(pair_times) if pair_times else None,
        "pair_seconds_max": max(pair_times, default=0),
        "detection_seconds_median": statistics.median(detection_times) if detection_times else None,
        "detection_seconds_max": max(detection_times, default=0),
        "latest_shard_elapsed": elapsed,
        "timing_scope": "Completed records only; excludes unfinished slow sources. Not a full-scan wall time."}
    (run / "live_progress.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int, required=True)
    args = parser.parse_args()
    report = summarize(args.run, args.expected_rows)
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ("slowest_saved", "latest_shard_elapsed")}))

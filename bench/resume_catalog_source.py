#!/usr/bin/env python3
"""Complete one audited missing source using the unchanged catalog detector."""
import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from search_mcule_retro import _worker_init, _search_one, _encode_records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
    args = parser.parse_args()
    audit = json.loads((args.prior_run / "progress_audit.json").read_text())
    source = audit["unfinished_sources"][args.index]
    prior_summary = sorted((args.prior_run / "parts").glob("*.summary.json"))[0]
    config = json.loads(prior_summary.read_text())
    output = args.output_dir / "parts" / f"part_{args.index}.jsonl.gz"
    checkpoint = args.output_dir / "parts" / f"part_{args.index}.checkpoint.json"
    if output.exists():
        raise FileExistsError(f"source result already saved: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    progress = {"source_id": source["source_id"], "row_index": source["row_index"],
        "phase": "running", "started_unix": time.time(), "workers": args.workers,
        "prior_run": str(args.prior_run.resolve()), "config": config["config"]}
    checkpoint.write_text(json.dumps(progress, indent=2) + "\n")
    print(json.dumps(progress), flush=True)
    _worker_init(config["target_smiles"], config["config"],
                 config["minimum_target_coverage_fraction"], True, None)
    counts, record = _search_one((source["row_index"], source["smiles"], source["source_id"]),
                                 seed_workers=max(1, args.workers - 1))
    with output.open("xb") as stream:
        stream.write(_encode_records([] if record is None else [record]))
    elapsed = time.perf_counter() - started
    summary = {**config, "output": str(output.resolve()), "counts": dict(counts),
        "shard_index": args.index, "shard_count": len(audit["unfinished_sources"]),
        "workers": args.workers, "scheduling": "resumed_one_source_per_node",
        "elapsed_seconds": elapsed, "source_id": source["source_id"],
        "resumed_from": str(args.prior_run.resolve()),
        "timing_scope": "fresh unfinished source only; not a fresh full-bank scan"}
    output.with_suffix(output.suffix + ".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    checkpoint.write_text(json.dumps({**progress, "phase": "saved", "elapsed_seconds": elapsed}, indent=2) + "\n")
    print(json.dumps({"source_id": source["source_id"], "elapsed_seconds": elapsed,
                      "counts": dict(counts), "output": str(output.resolve())}), flush=True)


if __name__ == "__main__":
    main()

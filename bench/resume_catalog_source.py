#!/usr/bin/env python3
"""Complete one audited missing source using the unchanged catalog detector."""
import argparse
import faulthandler
import gzip
import json
import os
from pathlib import Path
import signal
import pickle
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from search_mcule_retro import _worker_init, _search_one, _record_detection, _encode_records


def supervise(arguments, output_dir, index):
    """Wall-clock guard outside the detector, including native code and children."""
    with subprocess.Popen([sys.executable, __file__, *arguments, "--worker"],
                          start_new_session=True) as process:
        try:
            code = process.wait(timeout=300)
        except subprocess.TimeoutExpired:
            print(json.dumps({"event": "slow_run", "seconds": 300,
                "index": index, "scope": "precursor pipeline including AAM and post-processing"}), flush=True)
            try:
                code = process.wait(timeout=300)
            except subprocess.TimeoutExpired:
                # The isolated session contains the detector and its forked pool.
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                path = output_dir / "parts" / f"part_{index}.watchdog.json"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps({"phase": "timed_out", "limit_seconds": 600,
                    "result_complete": False, "index": index}) + "\n")
                raise SystemExit(124)
        if code:
            raise SystemExit(code)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior-run", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--resume-detection", action="store_true",
                        help="require and reuse the prior typed detection checkpoint")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not args.worker:
        return supervise(sys.argv[1:], args.output_dir, args.index)
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
    trace = (args.output_dir / 'parts' / f'part_{args.index}.stacks.txt').open('w')
    faulthandler.register(signal.SIGUSR1, file=trace)
    faulthandler.dump_traceback_later(60, repeat=True, file=trace)
    typed_directory = args.output_dir / 'checkpoints'
    typed_directory.mkdir(parents=True, exist_ok=True)
    _worker_init(config["target_smiles"], config["config"],
                 config["minimum_target_coverage_fraction"], True, None, typed_directory)
    row = source["row_index"], source["smiles"], source["source_id"]
    if args.resume_detection:
        saved = args.prior_run / 'checkpoints' / f'{row[0]}.detection.pkl.gz'
        with gzip.open(saved, 'rb') as stream:
            saved_row, detection = pickle.load(stream)
        if saved_row != row:
            raise ValueError('checkpoint source differs from audited source')
        counts, record = _record_detection(row, *detection)
    else:
        counts, record = _search_one(row, seed_workers=max(1, args.workers - 1))
    with output.open("xb") as stream:
        stream.write(_encode_records([] if record is None else [record]))
    elapsed = time.perf_counter() - started
    summary = {**config, "output": str(output.resolve()), "counts": dict(counts),
        "shard_index": args.index, "shard_count": len(audit["unfinished_sources"]),
        "workers": args.workers, "scheduling": "resumed_one_source_per_node",
        "elapsed_seconds": elapsed, "source_id": source["source_id"],
        "resumed_from": str(args.prior_run.resolve()),
        "reused_detection_checkpoint": args.resume_detection,
        "timing_scope": "fresh unfinished source only; not a fresh full-bank scan"}
    output.with_suffix(output.suffix + ".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    checkpoint.write_text(json.dumps({**progress, "phase": "saved", "elapsed_seconds": elapsed}, indent=2) + "\n")
    print(json.dumps({"source_id": source["source_id"], "elapsed_seconds": elapsed,
                      "counts": dict(counts), "output": str(output.resolve())}), flush=True)
    faulthandler.cancel_dump_traceback_later()
    faulthandler.unregister(signal.SIGUSR1)
    trace.close()


if __name__ == "__main__":
    main()

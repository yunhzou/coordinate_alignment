#!/usr/bin/env python3
"""Repack completed v4 detection shards into shared v6 archives without AAM."""
import argparse
import gzip
import json
from pathlib import Path

from rxn_core.fragment_matching.serialization import repack_fragment_detection_v4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError("choose a new output file; the original archive is preserved")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.input, "rt") as source, gzip.open(args.output, "xt") as target:
        for line in source:
            target.write(json.dumps(repack_fragment_detection_v4(json.loads(line)),
                                    separators=(",", ":")) + "\n")
    summary = args.input.with_suffix(args.input.suffix + ".summary.json")
    if summary.exists():
        record = json.loads(summary.read_text())
        record["output"] = str(args.output.resolve())
        record["repacked_from"] = str(args.input.resolve())
        args.output.with_suffix(args.output.suffix + ".summary.json").write_text(json.dumps(record) + "\n")


if __name__ == "__main__":
    main()

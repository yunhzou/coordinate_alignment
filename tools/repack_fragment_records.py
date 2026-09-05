#!/usr/bin/env python3
"""Explicitly migrate completed v4/v6 shards to shared v7 archives without AAM."""
import argparse
import gzip
import json
from pathlib import Path

from rxn_core.fragment_matching.serialization import (
    repack_fragment_detection_v4, repack_fragment_detection_v6)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--from-version", choices=("v4", "v6"), required=True)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError("choose a new output file; the original archive is preserved")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    repack = {"v4": repack_fragment_detection_v4, "v6": repack_fragment_detection_v6}[args.from_version]
    with gzip.open(args.input, "rt") as source, gzip.open(args.output, "xt") as target:
        for line in source:
            target.write(json.dumps(repack(json.loads(line)),
                                    separators=(",", ":")) + "\n")
    summary = args.input.with_suffix(args.input.suffix + ".summary.json")
    if summary.exists():
        record = json.loads(summary.read_text())
        record["output"] = str(args.output.resolve())
        record["repacked_from"] = str(args.input.resolve())
        args.output.with_suffix(args.output.suffix + ".summary.json").write_text(json.dumps(record) + "\n")


if __name__ == "__main__":
    main()

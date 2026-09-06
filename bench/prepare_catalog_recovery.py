#!/usr/bin/env python3
"""Prepare unchanged one-source recovery from specified interrupted shards."""
import argparse
import csv
import gzip
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--id-column", required=True)
    parser.add_argument("--shards", type=int, nargs="+", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary_path = next(path for path in sorted((args.run / "parts").glob("*.summary.json")) if path.exists())
    config = json.loads(summary_path.read_text())
    saved = set()
    for shard in args.shards:
        for path in (args.run / "logs").glob(f"*_{shard}.out"):
            for line in path.read_text().splitlines():
                if line.startswith('{"') and line.endswith('}'):
                    record = json.loads(line)
                    if "precursor_id" in record:
                        saved.add(record["precursor_id"])
    with gzip.open(args.catalog, "rt") as stream:
        missing = [{"source_id": row[args.id_column], "smiles": row["SMILES"], "row_index": index}
                   for index, row in enumerate(csv.DictReader(stream))
                   if index % config["shard_count"] in args.shards and row[args.id_column] not in saved]
    (args.output_dir / "parts").mkdir(parents=True)
    (args.output_dir / "parts/config.summary.json").write_text(json.dumps(config) + "\n")
    (args.output_dir / "checkpoints").symlink_to((args.run / "checkpoints").resolve(), target_is_directory=True)
    audit = {"scope": "Specified interrupted shards only, not a global bank audit",
             "shards": args.shards, "unfinished_sources": missing}
    (args.output_dir / "progress_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({**audit, "typed_checkpoints": [
        (args.output_dir / "checkpoints" / f"{row['row_index']}.detection.pkl.gz").exists()
        for row in missing]}), flush=True)


if __name__ == "__main__":
    main()

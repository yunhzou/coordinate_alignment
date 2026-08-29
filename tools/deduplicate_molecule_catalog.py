#!/usr/bin/env python3
"""Deduplicate a CSV molecule catalog by exact stereochemical structure."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi


def _open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def _split_values(value: str):
    return [part.strip() for part in value.split(";") if part.strip()]


def _append_unique(target: list[str], value: str):
    seen = set(target)
    for item in _split_values(value):
        if item not in seen:
            target.append(item)
            seen.add(item)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--smiles-column", default="SMILES")
    parser.add_argument("--compound-id-column", default="compound mcule ID")
    parser.add_argument("--product-id-column", default="product mcule IDs")
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    source = Path(args.input)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    records: dict[str, dict] = {}
    counts = Counter()

    with _open_text(source, "rt") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or args.smiles_column not in reader.fieldnames:
            raise ValueError(f"missing SMILES column: {args.smiles_column}")
        source_fields = list(reader.fieldnames)
        for row in reader:
            counts["rows"] += 1
            molecule = Chem.MolFromSmiles(row.get(args.smiles_column, ""))
            if molecule is None:
                counts["parse_errors"] += 1
                key = f"INVALID:{counts['rows']}"
                canonical = row.get(args.smiles_column, "")
            else:
                canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
                key = inchi.MolToInchiKey(molecule) or f"SMILES:{canonical}"
            if key not in records:
                records[key] = {
                    "row": dict(row),
                    "canonical_smiles": canonical,
                    "compound_ids": [],
                    "product_ids": [],
                    "source_rows": [],
                }
            else:
                counts["duplicate_rows"] += 1
            record = records[key]
            _append_unique(record["compound_ids"],
                           row.get(args.compound_id_column, ""))
            _append_unique(record["product_ids"],
                           row.get(args.product_id_column, ""))
            record["source_rows"].append(counts["rows"] - 1)

    extra_fields = ["canonical SMILES", "InChIKey", "duplicate count",
                    "all compound mcule IDs", "all product mcule IDs",
                    "source row indices"]
    with _open_text(output, "wt") as sink:
        writer = csv.DictWriter(sink, fieldnames=source_fields + extra_fields)
        writer.writeheader()
        for key, record in records.items():
            row = record["row"]
            row[args.smiles_column] = record["canonical_smiles"]
            row["canonical SMILES"] = record["canonical_smiles"]
            row["InChIKey"] = (key if not key.startswith(("INVALID:",
                                                           "SMILES:")) else "")
            row["duplicate count"] = len(record["source_rows"])
            row["all compound mcule IDs"] = ";".join(record["compound_ids"])
            row["all product mcule IDs"] = ";".join(record["product_ids"])
            row["source row indices"] = ";".join(
                str(index) for index in record["source_rows"])
            writer.writerow(row)

    summary = {
        "input": str(source),
        "output": str(output),
        "key": "full InChIKey (connectivity + stereochemistry)",
        "rows": counts["rows"],
        "unique_structures": len(records),
        "duplicate_rows_removed": counts["duplicate_rows"],
        "parse_errors": counts["parse_errors"],
    }
    output.with_suffix(output.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

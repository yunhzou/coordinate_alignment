#!/usr/bin/env python3
"""Build one stereochemically deduplicated fast-delivery inventory bank."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi


def _read_catalog(path, id_column):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
        for row in csv.DictReader(stream):
            yield row.get("SMILES", ""), row.get(id_column, "")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--class-1", required=True)
    parser.add_argument("--class-2", required=True)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    definitions = (
        ("fast_delivery_class_1", Path(args.class_1), "Mcule ID"),
        ("fast_delivery_class_2", Path(args.class_2), "Mcule ID"),
        ("inventory", Path(args.inventory), "Inventory ID"),
    )
    records = {}
    counts = Counter()
    for source, path, id_column in definitions:
        for smiles, compound_id in _read_catalog(path, id_column):
            counts[f"{source}_rows"] += 1
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                counts["parse_errors"] += 1
                continue
            canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
            key = inchi.MolToInchiKey(molecule) or f"SMILES:{canonical}"
            record = records.setdefault(key, {
                "smiles": canonical,
                "sources": [],
                "mcule_ids": [],
                "inventory_ids": [],
                "rows": 0,
            })
            record["rows"] += 1
            if source not in record["sources"]:
                record["sources"].append(source)
            ids = (record["inventory_ids"] if source == "inventory"
                   else record["mcule_ids"])
            if compound_id and compound_id not in ids:
                ids.append(compound_id)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "Bank ID", "SMILES", "Catalog Sources", "Mcule IDs",
        "Inventory IDs", "InChIKey", "Merged Row Count",
    )
    with gzip.open(output, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for key, record in records.items():
            identifiers = record["inventory_ids"] or record["mcule_ids"]
            writer.writerow({
                "Bank ID": identifiers[0] if identifiers else key,
                "SMILES": record["smiles"],
                "Catalog Sources": ";".join(record["sources"]),
                "Mcule IDs": ";".join(record["mcule_ids"]),
                "Inventory IDs": ";".join(record["inventory_ids"]),
                "InChIKey": key if not key.startswith("SMILES:") else "",
                "Merged Row Count": record["rows"],
            })

    memberships = Counter(
        ";".join(record["sources"]) for record in records.values())
    input_rows = {
        source: counts[f"{source}_rows"]
        for source, _path, _id_column in definitions
    }
    summary = {
        "schema": "rxn_core.merged_fast_delivery_with_inventory/v1",
        "inputs": {
            source: str(path) for source, path, _id_column in definitions
        },
        "input_rows": input_rows,
        "valid_input_rows": sum(input_rows.values()) - counts["parse_errors"],
        "unique_stereochemical_structures": len(records),
        "duplicate_rows_removed": (
            sum(input_rows.values()) - counts["parse_errors"] - len(records)),
        "source_membership": dict(sorted(memberships.items())),
        "parse_errors": counts["parse_errors"],
        "output": str(output),
    }
    output.with_suffix(output.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

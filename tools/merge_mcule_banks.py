#!/usr/bin/env python3
"""Merge deduplicated MCule catalogs with stereochemical deduplication."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import inchi


def _values(*fields):
    result = []
    seen = set()
    for field in fields:
        for value in str(field or "").split(";"):
            value = value.strip()
            if value and value not in seen:
                seen.add(value)
                result.append(value)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast-delivery", required=True)
    parser.add_argument("--natural-products", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")

    definitions = (
        ("fast_delivery", Path(args.fast_delivery), "Mcule ID",
         "all product mcule IDs"),
        ("natural_products", Path(args.natural_products),
         "compound mcule ID", "all product mcule IDs"),
    )
    records = {}
    counts = Counter()
    for source_name, path, id_column, product_column in definitions:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
            for row in csv.DictReader(stream):
                counts[f"{source_name}_rows"] += 1
                molecule = Chem.MolFromSmiles(row.get("SMILES", ""))
                if molecule is None:
                    counts["parse_errors"] += 1
                    continue
                canonical = Chem.MolToSmiles(
                    molecule, isomericSmiles=True)
                key = (row.get("InChIKey")
                       or inchi.MolToInchiKey(molecule)
                       or f"SMILES:{canonical}")
                record = records.setdefault(key, {
                    "smiles": canonical,
                    "ids": [],
                    "products": [],
                    "sources": [],
                    "rows": 0,
                })
                record["rows"] += 1
                record["ids"] = _values(
                    *record["ids"], row.get(id_column),
                    row.get("all compound mcule IDs"))
                record["products"] = _values(
                    *record["products"], row.get(product_column),
                    row.get("product mcule IDs"))
                if source_name not in record["sources"]:
                    record["sources"].append(source_name)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        "Mcule ID", "SMILES", "Catalog Sources", "Product Mcule IDs",
        "InChIKey", "Merged Row Count")
    with gzip.open(output, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for key, record in records.items():
            writer.writerow({
                "Mcule ID": record["ids"][0] if record["ids"] else key,
                "SMILES": record["smiles"],
                "Catalog Sources": ";".join(record["sources"]),
                "Product Mcule IDs": ";".join(record["products"]),
                "InChIKey": key if not key.startswith("SMILES:") else "",
                "Merged Row Count": record["rows"],
            })

    source_membership = Counter(
        ";".join(record["sources"]) for record in records.values())
    summary = {
        "schema": "rxn_core.merged_mcule_banks/v1",
        "inputs": {
            "fast_delivery": args.fast_delivery,
            "natural_products": args.natural_products,
        },
        "input_rows": {
            "fast_delivery": counts["fast_delivery_rows"],
            "natural_products": counts["natural_products_rows"],
        },
        "unique_stereochemical_structures": len(records),
        "source_membership": dict(source_membership),
        "cross_catalog_duplicates_removed": source_membership[
            "fast_delivery;natural_products"],
        "parse_errors": counts["parse_errors"],
        "output": str(output),
    }
    output.with_suffix(output.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

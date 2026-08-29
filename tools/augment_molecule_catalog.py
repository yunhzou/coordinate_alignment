#!/usr/bin/env python3
"""Add provenance-tracked molecules to a gzipped CSV catalog."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

from rdkit import Chem


def _structure(smiles: str):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    return (
        Chem.MolToSmiles(molecule, isomericSmiles=True),
        Chem.MolToInchiKey(molecule),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--additions", required=True,
                        help="Tab-separated SMILES, ID, and source label.")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    additions = []
    with Path(args.additions).open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError(
                    f"{args.additions}:{line_number}: expected 3 fields")
            smiles, compound_id, source = fields
            canonical, inchikey = _structure(smiles)
            additions.append({
                "smiles": smiles,
                "compound_id": compound_id,
                "source": source,
                "canonical": canonical,
                "inchikey": inchikey,
            })

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    existing_keys = set()
    added, duplicates = [], []
    row_count = 0
    with gzip.open(args.catalog, "rt", encoding="utf-8",
                   errors="replace") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames:
            raise ValueError("catalog has no CSV header")
        with gzip.open(output, "wt", encoding="utf-8") as sink:
            writer = csv.DictWriter(sink, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                writer.writerow(row)
                row_count += 1
                smiles = row.get("SMILES", "")
                molecule = Chem.MolFromSmiles(smiles)
                if molecule is not None:
                    existing_keys.add(Chem.MolToSmiles(
                        molecule, isomericSmiles=True))
            for addition in additions:
                if addition["canonical"] in existing_keys:
                    duplicates.append(addition)
                    continue
                row = {field: "" for field in reader.fieldnames}
                row["Mcule ID"] = addition["compound_id"]
                row["SMILES"] = addition["smiles"]
                if "canonical SMILES" in row:
                    row["canonical SMILES"] = addition["canonical"]
                if "InChIKey" in row:
                    row["InChIKey"] = addition["inchikey"]
                if "duplicate count" in row:
                    row["duplicate count"] = "1"
                if "all compound mcule IDs" in row:
                    row["all compound mcule IDs"] = addition["compound_id"]
                writer.writerow(row)
                existing_keys.add(addition["canonical"])
                added.append(addition)

    summary = {
        "schema": "rxn_core.augmented_molecule_catalog/v1",
        "base_catalog": args.catalog,
        "output": str(output),
        "base_rows": row_count,
        "output_rows": row_count + len(added),
        "added": added,
        "already_present": duplicates,
        "availability_note": (
            "Overlay entries retain blank fast-delivery availability fields."
        ),
    }
    summary_path = output.with_suffix(output.suffix + ".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

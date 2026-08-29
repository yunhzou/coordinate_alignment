#!/usr/bin/env python3
"""Write the exact-composition subset of a gzipped CSV molecule catalog."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem


def explicit_composition(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    return Counter(atom.GetSymbol() for atom in Chem.AddHs(molecule).GetAtoms())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--exclude-target-identity",
                        action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()

    target_composition = explicit_composition(args.target_smiles)
    if target_composition is None:
        raise ValueError("invalid target SMILES")
    target = Chem.MolFromSmiles(args.target_smiles)
    target_key = Chem.MolToSmiles(target, isomericSmiles=True)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    with gzip.open(args.catalog, "rt", encoding="utf-8", errors="replace") as source:
        reader = csv.DictReader(source)
        with gzip.open(output, "wt", encoding="utf-8") as sink:
            writer = csv.DictWriter(sink, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                counts["rows"] += 1
                composition = explicit_composition(row.get("SMILES", ""))
                if composition is None:
                    counts["parse_errors"] += 1
                elif (args.exclude_target_identity and
                      Chem.MolToSmiles(Chem.MolFromSmiles(row.get("SMILES", "")),
                                       isomericSmiles=True) == target_key):
                    counts["target_identity"] += 1
                elif composition == target_composition:
                    writer.writerow(row)
                    counts["selected"] += 1
                else:
                    counts["composition_rejected"] += 1
    summary = {
        "schema": "rxn_core.composition_catalog/v1",
        "source_catalog": args.catalog,
        "target_smiles": args.target_smiles,
        "target_composition": dict(sorted(target_composition.items())),
        "exclude_target_identity": args.exclude_target_identity,
        "counts": dict(counts),
    }
    output.with_suffix(output.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

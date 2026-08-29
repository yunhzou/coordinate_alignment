#!/usr/bin/env python3
"""Select strict target substructures without collapsing stereochemical variants."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--exclude-target-identity",
                        action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    target = Chem.MolFromSmiles(args.target_smiles)
    if target is None:
        raise ValueError("invalid target SMILES")
    target_key = Chem.MolToSmiles(target, isomericSmiles=True)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    coverages = []
    with gzip.open(args.catalog, "rt", encoding="utf-8",
                   errors="replace") as source:
        reader = csv.DictReader(source)
        with gzip.open(output, "wt", encoding="utf-8") as sink:
            writer = csv.DictWriter(sink, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                counts["rows"] += 1
                molecule = Chem.MolFromSmiles(row.get("SMILES", ""))
                if molecule is None:
                    counts["parse_errors"] += 1
                    continue
                if (args.exclude_target_identity and
                        Chem.MolToSmiles(molecule, isomericSmiles=True)
                        == target_key):
                    counts["target_identity"] += 1
                    continue
                if target.HasSubstructMatch(molecule, useChirality=False):
                    writer.writerow(row)
                    counts["selected"] += 1
                    coverages.append(
                        molecule.GetNumAtoms() / target.GetNumAtoms())
                else:
                    counts["rejected"] += 1
    summary = {
        "schema": "rxn_core.target_substructure_catalog/v1",
        "source_catalog": args.catalog,
        "target_smiles": args.target_smiles,
        "use_chirality": False,
        "stereochemical_variants_preserved": True,
        "exclude_target_identity": args.exclude_target_identity,
        "counts": dict(counts),
        "coverage_fraction": {
            "minimum": min(coverages) if coverages else None,
            "maximum": max(coverages) if coverages else None,
        },
    }
    output.with_suffix(output.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Select a target-similar percentile tier from a gzipped CSV catalog."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-smiles", required=True)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--retained-fraction", type=float)
    selection.add_argument("--selection-rule", choices=("sqrt",))
    parser.add_argument("--coverage-report")
    parser.add_argument("--exclude-catalog")
    args = parser.parse_args()
    if (args.retained_fraction is not None
            and not 0 < args.retained_fraction <= 1):
        raise ValueError("retained fraction must be in (0, 1]")

    target = Chem.MolFromSmiles(args.target_smiles)
    if target is None:
        raise ValueError("invalid target SMILES")
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048)
    additional = rdFingerprintGenerator.AdditionalOutput()
    additional.AllocateAtomToBits()
    target_fingerprint = generator.GetFingerprint(
        target, additionalOutput=additional)
    region_atoms = None
    if args.coverage_report:
        coverage = json.loads(Path(args.coverage_report).read_text())
        if coverage["target_smiles"] != args.target_smiles:
            raise ValueError("coverage report target does not match")
        region_atoms = {
            int(atom) for atom in coverage["union_uncovered_target_atoms"]
            if int(atom) < target.GetNumAtoms()
        }
        region_fingerprint = DataStructs.ExplicitBitVect(2048)
        atom_bits = additional.GetAtomToBits()
        for atom in region_atoms:
            for bit in atom_bits[atom]:
                region_fingerprint.SetBit(int(bit))
        target_fingerprint = region_fingerprint
    excluded_ids = set()
    if args.exclude_catalog:
        with gzip.open(args.exclude_catalog, "rt", encoding="utf-8") as stream:
            excluded_ids = {
                row["Inventory ID"] for row in csv.DictReader(stream)}
    rows = []
    parse_errors = 0
    with gzip.open(args.catalog, "rt", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        fieldnames = list(reader.fieldnames or ())
        for row_index, row in enumerate(reader):
            if row.get("Inventory ID") in excluded_ids:
                continue
            molecule = Chem.MolFromSmiles(row.get("SMILES", ""))
            if molecule is None:
                parse_errors += 1
                continue
            fingerprint = generator.GetFingerprint(molecule)
            similarity = DataStructs.TanimotoSimilarity(
                target_fingerprint, fingerprint)
            rows.append((
                -similarity,
                -molecule.GetNumHeavyAtoms(),
                row_index,
                row,
            ))
    rows.sort(key=lambda item: item[:3])
    retained_count = (
        math.ceil(math.sqrt(len(rows)))
        if args.selection_rule == "sqrt"
        else math.ceil(args.retained_fraction * len(rows)))
    selected = rows[:retained_count]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output_fields = fieldnames + ["Target Morgan Similarity"]
    with gzip.open(output, "wt", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=output_fields)
        writer.writeheader()
        for negative_similarity, _negative_size, _index, row in selected:
            writer.writerow(dict(
                row,
                **{"Target Morgan Similarity": -negative_similarity},
            ))

    summary = {
        "schema": "rxn_core.target_similarity_catalog/v1",
        "source_catalog": args.catalog,
        "target_smiles": args.target_smiles,
        "fingerprint": {"type": "Morgan", "radius": 2, "bits": 2048},
        "coverage_report": args.coverage_report,
        "region_target_atoms": sorted(region_atoms) if region_atoms else None,
        "excluded_catalog": args.exclude_catalog,
        "excluded_structures": len(excluded_ids),
        "selection_rule": args.selection_rule,
        "retained_fraction": retained_count / len(rows),
        "source_structures": len(rows),
        "selected_structures": len(selected),
        "parse_errors": parse_errors,
        "similarity_range": {
            "highest": -selected[0][0] if selected else None,
            "cutoff": -selected[-1][0] if selected else None,
        },
    }
    output.with_suffix(output.suffix + ".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

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
    parser.add_argument(
        "--coverage-region-field",
        choices=("best_partial_uncovered_target_atoms",
                 "union_uncovered_target_atoms"),
        default="union_uncovered_target_atoms")
    parser.add_argument("--exclude-catalog")
    parser.add_argument("--exclude-target-identity", action="store_true")
    parser.add_argument("--maximum-target-heavy-atom-fraction", type=float)
    args = parser.parse_args()
    if (args.retained_fraction is not None
            and not 0 < args.retained_fraction <= 1):
        raise ValueError("retained fraction must be in (0, 1]")
    if (args.maximum_target_heavy_atom_fraction is not None
            and not 0 < args.maximum_target_heavy_atom_fraction < 1):
        raise ValueError("maximum target heavy-atom fraction must be in (0, 1)")

    target = Chem.MolFromSmiles(args.target_smiles)
    if target is None:
        raise ValueError("invalid target SMILES")
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048)
    target_key = Chem.MolToSmiles(target, isomericSmiles=True)
    maximum_heavy_atoms = (
        math.floor(target.GetNumHeavyAtoms()
                   * args.maximum_target_heavy_atom_fraction)
        if args.maximum_target_heavy_atom_fraction is not None else None)
    additional = rdFingerprintGenerator.AdditionalOutput()
    additional.AllocateAtomToBits()
    target_fingerprint = generator.GetFingerprint(
        target, additionalOutput=additional)
    region_atoms = None
    requested_region_atoms = None
    if args.coverage_report:
        coverage = json.loads(Path(args.coverage_report).read_text())
        if coverage["target_smiles"] != args.target_smiles:
            raise ValueError("coverage report target does not match")
        requested_region_atoms = {
            int(atom) for atom in coverage[args.coverage_region_field]
        }
        explicit_target = Chem.AddHs(target)
        if any(atom < 0 or atom >= explicit_target.GetNumAtoms()
               for atom in requested_region_atoms):
            raise ValueError("coverage report contains invalid target atom")
        region_atoms = set()
        for atom in requested_region_atoms:
            if atom < target.GetNumAtoms():
                region_atoms.add(atom)
                continue
            region_atoms.update(
                neighbor.GetIdx()
                for neighbor in explicit_target.GetAtomWithIdx(
                    atom).GetNeighbors()
                if neighbor.GetIdx() < target.GetNumAtoms())
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
            if (args.exclude_target_identity
                    and Chem.MolToSmiles(
                        molecule, isomericSmiles=True) == target_key):
                continue
            if (maximum_heavy_atoms is not None
                    and molecule.GetNumHeavyAtoms() > maximum_heavy_atoms):
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
        "coverage_region_field": args.coverage_region_field,
        "requested_region_target_atoms": (
            sorted(requested_region_atoms)
            if requested_region_atoms is not None else None),
        "region_target_atoms": (
            sorted(region_atoms) if region_atoms is not None else None),
        "excluded_catalog": args.exclude_catalog,
        "excluded_target_identity": args.exclude_target_identity,
        "maximum_target_heavy_atom_fraction": (
            args.maximum_target_heavy_atom_fraction),
        "maximum_heavy_atoms": maximum_heavy_atoms,
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

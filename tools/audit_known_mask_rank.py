#!/usr/bin/env python3
"""Report known precursor positions within retained-target coverage masks."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from rdkit import Chem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parts", required=True)
    parser.add_argument("--known-id", action="append", required=True)
    args = parser.parse_args()
    known_ids = set(args.known_id)
    paths = sorted(Path(args.parts).glob("part_*.jsonl.gz"))
    known_rows = []
    for path in paths:
        with gzip.open(path, "rt") as stream:
            for line in stream:
                record = json.loads(line)
                if record["precursor_id"] not in known_ids:
                    continue
                molecule = Chem.MolFromSmiles(record["smiles"])
                explicit = Chem.AddHs(molecule)
                total = molecule.GetNumHeavyAtoms()
                for candidate in record["candidates"]:
                    mask = tuple(candidate["covered_target_atoms"])
                    retained = sum(
                        explicit.GetAtomWithIdx(atom).GetAtomicNum() > 1
                        for atom in candidate["retained_atoms"])
                    known_rows.append({
                        "precursor_id": record["precursor_id"],
                        "mask": mask,
                        "retained": retained,
                        "total": total,
                        "key": candidate_rank_key(
                            record, candidate, retained, total),
                        "better": 0,
                        "count": 0,
                    })

    masks = {item["mask"] for item in known_rows}
    for path in paths:
        with gzip.open(path, "rt") as stream:
            for line in stream:
                record = json.loads(line)
                molecule = Chem.MolFromSmiles(record["smiles"])
                if molecule is None:
                    continue
                total = molecule.GetNumHeavyAtoms()
                explicit = Chem.AddHs(molecule)
                for candidate in record["candidates"]:
                    mask = tuple(candidate["covered_target_atoms"])
                    if mask not in masks:
                        continue
                    retained = sum(
                        explicit.GetAtomWithIdx(atom).GetAtomicNum() > 1
                        for atom in candidate["retained_atoms"])
                    key = candidate_rank_key(record, candidate, retained, total)
                    for known in known_rows:
                        if known["mask"] != mask:
                            continue
                        known["count"] += 1
                        known["better"] += key < known["key"]

    output = []
    for known in known_rows:
        output.append({
            "precursor_id": known["precursor_id"],
            "covered_target_atom_count": len(known["mask"]),
            "retained_heavy_atoms": known["retained"],
            "total_heavy_atoms": known["total"],
            "mask_rank": known["better"] + 1,
            "mask_candidate_count": known["count"],
        })
    print(json.dumps(output, indent=2))


def candidate_rank_key(record, candidate, retained, total):
    return (
        record["status"] == "capped",
        -retained / total,
        len(candidate["boundary_bonds"]),
        sum(map(len, candidate["leftover_fragments"])),
        record["precursor_id"],
        candidate["mapping"],
    )


if __name__ == "__main__":
    main()

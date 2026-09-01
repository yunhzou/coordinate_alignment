#!/usr/bin/env python3
"""Find the best disjoint partial cover and the never-covered target region."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from rdkit import Chem


def _mask(atoms):
    value = 0
    for atom in atoms:
        value |= 1 << int(atom)
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", required=True, action="append")
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--maximum-precursors", type=int, default=12)
    parser.add_argument("--beam-width", type=int, default=20_000)
    args = parser.parse_args()

    target = Chem.AddHs(Chem.MolFromSmiles(args.target_smiles))
    records = []
    for path in args.records:
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            records.extend(json.loads(line) for line in stream)
    items = []
    union_mask = 0
    seen = set()
    for record in records:
        source = Chem.AddHs(Chem.MolFromSmiles(record["representation"]))
        for candidate in record["candidates"]:
            mask = _mask(candidate["covered_target_atoms"])
            union_mask |= mask
            key = (record["source_id"], mask)
            if key in seen:
                continue
            seen.add(key)
            items.append({
                "mask": mask,
                "source_id": record["source_id"],
                "smiles": record["representation"],
                "source_atom_count": source.GetNumAtoms(),
                "complete": record["complete"],
                "candidate": candidate,
            })

    def rank(selected, covered):
        retained = sum(item["mask"].bit_count() for item in selected)
        source_atoms = sum(item["source_atom_count"] for item in selected)
        retention = retained / source_atoms if source_atoms else 0
        return (
            -covered.bit_count(),
            len({item["source_id"] for item in selected}),
            -retention,
            sum(len(item["candidate"]["boundary_bonds"])
                for item in selected),
            tuple(item["source_id"] for item in selected),
        )

    states = {0: ()}
    for _depth in range(args.maximum_precursors):
        expanded = dict(states)
        for covered, selected in states.items():
            for item in items:
                if item["mask"] & covered:
                    continue
                new_covered = covered | item["mask"]
                new_selected = selected + (item,)
                current = expanded.get(new_covered)
                if (current is None
                        or rank(new_selected, new_covered)
                        < rank(current, new_covered)):
                    expanded[new_covered] = new_selected
        ranked = sorted(
            expanded.items(), key=lambda pair: rank(pair[1], pair[0]))
        states = dict(ranked[:args.beam_width])

    best_mask, best_items = min(
        states.items(), key=lambda pair: rank(pair[1], pair[0]))
    atom_count = target.GetNumAtoms()
    union_uncovered = [
        atom for atom in range(atom_count) if not union_mask & (1 << atom)]
    best_uncovered = [
        atom for atom in range(atom_count) if not best_mask & (1 << atom)]
    payload = {
        "schema": "rxn_core.partial_fragment_coverage/v1",
        "target_smiles": args.target_smiles,
        "target_atom_count": atom_count,
        "record_count": len(records),
        "candidate_count": len(items),
        "union_covered_target_atoms": [
            atom for atom in range(atom_count) if union_mask & (1 << atom)],
        "union_uncovered_target_atoms": union_uncovered,
        "best_partial_covered_target_atoms": [
            atom for atom in range(atom_count) if best_mask & (1 << atom)],
        "best_partial_uncovered_target_atoms": best_uncovered,
        "best_partial_precursors": [{
            "source_id": item["source_id"],
            "smiles": item["smiles"],
            "complete": item["complete"],
            "candidate": item["candidate"],
        } for item in best_items],
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "output": str(output),
        "union_coverage": atom_count - len(union_uncovered),
        "best_partial_coverage": atom_count - len(best_uncovered),
        "best_partial_precursors": [
            item["source_id"] for item in best_items],
    }, indent=2))


if __name__ == "__main__":
    main()

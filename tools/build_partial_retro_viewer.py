#!/usr/bin/env python3
"""Build the standard mapping viewer for one best partial target cover."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem

from build_retro_db_viewer import _html, _payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--partial-results", required=True)
    parser.add_argument("--scan-results", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--title", required=True)
    args = parser.parse_args()

    partial = json.loads(Path(args.partial_results).read_text())
    scan = json.loads(Path(args.scan_results).read_text())
    target = Chem.AddHs(Chem.MolFromSmiles(partial["target_smiles"]))
    raw_precursors = partial["best_partial_precursors"]
    precursors = []
    owner = {}
    relinquished = 0
    for module, item in enumerate(raw_precursors):
        candidate = item["candidate"]
        molecule = Chem.AddHs(Chem.MolFromSmiles(item["smiles"]))
        detected_coverage = candidate["covered_target_atoms"]
        covered = [atom for atom in detected_coverage if atom not in owner]
        relinquished += len(detected_coverage) - len(covered)
        covered_set = set(covered)
        retained_atoms = [
            source for source, target_atom in candidate["mapping"]
            if target_atom in covered_set
        ]
        mapping = [
            pair for pair in candidate["mapping"]
            if pair[1] in covered_set
        ]
        for atom in covered:
            owner[int(atom)] = module
        precursors.append({
            "precursor_id": item["source_id"],
            "smiles": item["smiles"],
            "complete": item["complete"],
            "covered_target_atoms": covered,
            "retained_atoms": retained_atoms,
            "leftover_fragments": candidate["leftover_fragments"],
            "boundary_bonds": candidate["boundary_bonds"],
            "attachment_atoms_target": candidate[
                "attachment_atoms_target"],
            "mapping": mapping,
            "retained_atom_count": len(retained_atoms),
            "total_atom_count": molecule.GetNumAtoms(),
        })
    formed = []
    for bond in target.GetBonds():
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if left in owner and right in owner and owner[left] != owner[right]:
            formed.append([left, right])
    retained = sum(item["retained_atom_count"] for item in precursors)
    total = sum(item["total_atom_count"] for item in precursors)
    score = {
        "chirality_violations": 0,
        "unique_precursor_structures": len({
            item["precursor_id"] for item in precursors}),
        "set_atom_retention": retained / total,
        "set_heavy_atom_retention": 0,
        "worst_heavy_atom_retention": 0,
        "mean_heavy_atom_retention": 0,
        "capped_precursors": sum(
            not item["complete"] for item in precursors),
        "broken_bonds": sum(
            len(item["boundary_bonds"]) for item in precursors),
        "leftover_atoms": sum(
            sum(map(len, item["leftover_fragments"]))
            for item in precursors),
        "formed_bonds": len(formed),
        "relinquished_overlap_atoms": relinquished,
    }
    report = {
        "target_smiles": partial["target_smiles"],
        "scan_counts": scan["scan_counts"],
        "construction_patterns": [{
            "pattern": 1,
            "coverage_atom_sets": [
                item["covered_target_atoms"] for item in precursors],
            "fragment_sizes": [
                len(item["covered_target_atoms"]) for item in precursors],
            "recommendation_ranks": [1],
            "best_score": score,
        }],
        "assemblies": [{
            "construction_pattern": 1,
            "precursors": precursors,
            "formed_bonds": formed,
            "score": score,
        }],
    }
    viewer_payload = _payload(
        report, 1, set(), {}, args.title,
        ("overlap-tolerant recursive cover"
         if partial.get("allow_overlap") else "partial geometric cover"))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_html(viewer_payload))
    print(output.resolve())


if __name__ == "__main__":
    main()

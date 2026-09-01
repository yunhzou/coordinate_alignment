#!/usr/bin/env python3
"""Run and persist the large repeated-arm geometric coverage experiment."""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path

from rdkit import Chem

from rxn_core.fragment_matching import FragmentDetectionConfig, detect_fragments
from rxn_core.fragment_matching.serialization import (
    fragment_detection_to_record,
)
from rxn_core.retrosynthesis import assemble_fragment_cover
from rxn_core.smiles import smiles_to_weighted_graph


CORE_SMILES = "Brc1cc(Br)cc(Br)c1"
ARM_SMILES = (
    "CC1(C)c2cc(B3OC(C)(C)C(C)(C)O3)ccc2-c2ccc"
    "(-c3ccc4c(c3)c3ccccc3n4-c3ccccc3)cc21"
)


def _coupling_scaffold(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    boron = next(
        atom for atom in molecule.GetAtoms() if atom.GetSymbol() == "B")
    attachment = next(
        atom for atom in boron.GetNeighbors()
        if atom.GetSymbol() == "C")
    attachment_index = attachment.GetIdx()
    editable = Chem.RWMol(molecule)
    editable.RemoveBond(boron.GetIdx(), attachment_index)
    atom_maps = []
    fragments = Chem.GetMolFrags(
        editable.GetMol(),
        asMols=True,
        sanitizeFrags=True,
        fragsMolAtomMapping=atom_maps,
    )
    return next(
        (fragment, indexes.index(attachment_index))
        for fragment, indexes in zip(fragments, atom_maps)
        if attachment_index in indexes
    )


def build_target_smiles():
    molecule = Chem.MolFromSmiles("c1ccccc1")
    for core_atom in (0, 2, 4):
        scaffold, attachment = _coupling_scaffold(ARM_SMILES)
        offset = molecule.GetNumAtoms()
        editable = Chem.RWMol(Chem.CombineMols(molecule, scaffold))
        editable.AddBond(
            core_atom, offset + attachment, Chem.BondType.SINGLE)
        molecule = editable.GetMol()
        Chem.SanitizeMol(molecule)
    return Chem.MolToSmiles(molecule)


def _graph(smiles):
    return smiles_to_weighted_graph(smiles, expand_hydrogens=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    target_smiles = build_target_smiles()
    target = _graph(target_smiles)
    config = FragmentDetectionConfig(
        minimum_fragment_size=1,
        iso_tolerance=0.5,
        branch_limit=100,
        candidate_limit=500,
    )
    sources = (
        ("1,3,5-tribromobenzene", CORE_SMILES),
        ("large OLED boronic ester", ARM_SMILES),
    )
    detections = []
    detection_records = []
    total_start = time.perf_counter()
    for row_index, (source_id, source_smiles) in enumerate(sources):
        start = time.perf_counter()
        detection = detect_fragments(
            _graph(source_smiles),
            target,
            source_id=source_id,
            config=config,
        )
        elapsed = time.perf_counter() - start
        detections.append(detection)
        record = fragment_detection_to_record(
            detection,
            row_index=row_index,
            representation=source_smiles,
        )
        record["elapsed_seconds"] = elapsed
        detection_records.append(record)
        print(
            source_id,
            f"{elapsed:.3f}s",
            f"candidates={len(detection.candidates)}",
            f"best={detection.best_fragment_size}",
            f"complete={detection.complete}",
            flush=True,
        )

    assembly_start = time.perf_counter()
    assembly_search = assemble_fragment_cover(
        target,
        tuple(
            candidate
            for detection in detections
            for candidate in detection.candidates
        ),
        maximum_precursors=4,
        assembly_limit=5_000,
        allow_repeated_precursors=True,
        require_attachment_bonds=False,
    )
    assembly_elapsed = time.perf_counter() - assembly_start
    assemblies = [
        {
            "precursor_ids": list(assembly.precursor_ids),
            "precursor_stoichiometry": dict(Counter(
                assembly.precursor_ids)),
            "formed_bonds": [list(bond) for bond in assembly.formed_bonds],
            "broken_bonds": [list(bond) for bond in assembly.broken_bonds],
            "candidate_mappings": [
                {
                    "source_id": candidate.source_id,
                    "mapping": [list(pair) for pair in candidate.mapping],
                    "covered_target_atoms": list(
                        candidate.covered_target_atoms),
                }
                for candidate in assembly.candidates
            ],
        }
        for assembly in assembly_search.assemblies
    ]
    payload = {
        "schema": "rxn_core.large_star_stress/v1",
        "target": {
            "smiles": target_smiles,
            "explicit_atom_count": len(target.nodes),
            "heavy_atom_count": Chem.MolFromSmiles(
                target_smiles).GetNumHeavyAtoms(),
        },
        "config": {
            "explicit_hydrogens": True,
            "isomorphism_tolerance": config.iso_tolerance,
            "branch_limit": config.branch_limit,
            "candidate_limit": config.candidate_limit,
            "maximum_precursors": 4,
        },
        "detections": detection_records,
        "assembly": {
            "status": assembly_search.status,
            "complete": assembly_search.complete,
            "elapsed_seconds": assembly_elapsed,
            "count": len(assemblies),
            "assemblies": assemblies,
        },
        "total_elapsed_seconds": time.perf_counter() - total_start,
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(output.resolve())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Benchmark one exact initial-fragment seed by its ordered position."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import time

from rdkit import Chem

from rxn_core.fragment_matching import FragmentDetectionConfig, prepare_fragment_target
from rxn_core.fragment_matching.detection import (
    _grow_initial_seed,
    _initial_seed_order,
    _prepare_fragment_detection,
)
from rxn_core.fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from rxn_core.matcher import _nauty_orbits


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--seed-position", type=int, required=True)
    parser.add_argument("--id-column", default="Inventory ID")
    args = parser.parse_args()

    with gzip.open(args.inventory, "rt", encoding="utf-8") as stream:
        row = next(
            row for row in csv.DictReader(stream)
            if row[args.id_column] == args.source_id)
    source_molecule = Chem.AddHs(Chem.MolFromSmiles(row["SMILES"]))
    target_molecule = Chem.AddHs(Chem.MolFromSmiles(args.target_smiles))
    config = FragmentDetectionConfig(
        iso_tolerance=0.5,
        branch_limit=100,
        candidate_limit=100,
        seed_mode="all",
    )
    target = prepare_fragment_target(
        molecule_to_weighted_graph(target_molecule), config=config)
    source, target_context, _ = _prepare_fragment_detection(
        molecule_to_weighted_graph(source_molecule), target, config)
    source_orbits = _nauty_orbits(source, wbo_tol=config.iso_tolerance)
    seed_order, _ = _initial_seed_order(source, config)
    seed = seed_order[args.seed_position]
    started = time.perf_counter()
    placements, capped, branch_count = _grow_initial_seed(
        source,
        target_context.graph,
        seed,
        config,
        source_orbits,
        target_context.atom_orbits,
    )
    print(json.dumps({
        "source_id": args.source_id,
        "seed_position": args.seed_position,
        "seed_atom": seed,
        "seconds": time.perf_counter() - started,
        "placement_count": len(placements),
        "capped": capped,
        "maximum_branch_count": branch_count,
    }))


if __name__ == "__main__":
    main()

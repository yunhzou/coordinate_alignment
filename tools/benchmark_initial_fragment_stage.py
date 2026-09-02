#!/usr/bin/env python3
"""Time exact parallel seed growth separately from family canonicalization."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import time

from rdkit import Chem

from rxn_core.fragment_matching import FragmentDetectionConfig, prepare_fragment_target
from rxn_core.fragment_matching.detection import _prepare_fragment_detection
from rxn_core.fragment_matching.parallel import _parallel_initial_fragment_placements
from rxn_core.fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from rxn_core.matcher.canonical import _PartialMappingCanonicalizer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--seed-limit", type=int)
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
        seed_limit=args.seed_limit,
        seed_mode="all",
    )
    target = prepare_fragment_target(
        molecule_to_weighted_graph(target_molecule), config=config)
    source, target_context, region = _prepare_fragment_detection(
        molecule_to_weighted_graph(source_molecule), target, config)

    original = _PartialMappingCanonicalizer.certificate
    certificate_calls = 0
    certificate_seconds = 0.0

    def measured(self, mapping):
        nonlocal certificate_calls, certificate_seconds
        started = time.perf_counter()
        result = original(self, mapping)
        certificate_calls += 1
        certificate_seconds += time.perf_counter() - started
        return result

    _PartialMappingCanonicalizer.certificate = measured
    started = time.perf_counter()
    result = _parallel_initial_fragment_placements(
        source,
        target_context.graph,
        config,
        target_context.atom_orbits,
        region,
        args.workers,
    )
    elapsed = time.perf_counter() - started
    _PartialMappingCanonicalizer.certificate = original
    print(json.dumps({
        "source_id": args.source_id,
        "workers": args.workers,
        "initial_stage_seconds": elapsed,
        "family_count": len(result[0]),
        "seed_attempt_count": result[5],
        "seed_pruned_count": result[6],
        "certificate_calls_in_parent": certificate_calls,
        "certificate_seconds_in_parent": certificate_seconds,
    }, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Benchmark one inventory precursor against one explicit-H target."""
from __future__ import annotations

import argparse
import cProfile
import csv
import gzip
import hashlib
import json
import pstats
import time

from rdkit import Chem

from rxn_core.fragment_matching import (
    FragmentDetectionConfig,
    FragmentDetectionExecution,
    detect_fragments,
    detect_fragments_parallel,
    prepare_fragment_target,
)
from rxn_core.fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from rxn_core.fragment_matching.serialization import fragment_candidate_to_record


def _inventory_row(path, source_id, row_index, id_column):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
        for index, row in enumerate(csv.DictReader(stream)):
            if ((source_id is not None and row[id_column] == source_id)
                    or (row_index is not None and index == row_index)):
                return row
    requested = source_id if source_id is not None else row_index
    raise ValueError(f"inventory source not found: {requested}")


def _result_digest(result):
    payload = {
        "status": result.status,
        "complete": result.complete,
        "maximum_branch_count": result.maximum_branch_count,
        "capped_seed_count": result.capped_seed_count,
        "best_fragment_size": result.best_fragment_size,
        "candidates": [
            fragment_candidate_to_record(candidate)
            for candidate in result.candidates
        ],
    }
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--inventory", required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--source-id")
    source.add_argument("--row-index", type=int)
    parser.add_argument("--id-column", default="Inventory ID")
    parser.add_argument("--seed-mode", choices=("all", "fragment_cover"),
                        default="all")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--profile")
    parser.add_argument("--profile-lines", type=int, default=40)
    parser.add_argument("--seed-workers", type=int, default=1)
    parser.add_argument("--seed-limit", type=int)
    parser.add_argument("--candidate-limit", type=int, default=100)
    args = parser.parse_args()
    if args.repeats < 1:
        raise ValueError("repeats must be positive")

    row = _inventory_row(
        args.inventory, args.source_id, args.row_index, args.id_column)
    source_id = row[args.id_column]
    source_molecule = Chem.AddHs(Chem.MolFromSmiles(row["SMILES"]))
    target_molecule = Chem.AddHs(Chem.MolFromSmiles(args.target_smiles))
    config = FragmentDetectionConfig(
        iso_tolerance=0.5,
        branch_limit=100,
        candidate_limit=args.candidate_limit,
        seed_limit=args.seed_limit,
        seed_mode=args.seed_mode,
        rough_retention_threshold=0.5,
    )
    source = molecule_to_weighted_graph(source_molecule)
    target_started = time.perf_counter()
    target = prepare_fragment_target(
        molecule_to_weighted_graph(target_molecule), config=config)
    target_seconds = time.perf_counter() - target_started

    profiler = cProfile.Profile() if args.profile else None
    durations = []
    results = []
    for _ in range(args.repeats):
        if profiler is not None:
            profiler.enable()
        started = time.perf_counter()
        if args.seed_workers == 1:
            result = detect_fragments(
                source,
                target,
                source_id=source_id,
                config=config,
            )
        else:
            result = detect_fragments_parallel(
                source,
                target,
                source_id=source_id,
                config=config,
                execution=FragmentDetectionExecution(
                    seed_workers=args.seed_workers),
            )
        durations.append(time.perf_counter() - started)
        if profiler is not None:
            profiler.disable()
        results.append(result)

    digests = [_result_digest(result) for result in results]
    if len(set(digests)) != 1:
        raise RuntimeError("benchmark repetitions produced different results")
    result = results[-1]
    print(json.dumps({
        "source_id": source_id,
        "source_explicit_atoms": source_molecule.GetNumAtoms(),
        "target_explicit_atoms": target_molecule.GetNumAtoms(),
        "seed_mode": args.seed_mode,
        "seed_workers": args.seed_workers,
        "target_preparation_seconds": target_seconds,
        "seconds": durations,
        "best_seconds": min(durations),
        "result_digest": digests[0],
        "status": result.status,
        "candidate_count": len(result.candidates),
        "best_fragment_size": result.best_fragment_size,
        "initial_placement_encounters": result.initial_placement_encounters,
        "initial_family_count": result.initial_family_count,
        "seed_attempt_count": result.seed_attempt_count,
        "seed_pruned_count": result.seed_pruned_count,
        "maximum_branch_count": result.maximum_branch_count,
    }, indent=2))
    if profiler is not None:
        profiler.dump_stats(args.profile)
        pstats.Stats(profiler).strip_dirs().sort_stats("cumulative").print_stats(
            args.profile_lines)


if __name__ == "__main__":
    main()

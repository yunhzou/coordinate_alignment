#!/usr/bin/env python3
"""Time exact source-orbit preparation for one inventory molecule."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import time

from rdkit import Chem

from rxn_core.fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from rxn_core.fragment_matching.detection import _coerce_graph
from rxn_core.matcher.orbits import (
    _nauty_colored_wbo_graph,
    _wbo_tolerance_bucket_lookup,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--id-column", default="Inventory ID")
    parser.add_argument("--iso-tolerance", type=float, default=0.5)
    args = parser.parse_args()

    with gzip.open(args.inventory, "rt", encoding="utf-8") as stream:
        row = next(
            row for row in csv.DictReader(stream)
            if row[args.id_column] == args.source_id)
    molecule = Chem.AddHs(Chem.MolFromSmiles(row["SMILES"]))
    started = time.perf_counter()
    graph = _coerce_graph(molecule_to_weighted_graph(molecule), 0.2)
    graph_seconds = time.perf_counter() - started

    started = time.perf_counter()
    pair_buckets, _ = _wbo_tolerance_bucket_lookup(
        graph, args.iso_tolerance)
    bucket_seconds = time.perf_counter() - started

    started = time.perf_counter()
    _, _, nauty_graph, _, _ = _nauty_colored_wbo_graph(
        graph, args.iso_tolerance)
    construction_seconds = time.perf_counter() - started

    import pynauty
    started = time.perf_counter()
    generators, _, _, orbits, _ = pynauty.autgrp(nauty_graph)
    automorphism_seconds = time.perf_counter() - started
    print(json.dumps({
        "source_id": args.source_id,
        "explicit_atoms": molecule.GetNumAtoms(),
        "graph_seconds": graph_seconds,
        "pair_bucket_seconds": bucket_seconds,
        "pair_bucket_count": len(pair_buckets),
        "nauty_graph_construction_seconds": construction_seconds,
        "nauty_automorphism_seconds": automorphism_seconds,
        "generator_count": len(generators),
        "orbit_count": len(set(orbits[:molecule.GetNumAtoms()])),
    }, indent=2))


if __name__ == "__main__":
    main()

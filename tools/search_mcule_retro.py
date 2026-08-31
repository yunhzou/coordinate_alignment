#!/usr/bin/env python3
"""Blind, sharded explicit-H retained-fragment search for one target SMILES."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import multiprocessing as mp
import os
import sys
import time
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

from rxn_core.fragment_matching import (
    FragmentDetectionConfig,
    detect_fragments,
)
from rxn_core.fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from rxn_core.fragment_matching.serialization import fragment_detection_to_record


_TARGET = None
_CONFIG = None
_MINIMUM_TARGET_COVERAGE_SIZE = 1


def _worker_init(target_smiles, config_record,
                 minimum_target_coverage_fraction):
    global _TARGET, _CONFIG, _MINIMUM_TARGET_COVERAGE_SIZE
    RDLogger.DisableLog("rdApp.*")
    target_implicit = Chem.MolFromSmiles(target_smiles)
    if target_implicit is None:
        raise ValueError(f"invalid target SMILES: {target_smiles!r}")
    target_molecule = Chem.AddHs(target_implicit)
    _TARGET = molecule_to_weighted_graph(target_molecule)
    _MINIMUM_TARGET_COVERAGE_SIZE = (
        max(1, math.ceil(minimum_target_coverage_fraction
                         * target_molecule.GetNumAtoms()))
        if minimum_target_coverage_fraction is not None else 1)
    _CONFIG = FragmentDetectionConfig(**config_record)


def _search_batch(batch):
    records = []
    counts = Counter(rows=len(batch))
    for row_index, smiles, precursor_id in batch:
        molecule_implicit = Chem.MolFromSmiles(smiles)
        if molecule_implicit is None:
            counts["parse_errors"] += 1
            continue
        molecule = Chem.AddHs(molecule_implicit)
        counts["searched"] += 1
        result = detect_fragments(
            molecule_to_weighted_graph(molecule),
            _TARGET,
            source_id=precursor_id,
            config=_CONFIG,
        )
        if result.status == "capped":
            counts["capped"] += 1
        candidates = [
            candidate for candidate in result.candidates
            if len(candidate.covered_target_atoms)
            >= _MINIMUM_TARGET_COVERAGE_SIZE
        ]
        if result.candidates and not candidates:
            counts["target_coverage_filtered"] += 1
        if not candidates:
            continue
        counts["matched_precursors"] += 1
        counts["fragment_candidates"] += len(candidates)
        records.append(fragment_detection_to_record(
            result,
            row_index=row_index,
            representation=smiles,
            candidates=candidates,
        ))
    return dict(counts), records


def _batches(catalog, shard_index, shard_count, batch_size, limit,
             catalog_format, id_column):
    batch = []
    accepted = 0
    with gzip.open(catalog, "rt", encoding="utf-8", errors="replace") as stream:
        if catalog_format == "csv":
            rows = (
                (row_index, row["SMILES"], row[id_column])
                for row_index, row in enumerate(csv.DictReader(stream))
            )
        else:
            def smi_rows():
                for row_index, line in enumerate(stream):
                    fields = line.rstrip("\n").split("\t", 1)
                    if len(fields) == 2:
                        yield row_index, fields[0], fields[1]
            rows = smi_rows()
        for row_index, smiles, precursor_id in rows:
            if row_index % shard_count != shard_index:
                continue
            if not smiles or not precursor_id:
                continue
            batch.append((row_index, smiles, precursor_id))
            accepted += 1
            if len(batch) >= batch_size:
                yield batch
                batch = []
            if limit is not None and accepted >= limit:
                break
    if batch:
        yield batch


def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--catalog-format", choices=("smi", "csv"),
                        default="smi")
    parser.add_argument("--id-column", default="Mcule ID")
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--minimum-fragment-size", type=int, default=1)
    parser.add_argument("--minimum-target-coverage-fraction", type=float)
    parser.add_argument("--iso-tolerance", type=float, default=0.5)
    parser.add_argument("--branch-limit", type=int, default=100)
    parser.add_argument("--candidate-limit", type=int, default=100)
    parser.add_argument("--maximum-boundary-bonds", type=int)
    parser.add_argument("--maximum-leftover-fragments", type=int)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("shard index must be in [0, shard count)")
    if (args.minimum_target_coverage_fraction is not None and
            not 0 < args.minimum_target_coverage_fraction <= 1):
        raise ValueError("target coverage fraction must be in (0, 1]")
    config_record = {
        "minimum_fragment_size": args.minimum_fragment_size,
        "iso_tolerance": args.iso_tolerance,
        "branch_limit": args.branch_limit,
        "candidate_limit": args.candidate_limit,
        "maximum_boundary_bonds": args.maximum_boundary_bonds,
        "maximum_leftover_fragments": args.maximum_leftover_fragments,
    }
    target_for_threshold = Chem.AddHs(Chem.MolFromSmiles(args.target_smiles))
    derived_minimum_fragment_size = (
        max(1, math.ceil(args.minimum_target_coverage_fraction
                         * target_for_threshold.GetNumAtoms()))
        if args.minimum_target_coverage_fraction is not None
        else args.minimum_fragment_size
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    totals = Counter()
    started = time.perf_counter()
    processed_batches = 0
    context = mp.get_context("fork")
    with gzip.open(output, "wt", encoding="utf-8") as sink:
        with context.Pool(
                processes=max(1, args.workers),
                initializer=_worker_init,
                initargs=(args.target_smiles, config_record,
                          args.minimum_target_coverage_fraction)) as pool:
            iterator = pool.imap_unordered(
                _search_batch,
                _batches(
                    args.catalog, args.shard_index, args.shard_count,
                    args.batch_size, args.limit, args.catalog_format,
                    args.id_column),
                chunksize=1,
            )
            for counts, records in iterator:
                totals.update(counts)
                for record in records:
                    sink.write(json.dumps(record, separators=(",", ":")) + "\n")
                processed_batches += 1
                if processed_batches % 100 == 0:
                    elapsed = time.perf_counter() - started
                    print(json.dumps({
                        "shard": args.shard_index,
                        "elapsed_seconds": round(elapsed, 3),
                        "rows": totals["rows"],
                        "searched": totals["searched"],
                        "matched": totals["matched_precursors"],
                        "capped": totals["capped"],
                        "rows_per_second": round(totals["rows"] / elapsed, 2),
                    }), flush=True)
    elapsed = time.perf_counter() - started
    count_record = {
        key: totals[key]
        for key in (
            "rows",
            "parse_errors",
            "searched",
            "capped",
            "target_coverage_filtered",
            "matched_precursors",
            "fragment_candidates",
        )
    }
    summary = {
        "schema": "rxn_core.retro_catalog_summary/v3",
        "target_smiles": args.target_smiles,
        "catalog": str(args.catalog),
        "output": str(output),
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "workers": args.workers,
        "config": config_record,
        "minimum_target_coverage_fraction": (
            args.minimum_target_coverage_fraction),
        "derived_minimum_fragment_size": derived_minimum_fragment_size,
        "explicit_hydrogens": True,
        "elapsed_seconds": elapsed,
        "counts": count_record,
    }
    summary_path = output.with_suffix(output.suffix + ".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

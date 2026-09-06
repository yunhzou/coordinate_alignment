#!/usr/bin/env python3
"""Blind, sharded explicit-H retained-fragment search for one target SMILES."""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import multiprocessing as mp
import os
import pickle
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path

from rdkit import Chem, RDLogger
from msgspec import json as fast_json

from rxn_core.fragment_matching import (
    FragmentDetectionConfig,
    FragmentDetectionExecution,
    detect_fragments,
    detect_fragments_parallel,
    prepare_fragment_target,
)
from rxn_core.fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from rxn_core.fragment_matching.serialization import fragment_detection_to_record
from rxn_core.retrosynthesis.config import DEFAULT_ISO_TOLERANCE


_TARGET = None
_CONFIG = None
_MINIMUM_TARGET_COVERAGE_SIZE = 1
_SAVE_ALL_RESULTS = False
_TARGET_REGION_ATOMS = None
_CHECKPOINT_DIRECTORY = None


def _worker_init(target_smiles, config_record,
                 minimum_target_coverage_fraction, save_all_results,
                 target_region_atoms, checkpoint_directory=None):
    global _TARGET, _CONFIG, _MINIMUM_TARGET_COVERAGE_SIZE
    global _SAVE_ALL_RESULTS
    global _TARGET_REGION_ATOMS
    global _CHECKPOINT_DIRECTORY
    RDLogger.DisableLog("rdApp.*")
    target_implicit = Chem.MolFromSmiles(target_smiles)
    if target_implicit is None:
        raise ValueError(f"invalid target SMILES: {target_smiles!r}")
    target_molecule = Chem.AddHs(target_implicit)
    _CONFIG = FragmentDetectionConfig(**config_record)
    _TARGET = prepare_fragment_target(
        molecule_to_weighted_graph(target_molecule), config=_CONFIG)
    _MINIMUM_TARGET_COVERAGE_SIZE = (
        max(1, math.ceil(minimum_target_coverage_fraction
                         * target_molecule.GetNumAtoms()))
        if minimum_target_coverage_fraction is not None else 1)
    _SAVE_ALL_RESULTS = bool(save_all_results)
    _TARGET_REGION_ATOMS = (
        frozenset(map(int, target_region_atoms))
        if target_region_atoms is not None else None)
    _CHECKPOINT_DIRECTORY = checkpoint_directory


def _detect_one(row, *, seed_workers=1):
    """Molecular search only; its typed result can be checkpointed and reused."""
    row_index, smiles, precursor_id = row
    counts = Counter(rows=1)
    molecule_implicit = Chem.MolFromSmiles(smiles)
    if molecule_implicit is None:
        counts["parse_errors"] += 1
        return counts, None, 0.0
    molecule = Chem.AddHs(molecule_implicit)
    counts["searched"] += 1
    source = molecule_to_weighted_graph(molecule)
    detection_started = time.perf_counter()
    if seed_workers == 1:
        result = detect_fragments(
            source,
            _TARGET,
            source_id=precursor_id,
            config=_CONFIG,
            target_region_atoms=_TARGET_REGION_ATOMS,
        )
    else:
        result = detect_fragments_parallel(
            source,
            _TARGET,
            source_id=precursor_id,
            config=_CONFIG,
            execution=FragmentDetectionExecution(
                seed_workers=seed_workers),
            target_region_atoms=_TARGET_REGION_ATOMS,
        )
    detection_seconds = time.perf_counter() - detection_started
    return counts, result, detection_seconds


def _record_detection(row, counts, result, detection_seconds):
    """Encode a completed detection without repeating any molecular search."""
    row_index, smiles, _precursor_id = row
    if result is None:
        return counts, None
    counts = Counter(counts)
    if result.status == "capped":
        counts["capped"] += 1
    counts["incomplete"] += not result.complete
    candidates = [
        candidate for candidate in result.candidates
        if len(candidate.covered_target_atoms)
        >= _MINIMUM_TARGET_COVERAGE_SIZE
    ]
    if len(candidates) != len(result.candidates):
        counts["target_coverage_filtered"] += 1
    if not candidates and not _SAVE_ALL_RESULTS:
        return counts, None
    if candidates:
        counts["matched_precursors"] += 1
        counts["fragment_candidates"] += len(candidates)
    record = fragment_detection_to_record(
        result,
        row_index=row_index,
        representation=smiles,
        candidates=candidates,
    )
    record["timing"] = {"detection_seconds": detection_seconds}
    return counts, record


def _search_one(row, *, seed_workers=1):
    detection = _detect_one(row, seed_workers=seed_workers)
    if _CHECKPOINT_DIRECTORY is not None:
        # Persist the complete typed result before archive construction. A
        # watchdog during later serialization must not require repeating AAM.
        path = Path(_CHECKPOINT_DIRECTORY) / f"{row[0]}.detection.pkl.gz"
        temporary = path.with_suffix('.partial')
        started = time.perf_counter()
        with gzip.open(temporary, 'xb', compresslevel=1) as stream:
            pickle.dump((row, detection), stream, protocol=5)
        temporary.replace(path)
        print(json.dumps({"stage": "detection_checkpoint", "source_id": row[2],
            "row_index": row[0], "detection_seconds": detection[2],
            "checkpoint_seconds": time.perf_counter() - started,
            "checkpoint": str(path)}), flush=True)
    return _record_detection(row, *detection)


def _search_batch(batch):
    records = []
    counts = Counter()
    for row in batch:
        row_counts, record = _search_one(row)
        counts.update(row_counts)
        if record is not None:
            records.append(record)
    return dict(counts), _encode_records(records)


def _encode_records(records):
    """Independent gzip members keep JSON encoding/compression on worker CPUs.

    Concatenated members form a standard gzip stream. The coordinator copies
    bytes, rather than unpickling and serially recompressing large AAM archives.
    """
    # Encode one record at a time; do not copy a whole batch into one huge string.
    output = io.BytesIO()
    with gzip.GzipFile(fileobj=output, mode="wb", compresslevel=1, mtime=0) as stream:
        for record in records:
            stream.write(fast_json.encode(record))
            stream.write(b"\n")
    return output.getvalue()


def _explicit_atom_count(row):
    molecule = Chem.MolFromSmiles(row[1])
    return 0 if molecule is None else Chem.AddHs(molecule).GetNumAtoms()


def _adaptive_partition(rows, workers):
    """Separate pairs whose seed work exceeds one worker's fair share."""
    weighted = [(row, _explicit_atom_count(row)) for row in rows]
    fair_share = sum(weight for _, weight in weighted) / max(1, workers)
    ordinary = [row for row, weight in weighted if weight <= fair_share]
    outliers = [
        row for row, weight in sorted(
            weighted, key=lambda item: (-item[1], item[0][0]))
        if weight > fair_share
    ]
    return ordinary, outliers, fair_share


def _outlier_worker_budgets(rows, workers):
    weights = [max(1, _explicit_atom_count(row)) for row in rows]
    remaining = max(0, workers - len(rows))
    total_weight = sum(weights)
    shares = [remaining * weight / total_weight for weight in weights]
    budgets = [1 + math.floor(share) for share in shares]
    unassigned = workers - sum(budgets)
    order = sorted(
        range(len(rows)),
        key=lambda index: (-(shares[index] - math.floor(shares[index])),
                           rows[index][0]),
    )
    for index in order[:unassigned]:
        budgets[index] += 1
    return tuple(budgets)


def _search_scheduled(row, worker_budget):
    started = time.perf_counter()
    seed_workers = 1 if worker_budget <= 2 else worker_budget - 1
    counts, record = _search_one(row, seed_workers=seed_workers)
    payload = _encode_records([] if record is None else [record])
    return row, counts, payload, time.perf_counter() - started, seed_workers


def _budgeted_results(executor, jobs, workers, function=_search_scheduled):
    """Admit independent jobs as CPU slots free; no phase-wide barriers."""
    pending = list(jobs)
    running = {}
    available = workers
    while pending or running:
        waiting = []
        for row, budget in pending:
            if not 1 <= budget <= workers:
                raise ValueError("job CPU budget outside allocation")
            if budget <= available:
                running[executor.submit(function, row, budget)] = budget
                available -= budget
            else:
                waiting.append((row, budget))
        pending = waiting
        completed, _ = wait(running, return_when=FIRST_COMPLETED)
        for future in completed:
            available += running.pop(future)
            yield future.result()


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
    parser.add_argument("--checkpoint-directory", type=Path,
                        help="save complete typed detections before archive construction")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--scheduling", choices=("precursor", "adaptive"),
        default="precursor")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--minimum-fragment-size", type=int, default=1)
    parser.add_argument("--minimum-target-coverage-fraction", type=float)
    parser.add_argument("--iso-tolerance", type=float, default=DEFAULT_ISO_TOLERANCE)
    parser.add_argument("--branch-limit", type=int, default=100)
    parser.add_argument("--candidate-limit", type=int)
    parser.add_argument("--seed-limit", type=int)
    parser.add_argument(
        "--seed-mode",
        choices=("all", "fragment_cover", "orbit_representatives"),
        default="all")
    parser.add_argument(
        "--rough-retention-threshold", type=float, default=0.5)
    parser.add_argument("--maximum-boundary-bonds", type=int)
    parser.add_argument("--maximum-leftover-fragments", type=int)
    parser.add_argument("--save-all-results", action=argparse.BooleanOptionalAction, default=True,
                        help="persist complete search evidence, including capped/no-match records")
    parser.add_argument("--target-region-report")
    parser.add_argument(
        "--target-region-field",
        choices=("best_partial_uncovered_target_atoms",
                 "union_uncovered_target_atoms"),
        default="best_partial_uncovered_target_atoms")
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
        "seed_limit": args.seed_limit,
        "seed_mode": args.seed_mode,
        "rough_retention_threshold": args.rough_retention_threshold,
        "maximum_boundary_bonds": args.maximum_boundary_bonds,
        "maximum_leftover_fragments": args.maximum_leftover_fragments,
    }
    target_region_atoms = None
    if args.target_region_report:
        region_report = json.loads(
            Path(args.target_region_report).read_text())
        if region_report["target_smiles"] != args.target_smiles:
            raise ValueError("target-region report target does not match")
        target_region_atoms = region_report[args.target_region_field]
    target_for_threshold = Chem.AddHs(Chem.MolFromSmiles(args.target_smiles))
    derived_minimum_fragment_size = (
        max(1, math.ceil(args.minimum_target_coverage_fraction
                         * target_for_threshold.GetNumAtoms()))
        if args.minimum_target_coverage_fraction is not None
        else args.minimum_fragment_size
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.checkpoint_directory is not None:
        args.checkpoint_directory.mkdir(parents=True, exist_ok=True)
    totals = Counter()
    started = time.perf_counter()
    processed_batches = 0
    outlier_count = 0
    fair_share_atoms = None
    raw_batches = _batches(
        args.catalog, args.shard_index, args.shard_count,
        args.batch_size, args.limit, args.catalog_format,
        args.id_column)
    outliers = []
    if args.scheduling == "adaptive":
        rows = [row for batch in raw_batches for row in batch]
        ordinary, outliers, fair_share_atoms = _adaptive_partition(
            rows, args.workers)
        outlier_count = len(outliers)
        batches = (
            ordinary[index:index + args.batch_size]
            for index in range(0, len(ordinary), args.batch_size)
        )
    else:
        batches = raw_batches

    def write_result(sink, counts, payload):
        nonlocal processed_batches
        totals.update(counts)
        sink.write(payload)
        sink.flush()
        processed_batches += 1

    context = mp.get_context("fork")
    with output.open("wb") as sink:
        # Also make the empty-bank output a valid gzip file.
        sink.write(_encode_records(()))
        if args.scheduling == "adaptive":
            print(json.dumps({
                "shard": args.shard_index,
                "scheduling": args.scheduling,
                "outlier_count": outlier_count,
                "fair_share_explicit_atoms": fair_share_atoms,
            }), flush=True)
        if args.scheduling == "adaptive":
            budgets = _outlier_worker_budgets(outliers, args.workers)
            jobs = list(zip(outliers, budgets)) + [(row, 1) for row in ordinary]
            with ProcessPoolExecutor(
                    max_workers=max(1, args.workers), mp_context=context,
                    initializer=_worker_init,
                    initargs=(args.target_smiles, config_record,
                              args.minimum_target_coverage_fraction,
                              args.save_all_results, target_region_atoms,
                              args.checkpoint_directory)) as executor:
                for row, counts, payload, pair_elapsed, seed_workers in _budgeted_results(
                        executor, jobs, max(1, args.workers)):
                    write_result(sink, counts, payload)
                    print(json.dumps({
                        "shard": args.shard_index, "rows": totals["rows"],
                        "precursor_id": row[2], "seed_workers": seed_workers,
                        "pair_elapsed_seconds": pair_elapsed,
                        "elapsed_seconds": time.perf_counter() - started,
                    }), flush=True)
        else:
            with context.Pool(
                processes=max(1, args.workers),
                initializer=_worker_init,
                initargs=(args.target_smiles, config_record,
                          args.minimum_target_coverage_fraction,
                          args.save_all_results, target_region_atoms,
                          args.checkpoint_directory),
            ) as pool:
                iterator = pool.imap_unordered(_search_batch, batches, chunksize=1)
                for counts, payload in iterator:
                    write_result(sink, counts, payload)
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
            "incomplete",
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
        "scheduling": args.scheduling,
        "adaptive_outlier_count": outlier_count,
        "adaptive_fair_share_explicit_atoms": fair_share_atoms,
        "config": config_record,
        "minimum_target_coverage_fraction": (
            args.minimum_target_coverage_fraction),
        "derived_minimum_fragment_size": derived_minimum_fragment_size,
        "explicit_hydrogens": True,
        "saved_all_results": args.save_all_results,
        "checkpoint_directory": str(args.checkpoint_directory) if args.checkpoint_directory else None,
        "target_region_report": args.target_region_report,
        "target_region_field": args.target_region_field,
        "target_region_atoms": target_region_atoms,
        "elapsed_seconds": elapsed,
        "counts": count_record,
    }
    summary_path = output.with_suffix(output.suffix + ".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

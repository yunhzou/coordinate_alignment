#!/usr/bin/env python3
"""Run an explicit-H one-edge cut sweep and report the best chiral result."""
from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

from rdkit import Chem

import search_mcule_rearrangement as search
from rxn_core.alignment.sweep import (
    cut_sweep, cut_sweep_items, run_cut_sweep_chunk)
from rxn_core.frag import classify_bonds


def _jsonable(value):
    """Convert the complete analytical AAM record to stable JSON values."""
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "item"):
        return value.item()
    return value


def _write_all_results(path, *, args, pool, evaluations, metrics):
    """Stream every mechanism and analytical branch without dropping data."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if output.suffix == ".gz" else open
    metadata = {
        "schema": "rxn_core.cut_sweep_aam/v1",
        "source_smiles": args.source_smiles,
        "target_smiles": args.target_smiles,
        "n_seeds_per_cut": args.n_seeds,
        "branch_limit": args.branch_limit,
        "only_cuts": args.only_cut,
        "selected_cuts": getattr(args, "selected_cuts", None),
        "heavy_cuts_only": args.heavy_cuts_only,
        "metrics": metrics,
    }
    with opener(output, "wt", encoding="utf-8") as handle:
        handle.write('{"metadata":')
        json.dump(_jsonable(metadata), handle, separators=(",", ":"))
        handle.write(',"mechanisms":[')
        for index, (signature, entry) in enumerate(pool.items()):
            if index:
                handle.write(",")
            record = {
                "signature": signature,
                "mapping": sorted(entry["mapping"].items()),
                "cuts": sorted(entry.get("cuts", ())),
                "has_no_cut": bool(entry.get("has_no_cut", False)),
                "dedup_count": int(entry.get("dedup_count", 1)),
                "branches": entry.get("branches", []),
                **evaluations[id(entry)],
            }
            json.dump(_jsonable(record), handle, separators=(",", ":"))
        handle.write("]}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-smiles", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--n-seeds", type=int, default=10)
    parser.add_argument("--branch-limit", type=int, default=100)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--only-cut", action="append", default=[], metavar="I-J",
        help="Run only this one-edge cut (repeatable), rather than the full sweep.")
    parser.add_argument("--cut-shard-index", type=int)
    parser.add_argument("--cut-shard-count", type=int)
    parser.add_argument("--heavy-cuts-only", action="store_true")
    parser.add_argument(
        "--best-result", help="Optional path for the full winning atom mapping.")
    parser.add_argument(
        "--all-results",
        help="Optional .json or .json.gz path for the complete AAM pool.")
    parser.add_argument(
        "--intermediate-dir",
        help="Persist raw and reduced signature buckets for restart/audit.")
    args = parser.parse_args()

    search._worker_init(
        args.target_smiles, 0, 0, 10_000, args.branch_limit, True, True)
    source = Chem.AddHs(Chem.MolFromSmiles(args.source_smiles))
    elements, wbo_source = search._mol_graph(source)
    started = time.perf_counter()
    common = dict(
        n_workers=args.workers, cut_floor=0.2, graph_floor=0.2,
        iso_tol=1.0, dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
        n_seeds=args.n_seeds, max_branches=args.branch_limit,
        chunksize=1, symmetry_repair=True,
        intermediate_dir=args.intermediate_dir)
    if args.cut_shard_index is not None:
        if not args.cut_shard_count or not (
                0 <= args.cut_shard_index < args.cut_shard_count):
            parser.error("cut shard index must be within cut shard count")
        all_cuts = cut_sweep_items(
            wbo_source, 0.2, elements=elements,
            heavy_only=args.heavy_cuts_only)
        cuts = all_cuts[args.cut_shard_index::args.cut_shard_count]
        args.selected_cuts = _jsonable(cuts)
        pool, metrics = run_cut_sweep_chunk(
            elements, wbo_source, search._TARGET_ELEMENTS,
            search._TARGET_WBO, cuts, return_metrics=True, **common)
        metrics.update({
            "all_cuts": len(all_cuts), "cuts": len(cuts),
            "seed_orders": len(cuts) * args.n_seeds,
            "shard_index": args.cut_shard_index,
            "shard_count": args.cut_shard_count,
        })
    elif args.only_cut:
        cuts = []
        for value in args.only_cut:
            left, right = value.split("-", 1)
            cuts.append(((int(left), int(right)),))
        pool, metrics = run_cut_sweep_chunk(
            elements, wbo_source, search._TARGET_ELEMENTS,
            search._TARGET_WBO, cuts, return_metrics=True, **common)
    else:
        common["heavy_cuts_only"] = args.heavy_cuts_only
        pool, metrics = cut_sweep(
            elements, wbo_source, search._TARGET_ELEMENTS,
            search._TARGET_WBO, return_metrics=True, **common)
    scored = []
    evaluations = {}
    for entry in pool.values():
        mapping = {int(a): int(b) for a, b in entry["mapping"].items()}
        broken, formed, _, _ = classify_bonds(
            mapping, wbo_source, search._TARGET_WBO,
            dwbo_threshold=0.5, elements_R=elements,
            elements_P=search._TARGET_ELEMENTS,
            metal_dwbo_threshold=0.3)
        chirality = search._chirality_violations(
            source, mapping, broken, formed)
        evaluations[id(entry)] = {
            "edit_count": len(broken) + len(formed),
            "chirality_violations": chirality,
            "broken": broken,
            "formed": formed,
        }
        if chirality == 0:
            scored.append((len(broken) + len(formed), broken, formed, entry))
    score, broken, formed, best = min(scored, key=lambda item: item[0])
    summary = {
        "n_seeds_per_cut": args.n_seeds,
        "branch_limit": args.branch_limit,
        "workers": args.workers,
        "elapsed_seconds": time.perf_counter() - started,
        "mechanism_count": len(pool),
        "zero_chirality_mechanism_count": len(scored),
        "best_zero_chirality_edits": score,
        "broken": len(broken),
        "formed": len(formed),
        "best_includes_uncut_search": bool(best.get("has_no_cut", False)),
        "best_discovery_cuts": repr(best.get("cuts")),
        "metrics": metrics,
    }
    if args.all_results:
        _write_all_results(
            args.all_results, args=args, pool=pool,
            evaluations=evaluations, metrics=metrics)
        summary["all_results"] = args.all_results
    if args.best_result:
        result = {
            "source_smiles": args.source_smiles,
            "target_smiles": args.target_smiles,
            "score": [score, 0],
            "chirality_violations": 0,
            "mapping": sorted([int(a), int(b)] for a, b in best["mapping"].items()),
            "broken": [list(event) for event in broken],
            "formed": [list(event) for event in formed],
            "has_no_cut": bool(best.get("has_no_cut", False)),
            "discovery_cuts": [list(edge) for edge in sorted(best.get("cuts", ()))],
        }
        output = Path(args.best_result)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2) + "\n")
        summary["best_result"] = str(output)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

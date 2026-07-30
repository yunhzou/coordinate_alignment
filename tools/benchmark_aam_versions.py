#!/usr/bin/env python3
"""Benchmark full AAM search and post-AAM finalization for one code tree."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import pickle
import resource
import subprocess
import sys
import time


def _digest(mapping):
    payload = json.dumps(
        sorted((int(r), int(p)) for r, p in mapping.items()),
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pool-output", type=Path)
    parser.add_argument("--pool-input", type=Path)
    args = parser.parse_args()

    source_root = args.source_root.resolve()
    sys.path.insert(0, str(source_root / "src"))
    from rxn_core import cut_sweep
    from rxn_core.alignment.index_chirality import fixed_mapping_aligned_rmsd
    from rxn_core.pipeline import (
        _rp_cut_kwargs,
        alignment_inputs_from_xyz,
        rp_stage_config,
        run_rp_stage_from_pool,
    )

    case_root = args.work_root.resolve() / args.case
    endpoint_root = case_root / "endpoints"
    inputs = alignment_inputs_from_xyz(
        endpoint_root / "R" / "reactant_combined.xyz",
        endpoint_root / "P" / "product_combined.xyz",
        reactant_workdir=endpoint_root / "R",
        product_workdir=endpoint_root / "P",
        xtb_mode="cache-only",
        name=args.case,
    )
    config = rp_stage_config()
    config.update({
        "index_chirality": "preserve",
        "search_mode": "full_cut_sweep",
        "n_seeds": 3,
        "max_branches": 100,
    })

    wall_started = time.perf_counter()
    aam_started = time.perf_counter()
    if args.pool_input:
        with args.pool_input.open("rb") as handle:
            pool = pickle.load(handle)
    else:
        pool = cut_sweep(
            inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
            n_workers=args.workers, **_rp_cut_kwargs(config))
    aam_seconds = time.perf_counter() - aam_started
    if args.pool_output:
        args.pool_output.parent.mkdir(parents=True, exist_ok=True)
        with args.pool_output.open("wb") as handle:
            pickle.dump(pool, handle, protocol=pickle.HIGHEST_PROTOCOL)
    post_started = time.perf_counter()
    result = run_rp_stage_from_pool(inputs, pool, config=config,
                                    elapsed=aam_seconds)
    post_seconds = time.perf_counter() - post_started

    mechanisms = []
    for mechanism in result.get("mechanisms") or ():
        mapping = {int(r): int(p)
                   for r, p in mechanism["mapping_RP"].items()}
        index = mechanism.get("index_chirality") or {}
        post = mechanism.get("post_aam") or {}
        mechanisms.append({
            "id": int(mechanism["id"]),
            "broken_bonds_R": mechanism.get("broken_bonds_R") or [],
            "formed_bonds_R": mechanism.get("formed_bonds_R") or [],
            "mapping_digest": _digest(mapping),
            "fixed_mapping_rmsd": fixed_mapping_aligned_rmsd(
                mapping, inputs.xyzR, inputs.xyzP),
            "index_chirality_violations": index.get(
                "selected_index_chirality_violation_count"),
            "rmsd_candidate_count": index.get("rmsd_candidate_count"),
            "rmsd_evaluated_leaf_count": index.get(
                "rmsd_evaluated_leaf_count"),
            "rmsd_pruned_leaf_count": index.get(
                "rmsd_pruned_leaf_count"),
            "rmsd_symmetry_factor_orders": index.get(
                "rmsd_symmetry_factor_orders"),
            "maximal_mapping_family_count": len(
                post.get("analytical_branches") or ()),
            "covered_path_counts": [
                branch.get("covered_path_count")
                for branch in post.get("analytical_branches") or ()
            ],
        })

    usage_self = resource.getrusage(resource.RUSAGE_SELF)
    usage_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        revision = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=source_root, text=True).strip()
    except Exception:
        revision = "unknown"
    record = {
        "source_root": str(source_root),
        "revision": revision,
        "case": args.case,
        "atom_count": len(inputs.elR),
        "workers": int(args.workers),
        "aam_seconds": aam_seconds,
        "post_aam_seconds": post_seconds,
        "total_seconds": time.perf_counter() - wall_started,
        "pool_mechanism_count": len(pool),
        "selected_mechanism_count": len(mechanisms),
        "mechanisms": mechanisms,
        "max_rss_self_kb": int(usage_self.ru_maxrss),
        "max_rss_children_kb": int(usage_children.ru_maxrss),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

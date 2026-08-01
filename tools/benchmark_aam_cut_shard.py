#!/usr/bin/env python3
"""Run or merge exact full-cut AAM shards for scaling experiments."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import time

from rxn_core.alignment.sweep import (
    cut_sweep_items,
    merge_cut_sweep_pools,
    run_cut_sweep_chunk,
)
from rxn_core.pipeline import alignment_inputs_from_xyz, rp_stage_config


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def run_shard(args):
    selection = read_json(args.selection_manifest)
    case = next(
        item for item in selection["cases"]
        if str(item["step_id"]) == str(args.step))
    run_root = args.run_root or selection.get("run_root")
    if run_root is None:
        raise ValueError("--run-root is required when absent from the manifest")
    endpoint_root = Path(run_root).resolve() / "work" / args.step / "endpoints"
    inputs = alignment_inputs_from_xyz(
        case["reactant_xyz"], case["product_xyz"], name=args.step,
        reactant_workdir=endpoint_root / "R",
        product_workdir=endpoint_root / "P", xtb_mode="cache-only")
    config = rp_stage_config()
    cuts = cut_sweep_items(inputs.wboR, float(config["cut_floor"]))
    selected = cuts[int(args.shard_index)::int(args.shard_count)]
    started = time.perf_counter()
    pool = run_cut_sweep_chunk(
        inputs.elR, inputs.wboR, inputs.elP, inputs.wboP, selected,
        n_workers=int(args.workers),
        cut_floor=float(config["cut_floor"]),
        graph_floor=float(config["graph_floor"]),
        iso_tol=float(config["iso_tol"]),
        dwbo_threshold=float(config["dwbo_threshold"]),
        metal_dwbo_threshold=config.get("metal_dwbo_threshold"),
        symmetry_wbo_tol=float(config["symmetry_wbo_tol"]),
        n_seeds=int(config["n_seeds"]),
        max_branches=int(config["max_branches"]),
        chunksize=int(config["chunksize"]),
        symmetry_repair=bool(config["symmetry_repair"]),
        symmetry_repair_min_changes=int(
            config["symmetry_repair_min_changes"]),
        symmetry_repair_max_evals=int(
            config["symmetry_repair_max_evals"]),
        anchor_map=config.get("anchor_map"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as handle:
        pickle.dump(pool, handle, protocol=pickle.HIGHEST_PROTOCOL)
    summary = {
        "step": args.step,
        "shard_index": int(args.shard_index),
        "shard_count": int(args.shard_count),
        "cut_count": len(selected),
        "workers": int(args.workers),
        "mechanism_count": len(pool),
        "branch_count": sum(
            len(entry.get("branches") or ()) for entry in pool.values()),
        "elapsed_seconds": time.perf_counter() - started,
        "output": str(args.output.resolve()),
    }
    write_json(args.summary, summary)
    print(json.dumps(summary), flush=True)


def merge_shards(args):
    started = time.perf_counter()
    pools = []
    for path in args.inputs:
        with path.open("rb") as handle:
            pools.append(pickle.load(handle))
    pool = merge_cut_sweep_pools(pools)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as handle:
        pickle.dump(pool, handle, protocol=pickle.HIGHEST_PROTOCOL)
    summary = {
        "input_count": len(pools),
        "mechanism_count": len(pool),
        "branch_count": sum(
            len(entry.get("branches") or ()) for entry in pool.values()),
        "elapsed_seconds": time.perf_counter() - started,
        "output": str(args.output.resolve()),
    }
    write_json(args.summary, summary)
    print(json.dumps(summary), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run")
    run.add_argument("--selection-manifest", type=Path, required=True)
    run.add_argument("--run-root", type=Path)
    run.add_argument("--step", required=True)
    run.add_argument("--shard-index", type=int, required=True)
    run.add_argument("--shard-count", type=int, required=True)
    run.add_argument("--workers", type=int, required=True)
    run.add_argument("--output", type=Path, required=True)
    run.add_argument("--summary", type=Path, required=True)
    merge = commands.add_parser("merge")
    merge.add_argument("--inputs", type=Path, nargs="+", required=True)
    merge.add_argument("--output", type=Path, required=True)
    merge.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "run":
        run_shard(args)
    else:
        merge_shards(args)


if __name__ == "__main__":
    main()

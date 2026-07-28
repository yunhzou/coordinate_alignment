#!/usr/bin/env python3
"""Run one selected case's full R/P cut sweep with structured tracing."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from rxn_core.alignment.sweep import cut_sweep_items, run_cut_sweep_chunk
from rxn_core.pipeline import alignment_inputs_from_xyz, rp_stage_config


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-manifest", type=Path, required=True)
    parser.add_argument("--step", required=True)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--workers", type=int, required=True)
    args = parser.parse_args()

    selection = read_json(args.selection_manifest.resolve())
    matches = [case for case in selection["cases"]
               if str(case["step_id"]) == args.step]
    if len(matches) != 1:
        raise ValueError(
            f"expected one selected case for {args.step!r}, found {len(matches)}")
    source = matches[0]
    endpoint_cache = (
        Path(selection["run_root"]).resolve()
        / "work" / args.step / "endpoints"
    )
    inputs = alignment_inputs_from_xyz(
        source["reactant_xyz"], source["product_xyz"], name=args.step,
        reactant_workdir=endpoint_cache / "R",
        product_workdir=endpoint_cache / "P",
        xtb_mode="cache-only",
    )
    config = rp_stage_config()
    cuts = cut_sweep_items(inputs.wboR, float(config["cut_floor"]))
    args.trace.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.trace.unlink(missing_ok=True)

    started = time.perf_counter()
    result = run_cut_sweep_chunk(
        inputs.elR, inputs.wboR, inputs.elP, inputs.wboP, cuts,
        n_workers=int(args.workers), trace_path=args.trace,
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
        symmetry_repair_min_changes=int(config["symmetry_repair_min_changes"]),
        symmetry_repair_max_evals=int(config["symmetry_repair_max_evals"]),
        anchor_map=config.get("anchor_map"),
    )
    summary = {
        "step": args.step,
        "atoms": len(inputs.elR),
        "cuts": len(cuts),
        "workers": int(args.workers),
        "mechanism_signatures": len(result),
        "elapsed_seconds": time.perf_counter() - started,
        "trace": str(args.trace.resolve()),
    }
    args.summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

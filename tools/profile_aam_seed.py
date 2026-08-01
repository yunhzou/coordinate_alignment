#!/usr/bin/env python3
"""Profile one exact cut/seed AAM growth path."""
from __future__ import annotations

import argparse
import cProfile
import io
import json
from pathlib import Path
import pstats
import time

from rxn_core.alignment.branch import _generate_seed_orders
from rxn_core.alignment.sweep import (
    _MechanismEventCanonicalizer,
    _nauty_orbits,
    _run_cut_work,
)
from rxn_core.frag import build_graph
from rxn_core.pipeline import alignment_inputs_from_xyz, rp_stage_config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-manifest", type=Path, required=True)
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--step", required=True)
    parser.add_argument("--cut", required=True, help="left,right")
    parser.add_argument("--seed-index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    selection = json.loads(args.selection_manifest.read_text())
    case = next(item for item in selection["cases"]
                if str(item["step_id"]) == args.step)
    run_root = args.run_root or selection.get("run_root")
    if run_root is None:
        raise ValueError("--run-root is required when absent from the manifest")
    endpoint_root = Path(run_root) / "work" / args.step / "endpoints"
    inputs = alignment_inputs_from_xyz(
        case["reactant_xyz"], case["product_xyz"], name=args.step,
        reactant_workdir=endpoint_root / "R",
        product_workdir=endpoint_root / "P", xtb_mode="cache-only")
    cfg = rp_stage_config()
    cfg["n_atoms"] = len(inputs.elR)
    graph_floor = float(cfg["graph_floor"])
    g_P = build_graph(inputs.elP, inputs.wboP, bond_cut=graph_floor)
    g_R_full = build_graph(inputs.elR, inputs.wboR, bond_cut=graph_floor)
    cut = (tuple(map(int, args.cut.split(","))),)
    g_R = g_R_full.copy()
    for left, right in cut:
        if g_R.has_edge(left, right):
            g_R.remove_edge(left, right)
    orders = _generate_seed_orders(g_R, int(cfg["n_seeds"]))
    p_orbits = _nauty_orbits(g_P, wbo_tol=float(cfg["iso_tol"]))
    r_orbits = _nauty_orbits(g_R_full, wbo_tol=float(cfg["iso_tol"]))
    # Construct once so its initialization appears outside the profile, as it
    # does in normal cut work after future worker caching.
    _MechanismEventCanonicalizer(
        g_R_full, wbo_tol=float(cfg["iso_tol"]))
    profiler = cProfile.Profile()
    started = time.perf_counter()
    profiler.enable()
    results, _events, _metrics = _run_cut_work(
        inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
        cfg, cut, [orders[int(args.seed_index)]], None,
        g_P, g_R_full, p_orbits, r_orbits,
        return_trace=False, collect_metrics=False)
    profiler.disable()
    elapsed = time.perf_counter() - started
    stream = io.StringIO()
    pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats(
        "cumulative").print_stats(80)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(stream.getvalue(), encoding="utf-8")
    print(json.dumps({
        "step": args.step, "cut": cut, "seed_index": args.seed_index,
        "elapsed_seconds": elapsed, "result_count": len(results),
        "profile": str(args.output.resolve()),
    }), flush=True)


if __name__ == "__main__":
    main()

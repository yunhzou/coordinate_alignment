#!/usr/bin/env python3
"""Benchmark full AAM search and post-AAM finalization for one code tree."""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import subprocess
import sys
import threading
import time


def _digest(mapping):
    payload = json.dumps(
        sorted((int(r), int(p)) for r, p in mapping.items()),
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def _process_children(pid):
    children = set()
    task_root = Path(f"/proc/{int(pid)}/task")
    try:
        tasks = list(task_root.iterdir())
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return children
    for task in tasks:
        try:
            raw = (task / "children").read_text(encoding="ascii").split()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        children.update(map(int, raw))
    return children


def _process_rss_kb(pid):
    try:
        for line in Path(f"/proc/{int(pid)}/status").read_text(
                encoding="ascii").splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1])
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        pass
    return 0


class _ProcessTreeMemorySampler:
    """Sample aggregate resident memory for this process and descendants."""

    def __init__(self, interval=0.1):
        self.interval = float(interval)
        self.peak_rss_kb = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _sample(self):
        pending = [os.getpid()]
        seen = set()
        while pending:
            pid = pending.pop()
            if pid in seen:
                continue
            seen.add(pid)
            pending.extend(_process_children(pid) - seen)
        total = sum(_process_rss_kb(pid) for pid in seen)
        self.peak_rss_kb = max(self.peak_rss_kb, total)

    def _run(self):
        while not self._stop.wait(self.interval):
            self._sample()

    def start(self):
        self._sample()
        self._thread.start()
        return self

    def stop(self):
        self._stop.set()
        self._thread.join()
        self._sample()
        return int(self.peak_rss_kb)


class _CallLedger:
    """Temporarily count and time selected post-AAM function calls."""

    def __init__(self):
        self.calls = {}
        self.seconds = {}
        self._originals = []

    def instrument(self, module, attribute, label=None):
        label = str(label or attribute)
        original = getattr(module, attribute)

        @functools.wraps(original)
        def measured(*args, **kwargs):
            started = time.perf_counter()
            self.calls[label] = self.calls.get(label, 0) + 1
            try:
                return original(*args, **kwargs)
            finally:
                self.seconds[label] = (
                    self.seconds.get(label, 0.0)
                    + time.perf_counter() - started)

        self._originals.append((module, attribute, original))
        setattr(module, attribute, measured)
        return self

    def restore(self):
        while self._originals:
            module, attribute, original = self._originals.pop()
            setattr(module, attribute, original)


def _pool_metrics(pool):
    entries = list(dict(pool or {}).values())
    branch_counts = [len(entry.get("branches") or ()) for entry in entries]
    branches = [branch for entry in entries
                for branch in entry.get("branches") or ()]
    fragments = [fragment for branch in branches
                 for fragment in ((branch.get("branch_symmetry") or {})
                                  .get("fragments") or ())]
    return {
        "pool_branch_count": sum(branch_counts),
        "max_branches_per_mechanism": max(branch_counts, default=0),
        "pool_fragment_record_count": len(fragments),
        "pool_symmetry_block_count": sum(
            len(((fragment.get("symmetry") or {}).get("blocks") or ()))
            for fragment in fragments),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pool-output", type=Path)
    parser.add_argument("--pool-input", type=Path)
    parser.add_argument("--contract", type=Path)
    parser.add_argument("--performance-profile")
    parser.add_argument("--fail-on-regression", action="store_true")
    args = parser.parse_args()

    source_root = args.source_root.resolve()
    sys.path.insert(0, str(source_root / "src"))
    from rxn_core import cut_sweep
    from rxn_core.benchmark_regression import evaluate_record, load_contract
    import rxn_core.alignment.sweep as sweep_module
    import rxn_core.pipeline as pipeline_module
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

    memory_sampler = _ProcessTreeMemorySampler().start()
    wall_started = time.perf_counter()
    aam_started = time.perf_counter()
    if args.pool_input:
        with args.pool_input.open("rb") as handle:
            pool = pickle.load(handle)
        aam_metrics = {}
    else:
        pool, aam_metrics = cut_sweep(
            inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
            n_workers=args.workers, return_metrics=True,
            **_rp_cut_kwargs(config))
    aam_seconds = time.perf_counter() - aam_started
    if args.pool_output:
        args.pool_output.parent.mkdir(parents=True, exist_ok=True)
        with args.pool_output.open("wb") as handle:
            pickle.dump(pool, handle, protocol=pickle.HIGHEST_PROTOCOL)
    post_started = time.perf_counter()
    ledger = _CallLedger()
    ledger.instrument(
        pipeline_module, "complete_chosen_automorphism_groups")
    ledger.instrument(
        pipeline_module, "_nauty_orbits", "pipeline._nauty_orbits")
    ledger.instrument(
        sweep_module, "_nauty_orbits", "sweep._nauty_orbits")
    ledger.instrument(
        sweep_module, "_nauty_atom_generators",
        "sweep._nauty_atom_generators")
    try:
        result = run_rp_stage_from_pool(inputs, pool, config=config,
                                        elapsed=aam_seconds)
    finally:
        ledger.restore()
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
    peak_process_tree_rss_kb = memory_sampler.stop()
    regression_metrics = {
        **aam_metrics,
        "configured_max_branches": int(config["max_branches"]),
        "stage_call_counts": {
            "aam_search": 1,
            "post_aam_finalize": 1,
        },
        "post_aam_call_counts": dict(sorted(ledger.calls.items())),
        "post_aam_call_seconds": dict(sorted(ledger.seconds.items())),
        **_pool_metrics(pool),
    }
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
        "stage_timing": result.get("timing") or {},
        "mechanisms": mechanisms,
        "peak_process_tree_rss_kb": peak_process_tree_rss_kb,
        "max_rss_self_kb": int(usage_self.ru_maxrss),
        "max_rss_children_kb": int(usage_children.ru_maxrss),
        "regression_metrics": regression_metrics,
    }
    if args.contract:
        record["regression_contract"] = evaluate_record(
            record, load_contract(args.contract),
            performance_profile=args.performance_profile)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record), flush=True)
    if (args.fail_on_regression
            and not (record.get("regression_contract") or {}).get(
                "passed", False)):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Stream diagnostics while running exactly one R/P no-cut work unit.

This tool wraps hot functions for observation only.  It does not change their
arguments, return values, matching order, limits, or exception behavior.
"""
from __future__ import annotations

import argparse
import faulthandler
import json
import os
from pathlib import Path
import threading
import time
import traceback

from rxn_core.pipeline import alignment_inputs_from_xyz, rp_stage_config
from rxn_core.alignment.sweep import run_cut_sweep_chunk


class StreamLog:
    def __init__(self, path: Path, heartbeat_seconds: float) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = path.open("a", encoding="utf-8", buffering=1)
        self.lock = threading.Lock()
        self.started = time.perf_counter()
        self.heartbeat_seconds = float(heartbeat_seconds)
        self.stop = threading.Event()
        self.state = {
            "phase": "initializing",
            "extend_calls": 0,
            "extend_seconds": 0.0,
            "boundary_calls": 0,
            "boundary_seconds": 0.0,
            "repair_calls": 0,
            "repair_seconds": 0.0,
        }

    def emit(self, event: str, **fields) -> None:
        record = {
            "event": event,
            "elapsed_sec": time.perf_counter() - self.started,
            **fields,
        }
        with self.lock:
            self.handle.write(json.dumps(record, sort_keys=True) + "\n")

    def set_phase(self, phase: str, **fields) -> None:
        with self.lock:
            self.state["phase"] = phase
            self.state["phase_fields"] = fields

    def add_timing(self, name: str, elapsed: float) -> None:
        with self.lock:
            self.state[f"{name}_calls"] += 1
            self.state[f"{name}_seconds"] += float(elapsed)

    def heartbeat(self) -> None:
        while not self.stop.wait(self.heartbeat_seconds):
            with self.lock:
                snapshot = dict(self.state)
            self.emit("heartbeat", **snapshot)

    def close(self) -> None:
        self.stop.set()
        self.handle.close()


def _automorph_domain_count(cands) -> int:
    return sum(
        len(getattr(cand, "automorph_blocks", ()) or ()) for cand in cands)


def install_observers(log: StreamLog) -> None:
    # Patch the names actually imported by each hot module.  Every wrapper is
    # observational: it forwards the original arguments and result unchanged.
    from rxn_core.growth import island
    from rxn_core.matcher import dedupe
    from rxn_core.alignment import branch, sweep

    original_extend = island._extend_sym_cands

    def observed_extend(cands, fragment, n, *args, **kwargs):
        call_no = log.state["extend_calls"] + 1
        fields = {
            "call": call_no,
            "fragment_size": len(fragment),
            "candidate_count": len(cands),
            "automorph_domain_count": _automorph_domain_count(cands),
            "new_r_atom": int(n),
        }
        log.set_phase("extend_sym_cands", **fields)
        log.emit("extend_start", **fields)
        started = time.perf_counter()
        try:
            result = original_extend(cands, fragment, n, *args, **kwargs)
        except BaseException as exc:
            log.emit("extend_error", **fields, error_type=type(exc).__name__,
                     error=str(exc), operation_sec=time.perf_counter() - started)
            raise
        elapsed = time.perf_counter() - started
        log.add_timing("extend", elapsed)
        log.emit(
            "extend_end", **fields, operation_sec=elapsed,
            result_count=len(result or ()),
            result_alternate_count=_alternate_count(result or ()),
        )
        return result

    island._extend_sym_cands = observed_extend

    def patch_boundary(module, label: str) -> None:
        if not hasattr(module, "_boundary_signature"):
            return
        original = module._boundary_signature

        def observed(*args, **kwargs):
            fields = {"site": label, "call": log.state["boundary_calls"] + 1}
            log.set_phase("boundary_signature", **fields)
            started = time.perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                elapsed = time.perf_counter() - started
                log.add_timing("boundary", elapsed)
                if elapsed >= 0.1:
                    log.emit("slow_boundary", **fields, operation_sec=elapsed)

        module._boundary_signature = observed

    patch_boundary(island, "growth.island")
    patch_boundary(dedupe, "matcher.dedupe")
    patch_boundary(branch, "alignment.branch")

    original_repair = sweep.symmetry_repair_mapping

    def observed_repair(mapping, *args, **kwargs):
        fields = {
            "call": log.state["repair_calls"] + 1,
            "mapped_atoms": len(mapping),
        }
        log.set_phase("symmetry_repair", **fields)
        log.emit("repair_start", **fields)
        started = time.perf_counter()
        try:
            result = original_repair(mapping, *args, **kwargs)
        except BaseException as exc:
            log.emit("repair_error", **fields, error_type=type(exc).__name__,
                     error=str(exc), operation_sec=time.perf_counter() - started)
            raise
        elapsed = time.perf_counter() - started
        log.add_timing("repair", elapsed)
        stats = result[1] if kwargs.get("return_stats") and isinstance(result, tuple) else None
        log.emit(
            "repair_end", **fields, operation_sec=elapsed,
            evaluated=(stats or {}).get("evaluated"),
            capped=(stats or {}).get("capped"),
        )
        return result

    sweep.symmetry_repair_mapping = observed_repair


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-manifest", type=Path, required=True)
    parser.add_argument("--step", required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--heartbeat-seconds", type=float, default=30.0)
    parser.add_argument("--stack-seconds", type=float, default=60.0)
    args = parser.parse_args()

    selection = read_json(args.selection_manifest.resolve())
    matches = [case for case in selection.get("cases", ())
               if str(case.get("step_id")) == args.step]
    if len(matches) != 1:
        raise ValueError(f"expected one selected case for {args.step!r}, found {len(matches)}")
    source = matches[0]
    source_run_root = Path(selection["run_root"]).resolve()
    endpoint_cache = source_run_root / "work" / args.step / "endpoints"

    log = StreamLog(args.log.resolve(), args.heartbeat_seconds)
    thread = threading.Thread(target=log.heartbeat, daemon=True)
    thread.start()
    faulthandler.enable()
    faulthandler.dump_traceback_later(args.stack_seconds, repeat=True)
    try:
        inputs = alignment_inputs_from_xyz(
            source["reactant_xyz"], source["product_xyz"], name=args.step,
            reactant_workdir=endpoint_cache / "R",
            product_workdir=endpoint_cache / "P",
            xtb_mode="cache-only",
        )
        config = rp_stage_config()
        log.emit(
            "case_start", step=args.step, atoms=len(inputs.elR), pid=os.getpid(),
            config={key: config.get(key) for key in (
                "cut_floor", "graph_floor", "iso_tol", "dwbo_threshold",
                "metal_dwbo_threshold", "n_seeds", "max_branches",
                "symmetry_repair", "symmetry_repair_min_changes",
                "symmetry_repair_max_evals",
            )},
        )
        install_observers(log)
        log.set_phase("run_no_cut")
        result = run_cut_sweep_chunk(
            inputs.elR, inputs.wboR, inputs.elP, inputs.wboP, [()],
            n_workers=1,
            trace_path=args.trace.resolve(),
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
        log.emit(
            "case_end", mechanisms=len(result),
            results=[{
                "signature": repr(signature),
                "mapping": {
                    str(int(r)): int(p)
                    for r, p in sorted((entry.get("mapping") or {}).items())
                },
            } for signature, entry in result.items()],
        )
        return 0
    except BaseException as exc:
        log.emit("case_error", error_type=type(exc).__name__, error=str(exc),
                 traceback=traceback.format_exc())
        raise
    finally:
        faulthandler.cancel_dump_traceback_later()
        log.close()


if __name__ == "__main__":
    raise SystemExit(main())

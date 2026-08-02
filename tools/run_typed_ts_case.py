#!/usr/bin/env python3
"""Run one cached benchmark case through the typed R/P and TS APIs."""
from __future__ import annotations

import argparse
import json
import os
import re
import time
from pathlib import Path

from rxn_core import (
    AAMProblem,
    AAMSearchConfig,
    MolecularEndpoint,
    TransitionStateTarget,
    VibrationalModes,
    align_reaction,
    analyze_transition_state,
    rp_record,
    ts_record,
    write_rp_bundle,
)
from rxn_core.chemistry_computations import parse_xyz
from rxn_core.chemistry_computations.xtb import read_wbo_file
from rxn_core.modes import parse_g98_modes


def _endpoint(cache, label):
    cache = Path(cache)
    xyz_files = sorted(
        path for path in cache.glob("*.xyz")
        if not path.name.startswith("xtb"))
    if len(xyz_files) != 1:
        raise RuntimeError(
            f"expected one source XYZ in {cache}, found {len(xyz_files)}")
    elements, coordinates = parse_xyz(xyz_files[0])
    wbo = read_wbo_file(cache / "wbo", len(elements))
    return MolecularEndpoint(
        tuple(elements), coordinates, wbo, label=label)


def _case_from_manifest(path, index):
    document = json.loads(Path(path).read_text())
    cases = document["cases"] if isinstance(document, dict) else document
    if index < 0 or index >= len(cases):
        raise IndexError(f"case index {index} outside 0..{len(cases)-1}")
    return cases[index]


def _iteration(path):
    match = re.search(r"iter(\d+)", str(path))
    return int(match.group(1)) if match else -1


def run_case(case, *, work_root, benchmark_root, output_root, workers,
             post_workers=None):
    case_id = str(case["step_id"])
    started = time.perf_counter()
    case_work = Path(work_root) / case_id
    case_output = Path(output_root) / "cases" / case_id
    case_output.mkdir(parents=True, exist_ok=True)
    config = AAMSearchConfig()

    reactant = _endpoint(case_work / "endpoints" / "R", "R")
    product = _endpoint(case_work / "endpoints" / "P", "P")
    rp_started = time.perf_counter()
    rp = align_reaction(
        AAMProblem(reactant, product, name=case_id),
        search_config=config,
        workers=max(1, int(workers)),
        post_workers=max(1, int(post_workers or workers)),
    )
    rp_seconds = time.perf_counter() - rp_started
    write_rp_bundle(rp, case_output / "rp")

    initial_guess = Path(benchmark_root) / case_id / "initial_guess"
    guess_files = {
        _iteration(path): path for path in initial_guess.glob("*.xyz")}
    target_records = []
    for iteration in range(1, 21):
        target_started = time.perf_counter()
        hessian_cache = case_work / "targets" / f"iter{iteration}_hess"
        missing = []
        if iteration not in guess_files:
            missing.append("initial_guess_xyz")
        if not (hessian_cache / "g98.out").is_file():
            missing.append("g98.out")
        if not (hessian_cache / "wbo").is_file():
            missing.append("wbo")
        if not list(hessian_cache.glob("*.xyz")):
            missing.append("cached_xyz")
        if missing:
            target_records.append({
                "iteration": iteration,
                "status": "missing",
                "missing": missing,
                "elapsed_seconds": time.perf_counter() - target_started,
            })
            continue
        molecule = _endpoint(hessian_cache, f"iter{iteration}")
        frequencies, modes = parse_g98_modes(hessian_cache / "g98.out")
        target = TransitionStateTarget(
            molecule, VibrationalModes(frequencies, modes), kind="initial_guess")
        result = analyze_transition_state(
            rp, target, search_config=config)
        record = ts_record(result)
        record.update({
            "iteration": iteration,
            "status": "ok",
            "initial_guess": str(guess_files[iteration]),
            "hessian_cache": str(hessian_cache),
            "elapsed_seconds": time.perf_counter() - target_started,
        })
        target_records.append(record)

    rankings = []
    for mechanism_index in range(len(rp.mechanisms)):
        scored = []
        for target in target_records:
            if target.get("status") != "ok":
                continue
            selected = target["mechanisms"][mechanism_index]["selected"]
            if selected is not None:
                scored.append({
                    "iteration": target["iteration"],
                    "score": selected["score"],
                    "overlap": selected["overlap"],
                    "wbo_progress": selected["wbo_progress"],
                    "frequency": selected["frequency"],
                    "sources": selected["sources"],
                })
        scored.sort(key=lambda item: (-item["score"], item["iteration"]))
        rankings.append({
            "mechanism_id": mechanism_index + 1,
            "ranked_initial_guesses": scored,
        })

    document = {
        "schema": "rxn_core.typed_ts_case/v1",
        "case": case_id,
        "atom_count": reactant.atom_count,
        "workers": int(workers),
        "post_workers": int(post_workers or workers),
        "rp_seconds": rp_seconds,
        "elapsed_seconds": time.perf_counter() - started,
        "rp": rp_record(rp),
        "targets": target_records,
        "rankings": rankings,
    }
    (case_output / "ts_scores.json").write_text(json.dumps(
        document, indent=2, default=float))
    summary = {
        "status": "ok",
        "case": case_id,
        "atom_count": reactant.atom_count,
        "mechanism_count": len(rp.mechanisms),
        "target_count": len(target_records),
        "scored_target_count": sum(
            target.get("status") == "ok" for target in target_records),
        "missing_target_count": sum(
            target.get("status") == "missing" for target in target_records),
        "rp_seconds": rp_seconds,
        "elapsed_seconds": document["elapsed_seconds"],
        "scores": str(case_output / "ts_scores.json"),
        "rp_view": str(case_output / "rp" / "view.html"),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
    }
    (case_output / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--index", type=int, default=None)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--benchmark-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--post-workers", type=int, default=None)
    args = parser.parse_args(argv)
    index = args.index
    if index is None:
        index = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    case = _case_from_manifest(args.manifest, index)
    run_case(
        case,
        work_root=args.work_root,
        benchmark_root=args.benchmark_root,
        output_root=args.output_root,
        workers=args.workers,
        post_workers=args.post_workers,
    )


if __name__ == "__main__":
    main()

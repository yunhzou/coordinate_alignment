#!/usr/bin/env python3
"""Run the current exact index-chirality R/P sweep over a tiered manifest."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import time
import traceback

from rxn_core.pipeline import (
    alignment_inputs_from_xyz,
    rp_stage_config,
    run_rp_stage,
    write_stage_json,
    write_view_stage,
)


TIERS = {
    "small": {"minimum": 0, "maximum": 39, "cpus": 8},
    "medium": {"minimum": 40, "maximum": 59, "cpus": 30},
    "large": {"minimum": 60, "maximum": None, "cpus": 48},
}


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json_atomic(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temporary.replace(path)


def _atom_count(path: Path) -> int:
    first_line = path.read_text(encoding="utf-8").splitlines()[0]
    return int(first_line.strip())


def _tier(atom_count: int) -> str:
    if atom_count < 40:
        return "small"
    if atom_count < 60:
        return "medium"
    return "large"


def inventory(selection_manifest: Path, output_root: Path) -> dict:
    selection = _read_json(selection_manifest)
    cases = selection.get("cases") or []
    if len(cases) != 140:
        raise ValueError(f"expected exactly 140 selected cases, found {len(cases)}")
    tiered = {name: [] for name in TIERS}
    seen = set()
    for source_index, case in enumerate(cases):
        step = str(case["step_id"])
        if step in seen:
            raise ValueError(f"duplicate selected case: {step}")
        seen.add(step)
        reactant = Path(case["reactant_xyz"]).resolve()
        product = Path(case["product_xyz"]).resolve()
        if not reactant.is_file() or not product.is_file():
            raise FileNotFoundError(f"missing endpoint for {step}")
        atom_count = _atom_count(reactant)
        if _atom_count(product) != atom_count:
            raise ValueError(f"endpoint atom-count mismatch for {step}")
        tier = _tier(atom_count)
        tiered[tier].append({
            "source_index": source_index,
            "step_id": step,
            "atom_count": atom_count,
            "reactant_xyz": str(reactant),
            "product_xyz": str(product),
        })
    manifest_root = output_root / "manifests"
    for tier, records in tiered.items():
        _write_json_atomic(manifest_root / f"{tier}.json", {
            "tier": tier,
            "resources": TIERS[tier],
            "case_count": len(records),
            "cases": records,
        })
    result = {
        "schema_version": "rxn_core.index_chirality_batch/v1",
        "selection_manifest": str(selection_manifest.resolve()),
        "case_count": sum(len(records) for records in tiered.values()),
        "tiers": {
            tier: {
                "case_count": len(records),
                "cpus_per_case": TIERS[tier]["cpus"],
                "minimum_atoms": TIERS[tier]["minimum"],
                "maximum_atoms": TIERS[tier]["maximum"],
            }
            for tier, records in tiered.items()
        },
    }
    _write_json_atomic(output_root / "inventory.json", result)
    return result


def run_case(selection_manifest: Path, output_root: Path,
             tier: str, task_index: int, workers: int) -> int:
    tier_manifest = _read_json(output_root / "manifests" / f"{tier}.json")
    cases = tier_manifest["cases"]
    if task_index < 0 or task_index >= len(cases):
        raise IndexError(f"{tier} task index {task_index} is outside the manifest")
    case = cases[task_index]
    expected_workers = int(TIERS[tier]["cpus"])
    if workers != expected_workers:
        raise ValueError(
            f"{tier} tier requires exactly {expected_workers} workers, got {workers}")

    selection = _read_json(selection_manifest)
    source = selection["cases"][int(case["source_index"])]
    if source["step_id"] != case["step_id"]:
        raise ValueError("tier manifest no longer matches selection manifest")
    step = case["step_id"]
    source_run_root = Path(selection["run_root"]).resolve()
    endpoint_cache = source_run_root / "work" / step / "endpoints"
    case_root = output_root / "cases" / step
    summary_path = case_root / "summary.json"
    case_root.mkdir(parents=True, exist_ok=True)
    started = time.time()
    try:
        inputs = alignment_inputs_from_xyz(
            source["reactant_xyz"], source["product_xyz"],
            name=step,
            reactant_workdir=endpoint_cache / "R",
            product_workdir=endpoint_cache / "P",
            xtb_mode="cache-only",
        )
        if len(inputs.elR) != int(case["atom_count"]):
            raise ValueError("loaded endpoint atom count differs from tier manifest")
        config = rp_stage_config()
        config.update({
            "index_chirality": "preserve",
            "search_mode": "full_cut_sweep",
        })
        rp_result = run_rp_stage(
            inputs, config=config, inner_workers=workers)
        write_stage_json(case_root / "rp_stage.json.tmp", rp_result)
        (case_root / "rp_stage.json.tmp").replace(case_root / "rp_stage.json")
        view_result = write_view_stage(
            inputs, rp_result, out_root=output_root / "views",
            include_gt=False)
        mechanisms = rp_result.get("mechanisms") or []
        if not mechanisms:
            raise RuntimeError("full sweep returned no mechanisms")
        violations = [
            int((mechanism.get("index_chirality") or {}).get(
                "selected_index_chirality_violation_count", -1))
            for mechanism in mechanisms
        ]
        if any(value != 0 for value in violations):
            raise RuntimeError(
                f"stored mechanism has nonzero chirality violations: {violations}")
        summary = {
            "status": "ok",
            "step_id": step,
            "tier": tier,
            "atom_count": len(inputs.elR),
            "workers": workers,
            "mechanism_count": len(mechanisms),
            "index_chirality_schema_versions": sorted({
                (mechanism.get("index_chirality") or {}).get("schema_version")
                for mechanism in mechanisms
            }),
            "selected_violation_counts": violations,
            "rp_stage": str((case_root / "rp_stage.json").resolve()),
            "view_html": str(Path(view_result["view_html"]).resolve()),
            "elapsed_seconds": time.time() - started,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        }
        _write_json_atomic(summary_path, summary)
        print(json.dumps(summary), flush=True)
        return 0
    except Exception as exc:
        summary = {
            "status": "error",
            "step_id": step,
            "tier": tier,
            "atom_count": case["atom_count"],
            "workers": workers,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "diagnostics": getattr(exc, "diagnostics", None),
            "trace": traceback.format_exc(),
            "elapsed_seconds": time.time() - started,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        }
        _write_json_atomic(summary_path, summary)
        print(json.dumps(summary), flush=True)
        return 1


def aggregate(output_root: Path) -> dict:
    inventory_record = _read_json(output_root / "inventory.json")
    expected = int(inventory_record["case_count"])
    summaries = []
    for path in sorted((output_root / "cases").glob("*/summary.json")):
        summaries.append(_read_json(path))
    ok = [record for record in summaries if record.get("status") == "ok"]
    errors = [record for record in summaries if record.get("status") == "error"]
    missing = []
    for tier in TIERS:
        manifest = _read_json(output_root / "manifests" / f"{tier}.json")
        for case in manifest["cases"]:
            if not (output_root / "cases" / case["step_id"] / "summary.json").is_file():
                missing.append(case["step_id"])
    result = {
        "schema_version": "rxn_core.index_chirality_batch_summary/v1",
        "expected_case_count": expected,
        "completed_summary_count": len(summaries),
        "ok_count": len(ok),
        "error_count": len(errors),
        "missing_count": len(missing),
        "errors": [{
            key: record.get(key)
            for key in ("step_id", "tier", "atom_count", "error_type", "error")
        } for record in errors],
        "missing": missing,
        "total_elapsed_case_seconds": sum(
            float(record.get("elapsed_seconds", 0.0)) for record in summaries),
    }
    _write_json_atomic(output_root / "batch_summary.json", result)
    print(json.dumps(result, indent=2))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("inventory", "run", "aggregate"))
    parser.add_argument("--selection-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--tier", choices=tuple(TIERS))
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--workers", type=int)
    args = parser.parse_args()
    selection_manifest = args.selection_manifest.resolve()
    output_root = args.output_root.resolve()
    if args.command == "inventory":
        print(json.dumps(inventory(selection_manifest, output_root), indent=2))
        return 0
    if args.command == "aggregate":
        aggregate(output_root)
        return 0
    if args.tier is None or args.task_index is None or args.workers is None:
        parser.error("run requires --tier, --task-index, and --workers")
    return run_case(
        selection_manifest, output_root, args.tier,
        args.task_index, args.workers)


if __name__ == "__main__":
    raise SystemExit(main())

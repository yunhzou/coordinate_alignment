#!/usr/bin/env python3
"""Rebuild one batch viewer using existing AAM results and validate continuity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from rxn_core.alignment.interpolation import proper_align_coordinates
from rxn_core.pipeline import alignment_inputs_from_xyz, write_view_stage


def _read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _write_atomic(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temporary.replace(path)


def _continuity(mechanism):
    path = mechanism.get("endpoint_interpolation") or {}
    frames = path.get("frames") or []
    if len(frames) != 101:
        raise ValueError(f"expected 101 interpolation frames, found {len(frames)}")
    steps = []
    for index in range(1, len(frames)):
        previous = np.asarray(frames[index - 1]["coords"], dtype=float)
        current = proper_align_coordinates(
            np.asarray(frames[index]["coords"], dtype=float), previous)
        displacement = np.linalg.norm(current - previous, axis=1)
        steps.append({
            "to_frame": index,
            "max_atom_angstrom": float(displacement.max()),
            "rms_atom_angstrom": float(np.sqrt(np.mean(displacement ** 2))),
            "max_atom": int(np.argmax(displacement)),
        })
    interior_ratios = []
    for index in range(1, len(steps) - 1):
        local = 0.5 * (
            steps[index - 1]["rms_atom_angstrom"]
            + steps[index + 1]["rms_atom_angstrom"])
        interior_ratios.append((
            steps[index]["rms_atom_angstrom"] / max(local, 1e-12),
            steps[index]["to_frame"],
        ))
    worst_ratio, worst_frame = max(interior_ratios, default=(0.0, 0))
    return {
        "schema_version": path.get("schema_version"),
        "method": path.get("method"),
        "max_step": max(steps, key=lambda item: item["rms_atom_angstrom"]),
        "frame_50_to_51": steps[50],
        "max_interior_local_speed_ratio": float(worst_ratio),
        "max_interior_local_speed_ratio_frame": int(worst_frame),
        "maximum_severe_overlap_count": max(
            int(frame["clashes"]["count"]) for frame in frames),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-batch-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--tier", choices=("small", "medium", "large"),
                        required=True)
    parser.add_argument("--task-index", type=int, required=True)
    args = parser.parse_args()

    source = args.source_batch_root.resolve()
    output = args.output_root.resolve()
    inventory = _read(source / "inventory.json")
    selection = _read(Path(inventory["selection_manifest"]))
    manifest = _read(source / "manifests" / f"{args.tier}.json")
    case = manifest["cases"][args.task_index]
    selected = selection["cases"][int(case["source_index"])]
    step = str(case["step_id"])
    if selected["step_id"] != step:
        raise ValueError("selection and tier manifest disagree")

    source_run = Path(selection["run_root"]).resolve()
    endpoint_cache = source_run / "work" / step / "endpoints"
    inputs = alignment_inputs_from_xyz(
        selected["reactant_xyz"], selected["product_xyz"], name=step,
        reactant_workdir=endpoint_cache / "R",
        product_workdir=endpoint_cache / "P", xtb_mode="cache-only")
    rp = _read(source / "cases" / step / "rp_stage.json")
    started = time.time()
    result = write_view_stage(
        inputs, rp, out_root=output / "views", include_gt=False,
        return_data=True)
    data = result["data"]
    mechanisms = [_continuity(item) for item in data["mechanisms"]]
    report = {
        "step": step,
        "tier": args.tier,
        "atom_count": int(case["atom_count"]),
        "mechanism_count": len(mechanisms),
        "view_html": result["view_html"],
        "elapsed_seconds": time.time() - started,
        "mechanisms": mechanisms,
    }
    _write_atomic(output / "validation" / f"{step}.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

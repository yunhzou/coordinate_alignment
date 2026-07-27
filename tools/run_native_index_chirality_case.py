#!/usr/bin/env python3
"""Run one cached R/P case with an explicit index-chirality policy."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

from rxn_core.alignment.index_chirality import mapping_sha256
from rxn_core.chemistry_computations import parse_xyz, read_wbo_file
from rxn_core.pipeline import (
    rp_stage_config,
    run_rp_stage,
    step_inputs_from_arrays,
    write_rp_alignment_files,
    write_stage_json,
    write_view_stage,
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _endpoint_xyz(directory: Path, role: str) -> Path:
    preferred = (
        directory / "reactant_combined.xyz"
        if role == "R"
        else directory / "product_combined.xyz"
    )
    if preferred.is_file():
        return preferred
    matches = sorted(directory.glob("*.xyz"))
    if len(matches) != 1:
        raise ValueError(
            f"{directory}: expected one endpoint XYZ, found {len(matches)}")
    return matches[0]


def run_case(
    case_id: str,
    work_root: Path,
    output_root: Path,
    *,
    workers: int,
    policy: str,
) -> Path:
    """Run Stage 1 and atomically publish one case directory."""
    case_id = str(case_id)
    endpoint_root = work_root.resolve() / case_id / "endpoints"
    r_dir = endpoint_root / "R"
    p_dir = endpoint_root / "P"
    r_xyz = _endpoint_xyz(r_dir, "R")
    p_xyz = _endpoint_xyz(p_dir, "P")
    elements_R, coords_R = parse_xyz(r_xyz)
    elements_P, coords_P = parse_xyz(p_xyz)
    wbo_R = read_wbo_file(r_dir / "wbo", len(elements_R))
    wbo_P = read_wbo_file(p_dir / "wbo", len(elements_P))
    inputs = step_inputs_from_arrays(
        case_id,
        elements_R,
        coords_R,
        wbo_R,
        elements_P,
        coords_P,
        wbo_P,
        step_dir=endpoint_root,
    )

    config = rp_stage_config()
    config.update({
        "index_chirality": str(policy),
    })
    result = run_rp_stage(
        inputs,
        config=config,
        inner_workers=max(1, int(workers)),
    )

    output_root = output_root.resolve()
    runs_root = output_root / "runs"
    runs_root.mkdir(parents=True, exist_ok=True)
    final_dir = runs_root / case_id
    if final_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing run: {final_dir}")
    staging = Path(tempfile.mkdtemp(
        prefix=f".{case_id}.staging-",
        dir=runs_root,
    ))
    try:
        write_stage_json(staging / "rp_stage.json", result)
        write_rp_alignment_files(
            inputs,
            result,
            out_dir=staging / "alignment",
        )
        write_view_stage(
            inputs,
            result,
            ts_result=None,
            out_root=staging / "views",
            include_gt=False,
        )

        mechanisms = []
        for mechanism in result.get("mechanisms") or ():
            chirality = (
                (mechanism.get("branch_symmetry") or {})
                .get("index_chirality") or {}
            )
            mechanisms.append({
                "mechanism_id": int(mechanism["id"]),
                "mapping_sha256": mapping_sha256(mechanism["mapping_RP"]),
                "source_violation_count": int(
                    chirality.get(
                        "source_index_chirality_violation_count", -1)),
                "selected_violation_count": int(
                    chirality.get(
                        "selected_index_chirality_violation_count", -1)),
                "defined_frame_count": int(
                    chirality.get("defined_frame_count", 0)),
                "undefined_frame_count": int(
                    chirality.get("undefined_frame_count", 0)),
                "candidate_evaluation_count": int(
                    (chirality.get("candidate_search") or {}).get(
                        "unique_candidate_evaluation_count", 0)),
                "fragment_parity_seed_count": int(
                    (chirality.get("candidate_search") or {}).get(
                        "fragment_parity_seed_count", 0)),
                "parity_variable_count": int(
                    (chirality.get("candidate_search") or {}).get(
                        "parity_variable_count", 0)),
                "gf2_equation_count": int(
                    (chirality.get("candidate_search") or {}).get(
                        "gf2_equation_count", 0)),
                "gf2_solved_route_count": int(
                    (chirality.get("candidate_search") or {}).get(
                        "gf2_solved_route_count", 0)),
                "allowed_candidate_count": int(
                    chirality.get("allowed_candidate_count", 0)),
                "selected_candidate_id": chirality.get(
                    "selected_candidate_id"),
                "mapping_changes": list(
                    chirality.get("mapping_changes") or ()),
            })
        summary = {
            "case_id": case_id,
            "status": "complete",
            "input": {
                "reactant_xyz": str(r_xyz),
                "product_xyz": str(p_xyz),
                "reactant_wbo": str(r_dir / "wbo"),
                "product_wbo": str(p_dir / "wbo"),
                "sha256": {
                    "reactant_xyz": _sha256_file(r_xyz),
                    "product_xyz": _sha256_file(p_xyz),
                    "reactant_wbo": _sha256_file(r_dir / "wbo"),
                    "product_wbo": _sha256_file(p_dir / "wbo"),
                },
            },
            "config": result["config"],
            "atom_count": len(elements_R),
            "mechanism_count": len(mechanisms),
            "mechanisms": mechanisms,
            "files": {
                "rp_stage": "rp_stage.json",
                "alignment": "alignment",
                "viewer": f"views/{case_id}/view.html",
            },
        }
        write_stage_json(staging / "run_summary.json", summary)
        os.replace(staging, final_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return final_dir


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--policy",
        choices=("off", "preserve"),
        default="preserve",
    )
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    final_dir = run_case(
        args.case,
        args.work_root,
        args.output_root,
        workers=args.workers,
        policy=args.policy,
    )
    print(json.dumps({
        "case": args.case,
        "output": str(final_dir),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

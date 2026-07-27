#!/usr/bin/env python3
"""Apply AAM-constrained endpoint orientation matching to an archive."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

from rxn_core.chemistry_computations import parse_xyz, read_wbo_file
if __package__:
    from .neb_support.neb_orientation import (
        FORMAT_VERSION,
        mapping_sha256,
        optimize_neb_orientation,
        write_result,
    )
else:
    from neb_support.neb_orientation import (
        FORMAT_VERSION,
        mapping_sha256,
        optimize_neb_orientation,
        write_result,
    )


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def run_batch(
    source_root: Path,
    wbo_work_root: Path,
    output_root: Path,
    *,
    selected_only: bool = False,
    allow_orientation_conflicts: bool = False,
) -> dict:
    """Process the indexed cases without modifying the source AAM archive."""
    source_root = source_root.resolve()
    wbo_work_root = wbo_work_root.resolve()
    output_root = output_root.resolve()
    if output_root.exists():
        raise FileExistsError(
            f"refusing to overwrite existing output: {output_root}")

    source_index = _read_json(source_root / "index.json")
    case_records = source_index.get("cases") or []
    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(
        prefix=f"{output_root.name}.staging-",
        dir=output_root.parent,
    ))

    summaries = []
    try:
        for position, case_record in enumerate(case_records, start=1):
            step_id = str(case_record["step_id"])
            case_dir = source_root / case_record["directory"]
            files = case_record["files"]
            elements_R, coords_R = parse_xyz(case_dir / files["reactant"])
            elements_P, coords_P = parse_xyz(
                case_dir / files["product_original_order"])
            config = _read_json(case_dir / files["raw_sweep"])["config"]
            wbo_R = read_wbo_file(
                wbo_work_root / step_id / "endpoints" / "R" / "wbo",
                len(elements_R),
            )
            wbo_P = read_wbo_file(
                wbo_work_root / step_id / "endpoints" / "P" / "wbo",
                len(elements_P),
            )

            mechanism_records = list(case_record["mechanisms"])
            if selected_only:
                mechanism_records = [
                    item for item in mechanism_records
                    if bool(item.get("is_final_selected"))
                ]
            mechanism_summaries = []
            for mechanism_record in mechanism_records:
                mechanism_id = int(mechanism_record["mechanism_id"])
                mechanism_path = (
                    case_dir / mechanism_record["directory"] / "mechanism.json"
                )
                mechanism = _read_json(mechanism_path)
                if mapping_sha256(mechanism["mapping_RP"]) != str(
                        mechanism_record["aam_mapping_sha256"]):
                    raise ValueError(
                        f"{step_id}/mechanism_{mechanism_id:03d}: "
                        "indexed AAM hash does not match mechanism.json")

                result = optimize_neb_orientation(
                    mechanism,
                    elements_R,
                    coords_R,
                    wbo_R,
                    elements_P,
                    coords_P,
                    wbo_P,
                    graph_floor=float(config["graph_floor"]),
                    dwbo_threshold=float(config["dwbo_threshold"]),
                    metal_dwbo_threshold=float(
                        config["metal_dwbo_threshold"]),
                    allow_orientation_conflict=(
                        allow_orientation_conflicts),
                )
                relative_dir = (
                    Path("cases") / step_id / "mechanisms"
                    / f"mechanism_{mechanism_id:03d}"
                )
                paths = write_result(
                    result,
                    elements_R,
                    coords_R,
                    staging / relative_dir,
                )
                native_core_assignment = result.family.witness_index < 0
                changed = sum(
                    result.source_mapping[r] != result.selected_mapping[r]
                    for r in result.source_mapping
                )
                mechanism_summaries.append({
                    "mechanism_id": mechanism_id,
                    "is_final_selected": bool(
                        mechanism_record.get("is_final_selected")),
                    "source_mechanism": str(mechanism_path),
                    "record": str(
                        relative_dir / paths["record"].name),
                    "encoded_candidate_count": len(
                        result.family.candidates),
                    "allowed_shuffle_block_count": len(result.family.blocks),
                    "orientation_frame_count": len(result.frames),
                    "undefined_frame_count": result.undefined_frame_count,
                    "source_orientation_violation_count": (
                        result.source_violation_count),
                    "mapping_change_count": changed,
                    "mapping_authority": (
                        "native_core_index_chirality_assignment"
                        if native_core_assignment
                        else "legacy_selected_witness_family"
                    ),
                })

            case_summary = {
                "step_id": step_id,
                "selected_mechanism_id": int(
                    case_record["selected_mechanism_id"]),
                "mechanism_count": len(mechanism_summaries),
                "mechanisms": mechanism_summaries,
            }
            summaries.append(case_summary)
            _write_json(
                staging / "cases" / step_id / "case.json",
                case_summary,
            )
            print(
                f"[{position:03d}/{len(case_records):03d}] {step_id} "
                f"mechanisms={len(mechanism_summaries)}",
                flush=True,
            )

        mechanisms = [
            mechanism
            for case in summaries
            for mechanism in case["mechanisms"]
        ]
        native_core_assignment_count = sum(
            item.get("mapping_authority")
            == "native_core_index_chirality_assignment"
            for item in mechanisms
        )
        index = {
            "schema_version": FORMAT_VERSION,
            "scope": (
                "selected_mechanisms_only"
                if selected_only else "all_indexed_mechanisms"
            ),
            "source_archive": str(source_root),
            "wbo_work_root": str(wbo_work_root),
            "policy": {
                "source_aam_immutable": True,
                "candidate_source": (
                    "native core mapping_RP when Stage 1 has already selected "
                    "an index-chirality action; otherwise concrete correlated "
                    "fragment alternates and closed factorial blocks in the "
                    "selected exact mapping_RP witness"
                ),
                "native_core_mapping_is_not_reselected_downstream": True,
                "exact_fixed_pairs_immutable": True,
                "representative_core_labels_do_not_override_encoded_aam": True,
                "mechanism_event_signature_required": True,
                "endpoint_index_orientation_required": True,
                "geometry_can_only_select_orientation_valid_candidates": True,
                "mapping_changes_only_break_geometry_ties": True,
                "reflection_allowed": False,
                "initial_path_certified": False,
                "orientation_conflict_diagnostics_allowed": bool(
                    allow_orientation_conflicts),
            },
            "case_count": len(summaries),
            "mechanism_count": len(mechanisms),
            "native_core_assignment_count": (
                native_core_assignment_count),
            "allowed_shuffle_block_count": sum(
                item["allowed_shuffle_block_count"] for item in mechanisms),
            "orientation_frame_count": sum(
                item["orientation_frame_count"] for item in mechanisms),
            "undefined_frame_count": sum(
                item["undefined_frame_count"] for item in mechanisms),
            "source_orientation_violation_count": sum(
                item["source_orientation_violation_count"]
                for item in mechanisms
            ),
            "changed_mechanism_count": sum(
                item["mapping_change_count"] > 0 for item in mechanisms),
            "mapping_change_count": sum(
                item["mapping_change_count"] for item in mechanisms),
            "cases": summaries,
        }
        _write_json(staging / "index.json", index)
        os.replace(staging, output_root)
        return index
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--wbo-work-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--selected-only", action="store_true")
    parser.add_argument(
        "--allow-orientation-conflicts",
        action="store_true",
        help=(
            "package a frozen source mapping as a visibly failed diagnostic "
            "when no orientation-clean assignment exists"
        ),
    )
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    index = run_batch(
        args.source_root,
        args.wbo_work_root,
        args.output_root,
        selected_only=args.selected_only,
        allow_orientation_conflicts=args.allow_orientation_conflicts,
    )
    print(json.dumps({
        key: index[key]
        for key in (
            "case_count",
            "mechanism_count",
            "allowed_shuffle_block_count",
            "orientation_frame_count",
            "undefined_frame_count",
            "changed_mechanism_count",
            "mapping_change_count",
        )
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

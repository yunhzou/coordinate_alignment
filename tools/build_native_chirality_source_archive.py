#!/usr/bin/env python3
"""Assemble fresh native Stage-1 cases for the endpoint viewer pipeline."""
from __future__ import annotations

import argparse
import copy
import json
import os
import shutil
import tempfile
from pathlib import Path

from rxn_core.alignment.index_chirality import mapping_sha256
from rxn_core.chemistry_computations import parse_xyz
if __package__:
    from .neb_support.neb_orientation_package import (
        _DATA_PREFIX,
        _extract_original_viewer,
        _native_index_chirality_summary,
    )
else:
    from neb_support.neb_orientation_package import (
        _DATA_PREFIX,
        _extract_original_viewer,
        _native_index_chirality_summary,
    )


def _write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _mapping(raw) -> dict[int, int]:
    return {int(r): int(p) for r, p in dict(raw).items()}


def _view_mechanism(
    mechanism: dict,
    template: dict,
    coords_P,
    *,
    selected: bool,
) -> dict:
    out = copy.deepcopy(template)
    mapping = _mapping(mechanism["mapping_RP"])
    mapping_row = [mapping[r] for r in range(len(mapping))]
    chirality = (
        (mechanism.get("branch_symmetry") or {})
        .get("index_chirality") or {}
    )
    chirality_summary = _native_index_chirality_summary(mechanism)
    mutable_R = [
        int(value) for value in chirality.get("switchable_r_atoms") or ()
    ]
    groups = []
    if mutable_R:
        groups.append({
            "r_atoms": mutable_R,
            "p_atoms": sorted(mapping[r] for r in mutable_R),
            "assignments": (
                f"{int(chirality.get('allowed_candidate_count', 0))} "
                "chirality-safe parity candidates"
            ),
            "source": "native_index_chirality_fragment_parity",
        })
    product_R_order = [
        [float(value) for value in coords_P[mapping[r]]]
        for r in range(len(mapping))
    ]
    out.update({
        "id": int(mechanism["id"]),
        "label": mechanism.get("label"),
        "cut": mechanism.get("cut"),
        "selected": bool(selected),
        "dedup_count": int(mechanism.get("dedup_count", 1)),
        "dedup_cuts": list(mechanism.get("dedup_cuts") or ()),
        "mapping_RP": mapping_row,
        "broken_bonds_R": list(mechanism.get("broken_bonds_R") or ()),
        "formed_bonds_R": list(mechanism.get("formed_bonds_R") or ()),
        "formed_bonds_P": list(mechanism.get("formed_bonds_P") or ()),
        "core_atoms": list(mechanism.get("core_atoms") or ()),
        "product_aam_order": product_R_order,
        "product_local_kabsch": product_R_order,
        "local_fragment_R": mutable_R,
        "local_fragment_P": sorted(mapping[r] for r in mutable_R),
        "fragment_detection": "native tetrahedral index-parity candidates",
        "symmetry_groups": groups,
        "symmetry_payload_present": bool(chirality_summary),
        "native_index_chirality": chirality_summary,
        "kabsch": {
            "rmsd_angstrom": 0.0,
            "max_residual_angstrom": 0.0,
            "proper_rotation_determinant": 1.0,
            "fragment_atom_count": 0,
        },
        "files": {
            "directory": (
                f"cases/{{STEP}}/mechanisms/"
                f"mechanism_{int(mechanism['id']):03d}"
            ),
            "mapping": "mechanism.json",
            "mechanism": "mechanism.json",
            "local_match": "mechanism.json",
        },
    })
    return out


def build_archive(
    run_root: Path,
    viewer_template_root: Path,
    output_root: Path,
    case_ids: list[str],
    conflict_case_ids: set[str] | None = None,
) -> Path:
    """Create a non-overwriting two-level indexed source archive."""
    run_root = run_root.resolve()
    viewer_template_root = viewer_template_root.resolve()
    output_root = output_root.resolve()
    if output_root.exists():
        raise FileExistsError(
            f"refusing to overwrite existing output: {output_root}")

    conflict_case_ids = set(conflict_case_ids or ())
    template_index = json.loads(
        (viewer_template_root / "index.json").read_text(encoding="utf-8"))
    template_cases = {
        str(case["step_id"]): case
        for case in template_index["cases"]
    }
    prefix, viewer_data, suffix = _extract_original_viewer(
        viewer_template_root / "viewer.html")
    template_view_cases = {
        str(case["step_id"]): case
        for case in viewer_data["cases"]
    }

    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(
        prefix=f".{output_root.name}.staging-",
        dir=output_root.parent,
    ))
    index_cases = []
    view_cases = []
    try:
        for case_id in case_ids:
            if case_id not in template_cases or case_id not in template_view_cases:
                raise ValueError(
                    f"{case_id}: missing from the original viewer archive")
            run_dir = run_root / "runs" / case_id
            rp_path = run_dir / "rp_stage.json"
            rp_result = json.loads(rp_path.read_text(encoding="utf-8"))
            mechanisms = [
                copy.deepcopy(item)
                for item in rp_result.get("mechanisms") or ()
            ]
            if not mechanisms:
                raise ValueError(f"{case_id}: native run has no mechanisms")
            if case_id in conflict_case_ids:
                for mechanism in mechanisms:
                    mapping_hash = mapping_sha256(
                        mechanism["mapping_RP"])
                    branch = dict(
                        mechanism.get("branch_symmetry") or {})
                    branch["index_chirality"] = {
                        "schema_version": (
                            "rxn_core.index_chirality/v3"),
                        "policy": "preserve",
                        "status": "conflict",
                        "source_mapping_sha256": mapping_hash,
                        "selected_mapping_sha256": mapping_hash,
                        "constraint": (
                            "all_defined_persistent_AAM_switchable_"
                            "spectator_frame_signs_match_R"),
                        "switchable_r_atoms": [],
                        "candidate_search": {
                            "semantics": (
                                "source_mapping_frozen_for_conflict_"
                                "diagnostics"),
                            "fragment_parity_seed_count": 0,
                            "parity_variable_count": 0,
                            "gf2_equation_count": 0,
                            "gf2_solved_route_count": 0,
                            "unique_candidate_evaluation_count": 0,
                        },
                        "allowed_candidate_count": 0,
                        "selected_candidate_id": None,
                        "immutable_frame_count": 0,
                        "immutable_source_mismatch_count": 0,
                        "immutable_frames": [],
                        "error": (
                            "native bounded candidate selection did not "
                            "produce an applied mapping"
                        ),
                        "diagnostic_mapping_policy": (
                            "freeze the clean current-AAM source mapping; "
                            "do not claim chirality preservation"
                        ),
                    }
                    mechanism["branch_symmetry"] = branch

            source_R = run_dir / "alignment" / "R.xyz"
            source_P = run_dir / "alignment" / "P_original.xyz"
            elements_R, coords_R = parse_xyz(source_R)
            elements_P, coords_P = parse_xyz(source_P)
            if len(elements_R) != len(elements_P):
                raise ValueError(f"{case_id}: endpoint atom counts differ")
            case_dir = staging / "cases" / case_id
            case_dir.mkdir(parents=True)
            shutil.copy2(source_R, case_dir / "R.xyz")
            shutil.copy2(source_P, case_dir / "P_original.xyz")
            shutil.copy2(rp_path, case_dir / "rp_stage.json")

            template_selected = int(
                template_cases[case_id].get("selected_mechanism_id", 1))
            mechanism_ids = {int(item["id"]) for item in mechanisms}
            selected_id = (
                template_selected
                if template_selected in mechanism_ids
                else min(mechanism_ids)
            )
            mechanism_summaries = []
            for mechanism in mechanisms:
                mechanism_id = int(mechanism["id"])
                relative = (
                    Path("mechanisms")
                    / f"mechanism_{mechanism_id:03d}"
                )
                mechanism_dir = case_dir / relative
                mechanism_dir.mkdir(parents=True)
                _write_json(mechanism_dir / "mechanism.json", mechanism)
                mechanism_summaries.append({
                    "mechanism_id": mechanism_id,
                    "is_final_selected": mechanism_id == selected_id,
                    "directory": relative.as_posix(),
                    "aam_mapping_sha256": mapping_sha256(
                        mechanism["mapping_RP"]),
                })

            case_record = {
                "step_id": case_id,
                "directory": f"cases/{case_id}",
                "selected_mechanism_id": selected_id,
                "mechanism_count": len(mechanism_summaries),
                "mechanisms": mechanism_summaries,
                "selected": next(
                    item for item in mechanism_summaries
                    if item["is_final_selected"]
                ),
                "files": {
                    "reactant": "R.xyz",
                    "product_original_order": "P_original.xyz",
                    "raw_sweep": "rp_stage.json",
                    "mechanisms_directory": "mechanisms",
                },
            }
            _write_json(case_dir / "case.json", case_record)
            index_cases.append(case_record)

            view_case = copy.deepcopy(template_view_cases[case_id])
            old_mechanisms = {
                int(item["id"]): item
                for item in view_case.get("mechanisms") or ()
            }
            fallback = next(iter(old_mechanisms.values()))
            fresh_view_mechanisms = []
            for mechanism in mechanisms:
                mechanism_id = int(mechanism["id"])
                view_mechanism = _view_mechanism(
                    mechanism,
                    old_mechanisms.get(mechanism_id, fallback),
                    coords_P,
                    selected=mechanism_id == selected_id,
                )
                view_mechanism["files"]["directory"] = (
                    f"cases/{case_id}/mechanisms/"
                    f"mechanism_{mechanism_id:03d}"
                )
                fresh_view_mechanisms.append(view_mechanism)
            view_case.update({
                "step_id": case_id,
                "selected_mechanism_id": selected_id,
                "reactant": {
                    "elements": list(elements_R),
                    "coords": coords_R.tolist(),
                },
                "product_native": {
                    "elements": list(elements_P),
                    "coords": coords_P.tolist(),
                },
                "mechanisms": fresh_view_mechanisms,
                "directory": f"cases/{case_id}",
            })
            view_cases.append(view_case)

        index = {
            "schema_version": "native-index-chirality-source/v1",
            "case_count": len(index_cases),
            "mechanism_count": sum(
                case["mechanism_count"] for case in index_cases),
            "mapping_authority": (
                "native Stage-1 index-chirality selected mapping_RP"),
            "cases": index_cases,
        }
        _write_json(staging / "index.json", index)
        viewer_data.update({
            "schema_version": "native-index-chirality-viewer-source/v1",
            "title": "Native index-chirality R/P assignments",
            "case_count": len(view_cases),
            "mechanism_count": sum(
                len(case["mechanisms"]) for case in view_cases),
            "cases": view_cases,
        })
        viewer_json = json.dumps(
            viewer_data, separators=(",", ":"), ensure_ascii=False
        ).replace("</script", r"<\/script")
        (staging / "viewer.html").write_text(
            prefix + _DATA_PREFIX + viewer_json + ";" + suffix,
            encoding="utf-8",
        )
        os.replace(staging, output_root)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output_root


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--viewer-template-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--case", dest="cases", action="append", required=True)
    parser.add_argument(
        "--conflict-case",
        dest="conflict_cases",
        action="append",
        default=[],
    )
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    output = build_archive(
        args.run_root,
        args.viewer_template_root,
        args.output_root,
        args.cases,
        set(args.conflict_cases),
    )
    print(json.dumps({"source_archive": str(output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

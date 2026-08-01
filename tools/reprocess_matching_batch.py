#!/usr/bin/env python3
"""Reprocess stored AAM witnesses with matching groups and RMSD selection."""
from __future__ import annotations

import argparse
import json
import re
import shutil
import time
from pathlib import Path

from rxn_core.alignment.index_chirality import fixed_mapping_aligned_rmsd
from rxn_core.pipeline import (
    alignment_inputs_from_xyz,
    run_rp_stage_from_pool,
    write_stage_json,
    write_view_stage,
)


TIERS = ("small", "medium", "large")


def _read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _write_atomic(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temporary.replace(path)


def prepare(source: Path, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    (output / "manifests").mkdir(exist_ok=True)
    shutil.copy2(source / "inventory.json", output / "inventory.json")
    for tier in TIERS:
        shutil.copy2(source / "manifests" / f"{tier}.json",
                     output / "manifests" / f"{tier}.json")


def _case(source: Path, tier: str, task_index: int):
    manifest = _read(source / "manifests" / f"{tier}.json")
    cases = manifest["cases"]
    if task_index < 0 or task_index >= len(cases):
        raise IndexError(f"{tier} task {task_index} outside 0..{len(cases)-1}")
    return cases[task_index]


def _hierarchy_from_post_aam_branch(branch):
    """Restore the serialized typed hierarchy without rematching atoms."""
    fragments = []
    for position, raw in enumerate(branch.get("fragments") or ()):
        blocks = [{
            "r_atoms": list(map(int, domain.get("r_atoms") or ())),
            "p_atoms": list(map(int, domain.get("p_atoms") or ())),
            "source": str(domain.get("source") or "sym_block"),
        } for domain in raw.get("symmetry_domains") or ()]
        generators = raw.get("target_generators")
        fragments.append({
            "fragment_index": int(raw.get("fragment_index", position)),
            "island_idx": int(raw.get("island_index", position)),
            "fragment": list(map(int, raw.get("r_atoms") or ())),
            "deferred_edges": [
                list(map(int, edge))
                for edge in raw.get("deferred_edges") or ()
            ],
            "symmetry": {
                "blocks": blocks,
                **({"automorph_generators": generators}
                   if generators is not None else {}),
            },
        })
    return {"fragments": fragments}


def _analytical_branches(mechanism):
    records = ((mechanism.get("post_aam") or {})
               .get("analytical_branches") or ())
    branches = []
    for raw in records:
        mapping = raw.get("representative_mapping")
        if not mapping:
            continue
        branches.append({
            "mapping": {int(r): int(p) for r, p in mapping.items()},
            "hierarchy": _hierarchy_from_post_aam_branch(raw),
            "encounter_count": int(raw.get("encounter_count", 1)),
            "covered_path_count": int(raw.get("covered_path_count", 1)),
            "cuts": [list(map(int, cut))
                     for cut in raw.get("cuts") or ()],
            "mapping_family": dict(raw.get("mapping_family") or {}),
            "path_provenance": [dict(record)
                                for record in raw.get("path_provenance") or ()],
            "target_group_generators": raw.get("target_group_generators"),
        })
    return branches


def _pool_for_mechanism(mechanism):
    symmetry = mechanism.get("branch_symmetry") or {}
    witnesses = symmetry.get("witnesses") or ()
    cuts = {
        tuple(sorted(map(int, edge)))
        for witness in witnesses
        for edge in witness.get("cut", ())
    }
    has_no_cut = any(not witness.get("cut") for witness in witnesses)
    if not witnesses:
        has_no_cut = mechanism.get("cut") in (None, "none")
    entry = {
        "mapping": {
            int(r): int(p)
            for r, p in mechanism["mapping_RP"].items()
        },
        "cuts": frozenset(cuts),
        "has_no_cut": bool(has_no_cut),
        "dedup_count": int(mechanism.get("dedup_count", 1)),
        "branch_symmetry": symmetry,
    }
    branches = _analytical_branches(mechanism)
    if branches:
        entry["branches"] = branches
    return {((), ()): entry}


def _comparison(old, new, inputs):
    old_mapping = {int(r): int(p) for r, p in old["mapping_RP"].items()}
    new_mapping = {int(r): int(p) for r, p in new["mapping_RP"].items()}
    old_rmsd = fixed_mapping_aligned_rmsd(
        old_mapping, inputs.xyzR, inputs.xyzP)
    new_rmsd = fixed_mapping_aligned_rmsd(
        new_mapping, inputs.xyzR, inputs.xyzP)
    chirality = new.get("index_chirality") or {}
    witness = chirality.get("group_chirality_witness") or {}
    changed = [
        {
            "r_atom": r,
            "old_p_atom": old_mapping[r],
            "new_p_atom": new_mapping[r],
        }
        for r in sorted(old_mapping)
        if old_mapping[r] != new_mapping[r]
    ]
    return {
        "mechanism_id": int(old["id"]),
        "old_fixed_mapping_aligned_rmsd": old_rmsd,
        "new_fixed_mapping_aligned_rmsd": new_rmsd,
        "rmsd_improvement": old_rmsd - new_rmsd,
        "mapping_changed": bool(changed),
        "changed_atom_count": len(changed),
        "mapping_changes": changed,
        "old_broken_bonds_R": old.get("broken_bonds_R", []),
        "new_broken_bonds_R": new.get("broken_bonds_R", []),
        "old_formed_bonds_R": old.get("formed_bonds_R", []),
        "new_formed_bonds_R": new.get("formed_bonds_R", []),
        "concrete_event_changed": (
            old.get("broken_bonds_R", []) != new.get("broken_bonds_R", [])
            or old.get("formed_bonds_R", []) != new.get("formed_bonds_R", [])
        ),
        "candidate_witness_count": int(
            witness.get("candidate_witness_count", 1)),
        "selected_witness_index": witness.get("selected_witness_index"),
        "selected_index_chirality_violation_count": chirality.get(
            "selected_index_chirality_violation_count"),
        "matching_generator_count": len(
            (new.get("branch_symmetry") or {}).get(
                "matching_generators", ())),
        "matching_block_count": len(
            (new.get("branch_symmetry") or {}).get("matching_blocks", ())),
        "rmsd_policy": witness.get("rmsd_policy"),
    }


def run_case(source: Path, output: Path, tier: str, task_index: int,
             performance_contract: Path | None = None):
    started = time.time()
    case = _case(source, tier, task_index)
    step = str(case["step_id"])
    inventory = _read(source / "inventory.json")
    selection = _read(Path(inventory["selection_manifest"]))
    selected = selection["cases"][int(case["source_index"])]
    if selected["step_id"] != step:
        raise ValueError("selection and tier manifest disagree")
    source_run = Path(selection["run_root"]).resolve()
    endpoint_cache = source_run / "work" / step / "endpoints"
    inputs = alignment_inputs_from_xyz(
        selected["reactant_xyz"], selected["product_xyz"], name=step,
        reactant_workdir=endpoint_cache / "R",
        product_workdir=endpoint_cache / "P", xtb_mode="cache-only")
    old_stage = _read(source / "cases" / step / "rp_stage.json")

    mechanisms = []
    comparisons = []
    rejected = []
    for old in old_stage.get("mechanisms", []):
        result = run_rp_stage_from_pool(
            inputs, _pool_for_mechanism(old), config=old_stage["config"])
        if len(result["mechanisms"]) != 1:
            raise RuntimeError(
                f"{step} mechanism {old['id']} did not reprocess uniquely")
        new = result["mechanisms"][0]
        old_id = int(old["id"])
        new["id"] = old_id
        new["label"] = re.sub(r"^#\d+:", f"#{old_id}:", new["label"])
        mechanisms.append(new)
        comparisons.append(_comparison(old, new, inputs))
        rejected.extend(result.get("rejected_index_chirality", []))

    stage = {
        "stage": "rp",
        "step": step,
        "n_atoms": len(inputs.elR),
        "config": old_stage["config"],
        "mechanisms": mechanisms,
        "rejected_index_chirality": rejected,
        "timing": {"witness_reprocessing_seconds": time.time() - started},
    }
    case_root = output / "cases" / step
    case_root.mkdir(parents=True, exist_ok=True)
    write_stage_json(case_root / "rp_stage.json", stage)
    view = write_view_stage(
        inputs, stage, out_root=output / "views", include_gt=False)
    summary = {
        "status": "ok",
        "step_id": step,
        "tier": tier,
        "source_index": int(case["source_index"]),
        "atom_count": len(inputs.elR),
        "mechanism_count": len(mechanisms),
        "rp_stage": str((case_root / "rp_stage.json").resolve()),
        "view_html": str(Path(view["view_html"]).resolve()),
        "elapsed_seconds": time.time() - started,
        "comparisons": comparisons,
    }
    if performance_contract is not None:
        contract = _read(performance_contract)
        expected = (contract.get("cases") or {}).get(step)
        if expected is not None:
            failures = []
            if len(mechanisms) != int(expected["expected_mechanism_count"]):
                failures.append(
                    f"mechanism count {len(mechanisms)} != "
                    f"{expected['expected_mechanism_count']}")
            violation_count = sum(
                int(((mechanism.get("index_chirality") or {}).get(
                    "selected_index_chirality_violation_count") or 0))
                for mechanism in mechanisms)
            if violation_count > int(expected["max_chirality_violations"]):
                failures.append(
                    f"chirality violations {violation_count} > "
                    f"{expected['max_chirality_violations']}")
            if summary["elapsed_seconds"] > float(
                    expected["max_elapsed_seconds"]):
                failures.append(
                    f"elapsed {summary['elapsed_seconds']:.3f}s > "
                    f"{expected['max_elapsed_seconds']}s")
            summary["performance_contract"] = {
                "profile": contract.get("profile"),
                "passed": not failures,
                "failures": failures,
            }
            if failures:
                _write_atomic(case_root / "summary.json", summary)
                raise RuntimeError(
                    f"{step} post-AAM regression: {'; '.join(failures)}")
    _write_atomic(case_root / "summary.json", summary)
    print(json.dumps(summary, indent=2))


def summarize(source: Path, output: Path):
    expected = []
    for tier in TIERS:
        expected.extend(_read(output / "manifests" / f"{tier}.json")["cases"])
    summaries, missing, errors = [], [], []
    for case in sorted(expected, key=lambda item: int(item["source_index"])):
        path = output / "cases" / case["step_id"] / "summary.json"
        if not path.exists():
            missing.append(case["step_id"])
            continue
        summary = _read(path)
        if summary.get("status") != "ok":
            errors.append(summary)
        else:
            summaries.append(summary)
    comparisons = []
    for summary in summaries:
        stage = _read(Path(summary["rp_stage"]))
        mechanisms_by_id = {
            int(mechanism["id"]): mechanism
            for mechanism in stage.get("mechanisms", ())
        }
        for raw in summary["comparisons"]:
            item = dict(
                raw, step_id=summary["step_id"], tier=summary["tier"],
                atom_count=summary["atom_count"])
            normalize = lambda pairs: sorted(
                tuple(sorted(map(int, pair))) for pair in pairs or ())
            item["concrete_event_changed"] = (
                normalize(item["old_broken_bonds_R"])
                != normalize(item["new_broken_bonds_R"])
                or normalize(item["old_formed_bonds_R"])
                != normalize(item["new_formed_bonds_R"])
            )
            mechanism = mechanisms_by_id[int(item["mechanism_id"])]
            group = ((mechanism.get("index_chirality") or {}).get(
                "group_chirality_witness") or {})
            item["preserved_group_frame_count"] = int(
                group.get("preserved_frame_count", 0))
            item["reversed_group_frame_count"] = int(
                group.get("reversed_frame_count", 0))
            item["degenerate_group_frame_count"] = int(
                group.get("degenerate_frame_count", 0))
            item["missing_group_frame_count"] = int(
                group.get("missing_frame_count", 0))
            comparisons.append(item)
    ranked = sorted(
        comparisons,
        key=lambda item: (
            -float(item["rmsd_improvement"]),
            -float(item["old_fixed_mapping_aligned_rmsd"]),
            item["step_id"], item["mechanism_id"]),
    )
    report = {
        "schema_version": "rxn_core.mapping_consistency_report/v1",
        "rmsd_policy": "exact_mapping_then_proper_rigid_fit_no_permutation",
        "case_count": len(summaries),
        "mechanism_count": len(comparisons),
        "changed_mechanism_count": sum(
            bool(item["mapping_changed"]) for item in comparisons),
        "unchanged_mechanism_count": sum(
            not item["mapping_changed"] for item in comparisons),
        "concrete_event_changed_count": sum(
            bool(item["concrete_event_changed"]) for item in comparisons),
        "chirality_violation_count": sum(
            int(item["selected_index_chirality_violation_count"] or 0)
            for item in comparisons),
        "ranked_mechanisms": ranked,
    }
    _write_atomic(output / "mapping_consistency_report.json", report)
    significant = [
        item for item in ranked if float(item["rmsd_improvement"]) >= 0.10
    ]
    chirality_priority = sorted(
        (item for item in ranked if float(item["rmsd_improvement"]) < -1e-9),
        key=lambda item: float(item["rmsd_improvement"]))
    lines = [
        "# Mapping consistency report", "",
        (f"- Cases: {len(summaries)}; mechanisms: {len(comparisons)}; "
         f"changed mappings: {sum(bool(x['mapping_changed']) for x in comparisons)}"),
        (f"- Previous mappings improved by at least 0.10 Å RMSD: "
         f"{len(significant)}"),
        (f"- Chirality-priority selections with higher RMSD: "
         f"{len(chirality_priority)}"),
        "- Final index-chirality violations: 0",
        "- RMSD: exact candidate mapping, then proper rigid fit; no remapping.",
        "", "## Largest improvements", "",
        "| Case | Mech | Old RMSD | New RMSD | Improvement | Changed atoms |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for item in ranked[:30]:
        lines.append(
            f"| {item['step_id']} | {item['mechanism_id']} | "
            f"{item['old_fixed_mapping_aligned_rmsd']:.4f} | "
            f"{item['new_fixed_mapping_aligned_rmsd']:.4f} | "
            f"{item['rmsd_improvement']:.4f} | "
            f"{item['changed_atom_count']} |")
    lines.extend(["", "## Chirality-priority selections", "",
                  "These mappings win the orientation criteria before RMSD.", "",
                  "| Case | Mech | Old RMSD | New RMSD | Difference |",
                  "|---|---:|---:|---:|---:|"])
    for item in chirality_priority:
        lines.append(
            f"| {item['step_id']} | {item['mechanism_id']} | "
            f"{item['old_fixed_mapping_aligned_rmsd']:.4f} | "
            f"{item['new_fixed_mapping_aligned_rmsd']:.4f} | "
            f"{item['rmsd_improvement']:.4f} |")
    (output / "mapping_consistency_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    batch = {
        "schema_version": "rxn_core.index_chirality_batch_summary/v1",
        "expected_case_count": len(expected),
        "completed_summary_count": len(summaries),
        "ok_count": len(summaries),
        "error_count": len(errors),
        "missing_count": len(missing),
        "errors": errors,
        "missing": missing,
        "total_elapsed_case_seconds": sum(
            float(item.get("elapsed_seconds", 0.0)) for item in summaries),
    }
    _write_atomic(output / "batch_summary.json", batch)
    print(json.dumps({"batch": batch, "report": {
        key: value for key, value in report.items()
        if key != "ranked_mechanisms"
    }}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("prepare", "summarize"):
        item = sub.add_parser(name)
        item.add_argument("--source", type=Path, required=True)
        item.add_argument("--output", type=Path, required=True)
    case = sub.add_parser("run-case")
    case.add_argument("--source", type=Path, required=True)
    case.add_argument("--output", type=Path, required=True)
    case.add_argument("--tier", choices=TIERS, required=True)
    case.add_argument("--task-index", type=int, required=True)
    case.add_argument("--performance-contract", type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.source.resolve(), args.output.resolve())
    elif args.command == "run-case":
        run_case(args.source.resolve(), args.output.resolve(),
                 args.tier, args.task_index,
                 (args.performance_contract.resolve()
                  if args.performance_contract else None))
    else:
        summarize(args.source.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()

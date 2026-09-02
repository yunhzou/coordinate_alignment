"""Strict serialization for fragment-detection records."""
from __future__ import annotations

from ..alignment.post_aam import AAMHierarchy
from .models import FragmentCandidate, FragmentDetectionResult


FRAGMENT_DETECTION_SCHEMA = "rxn_core.fragment_detection/v2"


def fragment_candidate_to_record(candidate: FragmentCandidate):
    return {
        "mapping": [list(item) for item in candidate.mapping],
        "retained_atoms": list(candidate.retained_atoms),
        "covered_target_atoms": list(candidate.covered_target_atoms),
        "leftover_fragments": [
            list(item) for item in candidate.leftover_fragments
        ],
        "boundary_bonds": [list(item) for item in candidate.boundary_bonds],
        "attachment_atoms_source": list(candidate.attachment_atoms_source),
        "attachment_atoms_target": list(candidate.attachment_atoms_target),
        "copied_residual_placements": [
            list(item) for item in candidate.copied_residual_placements
        ],
        "augmented_target_atom_count": candidate.augmented_target_atom_count,
        "retained_fragments": [
            list(item) for item in candidate.retained_fragments
        ],
        "aam_hierarchy": candidate.aam_hierarchy.to_record(),
    }


def fragment_candidate_from_record(record):
    return FragmentCandidate(
        source_id=str(record.get("source_id", "")),
        mapping=tuple(tuple(map(int, item))
                      for item in record.get("mapping") or ()),
        retained_atoms=tuple(map(int, record.get("retained_atoms") or ())),
        covered_target_atoms=tuple(map(
            int, record.get("covered_target_atoms") or ())),
        leftover_fragments=tuple(
            tuple(map(int, item))
            for item in record.get("leftover_fragments") or ()),
        boundary_bonds=tuple(
            tuple(map(int, item))
            for item in record.get("boundary_bonds") or ()),
        attachment_atoms_source=tuple(map(
            int, record.get("attachment_atoms_source") or ())),
        attachment_atoms_target=tuple(map(
            int, record.get("attachment_atoms_target") or ())),
        copied_residual_placements=tuple(
            tuple(map(int, item))
            for item in record.get("copied_residual_placements") or ()),
        augmented_target_atom_count=int(
            record.get("augmented_target_atom_count", 0)),
        retained_fragments=tuple(
            tuple(map(int, item))
            for item in record.get("retained_fragments") or ()),
        aam_hierarchy=AAMHierarchy.from_record(
            record.get("aam_hierarchy") or {}),
    )


def fragment_detection_to_record(
        result: FragmentDetectionResult, *, row_index, representation,
        candidates=None):
    selected = result.candidates if candidates is None else tuple(candidates)
    return {
        "schema": FRAGMENT_DETECTION_SCHEMA,
        "row_index": int(row_index),
        "source_id": result.source_id,
        "representation": representation,
        "status": result.status,
        "complete": result.complete,
        "branch_limit": result.branch_limit,
        "maximum_branch_count": result.maximum_branch_count,
        "capped_seed_count": result.capped_seed_count,
        "best_fragment_size": result.best_fragment_size,
        "initial_placement_encounters": result.initial_placement_encounters,
        "initial_family_count": result.initial_family_count,
        "best_initial_family_count": result.best_initial_family_count,
        "seed_attempt_count": result.seed_attempt_count,
        "seed_pruned_count": result.seed_pruned_count,
        "rough_stop_hit": result.rough_stop_hit,
        "candidates": [
            fragment_candidate_to_record(candidate) for candidate in selected
        ],
    }

"""Strict serialization for fragment-detection records."""
from __future__ import annotations

from .models import FragmentCandidate, FragmentDetectionResult


FRAGMENT_DETECTION_SCHEMA = "rxn_core.fragment_detection/v1"


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
    }


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
        "candidates": [
            fragment_candidate_to_record(candidate) for candidate in selected
        ],
    }

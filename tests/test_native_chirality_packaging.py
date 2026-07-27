from __future__ import annotations

import json

import numpy as np

from tools.build_native_chirality_source_archive import _view_mechanism
from tools.neb_support.neb_orientation_package import (
    _native_index_chirality_summary,
)


def _native_mechanism() -> dict:
    mapping = {0: 1, 1: 0, 2: 2, 3: 3, 4: 4}
    return {
        "id": 1,
        "mapping_RP": mapping,
        "branch_symmetry": {
            "index_chirality": {
                "schema_version": "rxn_core.index_chirality/v3",
                "policy": "preserve",
                "status": "applied",
                "source_mapping_sha256": "source",
                "selected_mapping_sha256": "selected",
                "source_index_chirality_violation_count": 2,
                "selected_index_chirality_violation_count": 0,
                "defined_frame_count": 4,
                "undefined_frame_count": 1,
                "switchable_r_atoms": [0, 1],
                "candidate_search": {
                    "semantics": "bounded_native_candidates",
                    "seed_route_count": 2,
                    "fragment_parity_seed_count": 2,
                    "parity_variable_count": 3,
                    "gf2_equation_count": 4,
                    "gf2_solved_route_count": 1,
                    "unique_candidate_evaluation_count": 4,
                    "candidate_evaluations": [
                        {"candidate_id": f"candidate:{index}"}
                        for index in range(4)
                    ],
                },
                "allowed_candidate_count": 1,
                "selected_candidate_id": "candidate:selected",
                "immutable_frame_count": 3,
                "immutable_source_mismatch_count": 1,
                "immutable_frames": [{
                    "id": "f:2:0-1-3-4",
                    "center_R": 2,
                    "neighbors_R_index_order": [0, 1, 3, 4],
                    "reason": (
                        "no_AAM_authorized_switchable_atom_in_frame"),
                    "source_index_chirality_mismatch": True,
                    "source_mismatch_details": [{
                        "reason": "index_orientation_reversed",
                    }],
                }],
                "mapping_changes": [{
                    "r_atom": 0,
                    "source_p_atom": 0,
                    "selected_p_atom": 1,
                }],
                "selection_rule": "bounded selection",
                "invariants": {
                    "product_automorphism_generation_used": False,
                },
            },
        },
    }


def test_v3_native_summary_uses_fragment_parity_metadata():
    summary = _native_index_chirality_summary(_native_mechanism())

    assert summary is not None
    assert summary["switchable_r_atoms"] == [0, 1]
    assert summary["candidate_search"] == {
        "semantics": "bounded_native_candidates",
        "seed_route_count": 2,
        "fragment_parity_seed_count": 2,
        "parity_variable_count": 3,
        "gf2_equation_count": 4,
        "gf2_solved_route_count": 1,
        "unique_candidate_evaluation_count": 4,
    }
    assert summary["allowed_candidate_count"] == 1
    assert summary["selected_candidate_id"] == "candidate:selected"
    assert summary["immutable_frame_count"] == 3
    assert summary["immutable_source_mismatch_count"] == 1
    assert summary["immutable_mismatch_frames"][0]["id"] == (
        "f:2:0-1-3-4")

    encoded = json.dumps(summary)
    assert "candidate_evaluations" not in encoded
    assert "symmetry_action" not in encoded
    assert "active_family" not in encoded
    assert "selected_action_id" not in encoded


def test_source_viewer_preserves_selected_mapping_without_reselection():
    mechanism = _native_mechanism()
    coords_P = np.array([
        [10.0, 0.0, 0.0],
        [20.0, 0.0, 0.0],
        [30.0, 0.0, 0.0],
        [40.0, 0.0, 0.0],
        [50.0, 0.0, 0.0],
    ])

    view = _view_mechanism(
        mechanism,
        {},
        coords_P,
        selected=True,
    )

    assert view["mapping_RP"] == [1, 0, 2, 3, 4]
    assert view["product_aam_order"] == [
        [20.0, 0.0, 0.0],
        [10.0, 0.0, 0.0],
        [30.0, 0.0, 0.0],
        [40.0, 0.0, 0.0],
        [50.0, 0.0, 0.0],
    ]
    assert view["local_fragment_R"] == [0, 1]
    assert view["native_index_chirality"][
        "selected_candidate_id"
    ] == "candidate:selected"
    assert view["symmetry_groups"] == [{
        "r_atoms": [0, 1],
        "p_atoms": [0, 1],
        "assignments": "1 chirality-safe parity candidates",
        "source": "native_index_chirality_fragment_parity",
    }]

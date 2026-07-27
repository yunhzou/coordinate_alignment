import numpy as np

import rxn_core.pipeline as pipeline
from rxn_core.alignment.index_chirality import (
    build_index_frames,
    select_index_chirality_assignment,
)


def _tetrahedral_case():
    elements = ["C", "F", "H", "H", "H"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [-0.9, -0.9, -0.9],
        [0.9, 0.9, -0.9],
        [0.9, -0.9, 0.9],
        [-0.9, 0.9, 0.9],
    ])
    wbo = np.zeros((5, 5))
    for neighbor in range(1, 5):
        wbo[0, neighbor] = wbo[neighbor, 0] = 1.0
    identity = {index: index for index in range(5)}
    odd = dict(identity)
    odd[3], odd[4] = odd[4], odd[3]
    branch_symmetry = {
        "blocks": [{
            "r_atoms": [3, 4],
            "p_atoms": [3, 4],
            "source": "alternate_witness",
        }],
    }
    return elements, coords, wbo, identity, odd, branch_symmetry


def test_index_chiral_selector_chooses_consistent_final_automorphism():
    elements, coords, wbo, identity, odd, symmetry = _tetrahedral_case()

    selection = select_index_chirality_assignment(
        odd, symmetry,
        elements, coords, wbo,
        elements, coords, wbo,
    )

    assert selection.selected_mapping == identity
    assert selection.metadata["source_index_chirality_violation_count"] == 1
    assert selection.metadata["selected_index_chirality_violation_count"] == 0
    assert selection.metadata["switchable_r_atoms"] == [3, 4]
    assert selection.metadata["candidate_count"] == 2


def test_index_frames_use_persistent_degree_four_shells():
    elements, coords, wbo, identity, _odd, _symmetry = _tetrahedral_case()

    frames, undefined = build_index_frames(
        identity, coords, coords, wbo, wbo)

    assert undefined == ()
    assert len(frames) == 1
    assert frames[0].center_R == 0
    assert frames[0].neighbors_R == (1, 2, 3, 4)


def test_immutable_orientation_mismatch_is_diagnostic_only():
    elements, coords_R, wbo, identity, _odd, _symmetry = _tetrahedral_case()
    coords_P = coords_R.copy()
    coords_P[[3, 4]] = coords_P[[4, 3]]

    selection = select_index_chirality_assignment(
        identity, {},
        elements, coords_R, wbo,
        elements, coords_P, wbo,
    )

    assert selection.selected_mapping == identity
    assert selection.metadata["defined_frame_count"] == 0
    assert selection.metadata["immutable_frame_count"] == 1
    assert selection.metadata["immutable_source_mismatch_count"] == 1


def test_rp_stage_index_chiral_mode_is_post_processing():
    elements, coords, wbo, identity, odd, symmetry = _tetrahedral_case()
    inputs = pipeline.step_inputs_from_arrays(
        "index_chiral", elements, coords, wbo,
        elements, coords, wbo)
    pool = {
        ((), ()): {
            "mapping": odd,
            "cuts": frozenset(),
            "has_no_cut": True,
            "dedup_count": 1,
            "branch_symmetry": symmetry,
        },
    }
    config = pipeline.rp_stage_config()
    config["index_chirality"] = "preserve"

    result = pipeline.run_rp_stage_from_pool(inputs, pool, config=config)
    mechanism = result["mechanisms"][0]

    assert mechanism["mapping_RP"] == identity
    assert mechanism["index_chirality"]["status"] == "applied"
    assert mechanism["broken_bonds_R"] == []
    assert mechanism["formed_bonds_R"] == []


def test_rp_stage_index_chiral_off_is_noop(monkeypatch):
    elements, coords, wbo, _identity, odd, symmetry = _tetrahedral_case()
    inputs = pipeline.step_inputs_from_arrays(
        "index_chiral_off", elements, coords, wbo,
        elements, coords, wbo)
    pool = {
        ((), ()): {
            "mapping": odd,
            "cuts": frozenset(),
            "has_no_cut": True,
            "dedup_count": 1,
            "branch_symmetry": symmetry,
        },
    }

    def forbidden(*_args, **_kwargs):
        raise AssertionError("index chirality ran while disabled")

    monkeypatch.setattr(pipeline, "select_index_chirality_assignment", forbidden)
    config = pipeline.rp_stage_config()
    config["index_chirality"] = "off"

    result = pipeline.run_rp_stage_from_pool(inputs, pool, config=config)

    assert result["mechanisms"][0]["mapping_RP"] == odd
    assert result["mechanisms"][0]["index_chirality"] is None

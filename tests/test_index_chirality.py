import numpy as np
import pytest

import rxn_core.pipeline as pipeline
from rxn_core.alignment.index_chirality import (
    IndexChiralityConflict,
    select_group_chiral_witness,
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
            "source": "exact_automorph_group",
        }],
    }
    return elements, coords, wbo, identity, odd, branch_symmetry


def _higher_coordinate_branch_case():
    elements = ["Sc", "O", "O", "O", "O", "O", "O"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.2, 0.1],
        [-0.4, 1.1, 0.3],
        [-0.6, -0.3, 1.2],
        [0.2, -1.2, -0.4],
        [0.7, 0.4, -1.0],
        [-1.0, 0.5, -0.7],
    ])
    wbo = np.zeros((7, 7))
    for neighbor in range(1, 7):
        wbo[0, neighbor] = wbo[neighbor, 0] = 0.5
    identity = {index: index for index in range(7)}
    reversed_witness = dict(identity)
    reversed_witness[1], reversed_witness[2] = 2, 1
    local_symmetry = {
        "fragments": [{"fragment_index": 0, "fragment": list(range(7))}],
        "blocks": [],
    }
    witnesses = [
        {"mapping": reversed_witness, "cut": [],
         "local_symmetry": local_symmetry},
        {"mapping": identity, "cut": [],
         "local_symmetry": local_symmetry},
    ]
    return elements, coords, wbo, identity, reversed_witness, witnesses


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
    assert selection.metadata["switchable_r_atoms"] == [2, 3, 4]
    assert selection.metadata["solver"] == (
        "pynauty_colored_relational_isomorphism")


def test_exact_symmetry_freedom_corrects_orientation_without_display_block():
    elements, coords_R, wbo, identity, _odd, _symmetry = _tetrahedral_case()
    coords_P = coords_R.copy()
    coords_P[[3, 4]] = coords_P[[4, 3]]

    selection = select_index_chirality_assignment(
        identity, {},
        elements, coords_R, wbo,
        elements, coords_P, wbo,
    )

    assert selection.selected_mapping != identity
    assert selection.selected_mapping[3] == 4
    assert selection.selected_mapping[4] == 3
    assert selection.metadata["defined_frame_count"] == 1
    assert selection.metadata["selected_index_chirality_violation_count"] == 0


def test_group_chirality_selects_between_distinct_aam_witnesses():
    elements, coords, wbo, identity, reversed_witness, witnesses = (
        _higher_coordinate_branch_case())

    selection = select_group_chiral_witness(
        reversed_witness, witnesses,
        elements, coords, wbo, elements, coords, wbo)

    assert selection.selected_mapping == identity
    assert selection.selected_witness_index == 1
    assert selection.metadata["preserved_frame_count"] > 0
    assert selection.metadata["reversed_frame_count"] == 0


def test_group_orientation_keeps_a_definite_near_planar_sign():
    elements = ["Sc", "O", "O", "O", "O", "O"]
    coords_R = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0], [-1.0, -0.2, 0.3], [0.2, -1.0, -0.4],
    ])
    coords_P = coords_R.copy()
    coords_P[3] = [0.5, 0.5, 0.01]
    wbo = np.zeros((6, 6))
    for neighbor in range(1, 6):
        wbo[0, neighbor] = wbo[neighbor, 0] = 0.5
    identity = {index: index for index in range(6)}
    witnesses = [{"mapping": identity}]

    selection = select_group_chiral_witness(
        identity, witnesses,
        elements, coords_R, wbo, elements, coords_P, wbo)

    frame = next(
        frame for frame in selection.preserved_frames
        if frame["center_R"] == 0
        and frame["neighbors_R_index_order"] == [1, 2, 3])
    assert frame["product_orientation_sign"] == 1
    assert 0 < frame["product_normalized_orientation"] < 0.1


def test_preserved_group_chirality_survives_final_automorphism():
    elements, coords, wbo, identity, reversed_witness, witnesses = (
        _higher_coordinate_branch_case())
    group = select_group_chiral_witness(
        reversed_witness, witnesses,
        elements, coords, wbo, elements, coords, wbo)

    selection = select_index_chirality_assignment(
        group.selected_mapping, witnesses[1]["local_symmetry"],
        elements, coords, wbo, elements, coords, wbo,
        group_chirality_frames=group.preserved_frames)

    assert selection.metadata["preserved_group_chirality_frame_count"] == (
        len(group.preserved_frames))
    assert selection.metadata["selected_index_chirality_violation_count"] == 0


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


def test_rp_stage_selects_group_chiral_branch_before_automorphism():
    elements, coords, wbo, identity, reversed_witness, witnesses = (
        _higher_coordinate_branch_case())
    inputs = pipeline.step_inputs_from_arrays(
        "group_chiral_branch", elements, coords, wbo,
        elements, coords, wbo)
    pool = {
        ((), ()): {
            "mapping": reversed_witness,
            "cuts": frozenset(),
            "has_no_cut": True,
            "dedup_count": 2,
            "branch_symmetry": {
                "witnesses": witnesses,
                "fragments": witnesses[0]["local_symmetry"]["fragments"],
            },
        },
    }
    config = pipeline.rp_stage_config()
    config["index_chirality"] = "preserve"

    result = pipeline.run_rp_stage_from_pool(inputs, pool, config=config)

    mechanism = result["mechanisms"][0]
    group = mechanism["index_chirality"]["group_chirality_witness"]
    assert group["selected_witness_index"] == 1
    assert group["reversed_frame_count"] == 0
    assert mechanism["index_chirality"][
        "preserved_group_chirality_frame_count"] > 0


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


def test_rp_stage_fails_only_after_every_mechanism_conflicts(monkeypatch):
    elements, coords, wbo, _identity, odd, symmetry = _tetrahedral_case()
    inputs = pipeline.step_inputs_from_arrays(
        "index_chiral_conflict", elements, coords, wbo,
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

    def conflict(*_args, **_kwargs):
        raise IndexChiralityConflict("no exact consensus")

    monkeypatch.setattr(
        pipeline, "select_index_chirality_assignment", conflict)
    config = pipeline.rp_stage_config()
    config["index_chirality"] = "preserve"

    with pytest.raises(
            IndexChiralityConflict,
            match="all 1 minimum-event mechanisms failed"):
        pipeline.run_rp_stage_from_pool(inputs, pool, config=config)


def test_rp_stage_evaluates_mechanisms_independently(monkeypatch):
    elements, coords, wbo, _identity, odd, symmetry = _tetrahedral_case()
    inputs = pipeline.step_inputs_from_arrays(
        "index_chiral_independent", elements, coords, wbo,
        elements, coords, wbo)
    entry = {
        "mapping": odd,
        "cuts": frozenset(),
        "has_no_cut": True,
        "dedup_count": 1,
        "branch_symmetry": symmetry,
    }
    pool = {
        (((0, 1),), ()): dict(entry),
        ((), ((0, 1),)): dict(entry),
    }
    actual_selector = pipeline.select_index_chirality_assignment
    calls = 0

    def first_conflicts(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise IndexChiralityConflict("mechanism-specific conflict")
        return actual_selector(*args, **kwargs)

    monkeypatch.setattr(
        pipeline, "select_index_chirality_assignment", first_conflicts)
    config = pipeline.rp_stage_config()
    config["index_chirality"] = "preserve"

    result = pipeline.run_rp_stage_from_pool(inputs, pool, config=config)

    assert calls == 2
    assert len(result["mechanisms"]) == 1
    rejected = result["rejected_index_chirality"]
    assert len(rejected) == 1
    assert rejected[0]["source_mechanism_id"] == 1
    assert rejected[0]["reason"] == "mechanism-specific conflict"
    assert rejected[0]["source_mapping_RP"] == odd


def test_index_direction_assignment_moves_a_coupled_whole_arm():
    elements = ["P", "F", "C", "H", "C", "H", "C", "H"]
    coords = np.array([
        [0.0, 0.0, 0.0], [0.0, 0.0, -1.2],
        [1.0, 0.0, 0.7], [1.7, 0.0, 1.2],
        [-0.5, 0.866, 0.7], [-0.85, 1.472, 1.2],
        [-0.5, -0.866, 0.7], [-0.85, -1.472, 1.2],
    ])
    wbo = np.zeros((8, 8))
    for left, right in [(0, 1), (0, 2), (2, 3),
                        (0, 4), (4, 5), (0, 6), (6, 7)]:
        wbo[left, right] = wbo[right, left] = 1.0
    source = {index: index for index in range(8)}
    source[2], source[4] = 4, 2
    source[3], source[5] = 5, 3

    selection = select_index_chirality_assignment(
        source, {}, elements, coords, wbo, elements, coords, wbo)

    assert selection.selected_mapping != source
    assert selection.metadata["selected_index_chirality_violation_count"] == 0
    changed = {item["r_atom"] for item in selection.metadata["mapping_changes"]}
    assert changed
    # Every arm is moved as one exact graph automorphism; a carbon can never
    # be detached from its hydrogen merely to repair the center orientation.
    for carbon, hydrogen in ((2, 3), (4, 5), (6, 7)):
        assert selection.selected_mapping[hydrogen] == (
            selection.selected_mapping[carbon] + 1)


def test_index_direction_is_coordination_number_independent():
    elements = ["C", "H", "H", "H", "H", "H"]
    coords_R = np.array([
        [0.0, 0.0, 0.0],
        [-0.913346, -0.402452, 0.061897],
        [-0.949395, 0.292956, -0.113251],
        [0.746586, 0.150256, -0.648100],
        [-0.839309, -0.530957, -0.116811],
        [0.258909, -0.548234, 0.795239],
    ])
    coords_P = coords_R.copy()
    coords_P[[2, 5]] = coords_P[[5, 2]]
    wbo = np.zeros((6, 6))
    for neighbor in range(1, 6):
        wbo[0, neighbor] = wbo[neighbor, 0] = 1.0
    identity = {index: index for index in range(6)}

    selection = select_index_chirality_assignment(
        identity, {}, elements, coords_R, wbo,
        elements, coords_P, wbo)

    assert selection.metadata["selected_index_chirality_violation_count"] == 0
    assert len(selection.metadata["active_frames"]) == 5
    assert selection.metadata["defined_frame_count"] >= 1
    assert selection.selected_mapping != identity


def test_planarized_endpoint_has_no_chirality_constraint():
    elements = ["C", "C", "H", "H"]
    coords_R = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 1.0],
        [0.0, -1.0, 1.0],
    ])
    coords_P = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ])
    wbo = np.zeros((4, 4))
    for neighbor in (1, 2, 3):
        wbo[0, neighbor] = wbo[neighbor, 0] = 1.0
    identity = {index: index for index in range(4)}

    selection = select_index_chirality_assignment(
        identity, {}, elements, coords_R, wbo,
        elements, coords_P, wbo)

    assert selection.metadata["selected_index_chirality_violation_count"] == 0
    assert selection.metadata["defined_frame_count"] == 0
    assert selection.metadata["nonstereogenic_frame_count"] == 1


def test_one_coplanar_high_coordinate_simplex_does_not_erase_the_others():
    elements = ["C", "H", "H", "H", "H", "F"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.2, 0.3, 1.0],
    ])
    wbo = np.zeros((6, 6))
    for neighbor in range(1, 6):
        wbo[0, neighbor] = wbo[neighbor, 0] = 1.0
    identity = {index: index for index in range(6)}

    selection = select_index_chirality_assignment(
        identity, {}, elements, coords, wbo,
        elements, coords, wbo)

    assert selection.metadata["nonstereogenic_frame_count"] == 1
    assert len(selection.metadata["active_frames"]) == 4
    assert {frame["center_R"] for frame in selection.metadata["active_frames"]} == {0}

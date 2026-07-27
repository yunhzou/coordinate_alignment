import copy
import inspect

import numpy as np
import pytest

from tools.neb_support.neb_orientation import (
    CandidateFamily,
    EncodedCandidate,
    build_orientation_frames,
    build_candidate_family,
    evaluate_orientations,
    mapping_sha256,
    normalized_orientation,
    optimize_neb_orientation,
    proper_kabsch,
    select_neb_mapping,
)
from tools.run_neb_orientation_batch import _parser


def _anchored_tetrahedral_endpoint(
    source_mapping=None,
    *,
    core_atoms=(),
    exact_fixed=(),
):
    """Return a methyl-like tetrahedral frame with non-coplanar fit anchors."""
    elements = ["C", "C", "H", "H", "H", "O", "N"]
    coords = np.array([
        [0.0, 0.0, 0.0],       # tetrahedral center
        [-0.9, -0.9, -0.9],    # persistent scaffold neighbor
        [0.9, 0.9, -0.9],      # equivalent H
        [0.9, -0.9, 0.9],      # equivalent H
        [-0.9, 0.9, 0.9],      # equivalent H
        [-2.0, -0.7, -0.8],    # fixed heavy-atom fit anchor
        [-0.7, -2.0, -0.8],    # fixed heavy-atom fit anchor
    ])
    wbo = np.zeros((len(elements), len(elements)))
    for left, right in (
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 5),
        (1, 6),
    ):
        wbo[left, right] = wbo[right, left] = 1.0

    source = (
        {index: index for index in range(len(elements))}
        if source_mapping is None
        else dict(source_mapping)
    )
    valid_block = {
        "source": "sym_block",
        "fragment_index": 0,
        "block_index": 7,
        "r_atoms": [2, 3, 4],
        "p_atoms": [2, 3, 4],
        "assignments": "3!",
        # This is a growth-state flag, not postprocessing authorization.
        "extendable": False,
    }
    nested_valid_block = {
        key: value for key, value in valid_block.items()
        if key != "source"
    }
    lossy_local_summary = {
        "source": "island_automorph",
        "fragment_index": 1,
        "block_index": 8,
        "r_atoms": [5, 6],
        "p_atoms": [5, 6],
        "assignments": "2!",
        "extendable": True,
    }
    mechanism = {
        "mapping_RP": source,
        "core_atoms": list(core_atoms),
        "branch_symmetry": {
            # Aggregate summaries must not enlarge the exact witness family.
            "blocks": [{
                "source": "sym_block",
                "r_atoms": [5, 6],
                "p_atoms": [5, 6],
                "assignments": "2!",
            }],
            "color_groups": [{
                "r_atoms": [5, 6],
                "p_atoms": [5, 6],
            }],
            "witnesses": [{
                "mapping": source,
                "local_symmetry": {
                    # This flattened list is display/provenance data only.
                    "blocks": [valid_block, lossy_local_summary],
                    "fragments": [{
                        "symmetry": {
                            "exact_fixed": list(exact_fixed),
                            "blocks": [nested_valid_block],
                        },
                    }, {
                        "symmetry": {
                            "exact_fixed": [],
                            "blocks": [lossy_local_summary],
                        },
                    }],
                },
            }],
        },
    }
    return elements, coords, wbo, mechanism


def _optimize(mechanism, elements, coords, wbo):
    return optimize_neb_orientation(
        mechanism,
        elements,
        coords,
        wbo,
        elements,
        coords.copy(),
        wbo.copy(),
    )


def test_odd_h_swap_is_repaired_without_mutating_source_mechanism():
    source = {index: index for index in range(7)}
    source[2], source[3] = source[3], source[2]
    elements, coords, wbo, mechanism = _anchored_tetrahedral_endpoint(
        source)
    mechanism_before = copy.deepcopy(mechanism)

    result = _optimize(mechanism, elements, coords, wbo)

    assert mechanism == mechanism_before
    assert result.source_mapping == source
    assert result.source_violation_count > 0
    assert result.selected_mapping == {index: index for index in range(7)}
    assert result.final_violation_count == 0
    assert result.transform.determinant == pytest.approx(1.0, abs=1e-12)


def test_even_cyclic_h_permutation_uses_zero_motion_tiebreak():
    source = {index: index for index in range(7)}
    source.update({2: 3, 3: 4, 4: 2})
    elements, coords, wbo, mechanism = _anchored_tetrahedral_endpoint(
        source)

    result = _optimize(mechanism, elements, coords, wbo)

    # The cyclic source assignment already has even local parity.  Geometry
    # resolves the equally oriented encoded alternatives in favor of zero
    # endpoint atom motion.
    assert result.source_violation_count == 0
    assert result.selected_mapping == {index: index for index in range(7)}
    assert result.final_violation_count == 0
    assert result.max_mutable_displacement == pytest.approx(0.0, abs=1e-12)
    assert result.mutable_rmsd == pytest.approx(0.0, abs=1e-12)


def test_native_core_index_chirality_assignment_is_not_remapped_downstream():
    source = {index: index for index in range(7)}
    source.update({2: 3, 3: 4, 4: 2})
    elements, coords, wbo, mechanism = _anchored_tetrahedral_endpoint(
        source)
    mechanism["branch_symmetry"]["index_chirality"] = {
        "schema_version": "rxn_core.index_chirality/v3",
        "status": "applied",
        "selected_mapping_sha256": mapping_sha256(source),
        "selected_candidate_id": "candidate:test",
        "frames": [{
            "center_R": 0,
            "neighbors_R_index_order": [1, 2, 3, 4],
            "reactant_normalized_orientation": normalized_orientation(
                coords, 1, (2, 3, 4)),
        }],
        "undefined_frame_count": 3,
    }

    result = _optimize(mechanism, elements, coords, wbo)

    # Without native-core metadata, the downstream geometry tie-break chooses
    # identity for this endpoint.  Once the core has selected an action, the
    # helper may fit/package coordinates but must not select a second mapping.
    assert result.selected_mapping == source
    assert result.family.witness_index == -1
    assert result.family.blocks == ()
    assert result.family.fixed_r_atoms == tuple(range(7))
    assert len(result.family.candidates) == 1
    assert len(result.frames) == 1
    assert result.frames[0].center == 0
    assert result.frames[0].neighbors == (1, 2, 3, 4)
    assert result.frames[0].to_dict()["orientation_model"] == (
        "affine_four_neighbor_tetrahedron")
    assert result.undefined_frame_count == 3
    assert result.selected_candidate.provenance_paths == ({
        "source": "native_core_index_chirality_assignment",
        "native_index_chirality_status": "applied",
        "selected_candidate_id": "candidate:test",
    },)
    assert result.to_dict()["invariants"][
        "selected_mapping_is_native_core_assignment"
    ] is True


def test_native_v2_does_not_reinvent_degree_three_constraints():
    elements = ["C", "H", "H", "H"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    wbo = np.zeros((4, 4))
    wbo[0, 1:] = 1.0
    wbo[1:, 0] = 1.0
    source = {0: 0, 1: 2, 2: 1, 3: 3}
    mechanism = {
        "mapping_RP": source,
        "branch_symmetry": {
            "index_chirality": {
                "schema_version": "rxn_core.index_chirality/v2",
                "status": "applied",
                "selected_mapping_sha256": mapping_sha256(source),
                "selected_candidate_id": "candidate:degree-three",
                # Native v2 intentionally has no hard frame at a
                # three-coordinate centre.
                "frames": [],
                "undefined_frame_count": 5,
            },
        },
    }
    family = build_candidate_family(mechanism, elements, elements)

    # The legacy downstream policy would invent one frame and reject the
    # already-certified native mapping.
    legacy_frames, _ = build_orientation_frames(
        family, coords, coords.copy(), wbo, wbo.copy())
    legacy_violations, _ = evaluate_orientations(
        source, legacy_frames, coords)
    assert len(legacy_frames) == 1
    assert legacy_violations == 1

    result = optimize_neb_orientation(
        mechanism,
        elements,
        coords,
        wbo,
        elements,
        coords.copy(),
        wbo.copy(),
    )

    assert result.selected_mapping == source
    assert result.frames == ()
    assert result.undefined_frame_count == 5
    assert result.source_violation_count == 0
    assert result.final_violation_count == 0


def test_core_label_does_not_override_encoded_aam_family():
    source = {index: index for index in range(7)}
    source[3], source[4] = source[4], source[3]
    elements, coords, wbo, mechanism = _anchored_tetrahedral_endpoint(
        source, core_atoms=[2])

    result = _optimize(mechanism, elements, coords, wbo)

    assert result.family.blocks[0].r_atoms == (2, 3, 4)
    assert 2 not in result.family.fixed_r_atoms
    assert result.selected_mapping == {index: index for index in range(7)}


def test_exact_fixed_pair_is_removed_from_encoded_shuffle():
    source = {index: index for index in range(7)}
    source[3], source[4] = source[4], source[3]
    elements, coords, wbo, mechanism = _anchored_tetrahedral_endpoint(
        source, exact_fixed=[2])

    result = _optimize(mechanism, elements, coords, wbo)

    assert result.family.blocks[0].r_atoms == (3, 4)
    assert 2 in result.family.fixed_r_atoms
    assert result.selected_mapping[2] == source[2]
    assert result.selected_mapping[3] == 3
    assert result.selected_mapping[4] == 4


def test_extendable_false_is_accepted_and_lossy_summaries_are_ignored():
    elements, _coords, _wbo, mechanism = _anchored_tetrahedral_endpoint()

    family = build_candidate_family(mechanism, elements, elements)

    assert len(family.blocks) == 1
    assert family.blocks[0].r_atoms == (2, 3, 4)
    assert family.blocks[0].p_atoms == (2, 3, 4)
    assert all(
        set(block.r_atoms) != {5, 6}
        for block in family.blocks
    )
    assert any(
        record["reason"] == "lossy_summary_is_not_candidate_authorization"
        and record["r_atoms"] == [5, 6]
        for record in family.discarded_blocks
    )


def test_flattened_factorial_does_not_complete_an_open_nested_state():
    elements, _coords, _wbo, mechanism = _anchored_tetrahedral_endpoint()
    nested = mechanism["branch_symmetry"]["witnesses"][0][
        "local_symmetry"]["fragments"][0]["symmetry"]["blocks"][0]
    nested.update({
        "r_atoms": [2],
        "p_atoms": [2, 3, 4],
        "assignments": "3",
        "open": True,
    })

    family = build_candidate_family(mechanism, elements, elements)

    assert family.blocks == ()
    assert any(
        record["reason"] == "open_symmetry_state_is_not_a_complete_shuffle"
        for record in family.discarded_blocks
    )


def test_correlated_fragment_alternate_is_used_without_core_override():
    source = {index: index for index in range(7)}
    source[2], source[3] = source[3], source[2]
    identity = {index: index for index in range(7)}
    elements, coords, wbo, mechanism = _anchored_tetrahedral_endpoint(
        source, core_atoms=[2])
    symmetry = mechanism["branch_symmetry"]["witnesses"][0][
        "local_symmetry"]["fragments"][0]["symmetry"]
    symmetry["witness"] = source
    symmetry["alternates"] = [{
        "witness": identity,
        "multiplicity": 1,
    }]
    symmetry["blocks"] = []

    result = _optimize(mechanism, elements, coords, wbo)

    assert result.family.blocks == ()
    assert len(result.family.candidates) == 2
    assert result.source_violation_count > 0
    assert result.selected_mapping == identity
    assert result.final_violation_count == 0
    assert any(
        choice["choice"] == "alternate"
        for path in result.selected_candidate.provenance_paths
        for choice in path["fragment_choices"]
    )


def test_fixed_core_center_still_constrains_mutable_ligand_orientation():
    source = {index: index for index in range(7)}
    source[2], source[3] = source[3], source[2]
    elements, coords, wbo, mechanism = _anchored_tetrahedral_endpoint(
        source, core_atoms=[0])

    result = _optimize(mechanism, elements, coords, wbo)

    assert result.source_violation_count > 0
    assert result.selected_mapping == {index: index for index in range(7)}
    assert result.selected_mapping[0] == source[0]


def test_proper_kabsch_does_not_fit_a_reflection():
    reference = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 3.0],
    ])
    reflection = np.diag([-1.0, 1.0, 1.0])
    mobile = reference @ reflection
    assert np.allclose(mobile @ reflection, reference)
    assert np.linalg.det(reflection) == pytest.approx(-1.0)

    rotation, translation, rmsd, rank = proper_kabsch(reference, mobile)

    assert rank == 3
    assert np.linalg.det(rotation) == pytest.approx(1.0, abs=1e-12)
    assert rmsd > 0.5
    assert not np.allclose(mobile @ rotation + translation, reference)


def test_orientation_defined_by_roundoff_not_empirical_volume_cutoff():
    for function in (
        build_orientation_frames,
        evaluate_orientations,
        select_neb_mapping,
        optimize_neb_orientation,
    ):
        assert "volume_tolerance" not in inspect.signature(
            function).parameters
    assert "--volume-tolerance" not in _parser()._option_string_actions

    source = {index: index for index in range(4)}
    family = CandidateFamily(
        source_mapping=source,
        witness_index=-1,
        candidates=(EncodedCandidate(source, ()),),
        blocks=(),
        fixed_r_atoms=tuple(source),
        discarded_blocks=(),
    )
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 1.0e-8],
    ])
    wbo = np.zeros((4, 4))
    wbo[0, 1:] = 1.0
    wbo[1:, 0] = 1.0

    frames, undefined = build_orientation_frames(
        family, coords, coords.copy(), wbo, wbo.copy())

    assert undefined == 0
    assert len(frames) == 1
    # This would have been discarded by the former 0.05 empirical cutoff,
    # despite its determinant being unambiguous at machine precision.
    assert 0.0 < abs(frames[0].reactant_orientation) < 0.05

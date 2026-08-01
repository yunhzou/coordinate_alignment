import numpy as np
import pytest

import rxn_core.alignment.index_chirality as index_chirality_module
import rxn_core.pipeline as pipeline
from rxn_core.alignment.index_chirality import (
    IndexChiralityConflict,
    _fixed_mappings_aligned_rmsd,
    _generated_atom_permutations,
    _masked_relation_data,
    _minimum_rmsd_group_action,
    analytical_family_static_context,
    compile_analytical_mapping_family,
    fixed_mapping_aligned_rmsd,
    select_group_chiral_witness,
    select_index_chirality_assignment,
)


def test_batched_fixed_mapping_rmsd_is_scalar_equivalent():
    rng = np.random.default_rng(20260729)
    coords_R = rng.normal(size=(9, 3))
    coords_P = rng.normal(size=(9, 3))
    permutations = [rng.permutation(9) for _ in range(12)]
    mappings = [
        {r: int(permutation[r]) for r in range(9)}
        for permutation in permutations
    ]
    scalar = np.asarray([
        fixed_mapping_aligned_rmsd(mapping, coords_R, coords_P)
        for mapping in mappings
    ])
    batched = _fixed_mappings_aligned_rmsd(
        mappings, coords_R, coords_P)
    np.testing.assert_allclose(batched, scalar, rtol=1e-13, atol=1e-13)


def test_symmetry_factor_rmsd_search_matches_exhaustive_group():
    rng = np.random.default_rng(41)
    coords_R = rng.normal(size=(8, 3))
    coords_P = rng.normal(size=(8, 3))
    identity = tuple(range(8))
    cycle_012 = (1, 2, 0, 3, 4, 5, 6, 7)
    swap_34 = (0, 1, 2, 4, 3, 5, 6, 7)
    swap_56 = (0, 1, 2, 3, 4, 6, 5, 7)
    generators = (cycle_012, swap_34, swap_56)
    canonical = {atom: atom for atom in range(8)}

    selected, rmsd, search = _minimum_rmsd_group_action(
        canonical, generators, coords_R, coords_P)
    exhaustive = []
    for action in _generated_atom_permutations(generators, 8):
        mapping = {r: int(action[r]) for r in range(8)}
        value = fixed_mapping_aligned_rmsd(mapping, coords_R, coords_P)
        exhaustive.append((round(value, 12), tuple(mapping.values()),
                           value, mapping))
    expected = min(exhaustive, key=lambda item: item[:2])

    assert selected == expected[3]
    assert rmsd == pytest.approx(expected[2], abs=1e-13)
    assert search["group_order"] == len(exhaustive) == 12
    assert (search["evaluated_leaf_count"]
            + search["pruned_leaf_count"]) == 12


def test_large_symmetry_branch_and_bound_matches_exhaustive_group():
    rng = np.random.default_rng(407088)
    atom_count = 30
    coords_R = rng.normal(size=(atom_count, 3))
    coords_P = coords_R + rng.normal(scale=0.03, size=(atom_count, 3))
    generators = []
    for left in range(4, atom_count, 2):
        permutation = list(range(atom_count))
        permutation[left], permutation[left + 1] = left + 1, left
        generators.append(tuple(permutation))
    canonical = {atom: atom for atom in range(atom_count)}

    selected, rmsd, search = _minimum_rmsd_group_action(
        canonical, generators, coords_R, coords_P)
    actions = _generated_atom_permutations(generators, atom_count)
    mappings = [{r: int(action[r]) for r in range(atom_count)}
                for action in actions]
    rmsds = _fixed_mappings_aligned_rmsd(mappings, coords_R, coords_P)
    exhaustive = min((
        round(float(value), 12), tuple(mapping.values()),
        float(value), mapping,
    ) for mapping, value in zip(mappings, rmsds))

    assert selected == exhaustive[3]
    assert rmsd == pytest.approx(exhaustive[2], abs=1e-13)
    assert search["group_order"] == 8192
    assert (search["evaluated_leaf_count"]
            + search["pruned_leaf_count"]) == 8192
    assert search["pruned_leaf_count"] > 0
    assert search["evaluated_leaf_count"] < 16
    assert search["search_method"] == "exact_covariance_action_ball_tree"


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


def test_finalized_aam_branch_family_is_used_without_reconstruction(
        monkeypatch):
    elements, coords, wbo, identity, odd, _symmetry = _tetrahedral_case()
    swap = list(range(len(elements)))
    swap[3], swap[4] = swap[4], swap[3]
    hierarchy = {"fragments": [{
        "fragment_index": 0,
        "fragment": list(range(len(elements))),
        "symmetry": {"automorph_generators": []},
    }]}

    def forbidden(*_args, **_kwargs):
        raise AssertionError(
            "finalized AAM groups must not rebuild a relational graph")

    monkeypatch.setattr(
        index_chirality_module, "_masked_relation_data", forbidden)
    selection = select_index_chirality_assignment(
        odd, hierarchy,
        elements, coords, wbo,
        elements, coords, wbo,
        branch_family_mappings=[odd, identity],
        aam_family_generators=[swap])

    assert selection.selected_mapping == identity
    assert selection.metadata["solver"] == "stored_AAM_generator_group"
    assert selection.metadata["selected_index_chirality_violation_count"] == 0


def test_entangled_group_reuses_compiled_aam_relation(monkeypatch):
    elements, coords, wbo, identity, odd, symmetry = _tetrahedral_case()
    family = compile_analytical_mapping_family(
        odd, symmetry, elements, wbo, elements, wbo)

    def forbidden(*_args, **_kwargs):
        raise AssertionError("compiled AAM relation must not be reconstructed")

    monkeypatch.setattr(
        index_chirality_module, "DIRECT_ACTION_COMPONENT_MAX_ORDER", 1)
    monkeypatch.setattr(
        index_chirality_module, "_masked_relation_data", forbidden)
    selection = select_index_chirality_assignment(
        odd, symmetry,
        elements, coords, wbo,
        elements, coords, wbo,
        aam_family_generators=family.target_generators,
        compiled_aam_family=family)

    assert selection.selected_mapping == identity
    assert selection.metadata["solver"] == (
        "compiled_AAM_relation_chirality_subgroup")


def test_compiled_family_cannot_drop_an_incompatible_ordinary_frame(
        monkeypatch):
    elements, coords, wbo, _identity, odd, symmetry = _tetrahedral_case()
    family = compile_analytical_mapping_family(
        odd, symmetry, elements, wbo, elements, wbo)
    base_vertices = family.graph_A.number_of_vertices
    canonical_isomorphism = index_chirality_module._canonical_isomorphism

    def reject_oriented_relation(graph_A, graph_B):
        # Model an exact AAM family whose base relation is valid but whose
        # orientation-colored subgroup is empty.
        if graph_A.number_of_vertices > base_vertices:
            return None
        return canonical_isomorphism(graph_A, graph_B)

    monkeypatch.setattr(
        index_chirality_module, "_canonical_isomorphism",
        reject_oriented_relation)

    with pytest.raises(IndexChiralityConflict) as caught:
        select_index_chirality_assignment(
            odd, symmetry,
            elements, coords, wbo,
            elements, coords, wbo,
            compiled_aam_family=family)

    diagnostics = caught.value.diagnostics
    assert diagnostics["constraint_model"] == (
        "hard_ordinary_affine_substituent_simplex")
    assert diagnostics["center_R"] == [0]


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


def test_group_chirality_rmsd_uses_each_fixed_candidate_mapping():
    elements = ["C"] * 4
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.1, 0.0, 0.0],
        [0.2, 1.7, 0.0],
        [0.3, 0.4, 2.3],
    ])
    wbo = np.zeros((4, 4))
    identity = {atom: atom for atom in range(4)}
    swapped = dict(identity)
    swapped[0], swapped[1] = 1, 0
    witnesses = [{"mapping": swapped}, {"mapping": identity}]

    selection = select_group_chiral_witness(
        swapped, witnesses,
        elements, coords, wbo, elements, coords, wbo)

    assert selection.selected_mapping == identity
    assert selection.selected_witness_index == 1
    assert selection.metadata["selected_fixed_mapping_aligned_rmsd"] == (
        pytest.approx(0.0, abs=1e-12))
    assert selection.metadata["rmsd_policy"] == (
        "exact_mapping_then_proper_rigid_fit_no_permutation")


def test_index_chirality_scores_every_valid_atom_action_before_rmsd_choice():
    elements = ["C"] * 4
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.1, 0.0, 0.0],
        [0.2, 1.7, 0.0],
        [0.3, 0.4, 2.3],
    ])
    wbo = np.zeros((4, 4))
    source = {0: 1, 1: 0, 2: 2, 3: 3}
    hierarchy = {"fragments": [{"fragment": [0, 1, 2, 3]}]}

    selection = select_index_chirality_assignment(
        source, hierarchy,
        elements, coords, wbo, elements, coords, wbo)

    assert selection.selected_mapping == {atom: atom for atom in range(4)}
    assert selection.metadata["chirality_valid_atom_bijection_count"] == 24
    assert selection.metadata["rmsd_candidate_count"] == 24
    assert selection.metadata["selected_fixed_mapping_aligned_rmsd"] == (
        pytest.approx(0.0, abs=1e-12))


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


def test_stored_group_keeps_maximal_feasible_high_coordinate_frame_basis():
    elements, coords_R, wbo, identity, _reversed, _witnesses = (
        _higher_coordinate_branch_case())
    coords_P = coords_R.copy()
    coords_P[:, 0] *= -1.0
    swap = list(range(len(elements)))
    swap[1], swap[2] = swap[2], swap[1]
    hierarchy = {"fragments": [{
        "fragment_index": 0,
        "fragment": list(range(len(elements))),
        "symmetry": {"automorph_generators": [swap]},
    }]}

    selection = select_index_chirality_assignment(
        identity, hierarchy,
        elements, coords_R, wbo,
        elements, coords_P, wbo)

    assert selection.metadata["solver"] == "stored_AAM_generator_group"
    assert selection.metadata["defined_frame_count"] == 8
    assert selection.metadata["reconfigured_frame_count"] == 7
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
            "branch_symmetry": witnesses[0]["local_symmetry"],
            "branches": [{
                "mapping": witness["mapping"],
                "cuts": witness["cut"],
                "encounter_count": 1,
                "hierarchy": witness["local_symmetry"],
            } for witness in witnesses],
        },
    }
    config = pipeline.rp_stage_config()
    config["index_chirality"] = "preserve"

    result = pipeline.run_rp_stage_from_pool(inputs, pool, config=config)

    mechanism = result["mechanisms"][0]
    group = mechanism["index_chirality"]["group_chirality_branch"]
    # Branch 0's representative is reversed, but its exact branch family can
    # repair the frame analytically; selection must not depend on choosing the
    # already-correct representative branch 1.
    assert mechanism["mapping_RP"] == identity
    assert mechanism["index_chirality"][
        "selected_analytical_branch_index"] == 0
    assert "selected_witness_index" not in group
    assert group["reversed_frame_count"] > 0
    assert mechanism["index_chirality"][
        "preserved_group_chirality_frame_count"] > 0


def test_rp_stage_preserves_preloaded_anchor_through_family_and_rmsd():
    elements = ["O", "C", "O"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[1, 2] = wbo[2, 1] = 1.0
    coords_R = np.array([[-1.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0],
                         [1.0, 0.0, 0.0]])
    coords_P = coords_R[::-1].copy()
    inputs = pipeline.step_inputs_from_arrays(
        "anchored_post_aam", elements, coords_R, wbo,
        elements, coords_P, wbo)
    config = pipeline.rp_stage_config()
    config.update({
        "anchor_map": {0: 2},
        "index_chirality": "preserve",
        "n_seeds": 1,
        "max_branches": 100,
    })
    pool = pipeline.cut_sweep(
        elements, wbo, elements, wbo,
        n_workers=0, **pipeline._rp_cut_kwargs(config))

    result = pipeline.run_rp_stage_from_pool(inputs, pool, config=config)

    assert len(result["mechanisms"]) == 1
    mapping = result["mechanisms"][0]["mapping_RP"]
    assert mapping == {0: 2, 1: 1, 2: 0}
    assert result["mechanisms"][0]["index_chirality"][
        "selected_index_chirality_violation_count"] == 0


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


def test_event_relation_uses_lossless_sparse_baseline_encoding():
    import pynauty

    elements = ["C"] * 4
    wbo_R = np.zeros((4, 4))
    wbo_P = np.zeros((4, 4))
    wbo_P[0, 1] = wbo_P[1, 0] = 0.8
    identity = {index: index for index in range(4)}

    relation, _persistent, _inverse, _fragments = _masked_relation_data(
        identity, {}, elements, wbo_R, elements, wbo_P,
        1.0, 0.2, 0.5, 0.3, None)
    graph = relation.graph("B")
    _generators, _m1, _m2, orbits, _count = pynauty.autgrp(graph)

    # Six complete atom pairs have two event colors.  Absence represents the
    # five-pair baseline exactly, leaving one exceptional relation vertex.
    assert graph.number_of_vertices == 5
    assert orbits[0] == orbits[1]
    assert orbits[2] == orbits[3]
    assert orbits[0] != orbits[2]


def test_cached_event_relation_is_identical_to_uncached_construction():
    rng = np.random.default_rng(73)
    elements = ["C", "C", "O", "H", "H", "H"]
    wbo_R = rng.uniform(0.0, 1.2, size=(6, 6))
    wbo_R = (wbo_R + wbo_R.T) / 2.0
    np.fill_diagonal(wbo_R, 0.0)
    wbo_P = wbo_R.copy()
    mapping = {index: index for index in range(6)}
    hierarchy = {"fragments": [
        {"fragment": [0, 1, 2]},
        {"fragment": [3, 4, 5]},
    ]}
    uncached = _masked_relation_data(
        mapping, hierarchy, elements, wbo_R, elements, wbo_P,
        0.2, 0.2, 0.5, 0.3, {})[0]
    context = analytical_family_static_context(
        elements, wbo_R, elements, wbo_P,
        graph_floor=0.2, dwbo_threshold=0.5,
        metal_dwbo_threshold=0.3)
    cached = _masked_relation_data(
        mapping, hierarchy, elements, wbo_R, elements, wbo_P,
        0.2, 0.2, 0.5, 0.3, {}, static_context=context)[0]

    assert cached.colors_A == uncached.colors_A
    assert cached.colors_B == uncached.colors_B
    assert cached.relation_records_A == uncached.relation_records_A
    assert cached.relation_records_B == uncached.relation_records_B

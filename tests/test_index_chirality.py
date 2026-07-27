import inspect

import numpy as np
import pytest

import rxn_core.alignment.index_chirality as index_chirality_module
import rxn_core.pipeline as pipeline_module
from rxn_core.alignment.index_chirality import (
    IndexChiralityConflict,
    mapping_event_signature,
    select_index_chirality_assignment,
)
from rxn_core.matcher import _SymBlock, _SymCand, _symmetry_state
from rxn_core.pipeline import (
    rp_stage_config,
    run_rp_stage,
    step_inputs_from_arrays,
)


def _star_wbo(atom_count=5):
    wbo = np.zeros((atom_count, atom_count), dtype=float)
    for neighbor in range(1, atom_count):
        wbo[0, neighbor] = wbo[neighbor, 0] = 1.0
    return wbo


def _tetrahedral_case():
    elements = ["C", "F", "H", "H", "H"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [-0.9, -0.9, -0.9],
        [0.9, 0.9, -0.9],
        [0.9, -0.9, 0.9],
        [-0.9, 0.9, 0.9],
    ])
    identity = {index: index for index in range(5)}
    odd = dict(identity)
    odd[3], odd[4] = odd[4], odd[3]
    return elements, coords, _star_wbo(), identity, odd


def _branch_with_candidate(source, candidate):
    return {
        "witnesses": [{
            "mapping": dict(source),
            "local_symmetry": {
                "fragments": [{
                    "symmetry": _symmetry_state(candidate),
                }],
            },
        }],
    }


def _mapping_tuple(mapping):
    return tuple(mapping[index] for index in sorted(mapping))


def _allowed_mapping_tuples(selection):
    return {
        _mapping_tuple(assignment.mapping)
        for assignment in selection.allowed_assignments
    }


def _selected_evaluation(selection):
    selected_id = selection.metadata["selected_candidate_id"]
    rows = selection.metadata["candidate_search"]["candidate_evaluations"]
    return next(row for row in rows if row["candidate_id"] == selected_id)


def _provenance_kinds(evaluation):
    return {
        str(record.get("kind", ""))
        for record in evaluation.get("provenance") or ()
    }


def test_selector_has_no_global_automorphism_enumerator():
    """The chirality option may only inspect choices retained by the AAM."""
    module_source = inspect.getsource(index_chirality_module)
    selector_source = inspect.getsource(
        index_chirality_module.select_index_chirality_assignment)
    pipeline_source = inspect.getsource(pipeline_module)

    assert "symmetry_action" not in module_source
    assert "_build_symmetry_action_family" not in module_source
    assert "isomorphisms_iter" not in module_source
    assert "itertools.permutations" not in selector_source
    assert "strictly_improving" not in module_source
    assert "while current" not in selector_source
    assert "_solve_gf2" in module_source
    selector_parameters = inspect.signature(
        index_chirality_module.select_index_chirality_assignment).parameters
    assert "volume_tolerance" not in selector_parameters
    assert "max_actions" not in selector_parameters
    assert "index_chirality_volume_tolerance" not in rp_stage_config()
    assert "index_chirality_max_actions" not in rp_stage_config()
    assert "--index-chirality-volume-tolerance" not in pipeline_source
    assert "--index-chirality-max-actions" not in pipeline_source


def test_index_chirality_off_is_a_true_pipeline_noop(monkeypatch):
    elements = ["S", "O", "O", "C", "N"]
    wbo_R = np.zeros((5, 5), dtype=float)
    wbo_P = np.zeros((5, 5), dtype=float)
    for neighbor, value in [
        (1, 1.70), (2, 1.05), (3, 0.75), (4, 0.99),
    ]:
        wbo_R[0, neighbor] = wbo_R[neighbor, 0] = value
    for neighbor, value in [
        (1, 1.61), (2, 1.61), (3, 0.75), (4, 0.90),
    ]:
        wbo_P[0, neighbor] = wbo_P[neighbor, 0] = value
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ])
    inputs = step_inputs_from_arrays(
        "index_chirality_off_noop",
        elements, coords, wbo_R,
        elements, coords, wbo_P,
    )

    def forbidden_selector(*_args, **_kwargs):
        raise AssertionError(
            "the native chirality selector ran while policy was off")

    monkeypatch.setattr(
        pipeline_module,
        "select_index_chirality_assignment",
        forbidden_selector,
    )
    config = rp_stage_config()
    config.update(index_chirality="off", n_seeds=1)
    result = run_rp_stage(inputs, config=config)

    assert result["config"]["index_chirality"] == "off"
    assert result["mechanisms"]
    assert all(
        "index_chirality" not in (mechanism.get("branch_symmetry") or {})
        for mechanism in result["mechanisms"]
    )


def test_primary_and_alternate_choose_clean_existing_ts01_style_mapping():
    elements, coords, wbo, identity, odd = _tetrahedral_case()
    candidate = _SymCand(odd).with_added_alternate(_SymCand(identity))
    branch_symmetry = _branch_with_candidate(odd, candidate)

    selection = select_index_chirality_assignment(
        odd,
        branch_symmetry,
        elements,
        coords,
        wbo,
        elements,
        coords,
        wbo,
    )

    assert selection.selected_mapping == identity
    assert _allowed_mapping_tuples(selection) == {
        _mapping_tuple(identity),
    }
    assert selection.metadata["source_index_chirality_violation_count"] > 0
    assert selection.metadata["selected_index_chirality_violation_count"] == 0
    assert selection.metadata["switchable_r_atoms"] == [3, 4]
    assert selection.metadata["defined_frame_count"] == 1
    assert selection.metadata["immutable_frame_count"] == 0
    assert (
        selection.metadata["candidate_search"][
            "unique_candidate_evaluation_count"]
        == 2
    )
    assert any(
        "alternate" in kind
        for kind in _provenance_kinds(_selected_evaluation(selection))
    )


def test_immutable_reversed_frame_is_diagnostic_and_does_not_mutate_mapping():
    """Endpoint differences outside AAM choices are reported, not repaired."""
    elements, coords_R, wbo, identity, _odd = _tetrahedral_case()
    coords_P = coords_R.copy()
    coords_P[[3, 4]] = coords_P[[4, 3]]
    fixed_candidate = _SymCand(identity)

    selection = select_index_chirality_assignment(
        identity,
        _branch_with_candidate(identity, fixed_candidate),
        elements,
        coords_R,
        wbo,
        elements,
        coords_P,
        wbo,
    )

    assert selection.selected_mapping == identity
    assert _allowed_mapping_tuples(selection) == {
        _mapping_tuple(identity),
    }
    assert selection.metadata["switchable_r_atoms"] == []
    assert selection.metadata["all_defined_frame_count"] == 1
    assert selection.metadata["defined_frame_count"] == 0
    assert selection.metadata["immutable_frame_count"] == 1
    assert selection.metadata["immutable_source_mismatch_count"] == 1
    assert len(selection.metadata["immutable_frames"]) == 1
    assert all(
        frame["reason"]
        == "no_AAM_authorized_switchable_atom_in_frame"
        and frame["source_index_chirality_mismatch"] is True
        and frame["source_mismatch_details"]
        for frame in selection.metadata["immutable_frames"]
    )
    assert selection.metadata["mapping_changes"] == []


def test_correlated_two_copy_alternate_is_never_cartesian_mixed():
    """One retained alternate row is one atomic choice, not per-index domains."""
    atom_count = 8
    elements = ["C"] * atom_count
    source = {index: index for index in range(atom_count)}
    alternate = {
        0: 3, 1: 2, 2: 1, 3: 0,
        4: 7, 5: 6, 6: 5, 7: 4,
    }
    mixed_first = {
        0: 3, 1: 2, 2: 1, 3: 0,
        4: 4, 5: 5, 6: 6, 7: 7,
    }
    mixed_second = {
        0: 0, 1: 1, 2: 2, 3: 3,
        4: 7, 5: 6, 6: 5, 7: 4,
    }
    wbo = np.zeros((atom_count, atom_count), dtype=float)
    for offset in (0, 4):
        for left in range(offset, offset + 3):
            wbo[left, left + 1] = wbo[left + 1, left] = 1.0
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [1.0, 3.0, 0.0],
        [2.0, 3.0, 0.0],
        [3.0, 3.0, 0.0],
    ])
    candidate = _SymCand(source).with_added_alternate(_SymCand(alternate))

    selection = select_index_chirality_assignment(
        source,
        _branch_with_candidate(source, candidate),
        elements,
        coords,
        wbo,
        elements,
        coords,
        wbo,
    )

    allowed = _allowed_mapping_tuples(selection)
    assert allowed == {
        _mapping_tuple(source),
        _mapping_tuple(alternate),
    }
    assert _mapping_tuple(mixed_first) not in allowed
    assert _mapping_tuple(mixed_second) not in allowed
    assert (
        selection.metadata["candidate_search"][
            "unique_candidate_evaluation_count"]
        == 2
    )


def test_exact_coplanar_frames_are_numerically_undefined_and_neutral():
    elements, _coords, wbo, identity, odd = _tetrahedral_case()
    planar_coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ])
    candidate = _SymCand(odd).with_added_alternate(_SymCand(identity))

    selection = select_index_chirality_assignment(
        odd,
        _branch_with_candidate(odd, candidate),
        elements,
        planar_coords,
        wbo,
        elements,
        planar_coords,
        wbo,
    )

    assert _allowed_mapping_tuples(selection) == {
        _mapping_tuple(odd),
        _mapping_tuple(identity),
    }
    assert selection.metadata["all_defined_frame_count"] == 0
    assert selection.metadata["defined_frame_count"] == 0
    assert selection.metadata["undefined_frame_count"] == 1
    assert all(
        frame["reason"] == "numerically_indeterminate_orientation"
        and frame["raw_determinant"] == 0.0
        and abs(frame["raw_determinant"])
        <= frame["determinant_error_bound"]
        for frame in selection.metadata["undefined_frames"]
    )


def test_closed_symmetry_block_uses_exact_parity_to_fix_orientation():
    elements, coords, wbo, _identity, odd = _tetrahedral_case()
    candidate = _SymCand(
        odd,
        blocks=(_SymBlock((2, 3, 4), (2, 3, 4)),),
    )

    selection = select_index_chirality_assignment(
        odd,
        _branch_with_candidate(odd, candidate),
        elements,
        coords,
        wbo,
        elements,
        coords,
        wbo,
    )

    changed = [
        r for r in odd
        if selection.selected_mapping[r] != odd[r]
    ]
    assert len(changed) == 2
    assert set(changed) <= {2, 3, 4}
    assert selection.metadata["selected_index_chirality_violation_count"] == 0
    # One atomic source row plus one exact GF(2) parity solution.
    assert (
        selection.metadata["candidate_search"][
            "unique_candidate_evaluation_count"]
        == 2
    )
    assert any(
        kind == "fragment_parity_solution"
        for kind in _provenance_kinds(_selected_evaluation(selection))
    )


def test_fixed_reaction_event_frames_are_diagnostic_not_blocking():
    elements, coords, wbo_R, _identity, odd = _tetrahedral_case()
    wbo_P = wbo_R.copy()
    # Still above the graph floor, but large enough to be a bond-order event.
    wbo_P[0, 1] = wbo_P[1, 0] = 0.4
    assert mapping_event_signature(
        odd, wbo_R, wbo_P, elements)[0] == ((0, 1),)

    selection = select_index_chirality_assignment(
        odd,
        _branch_with_candidate(odd, _SymCand(odd)),
        elements,
        coords,
        wbo_R,
        elements,
        coords,
        wbo_P,
    )

    assert selection.selected_mapping == odd
    assert _allowed_mapping_tuples(selection) == {
        _mapping_tuple(odd),
    }
    assert selection.metadata["defined_frame_count"] == 0
    assert selection.metadata["immutable_frame_count"] == 1
    assert selection.metadata["reaction_event_incident_frame_count"] == 1
    assert all(
        frame["reaction_event_incident"]
        for frame in selection.metadata["immutable_frames"]
    )
    assert selection.metadata["selected_index_chirality_violation_count"] == 0


def test_even_three_cycle_is_clean_under_one_tetrahedral_parity_frame():
    """An even ligand cycle must not hit the former triple-sign plateau."""
    elements = ["C", "H", "H", "H", "H"]
    coords = np.array([
        [0.0, 0.0, 0.0],
        [-0.9891213503478509, -0.3677866514678832, 1.2879252612892487],
        [0.1939744191326132, 0.9202308996398569, 0.5771037912572513],
        [-0.6364636463709805, 0.5419522204102933, -0.3165954511658161],
        [-0.32238911615896015, 0.09716731867045719, -1.5259304065189514],
    ])
    source = {0: 0, 1: 1, 2: 3, 3: 4, 4: 2}
    candidate = _SymCand(
        source,
        blocks=(_SymBlock((1, 2, 3, 4), (1, 2, 3, 4)),),
    )

    selection = select_index_chirality_assignment(
        source,
        _branch_with_candidate(source, candidate),
        elements,
        coords,
        _star_wbo(),
        elements,
        coords,
        _star_wbo(),
    )

    assert selection.selected_mapping == source
    assert selection.metadata["schema_version"] == (
        "rxn_core.index_chirality/v3")
    assert selection.metadata["all_defined_frame_count"] == 1
    assert selection.metadata["source_index_chirality_violation_count"] == 0
    assert (
        selection.metadata["candidate_search"]["gf2_solved_route_count"]
        == 1
    )


def test_anchor_removes_size_two_block_parity_without_false_conflict():
    elements, coords, wbo, _identity, odd = _tetrahedral_case()
    candidate = _SymCand(
        odd,
        blocks=(_SymBlock((3, 4), (3, 4)),),
    )

    selection = select_index_chirality_assignment(
        odd,
        _branch_with_candidate(odd, candidate),
        elements,
        coords,
        wbo,
        elements,
        coords,
        wbo,
        anchor_map={3: 4},
    )

    assert selection.selected_mapping == odd
    assert selection.metadata["switchable_r_atoms"] == []
    assert selection.metadata["defined_frame_count"] == 0
    assert selection.metadata["immutable_frame_count"] == 1
    assert selection.metadata["immutable_source_mismatch_count"] == 1
    assert any(
        record.get("reason")
        == "fewer_than_two_compatible_unanchored_atoms"
        for record in selection.metadata["candidate_search"][
            "fragment_state_diagnostics"]
    )


def _two_tetrahedral_centers():
    elements = [
        "C", "F", "H", "H", "H",
        "C", "F", "H", "H", "H",
    ]
    first = np.array([
        [0.0, 0.0, 0.0],
        [-0.9, -0.9, -0.9],
        [0.9, 0.9, -0.9],
        [0.9, -0.9, 0.9],
        [-0.9, 0.9, 0.9],
    ])
    coords = np.vstack((first, first + np.array([5.0, 0.0, 0.0])))
    wbo = np.zeros((10, 10), dtype=float)
    for center, neighbors in ((0, range(1, 5)), (5, range(6, 10))):
        for neighbor in neighbors:
            wbo[center, neighbor] = wbo[neighbor, center] = 1.0
    source = {index: index for index in range(10)}
    source[3], source[4] = source[4], source[3]
    source[8], source[9] = source[9], source[8]
    return elements, coords, wbo, source


def test_two_coexisting_size_three_blocks_are_solved_as_one_gf2_system():
    elements, coords, wbo, source = _two_tetrahedral_centers()
    candidate = _SymCand(
        source,
        blocks=(
            _SymBlock((2, 3, 4), (2, 3, 4)),
            _SymBlock((7, 8, 9), (7, 8, 9)),
        ),
    )

    selection = select_index_chirality_assignment(
        source,
        _branch_with_candidate(source, candidate),
        elements,
        coords,
        wbo,
        elements,
        coords,
        wbo,
    )

    assert selection.metadata["source_index_chirality_violation_count"] == 2
    assert selection.metadata["selected_index_chirality_violation_count"] == 0
    assert sum(
        selection.selected_mapping[r] != source[r] for r in source
    ) == 4
    search = selection.metadata["candidate_search"]
    assert search["parity_variable_count"] == 2
    assert search["gf2_equation_count"] == 2
    assert search["unique_candidate_evaluation_count"] == 2
    solved = next(
        record
        for record in search["fragment_state_diagnostics"]
        if record.get("status") == "solved"
    )
    assert solved["solution_bits"] == [1, 1]


def test_closed_blocks_from_different_fragments_are_never_composed():
    elements, coords, wbo, source = _two_tetrahedral_centers()
    state_left = _symmetry_state(_SymCand(
        source,
        blocks=(_SymBlock((2, 3, 4), (2, 3, 4)),),
    ))
    state_right = _symmetry_state(_SymCand(
        source,
        blocks=(_SymBlock((7, 8, 9), (7, 8, 9)),),
    ))
    branch_symmetry = {
        "witnesses": [{
            "mapping": source,
            "local_symmetry": {
                "fragments": [
                    {"symmetry": state_left},
                    {"symmetry": state_right},
                ],
            },
        }],
    }

    with pytest.raises(IndexChiralityConflict):
        select_index_chirality_assignment(
            source,
            branch_symmetry,
            elements,
            coords,
            wbo,
            elements,
            coords,
            wbo,
        )


def test_reaction_event_incidence_does_not_disable_parity_constraint():
    elements, coords, wbo_R, _identity, odd = _tetrahedral_case()
    wbo_P = wbo_R.copy()
    wbo_P[0, 1] = wbo_P[1, 0] = 0.4
    candidate = _SymCand(
        odd,
        blocks=(_SymBlock((2, 3, 4), (2, 3, 4)),),
    )

    selection = select_index_chirality_assignment(
        odd,
        _branch_with_candidate(odd, candidate),
        elements,
        coords,
        wbo_R,
        elements,
        coords,
        wbo_P,
    )

    assert selection.metadata["defined_frame_count"] == 1
    assert selection.metadata["frames"][0]["reaction_event_incident"] is True
    assert selection.metadata["selected_index_chirality_violation_count"] == 0
    assert selection.selected_mapping != odd


def test_nested_alternate_gets_only_its_own_fragment_parity_blocks():
    elements, coords, wbo, identity, odd = _tetrahedral_case()
    alternate_odd = dict(identity)
    alternate_odd[2], alternate_odd[3] = (
        alternate_odd[3], alternate_odd[2])
    block = _SymBlock((2, 3, 4), (2, 3, 4))
    candidate = _SymCand(
        odd,
        blocks=(block,),
    ).with_added_alternate(_SymCand(alternate_odd, blocks=(block,)))

    selection = select_index_chirality_assignment(
        odd,
        _branch_with_candidate(odd, candidate),
        elements,
        coords,
        wbo,
        elements,
        coords,
        wbo,
    )

    parity_records = [
        record
        for record in selection.metadata["candidate_search"][
            "fragment_state_diagnostics"]
        if record.get("status") == "solved"
    ]
    assert any(
        record["fragment_state"].get("kind")
        == "nested_alternate_fragment_parity_seed"
        and record["fragment_state"].get("alternate_index") == 0
        for record in parity_records
    )
    assert all(
        record["fragment_state"].get("fragment_index") == 0
        for record in parity_records
    )


def test_fragment_exact_fixed_atom_is_never_used_in_odd_representative():
    elements, coords, wbo, identity, odd = _tetrahedral_case()
    candidate = _SymCand(
        odd,
        blocks=(_SymBlock((2, 3, 4), (2, 3, 4)),),
        exact_fixed=(2,),
    )

    selection = select_index_chirality_assignment(
        odd,
        _branch_with_candidate(odd, candidate),
        elements,
        coords,
        wbo,
        elements,
        coords,
        wbo,
    )

    assert selection.selected_mapping == identity
    assert selection.selected_mapping[2] == odd[2]
    solved = next(
        record
        for record in selection.metadata["candidate_search"][
            "fragment_state_diagnostics"]
        if record.get("status") == "solved"
    )
    assert solved["exact_fixed_r_atoms"] == [2]
    assert solved["variable_blocks"][0]["canonical_odd_swap_R"] == [3, 4]


def test_nested_alternate_cannot_change_fragment_exact_fixed_atom():
    elements, coords, wbo, identity, odd = _tetrahedral_case()
    candidate = _SymCand(
        odd,
        exact_fixed=(3,),
    ).with_added_alternate(_SymCand(identity))

    selection = select_index_chirality_assignment(
        odd,
        _branch_with_candidate(odd, candidate),
        elements,
        coords,
        wbo,
        elements,
        coords,
        wbo,
    )

    # The identity alternate changes exact-fixed R3 from P4 to P3 and is
    # therefore not an authorized atomic route.
    assert all(
        assignment.mapping[3] == odd[3]
        for assignment in selection.allowed_assignments
    )
    assert all(
        provenance.get("kind") != "nested_alternate"
        for assignment in selection.allowed_assignments
        for provenance in assignment.provenance
    )

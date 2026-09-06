import numpy as np
import pytest
import rxn_core

from rxn_core.aam import search_aam
from rxn_core.analytical import compile_mechanism_families
from rxn_core.mechanisms import group_mechanisms
from rxn_core.rp import select_rp_mappings
from rxn_core.ts import analyze_transition_state
from rxn_core.domain import (
    AAMProblem,
    AAMSearchConfig,
    MolecularEndpoint,
    TransitionStateTarget,
    VibrationalModes,
)


def _endpoint(label, coordinates=None):
    coordinates = coordinates if coordinates is not None else [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
    wbo = np.array([[0.0, 1.0], [1.0, 0.0]])
    return MolecularEndpoint(
        elements=("H", "H"), coordinates=coordinates, wbo=wbo,
        label=label)


def test_endpoint_owns_immutable_validated_arrays():
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    endpoint = _endpoint("R", coordinates)
    coordinates[0, 0] = 9.0

    assert endpoint.coordinates[0, 0] == 0.0
    with pytest.raises(ValueError):
        endpoint.coordinates[0, 0] = 1.0


def test_package_root_exposes_only_typed_workflows():
    assert not hasattr(rxn_core, "run_rp_stage")
    assert not hasattr(rxn_core, "run_ts_stage")
    assert not hasattr(rxn_core, "cut_sweep")


def test_aam_problem_accepts_composition_mismatch():
    reactant = _endpoint("R")
    product = MolecularEndpoint(
        elements=("H", "C"),
        coordinates=np.zeros((2, 3)),
        wbo=np.zeros((2, 2)),
        label="P")

    problem = AAMProblem(reactant, product)
    assert not problem.balanced
    result = search_aam(problem, AAMSearchConfig(seed_count=1))
    assert result.graph.terminals
    assert max(len(result.graph.states[t].mapping) for t in result.graph.terminals) == 1


def test_partial_aam_keeps_unmatched_atoms_and_roundtrips(tmp_path):
    import json
    from rxn_core.artifacts import aam_record, aam_from_record
    source = MolecularEndpoint(('C','Cl'), np.zeros((2,3)), [[0,1],[1,0]])
    target = MolecularEndpoint(('C',), np.zeros((1,3)), [[0]])
    result = search_aam(AAMProblem(source,target), AAMSearchConfig(seed_count=1),intermediate_dir=tmp_path)
    restored = aam_from_record(aam_record(result))
    assert restored.problem.source_atom_count == 2
    assert restored.problem.target_atom_count == 1
    assert restored.graph == result.graph
    saved = aam_from_record(json.loads((tmp_path/'aam.json').read_text()))
    assert json.dumps(saved.graph.to_record(),sort_keys=True) == json.dumps(restored.graph.to_record(),sort_keys=True)
    assert not list(tmp_path.glob('*.tmp'))
    assert all(dict(result.graph.states[t].mapping) == {0:0} for t in result.graph.terminals)
    reverse = search_aam(AAMProblem(target, source), AAMSearchConfig(seed_count=1))
    assert reverse.graph.terminals
    assert any(dict(reverse.graph.states[t].mapping) == {0:0} for t in reverse.graph.terminals)


def test_search_aam_returns_complete_typed_hierarchy():
    problem = AAMProblem(_endpoint("R"), _endpoint("P"), name="h2")
    result = search_aam(
        problem,
        AAMSearchConfig(
            seed_count=1, branch_limit=100, symmetry_repair=False),
        workers=1)

    assert result.problem is problem
    assert result.graph.transitions
    assert result.branches
    grouped = group_mechanisms(result)
    assert grouped.minimum_event_mechanisms()
    mechanism = grouped.minimum_event_mechanisms()[0]
    assert mechanism.representative.degree == 2
    assert mechanism.branches
    assert mechanism.encounter_count >= 1
    assert all(branch.hierarchy.fragments for branch in mechanism.branches)
    fragment = mechanism.branches[0].hierarchy.fragments[0]
    assert fragment.representative_assignments
    assert fragment.multiplicity >= 1
    assert mechanism.branches[0].hierarchy.has_complete_exact_target_groups
    assert result.metrics.completed_group_requests >= 1
    assert result.metrics.completed_group_calculations >= 1
    assert result.metrics.retained_branch_count == len(result.graph.terminals)

    analytical = compile_mechanism_families(
        grouped, workers=1, minimum_events_only=True)
    assert analytical.mechanisms
    family = analytical.mechanisms[0].branches[0]
    assert family.family.contains(family.representative.as_dict())
    assert family.aam_branch.hierarchy.has_complete_exact_target_groups

    rp = select_rp_mappings(analytical)
    assert rp.mechanisms
    assert rp.mechanisms[0].mapping.degree == problem.atom_count
    assert rp.mechanisms[0].chirality[
        "selected_index_chirality_violation_count"] == 0


def test_nonempty_ts_processing_uses_exact_endpoint_consensus():
    elements = ("C", "N", "O", "H")
    coordinates = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0], [1.0, 1.0, 0.2],
    ])

    def bonds(entries):
        matrix = np.zeros((4, 4))
        for left, right, value in entries:
            matrix[left, right] = matrix[right, left] = value
        return matrix

    reactant = MolecularEndpoint(
        elements, coordinates,
        bonds(((0, 1, 1.0), (2, 3, 1.0))), label="R")
    product = MolecularEndpoint(
        elements, coordinates + 0.1,
        bonds(((0, 2, 1.0), (1, 3, 1.0))), label="P")
    target = MolecularEndpoint(
        elements, coordinates + 0.05,
        bonds(((0, 1, 0.5), (2, 3, 0.5),
               (0, 2, 0.5), (1, 3, 0.5))), label="TS")
    config = AAMSearchConfig(
        seed_count=1, branch_limit=100, symmetry_repair=False)
    rp = select_rp_mappings(compile_mechanism_families(
        group_mechanisms(search_aam(AAMProblem(reactant, product), config)),
        minimum_events_only=True))
    modes = np.zeros((1, 4, 3))
    modes[0, :, 0] = (1.0, -1.0, -1.0, 1.0)
    ts_target = TransitionStateTarget(
        target, VibrationalModes(np.array([-500.0]), modes))

    result = analyze_transition_state(rp, ts_target, search_config=config)

    assert len(result.mechanisms) == 1
    mechanism = result.mechanisms[0]
    assert mechanism.selected is not None
    assert mechanism.selected.sources == {"reactant", "product"}
    assert mechanism.selected.assignment.source_atoms == (
        mechanism.mechanism.core_atoms)
    assert mechanism.reactant_core_aam.branches
    assert mechanism.product_core_aam.branches


def test_no_event_mechanism_has_explicit_unscorable_ts_status():
    problem = AAMProblem(_endpoint("R"), _endpoint("P"), name="no_event")
    config = AAMSearchConfig(
        seed_count=1, branch_limit=100, symmetry_repair=False)
    rp = select_rp_mappings(compile_mechanism_families(
        group_mechanisms(search_aam(problem, config)), minimum_events_only=True))
    modes = np.zeros((1, 2, 3))
    target = TransitionStateTarget(
        _endpoint("TS"), VibrationalModes(np.array([-100.0]), modes))

    result = analyze_transition_state(rp, target, search_config=config)

    mechanism = result.mechanisms[0]
    assert mechanism.status == "no_reaction_core"
    assert mechanism.selected is None
    assert mechanism.reactant_core_aam is None
    assert mechanism.product_core_aam is None

import numpy as np
import pytest

from rxn_core.aam import search_aam
from rxn_core.analytical import compile_mapping_families
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


def test_aam_problem_rejects_composition_mismatch():
    reactant = _endpoint("R")
    product = MolecularEndpoint(
        elements=("H", "C"),
        coordinates=np.zeros((2, 3)),
        wbo=np.zeros((2, 2)),
        label="P")

    with pytest.raises(ValueError, match="compositions"):
        AAMProblem(reactant, product)


def test_search_aam_returns_complete_typed_hierarchy():
    problem = AAMProblem(_endpoint("R"), _endpoint("P"), name="h2")
    result = search_aam(
        problem,
        AAMSearchConfig(
            seed_count=1, branch_limit=100, symmetry_repair=False),
        workers=1)

    assert result.problem is problem
    assert result.mechanisms
    assert result.minimum_event_mechanisms()
    mechanism = result.minimum_event_mechanisms()[0]
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
    assert result.metrics.retained_branch_count == sum(
        len(item.branches) for item in result.mechanisms)

    analytical = compile_mapping_families(
        result, workers=1, minimum_events_only=True)
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
    rp = select_rp_mappings(compile_mapping_families(
        search_aam(AAMProblem(reactant, product), config),
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

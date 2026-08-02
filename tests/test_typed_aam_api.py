import numpy as np
import pytest

from rxn_core.aam import search_aam
from rxn_core.domain import (
    AAMProblem,
    AAMSearchConfig,
    MolecularEndpoint,
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
    assert result.metrics.retained_branch_count == sum(
        len(item.branches) for item in result.mechanisms)

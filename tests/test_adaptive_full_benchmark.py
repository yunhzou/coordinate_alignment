from types import SimpleNamespace

import numpy as np

from adaptive_full_benchmark import holdout_classes
from rxn_core import AAMProblem, AAMSearchConfig
from rxn_core.domain import MolecularEndpoint
from rxn_core.search_orientation import AAMSearchPlan


def test_holdout_comparison_normalizes_direction_and_counts_hydrogens():
    # Same heavy relation, but the hydrogen has moved from C to O: two events.
    elements = ['C', 'O', 'H']
    r = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], dtype=float)
    permutation = [1, 2, 0]
    p = np.zeros((3, 3))
    p[1, 2] = p[2, 1] = 1
    p[2, 0] = p[0, 2] = 1
    left = MolecularEndpoint(elements, np.zeros((3, 3)), r)
    right = MolecularEndpoint(['H', 'C', 'O'], np.zeros((3, 3)), p)
    problem = AAMProblem(left, right)
    raw = dict(reactant=dict(elements=elements, wbo=r),
               product=dict(elements=right.elements, wbo=p))
    observed = []
    for reverse in (False, True):
        mapping = dict(enumerate(permutation))
        if reverse:
            mapping = {p:r for r,p in mapping.items()}
        graph = SimpleNamespace(terminals=[0], states=[SimpleNamespace(mapping=tuple(mapping.items()))])
        plan = AAMSearchPlan(problem, AAMProblem(right, left) if reverse else problem,
                             AAMSearchConfig(), reverse)
        observed.append(holdout_classes(SimpleNamespace(graph=graph), plan, raw))
    assert observed[0] == observed[1]
    assert observed[0]['best_events'] == 2
    assert observed[0]['classes'][0]['counts'] == [1, 1, 0]


def test_holdout_missing_full_mapping_is_not_reported_as_zero_events():
    endpoint = MolecularEndpoint(['C'], np.zeros((1, 3)), np.zeros((1, 1)))
    problem = AAMProblem(endpoint, endpoint)
    plan = AAMSearchPlan(problem, problem, AAMSearchConfig(), False)
    raw = {side:dict(elements=['C'], wbo=[[0.]]) for side in ('reactant', 'product')}
    graph = SimpleNamespace(terminals=[0], states=[SimpleNamespace(mapping=())])
    result = holdout_classes(SimpleNamespace(graph=graph), plan, raw)
    assert result['valid_full_representatives'] == 0
    assert result['best_events'] is None

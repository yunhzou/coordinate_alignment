"""Tests of the TS comparison evaluator, not chemical accuracy evidence."""
from types import SimpleNamespace

import numpy as np
import z3

from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from score_real_ts_families import compile_slap
from publish_real_ts_comparison import scalar_events
from rxn_core.family_scoring import bond_events


def problem():
    endpoint = MolecularEndpoint(('H',) * 4, np.zeros((4, 3)), np.zeros((4, 4)))
    return AAMProblem(endpoint, endpoint)


def native(fingerprints):
    groups = {0: {'idxs': [0, 1]}, 1: {'idxs': [2, 3]}}
    graph = SimpleNamespace(labels=[1] * 4)
    return dict(lgp=(graph, graph), lap_sols={1: (
        [dict(fingerprint=np.array(fp)) for fp in fingerprints], (groups, groups))})


def test_lap_fingerprint_keeps_correlations_after_input_permutation():
    result = native([[[2, 0], [0, 2]]])
    solver, values, count = compile_slap(result, [2, 0, 3, 1], [1, 3, 0, 2], problem())
    assert count == 1 and solver.check() == z3.sat
    solver.push()
    solver.add(values[2][0] == 0)  # Final H label allows this; LAP fingerprint does not.
    assert solver.check() == z3.unsat
    solver.pop()
    solver.add(values[2][0] == 3)
    assert solver.check() == z3.sat


def test_all_returned_lap_fingerprints_are_kept():
    result = native([[[2, 0], [0, 2]], [[0, 2], [2, 0]]])
    solver, values, count = compile_slap(result, range(4), range(4), problem())
    assert count == 2
    for target in (0, 2):
        solver.push()
        solver.add(values[0][0] == target)
        assert solver.check() == z3.sat
        solver.pop()
    solver.add(values[0][0] == 0, values[1][0] == 2)
    assert solver.check() == z3.unsat  # No independently mixed fingerprints.


def test_scalar_full_atom_score_matches_event_counter():
    left = np.zeros((4, 4)); right = left.copy()
    left[0, 1] = left[1, 0] = 1.
    left[1, 2] = left[2, 1] = 1.8
    right[1, 2] = right[2, 1] = .8
    right[2, 3] = right[3, 2] = 1.
    endpoints = [MolecularEndpoint(('H',) * 4, np.zeros((4, 3)), w) for w in (left, right)]
    case = AAMProblem(*endpoints)
    expected = dict(broken=1, formed=1, order_changed=1, total=3)
    assert scalar_events(case, dict(enumerate(range(4)))) == expected
    assert bond_events(case, dict(enumerate(range(4)))) == expected

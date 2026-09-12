"""Bounded checks: no mapping-permutation or group-element enumeration."""
from collections import defaultdict
import math
from pathlib import Path
import random
import sys

import numpy as np
import pynauty
import pytest
from sympy.combinatorics import Permutation, PermutationGroup

from rxn_core import AAMProblem, MolecularEndpoint
from rxn_core.event_patterns import SignedEventIndex, _response_bounds, extract_path_events
from rxn_core.search_graph import AAMSearchGraph, SearchContext, SearchState, FragmentTransition, SearchStop

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bench'))
from reference_dense_events import DeltaPatterns


def problem(elements, r_edges, p_edges):
    def endpoint(edges):
        w = np.zeros((len(elements), len(elements)))
        for a, b, value in edges:
            w[a, b] = w[b, a] = value
        return MolecularEndpoint(tuple(elements), np.zeros((len(elements), 3)), w)
    return AAMProblem(endpoint(r_edges), endpoint(p_edges))


def raw_record(value):
    return {side: dict(elements=list(ep.elements), coordinates=ep.coordinates.tolist(), wbo=ep.wbo.tolist())
            for side, ep in (('reactant', value.reactant), ('product', value.product))}


def saved_path(value, generators=(), blocks=(), fragments=None):
    """Construct recorded path evidence directly; never run AAM search here."""
    n = value.source_atom_count
    identity = tuple((i, i) for i in range(n))
    parts = fragments or (tuple(range(n)),)
    states = [SearchState(0, 0, (), (), ())]
    transitions = []
    assigned = []
    for step, fragment in enumerate(parts):
        assigned.extend((i, i) for i in fragment)
        states.append(SearchState(step + 1, 0, tuple(sorted(assigned)), tuple((i, step) for i, _ in assigned), ()))
        symmetry = dict(witness={i: i for i in fragment}, blocks=list(blocks) if step == 0 else [],
                        exact_fixed=[], automorph_blocks=[],
                        automorph_generators=list(generators) if step == 0 else [])
        transitions.append(FragmentTransition(step, step, step + 1, fragment[0], (step, 0),
                                              dict(fragment=list(fragment), symmetry=symmetry)))
    graph = AAMSearchGraph((SearchContext(tuple(range(n)), tuple(range(n)), tuple(range(n)),
                                        iso_tolerance=1.0),), (0,), tuple(states), tuple(transitions),
                          (SearchStop(len(parts), 'objective_met'),))
    path = next(graph.paths())
    assert tuple(sorted(path.mapping.items())) == identity
    return path


@pytest.mark.parametrize('tau', [0.3, 0.5, 1.0])
def test_response_bounds_match_dense_float_comparisons(tau):
    values = np.array(sorted({x for v in (0., .1, .3, .5, .7, 1., 1.4, 1.8)
                             for x in (math.nextafter(v, -math.inf), v, math.nextafter(v, math.inf))}))
    for weight in values:
        a, b = _response_bounds(values, weight, tau)
        dense = np.where(values - weight <= -tau, -1, np.where(values - weight >= tau, 1, 0))
        assert np.array_equal(dense, [-1] * a + [0] * (b - a) + [1] * (len(values) - b))


def test_sparse_and_dense_groups_and_event_partitions_agree():
    rng = random.Random(117)
    for sample in range(8):
        elements = ('C', 'O', 'O', 'H', 'H') if sample % 2 else ('V', 'O', 'O', 'H', 'H')
        edges = [(a, b, rng.choice((0., .1, .3, .7, 1., 1.4, 1.8)))
                 for a in range(5) for b in range(a + 1, 5)]
        product = [(a, b, rng.choice((0., .1, .3, .7, 1., 1.4, 1.8))) for a, b, _ in edges]
        value = problem(elements, edges, product)
        legacy, index = DeltaPatterns(raw_record(value)), SignedEventIndex(value)
        dense_generators = [tuple(g[:5]) for g in pynauty.autgrp(legacy.base)[0]]
        dense_group = PermutationGroup([Permutation(g, size=5) for g in dense_generators] or [Permutation(list(range(5)))])
        sparse_group = PermutationGroup([Permutation(g, size=5) for g in index.source_generators] or [Permutation(list(range(5)))])
        # Generator containment proves equality; no closure/group enumeration.
        assert all(sparse_group.contains(g) for g in dense_group.generators)
        assert all(dense_group.contains(g) for g in sparse_group.generators)
        maps = [(0, 1, 2, 3, 4), (0, 2, 1, 3, 4), (0, 1, 2, 4, 3)]
        old, new = [legacy.describe(m) for m in maps], [index.describe(m) for m in maps]
        for a, b in zip(old, new):
            assert a['events'] == b['events'] and a['total'] == b['total']
        assert [[a['id'] == b['id'] for b in old] for a in old] == [[a['id'] == b['id'] for b in new] for a in new]
        assert index.base.number_of_vertices <= legacy.base.number_of_vertices


def test_equivalent_hydrogen_events_merge_but_signs_do_not():
    value = problem('CHH', [(0, 1, 1), (0, 2, 1)], [(0, 1, 1), (0, 2, .1)])
    index = SignedEventIndex(value)
    first, second = index.describe([0, 1, 2]), index.describe([0, 2, 1])
    assert first['events'] != second['events'] and first['id'] == second['id']
    value = problem('COO', [(0, 1, 1), (0, 2, 1.4)], [(0, 1, 1.4), (0, 2, 1.6)])
    assert SignedEventIndex(value).describe([0, 1, 2])['id'] != SignedEventIndex(value).describe([0, 2, 1])['id']


def test_window_retains_nonminimum_events_after_zero_event_seed():
    value = problem('COO', [(0, 1, 1), (0, 2, 1.4)], [(0, 1, 1.4), (0, 2, 1.6)])
    path = saved_path(value, [(0, 2, 1)])
    index = SignedEventIndex(value)
    result = extract_path_events(path, value, index, max_events=1)
    assert result['complete'], result
    assert {p['total'] for p in result['patterns']} == {0, 1}
    result = extract_path_events(path, value, index, max_events=0)
    assert result['complete'] and [p['total'] for p in result['patterns']] == [0]


def test_threshold_change_reuses_original_search_path():
    value = problem('COO', [(0, 1, 1), (0, 2, 1.4)], [(0, 1, 1.4), (0, 2, 1.6)])
    path = saved_path(value, [(0, 2, 1)])
    strict = extract_path_events(path, value, SignedEventIndex(value, threshold=.5))
    loose = extract_path_events(path, value, SignedEventIndex(value, threshold=.7))
    assert strict['complete'] and loose['complete']
    assert {p['total'] for p in strict['patterns']} == {0, 1}
    assert {p['total'] for p in loose['patterns']} == {0}
    assert path.context.iso_tolerance == 1.0


def test_raw_metal_events_remain_distinct_under_binary_search_constraints():
    raw = problem(('V', 'O', 'O'), [(0, 1, 2.2), (0, 2, .7)], [(0, 1, .8), (0, 2, 2.1)])
    search = problem(('V', 'O', 'O'), [(0, 1, 1), (0, 2, 1)], [(0, 1, 1), (0, 2, 1)])
    path = saved_path(search, [(0, 2, 1)])
    result = extract_path_events(path, search, SignedEventIndex(raw), max_events=2)
    assert result['complete'] and {p['total'] for p in result['patterns']} == {0, 2}
    assert np.array_equal(search.reactant.wbo, [[0, 1, 1], [1, 0, 0], [1, 0, 0]])


def test_cross_fragment_events_are_not_dropped():
    edges = [(0, 1, 1), (0, 2, 1.8), (1, 2, 1)]
    value = problem('COO', edges, edges)
    path = saved_path(value, [(0, 2, 1)], fragments=((1, 2), (0,)))
    result = extract_path_events(path, value, SignedEventIndex(value), max_events=2)
    assert result['complete'] and {p['total'] for p in result['patterns']} == {0, 2}


def test_generators_are_checked_on_all_responses_not_one_witness():
    edges = [(0, 1, 1), (0, 2, 1.4), (0, 3, 1.8)]
    value = problem('COOO', edges, edges)
    path = saved_path(value, [(0, 2, 1, 3), (0, 1, 3, 2)])
    index = SignedEventIndex(value)
    assert index.describe([0, 2, 1, 3])['total'] == index.describe([0, 1, 3, 2])['total'] == 0
    assert not index.invariant_path(path)
    result = extract_path_events(path, value, index)
    assert result['complete'] and any(p['total'] > 0 for p in result['patterns'])


def test_large_invariant_pool_never_enters_solver(monkeypatch):
    # The synthetic compressed pool represents 12! assignments. None are expanded.
    edges = [(0, i, 1) for i in range(1, 13)]
    value = problem(('C',) + ('H',) * 12, edges, edges)
    path = saved_path(value, blocks=[dict(r_atoms=list(range(1, 13)), p_atoms=list(range(1, 13)))])
    def forbidden(*args, **kwargs):
        raise AssertionError('invariant pool must remain compressed')
    monkeypatch.setattr('rxn_core.family_query.compile_path', forbidden)
    result = extract_path_events(path, value, SignedEventIndex(value))
    assert result['complete'] and result['solver_queries'] == 0 and len(result['patterns']) == 1


def test_budget_exhaustion_is_not_reported_as_complete():
    value = problem('COO', [(0, 1, 1), (0, 2, 1.4)], [(0, 1, 1.4), (0, 2, 1.6)])
    path = saved_path(value, [(0, 2, 1)])
    result = extract_path_events(path, value, SignedEventIndex(value), seconds=0)
    assert not result['complete'] and result['reason'] == 'time_budget'
    result = extract_path_events(path, value, SignedEventIndex(value), max_patterns=1)
    assert not result['complete'] and result['reason'] == 'pattern_budget'


def test_unlimited_class_extraction_checkpoints_discoveries():
    value = problem('COO', [(0, 1, 1), (0, 2, 1.4)], [(0, 1, 1.4), (0, 2, 1.6)])
    path = saved_path(value, [(0, 2, 1)])
    emitted = []
    result = extract_path_events(path, value, SignedEventIndex(value), max_patterns=None,
                                 on_pattern=lambda p: emitted.append(p['id']))
    assert result['complete']
    assert len(emitted) == len(set(emitted)) == 2
    assert set(emitted) == {p['id'] for p in result['patterns']}
    persisted = []
    class Interrupted(Exception):
        pass
    def interrupt_after_persisting(pattern):
        persisted.append(pattern)
        raise Interrupted
    with pytest.raises(Interrupted):
        extract_path_events(path, value, SignedEventIndex(value), max_patterns=None,
                            on_pattern=interrupt_after_persisting)
    assert len(persisted) == 1
    assert SignedEventIndex(value).describe(persisted[0]['mapping'])['id'] == persisted[0]['id']


def test_reject_partial_invalid_and_asymmetric_inputs():
    value = problem('COO', [], [])
    index = SignedEventIndex(value)
    assert index.counts([]).shape == (0, 2)
    for vector in ([], [0, 1], [0, 1, 1], [1, 0, 2], [0., 1., 2.]):
        with pytest.raises(ValueError):
            index.describe(vector)
    for tau in (0, -1, float('nan')):
        with pytest.raises(ValueError):
            SignedEventIndex(value, threshold=tau)
    unbalanced = AAMProblem(value.reactant, problem('CO', [], []).product)
    with pytest.raises(ValueError):
        SignedEventIndex(unbalanced)


def test_whole_event_orbit_is_excluded_in_one_query(monkeypatch):
    import rxn_core.event_patterns as module
    # Exercise the exact solver fallback separately from the faster group certificate.
    monkeypatch.setattr('rxn_core.event_certificates.EventFamilyCertificates.invariant_class',
                        lambda *args: False)
    edges = [(0, a, 1.) for a in range(1, 8)]
    value = problem(('C',) + ('H',) * 7, edges,
                    [(0, a, .4 if a == 1 else 1.) for a in range(1, 8)])
    path = saved_path(value, blocks=[dict(r_atoms=list(range(1, 8)), p_atoms=list(range(1, 8)))])
    result = extract_path_events(path, value, SignedEventIndex(value), max_patterns=1)
    assert result['complete'] and len(result['patterns']) == 1, result
    assert result['solver_queries'] == 1, result
    # Reproduce the previous prototype's repeated literal exclusions on this
    # tiny case, without enumerating the 7! complete atom assignments.
    monkeypatch.setattr(module, '_same_event_orbit', lambda t, o, p, i: module._same_event_pattern(t, o, p))
    old = extract_path_events(path, value, SignedEventIndex(value), max_patterns=8)
    assert old['complete'] and old['solver_queries'] == 7, old
    assert {p['id'] for p in old['patterns']} == {p['id'] for p in result['patterns']}


def test_correlated_actions_do_not_turn_into_independent_fragment_pools():
    edges = [(0, 1, 1.), (0, 2, 1.8), (3, 4, 1.), (3, 5, 1.8)]
    value = problem('COONOO', edges, edges)
    path = saved_path(value, [(0, 2, 1, 3, 5, 4)])
    result = extract_path_events(path, value, SignedEventIndex(value))
    assert result['complete'], result
    assert {p['total'] for p in result['patterns']} == {0, 4}
    # Swapping only one pair would produce two events, but is outside this
    # recorded two-element group and must never be invented by decoding.
    assert SignedEventIndex(value).describe([0, 2, 1, 3, 4, 5])['total'] == 2


def test_symbolic_orbit_equality_matches_canonical_identity():
    import z3
    from rxn_core.event_patterns import _same_event_orbit
    value = problem('CHHH', [(0, a, 1.) for a in (1, 2, 3)], [(0, 1, .1), (0, 2, 1.), (0, 3, 1.)])
    index = SignedEventIndex(value)
    pattern = index.describe([0, 1, 2, 3])
    for vector in ([0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 2, 1]):
        candidate = index.describe(vector)
        terms = {tuple((int(a), int(b))): 0 for a, b in zip(index.a, index.b)}
        for kind, sign in [('broken', -1), ('formed', 1)]:
            terms.update({tuple(pair): sign for pair in candidate['events'][kind]})
        solver = z3.Solver()
        solver.set(timeout=1000)
        solver.add(_same_event_orbit(terms, z3.IntVal(candidate['total']), pattern, index))
        assert (solver.check() == z3.sat) == (pattern['id'] == candidate['id'])
        # Same event location with reversed sign is a distinct output class.
        reversed_terms = {k: -v for k, v in terms.items()}
        solver = z3.Solver()
        solver.set(timeout=1000)
        solver.add(_same_event_orbit(reversed_terms, z3.IntVal(candidate['total']), pattern, index))
        assert solver.check() == z3.unsat


def test_saved_legacy_ids_and_scores_are_recomputed():
    from metal_binary_events import DeltaPatterns as ActiveIndex, recanonicalize_patterns
    value = problem('CHH', [(0, 1, 1.), (0, 2, 1.)], [(0, 1, .1), (0, 2, 1.)])
    raw = raw_record(value)
    old = DeltaPatterns(raw).describe([0, 1, 2])
    active = ActiveIndex(raw)
    saved = {old['id']: dict(old, total=99)}
    migrated = recanonicalize_patterns(active, saved)
    expected = active.describe([0, 1, 2])
    assert list(migrated) == [expected['id']]
    assert expected['id'].startswith(SignedEventIndex.schema)
    assert expected['id'] != old['id']
    assert migrated[expected['id']]['total'] == 1
    assert saved[old['id']]['total'] == 99


def test_comparison_never_claims_absence_outside_decoded_window():
    from validate_event_examples import compare_saved_slap
    slap = {'saved': dict(minimum=7, minimum_ids=['old'], patterns={'old': dict(total=7)})}
    result = compare_saved_slap(slap, {}, window=5, complete=True)
    assert result['saved']['unrecovered_status']['old'] == 'outside_event_window'
    result = compare_saved_slap(slap, {}, window=7, complete=False)
    assert result['saved']['unrecovered_status']['old'] == 'unresolved'
    result = compare_saved_slap(slap, {}, window=7, complete=True)
    assert result['saved']['unrecovered_status']['old'] == 'absent_from_this_saved_forward_graph'

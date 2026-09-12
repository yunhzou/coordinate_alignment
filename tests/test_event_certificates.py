"""Exact pruning and event-class symmetry; no group or mapping enumeration."""
import copy
from dataclasses import replace
import random

from rxn_core.event_certificates import EventFamilyCertificates
from rxn_core.event_patterns import SignedEventIndex, extract_path_events
from test_event_patterns import problem, saved_path


def test_source_symmetry_merges_events_without_entering_solver(monkeypatch):
    value = problem('CHHH', [(0,i,1) for i in (1,2,3)],
                    [(0,i,.4 if i == 1 else 1) for i in (1,2,3)])
    path = saved_path(value, blocks=[dict(r_atoms=[1,2,3], p_atoms=[1,2,3])])
    index = SignedEventIndex(value)
    assert not index.invariant_path(path)
    def forbidden(*args, **kwargs):
        raise AssertionError('source-equivalent shuffles need no solver')
    monkeypatch.setattr('rxn_core.family_query.compile_path', forbidden)
    result = extract_path_events(path, value, index, max_patterns=None)
    assert result['complete'] and len(result['patterns']) == 1
    assert result['solver_queries'] == 0 and result['reason'] == 'source_target_class_invariance'


def test_small_raw_weight_tolerance_does_not_preserve_event_count():
    value = problem('COO', [(0,1,1.), (0,2,.9)], [(0,1,.51), (0,2,.49)])
    path = saved_path(value, [(0,2,1)])
    index = SignedEventIndex(value)
    assert abs(value.product.wbo[0,1] - value.product.wbo[0,2]) < .2
    assert not index.family_certificates.invariant_class((0,1,2), index.family_certificates.actions(path))
    result = extract_path_events(path, value, index, max_patterns=None)
    assert result['complete'] and {p['total'] for p in result['patterns']} == {0,1}


def test_relaxed_bounds_preserve_all_classes_with_correlated_actions():
    rng = random.Random(732)
    for elements in ('COO', 'VOO'):
        for _ in range(4):
            r = [(0,i,rng.choice((.4,.8,1.2))) for i in (1,2)]
            p = [(0,i,rng.choice((.4,.8,1.2))) for i in (1,2)]
            value = problem(elements, r, p)
            path = saved_path(value, [(0,2,1)])
            index = SignedEventIndex(value)
            result = extract_path_events(path, value, index, max_patterns=None)
            assert result['complete']
            minimum = min(p['total'] for p in result['patterns'])
            certificate = index.family_certificates
            for window in range(3):
                assert certificate.lower_bound((0,1,2), certificate.actions(path), window) <= minimum
    value = problem('COO', [(0,1,1),(0,2,1.6)], [(0,1,1.6),(0,2,1)])
    path = saved_path(value, [(0,2,1)])
    index = SignedEventIndex(value)
    assert index.describe([0,1,2])['total'] == 2
    assert index.family_certificates.lower_bound((0,1,2), index.family_certificates.actions(path), 0) == 0


def test_source_target_symmetry_order_cannot_be_interchanged():
    value = problem('COOO', [(0,1,1),(0,2,1),(0,3,1.8)], [(0,1,1.8),(0,2,1),(0,3,1)])
    index = SignedEventIndex(value)
    source, target = (0,2,1,3), (0,1,3,2)
    safe = (('group',(source,)), ('group',(target,)))
    unsafe = tuple(reversed(safe))
    assert index.family_certificates.invariant_class((0,1,2,3), safe)
    assert not index.family_certificates.invariant_class((0,1,2,3), unsafe)
    initial = index.describe([0,1,2,3])
    assert index.describe([target[source[a]] for a in range(4)])['id'] == initial['id']
    assert index.describe([source[target[a]] for a in range(4)])['id'] != initial['id']


def test_edge_cache_respects_pool_locks_and_graph_identity():
    value = problem('CHHH', [(0,i,1) for i in (1,2,3)], [(0,i,1) for i in (1,2,3)])
    a = saved_path(value, blocks=[dict(r_atoms=[1,2,3],p_atoms=[1,2,3])])
    transition = a.graph.transitions[0]
    match = copy.deepcopy(transition.match)
    match['symmetry']['exact_fixed'] = [1]
    b = replace(a, graph=replace(a.graph, transitions=(replace(transition, match=match),)))
    certificate = EventFamilyCertificates(SignedEventIndex(value), cache_size=1)
    first = certificate.actions(a)
    assert certificate.actions(b) != first
    assert certificate.actions(a) == first


def test_factored_event_tables_equal_reference_for_all_symbolic_assignments():
    import time
    import z3
    from rxn_core.family_query import compile_path
    from rxn_core.event_patterns import _event_model
    from reference_event_model import _event_model as reference
    value = problem('COOO', [(0,1,1),(0,2,1.4),(0,3,1.8)], [(0,1,1.8),(0,2,1),(0,3,1.4)])
    path = saved_path(value, [(0,2,1,3),(0,1,3,2)])
    compiled = compile_path(path, value, {}, source_atoms=(), complete_reference=False)
    index = SignedEventIndex(value)
    old, old_total, old_metrics = reference(compiled, index, time.perf_counter()+3)
    new, new_total, new_metrics = _event_model(compiled, index, time.perf_counter()+3)
    assert new_metrics == old_metrics
    compiled.solver.add(z3.Or(old_total != new_total, *(old[k] != new[k] for k in old)))
    compiled.solver.set(timeout=1000)
    assert compiled.solver.check() == z3.unsat

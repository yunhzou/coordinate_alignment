"""Exact output comparisons for fine-grained and conditioned symmetry reuse."""
import numpy as np
import pytest

from rxn_core.alignment.branch import find_islands
from rxn_core.conditioned_symmetry import ConditionedSymmetryWorkspace
from rxn_core.cut_replay import FragmentRepair
from rxn_core.matcher import _nauty_orbits
from rxn_core.matcher.state import candidate_from_record
from rxn_core.search_symmetry import finalize_graph_symmetry, SymmetryWorkspace
from test_cut_replay import graph, raw


@pytest.mark.parametrize('tol', [.5, 1.])
@pytest.mark.parametrize('budget', [1, 2**20])
def test_extension_cache_preserves_changes_in_cuts_locks_and_boundaries(tol, budget):
    rng = np.random.default_rng(625018)
    for trial in range(4):
        elements = ['C'] * 6 + ['H'] * 4
        edges = [(i, i+1, float(rng.choice([.8, 1.2, 1.8]))) for i in range(5)]
        edges += [(0, 6, 1.), (2, 7, 1.), (4, 8, 1.), (4, 9, 1.)]
        if trial % 2:
            edges += [(0, 4, 1.)]
        source = graph(10, edges, elements)
        target = graph(10, edges[:-1] + [(1, 5, .8)], elements)
        orbits = _nauty_orbits(target, wbo_tol=tol)
        repair = FragmentRepair(source, target, orbits, cache_bytes=1, extension_cache_bytes=budget)
        for _ in range(70):
            mapping = [-1] * 10
            locked = list(map(int, rng.choice(6, 3, replace=False)))
            for a, b in zip(locked, rng.permutation([1, 3, 5])):
                mapping[a] = int(b)
            islands = [(r, 1+int(rng.integers(2))) for r in locked]
            cut = () if rng.integers(3) == 0 else (edges[int(rng.integers(len(edges)))][:2],)
            deferred = () if rng.integers(2) == 0 else (edges[int(rng.integers(len(edges)))][:2],)
            seed = int(rng.integers(10))
            options = dict(cap=int(rng.choice([1, 3, 100])), tol=tol, islands=islands, deferred=deferred)
            expected = raw(source, target, orbits, seed, mapping, cut, **options)
            assert raw(source, target, orbits, seed, mapping, cut, replay=repair, **options) == expected
            assert raw(source, target, orbits, seed, mapping, cut, replay=repair, **options) == expected
        if budget > 1:
            assert repair.stats()['extension_hits'] > 0


@pytest.mark.parametrize('tol', [.5, 1.])
def test_every_transition_and_ordered_generator_matches_reference(tol):
    source = graph(10, [(0,1),(1,2),(2,3),(3,4),(4,5),(0,6),(2,7),(4,8),(4,9)], ['C']*6+['H']*4)
    target = graph(10, [(0,1),(1,2),(2,3),(3,4),(4,5),(1,6),(3,7),(5,8),(5,9)], ['C']*6+['H']*4)
    orbits = _nauty_orbits(target, wbo_tol=tol)
    repair = FragmentRepair(source, target, orbits, extension_cache_bytes=2**20)
    workspace = ConditionedSymmetryWorkspace(target, tol, cache_bytes=2**20)
    for order in (list(range(10)), list(reversed(range(10))), [4,2,6,0,1,3,5,7,8,9]):
        for cut in ((), ((2,3),), ((0,6),), ((0,1),)):
            view = repair.for_cut(cut)
            opts = dict(cuts=cut, max_branches=100, p_orbits=orbits, iso_tol=tol)
            expected = find_islands(view.source, target, order, **opts)
            actual = find_islands(view.source, target, order, growth_replay=view, **opts)
            expected, _ = finalize_graph_symmetry(expected, target, iso_tolerance=tol)
            actual, _ = finalize_graph_symmetry(actual, target, iso_tolerance=tol, workspace=workspace)
            assert actual == expected


@pytest.mark.parametrize('budget', [1, 2**20])
def test_conditioned_generators_keep_correlated_ring_symmetry_and_sparse_indices(budget):
    import networkx as nx
    rng = np.random.default_rng(3911)
    target = graph(40, [(i,(i+1)%20) for i in range(20)] + [(i,i+20) for i in range(20)], ['C']*20+['H']*20)
    target = nx.relabel_nodes(target, {i: i*7 for i in range(40)})
    expanded = np.zeros((274, 274))
    expanded[np.ix_(np.arange(40)*7, np.arange(40)*7)] = target.graph['wbo_matrix']
    target.graph['wbo_matrix'] = expanded
    reference = SymmetryWorkspace(target, 1.)
    optimized = ConditionedSymmetryWorkspace(target, 1., cache_bytes=budget)
    for _ in range(50):
        locked = tuple((int(i), int(i*7)) for i in sorted(rng.choice(20, 4, replace=False)))
        state = dict(witness={0:0, 1:7}, blocks=[], exact_fixed=[], multiplicity=1,
                     automorph_blocks=[dict(r_atoms=[0,1], p_atoms=list(range(0,140,7)), extendable=False)])
        candidate = candidate_from_record(state)
        assert tuple(map(tuple, optimized.conditioned_generators(candidate, locked))) == reference.conditioned_generators(candidate, locked)
    if budget > 1:
        before = optimized.stats()['partition_hits']
        optimized.conditioned_generators(candidate, locked)
        assert optimized.stats()['partition_hits'] == before + 1

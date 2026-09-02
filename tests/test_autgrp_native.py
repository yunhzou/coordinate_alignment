"""Differential tests for the native ``_engine.AutGraph`` (native/src/autgrp.cpp).

``AutGraph(n, edges).generators(cells)`` must return exactly what
``pynauty.autgrp(graph)[0]`` returns for the same graph and vertex colouring:
the same generators, in the same (nauty emission) order.  The cells handed to
the native class are the pynauty ``vertex_coloring`` sets converted with
``list(cell)`` (``_pynauty_partition_cells``), which is the iteration order
nautywrap's C loop sees.

Skipped when the extension has not been built
(``.venv/bin/python native/build_engine.py``).
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

import pynauty
import pytest

from rxn_core.matcher import canonical
from rxn_core.matcher.canonical import (
    _CandidateAutomorphismCanonicalizer,
    _pynauty_partition_cells,
)

try:
    from rxn_core import _engine
except ImportError:  # pragma: no cover - depends on the build
    _engine = None

AutGraph = getattr(_engine, "AutGraph", None)

pytestmark = pytest.mark.skipif(
    AutGraph is None, reason="native AutGraph is not built")


# ---------------------------------------------------------------------------
# graph families
# ---------------------------------------------------------------------------

def _random_edges(rng, n, density):
    return [(a, b) for a in range(n) for b in range(a + 1, n)
            if rng.random() < density]


def _disjoint_copies(rng, copies, size):
    """``copies`` identical random components (large automorphism groups)."""
    component = _random_edges(rng, size, rng.uniform(0.2, 0.8))
    if not component:
        component = [(0, 1 % size)] if size > 1 else []
    return [(a + k * size, b + k * size)
            for k in range(copies) for a, b in component]


def _cycle(n):
    return [(i, (i + 1) % n) for i in range(n)]


def _complete(n):
    return [(a, b) for a in range(n) for b in range(a + 1, n)]


def _star(n):
    return [(0, i) for i in range(1, n)]


def _complete_bipartite(a, b):
    return [(i, a + j) for i in range(a) for j in range(b)]


def _graph_case(rng, index):
    """(n, edges) with a mix of random and highly symmetric graphs."""
    kind = index % 10
    if kind <= 3:
        n = rng.randint(5, 60)
        return n, _random_edges(rng, n, rng.choice([0.05, 0.15, 0.3, 0.5, 0.8]))
    if kind == 4:
        size = rng.randint(2, 6)
        copies = rng.randint(2, 60 // size)
        return size * copies, _disjoint_copies(rng, copies, size)
    if kind == 5:
        n = rng.randint(5, 60)
        return n, _cycle(n)
    if kind == 6:
        n = rng.randint(5, 30)
        return n, _complete(n)
    if kind == 7:
        n = rng.randint(5, 60)
        return n, _star(n)
    if kind == 8:
        # random graph plus isolated vertices (or entirely empty)
        core = rng.randint(0, 40)
        n = core + rng.randint(5, 20)
        return n, _random_edges(rng, core, rng.uniform(0.1, 0.6))
    a = rng.randint(2, 20)
    b = rng.randint(2, 20)
    return a + b, _complete_bipartite(a, b)


def _adjacency_dict(rng, n, edges):
    """pynauty adjacency dict; some entries list both directions / repeats,
    exactly the redundancy ``set_adjacency_dict`` tolerates."""
    adjacency = {}
    for a, b in edges:
        if rng.random() < 0.5:
            a, b = b, a
        adjacency.setdefault(a, []).append(b)
        if rng.random() < 0.3:
            adjacency.setdefault(b, []).append(a)
        if rng.random() < 0.1:
            adjacency[a].append(b)
    return adjacency


def _dict_pairs(adjacency):
    """Every (x, y) the C loop in ``_make_nygraph`` reads, in dict order."""
    return [(x, y) for x, ys in adjacency.items() for y in ys]


def _random_coloring(rng, n):
    """Random list of disjoint sets covering some or all of the vertices;
    within-cell order varies through the shuffled insertion order."""
    roll = rng.random()
    if roll < 0.08:
        return []                                    # no colouring at all
    n_cells = rng.randint(1, 6)
    vertices = list(range(n))
    rng.shuffle(vertices)
    if roll < 0.35:                                  # leave some uncoloured
        vertices = vertices[:rng.randint(0, n - 1)]
    cells = [[] for _ in range(n_cells)]
    for v in vertices:
        cells[rng.randrange(n_cells)].append(v)
    coloring = [set(cell) for cell in cells if cell or rng.random() < 0.1]
    return coloring


# ---------------------------------------------------------------------------
# (a) random coloured graphs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(4))
def test_generators_match_pynauty_on_random_graphs(seed):
    rng = random.Random(1000 + seed)
    graphs = compared = nontrivial = 0
    for index in range(150):
        n, edges = _graph_case(rng, index)
        adjacency = _adjacency_dict(rng, n, edges)
        native = AutGraph(n, _dict_pairs(adjacency))
        assert native.n_vertices == n
        pygraph = pynauty.Graph(n, directed=False, adjacency_dict=adjacency)
        graphs += 1
        # several colourings on the same (reused) objects, like production
        for _ in range(rng.randint(1, 4)):
            coloring = _random_coloring(rng, n)
            pygraph.set_vertex_coloring(coloring)
            cells = _pynauty_partition_cells(n, coloring)
            # the helper must hand over what pynauty stored, in its order
            assert cells == [list(cell) for cell in pygraph.vertex_coloring]
            expected = pynauty.autgrp(pygraph)[0]
            got = native.generators(cells)
            assert got == expected, (n, adjacency, coloring)
            assert all(type(x) is int for perm in got for x in perm)
            compared += 1
            if expected:
                nontrivial += 1
    assert graphs == 150
    assert compared >= 150
    assert nontrivial >= 100


def test_generators_match_pynauty_symmetric_families():
    """Hand-picked heavy groups and edge cases."""
    cases = []
    for n in (5, 12, 33, 60):
        cases.append((n, []))                               # all isolated
        cases.append((n, _complete(n)))
        cases.append((n, _cycle(n)))
        cases.append((n, _star(n)))
    cases.append((60, _disjoint_copies(random.Random(3), 20, 3)))
    cases.append((60, _disjoint_copies(random.Random(4), 6, 10)))
    cases.append((60, _disjoint_copies(random.Random(5), 30, 2)))
    cases.append((40, _complete_bipartite(20, 20)))
    cases.append((7, [(0, 1), (1, 2), (2, 0), (3, 4)]))   # loops-free, isolated 5,6
    rng = random.Random(77)
    for n, edges in cases:
        native = AutGraph(n, edges)
        pygraph = pynauty.Graph(
            n, directed=False,
            adjacency_dict={a: [b for x, b in edges if x == a] for a, _ in edges})
        colorings = [[], [set(range(n))], [set(range(0, n, 2))],
                     [set(range(n // 2)), set(range(n // 2, n))],
                     [{0}, {1}], _random_coloring(rng, n),
                     [set(), set(range(n))]]              # empty cell keeps 2 parts
        for coloring in colorings:
            pygraph.set_vertex_coloring(coloring)
            cells = _pynauty_partition_cells(n, coloring)
            assert cells == [list(cell) for cell in pygraph.vertex_coloring]
            assert native.generators(cells) == pynauty.autgrp(pygraph)[0], (
                n, coloring)


def test_partition_cells_follow_set_vertex_coloring_rules():
    n = 9
    # single cell covering everything -> dropped
    assert _pynauty_partition_cells(n, [set(range(n))]) == []
    # single partial cell -> completed to two cells
    cells = _pynauty_partition_cells(n, [{1, 5}])
    assert cells[0] == list({1, 5})
    assert sorted(cells[1]) == [0, 2, 3, 4, 6, 7, 8]
    assert cells == [list(c) for c in pynauty.Graph(
        n, vertex_coloring=[{1, 5}]).vertex_coloring]
    # invalid (overlapping / out of range) raises like pynauty
    with pytest.raises(ValueError):
        _pynauty_partition_cells(n, [{1, 2}, {2, 3}])
    with pytest.raises(ValueError):
        _pynauty_partition_cells(n, [{n}])


def test_autgraph_rejects_bad_input():
    graph = AutGraph(4, [(0, 1)])
    with pytest.raises(ValueError):
        graph.generators([[0, 1], [2]])            # does not cover vertex 3
    with pytest.raises(ValueError):
        graph.generators([[0, 1], [1, 2, 3]])      # repeated vertex
    with pytest.raises(ValueError):
        graph.generators([[0, 1, 2, 4]])           # out of range
    with pytest.raises(ValueError):
        AutGraph(4, [(0, 4)])
    with pytest.raises(ValueError):
        AutGraph(0, [])


def test_atom_generators_dispatch(monkeypatch):
    monkeypatch.delenv("RXN_CORE_NATIVE", raising=False)
    assert canonical._native_autgrp_available()
    monkeypatch.setenv("RXN_CORE_NATIVE", "0")
    assert not canonical._native_autgrp_available()


# ---------------------------------------------------------------------------
# (b) end to end through search_aam
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]


def _bench_cases():
    sys.path.insert(0, str(ROOT / "bench"))
    from cases import CASES
    return CASES


@pytest.mark.parametrize("case, minimum_calls", [("tempo", 500), ("tetraphenyl", 1)])
def test_search_aam_native_matches_pynauty(monkeypatch, case, minimum_calls):
    from rxn_core.aam import search_aam
    from rxn_core.domain import AAMProblem, AAMSearchConfig

    reactant, product = _bench_cases()[case]()
    problem = AAMProblem(reactant, product, name=case)
    compared = []

    def both(self, cand):
        colored = self._colored_vertices(cand, group_domains=True)
        raw_pynauty = self._raw_atom_generators_pynauty(colored)
        raw_native = self._raw_atom_generators_native(colored)
        assert raw_native == raw_pynauty          # elements and order
        expected = self._atom_generators_pynauty(cand)
        got = self._atom_generators_native(cand)
        assert got == expected
        compared.append(len(raw_pynauty))
        return got

    monkeypatch.setattr(
        _CandidateAutomorphismCanonicalizer, "atom_generators", both)
    result = search_aam(problem, AAMSearchConfig(), workers=1)
    assert result.mechanisms
    assert len(compared) >= minimum_calls
    print(f"\n{case}: compared {len(compared)} atom_generators calls "
          f"({sum(compared)} raw generators, "
          f"{sum(1 for c in compared if c)} non-trivial)")

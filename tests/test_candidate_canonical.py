import networkx as nx

from rxn_core.matcher.canonical import _CandidateAutomorphismCanonicalizer
from rxn_core.matcher.dedupe import _dedup_sym_cands
from rxn_core.matcher.orbits import _nauty_orbits
from rxn_core.matcher.state import _SymCand


def _carbon_path():
    graph = nx.Graph()
    graph.add_nodes_from((atom, {"element": "C"}) for atom in range(4))
    graph.add_edges_from(
        (left, right, {"wbo": 1.0})
        for left, right in ((0, 1), (1, 2), (2, 3)))
    return graph


def test_candidate_certificate_preserves_semantic_roles():
    canonicalizer = _CandidateAutomorphismCanonicalizer(_carbon_path())
    source = {10: 0, 11: 1}
    target = {10: 3, 11: 2}

    assert canonicalizer.certificate(source) == canonicalizer.certificate(target)


def test_candidate_certificate_distinguishes_wrong_role_assignment():
    canonicalizer = _CandidateAutomorphismCanonicalizer(_carbon_path())
    source = {10: 0, 11: 1}
    role_swapped_target = {10: 2, 11: 3}

    assert (canonicalizer.certificate(source)
            != canonicalizer.certificate(role_swapped_target))


def test_live_candidate_dedupe_never_computes_full_automorphism_group(
        monkeypatch):
    import pynauty

    def forbidden(*_args, **_kwargs):
        raise AssertionError("autgrp is forbidden in live candidate dedupe")

    monkeypatch.setattr(pynauty, "autgrp", forbidden)
    graph = _carbon_path()
    candidates = [
        _SymCand({10: 0, 11: 1}),
        _SymCand({10: 3, 11: 2}),
    ]

    deduped = _dedup_sym_cands(candidates, graph, graph)
    assert len(deduped) == 1


def test_candidate_canonicalizer_reuses_product_base(monkeypatch):
    from rxn_core.matcher import canonical

    graph = _carbon_path()
    orbits = _nauty_orbits(graph)
    original = canonical._nauty_colored_wbo_graph
    calls = []

    def counted(*args, **kwargs):
        calls.append((args, kwargs))
        return original(*args, **kwargs)

    monkeypatch.setattr(canonical, '_nauty_colored_wbo_graph', counted)
    first = _CandidateAutomorphismCanonicalizer(
        graph, p_orbits=orbits, locked_mapping={10: 0})
    second = _CandidateAutomorphismCanonicalizer(
        graph, p_orbits=orbits, locked_mapping={10: 1})

    assert len(calls) == 1
    assert first.adjacency is second.adjacency
    assert first.locked_roles != second.locked_roles


def test_candidate_canonicalizer_reuses_explicit_operation_cache(monkeypatch):
    from rxn_core.matcher import canonical

    graph = _carbon_path()
    original = canonical._nauty_colored_wbo_graph
    calls = []

    def counted(*args, **kwargs):
        calls.append((args, kwargs))
        return original(*args, **kwargs)

    monkeypatch.setattr(canonical, '_nauty_colored_wbo_graph', counted)
    cache = {}
    first = _CandidateAutomorphismCanonicalizer(
        graph, locked_mapping={10: 0}, base_cache=cache)
    second = _CandidateAutomorphismCanonicalizer(
        graph, locked_mapping={10: 1}, base_cache=cache)

    assert len(calls) == 1
    assert first.adjacency is second.adjacency

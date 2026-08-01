import networkx as nx

from rxn_core.matcher.canonical import _CandidateAutomorphismCanonicalizer
from rxn_core.matcher.dedupe import _dedup_sym_cands
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

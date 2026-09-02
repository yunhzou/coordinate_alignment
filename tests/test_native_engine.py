"""Differential tests: the compiled growth engine against the Python engine.

Every ``grow_island`` call made by a full ``search_aam`` run on the
tetraphenylmethane benchmark is executed by both engines and the results
(mappings, deferred edges, fragments, symmetry states, cap behaviour) must be
identical; the typed search results of the two modes must agree as well.
"""
import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "bench"))

from rxn_core.growth import native  # noqa: E402

pytestmark = pytest.mark.skipif(
    not native.available(), reason="native growth engine not built")


def _iso_key(isos):
    return [
        (tuple(sorted((int(k), int(v)) for k, v in dict(iso).items())),
         tuple(sorted(tuple(e) for e in iso.deferred_edges)),
         tuple(sorted(iso.fragment)),
         json.dumps(iso.symmetry, sort_keys=True, default=list))
        for iso in isos
    ]


def _search(case, workers=1):
    from cases import CASES
    from rxn_core.aam import search_aam
    from rxn_core.domain import AAMProblem, AAMSearchConfig

    R, P = CASES[case]()
    return search_aam(AAMProblem(R, P, name=case), AAMSearchConfig(),
                      workers=workers)


def _result_key(result):
    return [
        (mechanism.key,
         tuple(sorted(mechanism.representative.as_dict().items())),
         len(mechanism.branches),
         tuple(branch.hierarchy.to_record()["fragments"][0]["symmetry"].get(
             "automorph_generators") is not None
             for branch in mechanism.branches),
         tuple(tuple(sorted(branch.representative.as_dict().items()))
               for branch in mechanism.branches))
        for mechanism in result.mechanisms
    ]


def test_every_growth_call_agrees_with_python(monkeypatch):
    import rxn_core.alignment.branch as branch_mod
    import rxn_core.growth.island as island_mod
    from rxn_core.growth.result import IslandBranchLimitExceeded
    from rxn_core.matcher.policy import as_node_match_policy

    original = island_mod.grow_island
    stats = {"calls": 0, "native": 0}

    def both(g_R, g_P, seed, mapping, **kw):
        stats["calls"] += 1
        monkeypatch.setenv("RXN_CORE_NATIVE", "0")
        try:
            try:
                expected = ("ok", _iso_key(original(g_R, g_P, seed, dict(mapping), **kw)))
            except IslandBranchLimitExceeded as exc:
                expected = ("raised", (exc.count, exc.limit))
        finally:
            monkeypatch.setenv("RXN_CORE_NATIVE", "1")
        policy = as_node_match_policy(kw.get("node_policy"))
        if native.applicable(g_R, g_P, kw.get("p_orbits"), policy, kw.get("events")):
            stats["native"] += 1
            try:
                out = native.grow_island(
                    g_R, g_P, seed, dict(mapping),
                    graph_floor=kw.get("graph_floor", 0.2),
                    iso_tol=kw.get("iso_tol", 0.5),
                    min_lock_size=kw.get("min_lock_size", 1),
                    max_branches=kw.get("max_branches", 1_000_000),
                    islands_R=kw.get("islands_R"), p_orbits=kw.get("p_orbits"),
                    prior_deferred_edges=kw.get("prior_deferred_edges"),
                    allow_mapped_seed=kw.get("allow_mapped_seed", False),
                    profile=None, profile_context=None)
                got = ("ok", _iso_key(out)) if out is not None else None
            except IslandBranchLimitExceeded as exc:
                got = ("raised", (exc.count, exc.limit))
            assert got is None or got == expected, (seed, len(mapping))
        return original(g_R, g_P, seed, dict(mapping), **kw)

    monkeypatch.setattr(branch_mod, "grow_island", both)
    result = _search("tetraphenyl")
    assert result.mechanisms
    assert stats["calls"] > 300
    assert stats["native"] > 300


def test_search_results_agree_between_engines(monkeypatch):
    monkeypatch.setenv("RXN_CORE_NATIVE", "1")
    with_native = _result_key(_search("tetraphenyl"))
    monkeypatch.setenv("RXN_CORE_NATIVE", "0")
    with_python = _result_key(_search("tetraphenyl"))
    assert with_native == with_python


def test_reused_target_handles_source_elements_introduced_later():
    """Element codes are interned per process; a target built before a later
    source introduced new element symbols must not index past its table."""
    from rxn_core import _engine

    target = _engine.TargetGraph(
        ["NativeTargetC", "NativeTargetH"], [[0.0, 1.0], [1.0, 0.0]], 0.2,
        [(0, 1)], [0, 1], [(0, 1, 5)], 0)
    for index in range(8):
        source = _engine.SourceGraph(
            ["NativeTargetC", f"PreviouslyUnseenElement{index}"],
            [[0.0, 1.0], [1.0, 0.0]], 0.2, [(0, 1)])
        result = _engine.grow_island(
            source, target, 0, [-1, -1], 0.2, 0.5, 1, 100, None, [], False)
        assert result["capped"] is False
        # the unseen element has no compatible target atom: the island stops
        # at the seed with the C-X edge deferred
        assert result["isos"]
        assert all(dict(iso["mapping"]) == {0: 0} for iso in result["isos"])

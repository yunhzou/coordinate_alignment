import math

import numpy as np
import pytest

from rxn_core import classify_bonds
from rxn_core import build_graph
from rxn_core import bond_overlap_per_mode, rxn_overlap_per_mode
from rxn_core.alignment import (
    cut_sweep,
    match_wbo_graphs,
    select_min_mechanisms,
    _generate_seed_orders,
    symmetry_repair_mapping,
)
import rxn_core.alignment.branch as branch_mod
from rxn_core.matcher import (
    _SymBlock,
    _SymCand,
    _cand_map,
    _cand_possible_p_atoms,
    _dedup_sym_cands,
    _extend_sym_cands,
    _nauty_orbits,
    _support_witness_for_value,
)


def _represented_count(cand):
    count = getattr(cand, "multiplicity", 1)
    for block in getattr(cand, "blocks", ()):
        count *= math.factorial(len(block.p_atoms)) // math.factorial(
            len(block.p_atoms) - len(block.r_atoms))
    return count


def _tetramethyl_metal_graph():
    elements = ["Pd"]
    edges = []
    for i in range(4):
        c = 1 + 4 * i
        elements.append("C")
        edges.append((0, c, 0.8))
        for j in range(3):
            h = c + 1 + j
            elements.append("H")
            edges.append((c, h, 1.0))

    wbo = np.zeros((len(elements), len(elements)))
    for i, j, w in edges:
        wbo[i, j] = wbo[j, i] = w
    return elements, wbo


def test_symcand_reassigns_correlated_block_witness():
    wbo_p = np.zeros((4, 4))
    wbo_p[1, 3] = wbo_p[3, 1] = 1.0
    g_p = build_graph(["O", "O", "C", "C"], wbo_p, bond_cut=0.2)

    cand = _SymCand(
        {0: 0, 1: 1},
        (_SymBlock((0, 1), (0, 1), extendable=True),),
    )
    support = _support_witness_for_value(
        cand, 2, 3, bonded_in_frag=[0], r_wbos=[(0, 1.0)],
        g_P=g_p, iso_tol=0.1,
    )

    assert support == {0: 1}
    repaired = cand.with_witness(support)
    assert repaired.materialize()[0] == 1
    assert repaired.materialize()[1] == 0


def test_extension_collapses_correlated_orbit_duplicate_without_boundary():
    elements = ["C", "C"]
    wbo_r = np.zeros((2, 2))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    g_r = build_graph(elements, wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((4, 4))
    wbo_p[0, 2] = wbo_p[2, 0] = 1.0
    wbo_p[1, 3] = wbo_p[3, 1] = 1.0
    g_p = build_graph(["C", "C", "C", "C"], wbo_p, bond_cut=0.2)

    cand = _SymCand(
        {0: 0},
        (_SymBlock((0,), (0, 1), extendable=True),),
    )
    out = _extend_sym_cands(
        [cand], {0}, 1, g_r, g_p, {}, 0.1, None,
        p_orbits={0: 0, 1: 0, 2: 1, 3: 1},
        r_orbits={0: 0, 1: 0},
    )
    mappings = {tuple(sorted(_cand_map(c).items())) for c in out}

    assert mappings == {((0, 1), (1, 3))}


def test_extension_ignores_inactive_r_pairs():
    wbo_r = np.zeros((3, 3))
    wbo_r[1, 2] = wbo_r[2, 1] = 1.0
    wbo_r[0, 2] = wbo_r[2, 0] = 0.4
    g_r = build_graph(["C", "C", "O"], wbo_r, bond_cut=0.5)

    wbo_p = np.zeros((3, 3))
    wbo_p[1, 2] = wbo_p[2, 1] = 1.0
    g_p = build_graph(["C", "C", "O"], wbo_p, bond_cut=0.5)

    out = _extend_sym_cands(
        [_SymCand({0: 0, 1: 1})], {0, 1}, 2,
        g_r, g_p, {}, 0.1, None,
    )

    assert len(out) == 1


def test_extension_rejects_active_r_pair_mismatch():
    wbo_r = np.zeros((3, 3))
    wbo_r[1, 2] = wbo_r[2, 1] = 1.0
    wbo_r[0, 2] = wbo_r[2, 0] = 0.6
    g_r = build_graph(["C", "C", "O"], wbo_r, bond_cut=0.5)

    wbo_p = np.zeros((3, 3))
    wbo_p[1, 2] = wbo_p[2, 1] = 1.0
    g_p = build_graph(["C", "C", "O"], wbo_p, bond_cut=0.5)

    out = _extend_sym_cands(
        [_SymCand({0: 0, 1: 1})], {0, 1}, 2,
        g_r, g_p, {}, 0.1, None,
    )

    assert out == []


def test_unique_center_can_grow_to_many_compressed_targets():
    wbo_r = np.zeros((2, 2))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    g_r = build_graph(["C", "H"], wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((3, 3))
    wbo_p[0, 1] = wbo_p[1, 0] = 1.0
    wbo_p[0, 2] = wbo_p[2, 0] = 1.0
    g_p = build_graph(["C", "H", "H"], wbo_p, bond_cut=0.2)

    out = _extend_sym_cands(
        [_SymCand({0: 0})], {0}, 1, g_r, g_p, {}, 0.1, None,
        p_orbits={0: 0, 1: 1, 2: 1},
        r_orbits={0: 0, 1: 1},
    )

    assert len(out) == 1
    assert isinstance(out[0], _SymCand)
    assert out[0].blocks[0].r_atoms == (1,)
    assert out[0].blocks[0].p_atoms == (1, 2)


def test_planar_metal_tetramethyl_metal_carbon_branch_is_fourfold():
    elements, wbo = _tetramethyl_metal_graph()
    g_r = build_graph(elements, wbo, bond_cut=0.2)
    g_p = build_graph(elements, wbo, bond_cut=0.2)
    r_orbits = _nauty_orbits(g_r, wbo_tol=0.2)
    p_orbits = _nauty_orbits(g_p, wbo_tol=0.2)

    out = _extend_sym_cands(
        [_SymCand({0: 0})], {0}, 1,
        g_r, g_p, {}, 0.5, None,
        p_orbits=p_orbits, r_orbits=r_orbits,
        anchor_u=0, anchor_wbo=0.8,
    )

    assert len(out) == 1
    assert _represented_count(out[0]) == 4
    assert out[0].blocks[0].r_atoms == (1,)
    assert out[0].blocks[0].p_atoms == (1, 5, 9, 13)


def test_nauty_orbits_group_near_wbo_by_tolerance():
    pytest.importorskip("pynauty")
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[0, 2] = wbo[2, 0] = 1.1
    g = build_graph(["C", "H", "H"], wbo, bond_cut=0.2)

    loose = _nauty_orbits(g, wbo_tol=0.2)
    strict = _nauty_orbits(g, wbo_tol=0.05)

    assert loose[1] == loose[2]
    assert loose[0] != loose[1]
    assert strict[1] != strict[2]


def test_nauty_orbits_drop_into_single_step_extension():
    pytest.importorskip("pynauty")
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    wbo_r[0, 2] = wbo_r[2, 0] = 1.1
    g_r = build_graph(["C", "H", "H"], wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((3, 3))
    wbo_p[0, 1] = wbo_p[1, 0] = 1.0
    wbo_p[0, 2] = wbo_p[2, 0] = 1.1
    g_p = build_graph(["C", "H", "H"], wbo_p, bond_cut=0.2)
    r_orbits = _nauty_orbits(g_r, wbo_tol=0.2)
    p_orbits = _nauty_orbits(g_p, wbo_tol=0.2)

    hydrogen = _extend_sym_cands(
        [_SymCand({0: 0})], {0}, 1,
        g_r, g_p, {}, 0.2, None,
        p_orbits=p_orbits, r_orbits=r_orbits,
        anchor_u=0, anchor_wbo=1.0,
    )

    assert len(hydrogen) == 1
    assert hydrogen[0].blocks[0].r_atoms == (1,)
    assert hydrogen[0].blocks[0].p_atoms == (1, 2)


def test_planar_metal_tetramethyl_child_h_keeps_parent_degeneracy():
    elements, wbo = _tetramethyl_metal_graph()
    g_r = build_graph(elements, wbo, bond_cut=0.2)
    g_p = build_graph(elements, wbo, bond_cut=0.2)
    r_orbits = _nauty_orbits(g_r, wbo_tol=0.2)
    p_orbits = _nauty_orbits(g_p, wbo_tol=0.2)

    carbon = _extend_sym_cands(
        [_SymCand({0: 0})], {0}, 1,
        g_r, g_p, {}, 0.5, None,
        p_orbits=p_orbits, r_orbits=r_orbits,
        anchor_u=0, anchor_wbo=0.8,
    )
    hydrogen = _extend_sym_cands(
        carbon, {0, 1}, 2,
        g_r, g_p, {}, 0.5, None,
        p_orbits=p_orbits, r_orbits=r_orbits,
        anchor_u=1, anchor_wbo=1.0,
    )

    assert len(hydrogen) == 1
    assert hydrogen[0].multiplicity == 12
    assert _represented_count(hydrogen[0]) == 12


def test_planar_metal_tetramethyl_two_metal_carbon_edges_are_correlated():
    elements, wbo = _tetramethyl_metal_graph()
    g_r = build_graph(elements, wbo, bond_cut=0.2)
    g_p = build_graph(elements, wbo, bond_cut=0.2)
    r_orbits = _nauty_orbits(g_r, wbo_tol=0.2)
    p_orbits = _nauty_orbits(g_p, wbo_tol=0.2)

    first = _extend_sym_cands(
        [_SymCand({0: 0})], {0}, 1,
        g_r, g_p, {}, 0.5, None,
        p_orbits=p_orbits, r_orbits=r_orbits,
        anchor_u=0, anchor_wbo=0.8,
    )
    second = _extend_sym_cands(
        first, {0, 1}, 5,
        g_r, g_p, {}, 0.5, None,
        p_orbits=p_orbits, r_orbits=r_orbits,
        anchor_u=0, anchor_wbo=0.8,
    )

    assert len(second) == 1
    assert _represented_count(second[0]) == 12
    assert second[0].blocks[0].r_atoms == (1, 5)
    assert second[0].blocks[0].p_atoms == (1, 5, 9, 13)


def test_hidden_alternate_witness_can_extend_later_frontier_atom():
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 2] = wbo_r[2, 0] = 1.0
    g_r = build_graph(["C", "C", "O"], wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((3, 3))
    wbo_p[1, 2] = wbo_p[2, 1] = 1.0
    g_p = build_graph(["C", "C", "O"], wbo_p, bond_cut=0.2)

    compressed = _SymCand(
        {0: 0, 1: 1},
        multiplicity=2,
        alternates=[({0: 1, 1: 0}, 1)],
    )
    out = _extend_sym_cands(
        [compressed], {0, 1}, 2,
        g_r, g_p, {}, 0.5, None,
        anchor_u=0, anchor_wbo=1.0,
    )

    assert len(out) == 1
    assert _cand_map(out[0]) == {0: 1, 1: 0, 2: 2}


def test_weighted_tolerance_prevents_free_terminal_match_at_default_iso_tol():
    wbo_r = np.zeros((2, 2))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    g_r = build_graph(["C", "H"], wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((4, 4))
    wbo_p[0, 1] = wbo_p[1, 0] = 1.0
    wbo_p[2, 3] = wbo_p[3, 2] = 1.0
    g_p = build_graph(["C", "H", "C", "H"], wbo_p, bond_cut=0.2)

    out = _extend_sym_cands(
        [_SymCand({0: 0})], {0}, 1,
        g_r, g_p, {}, 0.5, None,
        anchor_u=0, anchor_wbo=1.0,
    )

    assert len(out) == 1
    assert _cand_possible_p_atoms(out[0]) == {0, 1}


def test_active_r_edge_requires_active_target_edge_even_with_loose_iso_tol():
    wbo_r = np.zeros((2, 2))
    wbo_r[0, 1] = wbo_r[1, 0] = 0.964
    g_r = build_graph(["C", "H"], wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((3, 3))
    wbo_p[0, 1] = wbo_p[1, 0] = 0.963
    g_p = build_graph(["C", "H", "H"], wbo_p, bond_cut=0.2)

    out = _extend_sym_cands(
        [_SymCand({0: 0})], {0}, 1,
        g_r, g_p, {}, 1.0, None,
        anchor_u=0, anchor_wbo=0.964,
    )

    assert len(out) == 1
    assert all(2 not in _cand_possible_p_atoms(c) for c in out)


def test_deferred_boundary_prevents_false_orbit_dedup():
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    wbo_r[0, 2] = wbo_r[2, 0] = 0.4
    g_r = build_graph(["C", "C", "O"], wbo_r, bond_cut=0.5)

    wbo_p = np.zeros((3, 3))
    wbo_p[0, 1] = wbo_p[1, 0] = 1.0
    wbo_p[0, 2] = wbo_p[2, 0] = 0.4
    g_p = build_graph(["C", "C", "O"], wbo_p, bond_cut=0.5)

    cands = [_SymCand({0: 0, 1: 1}), _SymCand({0: 1, 1: 0})]
    no_boundary = _dedup_sym_cands(
        cands, g_r, g_p,
        r_orbits={0: 0, 1: 0, 2: 1},
        p_orbits={0: 0, 1: 0, 2: 1},
        fragment={0, 1},
    )
    with_boundary = _dedup_sym_cands(
        cands, g_r, g_p,
        r_orbits={0: 0, 1: 0, 2: 1},
        p_orbits={0: 0, 1: 0, 2: 1},
        fragment={0, 1},
        deferred_edges={(0, 2)},
    )

    assert len(no_boundary) == 1
    assert len(with_boundary) == 2


def test_symmetry_repair_removes_false_orbit_swap_changes():
    elements = ["C", "O", "C", "O"]
    wbo = np.zeros((4, 4))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[2, 3] = wbo[3, 2] = 1.0
    g_r = build_graph(elements, wbo, bond_cut=0.2)
    g_p = build_graph(elements, wbo, bond_cut=0.2)

    bad = {0: 0, 1: 3, 2: 2, 3: 1}
    before = classify_bonds(bad, wbo, wbo)
    repaired, stats = symmetry_repair_mapping(
        bad, wbo, wbo, g_r, g_p,
        p_orbits={0: 0, 1: 1, 2: 2, 3: 1},
        min_changes=1,
        return_stats=True,
    )
    after = classify_bonds(repaired, wbo, wbo)

    assert (len(before[0]), len(before[1])) == (2, 2)
    assert (len(after[0]), len(after[1])) == (0, 0)
    assert stats["repaired"] is True
    assert stats["evaluated"] <= 4


def test_classify_bonds_uses_lower_metal_event_threshold():
    wbo_r = np.zeros((2, 2))
    wbo_p = np.zeros((2, 2))
    wbo_p[0, 1] = wbo_p[1, 0] = 0.35
    mapping = {0: 0, 1: 1}

    organic = classify_bonds(
        mapping, wbo_r, wbo_p,
        dwbo_threshold=0.5,
        elements_R=["C", "O"], elements_P=["C", "O"],
        metal_dwbo_threshold=0.3)
    metal = classify_bonds(
        mapping, wbo_r, wbo_p,
        dwbo_threshold=0.5,
        elements_R=["Pd", "O"], elements_P=["Pd", "O"],
        metal_dwbo_threshold=0.3)

    assert organic[1] == []
    assert metal[1] == [(0, 1, 0.0, 0.35)]


def test_generate_seed_orders_honors_trial_cap():
    elements = ["C", "C", "C", "C", "H", "H"]
    wbo = np.zeros((6, 6))
    for i, j in [(0, 1), (1, 2), (2, 3), (0, 4), (3, 5)]:
        wbo[i, j] = wbo[j, i] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)

    orders = _generate_seed_orders(g, n_trials=2)

    assert len(orders) == 2
    assert all(g.nodes[order[0]]["element"] != "H" for order in orders)


def test_match_wbo_graphs_uses_three_seed_contract():
    elements = ["C", "C", "H"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[1, 2] = wbo[2, 1] = 1.0

    result = match_wbo_graphs(
        elements, wbo, elements, wbo,
        n_seeds=3, graph_floor=0.2, iso_tol=0.5,
        repair_symmetry=False,
    )

    assert result.n_seeds == 3
    assert result.best is not None
    assert result.best.score == (0, 0)
    assert (len(result.best.broken), len(result.best.formed)) == (0, 0)


def test_find_islands_reuses_precomputed_orbits(monkeypatch):
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 1.0
    g_r = build_graph(["C", "H"], wbo, bond_cut=0.2)
    g_p = build_graph(["C", "H"], wbo, bond_cut=0.2)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("orbits should be supplied by caller")

    monkeypatch.setattr(branch_mod, "_nauty_orbits", fail_if_called)
    branches = branch_mod.find_islands(
        g_r, g_p, [0, 1],
        p_orbits={0: 0, 1: 1},
        r_orbits={0: 0, 1: 1},
    )

    assert branches
    assert branches[0].mapping == {0: 0, 1: 1}


def test_partial_witness_mode_scores_use_full_mode_norm():
    modes_R = np.zeros((1, 2, 3))
    modes_R[0, 0, 0] = 1.0
    full_mode_norms = np.array([math.sqrt(2.0)])

    V = np.zeros((2, 3))
    V[0, 0] = 1.0
    beta = bond_overlap_per_mode(modes_R, V, mode_norms=full_mode_norms)

    delta = np.zeros((2, 3))
    delta[0, 0] = 1.0
    rho = rxn_overlap_per_mode(modes_R, delta, [0],
                               mode_norms=full_mode_norms)

    assert np.allclose(beta, [1.0 / math.sqrt(2.0)])
    assert np.allclose(rho, [1.0 / math.sqrt(2.0)])


def test_cut_sweep_is_core_mechanism_discovery_api():
    elements = ["C", "C", "H"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[1, 2] = wbo[2, 1] = 1.0

    pool = cut_sweep(
        elements, wbo, elements, wbo,
        n_workers=0, n_seeds=1, max_branches=100,
        cut_floor=0.2, symmetry_repair=False,
    )
    minimal = select_min_mechanisms(pool)

    assert minimal
    assert all(len(sig[0]) + len(sig[1]) == 0 for sig in minimal)

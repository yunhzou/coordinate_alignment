import numpy as np

from rxn_core import classify_bonds
from rxn_core.pq import (
    _SymBlock,
    _SymCand,
    _cand_map,
    _cand_possible_p_atoms,
    _dedup_sym_cands,
    _extend_sym_cands,
    _generate_seed_orders,
    _support_witness_for_value,
    build_graph,
    symmetry_repair_mapping,
)


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
        [cand], {0}, 1, g_r, g_p, {}, {}, 0.1, None,
        p_orbits={0: 0, 1: 0, 2: 1, 3: 1},
        r_orbits={0: 0, 1: 0},
    )
    mappings = {tuple(sorted(_cand_map(c).items())) for c in out}

    assert mappings == {((0, 1), (1, 3))}


def test_extension_checks_full_wbo_vector_not_sparse_edges_only():
    wbo_r = np.zeros((3, 3))
    wbo_r[1, 2] = wbo_r[2, 1] = 1.0
    wbo_r[0, 2] = wbo_r[2, 0] = 0.4
    g_r = build_graph(["C", "C", "O"], wbo_r, bond_cut=0.5)

    wbo_p = np.zeros((3, 3))
    wbo_p[1, 2] = wbo_p[2, 1] = 1.0
    g_p = build_graph(["C", "C", "O"], wbo_p, bond_cut=0.5)

    out = _extend_sym_cands(
        [_SymCand({0: 0, 1: 1})], {0, 1}, 2,
        g_r, g_p, {}, {}, 0.1, None,
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
        [_SymCand({0: 0})], {0}, 1, g_r, g_p, {}, {}, 0.1, None,
        p_orbits={0: 0, 1: 1, 2: 1},
        r_orbits={0: 0, 1: 1},
    )

    assert len(out) == 1
    assert isinstance(out[0], _SymCand)
    assert out[0].blocks[0].r_atoms == (1,)
    assert out[0].blocks[0].p_atoms == (1, 2)


def test_growth_anchor_prevents_free_terminal_match_under_loose_iso_tol():
    wbo_r = np.zeros((2, 2))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    g_r = build_graph(["C", "H"], wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((4, 4))
    wbo_p[0, 1] = wbo_p[1, 0] = 1.0
    wbo_p[2, 3] = wbo_p[3, 2] = 1.0
    g_p = build_graph(["C", "H", "C", "H"], wbo_p, bond_cut=0.2)

    out = _extend_sym_cands(
        [_SymCand({0: 0})], {0}, 1,
        g_r, g_p, {}, {}, 1.0, None,
        anchor_u=0, anchor_wbo=1.0,
    )

    assert len(out) == 1
    assert _cand_possible_p_atoms(out[0]) == {0, 1}


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


def test_generate_seed_orders_honors_trial_cap():
    elements = ["C", "C", "C", "C", "H", "H"]
    wbo = np.zeros((6, 6))
    for i, j in [(0, 1), (1, 2), (2, 3), (0, 4), (3, 5)]:
        wbo[i, j] = wbo[j, i] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)

    orders = _generate_seed_orders(g, n_trials=2)

    assert len(orders) == 2
    assert all(g.nodes[order[0]]["element"] != "H" for order in orders)

import numpy as np

from rxn_core import classify_bonds
from rxn_core.pq import (
    _SymBlock,
    _SymCand,
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

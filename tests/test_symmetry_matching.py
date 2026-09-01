import json
import math

import numpy as np
import pytest

from rxn_core import classify_bonds
from rxn_core import build_graph
from rxn_core import bond_overlap_per_mode, rxn_overlap_per_mode
from rxn_core.alignment import (
    cut_sweep,
    match_wbo_graphs,
    run_cut_sweep_chunk,
    select_min_mechanisms,
    _generate_seed_orders,
    symmetry_repair_mapping,
)
from rxn_core.alignment.sweep import (
    attach_completed_candidate_groups,
    _branch_symmetry_record,
    _color_groups_from_blocks,
    _core_mapping_variants,
    _mechanism_signature,
    _orbit_bond_key,
    _pool_add,
    complete_chosen_automorphism_groups,
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
    _symmetry_state,
)
from rxn_core.matcher.orbits import _wbo_tolerance_bucket_lookup
from rxn_core.growth.result import _IsoResult
from rxn_core.growth import IslandBranchLimitExceeded, grow_island


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


def test_fragment_branch_cap_is_loud_and_allows_exact_limit():
    g_r = build_graph(["C"], np.zeros((1, 1)), bond_cut=0.2)
    wbo_two = np.zeros((4, 4))
    wbo_two[0, 2] = wbo_two[2, 0] = 1.0
    wbo_two[1, 3] = wbo_two[3, 1] = 1.0
    g_p_two = build_graph(["C", "C", "O", "N"], wbo_two, bond_cut=0.2)
    wbo_three = np.zeros((6, 6))
    for carbon, ligand in ((0, 3), (1, 4), (2, 5)):
        wbo_three[carbon, ligand] = wbo_three[ligand, carbon] = 1.0
    g_p_three = build_graph(
        ["C", "C", "C", "O", "N", "F"], wbo_three, bond_cut=0.2)

    assert len(grow_island(g_r, g_p_two, 0, {}, max_branches=2)) == 2
    with pytest.raises(IslandBranchLimitExceeded) as exc_info:
        grow_island(g_r, g_p_three, 0, {}, max_branches=2)
    assert exc_info.value.count == 3
    assert exc_info.value.limit == 2


def test_incremental_extension_enforces_cap_after_symmetry_deduplication():
    g_r = build_graph(
        ["Si", "N"], np.array([[0.0, 1.0], [1.0, 0.0]]), bond_cut=0.2)
    wbo = np.zeros((9, 9))
    for nitrogen, marker in ((1, 5), (2, 6), (3, 7), (4, 8)):
        wbo[0, nitrogen] = wbo[nitrogen, 0] = 1.0
        wbo[nitrogen, marker] = wbo[marker, nitrogen] = 1.0
    g_p = build_graph(
        ["Si", "N", "N", "N", "N", "O", "F", "Cl", "Br"],
        wbo,
        bond_cut=0.2,
    )

    with pytest.raises(IslandBranchLimitExceeded) as exc_info:
        grow_island(
            g_r,
            g_p,
            0,
            {},
            max_branches=2,
            r_orbits=_nauty_orbits(g_r),
            p_orbits=_nauty_orbits(g_p),
        )

    assert exc_info.value.count == 4
    assert exc_info.value.limit == 2


def test_cut_sweep_compact_metrics_are_opt_in_and_respect_branch_cap():
    elements = ["C", "C"]
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 1.0

    ordinary = cut_sweep(
        elements, wbo, elements, wbo,
        n_workers=0, n_seeds=1, max_branches=2,
        symmetry_repair=False)
    measured, metrics = cut_sweep(
        elements, wbo, elements, wbo,
        n_workers=0, n_seeds=1, max_branches=2,
        symmetry_repair=False, return_metrics=True)

    assert measured == ordinary
    assert metrics["configured_max_branches"] == 2
    assert metrics["cuts"] == 2
    assert metrics["seed_orders"] == 2
    assert metrics["growth_calls"] > 0
    assert metrics["max_live_branches"] <= 2
    assert metrics["max_growth_candidates"] <= 2

    parallel, parallel_metrics = cut_sweep(
        elements, wbo, elements, wbo,
        n_workers=2, n_seeds=1, max_branches=2,
        symmetry_repair=False, return_metrics=True)
    # Multiprocessing may return a different symmetry-equivalent witness;
    # mechanism-family keys and compact metrics are the invariant objects.
    assert set(parallel) == set(ordinary)
    assert parallel_metrics["cuts"] == metrics["cuts"]
    assert parallel_metrics["seed_orders"] == metrics["seed_orders"]
    assert parallel_metrics["max_live_branches"] <= 2


def test_heavy_only_cut_items_keep_hydrogens_but_do_not_cut_xh_edges():
    from rxn_core.alignment.sweep import cut_sweep_items

    elements = ["C", "C", "H", "H"]
    wbo = np.zeros((4, 4))
    for left, right in ((0, 1), (0, 2), (1, 3)):
        wbo[left, right] = wbo[right, left] = 1.0

    assert cut_sweep_items(
        wbo, elements=elements, heavy_only=True) == [(), ((0, 1),)]


def test_worker_local_cut_compression_is_exactly_serial_equivalent():
    elements = ["C", "O", "N"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[1, 2] = wbo[2, 1] = 1.0
    kwargs = {
        "n_seeds": 2,
        "max_branches": 100,
        "symmetry_repair": False,
    }

    serial = cut_sweep(
        elements, wbo, elements, wbo, n_workers=0, **kwargs)
    compressed_parallel = cut_sweep(
        elements, wbo, elements, wbo, n_workers=2, **kwargs)

    assert compressed_parallel == serial


def test_parallel_cut_sweep_persists_disjoint_reduction_buckets(tmp_path):
    elements = ["C", "O", "N"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[1, 2] = wbo[2, 1] = 1.0
    kwargs = {
        "n_seeds": 2,
        "max_branches": 100,
        "symmetry_repair": False,
    }

    serial = cut_sweep(
        elements, wbo, elements, wbo, n_workers=0, **kwargs)
    parallel = cut_sweep(
        elements, wbo, elements, wbo, n_workers=2,
        intermediate_dir=tmp_path, **kwargs)

    assert parallel == serial
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["mechanism_count"] == len(serial)
    assert manifest["worker_pool_entries"] >= len(serial)
    assert list(tmp_path.glob("raw_bucket_*.pkl"))
    assert list(tmp_path.glob("reduced_bucket_*.pkl"))


def test_one_cut_chunk_parallelizes_across_seed_orders(monkeypatch):
    import rxn_core.alignment.sweep as sweep_module

    observed = {}

    def fake_parallel(_el_r, _wbo_r, _el_t, _wbo_t, cfg, workers,
                      _core_r, cuts, **_kwargs):
        observed["workers"] = workers
        observed["seeds"] = cfg["n_seeds"]
        observed["cuts"] = cuts
        return {}

    monkeypatch.setattr(
        sweep_module, "_cut_sweep_chunk_parallel", fake_parallel)
    elements = ["C", "C"]
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 1.0

    result = run_cut_sweep_chunk(
        elements, wbo, elements, wbo, [((0, 1),)],
        n_workers=10, n_seeds=3, max_branches=10,
        symmetry_repair=False)

    assert result == {}
    assert observed == {
        "workers": 3,
        "seeds": 3,
        "cuts": [((0, 1),)],
    }


def test_parallel_cut_sweep_consumes_results_in_cut_seed_order(monkeypatch):
    import rxn_core.alignment.sweep as sweep_module

    consumed_work = []

    class OrderedPool:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def imap(self, _worker, work, chunksize=1):
            assert chunksize == 1
            consumed_work.extend(work)
            return iter({'pool': {}} for _item in work)

        def imap_unordered(self, *_args, **_kwargs):
            raise AssertionError(
                'completion-order merging changes retained representatives')

    monkeypatch.setattr(sweep_module.mp, 'Pool', OrderedPool)
    elements = ['C', 'O']
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 1.0
    cfg = {
        'cut_floor': 0.2,
        'max_branches': 100,
        'n_seeds': 3,
        'chunksize': 1,
    }

    result = sweep_module._cut_sweep_parallel(
        elements, wbo, elements, wbo, cfg, 2, (),
        collect_metrics=False)

    assert result == {}
    assert [seed for _cut, seed, *_rest in consumed_work] == [
        0, 1, 2, 0, 1, 2,
    ]


def test_live_branch_cap_discards_only_overflowing_parent_subtree(monkeypatch):
    g_r = build_graph(["C", "O"], np.zeros((2, 2)), bond_cut=0.2)
    g_p = build_graph(
        ["C", "C", "O", "O", "O"], np.zeros((5, 5)), bond_cut=0.2)

    def fake_grow(_g_r, _g_p, seed, mapping, **_kwargs):
        if seed == 0:
            return [
                _IsoResult({0: 0}, fragment={0}),
                _IsoResult({0: 1}, fragment={0}),
            ]
        if mapping[0] == 0:
            # This parent alone would create three children under a cap of
            # two, so its subtree must be removed atomically.
            return [
                _IsoResult({1: p}, fragment={1}) for p in (2, 3, 4)
            ]
        return [_IsoResult({1: 2}, fragment={1})]

    monkeypatch.setattr(branch_mod, "grow_island", fake_grow)
    monkeypatch.setattr(
        branch_mod, "_chemistry_orbit_signature",
        lambda mapping, *_args, **_kwargs: tuple(sorted(mapping.items())))

    branches = branch_mod.find_islands(
        g_r, g_p, [0, 1], orbit_dedup=False, max_branches=2)

    assert [branch.mapping for branch in branches] == [{0: 1, 1: 2}]


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


def test_mapping_variation_blocks_capture_branch_dedupe_pool():
    blocks = branch_mod._mapping_variation_blocks([
        {0: 11, 4: 9, 8: 0},
        {0: 9, 4: 11, 8: 0},
    ])

    assert blocks == [{
        "source": "interbranch",
        "r_atoms": [0, 4],
        "p_atoms": [9, 11],
        "extendable": False,
        "open": False,
        "assignments": 2,
    }]


def test_branch_symmetry_record_closes_open_pool_with_mapping_owner():
    branch = branch_mod._Branch()
    branch.commit(
        _IsoResult(
            {27: 27},
            fragment={27},
            symmetry={
                "witness": {27: 27},
                "blocks": [{
                    "r_atoms": [27],
                    "p_atoms": [26, 27],
                    "extendable": False,
                    "open": True,
                    "assignments": "2",
                }],
            },
        ),
        build_graph(["H"] * 28, np.zeros((28, 28)), bond_cut=0.2),
    )
    branch.mapping[22] = 26
    branch.islands_R[22] = 2
    branch.islands_P[26] = 2

    record = _branch_symmetry_record(branch)

    assert record["blocks"][0]["r_atoms"] == [22, 27]
    assert record["blocks"][0]["p_atoms"] == [26, 27]
    assert record["blocks"][0]["open"] is False
    assert record["blocks"][0]["assignments"] == "2!"


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


def test_core_alignment_expands_internal_branch_degeneracy():
    elements, wbo = _tetramethyl_metal_graph()

    pool = run_cut_sweep_chunk(
        elements, wbo, elements, wbo, [()],
        core_R=[1, 5], n_workers=1, n_seeds=1,
        max_branches=100, iso_tol=0.5, symmetry_wbo_tol=0.2)

    mappings = {
        tuple(sorted(info["mapping"].items()))
        for info in pool.values()
    }
    expected_targets = {1, 5, 9, 13}

    assert len(mappings) == 12
    assert mappings == {
        ((1, a), (5, b))
        for a in expected_targets
        for b in expected_targets
        if a != b
    }


def test_core_mapping_variants_do_not_invent_global_orbit_swaps():
    class Branch:
        mapping = {0: 0, 1: 1}
        symmetry_fragments = []

    wbo = np.zeros((2, 2))
    g_p = build_graph(["H", "H"], wbo, bond_cut=0.2)

    variants = _core_mapping_variants(
        Branch(), [0, 1], 100,
        g_P=g_p, p_orbits={0: 0, 1: 0})

    assert variants == [{0: 0, 1: 1}]


def test_core_alignment_does_not_filter_nonmatching_core_edges():
    el_r = ["N", "C", "O"]
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 1] = wbo_r[1, 0] = 0.79
    wbo_r[1, 2] = wbo_r[2, 1] = 1.10

    el_t = ["N", "C", "O"]
    wbo_t = np.zeros((3, 3))
    wbo_t[0, 1] = wbo_t[1, 0] = 0.0
    wbo_t[1, 2] = wbo_t[2, 1] = 1.67

    pool = run_cut_sweep_chunk(
        el_r, wbo_r, el_t, wbo_t, [()],
        core_R=[0, 1, 2], n_workers=1, n_seeds=1,
        graph_floor=0.2, iso_tol=0.1, max_branches=100)

    assert {
        tuple(sorted(info["mapping"].items()))
        for info in pool.values()
    } == {((0, 0), (1, 1), (2, 2))}


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


def test_nauty_loose_tolerance_preserves_edges_and_methyl_orbit():
    pytest.importorskip("pynauty")
    wbo = np.zeros((4, 4))
    for h, value in enumerate((0.982056, 0.980196, 0.960315), start=1):
        wbo[0, h] = wbo[h, 0] = value
    g = build_graph(["C", "H", "H", "H"], wbo, bond_cut=0.2)

    buckets, zero_bucket = _wbo_tolerance_bucket_lookup(g, 1.0)
    orbits = _nauty_orbits(g, wbo_tol=1.0)

    assert all(buckets[(0, h)] != zero_bucket for h in (1, 2, 3))
    assert buckets[(1, 2)] == zero_bucket
    assert orbits[1] == orbits[2] == orbits[3]
    assert orbits[0] != orbits[1]


def test_nauty_masks_positive_wbo_below_graph_floor():
    pytest.importorskip("pynauty")
    wbo = np.zeros((5, 5))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[0, 2] = wbo[2, 0] = 1.0
    wbo[1, 3] = wbo[3, 1] = 1.0
    wbo[2, 4] = wbo[4, 2] = 1.0
    wbo[1, 4] = wbo[4, 1] = 0.16
    graph = build_graph(["B", "O", "O", "H", "H"], wbo, bond_cut=0.2)

    buckets, zero_bucket = _wbo_tolerance_bucket_lookup(graph, 1.0)
    orbits = _nauty_orbits(graph, wbo_tol=1.0)

    assert buckets[(1, 4)] == zero_bucket
    assert orbits[1] == orbits[2]
    assert orbits[3] == orbits[4]


def test_chosen_automorphism_closure_recovers_early_singleton():
    pytest.importorskip("pynauty")
    wbo = np.zeros((4, 4))
    for h, value in enumerate((0.982056, 0.980196, 0.960315), start=1):
        wbo[0, h] = wbo[h, 0] = value
    graph = build_graph(["C", "H", "H", "H"], wbo, bond_cut=0.2)
    partial = {
        "fragments": [{"fragment_index": 0, "island_idx": 0,
                       "fragment": [0, 1, 2, 3]}],
        "blocks": [{
            "r_atoms": [1, 2], "p_atoms": [1, 2], "source": "sym_block",
        }],
    }

    complete = complete_chosen_automorphism_groups(
        partial, {0: 0, 1: 1, 2: 2, 3: 3}, graph, graph, 1.0)

    assert complete["color_groups"] == [{
        "r_atoms": [1, 2, 3],
        "p_atoms": [1, 2, 3],
        "sources": ["chosen_candidate_automorph"],
    }]


def test_stored_branch_group_preserves_correlated_local_orbits():
    wbo = np.zeros((5, 5))
    for atom in range(1, 5):
        wbo[0, atom] = wbo[atom, 0] = 1.0
    graph = build_graph(["C", "H", "H", "F", "F"], wbo, bond_cut=0.2)
    hierarchy = {"fragments": [{
        "fragment_index": 0,
        "island_idx": 1,
        "fragment": [0, 1, 2, 3, 4],
    }]}

    complete = complete_chosen_automorphism_groups(
        hierarchy, {atom: atom for atom in range(5)},
        graph, graph, 1.0,
        exact_target_generators=[
            [0, 2, 1, 3, 4],
            [0, 1, 2, 4, 3],
        ])

    groups = {
        (block.get("center_R"), tuple(block["r_atoms"]))
        for block in complete["blocks"]
    }
    assert (0, (1, 2)) in groups
    assert (0, (3, 4)) in groups
    assert all(block["source"] == "stored_AAM_branch_mapping_group"
               for block in complete["blocks"])


def test_cross_branch_assignment_relation_preserves_observed_shuffles():
    wbo = np.zeros((5, 5))
    for atom in range(1, 5):
        wbo[0, atom] = wbo[atom, 0] = 1.0
    graph = build_graph(["M", "H", "H", "F", "F"], wbo, bond_cut=0.2)
    hierarchy = {"fragments": [{
        "fragment_index": 0,
        "fragment": [0, 1, 2, 3, 4],
    }]}
    selected = {atom: atom for atom in range(5)}

    complete = complete_chosen_automorphism_groups(
        hierarchy, selected, graph, graph, 1.0,
        exact_target_generators=[],
        exact_branch_mappings=[
            selected,
            {0: 0, 1: 2, 2: 1, 3: 3, 4: 4},
            {0: 0, 1: 1, 2: 2, 3: 4, 4: 3},
        ])

    cross_branch = {
        (block["center_R"], tuple(block["r_atoms"]))
        for block in complete["blocks"]
        if block["source"] == "AAM_cross_branch_assignment"
    }
    assert cross_branch == {(0, (1, 2)), (0, (3, 4))}


def test_completed_candidate_groups_are_cached_after_branch_reduction(
        monkeypatch):
    from rxn_core.matcher.canonical import (
        _CandidateAutomorphismCanonicalizer,
    )

    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 1.0
    graph = build_graph(["H", "H"], wbo, bond_cut=0.2)
    state = {
        "witness": {0: 0, 1: 1},
        "blocks": [{
            "r_atoms": [0, 1],
            "p_atoms": [0, 1],
            "source": "exact_automorph_group",
        }],
        "automorph_blocks": [{
            "r_atoms": [0, 1],
            "p_atoms": [0, 1],
        }],
        "multiplicity": 2,
    }
    branch = {"hierarchy": {"fragments": [{
        "fragment_index": 0,
        "fragment": [0, 1],
        "symmetry": state,
    }]}}
    calls = 0
    original = _CandidateAutomorphismCanonicalizer.atom_generators

    def measured(self, candidate):
        nonlocal calls
        calls += 1
        return original(self, candidate)

    monkeypatch.setattr(
        _CandidateAutomorphismCanonicalizer, "atom_generators", measured)
    completed, metrics = attach_completed_candidate_groups(
        [branch, branch], graph, wbo_tol=1.0, return_metrics=True)

    assert calls == 1
    assert metrics == {
        "completed_candidate_group_requests": 2,
        "completed_candidate_group_calculations": 1,
        "completed_candidate_group_cache_hits": 1,
    }
    for item in completed:
        symmetry = item["hierarchy"]["fragments"][0]["symmetry"]
        assert symmetry["automorph_generators"] == [[1, 0]]
        assert symmetry["automorph_group_source"] == (
            "completed_candidate_after_branch_family_reduction")


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


def test_nonautomorphic_witnesses_are_not_compressed_before_later_frontier():
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 2] = wbo_r[2, 0] = 1.0
    g_r = build_graph(["C", "C", "O"], wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((3, 3))
    wbo_p[1, 2] = wbo_p[2, 1] = 1.0
    g_p = build_graph(["C", "C", "O"], wbo_p, bond_cut=0.2)

    cands = _dedup_sym_cands(
        [_SymCand({0: 0, 1: 1}), _SymCand({0: 1, 1: 0})],
        g_r, g_p,
    )
    assert len(cands) == 2

    out = _extend_sym_cands(
        cands, {0, 1}, 2,
        g_r, g_p, {}, 0.5, None,
        anchor_u=0, anchor_wbo=1.0,
    )

    assert len(out) == 1
    assert _cand_map(out[0]) == {0: 1, 1: 0, 2: 2}


def test_exact_certificate_preserves_coupled_group_action_not_just_orbits():
    wbo = np.zeros((4, 4))
    for a, b in [(0, 1), (1, 2), (2, 3), (3, 0)]:
        wbo[a, b] = wbo[b, a] = 1.0
    graph = build_graph(["C"] * 4, wbo, bond_cut=0.2)
    orbits = _nauty_orbits(graph, wbo_tol=0.2)

    # Every P atom belongs to one orbit, but a square automorphism cannot map
    # an adjacent ordered pair onto an opposite pair.  Independent orbit-pool
    # shuffling would merge these states incorrectly.
    adjacent = _SymCand({10: 0, 11: 1})
    rotated = _SymCand({10: 1, 11: 2})
    opposite = _SymCand({10: 0, 11: 2})
    deduped = _dedup_sym_cands(
        [adjacent, rotated, opposite], graph, graph,
        r_orbits=orbits, p_orbits=orbits,
    )

    assert len(deduped) == 2
    assert sorted(cand.multiplicity for cand in deduped) == [1, 2]


def test_successful_primary_extension_preserves_exact_group_evidence():
    elements = ["X", "C", "C", "H", "H"]
    wbo = np.zeros((5, 5))
    for carbon in [1, 2]:
        wbo[0, carbon] = wbo[carbon, 0] = 1.0
    wbo[1, 3] = wbo[3, 1] = 1.0
    wbo[2, 4] = wbo[4, 2] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)
    r_orbits = _nauty_orbits(g, wbo_tol=0.2)
    p_orbits = _nauty_orbits(g, wbo_tol=0.2)

    cands = [_SymCand({0: 0})]
    for fragment, atom, anchor in [
        ({0}, 1, 0),
        ({0, 1}, 2, 0),
        ({0, 1, 2}, 3, 1),
        ({0, 1, 2, 3}, 4, 2),
    ]:
        cands = _extend_sym_cands(
            cands, fragment, atom, g, g, {}, 0.5, None,
            p_orbits=p_orbits, r_orbits=r_orbits,
            anchor_u=anchor, anchor_wbo=1.0)

    assert len(cands) == 1
    assert cands[0].multiplicity == 2
    assert [(b.r_atoms, b.p_atoms) for b in cands[0].automorph_blocks] == [
        ((1, 2, 3, 4), (1, 2, 3, 4)),
    ]
    state = _symmetry_state(cands[0], r_orbits=r_orbits, p_orbits=p_orbits)
    assert state["automorph_blocks"][0]["source"] == "exact_automorph_group"


def test_deferred_boundary_dedupe_retains_concrete_witnesses_until_saturation():
    elements = ["X", "C", "C", "H", "H"]
    wbo = np.zeros((5, 5))
    for carbon in [1, 2]:
        wbo[0, carbon] = wbo[carbon, 0] = 1.0
    wbo[1, 3] = wbo[3, 1] = 1.0
    wbo[2, 4] = wbo[4, 2] = 1.0
    graph = build_graph(elements, wbo, bond_cut=0.2)
    r_orbits = _nauty_orbits(graph, wbo_tol=0.2)
    p_orbits = _nauty_orbits(graph, wbo_tol=0.2)

    cands = [_SymCand({0: 0})]
    for fragment, atom, anchor in [
        ({0}, 1, 0),
        ({0, 1}, 2, 0),
        ({0, 1, 2}, 3, 1),
        ({0, 1, 2, 3}, 4, 2),
    ]:
        cands = _extend_sym_cands(
            cands, fragment, atom, graph, graph, {}, 0.5, None,
            p_orbits=p_orbits, r_orbits=r_orbits,
            anchor_u=anchor, anchor_wbo=1.0,
            defer_boundary_dedupe=True)

    # No witness was quotient-merged during growth, so both correlated
    # concrete paths remain independently available to later atoms.
    assert len(cands) == 2
    assert {tuple(sorted(cand.mapping.items())) for cand in cands} == {
        tuple(sorted({0: 0, 1: 1, 2: 2, 3: 3, 4: 4}.items())),
        tuple(sorted({0: 0, 1: 2, 2: 1, 3: 4, 4: 3}.items())),
    }

    saturated = _dedup_sym_cands(
        cands, graph, graph, r_orbits=r_orbits, p_orbits=p_orbits,
        fragment=set(range(5)))
    assert len(saturated) == 1
    assert saturated[0].multiplicity == 2
    assert [(b.r_atoms, b.p_atoms) for b in saturated[0].automorph_blocks] == [
        ((1, 2, 3, 4), (1, 2, 3, 4)),
    ]


def test_block_join_refines_bucket_mismatch_before_later_frontier():
    wbo_r = np.zeros((4, 4))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.907
    wbo_r[0, 2] = wbo_r[2, 0] = 1.869
    wbo_r[0, 3] = wbo_r[3, 0] = 1.138
    g_r = build_graph(["W", "O", "O", "O"], wbo_r, bond_cut=0.2)

    wbo_p = np.zeros((4, 4))
    wbo_p[0, 1] = wbo_p[1, 0] = 1.910
    wbo_p[0, 2] = wbo_p[2, 0] = 1.875
    wbo_p[0, 3] = wbo_p[3, 0] = 1.127
    g_p = build_graph(["W", "O", "O", "O"], wbo_p, bond_cut=0.2)

    r_orbits = {0: 0, 1: 1, 2: 1, 3: 2}
    p_orbits = {0: 0, 1: 1, 2: 1, 3: 2}
    cand = _SymCand(
        {0: 0, 1: 1},
        (_SymBlock((1,), (1, 2), extendable=True),),
    )

    strong = _extend_sym_cands(
        [cand], {0, 1}, 2,
        g_r, g_p, {}, 1.0, None,
        p_orbits=p_orbits, r_orbits=r_orbits,
        anchor_u=0, anchor_wbo=1.869,
    )
    later = _extend_sym_cands(
        strong, {0, 1, 2}, 3,
        g_r, g_p, {}, 1.0, None,
        p_orbits=p_orbits, r_orbits=r_orbits,
        anchor_u=0, anchor_wbo=1.138,
    )

    assert any(
        _cand_map(c).get(2) in {1, 2} and _cand_map(c).get(3) == 3
        for c in later
    )


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


def test_symmetry_repair_rejects_non_automorphic_orbit_swap():
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
    assert (len(after[0]), len(after[1])) == (2, 2)
    assert repaired == bad
    assert stats["repaired"] is False


def test_symmetry_repair_accepts_strict_local_automorphism():
    elements = ["C", "H", "H"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 0.3
    wbo[0, 2] = wbo[2, 0] = 1.0
    g_r = build_graph(elements, wbo, bond_cut=0.2)
    g_p = build_graph(elements, wbo, bond_cut=0.2)
    p_orbits = _nauty_orbits(g_p, wbo_tol=1.0)

    bad = {0: 0, 1: 2, 2: 1}
    repaired, stats = symmetry_repair_mapping(
        bad, wbo, wbo, g_r, g_p, p_orbits,
        min_changes=1, return_stats=True)
    after = classify_bonds(repaired, wbo, wbo)

    assert repaired == {0: 0, 1: 1, 2: 2}
    assert (len(after[0]), len(after[1])) == (0, 0)
    assert stats["repaired"] is True


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


def test_mechanism_certificate_distinguishes_pair_orbits_not_vertex_orbits():
    # A distance-colored K5 is vertex-transitive, but cycle edges and
    # diagonals are different pair orbits.  Endpoint-orbit IDs alone conflate
    # these two one-bond events.
    n = 5
    elements = ["C"] * n
    wbo_r = np.full((n, n), 0.5)
    np.fill_diagonal(wbo_r, 0.0)
    for atom in range(n):
        neighbor = (atom + 1) % n
        wbo_r[atom, neighbor] = wbo_r[neighbor, atom] = 1.0
    wbo_cycle = wbo_r.copy()
    wbo_cycle[0, 1] = wbo_cycle[1, 0] = 0.0
    wbo_diagonal = wbo_r.copy()
    wbo_diagonal[0, 2] = wbo_diagonal[2, 0] = 0.0
    graph_r = build_graph(elements, wbo_r, bond_cut=0.2)
    r_orbits = _nauty_orbits(graph_r, wbo_tol=0.2)
    assert len(set(r_orbits.values())) == 1
    assert _orbit_bond_key([(0, 1)], r_orbits, "R") == (
        _orbit_bond_key([(0, 2)], r_orbits, "R"))
    mapping = {atom: atom for atom in range(n)}

    def signature(product_wbo):
        graph_p = build_graph(elements, product_wbo, bond_cut=0.2)
        return _mechanism_signature(
            mapping, wbo_r, product_wbo, r_orbits,
            _nauty_orbits(graph_p, wbo_tol=0.2),
            elements_R=elements, elements_P=elements,
            g_R_full=graph_r, symmetry_wbo_tol=0.2)

    cycle_signature = signature(wbo_cycle)
    diagonal_signature = signature(wbo_diagonal)
    assert len(cycle_signature[0]) == len(diagonal_signature[0]) == 1
    assert cycle_signature != diagonal_signature


def test_cut_sweep_pool_prefers_no_cut_representative():
    pool = {}
    sig = (("broken",), ("formed",))

    _pool_add(pool, sig, {0: 1}, ((0, 1),))
    _pool_add(pool, sig, {0: 0}, ())

    assert pool[sig]["mapping"] == {0: 0}
    assert pool[sig]["has_no_cut"] is True
    assert pool[sig]["cuts"] == frozenset({(0, 1)})
    assert pool[sig]["dedup_count"] == 2
    assert pool[sig]["branch_symmetry"]["analytical_branch_count"] == 2
    assert {
        tuple(sorted(branch["mapping"].items()))
        for branch in pool[sig]["branches"]
    } == {((0, 1),), ((0, 0),)}


def test_pool_keeps_distinct_completed_branches_under_one_mechanism():
    pool = {}
    sig = ((), ())

    _pool_add(pool, sig, {0: 0, 1: 1}, ())
    _pool_add(pool, sig, {0: 1, 1: 0}, ())

    assert pool[sig]["branch_symmetry"]["analytical_branch_count"] == 2
    assert len(pool[sig]["branches"]) == 2
    assert "matching_generators" not in pool[sig]["branch_symmetry"]
    assert "dedup_witnesses" not in pool[sig]


def test_display_symmetry_uses_only_selected_final_candidate_blocks():
    pool = {}
    sig = ((), ())
    selected_local = {
        "blocks": [
            {"r_atoms": [0, 1], "p_atoms": [10, 11]},
            {
                "r_atoms": [2, 3],
                "p_atoms": [12, 13],
                "source": "island_automorph",
            },
        ],
        "fragments": [],
    }
    other_local = {
        "blocks": [{"r_atoms": [4, 5], "p_atoms": [14, 15]}],
        "fragments": [],
    }

    _pool_add(pool, sig, {0: 10, 1: 11}, (), selected_local)
    _pool_add(pool, sig, {0: 11, 1: 10}, (), other_local)

    symmetry = pool[sig]["branch_symmetry"]
    assert symmetry["selected_branch_index"] == 0
    assert symmetry["analytical_branch_count"] == 2
    assert symmetry["color_groups"] == [{
        "r_atoms": [0, 1],
        "p_atoms": [10, 11],
        "sources": ["sym_block"],
    }]


def test_color_groups_do_not_transitively_merge_overlapping_blocks():
    groups = _color_groups_from_blocks([
        {
            "r_atoms": [1, 2],
            "p_atoms": [10, 11],
            "source": "sym_block",
        },
        {
            "r_atoms": [2, 3],
            "p_atoms": [11, 12],
            "source": "island_automorph",
        },
        {
            "r_atoms": [1, 2],
            "p_atoms": [10, 11],
            "source": "exact_automorph_group",
        },
    ])

    assert groups == [
        {
            "r_atoms": [1, 2],
            "p_atoms": [10, 11],
            "sources": ["exact_automorph_group", "sym_block"],
        },
        {
            "r_atoms": [2, 3],
            "p_atoms": [11, 12],
            "sources": ["island_automorph"],
        },
    ]


def test_generate_seed_orders_honors_trial_cap():
    elements = ["C", "C", "C", "C", "H", "H"]
    wbo = np.zeros((6, 6))
    for i, j in [(0, 1), (1, 2), (2, 3), (0, 4), (3, 5)]:
        wbo[i, j] = wbo[j, i] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)

    orders = _generate_seed_orders(g, n_trials=2)

    assert len(orders) == 2
    assert all(g.degree[order[0]] > 0 for order in orders)


def test_generate_seed_orders_deprioritizes_common_isolated_atoms():
    elements = ["C", "C", "C", "H", "H", "H", "H"]
    wbo = np.zeros((7, 7))
    for i, j in [(0, 1), (1, 2), (2, 3)]:
        wbo[i, j] = wbo[j, i] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)

    order = _generate_seed_orders(g, n_trials=1)[0]

    connected = {0, 1, 2, 3}
    isolated_common = {4, 5, 6}
    assert max(order.index(atom) for atom in connected) < min(
        order.index(atom) for atom in isolated_common)
    assert set(order) == set(range(7))


def test_core_atoms_do_not_reorder_seed_sequence(monkeypatch):
    elements = ["C", "C", "O"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[1, 2] = wbo[2, 1] = 1.0
    g = build_graph(elements, wbo, bond_cut=0.2)
    seen = []

    def fake_grow_island(_g_R, _g_P, seed, *_args, **_kwargs):
        seen.append(seed)
        return []

    monkeypatch.setattr(branch_mod, "grow_island", fake_grow_island)

    branch_mod.find_islands(
        g, g, [0, 1, 2],
        core_R=[2],
        stop_when_core_mapped=True,
        p_orbits={0: 0, 1: 1, 2: 2},
        r_orbits={0: 0, 1: 1, 2: 2},
    )

    assert seen == [0, 1, 2]


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


def test_cut_sweep_respects_hard_anchor_map():
    elements = ["O", "C", "O"]
    wbo = np.zeros((3, 3))
    wbo[0, 1] = wbo[1, 0] = 1.0
    wbo[1, 2] = wbo[2, 1] = 1.0

    pool = cut_sweep(
        elements, wbo, elements, wbo,
        n_workers=0, n_seeds=1, max_branches=100,
        cut_floor=0.2, symmetry_repair=True,
        anchor_map={0: 2},
    )
    minimal = select_min_mechanisms(pool)

    assert minimal
    assert all(info["mapping"][0] == 2 for info in minimal.values())
    assert all(info["mapping"][2] == 0 for info in minimal.values())

    incompatible = cut_sweep(
        elements, wbo, elements, wbo,
        n_workers=0, n_seeds=1, max_branches=100,
        cut_floor=0.2, symmetry_repair=True,
        anchor_map={0: 1},
    )

    assert incompatible == {}


def test_cut_sweep_preserves_branch_exact_automorph_group():
    elements = ["H", "H"]
    wbo = np.zeros((2, 2))
    wbo[0, 1] = wbo[1, 0] = 0.9

    pool = cut_sweep(
        elements, wbo, elements, wbo,
        n_workers=0, n_seeds=1, max_branches=100,
        cut_floor=0.2, symmetry_repair=False,
    )
    info = next(iter(pool.values()))
    branch_symmetry = info["branch_symmetry"]

    assert branch_symmetry["rule"] == "selected_analytical_branch"
    assert branch_symmetry["analytical_branch_count"] == 2
    assert any(
        block == {
            "fragment_index": 0,
            "block_index": 0,
            "island_idx": 1,
            "r_atoms": [0, 1],
            "p_atoms": [0, 1],
            "extendable": False,
            "open": False,
            "assignments": "exact_group",
            "source": "exact_automorph_group",
        }
        for block in branch_symmetry["blocks"]
    )
    assert branch_symmetry["color_groups"] == [{
        "r_atoms": [0, 1],
            "p_atoms": [0, 1],
            "sources": ["exact_automorph_group"],
    }]
    fragment_symmetry = branch_symmetry["fragments"][0]["symmetry"]
    assert fragment_symmetry["automorph_blocks"][0]["source"] == (
        "exact_automorph_group")
    assert "automorph_generators" not in fragment_symmetry


def test_anchored_noop_seed_does_not_keep_pass_loop_alive():
    elements = ["C", "C", "C"]
    wbo_r = np.zeros((3, 3))
    wbo_r[0, 1] = wbo_r[1, 0] = 1.0
    wbo_r[1, 2] = wbo_r[2, 1] = 1.0

    wbo_p = np.zeros((3, 3))
    wbo_p[1, 2] = wbo_p[2, 1] = 1.0

    g_r = build_graph(elements, wbo_r, bond_cut=0.2)
    g_p = build_graph(elements, wbo_p, bond_cut=0.2)
    events = []

    branches = branch_mod.find_islands(
        g_r, g_p, [1, 0],
        iso_tol=1.0,
        max_branches=100,
        anchor_map={0: 0},
        events=events,
    )

    assert branches
    assert all(branch.mapping[0] == 0 for branch in branches)
    assert [e["seed"] for e in events if e.get("type") == "seed_start"][:2] == [1, 0]
    assert sum(1 for e in events if e.get("type") == "pass_start") == 2

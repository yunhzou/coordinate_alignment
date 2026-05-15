"""Sweep-cut mechanism discovery for WBO graph alignment.

R-P sweep cut is part of the core alignment algorithm: mechanism discovery
tries the no-cut graph and every one-edge R cut above a WBO floor, then dedups
the resulting witnesses by symmetry-canonical broken/formed bond signatures.
"""
from __future__ import annotations

import multiprocessing as mp

from ..frag import build_graph, classify_bonds, expand_mapping
from ..matcher import _nauty_orbits
from .branch import (
    BranchLimitExceeded,
    _generate_seed_orders,
    find_islands,
    symmetry_repair_mapping,
)


def _canon_pair(a, b):
    return (a, b) if a <= b else (b, a)


def _strong_edges(wboR, cut_floor):
    n = int(wboR.shape[0])
    return [(i, j) for i in range(n) for j in range(i + 1, n)
            if float(wboR[i, j]) >= cut_floor]


def _orbit_bond_key(pairs, orbits, tag):
    return tuple(sorted(
        (tag, *_canon_pair(int(orbits[a]), int(orbits[b])))
        for a, b in pairs
    ))


def _core_mapping_key(mapping, core_R):
    return (
        tuple((int(r), int(mapping[r])) for r in sorted(core_R)),
        (),
    )


def _mechanism_signature(mapping, wboR, wboT, r_orbits, p_orbits,
                         dwbo_threshold=0.5,
                         elements_R=None, elements_P=None,
                         metal_dwbo_threshold=None):
    """Symmetry-canonical mechanism key for R-P discovery."""
    broken, formed, _, _ = classify_bonds(
        mapping, wboR, wboT, dwbo_threshold=dwbo_threshold,
        elements_R=elements_R, elements_P=elements_P,
        metal_dwbo_threshold=metal_dwbo_threshold)
    inv = {v: k for k, v in mapping.items()}
    br_pairs = [(a, b) for (a, b, _, _) in broken]
    fm_r_pairs = []
    fm_p_pairs = []
    for a, b, _, _ in formed:
        if a in inv and b in inv:
            fm_r_pairs.append((inv[a], inv[b]))
        else:
            fm_p_pairs.append((a, b))
    br = _orbit_bond_key(br_pairs, r_orbits, 'R')
    fm = (_orbit_bond_key(fm_r_pairs, r_orbits, 'R') +
          _orbit_bond_key(fm_p_pairs, p_orbits, 'P'))
    return br, tuple(sorted(fm))


def _pool_add(pool, sig, mapping, cuts):
    cuts = frozenset(cuts)
    no_cut = not cuts
    entry = pool.get(sig)
    if entry is None:
        pool[sig] = {
            'mapping': dict(mapping),
            'cuts': cuts,
            'has_no_cut': bool(no_cut),
            'dedup_count': 1,
        }
    else:
        if no_cut and not entry.get('has_no_cut', False):
            entry['mapping'] = dict(mapping)
        entry['has_no_cut'] = bool(entry.get('has_no_cut', False) or no_cut)
        entry['cuts'] = entry['cuts'] | cuts
        entry['dedup_count'] = entry.get('dedup_count', 1) + 1


def _run_find_islands_limited(g_R, g_P, order, core_R, cfg, *,
                              p_orbits=None, r_orbits=None):
    stop_on_core = bool(core_R)
    return find_islands(
        g_R, g_P, list(order),
        iso_tol=float(cfg['iso_tol']),
        max_branches=int(cfg['max_branches']),
        abort_on_branch_cap=bool(cfg.get('abort_on_branch_cap', False)),
        dwbo_threshold=float(cfg['dwbo_threshold']),
        metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'),
        symmetry_wbo_tol=float(cfg['symmetry_wbo_tol']),
        core_R=core_R,
        stop_when_core_mapped=stop_on_core,
        p_orbits=p_orbits,
        r_orbits=r_orbits,
    )


def _score_branch_mapping(mapping, g_R, g_P, wboR, wboT,
                          g_R_full, p_orbits, r_orbits, core_R, cfg,
                          elR=None, elT=None):
    if core_R:
        if not all(r in mapping for r in core_R):
            return None
        return _core_mapping_key(mapping, core_R), mapping

    if len(mapping) < int(cfg['n_atoms']) - 2:
        return None
    if cfg['symmetry_repair']:
        mapping = symmetry_repair_mapping(
            mapping, wboR, wboT, g_R_full, g_P, p_orbits,
            dwbo_threshold=float(cfg['dwbo_threshold']),
            metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'),
            min_changes=int(cfg['symmetry_repair_min_changes']),
            max_evals=int(cfg['symmetry_repair_max_evals']),
        )
    sig = _mechanism_signature(
        mapping, wboR, wboT, r_orbits, p_orbits,
        dwbo_threshold=float(cfg['dwbo_threshold']),
        elements_R=elR,
        elements_P=elT,
        metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'))
    return sig, mapping


_WORKER = {}


def _cs_winit(elR, wboR, elT, wboT, cfg):
    graph_floor = float(cfg['graph_floor'])
    _WORKER['elR'] = elR
    _WORKER['wboR'] = wboR
    _WORKER['elT'] = elT
    _WORKER['wboT'] = wboT
    _WORKER['cfg'] = dict(cfg)
    _WORKER['g_P'] = build_graph(elT, wboT, bond_cut=graph_floor)
    _WORKER['g_R_full'] = build_graph(elR, wboR, bond_cut=graph_floor)
    symmetry_wbo_tol = float(cfg['symmetry_wbo_tol'])
    _WORKER['p_orbits'] = _nauty_orbits(
        _WORKER['g_P'], wbo_tol=symmetry_wbo_tol)
    _WORKER['r_orbits'] = _nauty_orbits(
        _WORKER['g_R_full'], wbo_tol=symmetry_wbo_tol)


def _cs_wrun(args):
    cut, orders, core_R = args
    cfg = _WORKER['cfg']
    graph_floor = float(cfg['graph_floor'])
    g_R = build_graph(_WORKER['elR'], _WORKER['wboR'],
                      bond_cut=graph_floor)
    for i, j in cut:
        if g_R.has_edge(i, j):
            g_R.remove_edge(i, j)
    r_orbits_cut = _nauty_orbits(
        g_R, wbo_tol=float(cfg['symmetry_wbo_tol']))

    out = []
    try:
        for order in orders:
            try:
                branches = _run_find_islands_limited(
                    g_R, _WORKER['g_P'], order, core_R, cfg,
                    p_orbits=_WORKER['p_orbits'],
                    r_orbits=r_orbits_cut)
            except BranchLimitExceeded:
                raise
            except Exception:
                continue
            for branch in branches:
                mapping = expand_mapping(dict(branch.mapping), g_R, _WORKER['g_P'])
                scored = _score_branch_mapping(
                    mapping, g_R, _WORKER['g_P'],
                    _WORKER['wboR'], _WORKER['wboT'],
                    _WORKER['g_R_full'], _WORKER['p_orbits'],
                    _WORKER['r_orbits'], core_R, cfg,
                    _WORKER['elR'], _WORKER['elT'])
                if scored is None:
                    continue
                sig, mapping = scored
                out.append((sig, tuple(sorted(mapping.items())), cut))
    except BranchLimitExceeded:
        # Kill the whole cut.  Partial seed-order witnesses from this cut are
        # deliberately discarded because the cut has entered a pathological
        # outer-branch multiplication regime.
        return []
    return out


def cut_sweep_items(wboR, cut_floor=0.2):
    """Return the independent no-cut plus one-edge cut work items."""
    return [()] + [((int(i), int(j)),) for i, j in _strong_edges(
        wboR, float(cut_floor))]


def _cut_sweep_cfg(*, cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
                   dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
                   symmetry_wbo_tol=0.2,
                   n_seeds=3, max_branches=100,
                   chunksize=1,
                   symmetry_repair=True,
                   symmetry_repair_min_changes=1,
                   symmetry_repair_max_evals=20000,
                   n_atoms=0):
    return {
        'cut_floor': float(cut_floor),
        'graph_floor': float(graph_floor),
        'iso_tol': float(iso_tol),
        'dwbo_threshold': float(dwbo_threshold),
        'metal_dwbo_threshold': (
            None if metal_dwbo_threshold is None
            else float(metal_dwbo_threshold)
        ),
        'symmetry_wbo_tol': float(symmetry_wbo_tol),
        'n_seeds': int(n_seeds),
        'max_branches': int(max_branches),
        'chunksize': int(chunksize),
        'symmetry_repair': bool(symmetry_repair),
        'symmetry_repair_min_changes': int(symmetry_repair_min_changes),
        'symmetry_repair_max_evals': int(symmetry_repair_max_evals),
        'abort_on_branch_cap': True,
        'n_atoms': int(n_atoms),
    }


def _cut_sweep_chunk_serial(elR, wboR, elT, wboT, cfg, core_R, cuts):
    graph_floor = float(cfg['graph_floor'])
    g_P = build_graph(elT, wboT, bond_cut=graph_floor)
    g_R_full = build_graph(elR, wboR, bond_cut=graph_floor)
    symmetry_wbo_tol = float(cfg['symmetry_wbo_tol'])
    p_orbits = _nauty_orbits(g_P, wbo_tol=symmetry_wbo_tol)
    r_orbits = _nauty_orbits(g_R_full, wbo_tol=symmetry_wbo_tol)
    pool = {}

    def run(cuts):
        g_R = build_graph(elR, wboR, bond_cut=graph_floor)
        for i, j in cuts:
            if g_R.has_edge(i, j):
                g_R.remove_edge(i, j)
        r_orbits_cut = _nauty_orbits(g_R, wbo_tol=symmetry_wbo_tol)
        orders = _generate_seed_orders(g_R, n_trials=int(cfg['n_seeds']))
        cut_hits = []
        try:
            for order in orders:
                try:
                    branches = _run_find_islands_limited(
                        g_R, g_P, order, core_R, cfg,
                        p_orbits=p_orbits, r_orbits=r_orbits_cut)
                except BranchLimitExceeded:
                    raise
                except Exception:
                    continue
                for branch in branches:
                    mapping = expand_mapping(dict(branch.mapping), g_R, g_P)
                    scored = _score_branch_mapping(
                        mapping, g_R, g_P, wboR, wboT, g_R_full,
                        p_orbits, r_orbits, core_R, cfg, elR, elT)
                    if scored is None:
                        continue
                    sig, mapping = scored
                    cut_hits.append((sig, mapping))
        except BranchLimitExceeded:
            return
        for sig, mapping in cut_hits:
            _pool_add(pool, sig, mapping, cuts)

    for cut in cuts:
        run(tuple(tuple(pair) for pair in cut))
    return pool


def _cut_sweep_serial(elR, wboR, elT, wboT, cfg, core_R):
    return _cut_sweep_chunk_serial(
        elR, wboR, elT, wboT, cfg, core_R,
        cut_sweep_items(wboR, cfg['cut_floor']))


def _cut_sweep_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R):
    graph_floor = float(cfg['graph_floor'])
    g_R = build_graph(elR, wboR, bond_cut=graph_floor)
    orders = [tuple(order) for order in _generate_seed_orders(
        g_R, n_trials=int(cfg['n_seeds']))]
    cuts = cut_sweep_items(wboR, cfg['cut_floor'])
    work = [(cut, orders, core_R) for cut in cuts]
    pool = {}
    with mp.Pool(n_workers, initializer=_cs_winit,
                 initargs=(elR, wboR, elT, wboT, cfg)) as proc_pool:
        for results in proc_pool.imap_unordered(
                _cs_wrun, work, chunksize=max(1, int(cfg['chunksize']))):
            for sig, mapping_items, cut in results:
                _pool_add(pool, sig, dict(mapping_items), cut)
    return pool


def _cut_sweep_chunk_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R,
                              cuts):
    graph_floor = float(cfg['graph_floor'])
    g_R = build_graph(elR, wboR, bond_cut=graph_floor)
    orders = [tuple(order) for order in _generate_seed_orders(
        g_R, n_trials=int(cfg['n_seeds']))]
    work = [(cut, orders, core_R) for cut in cuts]
    pool = {}
    with mp.Pool(n_workers, initializer=_cs_winit,
                 initargs=(elR, wboR, elT, wboT, cfg)) as proc_pool:
        for results in proc_pool.imap_unordered(
                _cs_wrun, work, chunksize=max(1, int(cfg['chunksize']))):
            for sig, mapping_items, cut in results:
                _pool_add(pool, sig, dict(mapping_items), cut)
    return pool


def run_cut_sweep_chunk(elR, wboR, elT, wboT, cuts, *,
                        core_R=None,
                        n_workers=None,
                        cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
                        dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
                        symmetry_wbo_tol=0.2,
                        n_seeds=3, max_branches=100,
                        chunksize=1,
                        symmetry_repair=True,
                        symmetry_repair_min_changes=1,
                        symmetry_repair_max_evals=20000):
    """Run a chunk of independent cut-sweep work items.

    This is the Slurm-array friendly primitive.  The caller chooses which cut
    work items belong to the chunk; the merge step combines the returned pools.
    """
    cfg = _cut_sweep_cfg(
        cut_floor=cut_floor,
        graph_floor=graph_floor,
        iso_tol=iso_tol,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
        symmetry_wbo_tol=symmetry_wbo_tol,
        n_seeds=n_seeds,
        max_branches=max_branches,
        chunksize=chunksize,
        symmetry_repair=symmetry_repair,
        symmetry_repair_min_changes=symmetry_repair_min_changes,
        symmetry_repair_max_evals=symmetry_repair_max_evals,
        n_atoms=len(elR),
    )
    core_R = tuple(sorted(set(core_R or ())))
    normalized_cuts = [
        tuple(tuple(int(v) for v in pair) for pair in cut)
        for cut in cuts
    ]
    if n_workers and int(n_workers) > 1 and len(normalized_cuts) > 1:
        return _cut_sweep_chunk_parallel(
            elR, wboR, elT, wboT, cfg,
            min(int(n_workers), len(normalized_cuts)),
            core_R, normalized_cuts)
    return _cut_sweep_chunk_serial(
        elR, wboR, elT, wboT, cfg, core_R, normalized_cuts)


def merge_cut_sweep_pools(pools):
    """Merge partial cut-sweep pools produced by chunk tasks."""
    merged = {}
    for pool in pools:
        for sig, info in dict(pool or {}).items():
            cuts = frozenset(info.get('cuts', ()))
            no_cut = bool(info.get('has_no_cut', False))
            entry = merged.get(sig)
            if entry is None:
                merged[sig] = {
                    'mapping': dict(info['mapping']),
                    'cuts': cuts,
                    'has_no_cut': no_cut,
                    'dedup_count': int(info.get('dedup_count', 1)),
                }
                continue
            if no_cut and not entry.get('has_no_cut', False):
                entry['mapping'] = dict(info['mapping'])
            entry['cuts'] = frozenset(entry.get('cuts', ())) | cuts
            entry['has_no_cut'] = bool(
                entry.get('has_no_cut', False)
                or no_cut
            )
            entry['dedup_count'] = (
                int(entry.get('dedup_count', 1))
                + int(info.get('dedup_count', 1))
            )
    return merged


def cut_sweep(elR, wboR, elT, wboT, *,
              n_workers=None, core_R=None,
              cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
              dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
              symmetry_wbo_tol=0.2,
              n_seeds=3, max_branches=100,
              chunksize=1,
              symmetry_repair=True,
              symmetry_repair_min_changes=1,
              symmetry_repair_max_evals=20000):
    """Enumerate mechanism classes via no-cut plus one-edge R cuts.

    The returned pool maps a symmetry-canonical signature to:

    - `mapping`: representative symmetry-aware witness
    - `cuts`: set of R-edge cuts that led to that signature
    - `dedup_count`: number of witnesses collapsed into the signature

    With `core_R=None`, signatures are R-P mechanism signatures.  With
    `core_R` supplied, signatures are exact core mappings; this is useful for
    mechanism-local TS/IG scoring, but R-P mechanism discovery is the primary
    use.
    """
    cfg = _cut_sweep_cfg(
        cut_floor=cut_floor,
        graph_floor=graph_floor,
        iso_tol=iso_tol,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
        symmetry_wbo_tol=symmetry_wbo_tol,
        n_seeds=n_seeds,
        max_branches=max_branches,
        chunksize=chunksize,
        symmetry_repair=symmetry_repair,
        symmetry_repair_min_changes=symmetry_repair_min_changes,
        symmetry_repair_max_evals=symmetry_repair_max_evals,
        n_atoms=len(elR),
    )
    core_R = tuple(sorted(set(core_R or ())))
    if not n_workers or n_workers <= 1:
        return _cut_sweep_serial(elR, wboR, elT, wboT, cfg, core_R)
    return _cut_sweep_parallel(elR, wboR, elT, wboT, cfg,
                               int(n_workers), core_R)


def select_min_mechanisms(pool):
    """Keep only signatures with the fewest broken+formed orbit events."""
    if not pool:
        return {}
    best = min(len(sig[0]) + len(sig[1]) for sig in pool)
    return {
        sig: info for sig, info in pool.items()
        if len(sig[0]) + len(sig[1]) == best
    }

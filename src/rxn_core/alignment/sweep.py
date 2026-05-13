"""Sweep-cut mechanism discovery for WBO graph alignment.

R-P sweep cut is part of the core alignment algorithm: mechanism discovery
tries the no-cut graph and every one-edge R cut above a WBO floor, then dedups
the resulting witnesses by symmetry-canonical broken/formed bond signatures.
"""
from __future__ import annotations

import multiprocessing as mp
import signal

from ..frag import build_graph, classify_bonds, expand_mapping
from ..matcher import _nauty_orbits
from .branch import _generate_seed_orders, find_islands, symmetry_repair_mapping


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


def _mechanism_signature(mapping, wboR, wboT, r_orbits, p_orbits):
    """Symmetry-canonical mechanism key for R-P discovery."""
    broken, formed, _, _ = classify_bonds(mapping, wboR, wboT)
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
    entry = pool.get(sig)
    if entry is None:
        pool[sig] = {
            'mapping': dict(mapping),
            'cuts': cuts,
            'dedup_count': 1,
        }
    else:
        entry['cuts'] = entry['cuts'] | cuts
        entry['dedup_count'] = entry.get('dedup_count', 1) + 1


def _run_find_islands_limited(g_R, g_P, order, core_R, cfg):
    stop_on_core = bool(core_R)
    unit_timeout = float(cfg['unit_timeout'])
    kwargs = dict(
        iso_tol=float(cfg['iso_tol']),
        max_branches=int(cfg['max_branches']),
        core_R=core_R,
        stop_when_core_mapped=stop_on_core,
    )
    if unit_timeout <= 0 or not hasattr(signal, "SIGALRM"):
        return find_islands(g_R, g_P, list(order), **kwargs)

    def _raise_timeout(signum, frame):
        raise TimeoutError("cut_sweep work unit timed out")

    old_handler = signal.getsignal(signal.SIGALRM)
    signal.signal(signal.SIGALRM, _raise_timeout)
    signal.setitimer(signal.ITIMER_REAL, unit_timeout)
    try:
        return find_islands(g_R, g_P, list(order), **kwargs)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
        signal.signal(signal.SIGALRM, old_handler)


def _score_branch_mapping(mapping, g_R, g_P, wboR, wboT,
                          g_R_full, p_orbits, r_orbits, core_R, cfg):
    if core_R:
        if not all(r in mapping for r in core_R):
            return None
        return _core_mapping_key(mapping, core_R), mapping

    if len(mapping) < int(cfg['n_atoms']) - 2:
        return None
    if cfg['symmetry_repair']:
        mapping = symmetry_repair_mapping(
            mapping, wboR, wboT, g_R_full, g_P, p_orbits,
            min_changes=int(cfg['symmetry_repair_min_changes']),
            max_evals=int(cfg['symmetry_repair_max_evals']),
        )
    sig = _mechanism_signature(mapping, wboR, wboT, r_orbits, p_orbits)
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
    _WORKER['p_orbits'] = _nauty_orbits(_WORKER['g_P'], wbo_tol=0.2)
    _WORKER['r_orbits'] = _nauty_orbits(_WORKER['g_R_full'], wbo_tol=0.2)


def _cs_wrun(args):
    cut, order, core_R = args
    cfg = _WORKER['cfg']
    graph_floor = float(cfg['graph_floor'])
    g_R = build_graph(_WORKER['elR'], _WORKER['wboR'],
                      bond_cut=graph_floor)
    for i, j in cut:
        if g_R.has_edge(i, j):
            g_R.remove_edge(i, j)
    try:
        branches = _run_find_islands_limited(
            g_R, _WORKER['g_P'], order, core_R, cfg)
    except Exception:
        return []

    out = []
    for branch in branches:
        mapping = expand_mapping(dict(branch.mapping), g_R, _WORKER['g_P'])
        scored = _score_branch_mapping(
            mapping, g_R, _WORKER['g_P'],
            _WORKER['wboR'], _WORKER['wboT'],
            _WORKER['g_R_full'], _WORKER['p_orbits'],
            _WORKER['r_orbits'], core_R, cfg)
        if scored is None:
            continue
        sig, mapping = scored
        out.append((sig, tuple(sorted(mapping.items())), cut))
    return out


def _cut_sweep_serial(elR, wboR, elT, wboT, cfg, core_R):
    graph_floor = float(cfg['graph_floor'])
    g_P = build_graph(elT, wboT, bond_cut=graph_floor)
    g_R_full = build_graph(elR, wboR, bond_cut=graph_floor)
    p_orbits = _nauty_orbits(g_P, wbo_tol=0.2)
    r_orbits = _nauty_orbits(g_R_full, wbo_tol=0.2)
    strong = _strong_edges(wboR, float(cfg['cut_floor']))
    pool = {}

    def run(cuts):
        g_R = build_graph(elR, wboR, bond_cut=graph_floor)
        for i, j in cuts:
            if g_R.has_edge(i, j):
                g_R.remove_edge(i, j)
        orders = _generate_seed_orders(g_R, n_trials=int(cfg['n_seeds']))
        for order in orders:
            try:
                branches = _run_find_islands_limited(
                    g_R, g_P, order, core_R, cfg)
            except Exception:
                continue
            for branch in branches:
                mapping = expand_mapping(dict(branch.mapping), g_R, g_P)
                scored = _score_branch_mapping(
                    mapping, g_R, g_P, wboR, wboT, g_R_full,
                    p_orbits, r_orbits, core_R, cfg)
                if scored is None:
                    continue
                sig, mapping = scored
                _pool_add(pool, sig, mapping, cuts)

    run(())
    for edge in strong:
        run((edge,))
    return pool


def _cut_sweep_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R):
    graph_floor = float(cfg['graph_floor'])
    g_R = build_graph(elR, wboR, bond_cut=graph_floor)
    orders = [tuple(order) for order in _generate_seed_orders(
        g_R, n_trials=int(cfg['n_seeds']))]
    cuts = [()] + [((i, j),) for i, j in _strong_edges(
        wboR, float(cfg['cut_floor']))]
    work = [(cut, order, core_R) for order in orders for cut in cuts]
    pool = {}
    with mp.Pool(n_workers, initializer=_cs_winit,
                 initargs=(elR, wboR, elT, wboT, cfg)) as proc_pool:
        for results in proc_pool.imap_unordered(
                _cs_wrun, work, chunksize=max(1, int(cfg['chunksize']))):
            for sig, mapping_items, cut in results:
                _pool_add(pool, sig, dict(mapping_items), cut)
    return pool


def cut_sweep(elR, wboR, elT, wboT, *,
              n_workers=None, core_R=None,
              cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
              n_seeds=3, max_branches=1_000_000,
              chunksize=1, unit_timeout=0.0,
              symmetry_repair=True,
              symmetry_repair_min_changes=5,
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
    cfg = {
        'cut_floor': float(cut_floor),
        'graph_floor': float(graph_floor),
        'iso_tol': float(iso_tol),
        'n_seeds': int(n_seeds),
        'max_branches': int(max_branches),
        'chunksize': int(chunksize),
        'unit_timeout': float(unit_timeout),
        'symmetry_repair': bool(symmetry_repair),
        'symmetry_repair_min_changes': int(symmetry_repair_min_changes),
        'symmetry_repair_max_evals': int(symmetry_repair_max_evals),
        'n_atoms': len(elR),
    }
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

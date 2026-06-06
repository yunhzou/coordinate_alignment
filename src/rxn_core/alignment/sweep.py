"""Sweep-cut mechanism discovery for WBO graph alignment.

R-P sweep cut is part of the core alignment algorithm: mechanism discovery
tries the no-cut graph and every one-edge R cut above a WBO floor, then dedups
the resulting witnesses by symmetry-canonical broken/formed bond signatures.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import time
from pathlib import Path

from ..frag import build_graph, classify_bonds, expand_mapping
from ..matcher import _SymBlock, _nauty_orbits, _sym_block_assignment_expr
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


def _core_branch_record(branch, core_R):
    """Return the core-restricted compressed state for one aligned branch.

    This is intentionally not an exact core-map expansion.  It keeps the
    branch's set-to-set symmetry pools so TS/IG scoring can merge endpoint
    support before materializing automorphism variants.
    """
    core_R = tuple(sorted(set(int(r) for r in core_R)))
    core_set = set(core_R)
    mapping = {
        int(r): int(branch.mapping[r])
        for r in core_R
        if r in branch.mapping
    }
    if len(mapping) != len(core_R):
        return None

    blocks = {}
    for fragment in getattr(branch, 'symmetry_fragments', ()):
        symmetry = fragment.get('symmetry') or {}
        for block in symmetry.get('blocks') or ():
            r_atoms = tuple(sorted(
                int(r) for r in block.get('r_atoms', ())
                if int(r) in core_set
            ))
            if not r_atoms:
                continue
            p_atoms = tuple(sorted(int(p) for p in block.get('p_atoms', ())))
            if len(p_atoms) < len(r_atoms):
                continue
            blocks[(r_atoms, p_atoms)] = {
                'r_atoms': list(r_atoms),
                'p_atoms': list(p_atoms),
            }

    return {
        'mapping': mapping,
        'blocks': [blocks[key] for key in sorted(blocks)],
        'dedup_count': 1,
    }


def _core_branch_record_key(record, core_R):
    core_R = tuple(sorted(set(int(r) for r in core_R)))
    mapping = {int(r): int(t) for r, t in record.get('mapping', {}).items()}
    blocks = []
    block_r = set()
    for block in record.get('blocks', ()):
        r_atoms = tuple(sorted(int(r) for r in block.get('r_atoms', ())))
        p_atoms = tuple(sorted(int(p) for p in block.get('p_atoms', ())))
        if not r_atoms:
            continue
        blocks.append((r_atoms, p_atoms))
        block_r.update(r_atoms)
    fixed = tuple(
        (int(r), int(mapping[r]))
        for r in core_R
        if r in mapping and r not in block_r
    )
    return fixed, tuple(sorted(blocks))


def _core_mapping_variants(branch, core_R, max_variants, *,
                           g_P=None, p_orbits=None):
    """Enumerate exact core maps represented by one compressed branch.

    Fragment growth stores one concrete witness plus symmetry blocks.  For
    TS/IG verification the spectator witness is irrelevant, but core atoms
    inside those blocks must be shuffled before scoring.  This function does
    not impose any core-core WBO validity rule; normal island growth already
    decided what can be aligned, including separate islands and deferred
    nonmatching edges.
    """
    core_R = tuple(sorted(set(int(r) for r in core_R)))
    core_set = set(core_R)
    base = {
        int(r): int(branch.mapping[r])
        for r in core_R
        if r in branch.mapping
    }
    if len(base) != len(core_R):
        return []

    options = {r: {p} for r, p in base.items()}
    for fragment in getattr(branch, 'symmetry_fragments', ()):
        symmetry = fragment.get('symmetry') or {}
        for block in symmetry.get('blocks') or ():
            r_atoms = tuple(sorted(
                int(r) for r in block.get('r_atoms', ())
                if int(r) in core_set
            ))
            if not r_atoms:
                continue
            p_atoms = tuple(sorted(int(p) for p in block.get('p_atoms', ())))
            if len(p_atoms) < len(r_atoms):
                continue
            for r in r_atoms:
                options.setdefault(r, set()).update(p_atoms)

    order = sorted(core_R, key=lambda r: (len(options.get(r, ())), r))
    out = {}

    def rec(i, state, used):
        if len(out) > max_variants:
            raise BranchLimitExceeded(
                max_variants,
                branch_count=len(out),
                stage='core_symmetry_expansion',
            )
        if i == len(order):
            mapping = {int(r): int(state[r]) for r in core_R}
            out[_core_mapping_key(mapping, core_R)] = mapping
            if len(out) > max_variants:
                raise BranchLimitExceeded(
                    max_variants,
                    branch_count=len(out),
                    stage='core_symmetry_expansion',
                )
            return
        r = order[i]
        for p in sorted(options.get(r, {base[r]})):
            if p in used:
                continue
            state[r] = p
            used.add(p)
            rec(i + 1, state, used)
            used.remove(p)
            del state[r]

    rec(0, {}, set())
    return list(out.values())


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


def _branch_symmetry_record(branch):
    fragments = []
    blocks = []
    p_to_r = {
        int(p): int(r)
        for r, p in getattr(branch, 'mapping', {}).items()
    }
    for frag_index, fragment in enumerate(
            getattr(branch, 'symmetry_fragments', ())):
        record = {
            'fragment_index': int(frag_index),
            'island_idx': int(fragment.get('island_idx', frag_index)),
            'fragment': [int(r) for r in fragment.get('fragment', ())],
            'deferred_edges': [
                [int(a), int(b)]
                for a, b in fragment.get('deferred_edges', ())
            ],
            'symmetry': fragment.get('symmetry') or {},
        }
        fragments.append(record)
        symmetry = record['symmetry']
        for block_index, block in enumerate(symmetry.get('blocks') or ()):
            r_atoms = {
                int(r) for r in block.get('r_atoms', ())
            }
            p_atoms = [int(p) for p in block.get('p_atoms', ())]
            for p in p_atoms:
                if p in p_to_r:
                    r_atoms.add(p_to_r[p])
            r_atoms = sorted(r_atoms)
            if len(p_atoms) <= 1 and len(r_atoms) <= 1:
                continue
            normalized_block = _SymBlock(tuple(r_atoms), tuple(p_atoms),
                                         extendable=False)
            source = block.get('source')
            assignments = (
                block.get('assignments')
                if source in {'island_automorph', 'alternate_witness',
                              'sym_block', 'interbranch'}
                else _sym_block_assignment_expr(normalized_block)
            )
            item = {
                'fragment_index': int(frag_index),
                'block_index': int(block_index),
                'island_idx': record['island_idx'],
                'r_atoms': r_atoms,
                'p_atoms': p_atoms,
                'extendable': bool(block.get('extendable', False)),
                'open': bool(normalized_block.open),
                'assignments': assignments,
            }
            if source:
                item['source'] = source
            blocks.append(item)
        witness = {
            int(r): int(p)
            for r, p in dict(symmetry.get('witness') or {}).items()
        }
        for alt_index, alt in enumerate(symmetry.get('alternates') or ()):
            alt_witness = {
                int(r): int(p)
                for r, p in dict(alt.get('witness') or {}).items()
            }
            changed_r = sorted(
                r for r, p in alt_witness.items()
                if r in witness and int(witness[r]) != int(p)
            )
            p_atoms = sorted({
                int(witness[r])
                for r in changed_r
            } | {
                int(alt_witness[r])
                for r in changed_r
            })
            if len(p_atoms) <= 1 and len(changed_r) <= 1:
                continue
            blocks.append({
                'fragment_index': int(frag_index),
                'block_index': f"alt:{int(alt_index)}",
                'island_idx': record['island_idx'],
                'r_atoms': changed_r,
                'p_atoms': p_atoms,
                'extendable': False,
                'open': False,
                'assignments': alt.get('multiplicity'),
                'source': 'alternate_witness',
            })
    return {
        'rule': 'branch_symmetry_blocks',
        'fragments': fragments,
        'blocks': blocks,
    }


def _mapping_key(mapping):
    return tuple(
        (int(r), int(p))
        for r, p in sorted(dict(mapping).items())
    )


def _cut_record(cuts):
    return [
        [int(a), int(b)]
        for a, b in sorted(tuple(tuple(pair) for pair in cuts))
    ]


def _normalize_cut_record(cut):
    out = []
    for item in cut or ():
        if isinstance(item, str):
            continue
        try:
            a, b = item
        except (TypeError, ValueError):
            continue
        out.append([int(a), int(b)])
    return out


def _dedup_witness(mapping, cuts, branch_symmetry=None):
    return {
        'mapping': {int(r): int(p) for r, p in dict(mapping).items()},
        'cut': _cut_record(cuts),
        'local_symmetry': branch_symmetry or {
            'rule': 'branch_symmetry_blocks',
            'fragments': [],
            'blocks': [],
        },
    }


def _unique_witnesses(witnesses):
    out = []
    seen = set()
    for witness in witnesses or ():
        mapping = {
            int(r): int(p)
            for r, p in dict(witness.get('mapping') or {}).items()
        }
        if not mapping:
            continue
        key = _mapping_key(mapping)
        if key in seen:
            continue
        seen.add(key)
        out.append({
            'mapping': mapping,
            'cut': _normalize_cut_record(witness.get('cut', ())),
            'local_symmetry': witness.get('local_symmetry') or {
                'rule': 'branch_symmetry_blocks',
                'fragments': [],
                'blocks': [],
            },
        })
    return out


def _color_groups_from_blocks(blocks):
    groups = []
    by_key = {}
    for block in blocks:
        r_atoms = sorted(int(r) for r in block.get('r_atoms', ()))
        p_atoms = sorted(int(p) for p in block.get('p_atoms', ()))
        if len(r_atoms) <= 1 and len(p_atoms) <= 1:
            continue
        source = block.get('source') or 'sym_block'
        key = (tuple(r_atoms), tuple(p_atoms))
        if key in by_key:
            by_key[key]['sources'].add(source)
            continue
        group = {
            'r_atoms': r_atoms,
            'p_atoms': p_atoms,
            'sources': {source},
        }
        by_key[key] = group
        groups.append(group)
    for group in groups:
        group['sources'] = sorted(group['sources'])
    return groups


def combine_branch_symmetry_witnesses(witnesses):
    """Build the principled branch-degeneracy record for one mechanism.

    Displayed degeneracy is only the ambiguity carried by island/SymCand state:
    local `_SymCand` blocks/alternates and same-island candidate automorph
    blocks recorded during island growth. Final dedup witnesses remain as
    provenance but do not create extra coloring.
    """
    witnesses = _unique_witnesses(witnesses)
    fragments = []
    blocks = []
    block_keys = set()

    for witness_index, witness in enumerate(witnesses):
        local = witness.get('local_symmetry') or {}
        for fragment in local.get('fragments') or ():
            record = dict(fragment)
            record['witness_index'] = int(witness_index)
            fragments.append(record)
        for block in local.get('blocks') or ():
            r_atoms = tuple(sorted(int(r) for r in block.get('r_atoms', ())))
            p_atoms = tuple(sorted(int(p) for p in block.get('p_atoms', ())))
            source = block.get('source') or 'sym_block'
            key = (source, r_atoms, p_atoms)
            if key in block_keys:
                continue
            block_keys.add(key)
            record = dict(block)
            record['source'] = source
            record['r_atoms'] = list(r_atoms)
            record['p_atoms'] = list(p_atoms)
            record['witness_index'] = int(witness_index)
            blocks.append(record)

    return {
        'rule': 'mechanism_dedup_branch_symmetry',
        'fragments': fragments,
        'blocks': blocks,
        'color_groups': _color_groups_from_blocks(blocks),
        'witnesses': witnesses,
        'dedup_witness_count': len(witnesses),
    }


def _refresh_entry_branch_symmetry(entry):
    entry['dedup_witnesses'] = _unique_witnesses(entry.get('dedup_witnesses'))
    entry['branch_symmetry'] = combine_branch_symmetry_witnesses(
        entry['dedup_witnesses'])


def _pool_add(pool, sig, mapping, cuts, branch_symmetry=None):
    cuts = frozenset(cuts)
    no_cut = not cuts
    witness = _dedup_witness(mapping, cuts, branch_symmetry)
    entry = pool.get(sig)
    if entry is None:
        pool[sig] = {
            'mapping': dict(mapping),
            'cuts': cuts,
            'has_no_cut': bool(no_cut),
            'dedup_count': 1,
            'dedup_witnesses': [witness],
        }
        _refresh_entry_branch_symmetry(pool[sig])
    else:
        if no_cut and not entry.get('has_no_cut', False):
            entry['mapping'] = dict(mapping)
        entry['has_no_cut'] = bool(entry.get('has_no_cut', False) or no_cut)
        entry['cuts'] = entry['cuts'] | cuts
        entry['dedup_count'] = entry.get('dedup_count', 1) + 1
        entry.setdefault('dedup_witnesses', []).append(witness)
        _refresh_entry_branch_symmetry(entry)


def _anchor_mapping_ok(mapping, anchor_map):
    return all(
        int(mapping.get(int(r), -1)) == int(p)
        for r, p in dict(anchor_map or {}).items()
    )


def _run_find_islands_limited(g_R, g_P, order, core_R, cfg, *,
                              p_orbits=None, r_orbits=None,
                              profile=None):
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
        profile=profile,
        anchor_map=cfg.get('anchor_map'),
    )


def _score_branch_mapping(mapping, g_R, g_P, wboR, wboT,
                          g_R_full, p_orbits, r_orbits, core_R, cfg,
                          elR=None, elT=None,
                          return_repair_stats=False):
    anchor_map = cfg.get('anchor_map') or {}
    if core_R:
        if not all(r in mapping for r in core_R):
            return None
        if not _anchor_mapping_ok(mapping, anchor_map):
            return None
        scored = (_core_mapping_key(mapping, core_R), mapping)
        return (*scored, None) if return_repair_stats else scored

    if len(mapping) < int(cfg['n_atoms']) - 2:
        return None
    if not _anchor_mapping_ok(mapping, anchor_map):
        return None
    repair_stats = None
    if cfg['symmetry_repair']:
        base_mapping = dict(mapping)
        if return_repair_stats:
            mapping, repair_stats = symmetry_repair_mapping(
                mapping, wboR, wboT, g_R_full, g_P, p_orbits,
                dwbo_threshold=float(cfg['dwbo_threshold']),
                metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'),
                min_changes=int(cfg['symmetry_repair_min_changes']),
                max_evals=int(cfg['symmetry_repair_max_evals']),
                return_stats=True,
            )
        else:
            mapping = symmetry_repair_mapping(
                mapping, wboR, wboT, g_R_full, g_P, p_orbits,
                dwbo_threshold=float(cfg['dwbo_threshold']),
                metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'),
                min_changes=int(cfg['symmetry_repair_min_changes']),
                max_evals=int(cfg['symmetry_repair_max_evals']),
            )
        if not _anchor_mapping_ok(mapping, anchor_map):
            mapping = base_mapping
            if repair_stats is not None:
                repair_stats = dict(repair_stats)
                repair_stats['anchor_reverted'] = True
    sig = _mechanism_signature(
        mapping, wboR, wboT, r_orbits, p_orbits,
        dwbo_threshold=float(cfg['dwbo_threshold']),
        elements_R=elR,
        elements_P=elT,
        metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'))
    scored = (sig, mapping)
    return (*scored, repair_stats) if return_repair_stats else scored


def _cut_json(cut):
    return [list(map(int, pair)) for pair in cut]


def _repair_trace_stats(stats):
    if not stats:
        return None
    groups = stats.get('groups') or []
    return {
        'enabled': bool(stats.get('enabled', False)),
        'base_changes': stats.get('base_changes'),
        'best_changes': stats.get('best_changes'),
        'repaired': bool(stats.get('repaired', False)),
        'evaluated': int(stats.get('evaluated', 0) or 0),
        'capped': bool(stats.get('capped', False)),
        'n_groups': len(groups),
        'group_sizes': [int(group.get('size', 0) or 0)
                        for group in groups],
    }


def _growth_trace_summary(profile):
    if not profile:
        return {
            'calls': 0,
            'elapsed_sec': 0.0,
            'extend_elapsed_sec': 0.0,
            'heap_pops': 0,
            'extend_calls': 0,
            'commits': 0,
            'deferred': 0,
            'max_cands_before': 0,
            'max_cands_after': 0,
            'max_heap_len': 0,
            'max_fragment_size': 0,
            'slowest': [],
        }
    slowest = sorted(
        profile,
        key=lambda item: float(item.get('elapsed_sec', 0.0) or 0.0),
        reverse=True,
    )[:5]
    return {
        'calls': len(profile),
        'elapsed_sec': sum(float(item.get('elapsed_sec', 0.0) or 0.0)
                           for item in profile),
        'extend_elapsed_sec': sum(float(item.get('extend_elapsed_sec', 0.0) or 0.0)
                                  for item in profile),
        'heap_pops': sum(int(item.get('heap_pops', 0) or 0)
                         for item in profile),
        'extend_calls': sum(int(item.get('extend_calls', 0) or 0)
                            for item in profile),
        'commits': sum(int(item.get('commits', 0) or 0)
                       for item in profile),
        'deferred': sum(int(item.get('deferred', 0) or 0)
                        for item in profile),
        'max_cands_before': max(int(item.get('max_cands_before', 0) or 0)
                                for item in profile),
        'max_cands_after': max(int(item.get('max_cands_after', 0) or 0)
                               for item in profile),
        'max_heap_len': max(int(item.get('max_heap_len', 0) or 0)
                            for item in profile),
        'max_fragment_size': max(int(item.get('max_fragment_size', 0) or 0)
                                 for item in profile),
        'slowest': [
            {
                'seed': int(item.get('seed')),
                'pass': int(item.get('pass', 0) or 0),
                'branch_index': int(item.get('branch_index', 0) or 0),
                'mapped_before': int(item.get('mapped_before', 0) or 0),
                'result': item.get('result'),
                'elapsed_sec': float(item.get('elapsed_sec', 0.0) or 0.0),
                'extend_elapsed_sec': float(
                    item.get('extend_elapsed_sec', 0.0) or 0.0),
                'heap_pops': int(item.get('heap_pops', 0) or 0),
                'extend_calls': int(item.get('extend_calls', 0) or 0),
                'commits': int(item.get('commits', 0) or 0),
                'deferred': int(item.get('deferred', 0) or 0),
                'max_cands_before': int(item.get('max_cands_before', 0) or 0),
                'max_cands_after': int(item.get('max_cands_after', 0) or 0),
                'slowest_extend': item.get('slowest_extend'),
            }
            for item in slowest
        ],
    }


def _emit_trace(trace_path, events):
    if not trace_path or not events:
        return
    path = Path(trace_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a', encoding='utf-8') as handle:
        for event in events:
            handle.write(json.dumps(event, sort_keys=True) + '\n')


def _run_cut_work(elR, wboR, elT, wboT, cfg, cut, orders, core_R,
                  g_P, g_R_full, p_orbits, r_orbits_full,
                  *, return_trace=False):
    cut = tuple(tuple(int(v) for v in pair) for pair in cut)
    events = []
    out = []
    cut_t0 = time.perf_counter()

    graph_t0 = time.perf_counter()
    graph_floor = float(cfg['graph_floor'])
    g_R = build_graph(elR, wboR, bond_cut=graph_floor)
    for i, j in cut:
        if g_R.has_edge(i, j):
            g_R.remove_edge(i, j)
    r_orbits_cut = _nauty_orbits(
        g_R, wbo_tol=float(cfg['symmetry_wbo_tol']))
    graph_elapsed = time.perf_counter() - graph_t0
    if orders is None:
        orders = _generate_seed_orders(
            g_R, n_trials=int(cfg['n_seeds']))
    else:
        orders = list(orders)
    if return_trace:
        events.append({
            'event': 'cut_start',
            'cut': _cut_json(cut),
            'n_orders': len(orders),
            'max_branches': int(cfg['max_branches']),
            'symmetry_repair': bool(cfg['symmetry_repair']),
            'symmetry_repair_max_evals': int(
                cfg['symmetry_repair_max_evals']),
        })

    cut_status = 'completed'
    total_search_elapsed = 0.0
    total_score_elapsed = 0.0
    total_expand_elapsed = 0.0
    total_branches = 0
    total_accepted = 0
    total_repair_evals = 0
    total_repair_capped = 0
    try:
        for order_index, order in enumerate(orders):
            seed_t0 = time.perf_counter()
            seed_growth_profile = [] if return_trace else None
            try:
                branches = _run_find_islands_limited(
                    g_R, g_P, order, core_R, cfg,
                    p_orbits=p_orbits,
                    r_orbits=r_orbits_cut,
                    profile=seed_growth_profile)
            except BranchLimitExceeded as exc:
                cut_status = 'branch_cap'
                if return_trace:
                    for item in seed_growth_profile:
                        growth_event = dict(item)
                        growth_event['event'] = 'growth_call'
                        growth_event['cut'] = _cut_json(cut)
                        growth_event['seed_index'] = int(order_index)
                        events.append(growth_event)
                    events.append({
                        'event': 'seed_branch_cap',
                        'cut': _cut_json(cut),
                        'seed_index': int(order_index),
                        'elapsed_sec': time.perf_counter() - seed_t0,
                        'branch_count': exc.branch_count,
                        'max_branches': exc.max_branches,
                        'stage': exc.stage,
                        'growth': _growth_trace_summary(seed_growth_profile),
                    })
                raise
            except Exception as exc:
                if return_trace:
                    for item in seed_growth_profile:
                        growth_event = dict(item)
                        growth_event['event'] = 'growth_call'
                        growth_event['cut'] = _cut_json(cut)
                        growth_event['seed_index'] = int(order_index)
                        events.append(growth_event)
                    events.append({
                        'event': 'seed_error',
                        'cut': _cut_json(cut),
                        'seed_index': int(order_index),
                        'elapsed_sec': time.perf_counter() - seed_t0,
                        'error_type': type(exc).__name__,
                        'error': str(exc),
                        'growth': _growth_trace_summary(seed_growth_profile),
                    })
                continue

            search_elapsed = time.perf_counter() - seed_t0
            growth_summary = _growth_trace_summary(seed_growth_profile)
            total_search_elapsed += search_elapsed
            seed_score_elapsed = 0.0
            seed_expand_elapsed = 0.0
            seed_accepted = 0
            seed_repair_evals = 0
            seed_repair_capped = 0
            max_mapped = 0
            for branch_index, branch in enumerate(branches):
                branch_t0 = time.perf_counter()
                expand_elapsed = 0.0
                score_t0 = time.perf_counter()
                if core_R:
                    scored_items = []
                    for core_map in _core_mapping_variants(
                            branch, core_R, int(cfg['max_branches']),
                            g_P=g_P, p_orbits=p_orbits):
                        scored_items.append((
                            _core_mapping_key(core_map, core_R),
                            core_map,
                            None,
                        ))
                    mapping_for_stats = dict(branch.mapping)
                else:
                    expand_t0 = time.perf_counter()
                    mapping = expand_mapping(dict(branch.mapping), g_R, g_P)
                    expand_elapsed = time.perf_counter() - expand_t0
                    scored = _score_branch_mapping(
                        mapping, g_R, g_P, wboR, wboT, g_R_full,
                        p_orbits, r_orbits_full, core_R, cfg, elR, elT,
                        return_repair_stats=return_trace)
                    if scored is None:
                        scored_items = []
                    elif return_trace:
                        sig, repaired_mapping, repair_stats = scored
                        scored_items = [(sig, repaired_mapping, repair_stats)]
                    else:
                        sig, repaired_mapping = scored
                        scored_items = [(sig, repaired_mapping, None)]
                    mapping_for_stats = mapping
                score_elapsed = time.perf_counter() - score_t0
                branch_elapsed = time.perf_counter() - branch_t0
                seed_expand_elapsed += expand_elapsed
                seed_score_elapsed += score_elapsed
                max_mapped = max(max_mapped, len(mapping_for_stats))
                accepted = bool(scored_items)
                for sig, accepted_mapping, _repair_stats in scored_items:
                    out.append((
                        sig,
                        tuple(sorted(accepted_mapping.items())),
                        cut,
                        _branch_symmetry_record(branch),
                    ))
                    seed_accepted += 1
                repair_stats = scored_items[0][2] if scored_items else None
                repair_summary = _repair_trace_stats(repair_stats)
                if repair_summary:
                    seed_repair_evals += repair_summary['evaluated']
                    seed_repair_capped += int(repair_summary['capped'])
                if return_trace:
                    events.append({
                        'event': 'branch',
                        'cut': _cut_json(cut),
                        'seed_index': int(order_index),
                        'branch_index': int(branch_index),
                        'elapsed_sec': branch_elapsed,
                        'expand_elapsed_sec': expand_elapsed,
                        'score_elapsed_sec': score_elapsed,
                        'mapped_atoms': len(mapping_for_stats),
                        'accepted': bool(accepted),
                        'accepted_core_variants': len(scored_items),
                        'repair': repair_summary,
                    })
            n_branches = len(branches)
            total_branches += n_branches
            total_accepted += seed_accepted
            total_expand_elapsed += seed_expand_elapsed
            total_score_elapsed += seed_score_elapsed
            total_repair_evals += seed_repair_evals
            total_repair_capped += seed_repair_capped
            if return_trace:
                for item in seed_growth_profile:
                    growth_event = dict(item)
                    growth_event['event'] = 'growth_call'
                    growth_event['cut'] = _cut_json(cut)
                    growth_event['seed_index'] = int(order_index)
                    events.append(growth_event)
                events.append({
                    'event': 'seed_end',
                    'cut': _cut_json(cut),
                    'seed_index': int(order_index),
                    'search_elapsed_sec': search_elapsed,
                    'expand_elapsed_sec': seed_expand_elapsed,
                    'score_elapsed_sec': seed_score_elapsed,
                    'branches': n_branches,
                    'accepted': seed_accepted,
                    'max_mapped_atoms': max_mapped,
                    'repair_evals': seed_repair_evals,
                    'repair_capped_count': seed_repair_capped,
                    'growth': growth_summary,
                })
    except BranchLimitExceeded:
        # Kill the whole cut.  Partial seed-order witnesses from this cut are
        # deliberately discarded because the cut has entered a pathological
        # outer-branch multiplication regime.
        cut_status = 'branch_cap'
        out = []

    elapsed = time.perf_counter() - cut_t0
    if return_trace:
        events.append({
            'event': 'cut_end',
            'cut': _cut_json(cut),
            'status': cut_status,
            'elapsed_sec': elapsed,
            'graph_elapsed_sec': graph_elapsed,
            'search_elapsed_sec': total_search_elapsed,
            'expand_elapsed_sec': total_expand_elapsed,
            'score_elapsed_sec': total_score_elapsed,
            'branches': total_branches,
            'accepted': total_accepted,
            'repair_evals': total_repair_evals,
            'repair_capped_count': total_repair_capped,
            'hits': len(out),
            'avg_branch_elapsed_sec': (
                (total_expand_elapsed + total_score_elapsed) / total_branches
                if total_branches else 0.0
            ),
        })
    return out, events


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
    cut, orders, core_R, trace_enabled = args
    cfg = _WORKER['cfg']
    out, events = _run_cut_work(
        _WORKER['elR'], _WORKER['wboR'],
        _WORKER['elT'], _WORKER['wboT'],
        cfg, cut, orders, core_R,
        _WORKER['g_P'], _WORKER['g_R_full'],
        _WORKER['p_orbits'], _WORKER['r_orbits'],
        return_trace=trace_enabled)
    return {'results': out, 'events': events}


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
                   n_atoms=0,
                   anchor_map=None):
    anchor_map = {
        int(r): int(p)
        for r, p in dict(anchor_map or {}).items()
    }
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
        'anchor_map': anchor_map,
    }


def _cut_sweep_chunk_serial(elR, wboR, elT, wboT, cfg, core_R, cuts,
                            trace_path=None):
    graph_floor = float(cfg['graph_floor'])
    g_P = build_graph(elT, wboT, bond_cut=graph_floor)
    g_R_full = build_graph(elR, wboR, bond_cut=graph_floor)
    symmetry_wbo_tol = float(cfg['symmetry_wbo_tol'])
    p_orbits = _nauty_orbits(g_P, wbo_tol=symmetry_wbo_tol)
    r_orbits = _nauty_orbits(g_R_full, wbo_tol=symmetry_wbo_tol)
    pool = {}

    for cut in cuts:
        cut = tuple(tuple(pair) for pair in cut)
        results, events = _run_cut_work(
            elR, wboR, elT, wboT, cfg, cut, None, core_R,
            g_P, g_R_full, p_orbits, r_orbits,
            return_trace=bool(trace_path))
        _emit_trace(trace_path, events)
        for result in results:
            sig, mapping_items, _cut = result[:3]
            branch_symmetry = result[3] if len(result) > 3 else None
            mapping = dict(mapping_items)
            _pool_add(pool, sig, mapping, _cut, branch_symmetry)
    return pool


def _cut_sweep_serial(elR, wboR, elT, wboT, cfg, core_R):
    return _cut_sweep_chunk_serial(
        elR, wboR, elT, wboT, cfg, core_R,
        cut_sweep_items(wboR, cfg['cut_floor']))


def _cut_sweep_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R):
    cuts = cut_sweep_items(wboR, cfg['cut_floor'])
    work = [(cut, None, core_R, False) for cut in cuts]
    pool = {}
    with mp.Pool(n_workers, initializer=_cs_winit,
                 initargs=(elR, wboR, elT, wboT, cfg)) as proc_pool:
        for payload in proc_pool.imap_unordered(
                _cs_wrun, work, chunksize=max(1, int(cfg['chunksize']))):
            for result in payload['results']:
                sig, mapping_items, cut = result[:3]
                branch_symmetry = result[3] if len(result) > 3 else None
                _pool_add(pool, sig, dict(mapping_items), cut, branch_symmetry)
    return pool


def _cut_sweep_chunk_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R,
                              cuts, trace_path=None):
    work = [(cut, None, core_R, bool(trace_path)) for cut in cuts]
    pool = {}
    with mp.Pool(n_workers, initializer=_cs_winit,
                 initargs=(elR, wboR, elT, wboT, cfg)) as proc_pool:
        for payload in proc_pool.imap_unordered(
                _cs_wrun, work, chunksize=max(1, int(cfg['chunksize']))):
            _emit_trace(trace_path, payload.get('events', []))
            for result in payload['results']:
                sig, mapping_items, cut = result[:3]
                branch_symmetry = result[3] if len(result) > 3 else None
                _pool_add(pool, sig, dict(mapping_items), cut, branch_symmetry)
    return pool


def run_cut_sweep_chunk(elR, wboR, elT, wboT, cuts, *,
                        core_R=None,
                        n_workers=None,
                        trace_path=None,
                        cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
                        dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
                        symmetry_wbo_tol=0.2,
                        n_seeds=3, max_branches=100,
                        chunksize=1,
                        symmetry_repair=True,
                        symmetry_repair_min_changes=1,
                        symmetry_repair_max_evals=20000,
                        anchor_map=None):
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
        anchor_map=anchor_map,
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
            core_R, normalized_cuts, trace_path=trace_path)
    return _cut_sweep_chunk_serial(
        elR, wboR, elT, wboT, cfg, core_R, normalized_cuts,
        trace_path=trace_path)


def run_no_cut_core_branch_records(elS, wboS, elT, wboT, core_S, *,
                                   graph_floor=0.2, iso_tol=1.0,
                                   dwbo_threshold=0.5,
                                   metal_dwbo_threshold=0.3,
                                   symmetry_wbo_tol=0.2,
                                   n_seeds=3, max_branches=20000):
    """Return no-cut endpoint->target branch states for a mechanism core.

    Unlike ``run_cut_sweep_chunk(..., core_R=...)``, this keeps branch-level
    compressed symmetry records.  Callers can merge R-derived and P-derived
    endpoint support first, then expand only the automorphism variants inside
    each merged branch state for TS/IG scoring.
    """
    core_S = tuple(sorted(set(int(r) for r in core_S or ())))
    if not core_S:
        return []

    g_S = build_graph(elS, wboS, bond_cut=float(graph_floor))
    g_T = build_graph(elT, wboT, bond_cut=float(graph_floor))
    p_orbits = _nauty_orbits(g_T, wbo_tol=float(symmetry_wbo_tol))
    r_orbits = _nauty_orbits(g_S, wbo_tol=float(symmetry_wbo_tol))
    orders = _generate_seed_orders(g_S, n_trials=int(n_seeds))
    records = {}

    try:
        for order in orders:
            branches = find_islands(
                g_S, g_T, list(order),
                iso_tol=float(iso_tol),
                max_branches=int(max_branches),
                abort_on_branch_cap=True,
                dwbo_threshold=float(dwbo_threshold),
                metal_dwbo_threshold=metal_dwbo_threshold,
                symmetry_wbo_tol=float(symmetry_wbo_tol),
                core_R=core_S,
                stop_when_core_mapped=True,
                p_orbits=p_orbits,
                r_orbits=r_orbits,
            )
            for branch in branches:
                record = _core_branch_record(branch, core_S)
                if record is None:
                    continue
                key = _core_branch_record_key(record, core_S)
                existing = records.get(key)
                if existing is None:
                    records[key] = record
                else:
                    existing['dedup_count'] = (
                        int(existing.get('dedup_count', 1))
                        + int(record.get('dedup_count', 1))
                    )
    except BranchLimitExceeded:
        return []

    return list(records.values())


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
                    'dedup_witnesses': _unique_witnesses(
                        info.get('dedup_witnesses') or [
                            _dedup_witness(
                                info['mapping'], cuts,
                                info.get('branch_symmetry'))
                        ]),
                }
                _refresh_entry_branch_symmetry(merged[sig])
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
            entry.setdefault('dedup_witnesses', []).extend(
                _unique_witnesses(info.get('dedup_witnesses') or [
                    _dedup_witness(
                        info['mapping'], cuts, info.get('branch_symmetry'))
                ])
            )
            _refresh_entry_branch_symmetry(entry)
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
              symmetry_repair_max_evals=20000,
              anchor_map=None):
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
        anchor_map=anchor_map,
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

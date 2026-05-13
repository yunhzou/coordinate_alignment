"""Live priority-queue fragment growth loop."""
from __future__ import annotations

import heapq

from ..matcher import (
    _SymBlock,
    _SymCand,
    _boundary_signature,
    _cand_canon_signature,
    _cand_map,
    _cand_possible_p_atoms,
    _extend_sym_cands,
    _group_nodes_by_signature,
    _orbit_id,
    _symmetry_state,
)
from .frontier import _frontier_boundary_edges, _push_edges_from, _set_unique
from .result import _IsoResult
from .trace import (
    cand_possible_values,
    cands_pattern_sample,
    cands_sample,
    heap_snapshot,
    pool_by_frag_atom,
    represented_assignment_expr,
    why_extend_failed,
    why_merge_failed,
)


def grow_island(g_R, g_P, seed, mapping,
                graph_floor=0.2,
                iso_tol=0.5,
                min_lock_size=1,
                max_branches=1_000_000,
                events=None,
                islands_R=None,
                p_orbits=None,
                r_orbits=None,
                prior_deferred_edges=None):
    """
    Grow a fragment from `seed` using priority-queue propagation.

    Returns a list of isos:
      []            -- failed (no initial cands, or fragment too small)
      [single_iso]  -- locked successfully (set-unique or single cand)
      [iso_a, ...]  -- non-set-unique saturation; caller branches

    Optional `events` list receives diagnostic events (seed_start /
    commit / deferred / merge / seed_end) compatible with the
    existing trace_run.HTML viewer.
    """
    record = events is not None
    locked_p_atoms = set(mapping.values())
    if seed in mapping:
        return []
    seed_el = g_R.nodes[seed]['element']
    seed_targets = [v for v in g_P.nodes()
                    if v not in locked_p_atoms and g_P.nodes[v]['element'] == seed_el]
    if p_orbits is not None:
        seed_groups = _group_nodes_by_signature(
            seed_targets, lambda v: (g_P.nodes[v].get('element'),
                                     _orbit_id(p_orbits, v)))
    else:
        seed_groups = [(v,) for v in sorted(seed_targets)]
    cands = []
    for group in seed_groups:
        if len(group) > 1:
            cands.append(_SymCand({seed: group[0]},
                                  (_SymBlock((seed,), group,
                                             extendable=False),)))
        else:
            cands.append(_SymCand({seed: group[0]}))
    if not cands:
        if record:
            events.append({'type': 'seed_start', 'seed': int(seed),
                           'init_cands': 0, 'fragment': [int(seed)],
                           'p_atoms': []})
            events.append({'type': 'seed_end', 'result': 'no_initial_cands',
                           'final_cands': 0, 'fragment': [int(seed)],
                           'iso': None})
        return []
    fragment = {seed}
    distance = {seed: 0}
    used_edges = set()
    deferred_edges = {tuple(sorted(e)) for e in (prior_deferred_edges or ())}
    heap = []
    _push_edges_from(heap, used_edges, g_R, seed, fragment, graph_floor)

    if record:
        events.append({
            'type': 'seed_start',
            'seed': int(seed),
            'init_cands': len(cands),
            'represented_assignments': represented_assignment_expr(cands),
            'fragment': [int(seed)],
            'p_atoms': sorted({int(v) for c in cands
                               for v in _cand_possible_p_atoms(c)}),
            'cand_patterns': cands_pattern_sample(cands, 5),
        })

    while heap:
        neg_w, u, n = heapq.heappop(heap)
        wbo = -neg_w
        edge = frozenset({u, n})
        if edge in used_edges:
            continue
        used_edges.add(edge)
        if n in fragment:
            if record:
                events.append({
                    'type': 'pop_skip',
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3)},
                    'reason': 'ext_atom already in fragment',
                })
            continue

        n_in_mapping = n in mapping
        if record:
            events.append({
                'type': 'pop',
                'edge': {'frag_atom': int(u), 'ext_atom': int(n),
                         'wbo': round(wbo, 3),
                         'ext_element': g_R.nodes[n]['element']},
                'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                'island_id_at_ext': (int(islands_R[n]) if islands_R and n in islands_R else None),
                'island_image': (int(mapping[n]) if n_in_mapping else None),
                'pre_state': {
                    'fragment_size': len(fragment),
                    'fragment': sorted(int(x) for x in fragment),
                    'cands_count': len(cands),
                    'represented_assignments': represented_assignment_expr(cands),
                    'cands_sample': cands_sample(cands, 5),
                    'cands_pattern_sample': cands_pattern_sample(cands, 5),
                    'p_atoms_in_cands': sorted({int(v) for c in cands
                                                for v in _cand_possible_p_atoms(c)}),
                },
                'heap_top_after_pop': heap_snapshot(heap, used_edges, fragment, mapping, 8),
                'pool_by_frag_atom': pool_by_frag_atom(heap, used_edges, fragment, mapping, g_R),
            })

        # Tentatively add n to the shared fragment. If any weighted extension
        # survives, commit it; otherwise record this edge as deferred boundary.
        old_count = len(cands)
        old_fragment = set(fragment)
        candidate_fragment = fragment | {n}
        dedupe_fragment = set(candidate_fragment)
        if n_in_mapping and islands_R is not None and n in islands_R:
            target_iid = islands_R[n]
            dedupe_fragment |= {
                r for r, k in islands_R.items() if k == target_iid
            }
        dedupe_edges = (
            set(deferred_edges) |
            _frontier_boundary_edges(g_R, dedupe_fragment, graph_floor)
        )
        # Symmetry-compressed incremental extension.  It applies the same
        # element/WBO checks as the concrete incremental matcher, but groups
        # target atoms by local orbit/context before constructing children.
        new_cands = _extend_sym_cands(
            cands, fragment, n, g_R, g_P, mapping,
            iso_tol, islands_R, p_orbits=p_orbits, r_orbits=r_orbits,
            deferred_edges=deferred_edges, anchor_u=u, anchor_wbo=wbo,
            dedupe_edges=dedupe_edges)
        if new_cands:
            cands = new_cands
            ref_dist = distance.get(u, 0) + 1
            if n_in_mapping and islands_R is not None and n in islands_R:
                target_iid = islands_R[n]
                whole_island = [r for r, k in islands_R.items() if k == target_iid]
                candidate_fragment = candidate_fragment | set(whole_island)
                added_extra = [r for r in whole_island if r not in distance]
                for r in added_extra:
                    distance[r] = ref_dist
            fragment = candidate_fragment
            for r in fragment - old_fragment:
                if r not in distance:
                    distance[r] = ref_dist
                _push_edges_from(heap, used_edges, g_R, r, fragment, graph_floor)
            if record:
                added_atoms = sorted(int(r) for r in fragment - old_fragment)
                events.append({
                    'type': 'commit',
                    'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3)},
                    'added': int(n), 'element': g_R.nodes[n]['element'],
                    'atoms_added': added_atoms,
                    'island_size_absorbed': len(added_atoms) if n_in_mapping else None,
                    'island_image': int(mapping[n]) if n_in_mapping else None,
                    'cand_n_value_set': sorted({int(v) for c in cands
                                                for v in cand_possible_values(c, n)}),
                    'cands_before': old_count,
                    'cands_after': len(cands),
                    'represented_assignments_after': represented_assignment_expr(cands),
                    'cands_sample_after': cands_sample(cands, 5),
                    'cands_pattern_after': cands_pattern_sample(cands, 5),
                    'fragment': sorted(int(x) for x in fragment),
                    'p_atoms': sorted({int(v) for c in cands
                                       for v in _cand_possible_p_atoms(c)}),
                    'distance_from_seed': distance[n],
                    'bonds_to_fragment': [(int(u), round(wbo, 3))],
                    'heap_remaining': len(heap),
                    'heap_top': heap_snapshot(heap, used_edges, fragment, mapping, 8),
                    'pool_by_frag_atom': pool_by_frag_atom(heap, used_edges, fragment, mapping, g_R),
                })
        else:
            deferred_edges.add(tuple(sorted((u, n))))
            if record:
                events.append({
                    'type': 'consumed',
                    'deferred': True,
                    'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3),
                             'ext_element': g_R.nodes[n]['element']},
                    'reason': ('merge_failed' if n_in_mapping else 'cut_all_cands'),
                    'island_image': int(mapping[n]) if n_in_mapping else None,
                    'island_id': int(islands_R[n]) if islands_R and n in islands_R else None,
                    'why_per_cand': (why_merge_failed(cands, fragment, n, mapping, islands_R, g_R, g_P, iso_tol) if n_in_mapping
                                     else why_extend_failed(cands, fragment, n, u, wbo, g_R, g_P, mapping, iso_tol)),
                    'cands_count': len(cands),
                    'represented_assignments': represented_assignment_expr(cands),
                    'cands_sample': cands_sample(cands, 5),
                    'cands_pattern_sample': cands_pattern_sample(cands, 5),
                    'fragment': sorted(int(x) for x in fragment),
                    'heap_remaining': len(heap),
                    'heap_top': heap_snapshot(heap, used_edges, fragment, mapping, 8),
                    'pool_by_frag_atom': pool_by_frag_atom(heap, used_edges, fragment, mapping, g_R),
                })

    # heap empty
    if not cands or len(fragment) < min_lock_size:
        if record:
            events.append({
                'type': 'seed_end',
                'result': ('no_cands' if not cands else 'too_small'),
                'final_cands': len(cands),
                'fragment': sorted(int(x) for x in fragment),
                'iso': None,
            })
        return []
    if _set_unique(cands):
        if record:
            iso_map = _cand_map(cands[0])
            iso = {int(k): int(v) for k, v in iso_map.items()}
            events.append({
                'type': 'seed_end', 'result': 'success',
                'final_cands': len(cands),
                'fragment': sorted(int(x) for x in fragment),
                'iso': iso,
                'all_isos': [iso],
                'cand_patterns': cands_pattern_sample(cands, 5),
            })
        return [_IsoResult(_cand_map(cands[0]),
                           deferred_edges=deferred_edges,
                           fragment=fragment,
                           symmetry=_symmetry_state(cands[0]))]
    # Dedup by compressed structural signature.  Open symmetry blocks may
    # still contain many concrete witnesses; only one deterministic witness
    # is returned for each orbit/context-distinct saturation.
    by_set = {}
    for c in cands:
        if isinstance(c, _SymCand):
            key = (
                c.structural_signature(g_R, g_P, r_orbits, p_orbits),
                _boundary_signature(
                    c, g_R, g_P, fragment=fragment,
                    deferred_edges=deferred_edges, r_orbits=r_orbits,
                    p_orbits=p_orbits, locked_mapping=mapping),
            )
        elif p_orbits is not None:
            key = (
                _cand_canon_signature(c, p_orbits),
                _boundary_signature(
                    c, g_R, g_P, fragment=fragment,
                    deferred_edges=deferred_edges, r_orbits=r_orbits,
                    p_orbits=p_orbits, locked_mapping=mapping),
            )
        else:
            key = (
                tuple(sorted(c.items())),
                _boundary_signature(
                    c, g_R, g_P, fragment=fragment,
                    deferred_edges=deferred_edges, r_orbits=r_orbits,
                    p_orbits=p_orbits, locked_mapping=mapping),
            )
        if key not in by_set:
            by_set[key] = _IsoResult(
                _cand_map(c), deferred_edges=deferred_edges,
                fragment=fragment, symmetry=_symmetry_state(c))
    branches = list(by_set.values())[:max_branches]
    if record:
        events.append({
            'type': 'seed_end', 'result': 'branched',
            'final_cands': len(cands), 'n_branches': len(branches),
            'fragment': sorted(int(x) for x in fragment),
            'iso': {int(k): int(v) for k, v in branches[0].items()},
            'all_isos': [{int(k): int(v) for k, v in c.items()}
                         for c in branches],
            'cand_patterns': cands_pattern_sample(cands, 10),
        })
    return branches

"""Live priority-queue fragment growth loop."""
from __future__ import annotations

import heapq
import time

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
from ..matcher.policy import as_node_match_policy
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


def _cand_relation(cand):
    """Return R -> possible P atoms represented by one compressed candidate."""
    relation = {
        r: {p}
        for r, p in _cand_map(cand).items()
    }
    if isinstance(cand, _SymCand):
        for block in cand.blocks:
            for r in block.r_atoms:
                relation.setdefault(r, set()).update(block.p_atoms)
        for items, _mult in cand.alternates:
            for r, p in dict(items).items():
                relation.setdefault(r, set()).add(p)
    return relation


def _candidate_assignment(cand, r_atoms):
    mapping = _cand_map(cand)
    return tuple(
        (int(r), int(mapping[r]))
        for r in sorted(r_atoms)
        if r in mapping
    )


def _island_candidate_symmetry_blocks(cands, fragment):
    """Compressed mapping variation across same-island candidate states.

    Local `_SymCand.blocks` describe ambiguity inside one candidate.  During
    island growth, the live candidate pool can also contain complete correlated
    witnesses of the same local automorphism.  Store that same-island evidence
    as bipartite connected components, not as expanded concrete mappings.
    """
    cands = list(cands or ())
    if len(cands) <= 1:
        return []
    relations = [_cand_relation(c) for c in cands]
    graph = {}
    varied_r = set()
    for r in sorted(fragment):
        possible = set()
        for relation in relations:
            possible.update(relation.get(r, ()))
        if len(possible) <= 1:
            continue
        varied_r.add(r)
        r_node = ('r', int(r))
        graph.setdefault(r_node, set())
        for p in possible:
            p_node = ('p', int(p))
            graph.setdefault(p_node, set()).add(r_node)
            graph[r_node].add(p_node)

    blocks = []
    seen = set()
    block_keys = set()
    for start in list(graph):
        if start in seen:
            continue
        stack = [start]
        comp = set()
        seen.add(start)
        while stack:
            node = stack.pop()
            comp.add(node)
            for nb in graph.get(node, ()):
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        r_atoms = sorted(v for tag, v in comp if tag == 'r')
        p_atoms = sorted(v for tag, v in comp if tag == 'p')
        if not r_atoms or len(p_atoms) <= 1:
            continue
        assignments = {
            _candidate_assignment(c, r_atoms)
            for c in cands
        }
        assignments = {item for item in assignments if item}
        if len(assignments) <= 1 and all(r not in varied_r for r in r_atoms):
            continue
        key = (tuple(r_atoms), tuple(p_atoms))
        if key in block_keys:
            continue
        block_keys.add(key)
        blocks.append({
            'r_atoms': [int(r) for r in r_atoms],
            'p_atoms': [int(p) for p in p_atoms],
            'extendable': False,
            'open': len(r_atoms) < len(p_atoms),
            'assignments': len(assignments) if assignments else None,
            'source': 'island_automorph',
        })
    return blocks


def _remember_symmetry_blocks(accumulated, seen, blocks):
    for block in blocks or ():
        key = (
            tuple(block.get('r_atoms') or ()),
            tuple(block.get('p_atoms') or ()),
            block.get('source'),
        )
        if key in seen:
            continue
        seen.add(key)
        accumulated.append(dict(block))


def _symmetry_state_with_island_automorph(cand, island_blocks, *,
                                          r_orbits=None, p_orbits=None):
    state = _symmetry_state(cand, r_orbits=r_orbits, p_orbits=p_orbits)
    if island_blocks:
        state.setdefault('blocks', []).extend(dict(block)
                                              for block in island_blocks)
    return state


def grow_island(g_R, g_P, seed, mapping,
                graph_floor=0.2,
                iso_tol=0.5,
                min_lock_size=1,
                max_branches=1_000_000,
                events=None,
                islands_R=None,
                p_orbits=None,
                r_orbits=None,
                prior_deferred_edges=None,
                node_policy=None,
                allow_mapped_seed=False,
                profile=None,
                profile_context=None):
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
    node_policy = as_node_match_policy(node_policy)
    record = events is not None
    prof = None
    profile_t0 = None
    if profile is not None:
        profile_t0 = time.perf_counter()
        prof = {
            'seed': int(seed),
            'seed_targets': 0,
            'seed_groups': 0,
            'init_cands': 0,
            'initial_heap': 0,
            'heap_pops': 0,
            'stale_pops': 0,
            'fragment_skip_pops': 0,
            'extend_calls': 0,
            'extend_elapsed_sec': 0.0,
            'max_extend_elapsed_sec': 0.0,
            'max_cands_before': 0,
            'max_cands_after': 0,
            'max_fragment_size': 0,
            'max_heap_len': 0,
            'commits': 0,
            'deferred': 0,
            'merge_calls': 0,
            'free_extend_calls': 0,
            'slowest_extend': None,
        }
        if profile_context:
            prof.update(profile_context)

    def _finish_profile(result, cands_count=0, fragment_set=None, branches=0):
        if prof is None:
            return
        prof['result'] = result
        prof['final_cands'] = int(cands_count)
        prof['final_fragment_size'] = int(len(fragment_set or ()))
        prof['branches'] = int(branches)
        prof['elapsed_sec'] = time.perf_counter() - profile_t0
        profile.append(prof)

    locked_p_atoms = set(mapping.values())
    if seed in mapping:
        if not allow_mapped_seed:
            _finish_profile('already_mapped')
            return []
        seed_targets = [mapping[seed]]
        seed_groups = [(mapping[seed],)]
    else:
        seed_targets = [v for v in g_P.nodes()
                        if v not in locked_p_atoms
                        and node_policy.compatible(g_R, seed, g_P, v)]
        if p_orbits is not None:
            seed_groups = _group_nodes_by_signature(
                seed_targets, lambda v: (node_policy.key(g_P, v),
                                         _orbit_id(p_orbits, v)))
        else:
            seed_groups = [(v,) for v in sorted(seed_targets)]
    if prof is not None:
        prof['seed_targets'] = int(len(seed_targets))
        prof['seed_groups'] = int(len(seed_groups))
    cands = []
    for group in seed_groups:
        if seed in mapping:
            cands.append(_SymCand(
                {seed: mapping[seed]},
                exact_fixed=(seed,)))
        elif len(group) > 1:
            cands.append(_SymCand({seed: group[0]},
                                  (_SymBlock((seed,), group,
                                             extendable=False),)))
        else:
            cands.append(_SymCand({seed: group[0]}))
    if not cands:
        _finish_profile('no_initial_cands')
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
    island_symmetry_blocks = []
    island_symmetry_keys = set()
    heap = []
    _push_edges_from(heap, used_edges, g_R, seed, fragment, graph_floor)
    if prof is not None:
        prof['init_cands'] = int(len(cands))
        prof['initial_heap'] = int(len(heap))
        prof['max_cands_before'] = int(len(cands))
        prof['max_cands_after'] = int(len(cands))
        prof['max_fragment_size'] = int(len(fragment))
        prof['max_heap_len'] = int(len(heap))

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
        if prof is not None:
            prof['heap_pops'] += 1
        wbo = -neg_w
        edge = frozenset({u, n})
        if edge in used_edges:
            if prof is not None:
                prof['stale_pops'] += 1
            continue
        used_edges.add(edge)
        if n in fragment:
            if prof is not None:
                prof['fragment_skip_pops'] += 1
            if record:
                events.append({
                    'type': 'pop_skip',
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3)},
                    'reason': 'ext_atom already in fragment',
                })
            continue

        n_in_mapping = n in mapping
        if prof is not None:
            prof['extend_calls'] += 1
            if n_in_mapping:
                prof['merge_calls'] += 1
            else:
                prof['free_extend_calls'] += 1
        if record:
            events.append({
                'type': 'pop',
                'edge': {'frag_atom': int(u), 'ext_atom': int(n),
                         'wbo': round(wbo, 3),
                         'ext_element': g_R.nodes[n].get('element')},
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
        if prof is not None:
            prof['max_cands_before'] = max(
                int(prof['max_cands_before']), int(old_count))
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
        extend_t0 = time.perf_counter() if prof is not None else None
        new_cands = _extend_sym_cands(
            cands, fragment, n, g_R, g_P, mapping,
            iso_tol, islands_R, p_orbits=p_orbits, r_orbits=r_orbits,
            deferred_edges=deferred_edges, anchor_u=u, anchor_wbo=wbo,
            dedupe_edges=dedupe_edges, node_policy=node_policy)
        if prof is not None:
            extend_elapsed = time.perf_counter() - extend_t0
            prof['extend_elapsed_sec'] += extend_elapsed
            if extend_elapsed > prof['max_extend_elapsed_sec']:
                prof['max_extend_elapsed_sec'] = extend_elapsed
                prof['slowest_extend'] = {
                    'frag_atom': int(u),
                    'ext_atom': int(n),
                    'wbo': float(wbo),
                    'scenario': (
                        'merge_island' if n_in_mapping else 'extend_free'
                    ),
                    'elapsed_sec': extend_elapsed,
                    'cands_before': int(old_count),
                    'cands_after': int(len(new_cands or ())),
                    'fragment_size_before': int(len(fragment)),
                    'candidate_fragment_size': int(len(candidate_fragment)),
                }
        if new_cands:
            if prof is not None:
                prof['commits'] += 1
                prof['max_cands_after'] = max(
                    int(prof['max_cands_after']), int(len(new_cands)))
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
            if prof is not None:
                prof['max_fragment_size'] = max(
                    int(prof['max_fragment_size']), int(len(fragment)))
                prof['max_heap_len'] = max(
                    int(prof['max_heap_len']), int(len(heap)))
            if record:
                added_atoms = sorted(int(r) for r in fragment - old_fragment)
                events.append({
                    'type': 'commit',
                    'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3)},
                    'added': int(n), 'element': g_R.nodes[n].get('element'),
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
            _remember_symmetry_blocks(
                island_symmetry_blocks,
                island_symmetry_keys,
                _island_candidate_symmetry_blocks(cands, fragment),
            )
        else:
            if prof is not None:
                prof['deferred'] += 1
            deferred_edges.add(tuple(sorted((u, n))))
            if record:
                events.append({
                    'type': 'consumed',
                    'deferred': True,
                    'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3),
                             'ext_element': g_R.nodes[n].get('element')},
                    'reason': ('merge_failed' if n_in_mapping else 'cut_all_cands'),
                    'island_image': int(mapping[n]) if n_in_mapping else None,
                    'island_id': int(islands_R[n]) if islands_R and n in islands_R else None,
                    'why_per_cand': (why_merge_failed(cands, fragment, n, mapping, islands_R, g_R, g_P, iso_tol) if n_in_mapping
                                     else why_extend_failed(cands, fragment, n, u, wbo, g_R, g_P, mapping, iso_tol, node_policy=node_policy)),
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
        _finish_profile(
            'no_cands' if not cands else 'too_small',
            len(cands), fragment)
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
        _finish_profile('success', len(cands), fragment, 1)
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
                           symmetry=_symmetry_state_with_island_automorph(
                               cands[0], island_symmetry_blocks,
                               r_orbits=r_orbits, p_orbits=p_orbits))]
    # Dedup by compressed structural signature.  Open symmetry blocks may
    # still contain many concrete witnesses; only one deterministic witness
    # is returned for each orbit/context-distinct saturation.
    _remember_symmetry_blocks(
        island_symmetry_blocks,
        island_symmetry_keys,
        _island_candidate_symmetry_blocks(cands, fragment),
    )
    by_set = {}
    for c in cands:
        if isinstance(c, _SymCand):
            key = (
                c.structural_signature(g_R, g_P, r_orbits, p_orbits),
                _boundary_signature(
                    c, g_R, g_P, fragment=fragment,
                    deferred_edges=deferred_edges, r_orbits=r_orbits,
                    p_orbits=p_orbits, locked_mapping=mapping,
                    node_policy=node_policy),
            )
        elif p_orbits is not None:
            key = (
                _cand_canon_signature(c, p_orbits),
                _boundary_signature(
                    c, g_R, g_P, fragment=fragment,
                    deferred_edges=deferred_edges, r_orbits=r_orbits,
                    p_orbits=p_orbits, locked_mapping=mapping,
                    node_policy=node_policy),
            )
        else:
            key = (
                tuple(sorted(c.items())),
                _boundary_signature(
                    c, g_R, g_P, fragment=fragment,
                    deferred_edges=deferred_edges, r_orbits=r_orbits,
                    p_orbits=p_orbits, locked_mapping=mapping,
                    node_policy=node_policy),
            )
        if key not in by_set:
            by_set[key] = _IsoResult(
                _cand_map(c), deferred_edges=deferred_edges,
                fragment=fragment, symmetry=_symmetry_state_with_island_automorph(
                    c, island_symmetry_blocks, r_orbits=r_orbits,
                    p_orbits=p_orbits))
    branches = list(by_set.values())[:max_branches]
    _finish_profile('branched', len(cands), fragment, len(branches))
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

"""Branch scheduling and mechanism-state dedupe."""
from __future__ import annotations

import random
from collections import Counter, defaultdict

import numpy as np

from ..frag import bond_event_threshold
from ..growth import grow_island
from ..matcher import (
    _boundary_signature,
    _edge_wbo,
    _nauty_orbits,
    _orbit_id,
    _wbo_bucket,
    as_node_match_policy,
)

class BranchLimitExceeded(RuntimeError):
    """Raised when a seed order hits the live branch cap.

    Sweep-cut uses this as a cut-level abort signal: if any seed order for a
    cut reaches the cap, that cut is considered pathological and contributes
    no mechanism witnesses.
    """

    def __init__(self, max_branches, *, seed=None, pass_no=None,
                 branch_count=None, stage=None):
        self.max_branches = int(max_branches)
        self.seed = seed
        self.pass_no = pass_no
        self.branch_count = branch_count
        self.stage = stage
        msg = f"alignment branch cap hit: {branch_count}/{max_branches}"
        if seed is not None:
            msg += f" seed={seed}"
        if pass_no is not None:
            msg += f" pass={pass_no}"
        if stage:
            msg += f" stage={stage}"
        super().__init__(msg)


class _Branch:
    __slots__ = ('mapping', 'islands_R', 'islands_P', 'next_iid',
                 'deferred_edges', 'symmetry_fragments', '_signature_cache')
    def __init__(self):
        self.mapping = {}
        self.islands_R = {}
        self.islands_P = {}
        self.next_iid = 1
        self.deferred_edges = set()
        self.symmetry_fragments = []
        self._signature_cache = None

    @classmethod
    def from_anchor_map(cls, anchor_map):
        b = cls()
        for r, p in sorted((int(r), int(p))
                           for r, p in dict(anchor_map or {}).items()):
            iid = b.next_iid
            b.next_iid += 1
            b.mapping[r] = p
            b.islands_R[r] = iid
            b.islands_P[p] = iid
        return b

    def fork(self):
        b = _Branch()
        b.mapping = dict(self.mapping)
        b.islands_R = dict(self.islands_R)
        b.islands_P = dict(self.islands_P)
        b.next_iid = self.next_iid
        b.deferred_edges = set(self.deferred_edges)
        b.symmetry_fragments = list(self.symmetry_fragments)
        b._signature_cache = self._signature_cache
        return b

    def add_interbranch_symmetry(self, blocks):
        blocks = list(blocks or ())
        if not blocks:
            return
        r_atoms = sorted({
            int(r)
            for block in blocks
            for r in block.get('r_atoms', ())
        })
        self.symmetry_fragments.append({
            'island_idx': 0,
            'fragment': r_atoms,
            'deferred_edges': [],
            'symmetry': {
                'witness': {
                    int(r): int(self.mapping[r])
                    for r in r_atoms
                    if r in self.mapping
                },
                'blocks': blocks,
            },
        })

    def commit(self, iso, g_R, events=None):
        self._signature_cache = None
        touched = set()
        for r in iso:
            if r in self.islands_R:
                touched.add(self.islands_R[r])
        if touched:
            iid = min(touched)
        else:
            iid = self.next_iid
            self.next_iid += 1
        committed_new = []
        relabeled = []
        for r, p in iso.items():
            if r not in self.mapping:
                self.mapping[r] = p
                committed_new.append((int(r), int(p)))
            elif self.islands_R.get(r) != iid:
                relabeled.append((int(r), int(self.islands_R[r])))
            self.islands_R[r] = iid
            self.islands_P[p] = iid
        self.deferred_edges.update(getattr(iso, 'deferred_edges', ()))
        self.symmetry_fragments.append({
            'island_idx': int(iid),
            'fragment': sorted(int(r) for r in getattr(iso, 'fragment', ())),
            'deferred_edges': [list(map(int, e))
                               for e in sorted(getattr(iso, 'deferred_edges', ()))],
            'symmetry': getattr(iso, 'symmetry', {}),
        })
        for r, k in list(self.islands_R.items()):
            if k in touched and k != iid:
                relabeled.append((int(r), int(k)))
                self.islands_R[r] = iid
                self.islands_P[self.mapping[r]] = iid
        if events is not None:
            events.append({
                'type': 'island_locked',
                'island_idx': int(iid),
                'pairs': committed_new,
                'merged_with': sorted(int(t) for t in touched - {iid}),
                'relabeled': relabeled,
                'mapped_total': len(self.mapping),
            })


def _mapping_variation_blocks(mappings, source='interbranch'):
    """Compressed R/P components where equivalent mappings differ."""
    mappings = [
        {int(r): int(p) for r, p in dict(mapping or {}).items()}
        for mapping in mappings
    ]
    mappings = [mapping for mapping in mappings if mapping]
    if len(mappings) <= 1:
        return []

    graph = {}
    for r in sorted({r for mapping in mappings for r in mapping}):
        possible = {
            mapping[r]
            for mapping in mappings
            if r in mapping
        }
        if len(possible) <= 1:
            continue
        r_node = ('r', int(r))
        graph.setdefault(r_node, set())
        for p in possible:
            p_node = ('p', int(p))
            graph.setdefault(p_node, set()).add(r_node)
            graph[r_node].add(p_node)

    blocks = []
    seen = set()
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
            tuple(
                (r, mapping.get(r))
                for r in r_atoms
                if r in mapping
            )
            for mapping in mappings
        }
        assignments = {item for item in assignments if item}
        if len(assignments) <= 1:
            continue
        blocks.append({
            'source': source,
            'r_atoms': r_atoms,
            'p_atoms': p_atoms,
            'extendable': False,
            'open': len(r_atoms) < len(p_atoms),
            'assignments': len(assignments),
        })
    return blocks


def _orbit_pair(a, b, orbits):
    oa, ob = _orbit_id(orbits, a), _orbit_id(orbits, b)
    return (oa, ob) if oa <= ob else (ob, oa)


def _chemistry_orbit_signature(mapping, g_R, g_P, r_orbits=None, p_orbits=None,
                               dwbo_threshold=0.5,
                               metal_dwbo_threshold=None):
    """Broken/formed bond signature in joint R-orbit/P-orbit space."""
    br_pairs = []
    fm_pairs = []
    mapped = sorted(mapping)
    elements_R = [g_R.nodes[x].get('element') for x in range(len(g_R.nodes))]
    for i, u in enumerate(mapped):
        for v in mapped[i + 1:]:
            pu, pv = mapping[u], mapping[v]
            threshold = bond_event_threshold(
                elements_R, u, v,
                default_threshold=dwbo_threshold,
                metal_threshold=metal_dwbo_threshold)
            wR = _edge_wbo(g_R, u, v)
            wP = _edge_wbo(g_P, pu, pv)
            if wR - wP >= threshold:
                br_pairs.append((
                    _orbit_pair(u, v, r_orbits),
                    _orbit_pair(pu, pv, p_orbits),
                ))
            elif wP - wR >= threshold:
                fm_pairs.append((
                    _orbit_pair(u, v, r_orbits),
                    _orbit_pair(pu, pv, p_orbits),
                ))
    return tuple(sorted(br_pairs)), tuple(sorted(fm_pairs))


def _alignment_state_signature(mapping, deferred_edges, g_R, g_P,
                               r_orbits=None, p_orbits=None, core_R=(),
                               node_policy=None):
    core_R = tuple(sorted(core_R or ()))
    core_set = set(core_R)
    fixed = []
    for r, p in sorted(mapping.items()):
        if r in core_set:
            fixed.append((('exact', r), ('exact', p)))
        else:
            fixed.append((
                ('orbit', _orbit_id(r_orbits, r)),
                ('orbit', _orbit_id(p_orbits, p)),
            ))
    internal_pairs = []
    mapped = sorted(mapping)
    for i, r1 in enumerate(mapped):
        for r2 in mapped[i + 1:]:
            p1, p2 = mapping[r1], mapping[r2]
            internal_pairs.append((
                _orbit_pair(r1, r2, r_orbits),
                _orbit_pair(p1, p2, p_orbits),
                _wbo_bucket(_edge_wbo(g_R, r1, r2)),
                _wbo_bucket(_edge_wbo(g_P, p1, p2)),
            ))
    boundary = _boundary_signature(
        mapping, g_R, g_P, fragment=set(mapping),
        deferred_edges=deferred_edges, r_orbits=r_orbits,
        p_orbits=p_orbits, locked_mapping=mapping,
        node_policy=node_policy)
    return tuple(fixed), tuple(sorted(internal_pairs)), boundary


def _normalize_anchor_map(anchor_map, g_R, g_P):
    if not anchor_map:
        return {}
    out = {}
    used_p = set()
    for raw_r, raw_p in dict(anchor_map).items():
        r, p = int(raw_r), int(raw_p)
        if r not in g_R:
            raise ValueError(f"anchor R atom is not in query graph: {r}")
        if p not in g_P:
            raise ValueError(f"anchor P atom is not in target graph: {p}")
        if r in out and out[r] != p:
            raise ValueError(f"conflicting anchors for R atom {r}")
        if p in used_p:
            raise ValueError(f"anchor target atom is used more than once: {p}")
        out[r] = p
        used_p.add(p)
    return out


def _anchor_nodes_compatible(anchor_map, g_R, g_P, node_policy):
    return all(
        node_policy.compatible(g_R, r, g_P, p)
        for r, p in anchor_map.items()
    )


def find_islands(g_R, g_P, seed_order,
                 graph_floor=0.2, iso_tol=1.0,
                 dwbo_threshold=0.5, symmetry_wbo_tol=0.2,
                 metal_dwbo_threshold=None,
                 max_branches=1_000_000, events=None,
                 orbit_dedup=True, core_R=None,
                 stop_when_core_mapped=False,
                 p_orbits=None, r_orbits=None,
                 abort_on_branch_cap=False,
                 node_policy=None,
                 anchor_map=None,
                 profile=None):
    """Run growth over a single seed ordering, branching on
    non-set-unique locks. Returns list of _Branch.

    Optional `events` only records the FIRST (best-mapped) branch's
    trajectory — multi-branch traces would be confusing on a slider.

    orbit_dedup: when True (default), uses exact automorphism orbits on a
    tolerance-bucketed WBO-colored graph via pynauty as the compression key
    for orbit-equivalence. Callers that run many seed orders against the same
    graph can pass precomputed `p_orbits` / `r_orbits`; missing maps are
    computed here. Chemistry signatures and active R-pair extension checks
    remain the verifier.

    core_R: optional R atoms that define the scoring-relevant alignment.
    When supplied, branch dedup switches to exact core mapping as soon as
    every core atom is mapped.  With stop_when_core_mapped=True the search
    returns once all live branches have mapped the core; spectators remain
    represented only by the compressed witness state already discovered.

    anchor_map: optional exact R->P constraints.  Anchors are preloaded as
    locked single-atom islands, but their R atoms may still be used as growth
    seeds so the surrounding subgraph can be discovered from the constraint.
    """
    node_policy = as_node_match_policy(node_policy)
    anchor_map = _normalize_anchor_map(anchor_map, g_R, g_P)
    if not _anchor_nodes_compatible(anchor_map, g_R, g_P, node_policy):
        return []
    if orbit_dedup:
        if p_orbits is None:
            p_orbits = _nauty_orbits(
                g_P, wbo_tol=symmetry_wbo_tol,
                node_policy=node_policy)
        if r_orbits is None:
            r_orbits = _nauty_orbits(
                g_R, wbo_tol=symmetry_wbo_tol,
                node_policy=node_policy)
    else:
        p_orbits = None
        r_orbits = None
    core_R = tuple(sorted(set(core_R or ())))
    seed_order = list(dict.fromkeys(seed_order))
    branches = [_Branch.from_anchor_map(anchor_map)]
    progressed = True
    pass_no = 0

    def _core_complete(mapping):
        return core_R and all(r in mapping for r in core_R)

    def _mapped_anchor_can_seed(branch, seed):
        if seed not in anchor_map:
            return False
        seed_iid = branch.islands_R.get(seed)
        for nb in g_R.neighbors(seed):
            if g_R[seed][nb].get('wbo', 0.0) < graph_floor:
                continue
            if branch.islands_R.get(nb) != seed_iid:
                return True
        return False

    signature_cache = {}

    def _deferred_key(deferred_edges):
        return tuple(sorted(tuple(sorted(e)) for e in deferred_edges))

    def _mapping_signature(mapping, deferred_edges=()):
        # Full chemistry signatures are expensive for large near-complete
        # mappings.  The same branch state is checked repeatedly while seeds
        # are carried forward and cross-branch dedupe runs, so cache by exact
        # mapping/deferred state inside this find_islands invocation.
        cache_key = (
            tuple(sorted(mapping.items())),
            _deferred_key(deferred_edges),
        )
        cached = signature_cache.get(cache_key)
        if cached is not None:
            return cached
        core_key = tuple((r, mapping[r]) for r in core_R if r in mapping)
        if core_R and len(core_key) == len(core_R):
            sig = ('core_complete', core_key)
            signature_cache[cache_key] = sig
            return sig
        deferred_boundary = _boundary_signature(
            mapping, g_R, g_P, fragment=set(mapping),
            deferred_edges=deferred_edges, r_orbits=r_orbits,
            p_orbits=p_orbits, locked_mapping=mapping,
            node_policy=node_policy)
        sig = (
            'mechanism_state',
            core_key,
            _chemistry_orbit_signature(
                mapping, g_R, g_P, r_orbits, p_orbits,
                dwbo_threshold=dwbo_threshold,
                metal_dwbo_threshold=metal_dwbo_threshold),
            deferred_boundary,
        )
        signature_cache[cache_key] = sig
        return sig

    def _branch_signature(branch):
        if branch._signature_cache is None:
            branch._signature_cache = _mapping_signature(
                branch.mapping, branch.deferred_edges)
        return branch._signature_cache

    def _state_key(branch):
        return (
            tuple(sorted(branch.mapping.items())),
            tuple(sorted(branch.islands_R.items())),
            tuple(sorted(branch.islands_P.items())),
            _deferred_key(branch.deferred_edges),
        )

    def _hit_branch_cap(count, seed, stage):
        if not abort_on_branch_cap:
            return
        if events is not None:
            events.append({
                'type': 'done',
                'mapped': len(branches[0].mapping) if branches else 0,
                'stop_reason': 'branch_cap',
                'max_branches': int(max_branches),
                'branches': int(count),
                'seed': int(seed) if seed is not None else None,
                'stage': stage,
            })
        raise BranchLimitExceeded(
            max_branches,
            seed=int(seed) if seed is not None else None,
            pass_no=pass_no,
            branch_count=int(count),
            stage=stage,
        )

    while progressed:
        progressed = False
        pass_no += 1
        if events is not None:
            events.append({'type': 'pass_start', 'pass': pass_no,
                           'mapped': len(branches[0].mapping)})
        for seed in seed_order:
            new_branches = []
            pending_seen = set()

            def _append_pending(branch):
                sig = _branch_signature(branch)
                if sig in pending_seen:
                    return False
                pending_seen.add(sig)
                new_branches.append(branch)
                return True

            for bi, b in enumerate(branches):
                if len(new_branches) >= max_branches:
                    _hit_branch_cap(len(new_branches), seed, 'pre_branch_loop')
                    break
                if stop_when_core_mapped and _core_complete(b.mapping):
                    _append_pending(b)
                    continue
                seed_is_mapped_anchor = seed in b.mapping and seed in anchor_map
                if seed in b.mapping and not seed_is_mapped_anchor:
                    _append_pending(b)
                    continue
                if seed_is_mapped_anchor and not _mapped_anchor_can_seed(b, seed):
                    _append_pending(b)
                    continue
                # Only record events for branch 0 to keep trace linear
                ev_arg = events if (events is not None and bi == 0) else None
                isos = grow_island(g_R, g_P, seed, b.mapping,
                                    graph_floor=graph_floor, iso_tol=iso_tol,
                                    max_branches=max_branches,
                                    events=ev_arg,
                                    islands_R=b.islands_R,
                                    p_orbits=p_orbits,
                                    r_orbits=r_orbits,
                                    prior_deferred_edges=b.deferred_edges,
                                    node_policy=node_policy,
                                    allow_mapped_seed=seed_is_mapped_anchor,
                                    profile=profile,
                                    profile_context={
                                        'pass': int(pass_no),
                                        'branch_index': int(bi),
                                        'mapped_before': len(b.mapping),
                                    })
                if not isos:
                    _append_pending(b)
                    continue
                # Dedup island results by weighted alignment state before
                # forking branches.  This is still pre-mechanism dedupe: the
                # one-hop/deferred boundary is part of the key so future
                # distinguishability is not erased.
                seen_state = {}
                seen_state_mappings = defaultdict(list)
                for iso in isos:
                    full_m = dict(b.mapping); full_m.update(iso)
                    full_deferred = set(b.deferred_edges)
                    full_deferred.update(getattr(iso, 'deferred_edges', ()))
                    state_key = _mapping_signature(full_m, full_deferred)
                    if state_key not in seen_state:
                        seen_state[state_key] = iso
                    seen_state_mappings[state_key].append(full_m)
                deduped_isos = list(seen_state.values())
                for state_key, iso in seen_state.items():
                    blocks = _mapping_variation_blocks(
                        seen_state_mappings[state_key],
                        source='interbranch')
                    if blocks:
                        symmetry = dict(getattr(iso, 'symmetry', {}) or {})
                        symmetry['blocks'] = (
                            list(symmetry.get('blocks') or []) + blocks
                        )
                        iso.symmetry = symmetry
                for ii, iso in enumerate(deduped_isos):
                    if len(new_branches) >= max_branches:
                        _hit_branch_cap(len(new_branches), seed, 'island_fork')
                        break
                    before_state = _state_key(b)
                    b2 = b.fork()
                    b2.commit(iso, g_R,
                              events=events if (bi == 0 and ii == 0) else None)
                    after_state = _state_key(b2)
                    if after_state == before_state:
                        if _append_pending(b):
                            continue
                    elif _append_pending(b2):
                        progressed = True
            new_branches.sort(key=lambda b: -len(b.mapping))
            # Cross-branch dedup uses the same mechanism-state key.  It
            # collapses orbit-equivalent spectator permutations but keeps
            # states separated when their deferred one-hop boundary differs.
            seen = {}
            uniq = []
            for b in new_branches:
                state_sig = _branch_signature(b)
                kept = seen.get(state_sig)
                if kept is not None:
                    kept.add_interbranch_symmetry(
                        _mapping_variation_blocks(
                            [kept.mapping, b.mapping],
                            source='interbranch'))
                    continue
                seen[state_sig] = b
                uniq.append(b)
                if len(uniq) >= max_branches:
                    _hit_branch_cap(len(uniq), seed, 'cross_branch_dedupe')
                    break
            branches = uniq
            # Soft warning: pathological symmetry can blow this up. Default
            # cap is 1e6 so we should never hit it on real molecules; if we
            # do, surface it so we know.
            if len(branches) >= 10_000:
                import sys
                print(f"  [warn] alignment branch count = {len(branches)}  "
                      f"max_branches={max_branches}  "
                      f"new_branches_in={len(new_branches)}",
                      file=sys.stderr, flush=True)
            if (stop_when_core_mapped and core_R and branches and
                    all(_core_complete(b.mapping) for b in branches)):
                if events is not None:
                    events.append({'type': 'done',
                                   'mapped': len(branches[0].mapping),
                                   'stop_reason': 'core_mapped',
                                   'core_size': len(core_R)})
                return branches
    if events is not None:
        events.append({'type': 'done',
                       'mapped': len(branches[0].mapping)})
    return branches


def _ordered_seed_nodes(g_R, rng, common_element_threshold=3):
    nodes = list(g_R.nodes())
    element_counts = Counter(g_R.nodes[n].get('element') for n in nodes)
    threshold = max(1, int(common_element_threshold))

    def ambiguous_isolated(n):
        element = g_R.nodes[n].get('element')
        return int(g_R.degree[n]) == 0 and element_counts[element] > threshold

    contextual = [n for n in nodes if not ambiguous_isolated(n)]
    ambiguous = [n for n in nodes if ambiguous_isolated(n)]
    rng.shuffle(contextual)
    rng.shuffle(ambiguous)
    return contextual + ambiguous


def _generate_seed_orders(g_R, n_trials, rng_seed=42,
                          common_element_threshold=1):
    """Generate at most `n_trials` deterministic seed orderings.

    Seeds are deterministic random orderings, not chemistry/core-prioritized
    rankings.  The only special case is an isolated atom whose element occurs
    multiple times: that atom has no graph context, so it is retained but kept
    behind contextual atoms and is not used as an anchor unless no contextual
    anchor exists.
    """
    nodes = _ordered_seed_nodes(
        g_R,
        random.Random(rng_seed),
        common_element_threshold=common_element_threshold,
    )
    element_counts = Counter(g_R.nodes[n].get('element') for n in g_R.nodes())
    threshold = max(1, int(common_element_threshold))

    def ambiguous_isolated(n):
        element = g_R.nodes[n].get('element')
        return int(g_R.degree[n]) == 0 and element_counts[element] > threshold

    anchors = [n for n in nodes if not ambiguous_isolated(n)] or nodes
    orders = []
    n_trials = max(0, int(n_trials))
    if not anchors:
        return orders
    for idx in range(n_trials):
        anchor = anchors[idx % len(anchors)]
        rest = [
            n for n in _ordered_seed_nodes(
                g_R,
                random.Random(rng_seed + idx + 1),
                common_element_threshold=common_element_threshold,
            )
            if n != anchor
        ]
        orders.append([anchor] + rest)
    return orders

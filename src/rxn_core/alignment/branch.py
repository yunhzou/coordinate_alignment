"""Branch scheduling, mechanism-state dedupe, and symmetry repair."""
from __future__ import annotations

import random
from collections import Counter, defaultdict, deque

import numpy as np

from ..frag import classify_bonds, is_metal_element
from ..growth import IslandBranchLimitExceeded, grow_island
from ..matcher import (
    _boundary_signature,
    _edge_wbo,
    _nauty_atom_generators,
    _nauty_orbits,
    _orbit_id,
    _wbo_bucket,
    as_node_match_policy,
)

SYM_REPAIR_MAX_EVALS = 20000
# States scored per vectorised batch inside symmetry_repair_mapping; bounds the
# transient (states x local pairs) work arrays.
_SYM_REPAIR_SCORE_CHUNK = 512


class BranchLimitExceeded(RuntimeError):
    """Raised when an explicitly enumerated result set exceeds its cap."""

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
                 'deferred_edges', 'symmetry_paths', '_signature_cache')
    def __init__(self):
        self.mapping = {}
        self.islands_R = {}
        self.islands_P = {}
        self.next_iid = 1
        self.deferred_edges = set()
        self.symmetry_paths = [[]]
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
        b.symmetry_paths = [list(path) for path in self.symmetry_paths]
        b._signature_cache = self._signature_cache
        return b

    @property
    def symmetry_fragments(self):
        """Compatibility view of the first retained analytical path."""
        return self.symmetry_paths[0]

    def merge_exact_paths(self, other):
        """Union histories after exact equality of cumulative search state."""
        seen = {_freeze_branch_value(path) for path in self.symmetry_paths}
        for path in other.symmetry_paths:
            key = _freeze_branch_value(path)
            if key not in seen:
                seen.add(key)
                self.symmetry_paths.append(list(path))

    def add_interbranch_symmetry(self, blocks):
        blocks = list(blocks or ())
        if not blocks:
            return
        r_atoms = sorted({
            int(r)
            for block in blocks
            for r in block.get('r_atoms', ())
        })
        record = {
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
        }
        for path in self.symmetry_paths:
            path.append(record)

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
        record = {
            'island_idx': int(iid),
            'fragment': sorted(int(r) for r in getattr(iso, 'fragment', ())),
            'deferred_edges': [list(map(int, e))
                               for e in sorted(getattr(iso, 'deferred_edges', ()))],
            'symmetry': getattr(iso, 'symmetry', {}),
        }
        for path in self.symmetry_paths:
            path.append(record)
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


def _freeze_branch_value(value):
    if isinstance(value, dict):
        return tuple(sorted((str(key), _freeze_branch_value(item))
                            for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_branch_value(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return tuple(sorted((_freeze_branch_value(item) for item in value),
                            key=repr))
    return value


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


def _change_score_from_events(broken, formed):
    """Score of an already classified event set (former loop, unchanged)."""
    delta = 0.0
    for a, b, wR, wP in broken:
        delta += abs(float(wR or 0.0) - float(wP or 0.0))
    for a, b, wR, wP in formed:
        delta += abs(float(wP or 0.0) - float(wR or 0.0))
    return (len(broken) + len(formed), round(delta, 12))


def _mapping_change_score(mapping, wbo_R, wbo_P, dwbo_threshold=0.5,
                          elements_R=None, elements_P=None,
                          metal_dwbo_threshold=None):
    broken, formed, _, _ = classify_bonds(
        mapping, wbo_R, wbo_P, dwbo_threshold=dwbo_threshold,
        elements_R=elements_R, elements_P=elements_P,
        metal_dwbo_threshold=metal_dwbo_threshold)
    return _change_score_from_events(broken, formed)


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


def symmetry_repair_mapping(mapping, wbo_R, wbo_P, g_R, g_P, p_orbits,
                            dwbo_threshold=0.5, bond_floor=0.2,
                            metal_dwbo_threshold=None,
                            min_changes=1, full_permutation_size=6,
                            max_evals=SYM_REPAIR_MAX_EVALS,
                            return_stats=False):
    """Choose the best concrete realization inside product symmetry orbits.

    The matcher intentionally compresses high-symmetry choices and keeps one
    witness.  A witness is not chemistry: if changed-bond endpoints sit inside
    an unresolved P orbit, another realization of the same compressed match may
    remove fake broken/formed bonds.  Every tested reassignment is generated by
    pynauty while all atoms outside that local target set are fixed; sharing a
    coarse orbit alone is never treated as permission to swap.
    """
    if not mapping or p_orbits is None:
        return (mapping, {'enabled': False}) if return_stats else mapping
    mapping0 = dict(mapping)
    elements_R = [g_R.nodes[i].get('element') for i in range(wbo_R.shape[0])]
    elements_P = [g_P.nodes[i].get('element') for i in range(wbo_P.shape[0])]
    base_broken, base_formed, _, _ = classify_bonds(
        mapping0, wbo_R, wbo_P, dwbo_threshold=dwbo_threshold,
        elements_R=elements_R, elements_P=elements_P,
        metal_dwbo_threshold=metal_dwbo_threshold)
    base_changes = len(base_broken) + len(base_formed)
    stats = {
        'enabled': True,
        'base_changes': int(base_changes),
        'repaired': False,
        'groups': [],
        'evaluated': 0,
        'capped': False,
    }
    # This is a computational guard only.  Any nonzero changed-bond witness may
    # be an arbitrary representative artifact, so the default is to attempt
    # repair for every nonzero case.
    if base_changes < max(1, int(min_changes)):
        return (mapping0, stats) if return_stats else mapping0

    inv = {v: r for r, v in mapping0.items()}
    affected = set()
    for a, b, _, _ in base_broken:
        affected.add(a); affected.add(b)
    for a, b, _, _ in base_formed:
        if a in inv:
            affected.add(inv[a])
        if b in inv:
            affected.add(inv[b])
    if not affected:
        return (mapping0, stats) if return_stats else mapping0

    local = set(affected)
    for r in list(affected):
        if r not in g_R:
            continue
        for nb in g_R.neighbors(r):
            if g_R[r][nb].get('wbo', 0.0) >= bond_floor:
                local.add(nb)
    for p in {x for edge in base_formed for x in edge[:2]}:
        if p not in g_P:
            continue
        for q in g_P.neighbors(p):
            if q in inv and g_P[p][q].get('wbo', 0.0) >= bond_floor:
                local.add(inv[q])

    touched_orbits = set()
    for r in affected:
        if r in mapping0:
            p = mapping0[r]
            touched_orbits.add((g_R.nodes[r].get('element'), p_orbits[p]))
    groups = defaultdict(list)
    for r, p in mapping0.items():
        key = (g_R.nodes[r].get('element'), p_orbits[p])
        if key in touched_orbits:
            groups[key].append(r)
    groups = {
        key: sorted(rs) for key, rs in groups.items()
        if len(rs) > 1
    }
    if not groups:
        return (mapping0, stats) if return_stats else mapping0

    nR = wbo_R.shape[0]
    # Vectorised form of the former O(N^2) local_pairs loop.  Exact because
    # np.triu_indices enumerates the same i<j row-major pair order, the mask
    # applies the same "both mapped and at least one local" test, the WBO
    # values are the same float64 entries, and the threshold column applies
    # bond_event_threshold's rule (metal cutoff iff either endpoint is a
    # metal, else default) with the same float() conversions.
    mapped_mask = np.zeros(nR, dtype=bool)
    for r in mapping0:
        r = int(r)
        if 0 <= r < nR:
            mapped_mask[r] = True
    local_mask = np.zeros(nR, dtype=bool)
    for r in local:
        r = int(r)
        if 0 <= r < nR:
            local_mask[r] = True
    upper_i, upper_j = np.triu_indices(nR, 1)
    keep = (
        mapped_mask[upper_i] & mapped_mask[upper_j]
        & (local_mask[upper_i] | local_mask[upper_j]))
    pair_i = upper_i[keep].astype(np.intp)
    pair_j = upper_j[keep].astype(np.intp)
    pair_wbo_R = np.asarray(wbo_R)[pair_i, pair_j].astype(float)
    if elements_R is None or metal_dwbo_threshold is None:
        pair_threshold = np.full(pair_i.shape, float(dwbo_threshold))
    else:
        metal_mask = np.array(
            [is_metal_element(element) for element in elements_R],
            dtype=bool)
        pair_threshold = np.where(
            metal_mask[pair_i] | metal_mask[pair_j],
            float(metal_dwbo_threshold), float(dwbo_threshold))
    pair_r_active = pair_wbo_R >= float(bond_floor)
    current_images = np.full(nR, -1, dtype=np.intp)
    for r, p in mapping0.items():
        current_images[int(r)] = int(p)
    wbo_P_array = np.asarray(wbo_P)
    bond_floor_value = float(bond_floor)

    def local_scores(image_rows):
        """Exact former single-state score applied to every row at once.

        Each row undergoes the same elementwise IEEE operations the former
        per-state version used (subtract, abs, compare, select, scale by
        0.01), so a row's contribution vector is bit-identical to scoring that
        state alone; only the per-call interpreter overhead is shared.
        ``|d| >= t`` replaces ``(d >= t) | (-d >= t)``: the two agree for every
        double, including NaN, infinities and signed zero.  The former
        ``sum(contribution.tolist(), 0.0)`` was a strict left-to-right chain of
        double additions; ``np.add.accumulate`` along each row performs the
        identical chain (``out[i] = out[i-1] + x[i]``, no pairwise
        reassociation), and its last column is therefore the same double.  The
        leading ``0.0 +`` of ``sum`` is inert because contributions are
        ``abs`` values, ``abs * 0.01`` or the literal ``0.0`` (never ``-0.0``).
        The 12-place rounding tie-break is applied to that same double.
        """
        pair_wbo_P = wbo_P_array[
            image_rows[:, pair_i], image_rows[:, pair_j]]
        difference = pair_wbo_R - pair_wbo_P
        magnitude = np.abs(difference)
        changed_mask = magnitude >= pair_threshold
        contribution = np.where(
            changed_mask,
            magnitude,
            np.where(
                pair_r_active | (pair_wbo_P >= bond_floor_value),
                magnitude * 0.01,
                0.0,
            ),
        )
        if contribution.shape[1]:
            totals = np.add.accumulate(contribution, axis=1)[:, -1]
        else:
            totals = np.zeros(contribution.shape[0])
        return [
            (int(count), round(total, 12))
            for count, total in zip(
                np.count_nonzero(changed_mask, axis=1).tolist(),
                totals.tolist())
        ]

    current = dict(mapping0)
    current_score = local_scores(current_images[None, :])[0]
    repaired_any = False
    for key, rs in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
        stats['groups'].append({
            'element': key[0],
            'orbit': int(key[1]),
            'size': len(rs),
            'atoms': [int(r) for r in rs],
        })
    ordered_groups = [rs for _, rs in sorted(
        groups.items(), key=lambda item: (-len(item[1]), item[0]))]
    for group_index, rs in enumerate(ordered_groups):
        if stats['evaluated'] >= max_evals:
            stats['capped'] = True
            break
        targets = tuple(current[r] for r in rs)
        target_set = set(targets)
        atom_tags = {
            int(p): ('repair_fixed', int(p))
            for p in g_P.nodes()
            if int(p) not in target_set
        }
        for p in target_set:
            atom_tags[int(p)] = ('repair_group', int(group_index))
        generators = _nauty_atom_generators(
            g_P, wbo_tol=float(getattr(p_orbits, 'wbo_tol', 0.2)),
            atom_color_tags=atom_tags)

        states = {targets}
        queue = deque([targets])
        while queue and stats['evaluated'] < max_evals:
            state = queue.popleft()
            for generator in generators:
                # map(get, state, state) calls get(p, p) per atom exactly as
                # the former generator expression did, without a frame per
                # element.
                image = tuple(map(generator.get, state, state))
                if image in states:
                    continue
                states.add(image)
                queue.append(image)
                if len(states) + stats['evaluated'] >= max_evals:
                    stats['capped'] = bool(queue)
                    break

        best_score = current_score
        best_state = None
        rs_array = np.asarray(rs, dtype=np.intp)
        ordered_states = sorted(states)
        # The former loop checked ``evaluated >= max_evals`` before each state
        # in this sorted order, so it scored exactly the first ``budget``
        # states and flagged the cap iff a state was left over.  Scoring them
        # in batches reproduces the same scores in the same order; the
        # ``score < best_score`` scan below is unchanged.
        budget = max_evals - stats['evaluated']
        if len(ordered_states) > budget:
            stats['capped'] = True
        to_score = ordered_states[:max(budget, 0)]
        for start in range(0, len(to_score), _SYM_REPAIR_SCORE_CHUNK):
            chunk = to_score[start:start + _SYM_REPAIR_SCORE_CHUNK]
            image_rows = np.repeat(
                current_images[None, :], len(chunk), axis=0)
            image_rows[:, rs_array] = np.asarray(chunk, dtype=np.intp)
            for state, score in zip(chunk, local_scores(image_rows)):
                if score < best_score:
                    best_score = score
                    best_state = state
        stats['evaluated'] += len(to_score)
        if best_state is not None:
            current.update(zip(rs, best_state))
            current_images[rs_array] = best_state
            current_score = best_score
            repaired_any = True

    best = current
    # base_broken/base_formed were classified from mapping0 with exactly the
    # arguments _mapping_change_score(mapping0, ...) would pass, so scoring
    # them reproduces that call; when no group replaced an image, ``best``
    # holds the identical items as mapping0 and therefore the same score.
    base_score = _change_score_from_events(base_broken, base_formed)
    if repaired_any:
        best_score = _mapping_change_score(
            best, wbo_R, wbo_P,
            dwbo_threshold=dwbo_threshold,
            elements_R=elements_R,
            elements_P=elements_P,
            metal_dwbo_threshold=metal_dwbo_threshold)
    else:
        best_score = base_score
    stats['best_changes'] = int(best_score[0])
    if best_score < base_score:
        stats['repaired'] = True
        return (best, stats) if return_stats else best
    return (mapping0, stats) if return_stats else mapping0


def find_islands(g_R, g_P, seed_order,
                 graph_floor=0.2, iso_tol=1.0,
                 dwbo_threshold=0.5, symmetry_wbo_tol=0.2,
                 metal_dwbo_threshold=None,
                 max_branches=1_000_000, events=None,
                 orbit_dedup=True, core_R=None,
                 stop_when_core_mapped=False,
                 p_orbits=None, r_orbits=None,
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
    # Orbit equivalence and edge compatibility must describe the same
    # tolerance model.  Keep symmetry_wbo_tol in the signature for callers
    # written against older releases, but do not permit it to diverge.
    symmetry_wbo_tol = float(iso_tol)
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
    # Completed hierarchy/coset equivalence is richer than the live
    # mapping/island state.  Do not quotient distinct concrete live states by
    # endpoint automorphism here: two such states can carry different exact
    # fragment-assignment domains and therefore lead to different mechanisms.
    # Exact semantic family dedupe is performed after AAM has the complete
    # hierarchy available.
    branch_canonicalizer = None
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

    def _deferred_key(deferred_edges):
        return tuple(sorted(tuple(sorted(e)) for e in deferred_edges))

    def _progress_key(branch):
        return (
            tuple(sorted(branch.mapping.items())),
            tuple(sorted(branch.islands_R.items())),
            tuple(sorted(branch.islands_P.items())),
            _deferred_key(branch.deferred_edges),
        )

    def _island_partition(branch):
        groups = {}
        for atom, island in branch.islands_R.items():
            groups.setdefault(int(island), []).append(int(atom))
        return tuple(sorted(tuple(sorted(atoms))
                            for atoms in groups.values()))

    def _branch_signature(branch):
        if branch._signature_cache is not None:
            return branch._signature_cache
        signature = _progress_key(branch)
        branch._signature_cache = signature
        return signature

    def _merge_equivalent_paths(kept, other):
        if kept.mapping == other.mapping:
            kept.merge_exact_paths(other)
            return
        raise RuntimeError("nonidentical branches shared an exact state key")

    while progressed:
        progressed = False
        pass_no += 1
        if events is not None:
            events.append({'type': 'pass_start', 'pass': pass_no,
                           'mapped': len(branches[0].mapping)})
        for seed in seed_order:
            new_branches = []
            pending_seen = {}

            def _admit_subtree(subtree, *, made_progress=False,
                               source_branch=None):
                """Admit one parent branch's result atomically.

                A subtree that would take the live post-dedupe leaf count over
                ``max_branches`` is discarded by itself.  Sibling branches and
                other seed orders remain valid; exactly ``max_branches`` live
                leaves is allowed.
                """
                nonlocal progressed
                additions = []
                local_seen = {}
                for candidate in subtree:
                    sig = _branch_signature(candidate)
                    if sig in pending_seen:
                        _merge_equivalent_paths(
                            pending_seen[sig], candidate)
                        continue
                    if sig in local_seen:
                        _merge_equivalent_paths(local_seen[sig], candidate)
                        continue
                    local_seen[sig] = candidate
                    additions.append((sig, candidate))
                if len(new_branches) + len(additions) > max_branches:
                    if profile is not None:
                        profile.append({
                            'seed': int(seed),
                            'pass': int(pass_no),
                            'branch_index': (
                                None if source_branch is None
                                else int(source_branch)),
                            'mapped_before': (
                                0 if source_branch is None
                                else len(branches[source_branch].mapping)),
                            'result': 'subtree_branch_cap',
                            'stage': 'combined_live_leaves',
                            'subtree_branches': len(additions),
                            'live_branches_before': len(new_branches),
                            'max_branches': int(max_branches),
                            'elapsed_sec': 0.0,
                            'extend_elapsed_sec': 0.0,
                        })
                    return False
                for sig, candidate in additions:
                    pending_seen[sig] = candidate
                    new_branches.append(candidate)
                if additions and made_progress:
                    progressed = True
                return bool(additions)

            for bi, b in enumerate(branches):
                if stop_when_core_mapped and _core_complete(b.mapping):
                    _admit_subtree([b], source_branch=bi)
                    continue
                seed_is_mapped_anchor = seed in b.mapping and seed in anchor_map
                if seed in b.mapping and not seed_is_mapped_anchor:
                    _admit_subtree([b], source_branch=bi)
                    continue
                if seed_is_mapped_anchor and not _mapped_anchor_can_seed(b, seed):
                    _admit_subtree([b], source_branch=bi)
                    continue
                # Only record events for branch 0 to keep trace linear
                ev_arg = events if (events is not None and bi == 0) else None
                try:
                    isos = grow_island(
                        g_R, g_P, seed, b.mapping,
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
                except IslandBranchLimitExceeded:
                    # Only this parent branch's descendant subtree is
                    # pathological.  Remove it and continue its siblings.
                    continue
                if not isos:
                    _admit_subtree([b], source_branch=bi)
                    continue
                # Fragment candidates were already quotiented by an exact
                # pynauty transporter inside ``grow_island``.  Distinct
                # surviving candidates are distinct symbolic branch choices;
                # do not collapse them by a mechanism/event summary here.
                deduped_isos = list(isos)
                subtree = []
                subtree_progressed = False
                for ii, iso in enumerate(deduped_isos):
                    before_state = _progress_key(b)
                    b2 = b.fork()
                    b2.commit(iso, g_R,
                              events=events if (bi == 0 and ii == 0) else None)
                    after_state = _progress_key(b2)
                    if after_state == before_state:
                        subtree.append(b)
                    else:
                        subtree.append(b2)
                        subtree_progressed = True
                _admit_subtree(
                    subtree, made_progress=subtree_progressed,
                    source_branch=bi)
            new_branches.sort(key=lambda b: -len(b.mapping))
            # Cross-branch dedupe is intentionally exact.  Coupled
            # automorphic branches remain separate until their transporter is
            # represented analytically; an orbit/event signature is not a
            # proof that their mapping families are equal.
            seen = {}
            uniq = []
            for b in new_branches:
                state_sig = _branch_signature(b)
                kept = seen.get(state_sig)
                if kept is not None:
                    _merge_equivalent_paths(kept, b)
                    continue
                seen[state_sig] = b
                uniq.append(b)
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


def _chirality_violations(mapping, coords_R, coords_P,
                          broken, formed, elements_R, min_deg=4):
    """Count spectator-stereocenter chirality flips.
    Spectator = atom whose 4+ mapped neighbors are all preserved
    (no incident broken/formed bond). Sign = sign(det) of first-3
    neighbor displacement vectors against the 4th. Skip near-coplanar."""
    bad_atoms = set()
    for (i, j, _, _) in broken:
        bad_atoms.add(i); bad_atoms.add(j)
    inv = {v: k for k, v in mapping.items()}
    for (ip, jp, _, _) in formed:
        if ip in inv: bad_atoms.add(inv[ip])
        if jp in inv: bad_atoms.add(inv[jp])

    # Build neighbor lists from R-coords (assume bond-graph already filtered)
    # We use "any atom within reasonable bond distance" — use the bond-graph
    # caller passes us implicitly through coords ordering. For correctness
    # we re-derive from coords by looking up close neighbors.
    # Simpler: just take all atoms within 1.9 Å in R as "neighbors".
    n_R = len(elements_R)
    diff = coords_R[:, None, :] - coords_R[None, :, :]
    dist = np.linalg.norm(diff, axis=2)
    np.fill_diagonal(dist, np.inf)
    violations = 0
    checked = 0
    for u in range(n_R):
        if u in bad_atoms: continue
        if u not in mapping: continue
        nbrs = np.where(dist[u] < 1.9)[0].tolist()
        # all must be mapped + spectator
        if len(nbrs) < min_deg: continue
        if any(nb not in mapping for nb in nbrs): continue
        if any(nb in bad_atoms for nb in nbrs): continue
        # take first 4 neighbors
        nb4 = nbrs[:4]
        # signed volume in R-frame
        ru = coords_R[u]
        v1 = coords_R[nb4[1]] - coords_R[nb4[0]]
        v2 = coords_R[nb4[2]] - coords_R[nb4[0]]
        v3 = coords_R[nb4[3]] - coords_R[nb4[0]]
        det_R = float(np.linalg.det(np.stack([v1, v2, v3])))
        # corresponding atoms in P-frame
        pu = coords_P[mapping[u]]
        pn = [coords_P[mapping[nb]] for nb in nb4]
        v1p = pn[1] - pn[0]
        v2p = pn[2] - pn[0]
        v3p = pn[3] - pn[0]
        det_P = float(np.linalg.det(np.stack([v1p, v2p, v3p])))
        if abs(det_R) < 0.05 or abs(det_P) < 0.05:
            continue  # near-coplanar; sign unstable
        checked += 1
        if np.sign(det_R) != np.sign(det_P):
            violations += 1
    return violations


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

"""Branch scheduling, mechanism-state dedupe, and symmetry repair."""
from __future__ import annotations

import itertools
import random
from collections import defaultdict

import numpy as np

from ..frag import classify_bonds
from ..growth import grow_island
from ..matcher import (
    _boundary_signature,
    _edge_wbo,
    _nauty_orbits,
    _orbit_id,
    _wbo_bucket,
)

SYM_REPAIR_MAX_EVALS = 20000


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


def _orbit_pair(a, b, orbits):
    oa, ob = _orbit_id(orbits, a), _orbit_id(orbits, b)
    return (oa, ob) if oa <= ob else (ob, oa)


def _chemistry_orbit_signature(mapping, g_R, g_P, r_orbits=None, p_orbits=None,
                               dwbo_threshold=0.5):
    """Broken/formed bond signature in joint R-orbit/P-orbit space."""
    br_pairs = []
    fm_pairs = []
    mapped = sorted(mapping)
    for i, u in enumerate(mapped):
        for v in mapped[i + 1:]:
            pu, pv = mapping[u], mapping[v]
            wR = _edge_wbo(g_R, u, v)
            wP = _edge_wbo(g_P, pu, pv)
            if wR - wP >= dwbo_threshold:
                br_pairs.append((
                    _orbit_pair(u, v, r_orbits),
                    _orbit_pair(pu, pv, p_orbits),
                ))
            elif wP - wR >= dwbo_threshold:
                fm_pairs.append((
                    _orbit_pair(u, v, r_orbits),
                    _orbit_pair(pu, pv, p_orbits),
                ))
    return tuple(sorted(br_pairs)), tuple(sorted(fm_pairs))


def _alignment_state_signature(mapping, deferred_edges, g_R, g_P,
                               r_orbits=None, p_orbits=None, core_R=()):
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
        p_orbits=p_orbits, locked_mapping=mapping)
    return tuple(fixed), tuple(sorted(internal_pairs)), boundary


def _mapping_change_score(mapping, wbo_R, wbo_P, dwbo_threshold=0.5):
    broken, formed, _, _ = classify_bonds(
        mapping, wbo_R, wbo_P, dwbo_threshold=dwbo_threshold)
    delta = 0.0
    for a, b, wR, wP in broken:
        delta += abs(float(wR or 0.0) - float(wP or 0.0))
    for a, b, wR, wP in formed:
        delta += abs(float(wP or 0.0) - float(wR or 0.0))
    return (len(broken) + len(formed), round(delta, 12))


def symmetry_repair_mapping(mapping, wbo_R, wbo_P, g_R, g_P, p_orbits,
                            dwbo_threshold=0.5, bond_floor=0.2,
                            min_changes=5, full_permutation_size=6,
                            max_evals=SYM_REPAIR_MAX_EVALS,
                            return_stats=False):
    """Choose the best concrete realization inside product symmetry orbits.

    The matcher intentionally compresses high-symmetry choices and keeps one
    witness.  A witness is not chemistry: if changed-bond endpoints sit inside
    an unresolved P orbit, another realization of the same compressed match may
    remove fake broken/formed bonds.  This local search swaps
    only atoms that already occupy the same product orbit and element.
    """
    if not mapping or p_orbits is None:
        return (mapping, {'enabled': False}) if return_stats else mapping
    mapping0 = dict(mapping)
    base_broken, base_formed, _, _ = classify_bonds(
        mapping0, wbo_R, wbo_P, dwbo_threshold=dwbo_threshold)
    base_changes = len(base_broken) + len(base_formed)
    stats = {
        'enabled': True,
        'base_changes': int(base_changes),
        'repaired': False,
        'groups': [],
        'evaluated': 0,
        'capped': False,
    }
    if base_changes < min_changes:
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

    local_pairs = []
    nR = wbo_R.shape[0]
    for i in range(nR):
        for j in range(i + 1, nR):
            if i not in mapping0 or j not in mapping0:
                continue
            if i not in local and j not in local:
                continue
            local_pairs.append((i, j, float(wbo_R[i, j])))

    def local_score(m):
        changed = 0
        delta = 0.0
        for i, j, wR in local_pairs:
            wP = float(wbo_P[m[i], m[j]])
            if wR - wP >= dwbo_threshold:
                changed += 1
                delta += wR - wP
            elif wP - wR >= dwbo_threshold:
                changed += 1
                delta += wP - wR
            elif wR >= bond_floor or wP >= bond_floor:
                delta += abs(wR - wP) * 0.01
        return (changed, round(delta, 12))

    current = dict(mapping0)
    current_score = local_score(current)
    for key, rs in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
        stats['groups'].append({
            'element': key[0],
            'orbit': int(key[1]),
            'size': len(rs),
            'atoms': [int(r) for r in rs],
        })
    max_passes = max(4, 2 * sum(len(rs) for rs in groups.values()))
    ordered_groups = [rs for _, rs in sorted(
        groups.items(), key=lambda item: (-len(item[1]), item[0]))]
    for _ in range(max_passes):
        improved = False
        for rs in ordered_groups:
            if stats['evaluated'] >= max_evals:
                stats['capped'] = True
                break
            best_score = current_score
            best_map = None
            if len(rs) <= full_permutation_size:
                targets = [current[r] for r in rs]
                for perm in itertools.permutations(targets):
                    if stats['evaluated'] >= max_evals:
                        stats['capped'] = True
                        break
                    candidate = dict(current)
                    for r, p in zip(rs, perm):
                        candidate[r] = p
                    score = local_score(candidate)
                    stats['evaluated'] += 1
                    if score < best_score:
                        best_score = score
                        best_map = candidate
            else:
                for a, b in itertools.combinations(rs, 2):
                    if stats['evaluated'] >= max_evals:
                        stats['capped'] = True
                        break
                    candidate = dict(current)
                    candidate[a], candidate[b] = candidate[b], candidate[a]
                    score = local_score(candidate)
                    stats['evaluated'] += 1
                    if score < best_score:
                        best_score = score
                        best_map = candidate
            if best_map is not None:
                current = best_map
                current_score = best_score
                improved = True
        if stats['capped']:
            break
        if not improved:
            break

    best = current
    best_score = _mapping_change_score(best, wbo_R, wbo_P,
                                       dwbo_threshold=dwbo_threshold)
    base_score = _mapping_change_score(mapping0, wbo_R, wbo_P,
                                       dwbo_threshold=dwbo_threshold)
    stats['best_changes'] = int(best_score[0])
    if best_score < base_score:
        stats['repaired'] = True
        return (best, stats) if return_stats else best
    return (mapping0, stats) if return_stats else mapping0


def find_islands(g_R, g_P, seed_order,
                 graph_floor=0.2, iso_tol=1.0,
                 dwbo_threshold=0.5, symmetry_wbo_tol=0.2,
                 max_branches=1_000_000, events=None,
                 orbit_dedup=True, core_R=None,
                 stop_when_core_mapped=False,
                 p_orbits=None, r_orbits=None,
                 abort_on_branch_cap=False):
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
    """
    if orbit_dedup:
        if p_orbits is None:
            p_orbits = _nauty_orbits(g_P, wbo_tol=symmetry_wbo_tol)
        if r_orbits is None:
            r_orbits = _nauty_orbits(g_R, wbo_tol=symmetry_wbo_tol)
    else:
        p_orbits = None
        r_orbits = None
    core_R = tuple(sorted(set(core_R or ())))
    core_R_set = set(core_R)
    seed_order = list(seed_order)
    if core_R:
        seed_order = ([s for s in seed_order if s in core_R_set] +
                      [s for s in seed_order if s not in core_R_set])
    branches = [_Branch()]
    progressed = True
    pass_no = 0

    def _core_complete(mapping):
        return core_R and all(r in mapping for r in core_R)

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
            p_orbits=p_orbits, locked_mapping=mapping)
        sig = (
            'mechanism_state',
            core_key,
            _chemistry_orbit_signature(
                mapping, g_R, g_P, r_orbits, p_orbits,
                dwbo_threshold=dwbo_threshold),
            deferred_boundary,
        )
        signature_cache[cache_key] = sig
        return sig

    def _branch_signature(branch):
        if branch._signature_cache is None:
            branch._signature_cache = _mapping_signature(
                branch.mapping, branch.deferred_edges)
        return branch._signature_cache

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
                if seed in b.mapping:
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
                                    prior_deferred_edges=b.deferred_edges)
                if not isos:
                    _append_pending(b)
                    continue
                # Dedup island results by weighted alignment state before
                # forking branches.  This is still pre-mechanism dedupe: the
                # one-hop/deferred boundary is part of the key so future
                # distinguishability is not erased.
                seen_state = {}
                for iso in isos:
                    full_m = dict(b.mapping); full_m.update(iso)
                    full_deferred = set(b.deferred_edges)
                    full_deferred.update(getattr(iso, 'deferred_edges', ()))
                    state_key = _mapping_signature(full_m, full_deferred)
                    if state_key not in seen_state:
                        seen_state[state_key] = iso
                deduped_isos = list(seen_state.values())
                for ii, iso in enumerate(deduped_isos):
                    if len(new_branches) >= max_branches:
                        _hit_branch_cap(len(new_branches), seed, 'island_fork')
                        break
                    b2 = b.fork()
                    b2.commit(iso, g_R,
                              events=events if (bi == 0 and ii == 0) else None)
                    if _append_pending(b2):
                        progressed = True
            new_branches.sort(key=lambda b: -len(b.mapping))
            # Cross-branch dedup uses the same mechanism-state key.  It
            # collapses orbit-equivalent spectator permutations but keeps
            # states separated when their deferred one-hop boundary differs.
            seen = set()
            uniq = []
            for b in new_branches:
                state_sig = _branch_signature(b)
                if state_sig in seen:
                    continue
                seen.add(state_sig)
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


def _generate_seed_orders(g_R, n_trials, rng_seed=42):
    """Generate at most `n_trials` deterministic seed orderings.

    Heavy atoms are used as anchors first in graph order. Hydrogens are never
    used as explicit first seeds because they have few connectivity
    constraints. If more trials are requested than there are heavy atoms, pad
    with full random shuffles.
    """
    nodes = list(g_R.nodes())
    heavy = [n for n in nodes if g_R.nodes[n].get('element') != 'H']
    rng = random.Random(rng_seed)
    orders = []
    n_trials = max(0, int(n_trials))
    for h in heavy[:n_trials]:
        rest = [x for x in nodes if x != h]
        rng.shuffle(rest)
        orders.append([h] + rest)
    while len(orders) < n_trials:
        perm = list(nodes); rng.shuffle(perm); orders.append(perm)
    return orders

"""Mechanism-local R->TS/IG core mapping."""
from __future__ import annotations

from .sweep import _canon_pair, _core_mapping_key, _pool_add


def _core_edge_match_ok(wboR, wboT, r, t, u, v, *,
                        edge_floor=0.2, iso_tol=1.0):
    if wboR[r, u] < 0.2:
        return True
    wT = float(wboT[t, v])
    return wT >= edge_floor and abs(float(wboR[r, u]) - wT) <= iso_tol


def _core_edges_preserved(mapping, wboR, wboT, core_R, *,
                          edge_floor=0.2, iso_tol=1.0):
    """Validate a TS/IG core mapping against reactant-side core identity."""
    core = tuple(sorted(set(core_R or ())))
    for idx, a in enumerate(core):
        if a not in mapping:
            return False
        for b in core[idx + 1:]:
            if b not in mapping:
                return False
            if wboR[a, b] < 0.2:
                continue
            wT = float(wboT[mapping[a], mapping[b]])
            if wT < edge_floor:
                return False
            if abs(float(wboR[a, b]) - wT) > iso_tol:
                return False
    return True


def ts_core_pool(elR, wboR, elT, wboT, core_R, *,
                 broken_R=None, formed_R=None,
                 edge_floor=0.2, iso_tol=1.0,
                 max_candidates=20000):
    """Enumerate mechanism-local R->TS/IG core alternatives.

    Only mechanism core atoms are assigned.  Exact core mappings are preserved
    for scoring because a symmetry choice that touches the core can change
    beta/rho/kappa; spectator choices remain compressed/unmaterialized.
    """
    core_R = tuple(sorted(set(core_R or ())))
    if not core_R:
        return {}
    pool = {}
    core_set = set(core_R)
    core_edges = {
        r: [u for u in core_R if u != r and wboR[r, u] >= 0.2]
        for r in core_R
    }
    reactive_pairs = {
        _canon_pair(int(a), int(b))
        for a, b in list(broken_R or []) + list(formed_R or [])
        if a in core_set and b in core_set
    }

    domains = {
        r: [t for t, e in enumerate(elT) if e == elR[r]]
        for r in core_R
    }

    changed = True
    while changed:
        changed = False
        for r in core_R:
            kept = []
            for t in domains[r]:
                ok = True
                for u in core_edges[r]:
                    if not any(
                        v != t and _core_edge_match_ok(
                            wboR, wboT, r, t, u, v,
                            edge_floor=edge_floor, iso_tol=iso_tol)
                        for v in domains[u]
                    ):
                        ok = False
                        break
                if ok:
                    kept.append(t)
            if len(kept) != len(domains[r]):
                domains[r] = kept
                changed = True

    def reactive_hint(r, t, mapping):
        score = 0.0
        for u, v in mapping.items():
            if _canon_pair(r, u) in reactive_pairs:
                score += float(wboT[t, v])
        return score

    def feasible_values(r, mapping, used_T):
        vals = []
        for t in domains[r]:
            if t in used_T:
                continue
            ok = True
            for u, v in mapping.items():
                if not _core_edge_match_ok(
                        wboR, wboT, r, t, u, v,
                        edge_floor=edge_floor, iso_tol=iso_tol):
                    ok = False
                    break
            if not ok:
                continue
            for u in core_R:
                if u in mapping or u == r or wboR[r, u] < 0.2:
                    continue
                if not any(
                    v not in used_T
                    and v != t
                    and _core_edge_match_ok(
                        wboR, wboT, r, t, u, v,
                        edge_floor=edge_floor, iso_tol=iso_tol)
                    for v in domains[u]
                ):
                    ok = False
                    break
            if ok:
                vals.append(t)
        vals.sort(key=lambda t: (-reactive_hint(r, t, mapping), t))
        return vals

    def backtrack(mapping, used_T):
        if len(pool) >= max_candidates:
            return
        if len(mapping) == len(core_R):
            if _core_edges_preserved(
                    mapping, wboR, wboT, core_R,
                    edge_floor=edge_floor, iso_tol=iso_tol):
                _pool_add(pool, _core_mapping_key(mapping, core_R),
                          dict(mapping), ())
            return
        remaining = [r for r in core_R if r not in mapping]
        ranked = []
        for r in remaining:
            vals = feasible_values(r, mapping, used_T)
            ranked.append((len(vals), -len(core_edges[r]), r, vals))
        ranked.sort()
        n_vals, _, r, vals = ranked[0]
        if n_vals == 0:
            return
        for t in vals:
            mapping[r] = t
            used_T.add(t)
            backtrack(mapping, used_T)
            used_T.remove(t)
            del mapping[r]

    backtrack({}, set())
    return pool

"""Priority-queue frontier utilities for fragment growth."""
from __future__ import annotations

import heapq

from ..matcher import _cand_has_open_choice, _cand_map


def _set_unique(cands):
    """True iff there is at most one fully resolved distinct mapping.

    Earlier this used `frozenset(c.values())` as the equivalence key,
    which collapses every same-value-set witness to one, silently dropping
    symmetric variants like (R87->IG87, R88->IG88) vs
    (R87->IG88, R88->IG87).  The right equivalence for a resolved candidate is
    the witness pair set itself.  Compressed symmetry candidates with an open
    target pool are not lockable yet; they must either absorb enough context or
    reach saturation first."""
    if not cands:
        return False
    if any(_cand_has_open_choice(c) for c in cands):
        return False
    if len(cands) == 1:
        return True
    sig0 = tuple(sorted(_cand_map(cands[0]).items()))
    return all(tuple(sorted(_cand_map(c).items())) == sig0 for c in cands[1:])


def _push_edges_from(heap, used_edges, g_R, atom, fragment, graph_floor):
    """Push all not-yet-seen outgoing traversal edges of `atom` into the heap."""
    for nb in g_R.neighbors(atom):
        if nb in fragment:
            continue
        edge = frozenset({atom, nb})
        if edge in used_edges:
            continue
        w = g_R[atom][nb]['wbo']
        if w >= graph_floor:
            heapq.heappush(heap, (-w, atom, nb))


def _frontier_boundary_edges(g_R, fragment, graph_floor):
    """One-hop live frontier edges that may still distinguish symmetry."""
    fragment = set(fragment)
    edges = set()
    for atom in fragment:
        for nb in g_R.neighbors(atom):
            if nb in fragment:
                continue
            if g_R[atom][nb]['wbo'] >= graph_floor:
                edges.add(tuple(sorted((atom, nb))))
    return edges

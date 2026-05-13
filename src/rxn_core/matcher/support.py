"""Support-witness search inside compressed symmetry blocks."""
from __future__ import annotations

from collections import defaultdict

from .primitives import SYM_SUPPORT_MAX_STATES, _edge_wbo, _growth_edge_supported, _orbit_id, _wbo_bucket
from .state import _SymBlock, _SymCand, _sym_block_indexes


def _r_compatible_with_block(cand, block_idx, n, fragment, g_R, r_orbits):
    if not isinstance(cand, _SymCand):
        return False
    b = cand.blocks[block_idx]
    if not b.extendable or b.complete:
        return False
    if r_orbits is not None:
        n_orbit = r_orbits[n]
        if any(r_orbits[r] != n_orbit for r in b.r_atoms):
            return False
    block_r = set(b.r_atoms)
    outside = sorted(r for r in fragment if r not in block_r and r in cand.mapping)

    def rel_sig(r):
        out = []
        for x in outside:
            if g_R.has_edge(r, x):
                out.append((x, _wbo_bucket(_edge_wbo(g_R, r, x))))
        return tuple(out)

    n_sig = rel_sig(n)
    return all(rel_sig(r) == n_sig for r in b.r_atoms)


def _refine_sym_assignments(cand, assignments):
    """Freeze selected block assignments and shrink the remaining pools.

    This is used when a new atom is valid only under a particular assignment
    inside an existing symmetry block.  Keeping the old block would imply the
    new fixed atom is independent of that block, which is false for cases like
    R37=P38 <=> R104=P119 vs R37=P36 <=> R104=P123.
    """
    if not isinstance(cand, _SymCand):
        m = dict(cand)
        used = {p for r, p in m.items() if r not in assignments}
        for r, p in assignments.items():
            if p in used:
                return None
            m[r] = p
            used.add(p)
        return m

    assignments = dict(assignments)
    if len(set(assignments.values())) != len(assignments):
        return None

    r_to_block, _ = _sym_block_indexes(cand)
    block_r = set(r_to_block)
    m = {r: p for r, p in cand.mapping.items() if r not in block_r}
    used = set(m.values())
    for r, p in sorted(assignments.items()):
        if p in used and m.get(r) != p:
            return None
        m[r] = p
        used.add(p)

    new_blocks = []
    for block in cand.blocks:
        remaining_r = tuple(r for r in block.r_atoms if r not in assignments)
        remaining_p = tuple(p for p in block.p_atoms
                            if p not in assignments.values())
        if not remaining_r:
            continue
        if len(remaining_r) > len(remaining_p):
            return None
        for r in remaining_r:
            p = cand.mapping.get(r)
            if p in remaining_p and p not in used:
                m[r] = p
                used.add(p)
        new_blocks.append(_SymBlock(remaining_r, remaining_p,
                                    extendable=block.extendable))
    exact_fixed = set(cand.exact_fixed)
    try:
        return _SymCand(m, tuple(new_blocks), exact_fixed=exact_fixed,
                        multiplicity=cand.multiplicity,
                        alternates=cand.alternates)
    except ValueError:
        return None


def _support_witness_for_value(cand, n, v_n, bonded_in_frag, r_wbos,
                               g_P, iso_tol, join_block_idx=None,
                               strict_r_wbos=None,
                               max_states=SYM_SUPPORT_MAX_STATES):
    """Find a concrete witness that supports R[n] -> P[v_n].

    `_SymCand` blocks represent a pool of legal P atoms, but the stored
    mapping is only one witness.  Extension must therefore ask whether some
    assignment inside each unresolved block can satisfy the new edge pattern,
    not whether the current witness happens to satisfy it.
    """
    strict_by_r = dict(strict_r_wbos or ())
    graph_floor = float(g_P.graph.get("bond_cut", 0.2))

    def _pair_ok(r_atom, w_R, p_atom):
        w_P = _edge_wbo(g_P, p_atom, v_n)
        if r_atom in strict_by_r:
            return _growth_edge_supported(
                strict_by_r[r_atom], w_P, iso_tol, graph_floor)
        return _growth_edge_supported(w_R, w_P, iso_tol, graph_floor)

    if not isinstance(cand, _SymCand):
        used = set(cand.values())
        if v_n in used and cand.get(n) != v_n:
            return None
        support = {}
        for u, w in r_wbos:
            if u not in cand:
                return None
            p_u = cand[u]
            if p_u == v_n:
                return None
            if not _pair_ok(u, w, p_u):
                return None
        return support

    r_to_block, p_to_block = _sym_block_indexes(cand)
    block_r = set(r_to_block)
    fixed_used = {p for r, p in cand.mapping.items() if r not in block_r}
    if v_n in fixed_used and cand.get(n) != v_n:
        return None
    value_block = p_to_block.get(v_n)
    if value_block is not None and value_block != join_block_idx:
        return None
    if join_block_idx is not None:
        block = cand.blocks[join_block_idx]
        if v_n not in block.p_atoms or not block.open:
            return None

    by_block = defaultdict(list)
    support = {}
    w_by_r = {u: w for u, w in r_wbos}
    for u in bonded_in_frag:
        if u not in cand.mapping:
            return None
        bidx = r_to_block.get(u)
        if bidx is None:
            p_u = cand.mapping[u]
            if p_u == v_n:
                return None
            if not _pair_ok(u, w_by_r[u], p_u):
                return None
            continue
        by_block[bidx].append((u, w_by_r[u]))

    for bidx, items in by_block.items():
        block = cand.blocks[bidx]
        reserved = set()
        if bidx == join_block_idx:
            reserved.add(v_n)
        domains = []
        available_pool = [p for p in block.p_atoms if p not in reserved]
        for u, w in items:
            vals = []
            for p in available_pool:
                if _pair_ok(u, w, p):
                    vals.append(p)
            if not vals:
                return None
            domains.append((u, vals))
        if all(set(vals) == set(available_pool) for _, vals in domains):
            continue
        domains.sort(key=lambda item: (len(item[1]), item[0]))
        chosen = {}
        used = set(reserved)
        states = 0

        def backtrack(pos):
            nonlocal states
            states += 1
            if states > max_states:
                return False
            if pos == len(domains):
                return True
            u, vals = domains[pos]
            for p in vals:
                if p in used:
                    continue
                used.add(p)
                chosen[u] = p
                if backtrack(pos + 1):
                    return True
                chosen.pop(u, None)
                used.remove(p)
            return False

        if not backtrack(0):
            return None
        support.update(chosen)
    return support


def _force_sym_value(cand, r, p, fragment, g_R, r_orbits, p_orbits):
    if not isinstance(cand, _SymCand):
        if r in cand:
            return cand if cand[r] == p else None
        if p in cand.values():
            return None
        nc = dict(cand)
        nc[r] = p
        return nc
    r_to_block, p_to_block = _sym_block_indexes(cand)
    if r in cand.mapping:
        current = cand.mapping[r]
        if current == p:
            return cand
        block_idx = r_to_block.get(r)
        if block_idx is None or p not in cand.blocks[block_idx].p_atoms:
            return None
        assignment = {r: p}
        for other in cand.blocks[block_idx].r_atoms:
            if other != r and cand.mapping.get(other) == p:
                assignment[other] = current
                break
        return _refine_sym_assignments(cand, assignment)
    block_idx = p_to_block.get(p)
    if block_idx is not None:
        if not _r_compatible_with_block(cand, block_idx, r, fragment,
                                        g_R, r_orbits):
            return None
        nc = cand.with_extended_block(block_idx, r)
        return _refine_sym_assignments(nc, {r: p}) if nc is not None else None
    block_r = set(r_to_block)
    fixed_used = {v for rr, v in cand.mapping.items() if rr not in block_r}
    if p in fixed_used:
        return None
    return cand.with_fixed(r, p)

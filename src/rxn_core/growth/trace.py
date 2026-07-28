"""Trace/event helpers for priority-queue fragment growth.

These helpers format state for HTML/debug traces. They should not decide
whether a match is valid; validity lives in ``rxn_core.matcher`` and the
island-growth loop in ``growth.island``.
"""
from __future__ import annotations

from collections import Counter, defaultdict

from ..matcher import (
    _SymCand,
    _cand_map,
    _cand_possible_p_atoms,
    _edge_wbo,
    _growth_edge_supported,
    as_node_match_policy,
    _support_witness_for_value,
    _sym_block_assignment_expr,
    _sym_block_indexes,
)


def cands_sample(cands, k=10):
    return [{int(a): int(b) for a, b in _cand_map(c).items()}
            for c in cands[:k]]


def cand_possible_values(cand, r):
    if isinstance(cand, _SymCand):
        for block in cand.blocks:
            if r in block.r_atoms:
                return set(block.p_atoms)
    cm = _cand_map(cand)
    return {cm[r]} if r in cm else set()


def cand_assignment_expr(cand):
    if not isinstance(cand, _SymCand):
        return '1'
    factors = []
    if cand.multiplicity != 1:
        factors.append(str(cand.multiplicity))
    for block in cand.blocks:
        expr = _sym_block_assignment_expr(block)
        if expr != '1':
            factors.append(expr)
    return ' * '.join(factors) if factors else '1'


def represented_assignment_expr(cands, max_terms=6):
    counts = Counter(cand_assignment_expr(c) for c in cands)
    if not counts:
        return '0'
    terms = []
    for expr, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:max_terms]:
        if expr == '1':
            terms.append(str(count))
        elif count == 1:
            terms.append(expr)
        else:
            terms.append(f'{count}*({expr})')
    remaining = max(0, len(counts) - max_terms)
    if remaining:
        terms.append(f'... + {remaining} expression groups')
    return ' + '.join(terms)


def cands_pattern_sample(cands, k=5):
    out = []
    for cand in cands[:k]:
        item = {
            'witness': {int(a): int(b) for a, b in _cand_map(cand).items()},
            'blocks': [],
        }
        if isinstance(cand, _SymCand):
            item['multiplicity'] = int(cand.multiplicity)
            item['automorph_domains'] = len(cand.automorph_blocks)
            for block in cand.blocks:
                item['blocks'].append({
                    'r_atoms': [int(x) for x in block.r_atoms],
                    'p_atoms': [int(x) for x in block.p_atoms],
                    'extendable': bool(block.extendable),
                    'assignments': _sym_block_assignment_expr(block),
                })
        out.append(item)
    return out


def heap_snapshot(heap, used_edges, fragment, mapping, k=None):
    """Pending heap entries sorted by WBO desc, filtered to live entries."""
    peek = list(heap)
    peek.sort()
    live = [(w, uu, nn) for (w, uu, nn) in peek
            if frozenset({uu, nn}) not in used_edges and nn not in fragment]
    if k is not None:
        live = live[:k]
    return [{'frag_atom': int(uu), 'ext_atom': int(nn),
             'wbo': round(-w, 3),
             'ext_status': ('mapped' if nn in mapping else 'free')}
            for w, uu, nn in live]


def pool_by_frag_atom(heap, used_edges, fragment, mapping, g_R):
    """Live propagation pool grouped by fragment atom."""
    peek = list(heap)
    peek.sort()
    by_u = defaultdict(list)
    for w, uu, nn in peek:
        if frozenset({uu, nn}) in used_edges:
            continue
        if nn in fragment:
            continue
        by_u[int(uu)].append({
            'ext_atom': int(nn),
            'wbo': round(-w, 3),
            'ext_status': ('mapped' if nn in mapping else 'free'),
            'ext_element': g_R.nodes[nn]['element'],
        })
    return [{'frag_atom': int(u),
             'frag_element': g_R.nodes[u]['element'],
             'edges': sorted(by_u[u], key=lambda x: -x['wbo'])}
            for u in sorted(by_u.keys())]


def why_extend_failed(cands, fragment, n_atom, anchor_atom, anchor_wbo,
                      g_R, g_P, locked_mapping, iso_tol, node_policy=None):
    """Per-candidate explanation of why extension to n_atom failed."""
    node_policy = as_node_match_policy(node_policy)
    locked_p_atoms = set((locked_mapping or {}).values())
    graph_floor = float(g_P.graph.get("bond_cut", 0.2))
    bonded = sorted(u for u in fragment if g_R.has_edge(u, n_atom))
    r_wbos = [(u, _edge_wbo(g_R, u, n_atom)) for u in bonded]
    strict_r_wbos = {}
    if anchor_atom is not None and anchor_atom in fragment:
        strict_r_wbos[anchor_atom] = (
            _edge_wbo(g_R, anchor_atom, n_atom)
            if anchor_wbo is None else anchor_wbo
        )
    out = []
    for ci, raw_cand in enumerate(cands[:5]):
        cand = raw_cand if isinstance(raw_cand, _SymCand) else _SymCand(raw_cand)
        for vi, cand in enumerate((cand,)):
            cm = _cand_map(cand)
            used_p = _cand_possible_p_atoms(cand)
            v_set = {
                v for v in g_P.nodes()
                if v not in used_p and v not in locked_p_atoms
                and node_policy.compatible(g_R, n_atom, g_P, v)
            }
            tried = []
            for v in sorted(v_set)[:30]:
                why = []
                _, p_to_block = _sym_block_indexes(cand)
                join_idx = p_to_block.get(v)
                support = _support_witness_for_value(
                    cand, n_atom, v, bonded, r_wbos, g_P, iso_tol,
                    join_block_idx=join_idx,
                    strict_r_wbos=strict_r_wbos)
                if support is None:
                    bad = []
                    for u, w in r_wbos:
                        if u not in cm:
                            continue
                        wp = _edge_wbo(g_P, cm[u], v)
                        if u in strict_r_wbos:
                            if not _growth_edge_supported(
                                    strict_r_wbos[u], wp, iso_tol, graph_floor):
                                bad.append(f'R[{u}]: active edge needs target WBO>={graph_floor:.3f} and |{w:.3f}-{wp:.3f}|={abs(w-wp):.3f}<={iso_tol}')
                        elif not _growth_edge_supported(
                                w, wp, iso_tol, graph_floor):
                            bad.append(f'R[{u}]: active edge needs target WBO>={graph_floor:.3f} and |{w:.3f}-{wp:.3f}|={abs(w-wp):.3f}<={iso_tol}')
                    why.append('; '.join(bad[:5]) if bad else 'no block witness')
                tried.append({'v': int(v), 'rejected': bool(why),
                              'reason': '; '.join(why) if why else 'OK'})
            out.append({
                'cand_idx': ci,
                'variant_idx': vi,
                'variant_multiplicity': (
                    int(cand.multiplicity)
                    if isinstance(cand, _SymCand) else 1
                ),
                'cand_at_in_frag_neighbors': {int(u): int(cm[u]) for u in bonded},
                'candidate_v_count': len(v_set),
                'tried_v': tried,
            })
        if len(out) >= 12:
            break
    return out


def why_merge_failed(cands, fragment, n_atom, mapping, islands_R,
                     g_R, g_P, iso_tol):
    """Per-candidate explanation of why whole-island merge failed."""
    graph_floor = float(g_P.graph.get("bond_cut", 0.2))
    if islands_R is None or n_atom not in islands_R:
        island_atoms = [n_atom]
    else:
        target_iid = islands_R[n_atom]
        island_atoms = [r for r, k in islands_R.items() if k == target_iid]
    out = []
    for ci, cand in enumerate(cands[:5]):
        cm = _cand_map(cand)
        why = []
        used_p = set(cm.values())
        for r in island_atoms:
            p = mapping[r]
            if p in used_p and cm.get(r) != p:
                owner = [k for k, v in cm.items() if v == p][0]
                why.append(f'P[{p}] (image of R[{r}]) already used by R[{owner}]')
            if r in cm and cm[r] != p:
                why.append(f'R[{r}] in cand as P[{cm[r]}], conflicts with mapping P[{p}]')
        nc = dict(cm)
        for r, p in [(r, mapping[r]) for r in island_atoms if r not in nc]:
            nc[r] = p
        check_set = set(island_atoms) | fragment
        for r in island_atoms:
            for r2 in sorted(check_set):
                if r2 == r:
                    continue
                if r2 not in check_set or r2 not in nc:
                    continue
                if r >= r2 and r2 in island_atoms:
                    continue
                wR = _edge_wbo(g_R, r, r2)
                p, p2 = nc[r], nc[r2]
                wP = _edge_wbo(g_P, p, p2)
                if not _growth_edge_supported(wR, wP, iso_tol, graph_floor):
                    why.append(f'R[{r}]-R[{r2}]: target WBO must be >= {graph_floor:.3f}; |{wR:.3f}-{wP:.3f}|={abs(wR-wP):.3f} (P[{p}]-P[{p2}])')
        out.append({'cand_idx': ci, 'reasons': why[:8]})
    return out

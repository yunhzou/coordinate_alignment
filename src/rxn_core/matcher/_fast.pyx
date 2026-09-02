# cython: language_level=3, binding=True
"""Optional compiled kernels for the matcher's hottest leaf functions.

Every function here is a typed transcription of a pure-Python original in
``dedupe``, ``support`` and ``canonical``.  It returns exactly the objects the
original returns (same tuple structure, same ``int``/``float``/``bool`` types,
same ordering) so that grouping keys, sort orders and dedupe digests are
unchanged.  Inputs the transcription does not cover (non-integer state caps,
graphs without a WBO matrix on some paths) are delegated to the Python
original.  The module is imported only when ``RXN_CORE_FAST=1`` (see
``primitives._load_fast_kernels``); the pure-Python path stays the default.
"""
from libc.math cimport fabs

from cpython.dict cimport PyDict_GetItemWithError
from cpython.object cimport PyObject

from operator import itemgetter

from .policy import as_node_match_policy
from .primitives import SYM_SUPPORT_MAX_STATES, _wbo_bucket
from .state import _SymBlock, _SymCand, _cand_map


# --------------------------------------------------------------------------
# WBO access: the same list-of-lists view primitives._edge_wbo caches on the
# graph, built here when it does not exist yet.
# --------------------------------------------------------------------------

cdef object _wbo_rows(object g):
    """Return the ``_wbo_rows`` view of ``g`` or None when it has no matrix."""
    graph = g.graph
    mat = graph.get("wbo_matrix")
    if mat is None:
        return None
    cached = graph.get("_wbo_rows")
    if cached is None or cached[0] is not mat:
        rows = None
        if hasattr(mat, "tolist") and getattr(
                getattr(mat, "dtype", None), "kind", None) == "f":
            rows = mat.tolist()
        if rows is None:
            rows = [[float(value) for value in row] for row in mat]
        cached = (mat, rows)
        graph["_wbo_rows"] = cached
    return cached[1]


cdef inline object _edge_wbo(object rows, object g, object a, object b):
    """``primitives._edge_wbo`` given the graph's row view (or None)."""
    if a == b:
        return 0.0
    if rows is not None:
        return (<list>rows)[a][b]
    if g.has_edge(a, b):
        return float(g[a][b].get("wbo", 0.0))
    return 0.0


cdef inline object _orbit_bucket(bint has_buckets, object wbo_buckets,
                                 bint exact_dict, object zero_bucket,
                                 object rows, object g, object a, object b):
    """``orbits._orbit_wbo_bucket(orbits, a, b, _edge_wbo(g, a, b))``.

    The WBO value is only needed on the ``_wbo_bucket`` fallback, so it is
    read there; ``_wbo_rows`` has already been populated by the caller, which
    is the only side effect the eager read had.
    """
    cdef PyObject* hit
    if has_buckets:
        key = (a, b) if a <= b else (b, a)
        if exact_dict:
            hit = PyDict_GetItemWithError(<dict>wbo_buckets, key)
            if hit is not NULL:
                return <object>hit
        elif key in wbo_buckets:
            return wbo_buckets[key]
        if zero_bucket is not None:
            return zero_bucket
    return _wbo_bucket(_edge_wbo(rows, g, a, b))


# --------------------------------------------------------------------------
# (a) dedupe._p_relation_signature_from_parts
# --------------------------------------------------------------------------

def p_relation_signature_from_parts(cand, v, g_P, p_orbits, cm_items=None,
                                    blocks=None, node_policy=None,
                                    compact=False):
    """Compiled ``dedupe._p_relation_signature_from_parts``."""
    cdef list rel, block_rel, edge_wbos
    cdef dict r_by_p, neighbor_buckets
    cdef Py_ssize_t i
    cdef bint member
    if cm_items is None:
        cm_items = tuple(sorted(_cand_map(cand).items()))
    if blocks is None:
        blocks = cand.blocks if isinstance(cand, _SymCand) else ()
    node_policy = as_node_match_policy(node_policy)
    structural_zero = getattr(p_orbits, 'zero_bucket', None)

    rows = _wbo_rows(g_P)
    cdef bint has_buckets = hasattr(p_orbits, "wbo_buckets")
    wbo_buckets = p_orbits.wbo_buckets if has_buckets else None
    cdef bint exact_dict = has_buckets and type(wbo_buckets) is dict
    zero_bucket = structural_zero

    if compact and structural_zero is not None:
        active_cm = tuple((r, p) for r, p in cm_items if p != v)
        mapped_r = tuple(r for r, _ in active_cm)
        r_by_p = {p: r for r, p in active_cm}
        rel = []
        for p in g_P.neighbors(v):
            r = r_by_p.get(p)
            if r is None:
                continue
            bucket = _orbit_bucket(has_buckets, wbo_buckets, exact_dict,
                                   zero_bucket, rows, g_P, p, v)
            if bucket != structural_zero:
                rel.append((r, bucket))
        rel_out = ('sparse', mapped_r, tuple(sorted(rel)))

        neighbor_buckets = {}
        for p in g_P.neighbors(v):
            bucket = _orbit_bucket(has_buckets, wbo_buckets, exact_dict,
                                   zero_bucket, rows, g_P, p, v)
            if bucket != structural_zero:
                neighbor_buckets[p] = bucket
        block_rel = []
        i = 0
        for b in blocks:
            p_atoms = b.p_atoms
            member = v in p_atoms
            edge_wbos = []
            for p in p_atoms:
                if p != v and p in neighbor_buckets:
                    edge_wbos.append(neighbor_buckets[p])
            edge_wbos.sort()
            block_rel.append(
                (i, member, len(p_atoms) - (1 if member else 0),
                 tuple(edge_wbos)))
            i += 1
    else:
        rel = []
        for r, p in cm_items:
            if p == v:
                continue
            rel.append((r, _orbit_bucket(has_buckets, wbo_buckets, exact_dict,
                                         zero_bucket, rows, g_P, p, v)))
        rel_out = tuple(rel)
        block_rel = []
        i = 0
        for b in blocks:
            p_atoms = b.p_atoms
            edge_wbos = []
            for p in p_atoms:
                if p == v:
                    continue
                edge_wbos.append(_orbit_bucket(
                    has_buckets, wbo_buckets, exact_dict, zero_bucket,
                    rows, g_P, p, v))
            edge_wbos.sort()
            member = v in p_atoms
            block_rel.append((i, member, tuple(edge_wbos)))
            i += 1
    orbit_id = p_orbits[v] if p_orbits is not None else v
    return (node_policy.key(g_P, v), orbit_id, rel_out, tuple(block_rel))


# --------------------------------------------------------------------------
# (b) support._support_witness_for_value
# --------------------------------------------------------------------------

cdef inline bint _pair_ok(object rows, object g_P, object v_n,
                          dict strict_by_r, object r_atom, object w_R,
                          object p_atom, double iso_tol,
                          double graph_floor) except -1:
    """``support._pair_ok`` -> ``primitives._growth_edge_supported``."""
    cdef double w_P = _edge_wbo(rows, g_P, p_atom, v_n)
    cdef PyObject* hit = PyDict_GetItemWithError(strict_by_r, r_atom)
    if hit is not NULL:
        w_R = <object>hit
    cdef double wr = w_R
    return w_P >= graph_floor and fabs(wr - w_P) <= iso_tol


def _domain_sort_key(item):
    return (len(item[1]), item[0])


cdef bint _backtrack(list domains, Py_ssize_t pos, set used, dict chosen,
                     long long* states, long long max_states) except -1:
    states[0] += 1
    if states[0] > max_states:
        return False
    if pos == len(domains):
        return True
    u, vals = <tuple>domains[pos]
    for p in vals:
        if p in used:
            continue
        used.add(p)
        chosen[u] = p
        if _backtrack(domains, pos + 1, used, chosen, states, max_states):
            return True
        chosen.pop(u, None)
        used.remove(p)
    return False


def support_witness_for_value(cand, n, v_n, bonded_in_frag, r_wbos,
                              g_P, iso_tol, join_block_idx=None,
                              strict_r_wbos=None,
                              max_states=SYM_SUPPORT_MAX_STATES,
                              block_indexes=None):
    """Compiled ``support._support_witness_for_value``."""
    cdef long long max_states_c, states
    cdef dict strict_by_r, mapping, r_to_block, p_to_block, w_by_r
    cdef dict support, by_block, chosen
    cdef list order, items, domains, available_pool, vals
    cdef set used, reserved, pool_set
    cdef bint all_full
    cdef double tol, graph_floor
    if not isinstance(max_states, int):
        from .support import _support_witness_for_value_py
        return _support_witness_for_value_py(
            cand, n, v_n, bonded_in_frag, r_wbos, g_P, iso_tol,
            join_block_idx=join_block_idx, strict_r_wbos=strict_r_wbos,
            max_states=max_states, block_indexes=block_indexes)
    try:
        max_states_c = max_states
    except OverflowError:
        from .support import _support_witness_for_value_py
        return _support_witness_for_value_py(
            cand, n, v_n, bonded_in_frag, r_wbos, g_P, iso_tol,
            join_block_idx=join_block_idx, strict_r_wbos=strict_r_wbos,
            max_states=max_states, block_indexes=block_indexes)

    strict_by_r = dict(strict_r_wbos or ())
    graph_floor = float(g_P.graph.get("bond_cut", 0.2))
    tol = iso_tol
    rows = _wbo_rows(g_P)

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
            if not _pair_ok(rows, g_P, v_n, strict_by_r, u, w, p_u,
                            tol, graph_floor):
                return None
        return support

    if block_indexes is None:
        r_to_block = {}
        p_to_block = {}
        i = 0
        for block in cand.blocks:
            for r in block.r_atoms:
                r_to_block[r] = i
            for p in block.p_atoms:
                p_to_block[p] = i
            i += 1
    else:
        r_to_block, p_to_block = block_indexes
    mapping = cand.mapping
    # ``v_n in fixed_used`` with fixed_used = images of atoms outside blocks.
    for r, p in mapping.items():
        if p == v_n and r not in r_to_block:
            if mapping.get(n) != v_n:
                return None
            break
    value_block = p_to_block.get(v_n)
    if value_block is not None and value_block != join_block_idx:
        return None
    if join_block_idx is not None:
        block = cand.blocks[join_block_idx]
        if v_n not in block.p_atoms or not block.open:
            return None

    by_block = {}
    order = []
    support = {}
    w_by_r = {}
    for u, w in r_wbos:
        w_by_r[u] = w
    for u in bonded_in_frag:
        if u not in mapping:
            return None
        bidx = r_to_block.get(u)
        if bidx is None:
            p_u = mapping[u]
            if p_u == v_n:
                return None
            if not _pair_ok(rows, g_P, v_n, strict_by_r, u, w_by_r[u], p_u,
                            tol, graph_floor):
                return None
            continue
        items = <list>by_block.get(bidx)
        if items is None:
            items = []
            by_block[bidx] = items
            order.append(bidx)
        items.append((u, w_by_r[u]))

    for bidx in order:
        items = <list>by_block[bidx]
        block = cand.blocks[bidx]
        reserved = set()
        if bidx == join_block_idx:
            reserved.add(v_n)
        domains = []
        available_pool = [p for p in block.p_atoms if p not in reserved]
        for u, w in items:
            vals = []
            for p in available_pool:
                if _pair_ok(rows, g_P, v_n, strict_by_r, u, w, p,
                            tol, graph_floor):
                    vals.append(p)
            if not vals:
                return None
            domains.append((u, vals))
        pool_set = set(available_pool)
        all_full = True
        for _u, vals in domains:
            if set(vals) != pool_set:
                all_full = False
                break
        if all_full:
            continue
        domains.sort(key=_domain_sort_key)
        chosen = {}
        used = set(reserved)
        states = 0
        if not _backtrack(domains, 0, used, chosen, &states, max_states_c):
            return None
        support.update(chosen)
    return support


# --------------------------------------------------------------------------
# (c) canonical._CandidateAutomorphismCanonicalizer role kernels
# --------------------------------------------------------------------------

def candidate_roles(self, cand, *, group_domains=False):
    """Compiled ``_CandidateAutomorphismCanonicalizer._candidate_roles``."""
    cdef dict roles = {}
    cdef dict out
    cdef set block_r
    cdef list items
    cdef bint is_sym = isinstance(cand, _SymCand)
    if is_sym:
        mapping = cand.mapping
        blocks = cand.blocks
    else:
        mapping = dict(cand)
        blocks = ()

    block_r = set()
    for block in blocks:
        for r in block.r_atoms:
            block_r.add(r)
    if group_domains and is_sym:
        for block in cand.automorph_blocks:
            for r in block.r_atoms:
                block_r.add(r)
    # ``sorted(mapping.items())`` orders by the unique keys, so sorting the
    # keys alone visits the pairs in the same order.
    for r in sorted(mapping):
        if r not in block_r:
            p = mapping[r]
            key = int(p)
            items = <list>roles.get(key)
            if items is None:
                items = []
                roles[key] = items
            items.append(('mapped', int(r)))
    for block in blocks:
        block_role = (
            'pool', tuple(int(r) for r in block.r_atoms),
            bool(block.extendable),
        )
        for p in block.p_atoms:
            key = int(p)
            items = <list>roles.get(key)
            if items is None:
                items = []
                roles[key] = items
            items.append(block_role)
    if is_sym:
        for block in cand.automorph_blocks:
            group_role = (
                'automorph_domain', tuple(int(r) for r in block.r_atoms)
            )
            for p in block.p_atoms:
                key = int(p)
                items = <list>roles.get(key)
                if items is None:
                    items = []
                    roles[key] = items
                items.append(group_role)
    out = {}
    for p, items in roles.items():
        if len(items) == 1:
            out[p] = tuple(items)
        else:
            out[p] = tuple(sorted(items, key=repr))
    return out


def role_key_from_roles(self, candidate_roles, orbits):
    """Compiled ``_CandidateAutomorphismCanonicalizer.role_key_from_roles``."""
    cdef dict counts = {}
    cdef dict locked_roles = self.locked_roles
    cdef bint singleton = True
    cdef PyObject* hit
    sizes = getattr(orbits, '_orbit_sizes', None)
    if sizes is None:
        from collections import Counter
        sizes = Counter(orbits.values())
        orbits._orbit_sizes = sizes
    empty = ()
    for p, role in candidate_roles.items():
        orbit = orbits[p]
        if sizes[orbit] > 1:
            singleton = False
        hit = PyDict_GetItemWithError(locked_roles, p)
        locked = <object>hit if hit is not NULL else empty
        item = (orbit, locked, role)
        hit = PyDict_GetItemWithError(counts, item)
        if hit is NULL:
            counts[item] = 1
        else:
            counts[item] = <object>hit + 1
    return frozenset(counts.items()), singleton


def colored_vertices_from_roles(self, candidate_roles):
    """Compiled ``_CandidateAutomorphismCanonicalizer._colored_vertices_from_roles``."""
    cdef dict colors = {}
    cdef dict locked_roles = self.locked_roles
    cdef dict atom_base_color = self.atom_base_color
    cdef dict atom_index = self.atom_index
    cdef dict repr_cache = self.__dict__.setdefault('_color_repr_cache', {})
    cdef set vertices
    cdef list keyed
    cdef PyObject* hit
    empty = ()
    for p, vertex in atom_index.items():
        hit = PyDict_GetItemWithError(locked_roles, p)
        locked = <object>hit if hit is not NULL else empty
        role = (locked, candidate_roles.get(p, empty))
        color = ('atom', atom_base_color[vertex], role)
        hit = PyDict_GetItemWithError(colors, color)
        if hit is NULL:
            colors[color] = {vertex}
        else:
            (<set>hit).add(vertex)
    for color_index, edge_vertices in self.edge_color_classes:
        color = ('edge', color_index)
        hit = PyDict_GetItemWithError(colors, color)
        if hit is NULL:
            vertices = set()
            colors[color] = vertices
        else:
            vertices = <set>hit
        vertices.update(edge_vertices)
    keyed = []
    for color, vertices in colors.items():
        hit = PyDict_GetItemWithError(repr_cache, color)
        if hit is NULL:
            order_key = repr(color)
            repr_cache[color] = order_key
        else:
            order_key = <object>hit
        keyed.append((order_key, color, frozenset(vertices)))
    keyed.sort(key=itemgetter(0))
    return tuple((color, frozen) for _k, color, frozen in keyed)


# --------------------------------------------------------------------------
# (e) state._SymCand.__init__: witness completion inside symmetry blocks
# --------------------------------------------------------------------------

def symcand_init(self, mapping=None, blocks=(), exact_fixed=(),
                 multiplicity=1, automorph_blocks=()):
    """Compiled ``_SymCand.__init__``; same state, same ValueErrors."""
    cdef dict raw, m
    cdef set used, block_r, p_set, r_set
    cdef list available, missing, expanded
    blocks = tuple(blocks)
    block_r = set()
    for b in blocks:
        for r in b.r_atoms:
            block_r.add(r)
    raw = dict(mapping or {})
    m = {}
    for r, p in raw.items():
        if r not in block_r:
            m[r] = p
    used = set(m.values())
    for b in blocks:
        r_atoms = b.r_atoms
        p_atoms = b.p_atoms
        if len(r_atoms) > len(p_atoms):
            raise ValueError("symmetry block has more R atoms than P atoms")
        for r in r_atoms:
            if r not in raw:
                continue
            p = raw[r]
            if p not in p_atoms or p in used:
                raise ValueError(
                    "symmetry block witness conflicts with fixed atoms")
            m[r] = p
            used.add(p)
        available = [p for p in p_atoms if p not in used]
        missing = [r for r in r_atoms if r not in m]
        if len(available) < len(missing):
            raise ValueError(
                "symmetry block witness conflicts with fixed atoms")
        for r, p in zip(missing, available):
            m[r] = p
            used.add(p)
    self.mapping = m
    self.blocks = blocks
    self.exact_fixed = frozenset(exact_fixed)
    self.multiplicity = int(multiplicity)
    expanded = []
    for block in automorph_blocks:
        p_set = set(block.p_atoms)
        r_set = set(block.r_atoms)
        for r, p in m.items():
            if p in p_set:
                r_set.add(r)
        expanded.append(_SymBlock(tuple(r_set), tuple(p_set),
                                  extendable=False))
    self.automorph_blocks = tuple(expanded)


# --------------------------------------------------------------------------
# (d) dedupe._boundary_signature: per-pool target signature multiset
# --------------------------------------------------------------------------

def pool_target_signatures(pool, used_possible, inverse, mapped_r, blocks,
                           target_static, general):
    """Compiled ``dedupe._pool_target_signatures``."""
    cdef dict counts = {}
    cdef dict inverse_d = inverse
    cdef dict target_static_d = target_static
    cdef dict neighbor_buckets
    cdef list rel, block_rel, edge_wbos
    cdef Py_ssize_t i
    cdef bint member
    cdef PyObject* hit
    for v in pool:
        if v in used_possible:
            continue
        static = target_static_d.get(v)
        if static is None or v in inverse_d:
            sig = general(v)
        else:
            v_key, v_orbit, neighbours = static
            rel = []
            neighbor_buckets = {}
            for p, bucket in neighbours:
                neighbor_buckets[p] = bucket
                hit = PyDict_GetItemWithError(inverse_d, p)
                if hit is not NULL and <object>hit is not None:
                    rel.append((<object>hit, bucket))
            block_rel = []
            i = 0
            for b in blocks:
                p_atoms = b.p_atoms
                member = v in p_atoms
                edge_wbos = []
                for p in p_atoms:
                    if p != v:
                        hit = PyDict_GetItemWithError(neighbor_buckets, p)
                        if hit is not NULL:
                            edge_wbos.append(<object>hit)
                edge_wbos.sort()
                block_rel.append(
                    (i, member, len(p_atoms) - (1 if member else 0),
                     tuple(edge_wbos)))
                i += 1
            rel.sort()
            sig = (v_key, v_orbit, ('sparse', mapped_r, tuple(rel)),
                   tuple(block_rel))
        hit = PyDict_GetItemWithError(counts, sig)
        if hit is NULL:
            counts[sig] = 1
        else:
            counts[sig] = <object>hit + 1
    return frozenset(counts.items())

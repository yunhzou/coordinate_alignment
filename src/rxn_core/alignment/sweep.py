"""Sweep-cut mechanism discovery for WBO graph alignment.

R-P sweep cut is part of the core alignment algorithm: mechanism discovery
tries the no-cut graph and every one-edge R cut above a WBO floor, retains
unique analytical branch hierarchies, and groups them by exact canonical
broken/formed-event certificates.
"""
from __future__ import annotations

import copy
import json
import multiprocessing as mp
import time
from pathlib import Path

from ..frag import build_graph, classify_bonds, expand_mapping
from ..matcher import (
    _SymBlock,
    _atom_tuple_orbit,
    _nauty_atom_generators,
    _nauty_orbits,
    _sym_block_assignment_expr,
)
from ..matcher.orbits import _nauty_colored_wbo_graph
from .branch import (
    BranchLimitExceeded,
    _generate_seed_orders,
    _mapping_variation_blocks,
    find_islands,
    symmetry_repair_mapping,
)


def _canon_pair(a, b):
    return (a, b) if a <= b else (b, a)


def _strong_edges(wboR, cut_floor):
    n = int(wboR.shape[0])
    return [(i, j) for i in range(n) for j in range(i + 1, n)
            if float(wboR[i, j]) >= cut_floor]


def _orbit_bond_key(pairs, orbits, tag):
    return tuple(sorted(
        (tag, *_canon_pair(int(orbits[a]), int(orbits[b])))
        for a, b in pairs
    ))


def _core_mapping_key(mapping, core_R):
    return (
        tuple((int(r), int(mapping[r])) for r in sorted(core_R)),
        (),
    )


def _core_branch_record(branch, core_R, p_orbits=None):
    """Return the core-restricted compressed state for one aligned branch.

    This is intentionally not an exact core-map expansion.  It keeps the
    branch's set-to-set symmetry pools so TS/IG scoring can merge endpoint
    support before materializing automorphism variants.
    """
    core_R = tuple(sorted(set(int(r) for r in core_R)))
    core_set = set(core_R)
    mapping = {
        int(r): int(branch.mapping[r])
        for r in core_R
        if r in branch.mapping
    }
    if len(mapping) != len(core_R):
        return None

    # Core-map compression is based on exact target automorphism orbits from
    # pynauty, not on transient candidate-pool history.  Joint permutations
    # are validated later with strict target automorphism generators.
    orbit_members = {}
    for p, orbit in dict(p_orbits or {}).items():
        orbit_members.setdefault(int(orbit), []).append(int(p))
    core_by_orbit = {}
    for r, p in mapping.items():
        orbit = (p_orbits or {}).get(p)
        if orbit is not None:
            core_by_orbit.setdefault(int(orbit), []).append(int(r))
    blocks = []
    for orbit, r_atoms in sorted(core_by_orbit.items()):
        p_atoms = sorted(orbit_members.get(orbit, ()))
        if len(p_atoms) <= 1:
            continue
        blocks.append({
            'r_atoms': sorted(r_atoms),
            'p_atoms': p_atoms,
            'source': 'pynauty_target_orbit',
        })

    return {
        'mapping': mapping,
        'blocks': blocks,
        'dedup_count': 1,
    }


def _core_branch_record_key(record, core_R):
    core_R = tuple(sorted(set(int(r) for r in core_R)))
    mapping = {int(r): int(t) for r, t in record.get('mapping', {}).items()}
    blocks = []
    block_r = set()
    for block in record.get('blocks', ()):
        r_atoms = tuple(sorted(int(r) for r in block.get('r_atoms', ())))
        p_atoms = tuple(sorted(int(p) for p in block.get('p_atoms', ())))
        if not r_atoms:
            continue
        blocks.append((r_atoms, p_atoms))
        block_r.update(r_atoms)
    fixed = tuple(
        (int(r), int(mapping[r]))
        for r in core_R
        if r in mapping and r not in block_r
    )
    return fixed, tuple(sorted(blocks))


def _core_mapping_variants(branch, core_R, max_variants, *,
                           g_P=None, p_orbits=None):
    """Expand final candidate-carried core ambiguity with strict pynauty."""
    core_R = tuple(sorted(set(int(r) for r in core_R)))
    base = {
        int(r): int(branch.mapping[r])
        for r in core_R
        if r in branch.mapping
    }
    if len(base) != len(core_R):
        return []
    if g_P is None:
        return [base]

    blocks = []
    covered = set()
    final_symmetry = _branch_symmetry_record(branch)
    for block in final_symmetry.get('blocks', ()):
        source = block.get('source') or 'sym_block'
        if source not in {
                'sym_block', 'exact_automorph_group'}:
            continue
        r_atoms = tuple(sorted(
            int(r) for r in block.get('r_atoms', ())
            if int(r) in set(core_R) and int(r) not in covered
        ))
        if not r_atoms:
            continue
        p_atoms = tuple(sorted(int(p) for p in block.get('p_atoms', ())))
        if len(p_atoms) < len(r_atoms):
            continue
        blocks.append((r_atoms, frozenset(p_atoms)))
        covered.update(r_atoms)
    if not blocks:
        return [base]

    fixed = {r: base[r] for r in core_R if r not in covered}
    tag_parts = {}

    def add_tag(atom, tag):
        tag_parts.setdefault(int(atom), []).append(tag)

    for r, p in fixed.items():
        add_tag(p, ('fixed', int(r), int(p)))
    for block_index, (_r_atoms, p_atoms) in enumerate(blocks):
        for p in p_atoms:
            add_tag(p, ('block', int(block_index)))

    generators = _nauty_atom_generators(
        g_P,
        wbo_tol=float(getattr(p_orbits, 'wbo_tol', 0.2)),
        atom_color_tags={
            atom: tuple(tags) for atom, tags in tag_parts.items()
        },
    )
    seed = tuple(base[r] for r in core_R)
    states = _atom_tuple_orbit(seed, generators)
    out = []
    for state in states:
        mapping = {int(r): int(p) for r, p in zip(core_R, state)}
        if any(mapping[r] != p for r, p in fixed.items()):
            continue
        if any(mapping[r] not in p_atoms
               for r_atoms, p_atoms in blocks for r in r_atoms):
            continue
        out.append(mapping)
    if len(out) > int(max_variants):
        raise BranchLimitExceeded(
            max_variants,
            branch_count=len(out),
            stage='core_symmetry_expansion',
        )
    return out


class _MechanismEventCanonicalizer:
    """Reuse one full-R pynauty graph across concrete event sets."""

    def __init__(self, graph_R, wbo_tol):
        _nodes, self.atom_index, base, _buckets, _zero = (
            _nauty_colored_wbo_graph(
                graph_R, wbo_tol=float(wbo_tol)))
        self.base_adjacency = {
            int(vertex): set(neighbors)
            for vertex, neighbors in base._adjacency_dict.items()
        }
        self.base_colors = [set(cell) for cell in base._vertex_coloring]
        self.base_vertex_count = int(base.number_of_vertices)
        self.cache = {}

    def certificate(self, broken_pairs, formed_pairs):
        broken_pairs = tuple(sorted(tuple(sorted(map(int, pair)))
                                    for pair in broken_pairs))
        formed_pairs = tuple(sorted(tuple(sorted(map(int, pair)))
                                    for pair in formed_pairs))
        key = broken_pairs, formed_pairs
        cached = self.cache.get(key)
        if cached is not None:
            return cached
        import pynauty

        adjacency = {vertex: set(neighbors)
                     for vertex, neighbors in self.base_adjacency.items()}
        colors = [set(cell) for cell in self.base_colors]
        next_vertex = self.base_vertex_count
        type_markers = {}
        for event_name in ('broken', 'formed'):
            type_markers[event_name] = next_vertex
            adjacency[next_vertex] = set()
            colors.append({next_vertex})
            next_vertex += 1
        for event_name, pairs in (
                ('broken', broken_pairs), ('formed', formed_pairs)):
            event_vertices = set()
            for left, right in pairs:
                vertex = next_vertex
                next_vertex += 1
                a = self.atom_index[left]
                b = self.atom_index[right]
                adjacency.setdefault(a, set()).add(vertex)
                adjacency.setdefault(b, set()).add(vertex)
                marker = type_markers[event_name]
                adjacency[marker].add(vertex)
                adjacency[vertex] = {a, b, marker}
                event_vertices.add(vertex)
            if event_vertices:
                colors.append(event_vertices)
        graph = pynauty.Graph(
            next_vertex, directed=False,
            adjacency_dict={
                vertex: sorted(adjacency.get(vertex, ()))
                for vertex in range(next_vertex)
            },
            vertex_coloring=colors)
        certificate = pynauty.certificate(graph).hex()
        self.cache[key] = certificate
        return certificate


def _mechanism_signature(mapping, wboR, wboT, r_orbits, p_orbits,
                         dwbo_threshold=0.5,
                         elements_R=None, elements_P=None,
                         metal_dwbo_threshold=None,
                         g_R_full=None, symmetry_wbo_tol=0.2,
                         event_canonicalizer=None):
    """Exact symmetry-canonical mechanism-event certificate.

    Individual vertex-orbit IDs do not identify an edge orbit: two pairs can
    have endpoints in the same vertex orbits without one graph automorphism
    transporting the complete pair.  Canonicalize all event edges jointly on
    the full colored reactant graph instead.
    """
    broken, formed, _, _ = classify_bonds(
        mapping, wboR, wboT, dwbo_threshold=dwbo_threshold,
        elements_R=elements_R, elements_P=elements_P,
        metal_dwbo_threshold=metal_dwbo_threshold)
    inv = {v: k for k, v in mapping.items()}
    br_pairs = [(a, b) for (a, b, _, _) in broken]
    fm_r_pairs = []
    fm_p_pairs = []
    for a, b, _, _ in formed:
        if a in inv and b in inv:
            fm_r_pairs.append((inv[a], inv[b]))
        else:
            fm_p_pairs.append((a, b))
    if fm_p_pairs:
        # Complete R-P mappings should always pull formed edges into R.  Keep
        # an explicit failure instead of silently producing a weaker key.
        raise ValueError("mechanism signature received an unmapped formed edge")
    if g_R_full is None:
        # Compatibility for direct unit calls; production always supplies the
        # full reactant graph.
        br = _orbit_bond_key(br_pairs, r_orbits, 'R')
        fm = _orbit_bond_key(fm_r_pairs, r_orbits, 'R')
        return br, tuple(sorted(fm))

    event_canonicalizer = event_canonicalizer or (
        _MechanismEventCanonicalizer(
            g_R_full, wbo_tol=float(symmetry_wbo_tol)))
    certificate = event_canonicalizer.certificate(br_pairs, fm_r_pairs)
    return (
        tuple('broken' for _ in br_pairs),
        tuple('formed' for _ in fm_r_pairs),
        ('event_certificate_v1', certificate),
    )


def _branch_symmetry_record(branch, symmetry_fragments=None):
    fragments = []
    blocks = []
    p_to_r = {
        int(p): int(r)
        for r, p in getattr(branch, 'mapping', {}).items()
    }
    if symmetry_fragments is None:
        symmetry_fragments = getattr(branch, 'symmetry_fragments', ())
    for frag_index, fragment in enumerate(symmetry_fragments):
        record = {
            'fragment_index': int(frag_index),
            'island_idx': int(fragment.get('island_idx', frag_index)),
            'fragment': [int(r) for r in fragment.get('fragment', ())],
            'deferred_edges': [
                [int(a), int(b)]
                for a, b in fragment.get('deferred_edges', ())
            ],
            'symmetry': fragment.get('symmetry') or {},
        }
        fragments.append(record)
        symmetry = record['symmetry']
        for block_index, block in enumerate(symmetry.get('blocks') or ()):
            r_atoms = {
                int(r) for r in block.get('r_atoms', ())
            }
            p_atoms = [int(p) for p in block.get('p_atoms', ())]
            for p in p_atoms:
                if p in p_to_r:
                    r_atoms.add(p_to_r[p])
            r_atoms = sorted(r_atoms)
            if len(p_atoms) <= 1 and len(r_atoms) <= 1:
                continue
            normalized_block = _SymBlock(tuple(r_atoms), tuple(p_atoms),
                                         extendable=False)
            source = block.get('source')
            assignments = (
                block.get('assignments')
                if source in {'island_automorph', 'sym_block', 'interbranch',
                              'exact_automorph_group'}
                else _sym_block_assignment_expr(normalized_block)
            )
            item = {
                'fragment_index': int(frag_index),
                'block_index': int(block_index),
                'island_idx': record['island_idx'],
                'r_atoms': r_atoms,
                'p_atoms': p_atoms,
                'extendable': bool(block.get('extendable', False)),
                'open': bool(normalized_block.open),
                'assignments': assignments,
            }
            if source:
                item['source'] = source
            blocks.append(item)
    return {
        'rule': 'branch_symmetry_blocks',
        'fragments': fragments,
        'blocks': blocks,
    }


def _branch_analytical_derivations(branch, mapping):
    paths = getattr(branch, 'symmetry_paths', None)
    if not paths:
        return [(dict(mapping), _branch_symmetry_record(branch))]
    derivations = []
    for path in paths:
        symmetry = _branch_symmetry_record(
            branch, symmetry_fragments=path)
        derivations.append((dict(mapping), symmetry))
    return derivations


def _mapping_key(mapping):
    return tuple(
        (int(r), int(p))
        for r, p in sorted(dict(mapping).items())
    )


def _cut_record(cuts):
    return [
        [int(a), int(b)]
        for a, b in sorted(tuple(tuple(pair) for pair in cuts))
    ]


def _color_groups_from_blocks(blocks):
    groups = []
    by_key = {}
    for block in blocks:
        r_atoms = sorted(int(r) for r in block.get('r_atoms', ()))
        p_atoms = sorted(int(p) for p in block.get('p_atoms', ()))
        if len(r_atoms) <= 1 and len(p_atoms) <= 1:
            continue
        source = block.get('source') or 'sym_block'
        key = (tuple(r_atoms), tuple(p_atoms))
        if key in by_key:
            by_key[key]['sources'].add(source)
            continue
        group = {
            'r_atoms': r_atoms,
            'p_atoms': p_atoms,
            'sources': {source},
        }
        by_key[key] = group
        groups.append(group)
    for group in groups:
        group['sources'] = sorted(group['sources'])
    return groups


def _generator_orbit_map(g, generators):
    """Return exact atom orbits generated by pynauty permutations."""
    parent = {int(atom): int(atom) for atom in g.nodes()}

    def find(atom):
        atom = int(atom)
        while parent[atom] != atom:
            parent[atom] = parent[parent[atom]]
            atom = parent[atom]
        return atom

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for generator in generators:
        for atom, image in generator.items():
            union(atom, image)
    return {atom: find(atom) for atom in parent}


def complete_chosen_automorphism_groups(branch_symmetry, mapping, g_R, g_P,
                                        wbo_tol):
    """Complete chosen-candidate display blocks using exact stabilizers.

    Growth blocks are a compressed search trace and can omit an atom assigned
    earlier as a singleton.  For each mapped center, intersect the exact R and
    P neighbor orbits under automorphisms that fix that center.  The resulting
    sibling groups are therefore mutable at both endpoints and belong to the
    chosen mapping, rather than to discarded branches or orbit history.
    """
    mapping = {int(r): int(p) for r, p in dict(mapping).items()}
    result = dict(branch_symmetry or {})
    complete = []

    def stabilizer_orbits(g, center, cache):
        if center not in cache:
            generators = _nauty_atom_generators(
                g, wbo_tol=float(wbo_tol),
                atom_color_tags={int(center): ('fixed_center',)})
            cache[center] = _generator_orbit_map(g, generators)
        return cache[center]

    fragments = result.get('fragments') or [{
        'fragment_index': 0,
        'island_idx': 0,
        'fragment': sorted(mapping),
    }]
    for fragment_position, fragment in enumerate(fragments):
        fragment_R = {
            int(atom) for atom in fragment.get('fragment', ())
            if int(atom) in mapping and int(atom) in g_R
        }
        fragment_P = {mapping[atom] for atom in fragment_R}
        if not fragment_R or not fragment_P:
            continue
        # Final mutability belongs to this chosen fragment.  Edges crossing
        # its boundary are masked, just as they were for the candidate.
        sub_R = g_R.subgraph(fragment_R).copy()
        sub_P = g_P.subgraph(fragment_P).copy()
        r_cache = {}
        p_cache = {}
        fragment_group_start = len(complete)
        for center_R in sorted(fragment_R):
            center_P = mapping[center_R]
            if sub_R.degree(center_R) < 2 or sub_P.degree(center_P) < 2:
                continue
            r_orbits = stabilizer_orbits(sub_R, center_R, r_cache)
            p_orbits = stabilizer_orbits(sub_P, center_P, p_cache)
            groups = {}
            for atom_R in sub_R.neighbors(center_R):
                atom_P = mapping.get(int(atom_R))
                if atom_P is None or not sub_P.has_edge(center_P, atom_P):
                    continue
                key = (int(r_orbits[atom_R]), int(p_orbits[atom_P]))
                groups.setdefault(key, []).append((int(atom_R), int(atom_P)))
            for pairs in groups.values():
                if len(pairs) <= 1:
                    continue
                r_atoms = sorted(r for r, _ in pairs)
                p_atoms = sorted(p for _, p in pairs)
                complete.append({
                    'fragment_index': int(fragment.get(
                        'fragment_index', fragment_position)),
                    'island_idx': int(fragment.get(
                        'island_idx', fragment_position)),
                    'center_R': int(center_R),
                    'center_P': int(center_P),
                    'r_atoms': r_atoms,
                    'p_atoms': p_atoms,
                    'extendable': False,
                    'open': False,
                    'assignments': f"{len(r_atoms)}!",
                    'source': 'chosen_candidate_automorph',
                })

        # Some fragments have no degree-2 center (a symmetric diatomic or a
        # pair of equivalent disconnected atoms).  In that case the exact
        # whole-fragment orbit is the only meaningful fragment-level group.
        if len(complete) == fragment_group_start and len(fragment_R) > 1:
            r_orbits = _nauty_orbits(sub_R, wbo_tol=float(wbo_tol))
            p_orbits = _nauty_orbits(sub_P, wbo_tol=float(wbo_tol))
            groups = {}
            for atom_R in sorted(fragment_R):
                atom_P = mapping[atom_R]
                key = (int(r_orbits[atom_R]), int(p_orbits[atom_P]))
                groups.setdefault(key, []).append((atom_R, atom_P))
            for pairs in groups.values():
                if len(pairs) <= 1:
                    continue
                r_atoms = sorted(r for r, _ in pairs)
                p_atoms = sorted(p for _, p in pairs)
                complete.append({
                    'fragment_index': int(fragment.get(
                        'fragment_index', fragment_position)),
                    'island_idx': int(fragment.get(
                        'island_idx', fragment_position)),
                    'r_atoms': r_atoms,
                    'p_atoms': p_atoms,
                    'extendable': False,
                    'open': False,
                    'assignments': f"{len(r_atoms)}!",
                    'source': 'chosen_fragment_automorph',
                })

    unique_complete = {}
    for block in complete:
        key = (tuple(block['r_atoms']), tuple(block['p_atoms']))
        unique_complete.setdefault(key, block)
    complete = list(unique_complete.values())

    result['rule'] = 'chosen_analytical_branch_fragment_automorphisms'
    result.pop('witnesses', None)
    result.pop('matching_generators', None)
    result.pop('matching_blocks', None)
    result.pop('dedup_witness_count', None)
    result.pop('selected_witness_index', None)
    result['blocks'] = complete
    result['color_groups'] = _color_groups_from_blocks(complete)
    return result


def _freeze_analytical(value):
    if isinstance(value, dict):
        return tuple(sorted((str(key), _freeze_analytical(item))
                            for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_analytical(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return tuple(sorted((_freeze_analytical(item) for item in value),
                            key=repr))
    return value


def _analytical_branch(mapping, cuts, branch_symmetry=None, count=1):
    """One unique completed branch family retained below a mechanism."""
    return {
        'mapping': {int(r): int(p) for r, p in dict(mapping).items()},
        'cuts': _cut_record(cuts),
        'encounter_count': int(count),
        'hierarchy': copy.deepcopy(branch_symmetry or {
            'rule': 'branch_symmetry_blocks',
            'fragments': [],
            'blocks': [],
        }),
    }


def _analytical_branch_key(branch):
    # Cut/seed provenance is deliberately excluded.  Until a coupled
    # transporter is stored, only literal equality of the complete mapping
    # and hierarchy proves that two completed branch families are identical.
    return (
        _mapping_key(branch.get('mapping') or {}),
        _freeze_analytical(branch.get('hierarchy') or {}),
    )


def _merge_analytical_branch(branches, incoming, key_index=None):
    key = _analytical_branch_key(incoming)
    if key_index is not None:
        branch = key_index.get(key)
        candidates = () if branch is None else (branch,)
    else:
        candidates = (
            branch for branch in branches
            if _analytical_branch_key(branch) == key
        )
    for branch in candidates:
        branch['encounter_count'] = (
            int(branch.get('encounter_count', 1))
            + int(incoming.get('encounter_count', 1)))
        cuts = {
            tuple(map(int, cut))
            for cut in branch.get('cuts', ())
        } | {
            tuple(map(int, cut))
            for cut in incoming.get('cuts', ())
        }
        branch['cuts'] = _cut_record(cuts)
        return False
    branches.append(incoming)
    if key_index is not None:
        key_index[key] = incoming
    return True


def _refresh_entry_branch_symmetry(entry):
    branches = list(entry.get('branches') or ())
    representative_key = _mapping_key(entry.get('mapping') or {})
    selected = next((
        branch for branch in branches
        if _mapping_key(branch.get('mapping') or {}) == representative_key
    ), branches[0] if branches else None)
    hierarchy = copy.deepcopy(
        (selected or {}).get('hierarchy') or {
            'fragments': [], 'blocks': [],
        })
    hierarchy['rule'] = 'selected_analytical_branch'
    hierarchy['analytical_branch_count'] = len(branches)
    hierarchy['selected_branch_index'] = (
        branches.index(selected) if selected in branches else None)
    hierarchy['blocks'] = [
        block for block in hierarchy.get('blocks') or ()
        if (block.get('source') or 'sym_block') in {
            'sym_block', 'exact_automorph_group',
        }
    ]
    hierarchy['color_groups'] = _color_groups_from_blocks(hierarchy['blocks'])
    # Explicitly prevent legacy consumers from treating branch
    # representatives as a sampled permutation group.
    hierarchy.pop('witnesses', None)
    hierarchy.pop('matching_generators', None)
    hierarchy.pop('matching_blocks', None)
    hierarchy.pop('dedup_witness_count', None)
    entry['branch_symmetry'] = hierarchy


def _pool_add(pool, sig, mapping, cuts, branch_symmetry=None):
    cuts = frozenset(cuts)
    no_cut = not cuts
    branch = _analytical_branch(mapping, cuts, branch_symmetry)
    entry = pool.get(sig)
    if entry is None:
        pool[sig] = {
            'mapping': dict(mapping),
            'cuts': cuts,
            'has_no_cut': bool(no_cut),
            'dedup_count': 1,
            'branches': [branch],
            '_branch_key_index': {_analytical_branch_key(branch): branch},
        }
        _refresh_entry_branch_symmetry(pool[sig])
    else:
        representative_changed = bool(
            no_cut and not entry.get('has_no_cut', False))
        if no_cut and not entry.get('has_no_cut', False):
            entry['mapping'] = dict(mapping)
        entry['has_no_cut'] = bool(entry.get('has_no_cut', False) or no_cut)
        entry['cuts'] = entry['cuts'] | cuts
        entry['dedup_count'] = entry.get('dedup_count', 1) + 1
        branches = entry.setdefault('branches', [])
        key_index = entry.get('_branch_key_index')
        if key_index is None:
            key_index = {
                _analytical_branch_key(existing): existing
                for existing in branches
            }
            entry['_branch_key_index'] = key_index
        added = _merge_analytical_branch(
            branches, branch, key_index=key_index)
        if added or representative_changed:
            _refresh_entry_branch_symmetry(entry)


def _public_pool(pool):
    """Drop transient exact-key indexes before returning the public pool."""
    for entry in pool.values():
        entry.pop('_branch_key_index', None)
    return pool


def _anchor_mapping_ok(mapping, anchor_map):
    return all(
        int(mapping.get(int(r), -1)) == int(p)
        for r, p in dict(anchor_map or {}).items()
    )


def _run_find_islands_limited(g_R, g_P, order, core_R, cfg, *,
                              p_orbits=None, r_orbits=None,
                              profile=None):
    stop_on_core = bool(core_R)
    return find_islands(
        g_R, g_P, list(order),
        iso_tol=float(cfg['iso_tol']),
        max_branches=int(cfg['max_branches']),
        dwbo_threshold=float(cfg['dwbo_threshold']),
        metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'),
        symmetry_wbo_tol=float(cfg['symmetry_wbo_tol']),
        core_R=core_R,
        stop_when_core_mapped=stop_on_core,
        p_orbits=p_orbits,
        r_orbits=r_orbits,
        profile=profile,
        anchor_map=cfg.get('anchor_map'),
    )


def _score_branch_mapping(mapping, g_R, g_P, wboR, wboT,
                          g_R_full, p_orbits, r_orbits, core_R, cfg,
                          elR=None, elT=None,
                          event_canonicalizer=None,
                          return_repair_stats=False):
    anchor_map = cfg.get('anchor_map') or {}
    if core_R:
        if not all(r in mapping for r in core_R):
            return None
        if not _anchor_mapping_ok(mapping, anchor_map):
            return None
        scored = (_core_mapping_key(mapping, core_R), mapping)
        return (*scored, None) if return_repair_stats else scored

    if len(mapping) < int(cfg['n_atoms']) - 2:
        return None
    if not _anchor_mapping_ok(mapping, anchor_map):
        return None
    repair_stats = None
    if cfg['symmetry_repair']:
        base_mapping = dict(mapping)
        if return_repair_stats:
            mapping, repair_stats = symmetry_repair_mapping(
                mapping, wboR, wboT, g_R_full, g_P, p_orbits,
                dwbo_threshold=float(cfg['dwbo_threshold']),
                metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'),
                min_changes=int(cfg['symmetry_repair_min_changes']),
                max_evals=int(cfg['symmetry_repair_max_evals']),
                return_stats=True,
            )
        else:
            mapping = symmetry_repair_mapping(
                mapping, wboR, wboT, g_R_full, g_P, p_orbits,
                dwbo_threshold=float(cfg['dwbo_threshold']),
                metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'),
                min_changes=int(cfg['symmetry_repair_min_changes']),
                max_evals=int(cfg['symmetry_repair_max_evals']),
            )
        if not _anchor_mapping_ok(mapping, anchor_map):
            mapping = base_mapping
            if repair_stats is not None:
                repair_stats = dict(repair_stats)
                repair_stats['anchor_reverted'] = True
    sig = _mechanism_signature(
        mapping, wboR, wboT, r_orbits, p_orbits,
        dwbo_threshold=float(cfg['dwbo_threshold']),
        elements_R=elR,
        elements_P=elT,
        metal_dwbo_threshold=cfg.get('metal_dwbo_threshold'),
        g_R_full=g_R_full,
        symmetry_wbo_tol=float(cfg['symmetry_wbo_tol']),
        event_canonicalizer=event_canonicalizer)
    scored = (sig, mapping)
    return (*scored, repair_stats) if return_repair_stats else scored


def _cut_json(cut):
    return [list(map(int, pair)) for pair in cut]


def _repair_trace_stats(stats):
    if not stats:
        return None
    groups = stats.get('groups') or []
    return {
        'enabled': bool(stats.get('enabled', False)),
        'base_changes': stats.get('base_changes'),
        'best_changes': stats.get('best_changes'),
        'repaired': bool(stats.get('repaired', False)),
        'evaluated': int(stats.get('evaluated', 0) or 0),
        'capped': bool(stats.get('capped', False)),
        'n_groups': len(groups),
        'group_sizes': [int(group.get('size', 0) or 0)
                        for group in groups],
    }


def _growth_trace_summary(profile):
    if not profile:
        return {
            'calls': 0,
            'elapsed_sec': 0.0,
            'extend_elapsed_sec': 0.0,
            'heap_pops': 0,
            'extend_calls': 0,
            'commits': 0,
            'deferred': 0,
            'max_cands_before': 0,
            'max_cands_after': 0,
            'max_heap_len': 0,
            'max_fragment_size': 0,
            'slowest': [],
        }
    slowest = sorted(
        profile,
        key=lambda item: float(item.get('elapsed_sec', 0.0) or 0.0),
        reverse=True,
    )[:5]
    return {
        'calls': len(profile),
        'elapsed_sec': sum(float(item.get('elapsed_sec', 0.0) or 0.0)
                           for item in profile),
        'extend_elapsed_sec': sum(float(item.get('extend_elapsed_sec', 0.0) or 0.0)
                                  for item in profile),
        'heap_pops': sum(int(item.get('heap_pops', 0) or 0)
                         for item in profile),
        'extend_calls': sum(int(item.get('extend_calls', 0) or 0)
                            for item in profile),
        'commits': sum(int(item.get('commits', 0) or 0)
                       for item in profile),
        'deferred': sum(int(item.get('deferred', 0) or 0)
                        for item in profile),
        'max_cands_before': max(int(item.get('max_cands_before', 0) or 0)
                                for item in profile),
        'max_cands_after': max(int(item.get('max_cands_after', 0) or 0)
                               for item in profile),
        'max_heap_len': max(int(item.get('max_heap_len', 0) or 0)
                            for item in profile),
        'max_fragment_size': max(int(item.get('max_fragment_size', 0) or 0)
                                 for item in profile),
        'slowest': [
            {
                'seed': int(item.get('seed')),
                'pass': int(item.get('pass', 0) or 0),
                'branch_index': int(item.get('branch_index', 0) or 0),
                'mapped_before': int(item.get('mapped_before', 0) or 0),
                'result': item.get('result'),
                'elapsed_sec': float(item.get('elapsed_sec', 0.0) or 0.0),
                'extend_elapsed_sec': float(
                    item.get('extend_elapsed_sec', 0.0) or 0.0),
                'heap_pops': int(item.get('heap_pops', 0) or 0),
                'extend_calls': int(item.get('extend_calls', 0) or 0),
                'commits': int(item.get('commits', 0) or 0),
                'deferred': int(item.get('deferred', 0) or 0),
                'max_cands_before': int(item.get('max_cands_before', 0) or 0),
                'max_cands_after': int(item.get('max_cands_after', 0) or 0),
                'slowest_extend': item.get('slowest_extend'),
            }
            for item in slowest
        ],
    }


def _emit_trace(trace_path, events):
    if not trace_path or not events:
        return
    path = Path(trace_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a', encoding='utf-8') as handle:
        for event in events:
            handle.write(json.dumps(event, sort_keys=True) + '\n')


def _run_cut_work(elR, wboR, elT, wboT, cfg, cut, orders, core_R,
                  g_P, g_R_full, p_orbits, r_orbits_full,
                  *, return_trace=False):
    cut = tuple(tuple(int(v) for v in pair) for pair in cut)
    events = []
    out = []
    cut_t0 = time.perf_counter()

    graph_t0 = time.perf_counter()
    graph_floor = float(cfg['graph_floor'])
    event_canonicalizer = _MechanismEventCanonicalizer(
        g_R_full, wbo_tol=float(cfg['symmetry_wbo_tol']))
    g_R = build_graph(elR, wboR, bond_cut=graph_floor)
    for i, j in cut:
        if g_R.has_edge(i, j):
            g_R.remove_edge(i, j)
    r_orbits_cut = _nauty_orbits(
        g_R, wbo_tol=float(cfg['symmetry_wbo_tol']))
    graph_elapsed = time.perf_counter() - graph_t0
    if orders is None:
        orders = _generate_seed_orders(
            g_R, n_trials=int(cfg['n_seeds']))
    else:
        orders = list(orders)
    if return_trace:
        events.append({
            'event': 'cut_start',
            'cut': _cut_json(cut),
            'n_orders': len(orders),
            'max_branches': int(cfg['max_branches']),
            'symmetry_repair': bool(cfg['symmetry_repair']),
            'symmetry_repair_max_evals': int(
                cfg['symmetry_repair_max_evals']),
        })

    cut_status = 'completed'
    total_search_elapsed = 0.0
    total_score_elapsed = 0.0
    total_expand_elapsed = 0.0
    total_branches = 0
    total_accepted = 0
    total_repair_evals = 0
    total_repair_capped = 0
    try:
        for order_index, order in enumerate(orders):
            seed_t0 = time.perf_counter()
            seed_growth_profile = [] if return_trace else None
            try:
                branches = _run_find_islands_limited(
                    g_R, g_P, order, core_R, cfg,
                    p_orbits=p_orbits,
                    r_orbits=r_orbits_cut,
                    profile=seed_growth_profile)
            except BranchLimitExceeded as exc:
                cut_status = 'branch_cap'
                if return_trace:
                    for item in seed_growth_profile:
                        growth_event = dict(item)
                        growth_event['event'] = 'growth_call'
                        growth_event['cut'] = _cut_json(cut)
                        growth_event['seed_index'] = int(order_index)
                        events.append(growth_event)
                    events.append({
                        'event': 'seed_branch_cap',
                        'cut': _cut_json(cut),
                        'seed_index': int(order_index),
                        'elapsed_sec': time.perf_counter() - seed_t0,
                        'branch_count': exc.branch_count,
                        'max_branches': exc.max_branches,
                        'stage': exc.stage,
                        'growth': _growth_trace_summary(seed_growth_profile),
                    })
                raise
            except Exception as exc:
                if return_trace:
                    for item in seed_growth_profile:
                        growth_event = dict(item)
                        growth_event['event'] = 'growth_call'
                        growth_event['cut'] = _cut_json(cut)
                        growth_event['seed_index'] = int(order_index)
                        events.append(growth_event)
                    events.append({
                        'event': 'seed_error',
                        'cut': _cut_json(cut),
                        'seed_index': int(order_index),
                        'elapsed_sec': time.perf_counter() - seed_t0,
                        'error_type': type(exc).__name__,
                        'error': str(exc),
                        'growth': _growth_trace_summary(seed_growth_profile),
                    })
                raise

            search_elapsed = time.perf_counter() - seed_t0
            growth_summary = _growth_trace_summary(seed_growth_profile)
            total_search_elapsed += search_elapsed
            seed_score_elapsed = 0.0
            seed_expand_elapsed = 0.0
            seed_accepted = 0
            seed_repair_evals = 0
            seed_repair_capped = 0
            max_mapped = 0
            for branch_index, branch in enumerate(branches):
                branch_t0 = time.perf_counter()
                expand_elapsed = 0.0
                score_t0 = time.perf_counter()
                if core_R:
                    scored_items = []
                    try:
                        core_variants = _core_mapping_variants(
                            branch, core_R, int(cfg['max_branches']),
                            g_P=g_P, p_orbits=p_orbits)
                    except BranchLimitExceeded:
                        # This branch's explicit core expansion is the capped
                        # sub-result; sibling branches and seed orders remain.
                        core_variants = ()
                    for core_map in core_variants:
                        scored_items.append((
                            _core_mapping_key(core_map, core_R),
                            core_map,
                            None,
                        ))
                    mapping_for_stats = dict(branch.mapping)
                else:
                    expand_t0 = time.perf_counter()
                    mapping = expand_mapping(dict(branch.mapping), g_R, g_P)
                    expand_elapsed = time.perf_counter() - expand_t0
                    scored = _score_branch_mapping(
                        mapping, g_R, g_P, wboR, wboT, g_R_full,
                        p_orbits, r_orbits_full, core_R, cfg, elR, elT,
                        event_canonicalizer=event_canonicalizer,
                        return_repair_stats=return_trace)
                    if scored is None:
                        scored_items = []
                    elif return_trace:
                        sig, repaired_mapping, repair_stats = scored
                        scored_items = [(sig, repaired_mapping, repair_stats)]
                    else:
                        sig, repaired_mapping = scored
                        scored_items = [(sig, repaired_mapping, None)]
                    mapping_for_stats = mapping
                score_elapsed = time.perf_counter() - score_t0
                branch_elapsed = time.perf_counter() - branch_t0
                seed_expand_elapsed += expand_elapsed
                seed_score_elapsed += score_elapsed
                max_mapped = max(max_mapped, len(mapping_for_stats))
                accepted = bool(scored_items)
                for sig, accepted_mapping, _repair_stats in scored_items:
                    for derived_mapping, hierarchy in (
                            _branch_analytical_derivations(
                                branch, accepted_mapping)):
                        out.append((
                            sig,
                            tuple(sorted(derived_mapping.items())),
                            cut,
                            hierarchy,
                        ))
                        seed_accepted += 1
                repair_stats = scored_items[0][2] if scored_items else None
                repair_summary = _repair_trace_stats(repair_stats)
                if repair_summary:
                    seed_repair_evals += repair_summary['evaluated']
                    seed_repair_capped += int(repair_summary['capped'])
                if return_trace:
                    events.append({
                        'event': 'branch',
                        'cut': _cut_json(cut),
                        'seed_index': int(order_index),
                        'branch_index': int(branch_index),
                        'elapsed_sec': branch_elapsed,
                        'expand_elapsed_sec': expand_elapsed,
                        'score_elapsed_sec': score_elapsed,
                        'mapped_atoms': len(mapping_for_stats),
                        'accepted': bool(accepted),
                        'accepted_core_variants': len(scored_items),
                        'repair': repair_summary,
                    })
            n_branches = len(branches)
            total_branches += n_branches
            total_accepted += seed_accepted
            total_expand_elapsed += seed_expand_elapsed
            total_score_elapsed += seed_score_elapsed
            total_repair_evals += seed_repair_evals
            total_repair_capped += seed_repair_capped
            if return_trace:
                for item in seed_growth_profile:
                    growth_event = dict(item)
                    growth_event['event'] = 'growth_call'
                    growth_event['cut'] = _cut_json(cut)
                    growth_event['seed_index'] = int(order_index)
                    events.append(growth_event)
                events.append({
                    'event': 'seed_end',
                    'cut': _cut_json(cut),
                    'seed_index': int(order_index),
                    'search_elapsed_sec': search_elapsed,
                    'expand_elapsed_sec': seed_expand_elapsed,
                    'score_elapsed_sec': seed_score_elapsed,
                    'branches': n_branches,
                    'accepted': seed_accepted,
                    'max_mapped_atoms': max_mapped,
                    'repair_evals': seed_repair_evals,
                    'repair_capped_count': seed_repair_capped,
                    'growth': growth_summary,
                })
    except BranchLimitExceeded:
        # No broad cut-level fallback: a cap must be handled at the concrete
        # subtree that owns it.  Any uncaught cap is an implementation error.
        raise

    elapsed = time.perf_counter() - cut_t0
    if return_trace:
        events.append({
            'event': 'cut_end',
            'cut': _cut_json(cut),
            'status': cut_status,
            'elapsed_sec': elapsed,
            'graph_elapsed_sec': graph_elapsed,
            'search_elapsed_sec': total_search_elapsed,
            'expand_elapsed_sec': total_expand_elapsed,
            'score_elapsed_sec': total_score_elapsed,
            'branches': total_branches,
            'accepted': total_accepted,
            'repair_evals': total_repair_evals,
            'repair_capped_count': total_repair_capped,
            'hits': len(out),
            'avg_branch_elapsed_sec': (
                (total_expand_elapsed + total_score_elapsed) / total_branches
                if total_branches else 0.0
            ),
        })
    return out, events


_WORKER = {}


def _cs_winit(elR, wboR, elT, wboT, cfg):
    graph_floor = float(cfg['graph_floor'])
    _WORKER['elR'] = elR
    _WORKER['wboR'] = wboR
    _WORKER['elT'] = elT
    _WORKER['wboT'] = wboT
    _WORKER['cfg'] = dict(cfg)
    _WORKER['g_P'] = build_graph(elT, wboT, bond_cut=graph_floor)
    _WORKER['g_R_full'] = build_graph(elR, wboR, bond_cut=graph_floor)
    symmetry_wbo_tol = float(cfg['symmetry_wbo_tol'])
    _WORKER['p_orbits'] = _nauty_orbits(
        _WORKER['g_P'], wbo_tol=symmetry_wbo_tol)
    _WORKER['r_orbits'] = _nauty_orbits(
        _WORKER['g_R_full'], wbo_tol=symmetry_wbo_tol)


def _cs_wrun(args):
    cut, orders, core_R, trace_enabled = args
    cfg = _WORKER['cfg']
    out, events = _run_cut_work(
        _WORKER['elR'], _WORKER['wboR'],
        _WORKER['elT'], _WORKER['wboT'],
        cfg, cut, orders, core_R,
        _WORKER['g_P'], _WORKER['g_R_full'],
        _WORKER['p_orbits'], _WORKER['r_orbits'],
        return_trace=trace_enabled)
    return {'results': out, 'events': events}


def cut_sweep_items(wboR, cut_floor=0.2):
    """Return the independent no-cut plus one-edge cut work items."""
    return [()] + [((int(i), int(j)),) for i, j in _strong_edges(
        wboR, float(cut_floor))]


def _cut_sweep_cfg(*, cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
                   dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
                   symmetry_wbo_tol=0.2,
                   n_seeds=3, max_branches=100,
                   chunksize=1,
                   symmetry_repair=True,
                   symmetry_repair_min_changes=1,
                   symmetry_repair_max_evals=20000,
                   n_atoms=0,
                   anchor_map=None):
    anchor_map = {
        int(r): int(p)
        for r, p in dict(anchor_map or {}).items()
    }
    return {
        'cut_floor': float(cut_floor),
        'graph_floor': float(graph_floor),
        'iso_tol': float(iso_tol),
        'dwbo_threshold': float(dwbo_threshold),
        'metal_dwbo_threshold': (
            None if metal_dwbo_threshold is None
            else float(metal_dwbo_threshold)
        ),
        # One tolerance governs both edge verification and pynauty colors.
        'symmetry_wbo_tol': float(iso_tol),
        'n_seeds': int(n_seeds),
        'max_branches': int(max_branches),
        'chunksize': int(chunksize),
        'symmetry_repair': bool(symmetry_repair),
        'symmetry_repair_min_changes': int(symmetry_repair_min_changes),
        'symmetry_repair_max_evals': int(symmetry_repair_max_evals),
        'n_atoms': int(n_atoms),
        'anchor_map': anchor_map,
    }


def _cut_sweep_chunk_serial(elR, wboR, elT, wboT, cfg, core_R, cuts,
                            trace_path=None):
    graph_floor = float(cfg['graph_floor'])
    g_P = build_graph(elT, wboT, bond_cut=graph_floor)
    g_R_full = build_graph(elR, wboR, bond_cut=graph_floor)
    symmetry_wbo_tol = float(cfg['symmetry_wbo_tol'])
    p_orbits = _nauty_orbits(g_P, wbo_tol=symmetry_wbo_tol)
    r_orbits = _nauty_orbits(g_R_full, wbo_tol=symmetry_wbo_tol)
    pool = {}

    for cut in cuts:
        cut = tuple(tuple(pair) for pair in cut)
        results, events = _run_cut_work(
            elR, wboR, elT, wboT, cfg, cut, None, core_R,
            g_P, g_R_full, p_orbits, r_orbits,
            return_trace=bool(trace_path))
        _emit_trace(trace_path, events)
        for result in results:
            sig, mapping_items, _cut = result[:3]
            branch_symmetry = result[3] if len(result) > 3 else None
            mapping = dict(mapping_items)
            _pool_add(pool, sig, mapping, _cut, branch_symmetry)
    return pool


def _cut_sweep_serial(elR, wboR, elT, wboT, cfg, core_R):
    return _public_pool(_cut_sweep_chunk_serial(
        elR, wboR, elT, wboT, cfg, core_R,
        cut_sweep_items(wboR, cfg['cut_floor'])))


def _cut_sweep_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R):
    cuts = cut_sweep_items(wboR, cfg['cut_floor'])
    work = [(cut, None, core_R, False) for cut in cuts]
    pool = {}
    with mp.Pool(n_workers, initializer=_cs_winit,
                 initargs=(elR, wboR, elT, wboT, cfg)) as proc_pool:
        for payload in proc_pool.imap_unordered(
                _cs_wrun, work, chunksize=max(1, int(cfg['chunksize']))):
            for result in payload['results']:
                sig, mapping_items, cut = result[:3]
                branch_symmetry = result[3] if len(result) > 3 else None
                _pool_add(pool, sig, dict(mapping_items), cut, branch_symmetry)
    return _public_pool(pool)


def _cut_sweep_chunk_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R,
                              cuts, trace_path=None):
    work = [(cut, None, core_R, bool(trace_path)) for cut in cuts]
    pool = {}
    with mp.Pool(n_workers, initializer=_cs_winit,
                 initargs=(elR, wboR, elT, wboT, cfg)) as proc_pool:
        for payload in proc_pool.imap_unordered(
                _cs_wrun, work, chunksize=max(1, int(cfg['chunksize']))):
            _emit_trace(trace_path, payload.get('events', []))
            for result in payload['results']:
                sig, mapping_items, cut = result[:3]
                branch_symmetry = result[3] if len(result) > 3 else None
                _pool_add(pool, sig, dict(mapping_items), cut, branch_symmetry)
    return _public_pool(pool)


def run_cut_sweep_chunk(elR, wboR, elT, wboT, cuts, *,
                        core_R=None,
                        n_workers=None,
                        trace_path=None,
                        cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
                        dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
                        symmetry_wbo_tol=0.2,
                        n_seeds=3, max_branches=100,
                        chunksize=1,
                        symmetry_repair=True,
                        symmetry_repair_min_changes=1,
                        symmetry_repair_max_evals=20000,
                        anchor_map=None):
    """Run a chunk of independent cut-sweep work items.

    This is the Slurm-array friendly primitive.  The caller chooses which cut
    work items belong to the chunk; the merge step combines the returned pools.
    """
    cfg = _cut_sweep_cfg(
        cut_floor=cut_floor,
        graph_floor=graph_floor,
        iso_tol=iso_tol,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
        symmetry_wbo_tol=symmetry_wbo_tol,
        n_seeds=n_seeds,
        max_branches=max_branches,
        chunksize=chunksize,
        symmetry_repair=symmetry_repair,
        symmetry_repair_min_changes=symmetry_repair_min_changes,
        symmetry_repair_max_evals=symmetry_repair_max_evals,
        n_atoms=len(elR),
        anchor_map=anchor_map,
    )
    core_R = tuple(sorted(set(core_R or ())))
    normalized_cuts = [
        tuple(tuple(int(v) for v in pair) for pair in cut)
        for cut in cuts
    ]
    if n_workers and int(n_workers) > 1 and len(normalized_cuts) > 1:
        return _cut_sweep_chunk_parallel(
            elR, wboR, elT, wboT, cfg,
            min(int(n_workers), len(normalized_cuts)),
            core_R, normalized_cuts, trace_path=trace_path)
    return _public_pool(_cut_sweep_chunk_serial(
        elR, wboR, elT, wboT, cfg, core_R, normalized_cuts,
        trace_path=trace_path))


def run_no_cut_core_branch_records(elS, wboS, elT, wboT, core_S, *,
                                   graph_floor=0.2, iso_tol=1.0,
                                   dwbo_threshold=0.5,
                                   metal_dwbo_threshold=0.3,
                                   symmetry_wbo_tol=0.2,
                                   n_seeds=3, max_branches=20000):
    """Return no-cut endpoint->target branch states for a mechanism core.

    Unlike ``run_cut_sweep_chunk(..., core_R=...)``, this keeps branch-level
    compressed symmetry records.  Callers can merge R-derived and P-derived
    endpoint support first, then expand only the automorphism variants inside
    each merged branch state for TS/IG scoring.
    """
    core_S = tuple(sorted(set(int(r) for r in core_S or ())))
    if not core_S:
        return []

    g_S = build_graph(elS, wboS, bond_cut=float(graph_floor))
    g_T = build_graph(elT, wboT, bond_cut=float(graph_floor))
    symmetry_wbo_tol = float(iso_tol)
    p_orbits = _nauty_orbits(g_T, wbo_tol=symmetry_wbo_tol)
    r_orbits = _nauty_orbits(g_S, wbo_tol=symmetry_wbo_tol)
    orders = _generate_seed_orders(g_S, n_trials=int(n_seeds))
    records = {}

    try:
        for order in orders:
            branches = find_islands(
                g_S, g_T, list(order),
                iso_tol=float(iso_tol),
                max_branches=int(max_branches),
                dwbo_threshold=float(dwbo_threshold),
                metal_dwbo_threshold=metal_dwbo_threshold,
                symmetry_wbo_tol=float(symmetry_wbo_tol),
                core_R=core_S,
                stop_when_core_mapped=True,
                p_orbits=p_orbits,
                r_orbits=r_orbits,
            )
            for branch in branches:
                record = _core_branch_record(
                    branch, core_S, p_orbits=p_orbits)
                if record is None:
                    continue
                key = _core_branch_record_key(record, core_S)
                existing = records.get(key)
                if existing is None:
                    records[key] = record
                else:
                    existing['dedup_count'] = (
                        int(existing.get('dedup_count', 1))
                        + int(record.get('dedup_count', 1))
                    )
    except BranchLimitExceeded:
        return []

    return list(records.values())


def merge_cut_sweep_pools(pools):
    """Merge partial cut-sweep pools produced by chunk tasks."""
    def branches_from(info, cuts):
        if not info.get('branches'):
            raise ValueError(
                "cut-sweep chunk lacks analytical AAM branches; rerun it "
                "with the current branch schema")
        return [copy.deepcopy(branch) for branch in info['branches']]

    merged = {}
    for pool in pools:
        for sig, info in dict(pool or {}).items():
            cuts = frozenset(info.get('cuts', ()))
            no_cut = bool(info.get('has_no_cut', False))
            entry = merged.get(sig)
            if entry is None:
                branches = []
                for branch in branches_from(info, cuts):
                    _merge_analytical_branch(branches, branch)
                merged[sig] = {
                    'mapping': dict(info['mapping']),
                    'cuts': cuts,
                    'has_no_cut': no_cut,
                    'dedup_count': int(info.get('dedup_count', 1)),
                    'branches': branches,
                }
                _refresh_entry_branch_symmetry(merged[sig])
                continue
            if no_cut and not entry.get('has_no_cut', False):
                entry['mapping'] = dict(info['mapping'])
            entry['cuts'] = frozenset(entry.get('cuts', ())) | cuts
            entry['has_no_cut'] = bool(
                entry.get('has_no_cut', False)
                or no_cut
            )
            entry['dedup_count'] = (
                int(entry.get('dedup_count', 1))
                + int(info.get('dedup_count', 1))
            )
            for branch in branches_from(info, cuts):
                _merge_analytical_branch(
                    entry.setdefault('branches', []), branch)
            _refresh_entry_branch_symmetry(entry)
    return merged


def cut_sweep(elR, wboR, elT, wboT, *,
              n_workers=None, core_R=None,
              cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
              dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
              symmetry_wbo_tol=0.2,
              n_seeds=3, max_branches=100,
              chunksize=1,
              symmetry_repair=True,
              symmetry_repair_min_changes=1,
              symmetry_repair_max_evals=20000,
              anchor_map=None):
    """Enumerate mechanism classes via no-cut plus one-edge R cuts.

    The returned pool maps a symmetry-canonical signature to:

    - `mapping`: representative symmetry-aware witness
    - `cuts`: set of R-edge cuts that led to that signature
    - `dedup_count`: number of witnesses collapsed into the signature

    With `core_R=None`, signatures are R-P mechanism signatures.  With
    `core_R` supplied, signatures are exact core mappings; this is useful for
    mechanism-local TS/IG scoring, but R-P mechanism discovery is the primary
    use.
    """
    cfg = _cut_sweep_cfg(
        cut_floor=cut_floor,
        graph_floor=graph_floor,
        iso_tol=iso_tol,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold,
        symmetry_wbo_tol=symmetry_wbo_tol,
        n_seeds=n_seeds,
        max_branches=max_branches,
        chunksize=chunksize,
        symmetry_repair=symmetry_repair,
        symmetry_repair_min_changes=symmetry_repair_min_changes,
        symmetry_repair_max_evals=symmetry_repair_max_evals,
        n_atoms=len(elR),
        anchor_map=anchor_map,
    )
    core_R = tuple(sorted(set(core_R or ())))
    if not n_workers or n_workers <= 1:
        return _cut_sweep_serial(elR, wboR, elT, wboT, cfg, core_R)
    return _cut_sweep_parallel(elR, wboR, elT, wboT, cfg,
                               int(n_workers), core_R)


def select_min_mechanisms(pool):
    """Keep only signatures with the fewest broken+formed orbit events."""
    if not pool:
        return {}
    best = min(len(sig[0]) + len(sig[1]) for sig in pool)
    return {
        sig: info for sig, info in pool.items()
        if len(sig[0]) + len(sig[1]) == best
    }

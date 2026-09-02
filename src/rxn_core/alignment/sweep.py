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
import pickle
import time
from collections import deque
from pathlib import Path

from ..frag import build_graph, classify_bonds, expand_mapping
from ..matcher import (
    _SymBlock,
    _SymCand,
    _atom_tuple_orbit,
    _nauty_atom_generators,
    _nauty_orbits,
    _sym_block_assignment_expr,
)
from ..matcher.canonical import _CandidateAutomorphismCanonicalizer
from ..matcher.orbits import _nauty_colored_wbo_graph
from .branch import (
    BranchLimitExceeded,
    _generate_seed_order,
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


def _perm_compose(first, second):
    """Return ``second o first`` for dense image permutations."""
    return tuple(second[first[atom]] for atom in range(len(first)))


def _perm_inverse(permutation):
    inverse = [0] * len(permutation)
    for atom, image in enumerate(permutation):
        inverse[image] = atom
    return tuple(inverse)


def _point_stabilizer_generators(generators, point, degree):
    """Exact Schreier generators for the subgroup fixing ``point``."""
    identity = tuple(range(int(degree)))
    generators = tuple(dict.fromkeys(
        tuple(map(int, generator)) for generator in generators
        if tuple(map(int, generator)) != identity
    ))
    transversals = {int(point): identity}
    queue = deque([int(point)])
    while queue:
        current = queue.popleft()
        transversal = transversals[current]
        for generator in generators:
            image = generator[current]
            if image in transversals:
                continue
            transversals[image] = _perm_compose(transversal, generator)
            queue.append(image)
    stabilizers = []
    seen = set()
    for current, transversal in transversals.items():
        for generator in generators:
            image = generator[current]
            stabilizer = _perm_compose(
                _perm_compose(transversal, generator),
                _perm_inverse(transversals[image]))
            if stabilizer != identity and stabilizer not in seen:
                seen.add(stabilizer)
                stabilizers.append(stabilizer)
    return tuple(stabilizers)


def _stored_generator_orbit_map(atoms, generators):
    parent = {int(atom): int(atom) for atom in atoms}

    def find(atom):
        while parent[atom] != atom:
            parent[atom] = parent[parent[atom]]
            atom = parent[atom]
        return atom

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for generator in generators:
        for atom in tuple(parent):
            image = int(generator[atom])
            if image in parent:
                union(atom, image)
    return {atom: find(atom) for atom in parent}


def complete_chosen_automorphism_groups(branch_symmetry, mapping, g_R, g_P,
                                        wbo_tol, *,
                                        exact_target_generators=None,
                                        exact_branch_mappings=()):
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
    if exact_target_generators is not None:
        degree = len(mapping)
        exact_target_generators = tuple(dict.fromkeys(
            tuple(map(int, generator))
            for generator in exact_target_generators))
        if any(len(generator) != degree
               or set(generator) != set(range(degree))
               for generator in exact_target_generators):
            raise ValueError("stored AAM target generator is not bijective")
    exact_branch_mappings = tuple(
        {int(r): int(p) for r, p in dict(candidate).items()}
        for candidate in exact_branch_mappings or ())

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
        if exact_target_generators is not None:
            fragment_group_start = len(complete)
            for center_R in sorted(fragment_R):
                center_P = mapping[center_R]
                neighbors_R = [
                    int(atom) for atom in sub_R.neighbors(center_R)
                    if mapping.get(int(atom)) in fragment_P
                    and sub_P.has_edge(center_P, mapping[int(atom)])
                ]
                if len(neighbors_R) < 2:
                    continue
                stabilizers = _point_stabilizer_generators(
                    exact_target_generators, center_P, degree)
                orbit_map = _stored_generator_orbit_map(
                    fragment_P, stabilizers)
                groups = {}
                for atom_R in neighbors_R:
                    atom_P = mapping[atom_R]
                    groups.setdefault(orbit_map[atom_P], []).append(
                        (atom_R, atom_P))
                for pairs in groups.values():
                    if len(pairs) <= 1:
                        continue
                    complete.append({
                        'fragment_index': int(fragment.get(
                            'fragment_index', fragment_position)),
                        'island_idx': int(fragment.get(
                            'island_idx', fragment_position)),
                        'center_R': int(center_R),
                        'center_P': int(center_P),
                        'r_atoms': sorted(r for r, _ in pairs),
                        'p_atoms': sorted(p for _, p in pairs),
                        'extendable': False,
                        'open': False,
                        'assignments': f"{len(pairs)}!",
                        'source': 'stored_AAM_branch_mapping_group',
                    })
            if (len(complete) == fragment_group_start
                    and len(fragment_R) > 1):
                orbit_map = _stored_generator_orbit_map(
                    fragment_P, exact_target_generators)
                groups = {}
                for atom_R in sorted(fragment_R):
                    atom_P = mapping[atom_R]
                    groups.setdefault(orbit_map[atom_P], []).append(
                        (atom_R, atom_P))
                for pairs in groups.values():
                    if len(pairs) <= 1:
                        continue
                    complete.append({
                        'fragment_index': int(fragment.get(
                            'fragment_index', fragment_position)),
                        'island_idx': int(fragment.get(
                            'island_idx', fragment_position)),
                        'r_atoms': sorted(r for r, _ in pairs),
                        'p_atoms': sorted(p for _, p in pairs),
                        'extendable': False,
                        'open': False,
                        'assignments': f"{len(pairs)}!",
                        'source': 'stored_AAM_branch_mapping_group',
                    })
            # Separate exact branch families can carry correlated assignment
            # changes absent from every individual branch stabilizer.  Record
            # only permutations directly witnessed between deduplicated AAM
            # branch representatives with the same center and neighbor set.
            for center_R in sorted(fragment_R):
                center_P = mapping[center_R]
                neighbors_R = tuple(sorted(sub_R.neighbors(center_R)))
                selected_images = tuple(mapping[r] for r in neighbors_R)
                selected_by_image = {
                    mapping[r]: int(r) for r in neighbors_R}
                local_parent = {int(r): int(r) for r in neighbors_R}

                def local_find(atom):
                    while local_parent[atom] != atom:
                        local_parent[atom] = local_parent[local_parent[atom]]
                        atom = local_parent[atom]
                    return atom

                def local_union(left, right):
                    left, right = local_find(left), local_find(right)
                    if left != right:
                        local_parent[right] = left

                for alternative in exact_branch_mappings:
                    if alternative.get(center_R) != center_P:
                        continue
                    if any(r not in alternative for r in neighbors_R):
                        continue
                    alternative_images = tuple(
                        alternative[r] for r in neighbors_R)
                    if set(alternative_images) != set(selected_images):
                        continue
                    for r, image in zip(neighbors_R, alternative_images):
                        local_union(int(r), selected_by_image[image])
                groups = {}
                for atom_R in neighbors_R:
                    groups.setdefault(local_find(int(atom_R)), []).append(
                        int(atom_R))
                for r_atoms in groups.values():
                    if len(r_atoms) <= 1:
                        continue
                    complete.append({
                        'fragment_index': int(fragment.get(
                            'fragment_index', fragment_position)),
                        'island_idx': int(fragment.get(
                            'island_idx', fragment_position)),
                        'center_R': int(center_R),
                        'center_P': int(center_P),
                        'r_atoms': sorted(r_atoms),
                        'p_atoms': sorted(mapping[r] for r in r_atoms),
                        'extendable': False,
                        'open': False,
                        'assignments': 'exact_branch_relation',
                        'source': 'AAM_cross_branch_assignment',
                    })
            continue
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
    # A top-level copy is the minimal isolating copy: the pool only rebinds
    # branch-level keys ('encounter_count', 'cuts'), and every hierarchy
    # reader (_refresh_entry_branch_symmetry, _freeze_analytical,
    # complete_chosen_automorphism_groups, attach_completed_candidate_groups,
    # AAMHierarchy.from_record, index_chirality) copies before it writes, so
    # the nested fragment/block dicts (built fresh by _symmetry_state and
    # _branch_symmetry_record with int/str leaves, never written afterwards)
    # are shared instead of deep-copied.
    return {
        'mapping': {int(r): int(p) for r, p in dict(mapping).items()},
        'cuts': _cut_record(cuts),
        'encounter_count': int(count),
        'hierarchy': dict(branch_symmetry or {
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
    # Only top-level keys are rebound or popped below, and 'blocks' /
    # 'color_groups' are rebuilt as new lists, so a shallow copy isolates the
    # branch hierarchy; the shared nested dicts are never written by anyone.
    hierarchy = dict(
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
    # Branch representatives are provenance, never a sampled permutation
    # group.
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


def _merge_compressed_entry(target, sig, raw_entry, *, take_ownership=False,
                            refresh=True):
    """Merge one compressed entry, optionally consuming the incoming value.

    Without ``take_ownership`` (caller-owned pools from
    ``merge_cut_sweep_pools``) the deep copy stays: the merge rebinds entry
    and branch keys and appends to ``branches``, and the caller may keep
    mutating its own pool afterwards, which cannot be proven safe here.
    """
    entry = target.get(sig)
    if entry is None:
        entry = raw_entry if take_ownership else copy.deepcopy(raw_entry)
        branches = entry.setdefault('branches', [])
        entry['_branch_key_index'] = {
            _analytical_branch_key(branch): branch
            for branch in branches
        }
        target[sig] = entry
        return

    representative_changed = bool(
        raw_entry.get('has_no_cut', False)
        and not entry.get('has_no_cut', False))
    if representative_changed:
        entry['mapping'] = dict(raw_entry['mapping'])
    entry['has_no_cut'] = bool(
        entry.get('has_no_cut', False)
        or raw_entry.get('has_no_cut', False))
    entry['cuts'] = (
        frozenset(entry.get('cuts', ()))
        | frozenset(raw_entry.get('cuts', ())))
    entry['dedup_count'] = (
        int(entry.get('dedup_count', 1))
        + int(raw_entry.get('dedup_count', 1)))
    branches = entry.setdefault('branches', [])
    key_index = entry.get('_branch_key_index')
    if key_index is None:
        key_index = {
            _analytical_branch_key(branch): branch
            for branch in branches
        }
        entry['_branch_key_index'] = key_index
    changed = representative_changed
    for raw_branch in raw_entry.get('branches') or ():
        incoming_branch = (
            raw_branch if take_ownership else copy.deepcopy(raw_branch))
        changed = (_merge_analytical_branch(
            branches, incoming_branch, key_index=key_index) or changed)
    if changed:
        if refresh:
            _refresh_entry_branch_symmetry(entry)
        else:
            entry['_branch_symmetry_dirty'] = True


def _finalize_compressed_pool(pool):
    """Refresh each changed mechanism once after a bulk reduction."""
    for entry in pool.values():
        if entry.pop('_branch_symmetry_dirty', False):
            _refresh_entry_branch_symmetry(entry)
    return _public_pool(pool)


def _merge_compressed_pool(target, incoming, *, take_ownership=False,
                           refresh=True):
    """Merge an already deduplicated partial pool without losing counts."""
    for sig, raw_entry in dict(incoming or {}).items():
        _merge_compressed_entry(
            target, sig, raw_entry,
            take_ownership=take_ownership, refresh=refresh)
    return target


def _candidate_from_symmetry_state(state):
    """Reconstruct one completed compressed candidate from its AAM record."""
    state = dict(state or {})
    witness = {int(r): int(p)
               for r, p in dict(state.get('witness') or {}).items()}
    automorph_blocks = tuple(_SymBlock(
        tuple(map(int, block.get('r_atoms') or ())),
        tuple(map(int, block.get('p_atoms') or ())),
        extendable=False)
        for block in state.get('automorph_blocks') or ())
    blocks = tuple(_SymBlock(
        tuple(map(int, block.get('r_atoms') or ())),
        tuple(map(int, block.get('p_atoms') or ())),
        extendable=bool(block.get('extendable', False)))
        for block in state.get('blocks') or ()
        if str(block.get('source') or '') != 'exact_automorph_group')
    return _SymCand(
        witness, blocks,
        exact_fixed=tuple(map(int, state.get('exact_fixed') or ())),
        multiplicity=int(state.get('multiplicity', 1)),
        automorph_blocks=automorph_blocks)


def attach_completed_candidate_groups(branches, g_P, *, wbo_tol,
                                      node_policy=None, return_metrics=False):
    """Attach exact groups after completed branch-family reduction.

    The cache key is the complete locked prefix plus candidate state.  Live
    candidates and discarded growth paths never enter this bounded stage.
    """
    cache = {}
    canonical_base_cache = {}
    completed = []
    metrics = {
        'completed_candidate_group_requests': 0,
        'completed_candidate_group_calculations': 0,
        'completed_candidate_group_cache_hits': 0,
    }
    for raw_branch in branches:
        # Minimal structural copy: the loop below only rebinds
        # fragment['symmetry'] to a fresh dict, so copy the branch, its
        # hierarchy, the fragments container and each fragment dict; every
        # other nested object is shared read-only with the raw branch.
        branch = dict(raw_branch)
        hierarchy = branch.get('hierarchy') or {}
        raw_fragments = hierarchy.get('fragments') or ()
        if raw_fragments:
            hierarchy = dict(hierarchy)
            hierarchy['fragments'] = type(raw_fragments)(
                dict(fragment) for fragment in raw_fragments)
            branch['hierarchy'] = hierarchy
        locked = {}
        for fragment in hierarchy.get('fragments') or ():
            metrics['completed_candidate_group_requests'] += 1
            state = fragment.get('symmetry') or {}
            candidate = _candidate_from_symmetry_state(state)
            key = (tuple(sorted(locked.items())),
                   _freeze_analytical(state))
            generators = cache.get(key)
            if generators is None:
                metrics['completed_candidate_group_calculations'] += 1
                canonicalizer = _CandidateAutomorphismCanonicalizer(
                    g_P, locked_mapping=locked, node_policy=node_policy,
                    wbo_tol=float(wbo_tol),
                    base_cache=canonical_base_cache)
                generators = canonicalizer.atom_generators(candidate)
                cache[key] = generators
            else:
                metrics['completed_candidate_group_cache_hits'] += 1
            state = dict(state)
            state['automorph_generators'] = [
                list(generator) for generator in generators]
            state['automorph_group_source'] = (
                'completed_candidate_after_branch_family_reduction')
            fragment['symmetry'] = state
            for r, p in candidate.mapping.items():
                prior = locked.get(int(r))
                if prior is not None and prior != int(p):
                    raise RuntimeError(
                        "completed AAM fragment conflicts with locked prefix")
                locked[int(r)] = int(p)
        completed.append(branch)
    return (completed, metrics) if return_metrics else completed


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
                  *, return_trace=False, collect_metrics=False):
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
        orders = _generate_seed_orders(g_R, n_trials=int(cfg['n_seeds']))
    elif isinstance(orders, int):
        # Seed order ``idx`` depends only on the shared base ordering and its
        # own generator (rng_seed + idx + 1), so building it alone equals
        # ``_generate_seed_orders(g_R, n_seeds)[idx]``; ``range`` keeps the
        # original list-index semantics (negative index, IndexError).
        orders = [_generate_seed_order(
            g_R, range(int(cfg['n_seeds']))[int(orders)])]
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
    compact_metrics = {
        'cuts': 1,
        'seed_orders': 0,
        'growth_calls': 0,
        'subtree_branch_cap_count': 0,
        'max_live_branches': 0,
        'max_growth_candidates': 0,
    }
    try:
        for order_index, order in enumerate(orders):
            seed_t0 = time.perf_counter()
            seed_growth_profile = [] if (return_trace or collect_metrics) else None
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
            if collect_metrics:
                compact_metrics['seed_orders'] += 1
                compact_metrics['growth_calls'] += int(
                    growth_summary['calls'])
                compact_metrics['max_growth_candidates'] = max(
                    compact_metrics['max_growth_candidates'],
                    int(growth_summary['max_cands_before']),
                    int(growth_summary['max_cands_after']))
                compact_metrics['subtree_branch_cap_count'] += sum(
                    item.get('result') == 'subtree_branch_cap'
                    for item in seed_growth_profile)
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
            if collect_metrics:
                compact_metrics['max_live_branches'] = max(
                    compact_metrics['max_live_branches'], n_branches)
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
    compact_metrics['raw_result_count'] = len(out)
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
    return out, events, compact_metrics


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
    cut, orders, core_R, trace_enabled, metrics_enabled = args
    cfg = _WORKER['cfg']
    out, events, metrics = _run_cut_work(
        _WORKER['elR'], _WORKER['wboR'],
        _WORKER['elT'], _WORKER['wboT'],
        cfg, cut, orders, core_R,
        _WORKER['g_P'], _WORKER['g_R_full'],
        _WORKER['p_orbits'], _WORKER['r_orbits'],
        return_trace=trace_enabled, collect_metrics=metrics_enabled)
    local_pool = {}
    for result in out:
        sig, mapping_items, cut = result[:3]
        branch_symmetry = result[3] if len(result) > 3 else None
        _pool_add(
            local_pool, sig, dict(mapping_items), cut, branch_symmetry)
    local_pool = _public_pool(local_pool)
    metrics['worker_returned_branch_count'] = sum(
        len(entry.get('branches') or ())
        for entry in local_pool.values())
    return {'pool': local_pool, 'events': events, 'metrics': metrics}


def _merge_sweep_metrics(target, source):
    target['cuts'] += int(source.get('cuts', 0))
    target['seed_orders'] += int(source.get('seed_orders', 0))
    target['growth_calls'] += int(source.get('growth_calls', 0))
    target['subtree_branch_cap_count'] += int(
        source.get('subtree_branch_cap_count', 0))
    target['max_live_branches'] = max(
        target['max_live_branches'],
        int(source.get('max_live_branches', 0)))
    target['max_growth_candidates'] = max(
        target['max_growth_candidates'],
        int(source.get('max_growth_candidates', 0)))
    target['raw_result_count'] += int(source.get('raw_result_count', 0))
    target['worker_returned_branch_count'] += int(
        source.get('worker_returned_branch_count', 0))


def _new_sweep_metrics(max_branches):
    return {
        'configured_max_branches': int(max_branches),
        'cuts': 0,
        'seed_orders': 0,
        'growth_calls': 0,
        'subtree_branch_cap_count': 0,
        'max_live_branches': 0,
        'max_growth_candidates': 0,
        'raw_result_count': 0,
        'worker_returned_branch_count': 0,
        # Parent-side phases of the parallel sweep (see _stream_merge_pool):
        #   parent_route_seconds    incremental in-parent merge of worker pools
        #   parallel_reduce_seconds one-shot refresh of changed mechanisms
        #   parent_load_seconds     intermediate-file persistence I/O
        #   worker_stream_seconds   time spent waiting on worker results
        #   parent_merge_seconds    sum of the three parent-side phases
        'parent_merge_seconds': 0.0,
        'parent_route_seconds': 0.0,
        'parallel_reduce_seconds': 0.0,
        'parent_load_seconds': 0.0,
        'worker_stream_seconds': 0.0,
    }


def _stream_merge_pool(payloads, *, metrics, directory=None,
                       persistent=False, after_stream=None):
    """Merge ordered worker pools into the parent pool as they arrive.

    Results are consumed with ordered ``imap`` (chunksize 1), so merging each
    payload's entries in arrival order is exactly the serial ``_pool_add``
    order: a mechanism's representative and branch list depend only on the
    relative order of that mechanism's entries, which is the work order in
    both cases.  No bucket files, second process pool, or reload remain; with
    ``persistent`` the raw entry stream and the merged pool are still written
    next to a manifest for inspection.
    """
    pool = {}
    entry_count = 0
    raw_paths = []
    raw_handle = None
    if persistent:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        raw_paths = [directory / 'raw_bucket_000.pkl']
        raw_handle = raw_paths[0].open('xb')
    stream_started = time.perf_counter()
    merge_seconds = 0.0
    persist_seconds = 0.0
    try:
        for payload in payloads:
            _merge_sweep_metrics(metrics, payload.get('metrics') or {})
            entries = payload.get('pool') or {}
            if raw_handle is not None:
                persist_started = time.perf_counter()
                for sig, entry in entries.items():
                    pickle.dump((sig, entry), raw_handle,
                                protocol=pickle.HIGHEST_PROTOCOL)
                persist_seconds += time.perf_counter() - persist_started
            merge_started = time.perf_counter()
            for sig, entry in entries.items():
                _merge_compressed_entry(
                    pool, sig, entry, take_ownership=True, refresh=False)
                entry_count += 1
            merge_seconds += time.perf_counter() - merge_started
    finally:
        if raw_handle is not None:
            raw_handle.close()
    stream_elapsed = time.perf_counter() - stream_started
    if after_stream is not None:
        after_stream()

    finalize_started = time.perf_counter()
    _finalize_compressed_pool(pool)
    finalize_seconds = time.perf_counter() - finalize_started

    reduced_paths = []
    if persistent:
        persist_started = time.perf_counter()
        reduced_paths = [directory / 'reduced_bucket_000.pkl']
        with reduced_paths[0].open('xb') as handle:
            pickle.dump(pool, handle, protocol=pickle.HIGHEST_PROTOCOL)
        persist_seconds += time.perf_counter() - persist_started

    metrics['parent_route_seconds'] += merge_seconds
    metrics['worker_stream_seconds'] += max(
        0.0, stream_elapsed - merge_seconds - persist_seconds)
    metrics['parallel_reduce_seconds'] += finalize_seconds
    metrics['parent_load_seconds'] += persist_seconds
    metrics['parent_merge_seconds'] = (
        metrics['parent_route_seconds']
        + metrics['parallel_reduce_seconds']
        + metrics['parent_load_seconds'])

    if persistent:
        manifest = {
            'schema': 'rxn_core.cut_sweep_intermediate/v1',
            'bucket_count': 1,
            'worker_pool_entries': entry_count,
            'mechanism_count': len(pool),
            'raw_buckets': [path.name for path in raw_paths],
            'reduced_buckets': [path.name for path in reduced_paths],
            'metrics': dict(metrics),
        }
        (directory / 'manifest.json').write_text(
            json.dumps(manifest, indent=2) + '\n')
    return pool


def cut_sweep_items(wboR, cut_floor=0.2, *, elements=None,
                    heavy_only=False):
    """Return independent no-cut plus one-edge cut work items.

    ``heavy_only`` keeps explicit hydrogens in AAM mapping and final scoring;
    it only avoids using X-H bonds as artificial graph-search cuts.
    """
    edges = _strong_edges(wboR, float(cut_floor))
    if heavy_only:
        if elements is None:
            raise ValueError("elements are required for heavy-only cut sweep")
        edges = [
            (i, j) for i, j in edges
            if str(elements[i]) != 'H' and str(elements[j]) != 'H'
        ]
    return [()] + [((int(i), int(j)),) for i, j in edges]


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
                            trace_path=None, collect_metrics=False):
    graph_floor = float(cfg['graph_floor'])
    g_P = build_graph(elT, wboT, bond_cut=graph_floor)
    g_R_full = build_graph(elR, wboR, bond_cut=graph_floor)
    symmetry_wbo_tol = float(cfg['symmetry_wbo_tol'])
    p_orbits = _nauty_orbits(g_P, wbo_tol=symmetry_wbo_tol)
    r_orbits = _nauty_orbits(g_R_full, wbo_tol=symmetry_wbo_tol)
    pool = {}
    metrics = _new_sweep_metrics(cfg['max_branches'])

    for cut in cuts:
        cut = tuple(tuple(pair) for pair in cut)
        results, events, cut_metrics = _run_cut_work(
            elR, wboR, elT, wboT, cfg, cut, None, core_R,
            g_P, g_R_full, p_orbits, r_orbits,
            return_trace=bool(trace_path), collect_metrics=collect_metrics)
        if collect_metrics:
            _merge_sweep_metrics(metrics, cut_metrics)
        _emit_trace(trace_path, events)
        for result in results:
            sig, mapping_items, _cut = result[:3]
            branch_symmetry = result[3] if len(result) > 3 else None
            mapping = dict(mapping_items)
            _pool_add(pool, sig, mapping, _cut, branch_symmetry)
    return pool, metrics


def _cut_sweep_serial(elR, wboR, elT, wboT, cfg, core_R,
                      collect_metrics=False):
    pool, metrics = _cut_sweep_chunk_serial(
        elR, wboR, elT, wboT, cfg, core_R,
        cut_sweep_items(
            wboR, cfg['cut_floor'], elements=elR,
            heavy_only=cfg.get('heavy_cuts_only', False)),
        collect_metrics=collect_metrics)
    public = _public_pool(pool)
    return (public, metrics) if collect_metrics else public


def _cut_sweep_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R,
                        collect_metrics=False, intermediate_dir=None):
    cuts = cut_sweep_items(
        wboR, cfg['cut_floor'], elements=elR,
        heavy_only=cfg.get('heavy_cuts_only', False))
    work = [
        (cut, seed_index, core_R, False, collect_metrics)
        for cut in cuts
        for seed_index in range(int(cfg['n_seeds']))
    ]
    metrics = _new_sweep_metrics(cfg['max_branches'])
    proc_pool = mp.Pool(n_workers, initializer=_cs_winit,
                        initargs=(elR, wboR, elT, wboT, cfg))
    pool_closed = False

    def finish_search_pool():
        nonlocal pool_closed
        if not pool_closed:
            proc_pool.__exit__(None, None, None)
            pool_closed = True

    try:
        # Workers may finish in any order, but pool insertion order determines
        # the representative retained for analytically equivalent families.
        # Consume results in the explicit cut/seed work order so parallel
        # scheduling cannot change downstream chirality/RMSD representatives.
        payloads = proc_pool.imap(
            _cs_wrun, work, chunksize=max(1, int(cfg['chunksize'])))
        pool = _stream_merge_pool(
            payloads, metrics=metrics, directory=intermediate_dir,
            persistent=intermediate_dir is not None,
            after_stream=finish_search_pool)
    finally:
        finish_search_pool()
    metrics['cuts'] = len(cuts)
    public = _public_pool(pool)
    return (public, metrics) if collect_metrics else public


def _cut_sweep_chunk_parallel(elR, wboR, elT, wboT, cfg, n_workers, core_R,
                              cuts, trace_path=None, intermediate_dir=None,
                              collect_metrics=False):
    work = [
        (cut, seed_index, core_R, bool(trace_path), collect_metrics)
        for cut in cuts
        for seed_index in range(int(cfg['n_seeds']))
    ]
    metrics = _new_sweep_metrics(cfg['max_branches'])
    proc_pool = mp.Pool(n_workers, initializer=_cs_winit,
                        initargs=(elR, wboR, elT, wboT, cfg))
    pool_closed = False

    def finish_search_pool():
        nonlocal pool_closed
        if not pool_closed:
            proc_pool.__exit__(None, None, None)
            pool_closed = True

    try:
        # Preserve the same logical cut/seed order as serial execution.  The
        # work remains parallel; only parent-side result consumption is
        # deterministic.
        def payloads():
            for payload in proc_pool.imap(
                    _cs_wrun, work,
                    chunksize=max(1, int(cfg['chunksize']))):
                _emit_trace(trace_path, payload.get('events', []))
                yield payload

        pool = _stream_merge_pool(
            payloads(), metrics=metrics, directory=intermediate_dir,
            persistent=intermediate_dir is not None,
            after_stream=finish_search_pool)
    finally:
        finish_search_pool()
    metrics['cuts'] = len(cuts)
    metrics['seed_orders'] = len(work)
    public = _public_pool(pool)
    return (public, metrics) if collect_metrics else public


def run_cut_sweep_chunk(elR, wboR, elT, wboT, cuts, *,
                        core_R=None,
                        n_workers=None,
                        trace_path=None,
                        intermediate_dir=None,
                        return_metrics=False,
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
    work_count = len(normalized_cuts) * int(cfg['n_seeds'])
    if n_workers and int(n_workers) > 1 and work_count > 1:
        return _cut_sweep_chunk_parallel(
            elR, wboR, elT, wboT, cfg,
            min(int(n_workers), work_count),
            core_R, normalized_cuts, trace_path=trace_path,
            intermediate_dir=intermediate_dir,
            collect_metrics=bool(return_metrics))
    pool, metrics = _cut_sweep_chunk_serial(
        elR, wboR, elT, wboT, cfg, core_R, normalized_cuts,
        trace_path=trace_path, collect_metrics=bool(return_metrics))
    public = _public_pool(pool)
    return (public, metrics) if return_metrics else public


def merge_cut_sweep_pools(pools):
    """Merge partial cut-sweep pools produced by chunk tasks."""
    merged = {}
    for pool in pools:
        for info in dict(pool or {}).values():
            if not info.get('branches'):
                raise ValueError(
                    "cut-sweep chunk lacks analytical AAM branches; rerun it "
                    "with the current branch schema")
        _merge_compressed_pool(merged, pool)
    return _public_pool(merged)


def cut_sweep(elR, wboR, elT, wboT, *,
              n_workers=None, core_R=None,
              intermediate_dir=None,
              heavy_cuts_only=False,
              cut_floor=0.2, graph_floor=0.2, iso_tol=1.0,
              dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
              symmetry_wbo_tol=0.2,
              n_seeds=3, max_branches=100,
              chunksize=1,
              symmetry_repair=True,
              symmetry_repair_min_changes=1,
              symmetry_repair_max_evals=20000,
              anchor_map=None,
              return_metrics=False):
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
    cfg['heavy_cuts_only'] = bool(heavy_cuts_only)
    core_R = tuple(sorted(set(core_R or ())))
    if not n_workers or n_workers <= 1:
        return _cut_sweep_serial(
            elR, wboR, elT, wboT, cfg, core_R,
            collect_metrics=bool(return_metrics))
    return _cut_sweep_parallel(elR, wboR, elT, wboT, cfg,
                               int(n_workers), core_R,
                               collect_metrics=bool(return_metrics),
                               intermediate_dir=intermediate_dir)


def select_min_mechanisms(pool):
    """Keep only signatures with the fewest broken+formed orbit events."""
    if not pool:
        return {}
    best = min(len(sig[0]) + len(sig[1]) for sig in pool)
    return {
        sig: info for sig, info in pool.items()
        if len(sig[0]) + len(sig[1]) == best
    }

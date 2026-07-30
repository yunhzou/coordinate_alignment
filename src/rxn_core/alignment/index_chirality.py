"""Exact index-orientation consensus for one analytical AAM branch.

The selected AAM mapping is never recomputed here.  This module only composes
it with an automorphism shared by the selected, floor-masked R and P
fragments.  Local handedness is represented by affine substituent simplices,
so coordination number and display symmetry blocks play no role.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Mapping, Sequence

import numpy as np

from ..frag import bond_event_threshold, build_graph
from ..matcher.orbits import _wbo_tolerance_bucket_lookup


INDEX_CHIRALITY_SCHEMA = "rxn_core.index_chirality/v3"
ORIENTATION_DEGENERACY_TOL = 0.1
# Group orientation is topological: only a determinant indistinguishable from
# zero at working precision is undefined.  The local-center tolerance above
# remains appropriate for ordinary near-planar stereocenters.
GROUP_ORIENTATION_DEGENERACY_TOL = 0.0


class IndexChiralityError(ValueError):
    """Invalid index-orientation input or selected-fragment metadata."""


class IndexChiralityConflict(IndexChiralityError):
    """No exact allowed automorphism gives endpoint orientation consensus."""

    def __init__(self, message, *, diagnostics=None):
        super().__init__(message)
        self.diagnostics = diagnostics


@dataclass(frozen=True)
class _OrientedFrame:
    center_R: int
    neighbors_R: tuple[int, ...]
    sign_R: int
    sign_P_source: int
    normalized_R: float
    normalized_P_source: float

    @property
    def frame_id(self):
        shell = "-".join(str(value) for value in self.neighbors_R)
        return f"simplex:{self.center_R}:{shell}"


@dataclass(frozen=True)
class IndexChiralitySelection:
    source_mapping: dict[int, int]
    selected_mapping: dict[int, int]
    metadata: dict


@dataclass(frozen=True)
class GroupChiralityWitnessSelection:
    selected_mapping: dict[int, int]
    selected_witness_index: int | None
    preserved_frames: tuple[dict, ...]
    metadata: dict


@dataclass(frozen=True)
class GroupChiralityBranchAnalysis:
    mapping: dict[int, int]
    defined_frames: tuple[dict, ...]
    metadata: dict


@dataclass(frozen=True)
class _OrientationMeasure:
    normalized: float
    determinant: float
    determinant_error_bound: float
    sign: int
    zero_length: bool


class _RelationalGraph:
    """Two equally shaped colored graphs with shared incidence topology."""

    def __init__(self, atom_colors_A, atom_colors_B):
        if len(atom_colors_A) != len(atom_colors_B):
            raise IndexChiralityError("relational atom color lengths differ")
        self.colors_A = list(atom_colors_A)
        self.colors_B = list(atom_colors_B)
        self.adjacency = [set() for _ in atom_colors_A]
        self.atom_count = len(atom_colors_A)
        self.relation_records_A = []
        self.relation_records_B = []

    def _vertex(self, color_A, color_B=None):
        index = len(self.adjacency)
        self.adjacency.append(set())
        self.colors_A.append(color_A)
        self.colors_B.append(color_A if color_B is None else color_B)
        return index

    def _edge(self, left, right):
        self.adjacency[left].add(right)
        self.adjacency[right].add(left)

    def add_pair(self, left, right, color_A, color_B=None):
        color_B = color_A if color_B is None else color_B
        relation = self._vertex(color_A, color_B)
        self._edge(relation, int(left))
        self._edge(relation, int(right))
        atoms = tuple(sorted((int(left), int(right))))
        self.relation_records_A.append(("pair", repr(color_A), atoms))
        self.relation_records_B.append(("pair", repr(color_B), atoms))

    def add_ordered_relation(self, atoms, color_A, color_B):
        atoms = tuple(map(int, atoms))
        relation = self._vertex(color_A, color_B)
        for role, atom in enumerate(atoms):
            role_vertex = self._vertex(("orientation_role", role))
            self._edge(relation, role_vertex)
            self._edge(role_vertex, int(atom))
        self.relation_records_A.append(
            ("ordered", repr(color_A), atoms))
        self.relation_records_B.append(
            ("ordered", repr(color_B), atoms))

    def clone(self):
        result = _RelationalGraph([], [])
        result.colors_A = list(self.colors_A)
        result.colors_B = list(self.colors_B)
        result.adjacency = [set(neighbors) for neighbors in self.adjacency]
        result.atom_count = self.atom_count
        result.relation_records_A = list(self.relation_records_A)
        result.relation_records_B = list(self.relation_records_B)
        return result

    def graph(self, side, *, individualized=()):
        try:
            import pynauty
        except ImportError as exc:
            raise RuntimeError(
                "pynauty is required for index-orientation consensus"
            ) from exc
        colors = list(self.colors_A if side == "A" else self.colors_B)
        for rank, vertex in enumerate(individualized):
            vertex = int(vertex)
            colors[vertex] = (colors[vertex], "individualized", rank)
        classes = defaultdict(set)
        for vertex, color in enumerate(colors):
            classes[color].add(vertex)
        coloring = [
            set(vertices) for _, vertices in
            sorted(classes.items(), key=lambda item: repr(item[0]))
        ]
        return pynauty.Graph(
            len(self.adjacency), directed=False,
            adjacency_dict={
                vertex: sorted(neighbors)
                for vertex, neighbors in enumerate(self.adjacency)
            },
            vertex_coloring=coloring,
        )


def _int_mapping(mapping: Mapping[int, int]) -> dict[int, int]:
    return {int(r): int(p) for r, p in dict(mapping).items()}


def validate_mapping(mapping, elements_R, elements_P):
    mapping = _int_mapping(mapping)
    atom_count = len(elements_R)
    if len(elements_P) != atom_count:
        raise IndexChiralityError("endpoint atom counts differ")
    if set(mapping) != set(range(atom_count)):
        raise IndexChiralityError("mapping is not complete")
    if set(mapping.values()) != set(range(atom_count)):
        raise IndexChiralityError("mapping is not bijective")
    for r, p in mapping.items():
        if str(elements_R[r]) != str(elements_P[p]):
            raise IndexChiralityError(f"element mismatch at R{r}->P{p}")
    return mapping


def fixed_mapping_aligned_rmsd(mapping, coords_R, coords_P):
    """RMSD after a proper rigid fit of one immutable atom mapping.

    Product coordinates are first reindexed exactly once by ``mapping``.
    Kabsch removes only global translation and proper rotation; it performs no
    permutation, assignment, symmetry shuffle, or correspondence search.
    """
    mapping = _int_mapping(mapping)
    reactant = np.asarray(coords_R, dtype=float)
    product = np.asarray(
        [coords_P[mapping[r]] for r in range(len(reactant))], dtype=float)
    if product.shape != reactant.shape:
        raise IndexChiralityError("fixed-mapping RMSD coordinate mismatch")
    centered_P = product - product.mean(axis=0)
    centered_R = reactant - reactant.mean(axis=0)
    u, _singular, vt = np.linalg.svd(centered_P.T @ centered_R)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0.0:
        u[:, -1] *= -1.0
        rotation = u @ vt
    aligned = centered_P @ rotation
    return float(np.sqrt(np.mean(np.sum(
        (aligned - centered_R) ** 2, axis=1))))


def _fixed_mappings_aligned_rmsd(mappings, coords_R, coords_P):
    """Vectorized equivalent of :func:`fixed_mapping_aligned_rmsd`.

    Every product coordinate array is indexed by its supplied immutable
    mapping before a batched proper Kabsch fit.  This changes only scheduling:
    there is no correspondence search or symmetry rematching.
    """
    mappings = list(mappings)
    if not mappings:
        return np.empty(0, dtype=float)
    reactant = np.asarray(coords_R, dtype=float)
    product_xyz = np.asarray(coords_P, dtype=float)
    atom_count = len(reactant)
    indices = np.asarray([
        [int(mapping[r]) for r in range(atom_count)]
        for mapping in mappings
    ], dtype=int)
    products = product_xyz[indices]
    centered_P = products - products.mean(axis=1, keepdims=True)
    centered_R = reactant - reactant.mean(axis=0)
    covariance = np.einsum(
        'kni,nj->kij', centered_P, centered_R, optimize=True)
    u, _singular, vt = np.linalg.svd(covariance)
    rotations = u @ vt
    reflected = np.linalg.det(rotations) < 0.0
    if np.any(reflected):
        u = u.copy()
        u[reflected, :, -1] *= -1.0
        rotations = u @ vt
    aligned = centered_P @ rotations
    return np.sqrt(np.mean(np.sum(
        (aligned - centered_R[None, :, :]) ** 2, axis=2), axis=1))


def _orientation_measure(coords, origin, other_points, *,
                         degeneracy_tol=ORIENTATION_DEGENERACY_TOL):
    xyz = np.asarray(coords, dtype=np.longdouble)
    vectors = np.stack(
        [xyz[int(point)] - xyz[int(origin)] for point in other_points],
        axis=0,
    )
    squared_lengths = np.sum(vectors * vectors, axis=1)
    zero_length = bool(np.any(squared_lengths == 0))
    denominator = np.sqrt(np.prod(squared_lengths))
    a, b, c = vectors
    positive = (
        a[0] * b[1] * c[2],
        a[1] * b[2] * c[0],
        a[2] * b[0] * c[1],
    )
    negative = (
        a[2] * b[1] * c[0],
        a[1] * b[0] * c[2],
        a[0] * b[2] * c[1],
    )
    determinant = sum(positive) - sum(negative)
    permanent = sum(abs(term) for term in (*positive, *negative))
    eps = np.longdouble(np.finfo(np.longdouble).eps)
    gamma = 16 * eps / (1 - 16 * eps)
    error_bound = gamma * max(permanent, denominator)
    normalized = 0 if zero_length else determinant / denominator
    if (zero_length or abs(determinant) <= error_bound
            or abs(normalized) <= float(degeneracy_tol)):
        sign = 0
    else:
        sign = 1 if determinant > 0 else -1
    return _OrientationMeasure(
        normalized=float(normalized),
        determinant=float(determinant),
        determinant_error_bound=float(error_bound),
        sign=sign,
        zero_length=zero_length,
    )


def mapping_event_signature(mapping, wbo_R, wbo_P, elements_R, *,
                            dwbo_threshold=0.5,
                            metal_dwbo_threshold=0.3):
    """Exact broken/formed R-index pairs for one complete mapping."""
    mapping = _int_mapping(mapping)
    broken, formed = [], []
    for left in range(len(elements_R)):
        for right in range(left + 1, len(elements_R)):
            threshold = bond_event_threshold(
                elements_R, left, right,
                default_threshold=float(dwbo_threshold),
                metal_threshold=metal_dwbo_threshold)
            delta = float(
                wbo_R[left, right]
                - wbo_P[mapping[left], mapping[right]])
            if delta >= threshold:
                broken.append((left, right))
            elif -delta >= threshold:
                formed.append((left, right))
    return tuple(broken), tuple(formed)


def _event_class(r_wbo, p_wbo, threshold):
    delta = float(r_wbo) - float(p_wbo)
    if delta >= threshold:
        return 1
    if -delta >= threshold:
        return -1
    return 0


def _selected_fragments(branch_symmetry, source, fixed_singletons=()):
    raw = list(dict(branch_symmetry or {}).get("fragments") or ())
    if not raw:
        raw = [{"fragment_index": 0, "fragment": sorted(source)}]
    fragments, owner = [], {}
    for position, record in enumerate(raw):
        atoms = tuple(sorted(int(atom) for atom in record.get("fragment", ())))
        if not atoms:
            raise IndexChiralityError("selected AAM fragment is empty")
        for atom in atoms:
            if atom not in source:
                raise IndexChiralityError(
                    f"selected fragment contains unknown R atom {atom}")
            if atom in owner:
                raise IndexChiralityError(
                    f"selected fragments overlap at R atom {atom}")
            owner[atom] = position
        fragments.append(atoms)
    missing = sorted(set(source) - set(owner))
    fixed_singletons = set(map(int, fixed_singletons or ()))
    uncovered = sorted(set(missing) - fixed_singletons)
    if uncovered:
        raise IndexChiralityError(
            f"selected AAM fragments do not cover R atoms {uncovered}")
    # Hard anchors are preloaded into the branch before island growth.  If an
    # anchor never participates in a grown fragment, it is still a complete
    # individually fixed analytical fragment rather than missing hierarchy.
    for atom in missing:
        owner[atom] = len(fragments)
        fragments.append((atom,))
    return tuple(fragments), owner


def analytical_family_static_context(
        elements_R, wbo_R, elements_P, wbo_P, *, graph_floor=0.2,
        dwbo_threshold=0.5, metal_dwbo_threshold=0.3):
    """Precompute endpoint-only data shared by every analytical branch."""
    wbo_R = np.asarray(wbo_R)
    wbo_P = np.asarray(wbo_P)
    atom_count = len(elements_R)
    pair_groups_R = defaultdict(list)
    pair_groups_P = defaultdict(list)
    records_R, records_P = {}, {}
    for side, elements, wbo, groups, records in (
            ('R', elements_R, wbo_R, pair_groups_R, records_R),
            ('P', elements_P, wbo_P, pair_groups_P, records_P)):
        for left in range(atom_count):
            for right in range(left + 1, atom_count):
                element_pair = tuple(sorted((
                    str(elements[left]), str(elements[right]))))
                threshold = bond_event_threshold(
                    elements, left, right,
                    default_threshold=float(dwbo_threshold),
                    metal_threshold=metal_dwbo_threshold)
                group_key = (element_pair, float(threshold))
                value = float(wbo[left, right])
                groups[group_key].append(value)
                records[(left, right)] = (group_key, value)
    values_R = {
        key: tuple(sorted(set(values)))
        for key, values in pair_groups_R.items()
    }
    values_P = {
        key: tuple(sorted(set(values)))
        for key, values in pair_groups_P.items()
    }
    behavior_R = {}
    for pair, (group_key, value) in records_R.items():
        threshold = group_key[1]
        behavior_R[pair] = tuple(
            _event_class(value, other, threshold)
            for other in values_P[group_key])
    behavior_P = {}
    for pair, (group_key, value) in records_P.items():
        threshold = group_key[1]
        behavior_P[pair] = tuple(
            _event_class(other, value, threshold)
            for other in values_R[group_key])
    return {
        'graph_R': build_graph(elements_R, wbo_R, bond_cut=graph_floor),
        'graph_P': build_graph(elements_P, wbo_P, bond_cut=graph_floor),
        'pair_records_R': records_R,
        'pair_records_P': records_P,
        'behavior_R': behavior_R,
        'behavior_P': behavior_P,
    }


def _masked_relation_data(source, branch_symmetry, elements_R, wbo_R,
                          elements_P, wbo_P, graph_floor,
                          symmetry_wbo_tol, dwbo_threshold,
                          metal_dwbo_threshold, anchor_map,
                          static_context=None):
    atom_count = len(source)
    inverse = {p: r for r, p in source.items()}
    anchors = {int(r): int(p) for r, p in dict(anchor_map or {}).items()}
    for r, p in anchors.items():
        if r not in source or p not in inverse:
            raise IndexChiralityError(f"invalid anchor R{r}->P{p}")
        if str(elements_R[r]) != str(elements_P[p]):
            raise IndexChiralityError(f"anchor element mismatch R{r}->P{p}")
    fragments, owner_R = _selected_fragments(
        branch_symmetry, source, fixed_singletons=anchors)
    owner_P = {source[r]: owner_R[r] for r in source}

    source_anchor_tags = defaultdict(list)
    target_anchor_tags = defaultdict(list)
    for r, p in sorted(anchors.items()):
        source_anchor_tags[source[r]].append(r)
        target_anchor_tags[p].append(r)
    atom_colors_A, atom_colors_B = [], []
    for p in range(atom_count):
        common = ("atom", str(elements_P[p]), owner_P[p])
        atom_colors_A.append(common + (tuple(source_anchor_tags[p]),))
        atom_colors_B.append(common + (tuple(target_anchor_tags[p]),))
    relation = _RelationalGraph(atom_colors_A, atom_colors_B)

    static_context = static_context or analytical_family_static_context(
        elements_R, wbo_R, elements_P, wbo_P,
        graph_floor=graph_floor,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold)
    g_R = static_context['graph_R']
    g_P = static_context['graph_P']
    persistent_P = {p: set() for p in range(atom_count)}
    for fragment_id, fragment_R in enumerate(fragments):
        fragment_P = tuple(source[r] for r in fragment_R)
        sub_R = g_R.subgraph(fragment_R).copy()
        sub_P = g_P.subgraph(fragment_P).copy()
        buckets_R, zero_R = _wbo_tolerance_bucket_lookup(
            sub_R, float(symmetry_wbo_tol))
        buckets_P, zero_P = _wbo_tolerance_bucket_lookup(
            sub_P, float(symmetry_wbo_tol))
        for (left_R, right_R), bucket in sorted(buckets_R.items()):
            if bucket != zero_R:
                relation.add_pair(
                    source[left_R], source[right_R],
                    ("masked_R_wbo", fragment_id, bucket))
        for (left_P, right_P), bucket in sorted(buckets_P.items()):
            if bucket != zero_P:
                relation.add_pair(
                    left_P, right_P,
                    ("masked_P_wbo", fragment_id, bucket))
        for left_R in fragment_R:
            for right_R in fragment_R:
                if left_R >= right_R or not sub_R.has_edge(left_R, right_R):
                    continue
                left_P, right_P = source[left_R], source[right_R]
                if sub_P.has_edge(left_P, right_P):
                    persistent_P[left_P].add(right_P)
                    persistent_P[right_P].add(left_P)

    # These complete-pair colors are the coarsest endpoint scalar relations
    # needed to make the broken/formed event classification invariant.  They
    # split an iso-tolerance class only when crossing an actual event boundary.
    pair_records = []
    for left_P in range(atom_count):
        for right_P in range(left_P + 1, atom_count):
            left_R, right_R = inverse[left_P], inverse[right_P]
            pair_R = tuple(sorted((left_R, right_R)))
            pair_P = (left_P, right_P)
            group_key = static_context['pair_records_R'][pair_R][0]
            if static_context['pair_records_P'][pair_P][0] != group_key:
                raise IndexChiralityError(
                    "mapping changed an endpoint element-pair class")
            pair_records.append((left_P, right_P, group_key,
                                 static_context['behavior_R'][pair_R],
                                 static_context['behavior_P'][pair_P]))
    colored_pair_records = []
    color_counts = defaultdict(Counter)
    for left_P, right_P, group_key, r_behavior, p_behavior in pair_records:
        color = (
            "event_invariant_pair", group_key, r_behavior, p_behavior)
        colored_pair_records.append((left_P, right_P, group_key, color))
        color_counts[group_key][color] += 1

    # This is a complete edge-colored relation within each element/threshold
    # pair class.  Its most frequent color can be represented by absence:
    # atom colors already preserve the pair class, and any permutation that
    # preserves every exceptional colored pair necessarily maps the remaining
    # (baseline) pairs among themselves.  This is exactly equivalent to
    # materializing O(N^2) relation vertices, while molecular graphs normally
    # retain only O(E) exceptional event relations.
    baseline_color = {
        group_key: min(
            counts,
            key=lambda color: (-counts[color], repr(color)),
        )
        for group_key, counts in color_counts.items()
    }
    for left_P, right_P, group_key, color in colored_pair_records:
        if color == baseline_color[group_key]:
            continue
        relation.add_pair(
            left_P, right_P,
            color)

    return relation, persistent_P, inverse, fragments


def _canonical_isomorphism(graph_A, graph_B):
    import pynauty

    if pynauty.certificate(graph_A) != pynauty.certificate(graph_B):
        return None
    label_A = pynauty.canon_label(graph_A)
    label_B = pynauty.canon_label(graph_B)
    return {int(label_A[i]): int(label_B[i])
            for i in range(len(label_A))}


def _generated_atom_permutations(raw_generators, degree):
    """Enumerate the exact distinct atom action generated by pynauty.

    Relation vertices can permute while every atom stays fixed.  Restricting
    generators to the atom prefix before closure ensures such kernel elements
    never become RMSD candidates.
    """
    degree = int(degree)
    identity = tuple(range(degree))
    generators = tuple(dict.fromkeys(
        tuple(int(generator[atom]) for atom in range(degree))
        for generator in raw_generators
    ))
    generators = tuple(generator for generator in generators
                       if generator != identity)
    seen = {identity}
    queue = [identity]
    while queue:
        current = queue.pop()
        for generator in generators:
            image = tuple(generator[current[atom]]
                          for atom in range(degree))
            if image not in seen:
                seen.add(image)
                queue.append(image)
    return tuple(sorted(seen))


def _independent_atom_action_factors(raw_generators, degree):
    """Factor an atom action into exact disjoint-support subgroups.

    Generator supports that overlap belong to the same factor.  Distinct
    resulting factors have disjoint supports, commute, and intersect only in
    the identity, so their Cartesian product is the original atom action.
    Only each local subgroup is enumerated; the global product is not.
    """
    degree = int(degree)
    identity = tuple(range(degree))
    generators = tuple(dict.fromkeys(
        tuple(int(generator[atom]) for atom in range(degree))
        for generator in raw_generators
    ))
    generators = tuple(generator for generator in generators
                       if generator != identity)
    supports = [
        {atom for atom, image in enumerate(generator) if atom != image}
        for generator in generators
    ]
    parent = list(range(degree))

    def find(atom):
        while parent[atom] != atom:
            parent[atom] = parent[parent[atom]]
            atom = parent[atom]
        return atom

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for support in supports:
        support = tuple(sorted(support))
        for atom in support[1:]:
            union(support[0], atom)
    component_generators = defaultdict(list)
    for generator, support in zip(generators, supports):
        component_generators[find(min(support))].append(generator)
    factors = []
    for component in component_generators.values():
        local_actions = _generated_atom_permutations(component, degree)
        support = tuple(sorted({
            atom for action in local_actions
            for atom, image in enumerate(action) if atom != image
        }))
        factors.append((support, local_actions))
    return tuple(sorted(factors, key=lambda item: item[0]))


def _minimum_rmsd_group_action(canonical_mapping, raw_generators,
                               coords_R, coords_P):
    """Find the exact minimum-RMSD atom action without global enumeration.

    Independent symmetry factors form a search tree.  For every partial
    assignment, pair-distance disagreement among already fixed atoms gives a
    rigorous lower bound on any proper-fit RMSD below that node:

        RMSD >= sqrt(sum_ij (dR_ij - dP_ij)^2) / N.

    Therefore a complete remaining coset can be skipped when its bound is
    worse than the incumbent.  This is exact branch-and-bound, not a local or
    greedy selection; the greedy pass only supplies an initial incumbent.
    """
    canonical_mapping = dict(canonical_mapping)
    atom_count = len(canonical_mapping)
    factors = list(_independent_atom_action_factors(
        raw_generators, atom_count))
    group_order = 1
    for _support, actions in factors:
        group_order *= len(actions)
    if not factors:
        rmsd = fixed_mapping_aligned_rmsd(
            canonical_mapping, coords_R, coords_P)
        return canonical_mapping, rmsd, {
            'group_order': 1, 'evaluated_leaf_count': 1,
            'pruned_leaf_count': 0, 'factor_orders': [],
        }
    # Small groups are faster as one bounded vectorized batch.  The threshold
    # is explicit, so this path can never recreate an unbounded group closure.
    if group_order <= 4096:
        actions = _generated_atom_permutations(
            raw_generators, atom_count)
        candidates = [{
            r: int(action[canonical_mapping[r]])
            for r in canonical_mapping
        } for action in actions]
        rmsds = _fixed_mappings_aligned_rmsd(
            candidates, coords_R, coords_P)
        ranked = [(
            round(float(rmsd), 12),
            tuple(candidate[r] for r in range(atom_count)),
            float(rmsd), candidate,
        ) for candidate, rmsd in zip(candidates, rmsds)]
        _rounded, _key, rmsd, selected = min(
            ranked, key=lambda item: item[:2])
        return selected, rmsd, {
            'group_order': int(group_order),
            'evaluated_leaf_count': int(group_order),
            'pruned_leaf_count': 0,
            'factor_orders': [len(actions) for _support, actions in factors],
        }

    # Larger factors and factors touching more atoms are decided first so the
    # invariant distance bound becomes informative as early as possible.
    factors.sort(key=lambda item: (-len(item[0]), -len(item[1]), item[0]))
    coords_R = np.asarray(coords_R, dtype=float)
    coords_P = np.asarray(coords_P, dtype=float)
    distances_R = np.linalg.norm(
        coords_R[:, None, :] - coords_R[None, :, :], axis=2)
    distances_P = np.linalg.norm(
        coords_P[:, None, :] - coords_P[None, :, :], axis=2)
    movable_targets = set().union(*(set(support) for support, _ in factors))
    fixed_R = tuple(sorted(
        r for r, p in canonical_mapping.items() if p not in movable_targets))

    def distance_increment(mapping, new_atoms, decided_atoms):
        total = 0.0
        new_atoms = tuple(new_atoms)
        for offset, left in enumerate(new_atoms):
            for right in decided_atoms:
                delta = (distances_R[left, right]
                         - distances_P[mapping[left], mapping[right]])
                total += float(delta * delta)
            for right in new_atoms[:offset]:
                delta = (distances_R[left, right]
                         - distances_P[mapping[left], mapping[right]])
                total += float(delta * delta)
        return total

    initial_sum = distance_increment(canonical_mapping, fixed_R, ())
    factor_details = []
    for support, actions in factors:
        support_set = set(support)
        affected_R = tuple(sorted(
            r for r, p in canonical_mapping.items() if p in support_set))
        action_images = tuple(
            tuple(int(action[canonical_mapping[r]]) for r in affected_R)
            for action in actions)
        factor_details.append((affected_R, action_images))

    def cross_factor_cost(left_R, left_images, right_R, right_images):
        total = 0.0
        for left, left_P in zip(left_R, left_images):
            for right, right_P in zip(right_R, right_images):
                delta = (distances_R[left, right]
                         - distances_P[left_P, right_P])
                total += float(delta * delta)
        return total

    # Pairwise minima between every two still-undecided independent factors
    # are mutually relaxed (their minimizing actions need not agree), hence
    # their sum is a rigorous lower bound.  Including it prevents the search
    # from descending through millions of leaves before unavoidable geometric
    # disagreement becomes visible.
    pairwise_floor = {}
    for left_index, (left_R, left_actions) in enumerate(factor_details):
        for right_index in range(left_index + 1, len(factor_details)):
            right_R, right_actions = factor_details[right_index]
            pairwise_floor[(left_index, right_index)] = min(
                cross_factor_cost(left_R, left_images,
                                  right_R, right_images)
                for left_images in left_actions
                for right_images in right_actions
            )

    def action_mapping(base, affected_R, images):
        candidate = dict(base)
        candidate.update(zip(affected_R, images))
        return candidate

    def optimistic_remaining(index, mapping, decided):
        floor = 0.0
        for factor_index in range(index, len(factor_details)):
            affected_R, action_images = factor_details[factor_index]
            floor += min(
                distance_increment(
                    action_mapping(mapping, affected_R, images),
                    affected_R, decided)
                for images in action_images
            )
        for left_index in range(index, len(factor_details)):
            for right_index in range(left_index + 1, len(factor_details)):
                floor += pairwise_floor[(left_index, right_index)]
        return floor

    # A deterministic greedy descent supplies a strong incumbent but never
    # removes alternatives from the exact search below.
    greedy = dict(canonical_mapping)
    for affected_R, action_images in factor_details:
        trials = []
        for images in action_images:
            candidate = action_mapping(greedy, affected_R, images)
            rmsd = fixed_mapping_aligned_rmsd(candidate, coords_R, coords_P)
            trials.append((round(rmsd, 12), tuple(candidate[r]
                                                  for r in range(atom_count)),
                           rmsd, candidate))
        greedy = min(trials, key=lambda item: item[:2])[3]
    best_mapping = greedy
    best_rmsd = fixed_mapping_aligned_rmsd(greedy, coords_R, coords_P)
    best_rank = (round(best_rmsd, 12),
                 tuple(greedy[r] for r in range(atom_count)))
    evaluated = 0
    pruned = 0
    suffix_orders = [1] * (len(factors) + 1)
    for index in range(len(factors) - 1, -1, -1):
        suffix_orders[index] = (
            suffix_orders[index + 1] * len(factors[index][1]))

    def search(index, mapping, decided, distance_sum):
        nonlocal best_mapping, best_rmsd, best_rank, evaluated, pruned
        relaxed_sum = distance_sum + optimistic_remaining(
            index, mapping, decided)
        lower_bound = np.sqrt(max(relaxed_sum, 0.0)) / atom_count
        if lower_bound > best_rmsd + 1e-12:
            pruned += suffix_orders[index]
            return
        if index == len(factors):
            evaluated += 1
            rmsd = fixed_mapping_aligned_rmsd(mapping, coords_R, coords_P)
            rank = (round(rmsd, 12),
                    tuple(mapping[r] for r in range(atom_count)))
            if rank < best_rank:
                best_mapping = dict(mapping)
                best_rmsd = float(rmsd)
                best_rank = rank
            return
        affected_R, action_images = factor_details[index]
        children = []
        for images in action_images:
            child = action_mapping(mapping, affected_R, images)
            increment = distance_increment(child, affected_R, decided)
            children.append((distance_sum + increment,
                             tuple(child[r] for r in affected_R), child))
        for child_sum, _key, child in sorted(children, key=lambda item: item[:2]):
            search(index + 1, child, decided + affected_R, child_sum)

    search(0, dict(canonical_mapping), fixed_R, initial_sum)
    return best_mapping, best_rmsd, {
        'group_order': int(group_order),
        'evaluated_leaf_count': int(evaluated),
        'pruned_leaf_count': int(pruned),
        'factor_orders': [len(actions) for _support, actions in factors],
    }


class AnalyticalMappingFamily:
    """Exact isomorphism coset represented by one completed AAM branch.

    The branch relation defines ``Iso(A, B)``.  Restricting those
    isomorphisms to atom vertices and composing with the branch's R->P source
    mapping gives every concrete bijection represented by the branch.  Group
    equality is proven by generator containment; no elements are enumerated.
    """

    def __init__(self, source_mapping, relation):
        import pynauty

        self.source_mapping = dict(source_mapping)
        self.degree = len(self.source_mapping)
        self.relation = relation
        self.graph_A = relation.graph("A")
        self.graph_B = relation.graph("B")
        isomorphism = _canonical_isomorphism(self.graph_A, self.graph_B)
        if isomorphism is None:
            raise IndexChiralityConflict(
                "analytical AAM branch relation has no endpoint isomorphism")
        self.representative_mapping = {
            r: int(isomorphism[self.source_mapping[r]])
            for r in sorted(self.source_mapping)
        }
        (raw_generators, mantissa, exponent,
         raw_orbits, _orbit_count) = pynauty.autgrp(self.graph_B)
        identity = tuple(range(self.degree))
        self.target_generators = tuple(dict.fromkeys(
            tuple(int(generator[atom]) for atom in range(self.degree))
            for generator in raw_generators
            if tuple(int(generator[atom])
                     for atom in range(self.degree)) != identity
        ))
        orbit_members = {}
        for atom in range(self.degree):
            orbit_members.setdefault(int(raw_orbits[atom]), []).append(atom)
        self.target_orbits = tuple(sorted(
            (tuple(members) for members in orbit_members.values()),
            key=lambda members: members[0]))
        self.group_order = (round(float(mantissa), 12), int(exponent))
        self._membership_cache = {}
        self._relation_record_counts_B = Counter(
            self.relation.relation_records_B)

    @property
    def invariant(self):
        return self.degree, self.group_order, self.target_orbits

    def contains(self, mapping):
        mapping = validate_mapping(
            mapping,
            ["X"] * self.degree,
            ["X"] * self.degree)
        key = tuple(mapping[r] for r in range(self.degree))
        cached = self._membership_cache.get(key)
        if cached is not None:
            return cached
        source_order = tuple(self.source_mapping[r]
                             for r in range(self.degree))
        target_order = tuple(mapping[r] for r in range(self.degree))
        sigma = {
            int(source_atom): int(target_atom)
            for source_atom, target_atom in zip(source_order, target_order)
        }
        atom_colors_match = all(
            self.relation.colors_A[atom]
            == self.relation.colors_B[sigma[atom]]
            for atom in range(self.degree)
        )

        def transported_record_counts():
            records = Counter()
            for kind, color, atoms in self.relation.relation_records_A:
                images = tuple(sigma[atom] for atom in atoms)
                if kind == "pair":
                    images = tuple(sorted(images))
                record = (kind, color, images)
                if record not in self._relation_record_counts_B:
                    return None
                records[record] += 1
            return records

        transported = transported_record_counts() if atom_colors_match else None
        result = (transported is not None
                  and transported == self._relation_record_counts_B)
        self._membership_cache[key] = bool(result)
        return bool(result)

    @staticmethod
    def _left_act(generator, mapping):
        return {
            r: int(generator[mapping[r]])
            for r in mapping
        }

    def equivalent(self, other):
        if not isinstance(other, AnalyticalMappingFamily):
            return False
        return self.is_subset_of(other) and other.is_subset_of(self)

    def is_subset_of(self, other):
        """Prove coset inclusion by representative/generator membership."""
        if not isinstance(other, AnalyticalMappingFamily):
            return False
        if self.degree != other.degree:
            return False
        self_log_order = (np.log10(self.group_order[0])
                          + self.group_order[1])
        other_log_order = (np.log10(other.group_order[0])
                           + other.group_order[1])
        if self_log_order > other_log_order + 1e-10:
            return False
        other_orbit = {
            atom: orbit_index
            for orbit_index, orbit in enumerate(other.target_orbits)
            for atom in orbit
        }
        # A subgroup orbit cannot cross two orbits of the containing group.
        # This is a necessary exact group-inclusion condition and rejects most
        # unrelated family pairs before any dense relation transport.
        if any(len({other_orbit[atom] for atom in orbit}) > 1
               for orbit in self.target_orbits):
            return False
        if not other.contains(self.representative_mapping):
            return False
        return not any(not other.contains(self._left_act(
            generator, self.representative_mapping))
            for generator in self.target_generators)

    def record(self):
        return {
            "representative_mapping": dict(self.representative_mapping),
            "target_generators": [list(generator)
                                  for generator in self.target_generators],
            "target_orbits": [list(orbit) for orbit in self.target_orbits],
            "group_order": {
                "mantissa": self.group_order[0],
                "decimal_exponent": self.group_order[1],
            },
        }


def compile_analytical_mapping_family(
        source_mapping, branch_symmetry,
        elements_R: Sequence[str], wbo_R,
        elements_P: Sequence[str], wbo_P, *,
        graph_floor=0.2, symmetry_wbo_tol=0.2,
        dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
        anchor_map=None, static_context=None):
    """Compile one completed branch into its exact pre-chirality coset."""
    source = validate_mapping(source_mapping, elements_R, elements_P)
    relation, _persistent, _inverse, _fragments = _masked_relation_data(
        source, branch_symmetry, elements_R, np.asarray(wbo_R),
        elements_P, np.asarray(wbo_P), graph_floor,
        symmetry_wbo_tol, dwbo_threshold,
        metal_dwbo_threshold, anchor_map or {},
        static_context=static_context)
    return AnalyticalMappingFamily(source, relation)


def _simplex_measure(coords, center, neighbors, degeneracy_tol):
    """Affine orientation of one local substituent simplex."""
    neighbors = tuple(int(atom) for atom in neighbors)
    if len(neighbors) == 3:
        origin, points = int(center), neighbors
    elif len(neighbors) == 4:
        origin, points = neighbors[0], neighbors[1:]
    else:
        raise IndexChiralityError(
            "an oriented substituent simplex needs three or four neighbors")
    return _orientation_measure(
        coords, origin, points, degeneracy_tol=degeneracy_tol)


def select_group_chiral_witness(
        source_mapping, witnesses,
        elements_R: Sequence[str], coords_R, wbo_R,
        elements_P: Sequence[str], coords_P, wbo_P, *,
        graph_floor=0.2,
        orientation_degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL):
    """Choose the same-event AAM witness preserving most group orientation.

    Higher-coordinate chirality is a relationship between separate AAM
    isomorphisms, so it must be resolved before one witness is collapsed into
    its own automorphism group.  For every R center that is higher-coordinate
    at either endpoint under any witness, the fixed reference set is all
    defined center-relative triples of its R neighbors.  Missing, planar, and
    reversed P relations are distinguished; the witness preserving the most
    defined relations wins deterministically.
    """
    source = validate_mapping(source_mapping, elements_R, elements_P)
    raw_witnesses = list(witnesses or ())
    candidate_records = []
    seen = set()
    for witness_index, witness in enumerate(raw_witnesses):
        mapping = validate_mapping(
            witness.get("mapping") or {}, elements_R, elements_P)
        key = tuple(mapping[r] for r in sorted(mapping))
        if key in seen:
            continue
        seen.add(key)
        candidate_records.append((witness_index, mapping))
    source_key = tuple(source[r] for r in sorted(source))
    if source_key not in seen:
        candidate_records.insert(0, (None, source))
    if not candidate_records:
        candidate_records = [(None, source)]

    graph_R = build_graph(elements_R, wbo_R, bond_cut=graph_floor)
    graph_P = build_graph(elements_P, wbo_P, bond_cut=graph_floor)
    high_coordinate_centers = {
        int(center_R)
        for center_R in graph_R.nodes()
        if graph_R.degree(center_R) > 4
        or any(graph_P.degree(mapping[center_R]) > 4
               for _index, mapping in candidate_records)
    }
    reference_frames = []
    for center_R in sorted(high_coordinate_centers):
        for neighbors_R in combinations(
                sorted(graph_R.neighbors(center_R)), 3):
            measure_R = _orientation_measure(
                coords_R, center_R, neighbors_R,
                degeneracy_tol=orientation_degeneracy_tol)
            if measure_R.sign:
                reference_frames.append((
                    center_R, tuple(neighbors_R), measure_R))

    evaluated = []
    for witness_index, mapping in candidate_records:
        preserved, reversed_frames = [], []
        missing_count = 0
        degenerate_count = 0
        for center_R, neighbors_R, measure_R in reference_frames:
            center_P = mapping[center_R]
            neighbors_P = tuple(mapping[r] for r in neighbors_R)
            if not all(graph_P.has_edge(center_P, p) for p in neighbors_P):
                missing_count += 1
                continue
            measure_P = _orientation_measure(
                coords_P, center_P, neighbors_P,
                degeneracy_tol=orientation_degeneracy_tol)
            record = {
                "center_R": int(center_R),
                "neighbors_R_index_order": list(neighbors_R),
                "reactant_orientation_sign": int(measure_R.sign),
                "product_orientation_sign": int(measure_P.sign),
                "reactant_normalized_orientation": measure_R.normalized,
                "product_normalized_orientation": measure_P.normalized,
            }
            if not measure_P.sign:
                degenerate_count += 1
            elif measure_P.sign == measure_R.sign:
                preserved.append(record)
            else:
                reversed_frames.append(record)
        mapping_key = tuple(mapping[r] for r in sorted(mapping))
        fixed_rmsd = fixed_mapping_aligned_rmsd(
            mapping, coords_R, coords_P)
        rank = (
            -len(preserved),
            len(reversed_frames),
            degenerate_count,
            missing_count,
            fixed_rmsd,
            mapping_key,
        )
        evaluated.append((
            rank, witness_index, mapping, preserved, reversed_frames,
            degenerate_count, missing_count, fixed_rmsd))
    evaluated.sort(key=lambda item: item[0])
    (_rank, witness_index, selected, preserved, reversed_frames,
     degenerate_count, missing_count, selected_rmsd) = evaluated[0]
    metadata = {
        "schema_version": "rxn_core.group_chirality_witness/v1",
        "policy": (
            "maximize_preserved_orientation_then_minimize_fixed_mapping_rmsd"),
        "candidate_witness_count": len(candidate_records),
        "high_coordinate_centers_R": sorted(high_coordinate_centers),
        "reference_frame_count": len(reference_frames),
        "preserved_frame_count": len(preserved),
        "reversed_frame_count": len(reversed_frames),
        "degenerate_frame_count": int(degenerate_count),
        "missing_frame_count": int(missing_count),
        "selected_witness_index": witness_index,
        "selected_mapping_changed": selected != source,
        "selected_fixed_mapping_aligned_rmsd": float(selected_rmsd),
        "rmsd_policy": (
            "exact_mapping_then_proper_rigid_fit_no_permutation"),
        "selected_reversed_frames": reversed_frames,
    }
    return GroupChiralityWitnessSelection(
        selected_mapping=selected,
        selected_witness_index=witness_index,
        preserved_frames=tuple(preserved),
        metadata=metadata,
    )


def analyze_group_chirality_branch(
        mapping,
        elements_R: Sequence[str], coords_R, wbo_R,
        elements_P: Sequence[str], coords_P, wbo_P, *,
        graph_floor=0.2,
        orientation_degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL):
    """Analyze all defined persistent coordination frames in one branch.

    This is the post-AAM branch API.  It has exactly one immutable branch
    representative and performs no candidate ranking or witness sampling.
    Both currently preserved and reversed frames are returned as constraints;
    the relational solver decides whether the branch family can satisfy them.
    """
    legacy = select_group_chiral_witness(
        mapping, (), elements_R, coords_R, wbo_R,
        elements_P, coords_P, wbo_P,
        graph_floor=graph_floor,
        orientation_degeneracy_tol=orientation_degeneracy_tol)
    metadata = dict(legacy.metadata)
    reversed_frames = list(metadata.pop('selected_reversed_frames', ()))
    metadata.pop('selected_witness_index', None)
    metadata.pop('candidate_witness_count', None)
    metadata['schema_version'] = "rxn_core.group_chirality_branch/v1"
    metadata['policy'] = "analytical_branch_constraint_construction"
    frames = tuple(list(legacy.preserved_frames) + reversed_frames)
    metadata['defined_constraint_frame_count'] = len(frames)
    return GroupChiralityBranchAnalysis(
        mapping=dict(legacy.selected_mapping),
        defined_frames=frames,
        metadata=metadata)


def _frame_measure_for_mapping(frame, mapping, coords_P, degeneracy_tol):
    return _simplex_measure(
        coords_P, mapping[frame.center_R],
        tuple(mapping[r] for r in frame.neighbors_R), degeneracy_tol)


def select_index_chirality_assignment(
        source_mapping, branch_symmetry,
        elements_R: Sequence[str], coords_R, wbo_R,
        elements_P: Sequence[str], coords_P, wbo_P, *,
        graph_floor=0.2, symmetry_wbo_tol=0.2,
        dwbo_threshold=0.5, metal_dwbo_threshold=0.3,
        orientation_degeneracy_tol=ORIENTATION_DEGENERACY_TOL,
        anchor_map=None, group_chirality_frames=(), static_context=None):
    """Compose one AAM witness with one exact orientation-preserving automorph.

    The constraint is solved simultaneously as a colored relational-graph
    isomorphism.  There is no degree-specific rule, greedy generator sequence,
    capped orbit enumeration, or alternate matching path.
    """
    source = validate_mapping(source_mapping, elements_R, elements_P)
    relation, persistent_P, inverse, fragments = _masked_relation_data(
        source, branch_symmetry, elements_R, np.asarray(wbo_R),
        elements_P, np.asarray(wbo_P), graph_floor,
        symmetry_wbo_tol, dwbo_threshold, metal_dwbo_threshold,
        anchor_map, static_context=static_context)

    group_frame_candidates = []
    for raw_frame in group_chirality_frames or ():
        center_R = int(raw_frame["center_R"])
        neighbors_R = tuple(
            int(r) for r in raw_frame["neighbors_R_index_order"])
        if len(neighbors_R) != 3:
            raise IndexChiralityError(
                "a group-chirality frame requires exactly three ligands")
        center_P = source[center_R]
        neighbors_P = tuple(source[r] for r in neighbors_R)
        if not all(
                float(wbo_R[center_R, r]) >= float(graph_floor)
                and float(wbo_P[center_P, p]) >= float(graph_floor)
                for r, p in zip(neighbors_R, neighbors_P)):
            raise IndexChiralityError(
                "selected group-chirality frame is not persistent")
        measure_R = _orientation_measure(
            coords_R, center_R, neighbors_R,
            degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL)
        measure_P = _orientation_measure(
            coords_P, center_P, neighbors_P,
            degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL)
        if not measure_R.sign or not measure_P.sign:
            raise IndexChiralityError(
                "group-chirality frame has undefined endpoint orientation")
        group_frame_candidates.append((
            min(abs(measure_R.normalized), abs(measure_P.normalized)),
            center_R, neighbors_R, measure_R.sign, center_P,
        ))

    # Higher-coordinate geometries can continuously reconfigure so that one
    # ligand triple crosses coplanarity while the overall ligand assignment
    # remains orientation-consistent.  Build a maximal feasible signed-frame
    # basis, strongest simplices first, rather than making one conflicting
    # triple reject the entire exact mapping family.
    group_frames = []
    reconfigured_group_frames = []
    for robustness, center_R, neighbors_R, sign_R, center_P in sorted(
            group_frame_candidates,
            key=lambda item: (-item[0], item[1], item[2])):
        trial = relation.clone()
        for ordered_R in permutations(neighbors_R):
            ordered_P = tuple(source[r] for r in ordered_R)
            ordered_sign_R = _orientation_measure(
                coords_R, center_R, ordered_R,
                degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL).sign
            sign_P = _orientation_measure(
                coords_P, center_P, ordered_P,
                degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL).sign
            trial.add_ordered_relation(
                (center_P, *ordered_P),
                ("oriented_coordination_triple", ordered_sign_R),
                ("oriented_coordination_triple", sign_P),
            )
        if _canonical_isomorphism(trial.graph("A"), trial.graph("B")) is None:
            reconfigured_group_frames.append({
                "center_R": int(center_R),
                "neighbors_R_index_order": list(neighbors_R),
                "reactant_orientation_sign": int(sign_R),
                "normalized_orientation_robustness": float(robustness),
                "reason": "incompatible_with_stronger_group_frame_basis",
            })
            continue
        relation = trial
        group_frames.append((center_R, neighbors_R, sign_R))

    base_A, base_B = relation.graph("A"), relation.graph("B")
    if _canonical_isomorphism(base_A, base_B) is None:
        raise IndexChiralityConflict(
            "selected AAM fragments admit no automorphism satisfying the "
            "event and anchor constraints")

    import pynauty
    (_generators, group_mantissa, group_exponent,
     raw_orbits, _orbit_count) = pynauty.autgrp(base_B)
    atom_orbit_sizes = Counter(raw_orbits[p] for p in range(len(source)))
    movable_P = {
        p for p in range(len(source))
        if atom_orbit_sizes[raw_orbits[p]] > 1
    }
    switchable_R = sorted(r for r, p in source.items() if p in movable_P)

    locally_mutable_centers = set()
    for center_P, neighbors_P in persistent_P.items():
        if len(neighbors_P) < 3:
            continue
        center_fixed_B = relation.graph("B", individualized=(center_P,))
        (_center_generators, _center_mantissa, _center_exponent,
         center_orbits, _center_orbit_count) = pynauty.autgrp(center_fixed_B)
        neighbor_orbit_sizes = Counter(
            center_orbits[p] for p in neighbors_P)
        if any(neighbor_orbit_sizes[center_orbits[p]] > 1
               for p in neighbors_P):
            locally_mutable_centers.add(center_P)

    frame_records = {}
    for center_P in range(len(source)):
        if center_P not in locally_mutable_centers:
            continue
        neighbor_count = len(persistent_P[center_P])
        simplex_size = 3 if neighbor_count == 3 else 4
        for neighbors_P in combinations(
                sorted(persistent_P[center_P]), simplex_size):
            center_R = inverse[center_P]
            neighbors_R = tuple(inverse[p] for p in neighbors_P)
            measure_R = _simplex_measure(
                coords_R, center_R, neighbors_R,
                orientation_degeneracy_tol)
            measure_P = _simplex_measure(
                coords_P, center_P, neighbors_P,
                orientation_degeneracy_tol)
            frame = _OrientedFrame(
                center_R=center_R,
                neighbors_R=neighbors_R,
                sign_R=measure_R.sign,
                sign_P_source=measure_P.sign,
                normalized_R=measure_R.normalized,
                normalized_P_source=measure_P.normalized,
            )
            frame_records[(center_P, neighbors_P)] = (
                frame, measure_R, measure_P)

    # A coplanar simplex orbit carries no handedness, but does not erase the
    # handedness of other simplices at a high-coordinate center.  First form
    # exact simplex orbits, remove only zero-containing orbits, then test the
    # remaining simplices jointly for each symmetry-equivalent center orbit.
    parent = {key: key for key in frame_records}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for generator in _generators:
        for center_P, neighbors_P in frame_records:
            image = (
                int(generator[center_P]),
                tuple(sorted(int(generator[p]) for p in neighbors_P)),
            )
            if image in frame_records:
                union((center_P, neighbors_P), image)

    simplex_orbit_keys = defaultdict(list)
    for key in frame_records:
        simplex_orbit_keys[find(key)].append(key)
    simplex_orbit_is_stereogenic = {}
    for root, keys in simplex_orbit_keys.items():
        simplex_orbit_is_stereogenic[root] = not any(
            frame_records[key][1].sign == 0
            or frame_records[key][2].sign == 0
            for key in keys)

    center_orbit_keys = defaultdict(list)
    for simplex_root, keys in simplex_orbit_keys.items():
        if not simplex_orbit_is_stereogenic[simplex_root]:
            continue
        for key in keys:
            center_orbit_keys[int(raw_orbits[key[0]])].append(key)

    def add_frame_relations(target, keys):
        for center_P, neighbors_P in keys:
            center_R = inverse[center_P]
            for ordered_P in permutations(neighbors_P):
                ordered_R = tuple(inverse[p] for p in ordered_P)
                sign_R = _simplex_measure(
                    coords_R, center_R, ordered_R,
                    orientation_degeneracy_tol).sign
                sign_P = _simplex_measure(
                    coords_P, center_P, ordered_P,
                    orientation_degeneracy_tol).sign
                target.add_ordered_relation(
                    (center_P, *ordered_P),
                    ("oriented_substituent_simplex", len(ordered_P), sign_R),
                    ("oriented_substituent_simplex", len(ordered_P), sign_P),
                )

    nonstereogenic_frames = []
    reconfigured_frames = []
    for root, keys in simplex_orbit_keys.items():
        if not simplex_orbit_is_stereogenic[root]:
            for key in keys:
                frame, measure_R, measure_P = frame_records[key]
                nonstereogenic_frames.append({
                    "id": frame.frame_id,
                    "center_R": frame.center_R,
                    "neighbors_R_index_order": list(frame.neighbors_R),
                    "reactant_orientation_sign": measure_R.sign,
                    "source_product_orientation_sign": measure_P.sign,
                    "reason": (
                        "automorphism_orbit_contains_nonstereogenic_endpoint"),
                })

    eligible_center_roots = set()
    for root, keys in center_orbit_keys.items():
        local_relation = relation.clone()
        add_frame_relations(local_relation, keys)
        if _canonical_isomorphism(
                local_relation.graph("A"), local_relation.graph("B")) is None:
            for key in keys:
                frame, measure_R, measure_P = frame_records[key]
                reconfigured_frames.append({
                    "id": frame.frame_id,
                    "center_R": frame.center_R,
                    "neighbors_R_index_order": list(frame.neighbors_R),
                    "reactant_orientation_sign": measure_R.sign,
                    "source_product_orientation_sign": measure_P.sign,
                    "reason": (
                        "local_coordination_geometry_is_not_"
                        "stereochemically_equivalent"),
                })
            continue
        eligible_center_roots.add(root)

    active_frames = [
        frame_records[key][0]
        for root in eligible_center_roots for key in center_orbit_keys[root]
    ]
    oriented_relation = relation.clone()
    for root in eligible_center_roots:
        add_frame_relations(oriented_relation, center_orbit_keys[root])

    oriented_A = oriented_relation.graph("A")
    oriented_B = oriented_relation.graph("B")
    isomorphism = _canonical_isomorphism(oriented_A, oriented_B)
    if isomorphism is None:
        source_mismatch_frames = [
            frame for frame in active_frames
            if frame.sign_R != frame.sign_P_source
        ]
        raise IndexChiralityConflict(
            "no exact selected-fragment automorphism satisfies all signed "
            f"substituent simplices simultaneously; source mismatches="
            f"{len(source_mismatch_frames)}",
            diagnostics={
                "constraint_model": (
                    "simultaneous_affine_substituent_simplices"),
                "orientation_degeneracy_tolerance": float(
                    orientation_degeneracy_tol),
                "switchable_r_atoms": switchable_R,
                "active_frame_count": len(active_frames),
                "source_mismatch_frames": [{
                    "id": frame.frame_id,
                    "center_R": frame.center_R,
                    "neighbors_R_index_order": list(frame.neighbors_R),
                    "reactant_orientation_sign": frame.sign_R,
                    "source_product_orientation_sign": frame.sign_P_source,
                    "reactant_normalized_orientation": frame.normalized_R,
                    "source_product_normalized_orientation": (
                        frame.normalized_P_source),
                } for frame in source_mismatch_frames],
            })
    import pynauty
    (oriented_generators, oriented_mantissa, oriented_exponent,
     _oriented_orbits, _oriented_orbit_count) = pynauty.autgrp(oriented_B)
    canonical_mapping = {
        r: int(isomorphism[source[r]]) for r in source
    }
    selected, selected_rmsd, rmsd_search = _minimum_rmsd_group_action(
        canonical_mapping, oriented_generators, coords_R, coords_P)
    validate_mapping(selected, elements_R, elements_P)

    source_signature = mapping_event_signature(
        source, wbo_R, wbo_P, elements_R,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold)
    selected_signature = mapping_event_signature(
        selected, wbo_R, wbo_P, elements_R,
        dwbo_threshold=dwbo_threshold,
        metal_dwbo_threshold=metal_dwbo_threshold)
    if selected_signature != source_signature:
        raise IndexChiralityConflict(
            "internal invariant failure: constrained automorphism changed the "
            "selected AAM bond event")

    source_violations = [
        frame for frame in active_frames
        if frame.sign_R != frame.sign_P_source
    ]
    selected_violations = []
    for frame in active_frames:
        if _frame_measure_for_mapping(
                frame, selected, coords_P,
                orientation_degeneracy_tol).sign != frame.sign_R:
            selected_violations.append(frame)
    if selected_violations:
        raise IndexChiralityConflict(
            "internal invariant failure: relational isomorphism did not "
            "preserve every affine substituent simplex")
    selected_group_violations = []
    for center_R, neighbors_R, sign_R in group_frames:
        sign_P = _orientation_measure(
            coords_P, selected[center_R],
            tuple(selected[r] for r in neighbors_R),
            degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL).sign
        if sign_P != sign_R:
            selected_group_violations.append((center_R, neighbors_R))
    if selected_group_violations:
        raise IndexChiralityConflict(
            "internal invariant failure: relational isomorphism changed "
            "preserved group-level coordination orientation")
    metadata = {
        "schema_version": INDEX_CHIRALITY_SCHEMA,
        "policy": "preserve",
        "status": "applied" if selected != source else "already_consistent",
        "solver": "pynauty_colored_relational_isomorphism",
        "constraint_model": "simultaneous_affine_substituent_simplices",
        "orientation_degeneracy_tolerance": float(
            orientation_degeneracy_tol),
        "candidate_source": "selected_AAM_masked_fragment_automorphism",
        "allowed_automorphism_group_order": {
            "mantissa": float(group_mantissa),
            "decimal_exponent": int(group_exponent),
        },
        "chirality_valid_atom_bijection_count": rmsd_search['group_order'],
        "chirality_relation_automorphism_group_order": {
            "mantissa": float(oriented_mantissa),
            "decimal_exponent": int(oriented_exponent),
        },
        "rmsd_candidate_count": rmsd_search['group_order'],
        "rmsd_evaluated_leaf_count": rmsd_search['evaluated_leaf_count'],
        "rmsd_pruned_leaf_count": rmsd_search['pruned_leaf_count'],
        "rmsd_symmetry_factor_orders": rmsd_search['factor_orders'],
        "selected_fixed_mapping_aligned_rmsd": float(selected_rmsd),
        "rmsd_policy": (
            "exact_symmetry_factor_branch_and_bound_then_fixed_mapping_"
            "proper_fit_no_remapping"),
        "selected_fragment_count": len(fragments),
        "preserved_group_chirality_frame_count": len(group_frames),
        "reconfigured_group_chirality_frame_count": len(
            reconfigured_group_frames),
        "switchable_r_atoms": switchable_R,
        "defined_frame_count": sum(
            frame.sign_R != 0 and frame.sign_P_source != 0
            for frame in active_frames),
        "nonstereogenic_frame_count": len(nonstereogenic_frames),
        "reconfigured_frame_count": len(reconfigured_frames),
        "source_index_chirality_violation_count": len(source_violations),
        "selected_index_chirality_violation_count": 0,
        "mapping_changes": [{
            "r_atom": r,
            "source_p_atom": source[r],
            "selected_p_atom": selected[r],
        } for r in sorted(source) if source[r] != selected[r]],
        "active_frames": [{
            "id": frame.frame_id,
            "center_R": frame.center_R,
            "neighbors_R_index_order": list(frame.neighbors_R),
            "reactant_orientation_sign": frame.sign_R,
            "source_product_orientation_sign": frame.sign_P_source,
            "reactant_normalized_orientation": frame.normalized_R,
            "source_product_normalized_orientation": (
                frame.normalized_P_source),
        } for frame in active_frames],
        "preserved_group_chirality_frames": [{
            "center_R": center_R,
            "neighbors_R_index_order": list(neighbors_R),
            "reactant_orientation_sign": sign_R,
        } for center_R, neighbors_R, sign_R in group_frames],
        "reconfigured_group_chirality_frames": reconfigured_group_frames,
        "nonstereogenic_frames": nonstereogenic_frames,
        "reconfigured_frames": reconfigured_frames,
        "event_signature_unchanged": True,
    }
    return IndexChiralitySelection(
        source_mapping=source,
        selected_mapping=selected,
        metadata=metadata,
    )

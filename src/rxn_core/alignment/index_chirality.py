"""Exact index-orientation consensus for one selected AAM witness.

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
        relation = self._vertex(color_A, color_B)
        self._edge(relation, int(left))
        self._edge(relation, int(right))

    def add_ordered_relation(self, atoms, color_A, color_B):
        relation = self._vertex(color_A, color_B)
        for role, atom in enumerate(atoms):
            role_vertex = self._vertex(("orientation_role", role))
            self._edge(relation, role_vertex)
            self._edge(role_vertex, int(atom))

    def clone(self):
        result = _RelationalGraph([], [])
        result.colors_A = list(self.colors_A)
        result.colors_B = list(self.colors_B)
        result.adjacency = [set(neighbors) for neighbors in self.adjacency]
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


def _selected_fragments(branch_symmetry, source):
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
    if missing:
        raise IndexChiralityError(
            f"selected AAM fragments do not cover R atoms {missing}")
    return tuple(fragments), owner


def _masked_relation_data(source, branch_symmetry, elements_R, wbo_R,
                          elements_P, wbo_P, graph_floor,
                          symmetry_wbo_tol, dwbo_threshold,
                          metal_dwbo_threshold, anchor_map):
    atom_count = len(source)
    inverse = {p: r for r, p in source.items()}
    fragments, owner_R = _selected_fragments(branch_symmetry, source)
    owner_P = {source[r]: owner_R[r] for r in source}
    anchors = {int(r): int(p) for r, p in dict(anchor_map or {}).items()}
    for r, p in anchors.items():
        if r not in source or p not in inverse:
            raise IndexChiralityError(f"invalid anchor R{r}->P{p}")
        if str(elements_R[r]) != str(elements_P[p]):
            raise IndexChiralityError(f"anchor element mismatch R{r}->P{p}")

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

    g_R = build_graph(elements_R, wbo_R, bond_cut=graph_floor)
    g_P = build_graph(elements_P, wbo_P, bond_cut=graph_floor)
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
    pair_groups_R = defaultdict(list)
    pair_groups_P = defaultdict(list)
    pair_records = []
    for left_P in range(atom_count):
        for right_P in range(left_P + 1, atom_count):
            left_R, right_R = inverse[left_P], inverse[right_P]
            element_pair = tuple(sorted((
                str(elements_R[left_R]), str(elements_R[right_R]))))
            threshold = bond_event_threshold(
                elements_R, left_R, right_R,
                default_threshold=float(dwbo_threshold),
                metal_threshold=metal_dwbo_threshold)
            r_value = float(wbo_R[left_R, right_R])
            p_value = float(wbo_P[left_P, right_P])
            group_key = (element_pair, float(threshold))
            pair_groups_R[group_key].append(r_value)
            pair_groups_P[group_key].append(p_value)
            pair_records.append((
                left_P, right_P, group_key, r_value, p_value))
    values_R = {
        key: tuple(sorted(set(values))) for key, values in pair_groups_R.items()
    }
    values_P = {
        key: tuple(sorted(set(values))) for key, values in pair_groups_P.items()
    }
    for left_P, right_P, group_key, r_value, p_value in pair_records:
        threshold = group_key[1]
        r_behavior = tuple(
            _event_class(r_value, other, threshold)
            for other in values_P[group_key])
        p_behavior = tuple(
            _event_class(other, p_value, threshold)
            for other in values_R[group_key])
        relation.add_pair(
            left_P, right_P,
            ("event_invariant_pair", group_key, r_behavior, p_behavior))

    return relation, persistent_P, inverse, fragments


def _canonical_isomorphism(graph_A, graph_B):
    import pynauty

    if pynauty.certificate(graph_A) != pynauty.certificate(graph_B):
        return None
    label_A = pynauty.canon_label(graph_A)
    label_B = pynauty.canon_label(graph_B)
    return {int(label_A[i]): int(label_B[i])
            for i in range(len(label_A))}


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
        rank = (
            -len(preserved),
            len(reversed_frames),
            degenerate_count,
            missing_count,
            mapping_key != source_key,
            mapping_key,
        )
        evaluated.append((
            rank, witness_index, mapping, preserved, reversed_frames,
            degenerate_count, missing_count))
    evaluated.sort(key=lambda item: item[0])
    (_rank, witness_index, selected, preserved, reversed_frames,
     degenerate_count, missing_count) = evaluated[0]
    metadata = {
        "schema_version": "rxn_core.group_chirality_witness/v1",
        "policy": "maximize_preserved_high_coordinate_orientation",
        "candidate_witness_count": len(candidate_records),
        "high_coordinate_centers_R": sorted(high_coordinate_centers),
        "reference_frame_count": len(reference_frames),
        "preserved_frame_count": len(preserved),
        "reversed_frame_count": len(reversed_frames),
        "degenerate_frame_count": int(degenerate_count),
        "missing_frame_count": int(missing_count),
        "selected_witness_index": witness_index,
        "selected_mapping_changed": selected != source,
        "selected_reversed_frames": reversed_frames,
    }
    return GroupChiralityWitnessSelection(
        selected_mapping=selected,
        selected_witness_index=witness_index,
        preserved_frames=tuple(preserved),
        metadata=metadata,
    )


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
        anchor_map=None, group_chirality_frames=()):
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
        anchor_map)

    group_frames = []
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
        if not measure_R.sign or measure_R.sign != measure_P.sign:
            raise IndexChiralityError(
                "selected group-chirality frame is not orientation-preserving")
        group_frames.append((center_R, neighbors_R, measure_R.sign))
        for ordered_R in permutations(neighbors_R):
            ordered_P = tuple(source[r] for r in ordered_R)
            sign_R = _orientation_measure(
                coords_R, center_R, ordered_R,
                degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL).sign
            sign_P = _orientation_measure(
                coords_P, center_P, ordered_P,
                degeneracy_tol=GROUP_ORIENTATION_DEGENERACY_TOL).sign
            relation.add_ordered_relation(
                (center_P, *ordered_P),
                ("oriented_coordination_triple", sign_R),
                ("oriented_coordination_triple", sign_P),
            )

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
    selected = {r: int(isomorphism[source[r]]) for r in source}
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
        "selected_fragment_count": len(fragments),
        "preserved_group_chirality_frame_count": len(group_frames),
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
        "nonstereogenic_frames": nonstereogenic_frames,
        "reconfigured_frames": reconfigured_frames,
        "event_signature_unchanged": True,
    }
    return IndexChiralitySelection(
        source_mapping=source,
        selected_mapping=selected,
        metadata=metadata,
    )

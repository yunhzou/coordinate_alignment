"""Direct mapped bond-angle-dihedral interpolation for viewer diagnostics."""
from __future__ import annotations

from collections import deque
from itertools import combinations

import numpy as np


_COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66,
    "F": 0.57, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20,
    "I": 1.39, "Pd": 1.39, "Pt": 1.36, "Ni": 1.24, "Cu": 1.32,
    "Fe": 1.32, "Co": 1.26, "Rh": 1.42, "Ir": 1.41,
}


def proper_align_coordinates(coords, reference):
    coords = np.asarray(coords, dtype=float)
    reference = np.asarray(reference, dtype=float)
    x_mean, y_mean = coords.mean(axis=0), reference.mean(axis=0)
    x, y = coords - x_mean, reference - y_mean
    u, _s, vt = np.linalg.svd(x.T @ y)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return x @ rotation + y_mean


def _wrap_angle(value):
    return (float(value) + np.pi) % (2.0 * np.pi) - np.pi


def _unit(vector):
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if norm < 1e-12:
        return None
    return vector / norm


def _rotation_between(left, right):
    """Shortest proper rotation taking one unit vector to another."""
    left, right = _unit(left), _unit(right)
    if left is None or right is None:
        return np.eye(3)
    cross = np.cross(left, right)
    sine = float(np.linalg.norm(cross))
    cosine = float(np.clip(np.dot(left, right), -1.0, 1.0))
    if sine < 1e-12:
        if cosine > 0:
            return np.eye(3)
        axis = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(axis, left))) > 0.8:
            axis = np.array([0.0, 1.0, 0.0])
        axis = _unit(np.cross(left, axis))
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    axis = cross / sine
    skew = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return np.eye(3) + sine * skew + (1.0 - cosine) * (skew @ skew)


def _slerp_vector(left, right, t):
    left, right = _unit(left), _unit(right)
    if left is None or right is None:
        return left if right is None else right
    dot = float(np.clip(np.dot(left, right), -1.0, 1.0))
    if dot > 1.0 - 1e-10:
        return _unit((1.0 - t) * left + t * right)
    if dot < -1.0 + 1e-10:
        axis = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(axis, left))) > 0.8:
            axis = np.array([0.0, 1.0, 0.0])
        middle = _unit(np.cross(left, axis))
        return (_slerp_vector(left, middle, 2.0 * t) if t <= 0.5
                else _slerp_vector(middle, right, 2.0 * t - 1.0))
    angle = float(np.arccos(dot))
    return (
        np.sin((1.0 - t) * angle) * left
        + np.sin(t * angle) * right
    ) / np.sin(angle)


def _anchor_frame(coords, center, first, second):
    e1 = _unit(coords[first] - coords[center])
    if e1 is None:
        return None
    second_vector = coords[second] - coords[center]
    e2 = _unit(second_vector - e1 * np.dot(second_vector, e1))
    if e2 is None:
        return None
    e3 = _unit(np.cross(e1, e2))
    return np.column_stack((e1, e2, e3))


def _interpolate_frame(frame_R, frame_P, t):
    relative = frame_P @ frame_R.T
    cosine = float(np.clip((np.trace(relative) - 1.0) / 2.0, -1.0, 1.0))
    angle = float(np.arccos(cosine))
    if angle < 1e-12:
        return frame_R.copy()
    axis = _unit(np.array([
        relative[2, 1] - relative[1, 2],
        relative[0, 2] - relative[2, 0],
        relative[1, 0] - relative[0, 1],
    ]))
    if axis is None:
        values, vectors = np.linalg.eig(relative)
        axis = _unit(np.real(vectors[:, np.argmin(abs(values - 1.0))]))
    theta = t * angle
    skew = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    rotation = np.eye(3) + np.sin(theta) * skew + (
        1.0 - np.cos(theta)) * (skew @ skew)
    return rotation @ frame_R


def _components_without_center(n_atoms, bonded_pairs, center):
    adjacency = [set() for _ in range(n_atoms)]
    for left, right in bonded_pairs:
        left, right = int(left), int(right)
        if center in (left, right):
            continue
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(range(n_atoms)) - {int(center)}
    components = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        component = {root}
        queue = deque([root])
        while queue:
            atom = queue.popleft()
            for neighbor in adjacency[atom]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.add(neighbor)
                    queue.append(neighbor)
        components.append(component)
    return components


def _coordination_constraints(coords_R, coords_P, bonded_pairs,
                              persistent_bonded_pairs):
    n_atoms = len(coords_R)
    union_neighbors = [set() for _ in range(n_atoms)]
    persistent_neighbors = [set() for _ in range(n_atoms)]
    for left, right in bonded_pairs:
        left, right = int(left), int(right)
        union_neighbors[left].add(right)
        union_neighbors[right].add(left)
    for left, right in persistent_bonded_pairs:
        left, right = int(left), int(right)
        persistent_neighbors[left].add(right)
        persistent_neighbors[right].add(left)
    constraints = []
    for center in range(n_atoms):
        neighbors = sorted(persistent_neighbors[center])
        if len(union_neighbors[center]) <= 4 or len(neighbors) < 3:
            continue
        components = _components_without_center(
            n_atoms, bonded_pairs, center)
        owner = {
            atom: component_index
            for component_index, component in enumerate(components)
            for atom in component
        }
        anchor_candidates = []
        for first, second in combinations(neighbors, 2):
            frame_R = _anchor_frame(coords_R, center, first, second)
            frame_P = _anchor_frame(coords_P, center, first, second)
            if frame_R is None or frame_P is None:
                continue
            same_component = owner[first] == owner[second]
            cross_R = np.linalg.norm(np.cross(
                _unit(coords_R[first] - coords_R[center]),
                _unit(coords_R[second] - coords_R[center])))
            cross_P = np.linalg.norm(np.cross(
                _unit(coords_P[first] - coords_P[center]),
                _unit(coords_P[second] - coords_P[center])))
            anchor_candidates.append((
                same_component, min(cross_R, cross_P),
                -first, -second, first, second, frame_R, frame_P))
        if not anchor_candidates:
            continue
        (*_score, first, second, frame_R, frame_P) = max(anchor_candidates)
        attached = []
        for component_index, component in enumerate(components):
            donors = sorted(set(neighbors).intersection(component))
            if donors:
                attached.append((tuple(sorted(component)), tuple(donors)))
        constraints.append({
            "center": center,
            "anchors": (first, second),
            "frame_R": frame_R,
            "frame_P": frame_P,
            "neighbors": tuple(neighbors),
            "attached": tuple(attached),
        })
    return tuple(constraints)


def _apply_coordination_constraint(coords, constraint, coords_R, coords_P, t):
    center = constraint["center"]
    center_xyz = coords[center].copy()
    frame_R, frame_P = constraint["frame_R"], constraint["frame_P"]
    frame_t = _interpolate_frame(frame_R, frame_P, t)
    targets = {}
    for donor in constraint["neighbors"]:
        vector_R = coords_R[donor] - coords_R[center]
        vector_P = coords_P[donor] - coords_P[center]
        distance = ((1.0 - t) * np.linalg.norm(vector_R)
                    + t * np.linalg.norm(vector_P))
        local_R = frame_R.T @ _unit(vector_R)
        local_P = frame_P.T @ _unit(vector_P)
        direction = frame_t @ _slerp_vector(local_R, local_P, t)
        targets[donor] = center_xyz + distance * direction
    for component, donors in constraint["attached"]:
        donor = donors[0]
        current_vector = coords[donor] - center_xyz
        target_vector = targets[donor] - center_xyz
        rotation = _rotation_between(current_vector, target_vector)
        atom_indices = np.asarray(component, dtype=int)
        moved = (coords[atom_indices] - center_xyz) @ rotation.T + center_xyz
        if len(donors) > 1:
            current_directions = np.stack([
                _unit(coords[d] - center_xyz) for d in donors])
            target_directions = np.stack([
                _unit(targets[d] - center_xyz) for d in donors])
            u, _s, vt = np.linalg.svd(
                current_directions.T @ target_directions)
            rotation = u @ vt
            if np.linalg.det(rotation) < 0:
                u[:, -1] *= -1
                rotation = u @ vt
            moved = (
                (coords[atom_indices] - center_xyz) @ rotation + center_xyz)
        donor_positions = {
            int(atom): moved[list(component).index(atom)] for atom in donors
        }
        correction = np.mean([
            targets[atom] - donor_positions[atom] for atom in donors], axis=0)
        coords[atom_indices] = moved + correction


def _basis(coords, origin, axis_atom, plane_atom):
    e1 = np.asarray(coords[axis_atom] - coords[origin], dtype=float)
    norm = np.linalg.norm(e1)
    if norm < 1e-10:
        e1 = np.array([1.0, 0.0, 0.0])
    else:
        e1 /= norm
    plane = np.asarray(coords[plane_atom] - coords[axis_atom], dtype=float)
    e3 = np.cross(e1, plane)
    norm = np.linalg.norm(e3)
    if norm < 1e-10:
        candidate = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(candidate, e1))) > 0.8:
            candidate = np.array([0.0, 1.0, 0.0])
        e3 = np.cross(e1, candidate)
        norm = np.linalg.norm(e3)
    e3 /= norm
    e2 = np.cross(e3, e1)
    return e1, e2, e3


def _local_internal(coords, atom, refs):
    origin, axis_atom, plane_atom = refs
    vector = np.asarray(coords[atom] - coords[origin], dtype=float)
    distance = float(np.linalg.norm(vector))
    if distance < 1e-12:
        return 0.0, 0.0, 0.0
    unit = vector / distance
    e1, e2, e3 = _basis(coords, origin, axis_atom, plane_atom)
    theta = float(np.arccos(np.clip(np.dot(unit, e1), -1.0, 1.0)))
    phi = float(np.arctan2(np.dot(unit, e3), np.dot(unit, e2)))
    return distance, theta, phi


def _place_from_internal(coords, refs, internal):
    origin, axis_atom, plane_atom = refs
    distance, theta, phi = internal
    e1, e2, e3 = _basis(coords, origin, axis_atom, plane_atom)
    direction = (
        np.cos(theta) * e1
        + np.sin(theta) * (np.cos(phi) * e2 + np.sin(phi) * e3)
    )
    return coords[origin] + distance * direction


def _zmatrix_plan(n_atoms, bonded_pairs, persistent_bonded_pairs):
    union = [set() for _ in range(n_atoms)]
    persistent = {tuple(sorted(map(int, pair)))
                  for pair in persistent_bonded_pairs}
    for left, right in {tuple(sorted(map(int, pair))) for pair in bonded_pairs}:
        union[left].add(right)
        union[right].add(left)
    unseen = set(range(n_atoms))
    components = []
    while unseen:
        root = min(unseen)
        parent = {root: None}
        order = []
        queue = deque([root])
        unseen.remove(root)
        while queue:
            atom = queue.popleft()
            order.append(atom)
            neighbors = sorted(
                union[atom],
                key=lambda other: (
                    tuple(sorted((atom, other))) not in persistent, other))
            for neighbor in neighbors:
                if neighbor not in unseen:
                    continue
                unseen.remove(neighbor)
                parent[neighbor] = atom
                queue.append(neighbor)
        components.append((order, parent))

    anchors = []
    entries = []
    for order, parent in components:
        component_anchors = order[:min(3, len(order))]
        anchors.extend(component_anchors)
        placed = list(component_anchors)
        for atom in order[len(component_anchors):]:
            origin = parent.get(atom)
            if origin is None or origin not in placed:
                origin = placed[0]
            candidates = []
            parent_origin = parent.get(origin)
            if parent_origin is not None:
                candidates.append(parent_origin)
            candidates.extend(reversed(placed))
            axis_atom = next(item for item in candidates if item != origin)
            candidates = []
            parent_axis = parent.get(axis_atom)
            if parent_axis is not None:
                candidates.append(parent_axis)
            candidates.extend(reversed(placed))
            plane_atom = next(
                item for item in candidates
                if item not in {origin, axis_atom})
            entries.append((atom, (origin, axis_atom, plane_atom)))
            placed.append(atom)
    return tuple(anchors), tuple(entries)


def _clash_report(coords, elements, bonded_pairs, threshold=0.70):
    coords = np.asarray(coords, dtype=float)
    bonded = {tuple(sorted(map(int, pair))) for pair in bonded_pairs}
    clashes, minimum = [], float("inf")
    for left in range(len(elements)):
        for right in range(left + 1, len(elements)):
            if (left, right) in bonded:
                continue
            radius_sum = (_COVALENT_RADII.get(str(elements[left]), 0.77)
                          + _COVALENT_RADII.get(str(elements[right]), 0.77))
            distance = float(np.linalg.norm(coords[left] - coords[right]))
            ratio = distance / radius_sum
            minimum = min(minimum, ratio)
            if ratio < float(threshold):
                clashes.append({"atoms": [left, right], "distance": distance,
                                "radius_ratio": ratio})
    clashes.sort(key=lambda item: item["radius_ratio"])
    return {"count": len(clashes),
            "minimum_radius_ratio": minimum if np.isfinite(minimum) else None,
            "pairs": clashes[:12], "threshold": float(threshold)}


def internal_coordinate_interpolation(
        coords_R, coords_P, elements, *, bonded_pairs=(),
        persistent_bonded_pairs=None, n_frames=101, clash_threshold=0.70):
    """Interpolate mapped local distance, angle, and signed torsion values."""
    coords_R = np.asarray(coords_R, dtype=float)
    coords_P = np.asarray(coords_P, dtype=float)
    if coords_R.shape != coords_P.shape or coords_R.ndim != 2:
        raise ValueError("R and P coordinates must have the same (N, 3) shape")
    if coords_R.shape[1] != 3 or len(elements) != len(coords_R):
        raise ValueError("coordinates and elements are inconsistent")
    persistent_bonded_pairs = (bonded_pairs if persistent_bonded_pairs is None
                               else persistent_bonded_pairs)
    product = proper_align_coordinates(coords_P, coords_R)
    anchors, entries = _zmatrix_plan(
        len(coords_R), bonded_pairs, persistent_bonded_pairs)
    coordination_constraints = _coordination_constraints(
        coords_R, product, bonded_pairs, persistent_bonded_pairs)
    endpoint_internals = [
        (_local_internal(coords_R, atom, refs),
         _local_internal(product, atom, refs))
        for atom, refs in entries
    ]
    frames = []
    for t in np.linspace(0.0, 1.0, max(2, int(n_frames))):
        coords = np.zeros_like(coords_R)
        reference = (1.0 - t) * coords_R + t * product
        coords[list(anchors)] = reference[list(anchors)]
        for (atom, refs), (internal_R, internal_P) in zip(
                entries, endpoint_internals):
            distance = (1.0 - t) * internal_R[0] + t * internal_P[0]
            theta = (1.0 - t) * internal_R[1] + t * internal_P[1]
            phi = internal_R[2] + t * _wrap_angle(
                internal_P[2] - internal_R[2])
            coords[atom] = _place_from_internal(
                coords, refs, (distance, theta, phi))
        for constraint in coordination_constraints:
            _apply_coordination_constraint(
                coords, constraint, coords_R, product, float(t))
        frames.append({"t": float(t), "coords": coords.tolist(),
                       "clashes": _clash_report(
                           coords, elements, bonded_pairs, clash_threshold)})
    return {
        "schema_version": "rxn_core.zmatrix_interpolation/v2",
        "method": (
            "mapped_internal_coordinates_with_spherical_coordination"),
        "primitive_counts": {"anchors": len(anchors),
                             "bond_angle_torsion": len(entries),
                             "coordination_centers": len(
                                 coordination_constraints)},
        "coordination_constraints": [{
            "center": constraint["center"],
            "anchors": list(constraint["anchors"]),
            "persistent_donors": list(constraint["neighbors"]),
            "anchor_rule": (
                "persistent_chelating_pair_then_maximum_angular_condition"),
        } for constraint in coordination_constraints],
        "n_frames": len(frames), "clash_threshold": float(clash_threshold),
        "frames": frames,
    }


# Backward-compatible name for callers from the earlier viewer prototype.
internal_distance_interpolation = internal_coordinate_interpolation

# Compatibility for tests and callers of the initial viewer prototype.
_proper_align = proper_align_coordinates

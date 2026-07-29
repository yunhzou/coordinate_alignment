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


def _pairwise_distances(coords):
    delta = np.asarray(coords, dtype=float)[:, None, :] - np.asarray(
        coords, dtype=float)[None, :, :]
    return np.linalg.norm(delta, axis=2)


def _bond_angle(coords, left, center, right):
    a, b = coords[left] - coords[center], coords[right] - coords[center]
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denominator < 1e-12:
        return 0.0
    return float(np.arccos(np.clip(np.dot(a, b) / denominator, -1.0, 1.0)))


def _dihedral(coords, first, second, third, fourth):
    b0 = coords[first] - coords[second]
    b1 = coords[third] - coords[second]
    b2 = coords[fourth] - coords[third]
    axis = _unit(b1)
    if axis is None:
        return 0.0
    v = b0 - np.dot(b0, axis) * axis
    w = b2 - np.dot(b2, axis) * axis
    if _unit(v) is None or _unit(w) is None:
        return 0.0
    return float(np.arctan2(np.dot(np.cross(axis, v), w), np.dot(v, w)))


def _proper_rotation(left, right):
    """Least-squares proper row-vector rotation from left to right."""
    u, _singular, vt = np.linalg.svd(
        np.asarray(left, dtype=float).T @ np.asarray(right, dtype=float))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0.0:
        u[:, -1] *= -1.0
        rotation = u @ vt
    return rotation


def _rotation_fraction(rotation, t):
    """Shortest SO(3) fraction of a proper rotation, including 180 degrees."""
    rotation = np.asarray(rotation, dtype=float)
    cosine = float(np.clip((np.trace(rotation) - 1.0) / 2.0, -1.0, 1.0))
    angle = float(np.arccos(cosine))
    if angle < 1e-12:
        return np.eye(3)
    values, vectors = np.linalg.eig(rotation.T)
    axis = _unit(np.real(vectors[:, np.argmin(np.abs(values - 1.0))]))
    if axis is None:
        return np.eye(3)
    pivot = int(np.argmax(np.abs(axis)))
    if axis[pivot] < 0.0:
        axis = -axis
    column_rotation = rotation.T
    sine_axis = np.array([
        column_rotation[2, 1] - column_rotation[1, 2],
        column_rotation[0, 2] - column_rotation[2, 0],
        column_rotation[1, 0] - column_rotation[0, 1],
    ]) / 2.0
    signed_angle = (-angle if np.dot(sine_axis, axis) < -1e-12 else angle)
    theta = float(t) * signed_angle
    skew = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    # Transpose the column-vector Rodrigues rotation for row vectors.
    return (np.eye(3) + np.sin(theta) * skew
            + (1.0 - np.cos(theta)) * (skew @ skew)).T


def _slerp_direction(left, right, t):
    left, right = _unit(left), _unit(right)
    if left is None or right is None:
        return left if right is None else right
    dot = float(np.clip(np.dot(left, right), -1.0, 1.0))
    if dot > 1.0 - 1e-10:
        return _unit((1.0 - t) * left + t * right)
    if dot < -1.0 + 1e-10:
        axis = np.eye(3)[int(np.argmin(np.abs(left)))]
        middle = _unit(np.cross(left, axis))
        return (_slerp_direction(left, middle, 2.0 * t) if t <= 0.5
                else _slerp_direction(middle, right, 2.0 * t - 1.0))
    angle = float(np.arccos(dot))
    return ((np.sin((1.0 - t) * angle) * left
             + np.sin(t * angle) * right) / np.sin(angle))


def _basis(coords, origin, axis_atom, plane_atom):
    e1 = _unit(coords[axis_atom] - coords[origin])
    if e1 is None:
        e1 = np.array([1.0, 0.0, 0.0])
    plane = coords[plane_atom] - coords[axis_atom]
    e3 = _unit(np.cross(e1, plane))
    if e3 is None:
        candidate = np.eye(3)[int(np.argmin(np.abs(e1)))]
        e3 = _unit(np.cross(e1, candidate))
    e2 = np.cross(e3, e1)
    return e1, e2, e3


def _local_internal(coords, atom, refs):
    origin, axis_atom, plane_atom = refs
    vector = coords[atom] - coords[origin]
    distance = float(np.linalg.norm(vector))
    if distance < 1e-12:
        return 0.0, 0.0, 0.0
    direction = vector / distance
    e1, e2, e3 = _basis(coords, origin, axis_atom, plane_atom)
    theta = float(np.arccos(np.clip(np.dot(direction, e1), -1.0, 1.0)))
    phi = float(np.arctan2(np.dot(direction, e3), np.dot(direction, e2)))
    return distance, theta, phi


def _place_internal(coords, refs, internal):
    origin, axis_atom, plane_atom = refs
    distance, theta, phi = internal
    e1, e2, e3 = _basis(coords, origin, axis_atom, plane_atom)
    direction = (np.cos(theta) * e1
                 + np.sin(theta) * (np.cos(phi) * e2 + np.sin(phi) * e3))
    return coords[origin] + distance * direction


def _continuous_basis(coords, refs, basis_R, basis_P, frame_rotation, t):
    origin, axis_atom, _plane_atom = refs
    current_axis = _unit(coords[axis_atom] - coords[origin])
    canonical = np.asarray(basis_R) @ _rotation_fraction(frame_rotation, t)
    e1 = current_axis
    e2 = canonical[1]
    e2 = _unit(e2 - e1 * np.dot(e2, e1))
    if e2 is None:
        e2 = canonical[2]
        e2 = _unit(e2 - e1 * np.dot(e2, e1))
    if e2 is None:
        e2 = _unit(np.cross(np.eye(3)[int(np.argmin(np.abs(e1)))], e1))
    e3 = _unit(np.cross(e1, e2))
    return e1, e2, e3


def _place_internal_continuous(coords, refs, internal, basis):
    origin, _axis_atom, _plane_atom = refs
    distance, theta, phi = internal
    e1, e2, e3 = basis
    direction = (np.cos(theta) * e1
                 + np.sin(theta) * (np.cos(phi) * e2 + np.sin(phi) * e3))
    return coords[origin] + distance * direction


def _analytic_internal_plans(coords_R, coords_P, persistent_bonded_pairs):
    """Deterministic spanning internal coordinates; no frame-time solving."""
    n_atoms = len(coords_R)
    adjacency = [set() for _ in range(n_atoms)]
    for left, right in persistent_bonded_pairs:
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(range(n_atoms))
    plans = []
    while unseen:
        component = set()
        seed = min(unseen)
        queue = deque([seed])
        unseen.remove(seed)
        while queue:
            atom = queue.popleft()
            component.add(atom)
            for neighbor in adjacency[atom]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    queue.append(neighbor)
        if len(component) < 4:
            continue
        root = min(component, key=lambda atom: (-len(adjacency[atom]), atom))
        parent = {root: None}
        order = []
        queue = deque([root])
        while queue:
            atom = queue.popleft()
            order.append(atom)
            for neighbor in sorted(
                    adjacency[atom], key=lambda item: (-len(adjacency[item]), item)):
                if neighbor not in component or neighbor in parent:
                    continue
                parent[neighbor] = atom
                queue.append(neighbor)
        second = order[1]
        third_candidates = [
            atom for atom in order[2:]
            if np.linalg.norm(np.cross(
                coords_R[second] - coords_R[root],
                coords_R[atom] - coords_R[root])) > 1e-5]
        third = third_candidates[0] if third_candidates else order[2]
        anchors = (root, second, third)
        placed = list(anchors)
        entries = []
        for atom in order:
            if atom in anchors:
                continue
            origin = parent[atom]
            if origin not in placed:
                # The selected third anchor may have displaced one BFS item;
                # postpone until its parent has been reconstructed.
                continue
            if origin == root:
                axis_atom, plane_atom = second, third
            elif parent[origin] == root:
                axis_atom = root
                plane_atom = third if origin != third else second
            else:
                axis_atom = parent[origin]
                plane_atom = parent[axis_atom]
            refs = (origin, axis_atom, plane_atom)
            internal_R = _local_internal(coords_R, atom, refs)
            internal_P = _local_internal(coords_P, atom, refs)
            basis_R = np.stack(_basis(coords_R, *refs))
            basis_P = np.stack(_basis(coords_P, *refs))
            entries.append((atom, refs, internal_R, internal_P,
                            basis_R, basis_P,
                            _proper_rotation(basis_R, basis_P)))
            placed.append(atom)
        star_groups = []
        for center in order:
            if center == root:
                continue
            children = [atom for atom in order
                        if parent.get(atom) == center and atom not in anchors]
            if not children:
                continue
            axis_atom = parent[center]
            if parent.get(axis_atom) is not None:
                plane_atom = parent[axis_atom]
            else:
                plane_atom = next(atom for atom in anchors
                                  if atom not in {center, axis_atom})
            refs = (center, axis_atom, plane_atom)
            basis_R = np.column_stack(_basis(coords_R, *refs))
            basis_P = np.column_stack(_basis(coords_P, *refs))
            local_R, local_P, distances = [], [], []
            for child in children:
                vector_R = coords_R[child] - coords_R[center]
                vector_P = coords_P[child] - coords_P[center]
                local_R.append(_unit(vector_R) @ basis_R)
                local_P.append(_unit(vector_P) @ basis_P)
                distances.append((float(np.linalg.norm(vector_R)),
                                  float(np.linalg.norm(vector_P))))
            local_R, local_P = np.stack(local_R), np.stack(local_P)
            # Both endpoint bases make the incoming parent bond +x.  The
            # collective freedom of the child star is therefore a signed
            # twist about x; allowing a general rotation would incorrectly
            # move children through the parent plane.
            cosine_sum = float(np.sum(
                local_R[:, 1] * local_P[:, 1]
                + local_R[:, 2] * local_P[:, 2]))
            sine_sum = float(np.sum(
                local_R[:, 1] * local_P[:, 2]
                - local_R[:, 2] * local_P[:, 1]))
            twist = float(np.arctan2(sine_sum, cosine_sum))
            cosine, sine = np.cos(twist), np.sin(twist)
            rotation = np.array([
                [1.0, 0.0, 0.0],
                [0.0, cosine, sine],
                [0.0, -sine, cosine],
            ])
            star_groups.append({
                "center": center,
                "axis_atom": axis_atom,
                "plane_atom": plane_atom,
                "children": tuple(children),
                "local_R": local_R,
                "local_P_in_R_frame": local_P @ rotation.T,
                "distances": tuple(distances),
                "rotation": rotation,
                "basis_R": np.stack(_basis(coords_R, *refs)),
                "basis_P": np.stack(_basis(coords_P, *refs)),
                "basis_rotation": _proper_rotation(
                    np.stack(_basis(coords_R, *refs)),
                    np.stack(_basis(coords_P, *refs))),
            })
        anchor_distances = tuple(
            (float(np.linalg.norm(coords_R[left] - coords_R[right])),
             float(np.linalg.norm(coords_P[left] - coords_P[right])))
            for left, right in ((root, second), (root, third), (second, third)))
        anchor_basis_R = np.stack(_basis(
            coords_R, root, second, third))
        anchor_basis_P = np.stack(_basis(
            coords_P, root, second, third))
        plans.append({"anchors": anchors,
                      "anchor_distances": anchor_distances,
                      "anchor_basis_R": anchor_basis_R,
                      "anchor_basis_P": anchor_basis_P,
                      "anchor_frame_rotation": _proper_rotation(
                          anchor_basis_R, anchor_basis_P),
                      "entries": tuple(entries),
                      "star_groups": tuple(star_groups)})
    return tuple(plans)


def _apply_analytic_internal_plans(coords, plans, t, basis_state):
    coords = np.asarray(coords, dtype=float).copy()
    for plan_index, plan in enumerate(plans):
        root, second, third = plan["anchors"]
        interpolated = [
            (1.0 - t) * pair[0] + t * pair[1]
            for pair in plan["anchor_distances"]]
        d_root_second, d_root_third, d_second_third = interpolated
        anchor_basis = (plan["anchor_basis_R"]
                        @ _rotation_fraction(
                            plan["anchor_frame_rotation"], t))
        axis = _unit(anchor_basis[0])
        if axis is not None:
            perpendicular = _unit(anchor_basis[1])
            coords[second] = coords[root] + d_root_second * axis
            along = ((d_root_third ** 2 - d_second_third ** 2
                      + d_root_second ** 2) / (2.0 * d_root_second))
            radial = np.sqrt(max(0.0, d_root_third ** 2 - along ** 2))
            coords[third] = (
                coords[root] + along * axis + radial * perpendicular)
        for (atom, refs, internal_R, internal_P, basis_R, basis_P,
             frame_rotation) in plan["entries"]:
            internal = (
                (1.0 - t) * internal_R[0] + t * internal_P[0],
                (1.0 - t) * internal_R[1] + t * internal_P[1],
                internal_R[2] + t * _wrap_angle(internal_P[2] - internal_R[2]),
            )
            basis = list(_basis(coords, *refs))
            key = (plan_index, "entry", atom)
            previous_e2 = basis_state.get(key, basis_R[1])
            if np.dot(basis[1], previous_e2) < 0.0:
                basis[1], basis[2] = -basis[1], -basis[2]
            basis_state[key] = basis[1].copy()
            coords[atom] = _place_internal_continuous(
                coords, refs, internal, basis)
        for group in plan["star_groups"]:
            refs = (group["center"], group["axis_atom"], group["plane_atom"])
            basis_rows = list(_basis(coords, *refs))
            key = (plan_index, "star", group["center"])
            previous_e2 = basis_state.get(key, group["basis_R"][1])
            if np.dot(basis_rows[1], previous_e2) < 0.0:
                basis_rows[1], basis_rows[2] = -basis_rows[1], -basis_rows[2]
            basis_state[key] = basis_rows[1].copy()
            basis = np.column_stack(basis_rows)
            rotation = _rotation_fraction(group["rotation"], t)
            for index, child in enumerate(group["children"]):
                local = _slerp_direction(
                    group["local_R"][index],
                    group["local_P_in_R_frame"][index], t)
                direction = _unit(local @ rotation @ basis.T)
                distance = ((1.0 - t) * group["distances"][index][0]
                            + t * group["distances"][index][1])
                coords[child] = coords[group["center"]] + distance * direction
    return coords


def _redundant_internal_constraints(coords_R, coords_P, bonded_pairs,
                                    persistent_bonded_pairs,
                                    reactant_bonded_pairs,
                                    product_bonded_pairs):
    """Build one consistent, redundant distance/angle constraint network.

    A Z-matrix only constrains a spanning tree.  That is intrinsically wrong
    for rings, chelates, and several fragments attached to one centre.  Here
    every atom pair is an internal coordinate, with stronger graph-derived
    constraints for bonds and angles.  Multiple constraints on the same pair
    are combined exactly as a weighted least-squares term.
    """
    n_atoms = len(coords_R)
    distances_R = _pairwise_distances(coords_R)
    distances_P = _pairwise_distances(coords_P)
    bonded = {tuple(sorted(map(int, pair))) for pair in bonded_pairs}
    persistent = {
        tuple(sorted(map(int, pair))) for pair in persistent_bonded_pairs}
    # Each tuple is (i, j, weight, endpoint_R, endpoint_P, kind).  The weak
    # complete graph is deliberate: it is the redundant all-atom internal
    # coordinate requested by the viewer, not a collision-repulsion term.
    terms = []
    for left in range(n_atoms):
        for right in range(left + 1, n_atoms):
            terms.append((left, right, 20.0, distances_R[left, right],
                          distances_P[left, right], "all_pair"))
    for left, right in sorted(bonded):
        is_persistent = (left, right) in persistent
        weight = 400.0 if is_persistent else 20.0
        terms.append((left, right, weight, distances_R[left, right],
                      distances_P[left, right],
                      ("bond", is_persistent)))

    # The distance between two neighbours plus the two bond lengths is an
    # exact representation of their bond angle.  Interpolate the angle itself
    # and convert it back to its opposite distance for the solver.
    angle_terms = []
    persistent_neighbors = [set() for _ in range(n_atoms)]
    union_neighbors = [set() for _ in range(n_atoms)]
    for left, right in bonded:
        union_neighbors[left].add(right)
        union_neighbors[right].add(left)
    for left, right in persistent:
        persistent_neighbors[left].add(right)
        persistent_neighbors[right].add(left)
    for center in range(n_atoms):
        # Only persistent legs define a molecular angle throughout the path.
        # An angle containing a forming/breaking leg does not exist at one
        # endpoint and must not distort the persistent scaffold.
        for left, right in combinations(
                sorted(persistent_neighbors[center]), 2):
            endpoint_targets = []
            endpoint_angles = []
            for coords, distances in ((coords_R, distances_R),
                                      (coords_P, distances_P)):
                a = coords[left] - coords[center]
                b = coords[right] - coords[center]
                denom = float(np.linalg.norm(a) * np.linalg.norm(b))
                cosine = (1.0 if denom < 1e-12 else
                          float(np.clip(np.dot(a, b) / denom, -1.0, 1.0)))
                endpoint_angles.append(float(np.arccos(cosine)))
                endpoint_targets.append(distances[left, right])
            angle_terms.append((center, left, right, tuple(endpoint_angles)))
            terms.append((left, right, 200.0, endpoint_angles[0],
                          endpoint_angles[1], ("angle", center)))
    # A changing leg attached beside a persistent leg needs a directional
    # guide; distances alone can otherwise send a transferring atom straight
    # through the persistent neighbour.  It remains soft because the changing
    # leg is not a molecular angle at both endpoints.
    for center in range(n_atoms):
        for left, right in combinations(sorted(union_neighbors[center]), 2):
            left_edge = tuple(sorted((center, left)))
            right_edge = tuple(sorted((center, right)))
            persistent_count = int(left_edge in persistent) + int(
                right_edge in persistent)
            if persistent_count != 1:
                continue
            theta_R = _bond_angle(coords_R, left, center, right)
            theta_P = _bond_angle(coords_P, left, center, right)
            terms.append((left, right, 80.0, theta_R, theta_P,
                          ("soft_angle", center)))
    rotatable_constraints = []
    for second, third in sorted(persistent):
        # A graph bridge is a genuine torsional freedom: removing it separates
        # two rigidly connected sides.  Ring bonds are deliberately excluded
        # because rotating one side would break ring closure.
        adjacency = [set(items) for items in persistent_neighbors]
        adjacency[second].discard(third)
        adjacency[third].discard(second)
        moving = {third}
        queue = deque([third])
        while queue:
            atom = queue.popleft()
            for neighbor in adjacency[atom]:
                if neighbor not in moving:
                    moving.add(neighbor)
                    queue.append(neighbor)
        if second in moving:
            continue
        candidates = []
        for first in sorted(persistent_neighbors[second] - {third}):
            for fourth in sorted(persistent_neighbors[third] - {second}):
                condition = []
                for coords in (coords_R, coords_P):
                    axis = _unit(coords[third] - coords[second])
                    left = _unit(coords[first] - coords[second])
                    right = _unit(coords[fourth] - coords[third])
                    if axis is None or left is None or right is None:
                        condition.append(0.0)
                    else:
                        condition.append(min(np.linalg.norm(np.cross(axis, left)),
                                             np.linalg.norm(np.cross(axis, right))))
                candidates.append((min(condition), -first, -fourth,
                                   first, fourth))
        if not candidates or max(candidates)[0] < 1e-5:
            continue
        _score, _neg_first, _neg_fourth, first, fourth = max(candidates)
        atoms = (first, second, third, fourth)
        phi_R = _dihedral(coords_R, *atoms)
        phi_P = _dihedral(coords_P, *atoms)
        rotatable_constraints.append({
            "atoms": atoms,
            "moving": tuple(sorted(moving)),
            "phi_R": phi_R,
            "phi_delta": _wrap_angle(phi_P - phi_R),
        })
    frame_constraints = []
    for center in range(n_atoms):
        neighbors = sorted(persistent_neighbors[center])
        if len(neighbors) < 3:
            continue
        directions_R = np.stack([
            _unit(coords_R[atom] - coords_R[center]) for atom in neighbors])
        directions_P = np.stack([
            _unit(coords_P[atom] - coords_P[center]) for atom in neighbors])
        if np.linalg.matrix_rank(directions_R, tol=1e-6) < 2:
            continue
        rotation = _proper_rotation(directions_R, directions_P)
        residual_P = directions_P @ rotation.T

        adjacency = [set(items) for items in persistent_neighbors]
        for neighbor in neighbors:
            adjacency[neighbor].discard(center)
        components = []
        unseen = set(range(n_atoms)) - {center}
        while unseen:
            root = min(unseen)
            unseen.remove(root)
            members = {root}
            queue = deque([root])
            while queue:
                atom = queue.popleft()
                for adjacent in adjacency[atom]:
                    if adjacent != center and adjacent in unseen:
                        unseen.remove(adjacent)
                        members.add(adjacent)
                        queue.append(adjacent)
            components.append(members)
        owned = []
        for members in components:
            donors = sorted(set(neighbors).intersection(members))
            if len(donors) == 1:
                owned.append((len(members), donors[0], tuple(sorted(members))))
        if len(owned) < 2:
            continue
        # The largest branch defines the local reference; all other graph-
        # separable branches rotate around the centre relative to it.
        fixed = max(owned, key=lambda item: (item[0], -item[1]))
        movable = []
        neighbor_index = {atom: index for index, atom in enumerate(neighbors)}
        for _size, donor, members in owned:
            if donor == fixed[1]:
                continue
            index = neighbor_index[donor]
            movable.append({
                "donor": donor,
                "atoms": members,
                "direction_R": directions_R[index],
                "direction_P_in_R_frame": residual_P[index],
            })
        if movable:
            fixed_index = neighbor_index[fixed[1]]
            frame_constraints.append({
                "center": center,
                "rotation": rotation,
                "fixed_donor": fixed[1],
                "fixed_direction_R": directions_R[fixed_index],
                "fixed_direction_P_in_R_frame": residual_P[fixed_index],
                "movable": tuple(movable),
            })
    component_pose_constraints = []
    unseen = set(range(n_atoms))
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        members = {root}
        queue = deque([root])
        while queue:
            atom = queue.popleft()
            for neighbor in persistent_neighbors[atom]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    members.add(neighbor)
                    queue.append(neighbor)
        atoms = np.asarray(sorted(members), dtype=int)
        if len(atoms) < 2:
            continue
        endpoint_R = coords_R[atoms]
        endpoint_P = coords_P[atoms]
        centroid_R = endpoint_R.mean(axis=0)
        centroid_P = endpoint_P.mean(axis=0)
        centered_R = endpoint_R - centroid_R
        centered_P = endpoint_P - centroid_P
        if len(atoms) >= 3 and np.linalg.matrix_rank(centered_R) >= 2:
            rotation = _proper_rotation(centered_R, centered_P)
        else:
            rotation = _row_rotation_between(
                centered_R[-1] - centered_R[0],
                centered_P[-1] - centered_P[0])
        component_pose_constraints.append({
            "atoms": tuple(map(int, atoms)),
            "centroid_R": centroid_R,
            "centroid_P": centroid_P,
            "centered_R": centered_R,
            "centered_P_in_R_frame": centered_P @ rotation.T,
            "rotation": rotation,
        })
    reactant_neighbors = [set() for _ in range(n_atoms)]
    product_neighbors = [set() for _ in range(n_atoms)]
    for left, right in reactant_bonded_pairs:
        reactant_neighbors[left].add(right)
        reactant_neighbors[right].add(left)
    for left, right in product_bonded_pairs:
        product_neighbors[left].add(right)
        product_neighbors[right].add(left)
    transfer_constraints = []
    for atom in range(n_atoms):
        if persistent_neighbors[atom]:
            continue
        if len(reactant_neighbors[atom]) != 1 or len(product_neighbors[atom]) != 1:
            continue
        donor = next(iter(reactant_neighbors[atom]))
        acceptor = next(iter(product_neighbors[atom]))
        if donor == acceptor:
            continue
        endpoint_frames = []
        valid = True
        for coords in (coords_R, coords_P):
            axis = _unit(coords[acceptor] - coords[donor])
            if axis is None:
                valid = False
                break
            atom_vector = coords[atom] - coords[donor]
            perpendicular = _unit(
                atom_vector - axis * np.dot(atom_vector, axis))
            if perpendicular is None:
                valid = False
                break
            endpoint_frames.append((
                axis, perpendicular, _unit(np.cross(axis, perpendicular))))
        if not valid:
            continue
        frame_R = np.stack(endpoint_frames[0])
        frame_P = np.stack(endpoint_frames[1])
        transfer_constraints.append({
            "atom": atom,
            "donor": donor,
            "acceptor": acceptor,
            "distance_donor_R": float(np.linalg.norm(
                coords_R[atom] - coords_R[donor])),
            "distance_donor_P": float(np.linalg.norm(
                coords_P[atom] - coords_P[donor])),
            "distance_acceptor_R": float(np.linalg.norm(
                coords_R[atom] - coords_R[acceptor])),
            "distance_acceptor_P": float(np.linalg.norm(
                coords_P[atom] - coords_P[acceptor])),
            "axis_R": endpoint_frames[0][0],
            "perpendicular_R": endpoint_frames[0][1],
            "frame_rotation": _proper_rotation(frame_R, frame_P),
        })
    return (terms, tuple(angle_terms), tuple(rotatable_constraints),
            tuple(frame_constraints), tuple(component_pose_constraints),
            tuple(transfer_constraints))


def _constraint_matrices(terms, t, n_atoms, angle_terms, coords_R, coords_P):
    weights = np.zeros((n_atoms, n_atoms), dtype=float)
    weighted_targets = np.zeros((n_atoms, n_atoms), dtype=float)
    distances_R = _pairwise_distances(coords_R)
    distances_P = _pairwise_distances(coords_P)
    for left, right, weight, target_R, target_P, kind in terms:
        if (isinstance(kind, tuple)
                and kind[0] in {"angle", "soft_angle"}):
            center = int(kind[1])
            bond_left = ((1.0 - t) * distances_R[center, left]
                         + t * distances_P[center, left])
            bond_right = ((1.0 - t) * distances_R[center, right]
                          + t * distances_P[center, right])
            theta = (1.0 - t) * target_R + t * target_P
            target = np.sqrt(max(
                0.0, bond_left * bond_left + bond_right * bond_right
                - 2.0 * bond_left * bond_right * np.cos(theta)))
        else:
            target = (1.0 - t) * target_R + t * target_P
        weights[left, right] += weight
        weights[right, left] += weight
        weighted_targets[left, right] += weight * target
        weighted_targets[right, left] += weight * target
    targets = np.divide(weighted_targets, weights, out=np.zeros_like(weights),
                        where=weights > 0.0)
    return weights, targets


def _stress_solve(initial, reference, weights, targets, *, tether=0.05,
                  system_inverse=None, max_iterations=120, tolerance=1e-9):
    """Weighted SMACOF with a small continuation tether.

    The previous frame is the reference.  Consequently the solver follows a
    continuous proper branch instead of independently embedding each frame
    (which can arbitrarily reflect a stereocentre).
    """
    coords = np.asarray(initial, dtype=float).copy()
    reference = np.asarray(reference, dtype=float)
    if system_inverse is None:
        laplacian = np.diag(weights.sum(axis=1)) - weights
        system = laplacian + float(tether) * np.eye(len(coords))
        system_inverse = np.linalg.inv(system)
    for _iteration in range(int(max_iterations)):
        distances = _pairwise_distances(coords)
        ratios = np.divide(weights * targets, distances,
                           out=np.zeros_like(weights), where=distances > 1e-12)
        b_matrix = -ratios
        np.fill_diagonal(b_matrix, ratios.sum(axis=1))
        updated = system_inverse @ (
            b_matrix @ coords + float(tether) * reference)
        if float(np.max(np.linalg.norm(updated - coords, axis=1))) < tolerance:
            coords = updated
            break
        coords = updated
    return coords


def _hard_local_targets(terms, t, coords_R, coords_P):
    """Return persistent bond and angle-side distances for projection."""
    distances_R = _pairwise_distances(coords_R)
    distances_P = _pairwise_distances(coords_P)
    targets = []
    for left, right, _weight, target_R, target_P, kind in terms:
        if isinstance(kind, tuple) and kind == ("bond", True):
            targets.append((left, right,
                            (1.0 - t) * target_R + t * target_P))
        elif isinstance(kind, tuple) and kind[0] == "angle":
            center = int(kind[1])
            bond_left = ((1.0 - t) * distances_R[center, left]
                         + t * distances_P[center, left])
            bond_right = ((1.0 - t) * distances_R[center, right]
                          + t * distances_P[center, right])
            theta = (1.0 - t) * target_R + t * target_P
            opposite = np.sqrt(max(
                0.0, bond_left * bond_left + bond_right * bond_right
                - 2.0 * bond_left * bond_right * np.cos(theta)))
            targets.append((left, right, opposite))
    return tuple(targets)


def _project_distance_constraints(coords, targets, *, max_iterations=1,
                                  tolerance=2e-7, relaxation=1.0):
    """SHAKE-style projection of the persistent molecular scaffold."""
    coords = np.asarray(coords, dtype=float).copy()
    if not targets:
        return coords
    left = np.asarray([item[0] for item in targets], dtype=int)
    right = np.asarray([item[1] for item in targets], dtype=int)
    desired = np.asarray([item[2] for item in targets], dtype=float)
    counts = np.bincount(
        np.concatenate((left, right)), minlength=len(coords)).astype(float)
    counts[counts == 0.0] = 1.0
    for _iteration in range(int(max_iterations)):
        vectors = coords[right] - coords[left]
        distances = np.linalg.norm(vectors, axis=1)
        errors = distances - desired
        if float(np.max(np.abs(errors))) < tolerance:
            break
        corrections = np.divide(
            0.5 * errors[:, None] * vectors, distances[:, None],
            out=np.zeros_like(vectors), where=distances[:, None] > 1e-12)
        displacement = np.zeros_like(coords)
        np.add.at(displacement, left, corrections)
        np.add.at(displacement, right, -corrections)
        coords += float(relaxation) * displacement / counts[:, None]
    return coords


def _rotate_about_axis(points, origin, axis, angle):
    axis = _unit(axis)
    if axis is None:
        return np.asarray(points, dtype=float).copy()
    points = np.asarray(points, dtype=float) - origin
    cosine, sine = np.cos(angle), np.sin(angle)
    return (points * cosine
            + np.cross(axis, points) * sine
            + np.outer(points @ axis, axis) * (1.0 - cosine)) + origin


def _row_rotation_between(left, right):
    left, right = _unit(left), _unit(right)
    if left is None or right is None:
        return np.eye(3)
    cosine = float(np.clip(np.dot(left, right), -1.0, 1.0))
    axis = _unit(np.cross(left, right))
    if axis is None:
        if cosine > 0.0:
            return np.eye(3)
        basis = np.eye(3)[int(np.argmin(np.abs(left)))]
        axis = _unit(np.cross(left, basis))
    angle = float(np.arccos(cosine))
    skew = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    column = (np.eye(3) + np.sin(angle) * skew
              + (1.0 - np.cos(angle)) * (skew @ skew))
    return column.T


def _apply_coordination_frame_constraints(coords, constraints, t):
    """Rotate graph-separable branches along signed local SO(3) paths."""
    coords = np.asarray(coords, dtype=float).copy()
    for constraint in constraints:
        center = constraint["center"]
        frame_rotation = _rotation_fraction(constraint["rotation"], t)
        for branch in constraint["movable"]:
            local = _slerp_direction(
                branch["direction_R"],
                branch["direction_P_in_R_frame"], t)
            target = _unit(local @ frame_rotation)
            current = _unit(coords[branch["donor"]] - coords[center])
            if target is None or current is None:
                continue
            cosine = float(np.clip(np.dot(current, target), -1.0, 1.0))
            axis = _unit(np.cross(current, target))
            if axis is None:
                if cosine > 0.0:
                    continue
                basis = np.eye(3)[int(np.argmin(np.abs(current)))]
                axis = _unit(np.cross(current, basis))
            angle = float(np.arccos(cosine))
            atoms = np.asarray(branch["atoms"], dtype=int)
            coords[atoms] = _rotate_about_axis(
                coords[atoms], coords[center], axis, angle)
    return coords


def _apply_component_pose_constraints(coords, constraints, t):
    """Follow each persistent component's shortest proper endpoint pose."""
    coords = np.asarray(coords, dtype=float).copy()
    for constraint in constraints:
        atoms = np.asarray(constraint["atoms"], dtype=int)
        rotation = _rotation_fraction(constraint["rotation"], t)
        residual_shape = (
            (1.0 - t) * constraint["centered_R"]
            + t * constraint["centered_P_in_R_frame"])
        target_shape = residual_shape @ rotation
        target_centroid = (
            (1.0 - t) * constraint["centroid_R"]
            + t * constraint["centroid_P"])
        coords[atoms] = target_shape + target_centroid
    return coords


def _persistent_fragment_plan(coords_R, coords_P, persistent_bonded_pairs):
    """Build rigid-fragment trees joined by graph-proven bridge bonds.

    A persistent bridge is a genuine inter-fragment freedom only when both
    sides contain at least two atoms.  Terminal bonds stay in their parent
    fragment because they cannot define a relative torsion by themselves.
    The construction is graph-only; no element or chemical motif rules enter.
    """
    n_atoms = len(coords_R)
    adjacency = [set() for _ in range(n_atoms)]
    edges = tuple(sorted({tuple(sorted(map(int, edge)))
                          for edge in persistent_bonded_pairs}))
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)

    cut_edges = []

    def split_component(members):
        # A fragment needs at least three atoms to carry an unambiguous 3-D
        # local pose.  Choose the most balanced admissible bridge, then apply
        # the same graph rule recursively to its two sides.
        candidates = []
        for left, right in edges:
            if left not in members or right not in members:
                continue
            side = {left}
            queue = deque([left])
            while queue:
                atom = queue.popleft()
                for neighbor in adjacency[atom]:
                    if neighbor not in members:
                        continue
                    if ((atom == left and neighbor == right)
                            or (atom == right and neighbor == left)):
                        continue
                    if neighbor not in side:
                        side.add(neighbor)
                        queue.append(neighbor)
            other = set(members) - side
            if right in side or min(len(side), len(other)) < 3:
                continue
            candidates.append((min(len(side), len(other)),
                               -abs(len(side) - len(other)),
                               -left, -right, (left, right), side, other))
        if not candidates:
            return
        _small, _balance, _left, _right, edge, side, other = max(candidates)
        cut_edges.append(edge)
        split_component(side)
        split_component(other)

    unseen_components = set(range(n_atoms))
    while unseen_components:
        seed = min(unseen_components)
        members = {seed}
        unseen_components.remove(seed)
        queue = deque([seed])
        while queue:
            atom = queue.popleft()
            for neighbor in adjacency[atom]:
                if neighbor in unseen_components:
                    unseen_components.remove(neighbor)
                    members.add(neighbor)
                    queue.append(neighbor)
        split_component(members)

    fragment_adjacency = [set(items) for items in adjacency]
    for left, right in cut_edges:
        fragment_adjacency[left].discard(right)
        fragment_adjacency[right].discard(left)
    atom_to_fragment = {}
    fragments = []
    unseen = set(range(n_atoms))
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        members = {seed}
        queue = deque([seed])
        while queue:
            atom = queue.popleft()
            for neighbor in fragment_adjacency[atom]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    members.add(neighbor)
                    queue.append(neighbor)
        fragment_index = len(fragments)
        atoms = np.asarray(sorted(members), dtype=int)
        for atom in atoms:
            atom_to_fragment[int(atom)] = fragment_index
        endpoint_R = coords_R[atoms]
        endpoint_P = coords_P[atoms]
        centroid_R = endpoint_R.mean(axis=0)
        centroid_P = endpoint_P.mean(axis=0)
        centered_R = endpoint_R - centroid_R
        centered_P = endpoint_P - centroid_P
        if len(atoms) >= 3 and np.linalg.matrix_rank(centered_R) >= 2:
            rotation_P = _proper_rotation(centered_R, centered_P)
        elif len(atoms) >= 2:
            rotation_P = _row_rotation_between(
                centered_R[-1] - centered_R[0],
                centered_P[-1] - centered_P[0])
        else:
            rotation_P = np.eye(3)
        atom_index = {int(atom): index for index, atom in enumerate(atoms)}
        local_neighbors = {int(atom): set() for atom in atoms}
        local_bonds = []
        for left, right in edges:
            if left not in members or right not in members:
                continue
            local_neighbors[left].add(right)
            local_neighbors[right].add(left)
            local_bonds.append((atom_index[left], atom_index[right]))
        local_angles = []
        for center in atoms:
            for left, right in combinations(
                    sorted(local_neighbors[int(center)]), 2):
                local_angles.append((atom_index[left], atom_index[int(center)],
                                     atom_index[right]))
        fragments.append({
            "atoms": atoms,
            "atom_index": atom_index,
            "centroid_R": centroid_R,
            "centroid_P": centroid_P,
            "local_R": centered_R,
            "local_P": centered_P @ rotation_P.T,
            "rotation_P": rotation_P,
            "distances_R": _pairwise_distances(centered_R),
            "distances_P": _pairwise_distances(
                centered_P @ rotation_P.T),
            "local_bonds": tuple(local_bonds),
            "local_angles": tuple(local_angles),
            "internal_plans": _analytic_internal_plans(
                centered_R, centered_P @ rotation_P.T, local_bonds),
        })

    joints = []
    fragment_neighbors = [set() for _ in fragments]
    for left, right in cut_edges:
        left_fragment = atom_to_fragment[left]
        right_fragment = atom_to_fragment[right]
        joint_index = len(joints)
        joints.append({"atoms": (left, right),
                       "fragments": (left_fragment, right_fragment)})
        fragment_neighbors[left_fragment].add((right_fragment, joint_index))
        fragment_neighbors[right_fragment].add((left_fragment, joint_index))

    roots = []
    parent = {}
    order = []
    unvisited = set(range(len(fragments)))
    while unvisited:
        connected = set()
        seed = min(unvisited)
        queue = deque([seed])
        connected.add(seed)
        while queue:
            fragment = queue.popleft()
            for neighbor, _joint in fragment_neighbors[fragment]:
                if neighbor not in connected:
                    connected.add(neighbor)
                    queue.append(neighbor)
        root = max(connected, key=lambda index: (
            len(fragments[index]["atoms"]), -index))
        roots.append(root)
        unvisited.difference_update(connected)
        parent[root] = None
        queue = deque([root])
        while queue:
            fragment = queue.popleft()
            order.append(fragment)
            for neighbor, joint_index in sorted(fragment_neighbors[fragment]):
                if neighbor in parent:
                    continue
                parent[neighbor] = (fragment, joint_index)
                queue.append(neighbor)

    for child, relation in parent.items():
        if relation is None:
            continue
        parent_fragment, joint_index = relation
        joint = joints[joint_index]
        left, right = joint["atoms"]
        left_fragment, right_fragment = joint["fragments"]
        if left_fragment == parent_fragment:
            parent_atom, child_atom = left, right
        else:
            parent_atom, child_atom = right, left
        parent_model = fragments[parent_fragment]
        child_model = fragments[child]
        vector_R = coords_R[child_atom] - coords_R[parent_atom]
        vector_P_local = ((coords_P[child_atom] - coords_P[parent_atom])
                          @ parent_model["rotation_P"].T)
        joint.update({
            "parent_fragment": parent_fragment,
            "child_fragment": child,
            "parent_atom": parent_atom,
            "child_atom": child_atom,
            "length_R": float(np.linalg.norm(vector_R)),
            "length_P": float(np.linalg.norm(vector_P_local)),
            "direction_R": _unit(vector_R),
            "direction_P_local": _unit(vector_P_local),
            "relative_rotation_P": (
                child_model["rotation_P"]
                @ parent_model["rotation_P"].T),
        })
    return {"fragments": tuple(fragments), "joints": tuple(joints),
            "roots": tuple(roots), "parent": parent, "order": tuple(order)}


def _apply_persistent_fragment_plan(plan, t, n_atoms, local_state=None):
    """Reconstruct one continuous hierarchical rigid-fragment frame."""
    coords = np.zeros((n_atoms, 3), dtype=float)
    rotations = {}
    local_state = {} if local_state is None else local_state
    for fragment_index in plan["order"]:
        model = plan["fragments"][fragment_index]
        raw_local = ((1.0 - t) * model["local_R"]
                     + t * model["local_P"])
        if t <= 0.0:
            local = model["local_R"].copy()
        elif t >= 1.0:
            local = model["local_P"].copy()
        elif len(model["atoms"]) < 2:
            local = raw_local
        else:
            # A deterministic local Z-tree carries bond lengths, angles, and
            # signed torsions.  It has no iterative embedding and therefore
            # no reflected solution or frame-time branch choice.
            basis_state = local_state.setdefault(
                ("basis", fragment_index), {})
            local = _apply_analytic_internal_plans(
                raw_local, model["internal_plans"], t, basis_state)
            local = proper_align_coordinates(local, raw_local)
        local_state[fragment_index] = local.copy()
        relation = plan["parent"][fragment_index]
        if relation is None:
            rotation = _rotation_fraction(model["rotation_P"], t)
            centroid = ((1.0 - t) * model["centroid_R"]
                        + t * model["centroid_P"])
        else:
            parent_fragment, joint_index = relation
            joint = plan["joints"][joint_index]
            parent_rotation = rotations[parent_fragment]
            relative = _rotation_fraction(
                joint["relative_rotation_P"], t)
            rotation = relative @ parent_rotation
            direction = _slerp_direction(
                joint["direction_R"], joint["direction_P_local"], t)
            length = ((1.0 - t) * joint["length_R"]
                      + t * joint["length_P"])
            parent_atom = joint["parent_atom"]
            child_atom = joint["child_atom"]
            child_local_index = model["atom_index"][child_atom]
            child_target = (coords[parent_atom]
                            + length * (direction @ parent_rotation))
            centroid = child_target - local[child_local_index] @ rotation
        rotations[fragment_index] = rotation
        coords[model["atoms"]] = local @ rotation + centroid
    return coords


def _apply_transfer_constraints(coords, constraints, t):
    """Move a changing one-bond atom around, rather than through, its axis."""
    coords = np.asarray(coords, dtype=float).copy()
    for constraint in constraints:
        donor = constraint["donor"]
        acceptor = constraint["acceptor"]
        axis = _unit(coords[acceptor] - coords[donor])
        if axis is None:
            continue
        canonical_rotation = _rotation_fraction(
            constraint["frame_rotation"], t)
        canonical_axis = _unit(constraint["axis_R"] @ canonical_rotation)
        perpendicular = _unit(
            constraint["perpendicular_R"] @ canonical_rotation)
        adjustment = _row_rotation_between(canonical_axis, axis)
        perpendicular = _unit(perpendicular @ adjustment)
        perpendicular = _unit(
            perpendicular - axis * np.dot(perpendicular, axis))
        if perpendicular is None:
            continue
        distance_donor = (
            (1.0 - t) * constraint["distance_donor_R"]
            + t * constraint["distance_donor_P"])
        distance_acceptor = (
            (1.0 - t) * constraint["distance_acceptor_R"]
            + t * constraint["distance_acceptor_P"])
        axis_length = float(np.linalg.norm(
            coords[acceptor] - coords[donor]))
        along = ((distance_donor ** 2 - distance_acceptor ** 2
                  + axis_length ** 2) / (2.0 * axis_length))
        radial = np.sqrt(max(0.0, distance_donor ** 2 - along ** 2))
        coords[constraint["atom"]] = (
            coords[donor] + along * axis + radial * perpendicular)
    return coords


def _refine_global_path(frames, grid, terms, coords_R, coords_P,
                        *, iterations=1000, learning_rate=0.020,
                        pose_weight=0.01, temporal_weight=1.0):
    """Optimize one continuous path against all internal coordinates."""
    initial = np.asarray(frames, dtype=float)
    path = initial.copy()
    n_atoms = path.shape[1]
    left, right = np.triu_indices(n_atoms, 1)
    weights, _unused = _constraint_matrices(
        terms, 0.0, n_atoms, (), coords_R, coords_P)
    pair_weights = weights[left, right]
    target_stack = []
    for t in grid:
        _weights, targets = _constraint_matrices(
            terms, float(t), n_atoms, (), coords_R, coords_P)
        target_stack.append(targets[left, right])
    target_stack = np.asarray(target_stack)
    atom_weight = np.bincount(
        np.concatenate((left, right)),
        weights=np.concatenate((pair_weights, pair_weights)),
        minlength=n_atoms)
    atom_weight = np.maximum(atom_weight, 1.0)
    first_moment = np.zeros_like(path)
    second_moment = np.zeros_like(path)
    beta1, beta2 = 0.9, 0.999
    for iteration in range(1, int(iterations) + 1):
        gradient = np.zeros_like(path)
        for frame in range(1, len(path) - 1):
            vectors = path[frame, left] - path[frame, right]
            distances = np.linalg.norm(vectors, axis=1)
            coefficient = np.divide(
                2.0 * pair_weights
                * (distances - target_stack[frame]), distances,
                out=np.zeros_like(distances), where=distances > 1e-12)
            pair_gradient = coefficient[:, None] * vectors
            np.add.at(gradient[frame], left, pair_gradient)
            np.add.at(gradient[frame], right, -pair_gradient)
        gradient /= atom_weight[None, :, None]
        gradient += float(pose_weight) * (path - initial)
        acceleration = path[2:] - 2.0 * path[1:-1] + path[:-2]
        gradient[:-2] += float(temporal_weight) * acceleration
        gradient[1:-1] -= 2.0 * float(temporal_weight) * acceleration
        gradient[2:] += float(temporal_weight) * acceleration
        gradient[[0, -1]] = 0.0
        first_moment = beta1 * first_moment + (1.0 - beta1) * gradient
        second_moment = (beta2 * second_moment
                         + (1.0 - beta2) * gradient * gradient)
        corrected_first = first_moment / (1.0 - beta1 ** iteration)
        corrected_second = second_moment / (1.0 - beta2 ** iteration)
        path -= float(learning_rate) * corrected_first / (
            np.sqrt(corrected_second) + 1e-8)
        path[0], path[-1] = coords_R, coords_P
    return path


def _remove_numeric_path_kinks(frames, *, ratio_limit=1.45,
                               max_passes=120):
    """Diffuse only numerical frame kinks; linear motion is unchanged."""
    path = np.asarray(frames, dtype=float).copy()
    for _pass in range(int(max_passes)):
        step = np.sqrt(np.mean((path[1:] - path[:-1]) ** 2, axis=(1, 2)))
        local = 0.5 * (step[:-2] + step[2:])
        ratios = np.divide(step[1:-1], np.maximum(local, 1e-12))
        bad = np.flatnonzero(ratios > float(ratio_limit)) + 1
        if not len(bad):
            break
        mask = np.zeros(len(path), dtype=bool)
        for index in bad:
            mask[max(1, index - 6):min(len(path) - 1, index + 7)] = True
        averaged = 0.25 * path[:-2] + 0.50 * path[1:-1] + 0.25 * path[2:]
        interior_mask = mask[1:-1]
        path[1:-1][interior_mask] = averaged[interior_mask]
    path[0], path[-1] = frames[0], frames[-1]
    return path


def _apply_rotatable_constraints(coords, constraints, t, *, sweeps=2):
    """Interpolate only graph-proven torsional freedoms by rigid rotation."""
    coords = np.asarray(coords, dtype=float).copy()
    epsilon = 1e-4
    for _sweep in range(int(sweeps)):
        for constraint in constraints:
            first, second, third, fourth = constraint["atoms"]
            target = constraint["phi_R"] + t * constraint["phi_delta"]
            current = _dihedral(coords, first, second, third, fourth)
            error = _wrap_angle(target - current)
            if abs(error) < 1e-10:
                continue
            axis = coords[third] - coords[second]
            trial = coords.copy()
            moving = np.asarray(constraint["moving"], dtype=int)
            trial[moving] = _rotate_about_axis(
                trial[moving], trial[second], axis, epsilon)
            slope = _wrap_angle(
                _dihedral(trial, first, second, third, fourth) - current)
            slope = slope / epsilon
            if abs(slope) < 0.5:
                continue
            coords[moving] = _rotate_about_axis(
                coords[moving], coords[second], axis, error / slope)
    return coords


def _constraint_report(coords, coords_R, coords_P, t, bonded_pairs,
                       persistent_bonded_pairs, angle_terms,
                       rotatable_constraints):
    distances = _pairwise_distances(coords)
    distances_R = _pairwise_distances(coords_R)
    distances_P = _pairwise_distances(coords_P)
    persistent = set(persistent_bonded_pairs)
    bond_errors, persistent_bond_errors = [], []
    for left, right in bonded_pairs:
        target = ((1.0 - t) * distances_R[left, right]
                  + t * distances_P[left, right])
        if target > 1e-12:
            error = abs(distances[left, right] - target) / target
            bond_errors.append(error)
            if (left, right) in persistent:
                persistent_bond_errors.append(error)
    angle_errors = []
    for center, left, right, endpoint_angles in angle_terms:
        a, b = coords[left] - coords[center], coords[right] - coords[center]
        denom = float(np.linalg.norm(a) * np.linalg.norm(b))
        if denom < 1e-12:
            continue
        actual = float(np.arccos(np.clip(np.dot(a, b) / denom, -1.0, 1.0)))
        target = ((1.0 - t) * endpoint_angles[0] + t * endpoint_angles[1])
        angle_errors.append(abs(actual - target))
    torsion_errors = []
    for constraint in rotatable_constraints:
        actual = _dihedral(coords, *constraint["atoms"])
        target = constraint["phi_R"] + t * constraint["phi_delta"]
        torsion_errors.append(abs(_wrap_angle(actual - target)))
    return {
        "max_bond_relative_error": max(bond_errors, default=0.0),
        "max_persistent_bond_relative_error": max(
            persistent_bond_errors, default=0.0),
        "max_angle_error_degrees": np.degrees(
            max(angle_errors, default=0.0)),
        "max_torsion_error_degrees": np.degrees(
            max(torsion_errors, default=0.0)),
    }


def internal_coordinate_interpolation(
        coords_R, coords_P, elements, *, bonded_pairs=(),
        persistent_bonded_pairs=None, reactant_bonded_pairs=None,
        product_bonded_pairs=None, n_frames=101, clash_threshold=0.70):
    """Interpolate a mapped redundant all-atom internal-coordinate network."""
    coords_R = np.asarray(coords_R, dtype=float)
    coords_P = np.asarray(coords_P, dtype=float)
    if coords_R.shape != coords_P.shape or coords_R.ndim != 2:
        raise ValueError("R and P coordinates must have the same (N, 3) shape")
    if coords_R.shape[1] != 3 or len(elements) != len(coords_R):
        raise ValueError("coordinates and elements are inconsistent")
    persistent_bonded_pairs = (bonded_pairs if persistent_bonded_pairs is None
                               else persistent_bonded_pairs)
    product = proper_align_coordinates(coords_P, coords_R)
    bonded_pairs = tuple(sorted({tuple(sorted(map(int, pair)))
                                 for pair in bonded_pairs}))
    persistent_bonded_pairs = tuple(sorted({tuple(sorted(map(int, pair)))
                                            for pair in
                                            persistent_bonded_pairs}))
    reactant_bonded_pairs = tuple(sorted({tuple(sorted(map(int, pair)))
                                          for pair in (
                                              bonded_pairs if
                                              reactant_bonded_pairs is None
                                              else reactant_bonded_pairs)}))
    product_bonded_pairs = tuple(sorted({tuple(sorted(map(int, pair)))
                                         for pair in (
                                             bonded_pairs if
                                             product_bonded_pairs is None
                                             else product_bonded_pairs)}))
    (terms, angle_terms, rotatable_constraints, frame_constraints,
     component_pose_constraints,
     transfer_constraints) = _redundant_internal_constraints(
         coords_R, product, bonded_pairs, persistent_bonded_pairs,
         reactant_bonded_pairs, product_bonded_pairs)
    fragment_plan = _persistent_fragment_plan(
        coords_R, product, persistent_bonded_pairs)
    grid = np.linspace(0.0, 1.0, max(2, int(n_frames)))
    initial_path = []
    fragment_local_state = {}
    for frame_index, t in enumerate(grid):
        if frame_index == 0:
            coords = coords_R.copy()
        elif frame_index == len(grid) - 1:
            coords = product.copy()
        else:
            # Each persistent fragment follows one local co-rotating path;
            # bridge joints carry the complete relative fragment pose.  No
            # frame is independently embedded or allowed to switch branches.
            coords = _apply_persistent_fragment_plan(
                fragment_plan, float(t), len(coords_R), fragment_local_state)
            coords = _apply_transfer_constraints(
                coords, transfer_constraints, float(t))
            coords = _apply_rotatable_constraints(
                coords, rotatable_constraints, float(t), sweeps=2)
        initial_path.append(coords)
    frames = []
    for t, coords in zip(grid, initial_path):
        residuals = _constraint_report(
            coords, coords_R, product, float(t), bonded_pairs,
            persistent_bonded_pairs, angle_terms, rotatable_constraints)
        frames.append({"t": float(t), "coords": coords.tolist(),
                       "constraint_residuals": residuals,
                       "clashes": _clash_report(
                           coords, elements, bonded_pairs, clash_threshold)})
    return {
        "schema_version": "rxn_core.fragment_kinematic_interpolation/v8",
        "method": "persistent_fragment_local_shape_and_joint_so3_path",
        "primitive_counts": {
            "all_atom_distances": len(coords_R) * (len(coords_R) - 1) // 2,
            "bonds": len(bonded_pairs),
            "angles": len(angle_terms),
            "persistent_fragments": len(fragment_plan["fragments"]),
            "interfragment_joints": len(fragment_plan["joints"]),
            "signed_transfer_paths": len(transfer_constraints),
        },
        "continuation": "analytic_no_iterations_no_frame_state",
        "n_frames": len(frames), "clash_threshold": float(clash_threshold),
        "frames": frames,
    }


# Backward-compatible name for callers from the earlier viewer prototype.
internal_distance_interpolation = internal_coordinate_interpolation

# Compatibility for tests and callers of the initial viewer prototype.
_proper_align = proper_align_coordinates

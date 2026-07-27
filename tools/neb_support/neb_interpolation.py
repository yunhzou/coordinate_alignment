"""Generate and audit endpoint-inclusive reaction-path previews.

The legacy Cartesian preview keeps the packaged reactant row order throughout:

``X_i(t) = R_i + t * (P_final_i - R_i)``, for ``0 <= t <= 1``.

The preferred preview uses ASE's image-dependent pair potential (IDPP), whose
targets are interpolated all-pair distances.  Native signed tetrahedral frames
select the correct mirror branch after IDPP.  Both paths are initial-path
diagnostics, not optimized NEB trajectories.
"""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from rxn_core.chemistry_computations import write_xyz_str


INTERPOLATION_VERSION = "rxn_core.neb_linear_interpolation/v1"
INTERNAL_INTERPOLATION_VERSION = (
    "rxn_core.internal_coordinate_interpolation/v1"
)
INTERNAL_PATH_MODEL = (
    "idpp_pair_distance_with_native_signed_frame_branch"
)

# Single-bond covalent radii in angstrom, rounded from Cordero et al.,
# Dalton Transactions 2008, DOI: 10.1039/B801115J.
COVALENT_RADII_ANGSTROM = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Al": 1.21,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Sc": 1.70,
    "Ti": 1.60,
    "V": 1.53,
    "Fe": 1.32,
    "Ni": 1.24,
    "Cu": 1.32,
    "Br": 1.20,
    "Ru": 1.46,
    "Rh": 1.42,
    "Pd": 1.39,
    "I": 1.39,
    "Pt": 1.36,
    "Au": 1.36,
}


class InterpolationError(ValueError):
    """The endpoint pair cannot form a trustworthy row-wise preview."""


def _validated_endpoints(
    elements: Sequence[str],
    coords_R: np.ndarray,
    coords_P_final: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    reactant = np.asarray(coords_R, dtype=float)
    product = np.asarray(coords_P_final, dtype=float)
    expected = (len(elements), 3)
    if reactant.shape != expected or product.shape != expected:
        raise InterpolationError(
            "R and P_final must match elements and have equal (n, 3) shape"
        )
    if not np.all(np.isfinite(reactant)) or not np.all(np.isfinite(product)):
        raise InterpolationError("endpoint coordinates contain nonfinite data")
    return reactant, product


def linear_images(
    coords_R: np.ndarray,
    coords_P_final: np.ndarray,
    image_count: int = 21,
) -> np.ndarray:
    """Return endpoint-inclusive Cartesian interpolation images."""
    reactant = np.asarray(coords_R, dtype=float)
    product = np.asarray(coords_P_final, dtype=float)
    if reactant.shape != product.shape or (
            reactant.ndim != 2 or reactant.shape[1] != 3):
        raise InterpolationError("R and P_final must have equal (n, 3) shape")
    if image_count < 2:
        raise InterpolationError("image_count must be at least two")
    fractions = np.linspace(0.0, 1.0, int(image_count))
    images = (
        reactant[None, :, :]
        + fractions[:, None, None] * (product - reactant)[None, :, :]
    )
    # Preserve the exact input arrays at both ends.
    images[0] = reactant
    images[-1] = product
    return images


def _load_ase_idpp():
    """Import the optional path dependency only when it is requested."""
    try:
        import ase
        from ase import Atoms
        from ase.mep import NEB
    except ImportError as exc:
        raise InterpolationError(
            "internal-coordinate interpolation requires ASE; install the "
            "optional dependency with `pip install 'rxn_core[path]'`"
        ) from exc
    return ase, Atoms, NEB


def _ase_idpp_images(
    elements: Sequence[str],
    reactant: np.ndarray,
    product: np.ndarray,
    image_count: int,
) -> tuple[np.ndarray, str]:
    """Run standard all-pair-distance IDPP through ASE."""
    ase, Atoms, NEB = _load_ase_idpp()
    atoms = [Atoms(symbols=list(elements), positions=reactant)]
    atoms.extend(
        Atoms(symbols=list(elements), positions=reactant)
        for _ in range(image_count - 2)
    )
    atoms.append(Atoms(symbols=list(elements), positions=product))
    # Use ASE's established default IDPP/NEB materialization.  The native
    # signed-frame audit below remains authoritative for handedness.
    band = NEB(atoms)
    band.interpolate(method="idpp")
    images = np.asarray([image.get_positions() for image in atoms], dtype=float)
    return images, str(getattr(ase, "__version__", "unknown"))


def _frame_points(frame: Mapping) -> tuple[int, ...]:
    center = int(frame["center"])
    neighbors = tuple(
        int(value) for value in frame["neighbors_R_index_order"]
    )
    if len(neighbors) == 4:
        return neighbors
    if len(neighbors) == 3:
        return (center, *neighbors)
    raise InterpolationError(
        "orientation frame must have three legacy neighbors or four "
        "native tetrahedral neighbors"
    )


def _orientation_measure(
    coords: np.ndarray,
    frame: Mapping,
) -> dict:
    """Evaluate one signed frame using only a machine-roundoff bound."""
    points = _frame_points(frame)
    xyz = np.asarray(coords, dtype=np.longdouble)
    origin = points[0]
    vectors = np.stack(
        [xyz[index] - xyz[origin] for index in points[1:]], axis=0
    )
    squared_lengths = np.sum(vectors * vectors, axis=1)
    zero_length = bool(np.any(squared_lengths == 0))
    denominator = np.sqrt(np.prod(squared_lengths))
    a, b, c = vectors
    positive_terms = (
        a[0] * b[1] * c[2],
        a[1] * b[2] * c[0],
        a[2] * b[0] * c[1],
    )
    negative_terms = (
        a[2] * b[1] * c[0],
        a[1] * b[0] * c[2],
        a[0] * b[2] * c[1],
    )
    determinant = sum(positive_terms) - sum(negative_terms)
    permanent = sum(
        abs(term) for term in (*positive_terms, *negative_terms)
    )
    eps = np.longdouble(np.finfo(np.longdouble).eps)
    operation_count = np.longdouble(16)
    gamma = operation_count * eps / (
        np.longdouble(1) - operation_count * eps
    )
    error_bound = gamma * max(permanent, denominator)
    defined = bool(
        not zero_length and abs(determinant) > error_bound
    )
    normalized = (
        np.longdouble(0)
        if zero_length
        else determinant / denominator
    )
    sign = 0 if not defined else (1 if determinant > 0 else -1)
    return {
        "defined": defined,
        "sign": sign,
        "normalized": float(normalized),
        "determinant": float(determinant),
        "determinant_error_bound": float(error_bound),
        "zero_length": zero_length,
    }


def _frame_reference_signs(
    reactant: np.ndarray,
    product: np.ndarray,
    frames: Sequence[Mapping],
) -> tuple[int, ...]:
    signs = []
    for frame_index, frame in enumerate(frames):
        measure_R = _orientation_measure(reactant, frame)
        measure_P = _orientation_measure(product, frame)
        if not measure_R["defined"] or not measure_P["defined"]:
            raise InterpolationError(
                f"orientation frame {frame_index} is undefined at an endpoint"
            )
        if measure_R["sign"] != measure_P["sign"]:
            raise InterpolationError(
                f"orientation frame {frame_index} has opposite endpoint signs"
            )
        signs.append(int(measure_R["sign"]))
    return tuple(signs)


def _frame_violations(
    coords: np.ndarray,
    frames: Sequence[Mapping],
    reference_signs: Sequence[int],
) -> tuple[int, ...]:
    violations = []
    for frame_index, (frame, reference_sign) in enumerate(
        zip(frames, reference_signs)
    ):
        measure = _orientation_measure(coords, frame)
        if not measure["defined"] or measure["sign"] != reference_sign:
            violations.append(frame_index)
    return tuple(violations)


def _reflect_point_across_plane(
    coords: np.ndarray,
    point: int,
    plane_points: Sequence[int],
) -> np.ndarray:
    """Reflect one tetrahedron point in the plane through the other three."""
    if len(plane_points) != 3:
        raise InterpolationError("a reflection plane needs exactly three points")
    out = np.asarray(coords, dtype=float).copy()
    xyz = np.asarray(coords, dtype=np.longdouble)
    first, second, third = (
        xyz[int(index)] for index in plane_points
    )
    normal = np.cross(second - first, third - first)
    squared_norm = np.dot(normal, normal)
    if squared_norm == 0:
        raise InterpolationError(
            "cannot repair a tetrahedral branch through a collinear plane"
        )
    displacement = (
        np.longdouble(2)
        * np.dot(xyz[int(point)] - first, normal)
        / squared_norm
    )
    out[int(point)] = np.asarray(
        xyz[int(point)] - displacement * normal,
        dtype=float,
    )
    return out


def _pair_distances(coords: np.ndarray) -> np.ndarray:
    xyz = np.asarray(coords, dtype=float)
    return np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)


def _relative_pair_distance_stress(
    coords: np.ndarray,
    target_distances: np.ndarray,
) -> float:
    distances = _pair_distances(coords)
    upper = np.triu_indices(len(distances), k=1)
    targets = np.asarray(target_distances, dtype=float)[upper]
    if np.any(targets <= 0.0):
        raise InterpolationError(
            "distinct endpoint atoms produce a nonpositive target distance"
        )
    residuals = (distances[upper] - targets) / targets
    return float(np.mean(residuals * residuals)) if len(residuals) else 0.0


def _repair_tetrahedral_branches(
    coords: np.ndarray,
    target_distances: np.ndarray,
    frames: Sequence[Mapping],
    reference_signs: Sequence[int],
) -> tuple[np.ndarray, dict]:
    """Greedily choose deterministic distance-preserving mirror branches."""
    current = np.asarray(coords, dtype=float).copy()
    initial = _frame_violations(current, frames, reference_signs)
    reflections = []
    while True:
        violations = _frame_violations(current, frames, reference_signs)
        if not violations:
            break
        candidates = []
        for frame_index in violations:
            frame = frames[frame_index]
            points = _frame_points(frame)
            frame_id = (int(frame["center"]), points, int(frame_index))
            for point in points:
                plane = tuple(value for value in points if value != point)
                candidate = _reflect_point_across_plane(
                    current, point, plane
                )
                candidate_violations = _frame_violations(
                    candidate, frames, reference_signs
                )
                stress = _relative_pair_distance_stress(
                    candidate, target_distances
                )
                candidates.append((
                    len(candidate_violations),
                    stress,
                    frame_id,
                    int(point),
                    candidate_violations,
                    candidate,
                ))
        best = min(candidates, key=lambda item: item[:4])
        if best[0] >= len(violations):
            raise InterpolationError(
                "local tetrahedral mirror repair could not reduce the hard-"
                "frame violation count"
            )
        current = best[5]
        reflections.append({
            "frame_index": int(best[2][2]),
            "center_R": int(best[2][0]),
            "reflected_atom_R": int(best[3]),
            "violation_count_before": len(violations),
            "violation_count_after": int(best[0]),
            "target_relative_pair_distance_stress": float(best[1]),
        })
        if len(reflections) > len(frames):
            raise InterpolationError(
                "local tetrahedral mirror repair did not terminate"
            )
    return current, {
        "initial_violation_count": len(initial),
        "final_violation_count": 0,
        "reflection_count": len(reflections),
        "reflections": reflections,
    }


def _proper_kabsch_to_reference(
    coords: np.ndarray,
    reference: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Apply a proper whole-image rigid gauge and return its determinant."""
    moving = np.asarray(coords, dtype=float)
    target = np.asarray(reference, dtype=float)
    if not len(moving):
        return moving.copy(), 1.0
    moving_centroid = np.mean(moving, axis=0)
    target_centroid = np.mean(target, axis=0)
    covariance = (
        (moving - moving_centroid).T
        @ (target - target_centroid)
    )
    left, _singular_values, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0.0:
        left[:, -1] *= -1.0
        rotation = left @ right_t
    determinant = float(np.linalg.det(rotation))
    aligned = (
        (moving - moving_centroid) @ rotation + target_centroid
    )
    return aligned, determinant


def _contiguous_runs(indices: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    runs = []
    for value in sorted(int(index) for index in indices):
        if not runs or value != runs[-1][-1] + 1:
            runs.append([value])
        else:
            runs[-1].append(value)
    return tuple(tuple(run) for run in runs)


def _path_hard_frame_violation_count(
    images: np.ndarray,
    frames: Sequence[Mapping],
    reference_signs: Sequence[int],
) -> tuple[int, int, int]:
    """Count invalid knots and exact roots over the whole materialized path."""
    knot_count = sum(
        len(_frame_violations(image, frames, reference_signs))
        for image in images
    )
    root_count = 0
    for start, end in zip(images[:-1], images[1:]):
        for frame in frames:
            neighbors = tuple(
                int(value)
                for value in frame["neighbors_R_index_order"]
            )
            coefficients = _volume_coefficients_exact(
                start, end, int(frame["center"]), neighbors
            )
            root_count += _exact_open_unit_root_count(coefficients)
    return knot_count + root_count, knot_count, root_count


def _repair_path_tetrahedral_branches(
    images: np.ndarray,
    targets: Sequence[np.ndarray],
    guides: Sequence[np.ndarray],
    frames: Sequence[Mapping],
    reference_signs: Sequence[int],
) -> tuple[np.ndarray, list[dict]]:
    """Repair each contiguous inverted run with one consistent point branch."""
    current = np.asarray(images, dtype=float).copy()
    records = [
        {
            "image_index": int(index),
            "reflection_count": 0,
            "reflections": [],
        }
        for index in range(len(current))
    ]
    while True:
        current_counts = _path_hard_frame_violation_count(
            current, frames, reference_signs
        )
        if current_counts[0] == 0:
            break
        candidates = []
        for frame_index, frame in enumerate(frames):
            violating_images = [
                image_index
                for image_index in range(1, len(current) - 1)
                if frame_index in _frame_violations(
                    current[image_index], frames, reference_signs
                )
            ]
            for run in _contiguous_runs(violating_images):
                points = _frame_points(frame)
                frame_id = (
                    int(frame["center"]), points, int(frame_index)
                )
                for point in points:
                    candidate = current.copy()
                    for image_index in run:
                        plane = tuple(
                            value for value in points if value != point
                        )
                        reflected = _reflect_point_across_plane(
                            candidate[image_index], point, plane
                        )
                        candidate[image_index], _determinant = (
                            _proper_kabsch_to_reference(
                                reflected, guides[image_index]
                            )
                        )
                    hard_counts = _path_hard_frame_violation_count(
                        candidate, frames, reference_signs
                    )
                    stress = sum(
                        _relative_pair_distance_stress(
                            candidate[index], targets[index]
                        )
                        for index in range(1, len(candidate) - 1)
                    )
                    candidates.append((
                        hard_counts[0],
                        stress,
                        frame_id,
                        int(run[0]),
                        int(point),
                        run,
                        hard_counts,
                        candidate,
                    ))
        if not candidates:
            break
        best = min(candidates, key=lambda item: item[:5])
        if best[0] >= current_counts[0]:
            break
        current = best[7]
        for image_index in best[5]:
            item = {
                "frame_index": int(best[2][2]),
                "center_R": int(best[2][0]),
                "reflected_atom_R": int(best[4]),
                "contiguous_run": [int(value) for value in best[5]],
                "hard_violation_count_before": int(current_counts[0]),
                "hard_violation_count_after": int(best[0]),
                "knot_violation_count_after": int(best[6][1]),
                "piecewise_root_count_after": int(best[6][2]),
                "target_relative_pair_distance_stress_sum": float(best[1]),
            }
            records[image_index]["reflection_count"] += 1
            records[image_index]["reflections"].append(item)
    return current, records


def internal_coordinate_images(
    elements: Sequence[str],
    coords_R: np.ndarray,
    coords_P_final: np.ndarray,
    orientation_frames: Sequence[Mapping],
    image_count: int = 21,
) -> tuple[np.ndarray, dict]:
    """Generate IDPP images on the native signed tetrahedral branch.

    IDPP determines shape from interpolated all-pair distances.  Those
    distances cannot distinguish mirror branches, so a violating local
    tetrahedron is repaired by reflecting one of its four points through the
    plane of the other three.  That operation leaves all six distances within
    the tetrahedron unchanged.  A proper Kabsch transform then chooses only
    the global rigid gauge of each image.
    """
    reactant, product = _validated_endpoints(
        elements, coords_R, coords_P_final
    )
    if image_count < 2:
        raise InterpolationError("image_count must be at least two")
    frames = tuple(orientation_frames)
    reference_signs = _frame_reference_signs(
        reactant, product, frames
    )
    if image_count == 2:
        raw = np.stack([reactant, product])
        ase_version = None
    else:
        raw, ase_version = _ase_idpp_images(
            elements, reactant, product, int(image_count)
        )
    if raw.shape != (int(image_count), len(elements), 3):
        raise InterpolationError("ASE IDPP returned an unexpected path shape")
    if not np.all(np.isfinite(raw)):
        raise InterpolationError("ASE IDPP returned nonfinite coordinates")
    raw = np.asarray(raw, dtype=float).copy()
    raw[0] = reactant
    raw[-1] = product

    distances_R = _pair_distances(reactant)
    distances_P = _pair_distances(product)
    fractions = np.linspace(0.0, 1.0, int(image_count))
    images = raw.copy()
    targets = []
    guides = []
    gauge_determinants = [1.0] * int(image_count)
    for fraction in fractions:
        value = float(fraction)
        targets.append(
            (1.0 - value) * distances_R + value * distances_P
        )
        guides.append(
            (1.0 - value) * reactant + value * product
        )
    for image_index in range(1, int(image_count) - 1):
        images[image_index], gauge_determinants[image_index] = (
            _proper_kabsch_to_reference(
                images[image_index], guides[image_index]
            )
        )
    images, repair_records = _repair_path_tetrahedral_branches(
        images, targets, guides, frames, reference_signs
    )
    image_records = []
    for image_index in range(1, int(image_count) - 1):
        repair = repair_records[image_index]
        repair.update({
            "image_index": int(image_index),
            "fraction": float(fractions[image_index]),
            "proper_rotation_determinant": float(
                gauge_determinants[image_index]
            ),
            "target_relative_pair_distance_stress": (
                _relative_pair_distance_stress(
                    images[image_index], targets[image_index]
                )
            ),
        })
        image_records.append(repair)

    # Preserve the caller's endpoint arrays exactly, after every transformation.
    images[0] = reactant
    images[-1] = product
    final_counts = _path_hard_frame_violation_count(
        images, frames, reference_signs
    )
    reflection_count = sum(
        int(record["reflection_count"]) for record in image_records
    )
    diagnostics = {
        "schema_version": INTERNAL_INTERPOLATION_VERSION,
        "path_model": INTERNAL_PATH_MODEL,
        "generator": "ASE standard IDPP with deterministic mirror-branch repair",
        "ase_version": ase_version,
        "image_count": int(image_count),
        "all_pair_distance_targets": True,
        "target_distance_formula": (
            "d_ij(t) = (1-t) d_ij(R) + t d_ij(P_final)"
        ),
        "orientation_frame_count": len(frames),
        "hard_frame_definition": (
            "long-double signed determinant with forward error derived "
            "only from machine epsilon"
        ),
        "empirical_volume_tolerance_used": False,
        "local_reflection_preserves_six_tetrahedron_distances": True,
        "local_reflection_count": reflection_count,
        "repaired_image_count": sum(
            bool(record["reflection_count"]) for record in image_records
        ),
        "final_hard_frame_violation_count": int(final_counts[0]),
        "final_hard_frame_knot_violation_count": int(final_counts[1]),
        "final_hard_frame_piecewise_root_count": int(final_counts[2]),
        "continuous_hard_frame_status": (
            "pass" if final_counts[0] == 0 else "fail"
        ),
        "proper_global_gauge_only": True,
        "endpoint_exact": bool(
            np.array_equal(images[0], reactant)
            and np.array_equal(images[-1], product)
        ),
        "images": image_records,
    }
    return images, diagnostics


def write_interpolation_xyz(
    path: str | Path,
    elements: Sequence[str],
    images: np.ndarray,
    *,
    step_id: str,
    mechanism_id: int,
) -> None:
    """Write an endpoint-inclusive internal-coordinate preview."""
    frames = np.asarray(images, dtype=float)
    if frames.ndim != 3 or frames.shape[1:] != (len(elements), 3):
        raise InterpolationError("trajectory shape does not match elements")
    if len(frames) < 2:
        raise InterpolationError("trajectory must contain at least two images")
    last = len(frames) - 1
    blocks = []
    for index, coords in enumerate(frames):
        fraction = index / last
        blocks.append(write_xyz_str(
            elements,
            coords,
            f"{step_id} mechanism {mechanism_id:03d}; "
            f"IDPP internal-coordinate preview {index}/{last}; "
            f"t={fraction:.8f}; not an optimized NEB path",
        ))
    Path(path).write_text("".join(blocks), encoding="utf-8")


def write_linear_xyz(
    path: str | Path,
    elements: Sequence[str],
    images: np.ndarray,
    *,
    step_id: str,
    mechanism_id: int,
) -> None:
    """Write an endpoint-inclusive, multi-frame XYZ preview."""
    frames = np.asarray(images, dtype=float)
    if frames.ndim != 3 or frames.shape[1:] != (len(elements), 3):
        raise InterpolationError("trajectory shape does not match elements")
    last = len(frames) - 1
    blocks = []
    for index, coords in enumerate(frames):
        fraction = index / last
        blocks.append(write_xyz_str(
            elements,
            coords,
            f"{step_id} mechanism {mechanism_id:03d}; "
            f"linear endpoint preview {index}/{last}; t={fraction:.8f}; "
            "not an optimized NEB path",
        ))
    Path(path).write_text("".join(blocks), encoding="utf-8")


def _fraction(value: float) -> Fraction:
    # Packaged XYZ coordinates have finite decimal representations. Converting
    # through str recovers that decimal, instead of the binary float fraction.
    return Fraction(str(float(value)))


def _determinant3(
    first: Sequence[Fraction],
    second: Sequence[Fraction],
    third: Sequence[Fraction],
) -> Fraction:
    return (
        first[0] * (second[1] * third[2] - second[2] * third[1])
        - first[1] * (second[0] * third[2] - second[2] * third[0])
        + first[2] * (second[0] * third[1] - second[1] * third[0])
    )


def _volume_coefficients_exact(
    coords_R: np.ndarray,
    coords_P: np.ndarray,
    center: int,
    neighbors: Sequence[int],
) -> tuple[Fraction, ...]:
    """Exact cubic coefficients for a linearly interpolated signed volume."""
    neighbors = tuple(int(value) for value in neighbors)
    if len(neighbors) == 4:
        # Native v3 affine tetrahedron: the first index-ordered neighbor is
        # the origin and the other three are its edge endpoints.
        origin = neighbors[0]
        endpoints = neighbors[1:]
    elif len(neighbors) == 3:
        # Legacy non-native/v2 center-to-three-neighbor frame.
        origin = int(center)
        endpoints = neighbors
    else:
        raise InterpolationError(
            "orientation frame must have three legacy neighbors or four "
            "native tetrahedral neighbors")
    reactant = [
        [_fraction(coords_R[index, axis]) for axis in range(3)]
        for index in (origin, *endpoints)
    ]
    product = [
        [_fraction(coords_P[index, axis]) for axis in range(3)]
        for index in (origin, *endpoints)
    ]
    r_vectors = []
    deltas = []
    for offset in range(1, 4):
        r_vector = [
            reactant[offset][axis] - reactant[0][axis]
            for axis in range(3)
        ]
        p_vector = [
            product[offset][axis] - product[0][axis]
            for axis in range(3)
        ]
        r_vectors.append(r_vector)
        deltas.append([
            p_vector[axis] - r_vector[axis] for axis in range(3)
        ])
    r1, r2, r3 = r_vectors
    d1, d2, d3 = deltas
    return (
        _determinant3(r1, r2, r3),
        (
            _determinant3(d1, r2, r3)
            + _determinant3(r1, d2, r3)
            + _determinant3(r1, r2, d3)
        ),
        (
            _determinant3(d1, d2, r3)
            + _determinant3(d1, r2, d3)
            + _determinant3(r1, d2, d3)
        ),
        _determinant3(d1, d2, d3),
    )


def _trim_polynomial(
    coefficients: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    values = list(coefficients)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def _polynomial_derivative(
    coefficients: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    return _trim_polynomial(tuple(
        Fraction(index) * value
        for index, value in enumerate(coefficients)
        if index
    ) or (Fraction(0),))


def _polynomial_divmod(
    dividend: Sequence[Fraction],
    divisor: Sequence[Fraction],
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    numerator = list(_trim_polynomial(dividend))
    denominator = _trim_polynomial(divisor)
    if len(denominator) == 1 and denominator[0] == 0:
        raise ZeroDivisionError("zero polynomial divisor")
    if len(numerator) < len(denominator):
        return (Fraction(0),), tuple(numerator)
    quotient = [Fraction(0)] * (len(numerator) - len(denominator) + 1)
    while len(numerator) >= len(denominator) and any(numerator):
        degree = len(numerator) - len(denominator)
        factor = numerator[-1] / denominator[-1]
        quotient[degree] = factor
        for index, value in enumerate(denominator):
            numerator[index + degree] -= factor * value
        numerator = list(_trim_polynomial(numerator))
    return _trim_polynomial(quotient), _trim_polynomial(numerator)


def _sturm_sequence(
    coefficients: Sequence[Fraction],
) -> tuple[tuple[Fraction, ...], ...]:
    first = _trim_polynomial(coefficients)
    if len(first) <= 1:
        return (first,)
    sequence = [first, _polynomial_derivative(first)]
    while not (
        len(sequence[-1]) == 1 and sequence[-1][0] == 0
    ):
        _quotient, remainder = _polynomial_divmod(
            sequence[-2], sequence[-1])
        if len(remainder) == 1 and remainder[0] == 0:
            break
        negated = tuple(-value for value in remainder)
        # Positive normalization controls Fraction growth without changing
        # signs or the Sturm variation count.
        scale = abs(negated[-1])
        sequence.append(tuple(value / scale for value in negated))
    return tuple(sequence)


def _polynomial_value(
    coefficients: Sequence[Fraction],
    value: Fraction,
) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def _sign_variations(
    sequence: Sequence[Sequence[Fraction]],
    value: Fraction,
) -> int:
    signs = []
    for polynomial in sequence:
        evaluated = _polynomial_value(polynomial, value)
        if evaluated:
            signs.append(1 if evaluated > 0 else -1)
    return sum(
        signs[index] != signs[index - 1]
        for index in range(1, len(signs))
    )


def _exact_open_unit_root_count(
    coefficients: Sequence[Fraction],
) -> int:
    polynomial = _trim_polynomial(coefficients)
    if len(polynomial) <= 1:
        return 0
    if _polynomial_value(polynomial, Fraction(0)) == 0 or (
            _polynomial_value(polynomial, Fraction(1)) == 0):
        raise InterpolationError(
            "defined orientation frame has a zero-volume endpoint")
    sequence = _sturm_sequence(polynomial)
    return (
        _sign_variations(sequence, Fraction(0))
        - _sign_variations(sequence, Fraction(1))
    )


def _numerical_open_unit_roots(
    coefficients: Sequence[Fraction],
) -> list[float]:
    polynomial = _trim_polynomial(coefficients)
    if len(polynomial) <= 1:
        return []
    roots = np.polynomial.polynomial.polyroots(
        np.asarray([float(value) for value in polynomial]))
    values = sorted(
        float(root.real)
        for root in roots
        if abs(float(root.imag)) <= 1.0e-8
        and 1.0e-10 < float(root.real) < 1.0 - 1.0e-10
    )
    unique = []
    for value in values:
        if not unique or abs(value - unique[-1]) > 1.0e-7:
            unique.append(value)
    return unique


def _audit_orientation_frames(
    coords_R: np.ndarray,
    coords_P: np.ndarray,
    frames: Sequence[Mapping],
    *,
    margin_samples: int,
    image_count: int,
) -> dict:
    failures = []
    worst_margin = None
    fractions = np.linspace(0.0, 1.0, margin_samples)
    for frame_index, frame in enumerate(frames):
        center = int(frame["center"])
        neighbors = tuple(
            int(value) for value in frame["neighbors_R_index_order"])
        if len(neighbors) not in {3, 4}:
            raise InterpolationError(
                f"orientation frame {frame_index} must have three legacy "
                "neighbors or four native tetrahedral neighbors")
        coefficients = _volume_coefficients_exact(
            coords_R, coords_P, center, neighbors)
        exact_count = _exact_open_unit_root_count(coefficients)
        numerical_roots = _numerical_open_unit_roots(coefficients)
        if exact_count != len(numerical_roots):
            raise InterpolationError(
                f"orientation frame {frame_index}: exact root count "
                f"{exact_count} != numerical locations "
                f"{len(numerical_roots)}")
        if len(neighbors) == 4:
            origin = neighbors[0]
            endpoints = neighbors[1:]
            orientation_model = "affine_four_neighbor_tetrahedron"
        else:
            origin = center
            endpoints = neighbors
            orientation_model = "center_to_three_neighbor_vectors"
        reactant_vectors = np.stack([
            coords_R[neighbor] - coords_R[origin]
            for neighbor in endpoints
        ])
        product_vectors = np.stack([
            coords_P[neighbor] - coords_P[origin]
            for neighbor in endpoints
        ])
        vectors = (
            reactant_vectors[None, :, :]
            + fractions[:, None, None]
            * (product_vectors - reactant_vectors)[None, :, :]
        )
        denominators = np.prod(
            np.linalg.norm(vectors, axis=2), axis=1)
        numerators = np.linalg.det(vectors)
        values = np.divide(
            numerators,
            denominators,
            out=np.zeros_like(numerators),
            where=denominators != 0.0,
        )
        reference_sign = 1.0 if values[0] > 0.0 else -1.0
        signed_values = reference_sign * values
        minimum_index = int(np.argmin(signed_values))
        minimum_signed = float(signed_values[minimum_index])
        minimum_t = float(fractions[minimum_index])
        frame_summary = {
            "frame_index": frame_index,
            "center_R": center,
            "neighbors_R_index_order": list(neighbors),
            "orientation_model": orientation_model,
            "reactant_normalized_orientation": float(values[0]),
            "product_normalized_orientation": float(values[-1]),
            "exact_interior_zero_count": exact_count,
            "interior_zero_t": numerical_roots,
            "minimum_sampled_signed_normalized_orientation": minimum_signed,
            "minimum_sampled_t": minimum_t,
            "nearest_preview_images": sorted({
                int(round(value * (image_count - 1)))
                for value in numerical_roots
            }),
        }
        if worst_margin is None or minimum_signed < (
                worst_margin[
                    "minimum_sampled_signed_normalized_orientation"]):
            worst_margin = frame_summary
        if exact_count:
            failures.append(frame_summary)
    return {
        "method": (
            "exact rational Sturm root count for the affine signed-volume "
            "cubic; numerical root locations; sampled normalized margin"
        ),
        "frame_policy": (
            "native four-neighbor tetrahedral hard frames recorded by core; "
            "legacy three-neighbor frames accepted for older inputs"
        ),
        "frame_count": len(frames),
        "frames_with_interior_zero_count": len(failures),
        "interior_zero_count": sum(
            item["exact_interior_zero_count"] for item in failures),
        "normalized_margin_sample_count": int(margin_samples),
        "continuous_orientation_status": (
            "fail" if failures else "pass"
        ),
        "zero_crossing_frames": failures,
        "worst_sampled_frame": worst_margin,
    }


def _pair_closest_approach(
    coords_R: np.ndarray,
    coords_P: np.ndarray,
    left: int,
    right: int,
) -> tuple[float, float, float, float]:
    initial = coords_R[left] - coords_R[right]
    final = coords_P[left] - coords_P[right]
    change = final - initial
    denominator = float(change @ change)
    if denominator <= 1.0e-24:
        fraction = 0.0
    else:
        fraction = float(np.clip(
            -(initial @ change) / denominator, 0.0, 1.0))
    minimum = float(np.linalg.norm(initial + fraction * change))
    return (
        fraction,
        minimum,
        float(np.linalg.norm(initial)),
        float(np.linalg.norm(final)),
    )


def _audit_pair_distances(
    elements: Sequence[str],
    coords_R: np.ndarray,
    coords_P: np.ndarray,
) -> dict:
    hard_failures = []
    warnings = []
    closest = None
    unavailable = sorted({
        str(element) for element in elements
        if str(element) not in COVALENT_RADII_ANGSTROM
    })
    for left in range(len(elements)):
        for right in range(left + 1, len(elements)):
            fraction, minimum, distance_R, distance_P = (
                _pair_closest_approach(
                    coords_R, coords_P, left, right))
            radii = (
                COVALENT_RADII_ANGSTROM.get(str(elements[left])),
                COVALENT_RADII_ANGSTROM.get(str(elements[right])),
            )
            radius_sum = (
                None if None in radii else float(radii[0] + radii[1]))
            hard_threshold = max(
                0.35,
                0.0 if radius_sum is None else 0.50 * radius_sum,
            )
            warning_threshold = max(
                0.50,
                0.0 if radius_sum is None else 0.70 * radius_sum,
            )
            record = {
                "atoms_R": [left, right],
                "elements": [
                    str(elements[left]), str(elements[right])],
                "t": fraction,
                "minimum_distance_angstrom": minimum,
                "reactant_distance_angstrom": distance_R,
                "product_distance_angstrom": distance_P,
                "covalent_radius_sum_angstrom": radius_sum,
                "hard_threshold_angstrom": hard_threshold,
                "warning_threshold_angstrom": warning_threshold,
            }
            if closest is None or minimum < (
                    closest["minimum_distance_angstrom"]):
                closest = record
            interior = 1.0e-10 < fraction < 1.0 - 1.0e-10
            if interior and minimum < hard_threshold:
                hard_failures.append(record)
            elif (
                interior
                and minimum < warning_threshold
                and minimum < 0.8 * min(distance_R, distance_P)
            ):
                warnings.append(record)
    hard_failures.sort(key=lambda item: item[
        "minimum_distance_angstrom"])
    warnings.sort(key=lambda item: item["minimum_distance_angstrom"])
    return {
        "method": (
            "analytic minimum distance for every pair of linear atom paths"
        ),
        "threshold_policy": {
            "hard": (
                "max(0.35 angstrom, 0.50 * covalent-radius sum)"
            ),
            "warning": (
                "max(0.50 angstrom, 0.70 * covalent-radius sum), "
                "and at least 20% below both endpoint distances"
            ),
            "interpretation": (
                "conservative engineering screen, not a universal "
                "chemical-law threshold"
            ),
        },
        "radius_source": (
            "Cordero et al., Dalton Trans. 2008, DOI 10.1039/B801115J"
        ),
        "elements_without_radius": unavailable,
        "hard_collision_count": len(hard_failures),
        "close_approach_warning_count": len(warnings),
        "continuous_collision_status": (
            "fail" if hard_failures
            else "review" if warnings or unavailable
            else "pass"
        ),
        "closest_pair": closest,
        "hard_collisions": hard_failures,
        "close_approach_warnings": warnings,
    }


def _audit_displacements(
    elements: Sequence[str],
    coords_R: np.ndarray,
    coords_P: np.ndarray,
    image_count: int,
) -> dict:
    distances = np.linalg.norm(coords_P - coords_R, axis=1)
    heavy = np.asarray([
        str(element) != "H" for element in elements], dtype=bool)

    def statistics(mask: np.ndarray) -> dict | None:
        values = distances[mask]
        if not len(values):
            return None
        indices = np.flatnonzero(mask)
        maximum_position = int(np.argmax(values))
        return {
            "atom_count": int(len(values)),
            "rms_angstrom": float(np.sqrt(np.mean(values ** 2))),
            "median_angstrom": float(np.median(values)),
            "maximum_angstrom": float(values[maximum_position]),
            "maximum_atom_R": int(indices[maximum_position]),
        }

    top_indices = sorted(
        range(len(distances)),
        key=lambda index: float(distances[index]),
        reverse=True,
    )[:10]
    maximum = float(np.max(distances)) if len(distances) else 0.0
    return {
        "interpretation": (
            "informational/manual-review metric; displacement alone does "
            "not invalidate a reaction path"
        ),
        "all_atoms": statistics(np.ones(len(elements), dtype=bool)),
        "heavy_atoms": statistics(heavy),
        "hydrogen_atoms": statistics(~heavy),
        "maximum_per_preview_interval_angstrom": (
            maximum / (image_count - 1)
        ),
        "above_five_angstrom_count": int(np.sum(distances > 5.0)),
        "status": "review" if maximum > 5.0 else "pass",
        "largest_atom_displacements": [
            {
                "atom_R": int(index),
                "element": str(elements[index]),
                "distance_angstrom": float(distances[index]),
            }
            for index in top_indices
        ],
    }


def _audit_piecewise_orientation_frames(
    images: np.ndarray,
    frames: Sequence[Mapping],
    *,
    margin_samples_per_segment: int,
) -> dict:
    """Audit every linear segment joining consecutive generated images."""
    if margin_samples_per_segment < 2:
        raise InterpolationError(
            "margin_samples_per_segment must be at least two"
        )
    path = np.asarray(images, dtype=float)
    segment_count = len(path) - 1
    local_fractions = np.linspace(
        0.0, 1.0, int(margin_samples_per_segment)
    )
    failures = []
    worst_margin = None
    total_roots = 0
    for frame_index, frame in enumerate(frames):
        center = int(frame["center"])
        neighbors = tuple(
            int(value) for value in frame["neighbors_R_index_order"]
        )
        reference = _orientation_measure(path[0], frame)
        if not reference["defined"]:
            raise InterpolationError(
                f"orientation frame {frame_index} is undefined at R"
            )
        reference_sign = int(reference["sign"])
        knot_failures = []
        for image_index, coords in enumerate(path):
            measure = _orientation_measure(coords, frame)
            if not measure["defined"] or measure["sign"] != reference_sign:
                knot_failures.append({
                    "image_index": int(image_index),
                    "path_t": float(image_index / segment_count),
                    "defined": bool(measure["defined"]),
                    "sign": int(measure["sign"]),
                    "normalized_orientation": float(
                        measure["normalized"]
                    ),
                    "determinant": float(measure["determinant"]),
                    "determinant_error_bound": float(
                        measure["determinant_error_bound"]
                    ),
                })

        roots = []
        minimum_signed = None
        minimum_t = None
        for segment_index in range(segment_count):
            start = path[segment_index]
            end = path[segment_index + 1]
            coefficients = _volume_coefficients_exact(
                start, end, center, neighbors
            )
            exact_count = _exact_open_unit_root_count(coefficients)
            numerical_roots = _numerical_open_unit_roots(coefficients)
            if exact_count != len(numerical_roots):
                raise InterpolationError(
                    f"orientation frame {frame_index}, segment "
                    f"{segment_index}: exact root count {exact_count} != "
                    f"numerical locations {len(numerical_roots)}"
                )
            for local_t in numerical_roots:
                roots.append({
                    "segment_index": int(segment_index),
                    "segment_t": float(local_t),
                    "path_t": float(
                        (segment_index + local_t) / segment_count
                    ),
                    "nearest_preview_image": int(round(
                        segment_index + local_t
                    )),
                })
            for local_t in local_fractions:
                coords = start + float(local_t) * (end - start)
                measure = _orientation_measure(coords, frame)
                signed = reference_sign * float(measure["normalized"])
                path_t = float(
                    (segment_index + float(local_t)) / segment_count
                )
                if minimum_signed is None or signed < minimum_signed:
                    minimum_signed = signed
                    minimum_t = path_t

        frame_summary = {
            "frame_index": int(frame_index),
            "center_R": center,
            "neighbors_R_index_order": list(neighbors),
            "orientation_model": (
                "affine_four_neighbor_tetrahedron"
                if len(neighbors) == 4
                else "center_to_three_neighbor_vectors"
            ),
            "reactant_normalized_orientation": float(
                reference["normalized"]
            ),
            "piecewise_interior_zero_count": len(roots),
            "interior_zero_t": [item["path_t"] for item in roots],
            "interior_zero_records": roots,
            "knot_violation_count": len(knot_failures),
            "knot_violations": knot_failures,
            "minimum_sampled_signed_normalized_orientation": float(
                minimum_signed
            ),
            "minimum_sampled_t": float(minimum_t),
        }
        if worst_margin is None or frame_summary[
            "minimum_sampled_signed_normalized_orientation"
        ] < worst_margin[
            "minimum_sampled_signed_normalized_orientation"
        ]:
            worst_margin = frame_summary
        if roots or knot_failures:
            failures.append(frame_summary)
        total_roots += len(roots)
    return {
        "method": (
            "exact rational Sturm root count for the signed-volume cubic "
            "on every adjacent image segment; long-double machine-roundoff "
            "classification at every knot"
        ),
        "frame_policy": (
            "native four-neighbor tetrahedral hard frames recorded by core; "
            "legacy three-neighbor frames accepted for older inputs"
        ),
        "frame_count": len(frames),
        "segment_count": segment_count,
        "frames_with_interior_zero_count": sum(
            bool(item["piecewise_interior_zero_count"])
            for item in failures
        ),
        "interior_zero_count": total_roots,
        "frames_with_knot_violation_count": sum(
            bool(item["knot_violation_count"]) for item in failures
        ),
        "margin_samples_per_segment": int(margin_samples_per_segment),
        "continuous_orientation_status": (
            "fail" if failures else "pass"
        ),
        "zero_crossing_frames": failures,
        "worst_sampled_frame": worst_margin,
        "empirical_volume_tolerance_used": False,
    }


def _audit_piecewise_pair_distances(
    elements: Sequence[str],
    images: np.ndarray,
) -> dict:
    """Find the analytic closest approach of every pair over every segment."""
    path = np.asarray(images, dtype=float)
    segment_count = len(path) - 1
    unavailable = sorted({
        str(element) for element in elements
        if str(element) not in COVALENT_RADII_ANGSTROM
    })
    pair_minima = []
    for left in range(len(elements)):
        for right in range(left + 1, len(elements)):
            best = None
            for segment_index, (start, end) in enumerate(
                zip(path[:-1], path[1:])
            ):
                local_t, minimum, _start_distance, _end_distance = (
                    _pair_closest_approach(
                        start, end, left, right
                    )
                )
                path_t = (segment_index + local_t) / segment_count
                key = (minimum, path_t, segment_index)
                if best is None or key < best[0]:
                    best = (key, segment_index, local_t, path_t, minimum)
            assert best is not None
            radii = (
                COVALENT_RADII_ANGSTROM.get(str(elements[left])),
                COVALENT_RADII_ANGSTROM.get(str(elements[right])),
            )
            radius_sum = (
                None if None in radii else float(radii[0] + radii[1])
            )
            hard_threshold = max(
                0.35,
                0.0 if radius_sum is None else 0.50 * radius_sum,
            )
            warning_threshold = max(
                0.50,
                0.0 if radius_sum is None else 0.70 * radius_sum,
            )
            pair_minima.append({
                "atoms_R": [left, right],
                "elements": [
                    str(elements[left]), str(elements[right])
                ],
                "t": float(best[3]),
                "segment_index": int(best[1]),
                "segment_t": float(best[2]),
                "minimum_distance_angstrom": float(best[4]),
                "reactant_distance_angstrom": float(np.linalg.norm(
                    path[0, left] - path[0, right]
                )),
                "product_distance_angstrom": float(np.linalg.norm(
                    path[-1, left] - path[-1, right]
                )),
                "covalent_radius_sum_angstrom": radius_sum,
                "hard_threshold_angstrom": hard_threshold,
                "warning_threshold_angstrom": warning_threshold,
            })

    closest = min(
        pair_minima,
        key=lambda item: (
            item["minimum_distance_angstrom"], item["atoms_R"]
        ),
        default=None,
    )
    hard_failures = []
    warnings = []
    for record in pair_minima:
        path_t = float(record["t"])
        interior = 1.0e-10 < path_t < 1.0 - 1.0e-10
        minimum = float(record["minimum_distance_angstrom"])
        distance_R = float(record["reactant_distance_angstrom"])
        distance_P = float(record["product_distance_angstrom"])
        if interior and minimum < record["hard_threshold_angstrom"]:
            hard_failures.append(record)
        elif (
            interior
            and minimum < record["warning_threshold_angstrom"]
            and minimum < 0.8 * min(distance_R, distance_P)
        ):
            warnings.append(record)
    hard_failures.sort(key=lambda item: (
        item["minimum_distance_angstrom"], item["atoms_R"]
    ))
    warnings.sort(key=lambda item: (
        item["minimum_distance_angstrom"], item["atoms_R"]
    ))
    return {
        "method": (
            "analytic minimum distance for every atom pair on every "
            "piecewise-linear segment joining generated images"
        ),
        "threshold_policy": {
            "hard": "max(0.35 angstrom, 0.50 * covalent-radius sum)",
            "warning": (
                "max(0.50 angstrom, 0.70 * covalent-radius sum), and at "
                "least 20% below both endpoint distances"
            ),
            "interpretation": (
                "conservative engineering screen, not a universal "
                "chemical-law threshold"
            ),
        },
        "radius_source": (
            "Cordero et al., Dalton Trans. 2008, DOI 10.1039/B801115J"
        ),
        "elements_without_radius": unavailable,
        "segment_count": segment_count,
        "hard_collision_count": len(hard_failures),
        "close_approach_warning_count": len(warnings),
        "continuous_collision_status": (
            "fail" if hard_failures
            else "review" if warnings or unavailable
            else "pass"
        ),
        "closest_pair": closest,
        "hard_collisions": hard_failures,
        "close_approach_warnings": warnings,
    }


def _audit_piecewise_displacements(
    elements: Sequence[str],
    images: np.ndarray,
) -> dict:
    path = np.asarray(images, dtype=float)
    report = _audit_displacements(
        elements, path[0], path[-1], len(path)
    )
    step_distances = np.linalg.norm(path[1:] - path[:-1], axis=2)
    per_atom_travel = np.sum(step_distances, axis=0)
    report["maximum_per_preview_interval_angstrom"] = (
        float(np.max(step_distances)) if step_distances.size else 0.0
    )
    report["maximum_total_atom_travel_angstrom"] = (
        float(np.max(per_atom_travel)) if len(per_atom_travel) else 0.0
    )
    report["rms_total_atom_travel_angstrom"] = (
        float(np.sqrt(np.mean(per_atom_travel ** 2)))
        if len(per_atom_travel) else 0.0
    )
    report["segment_count"] = len(path) - 1
    return report


def audit_internal_coordinate_interpolation(
    elements: Sequence[str],
    images: np.ndarray,
    orientation_frames: Sequence[Mapping],
    generation: Mapping | None = None,
    *,
    margin_samples_per_segment: int = 51,
) -> dict:
    """Audit a generated internal-coordinate path and return a JSON record."""
    path = np.asarray(images, dtype=float)
    if path.ndim != 3 or path.shape[1:] != (len(elements), 3):
        raise InterpolationError(
            "trajectory coordinates/elements differ"
        )
    if len(path) < 2:
        raise InterpolationError("trajectory must contain at least two images")
    if not np.all(np.isfinite(path)):
        raise InterpolationError("trajectory contains nonfinite coordinates")
    orientation = _audit_piecewise_orientation_frames(
        path,
        orientation_frames,
        margin_samples_per_segment=margin_samples_per_segment,
    )
    pairs = _audit_piecewise_pair_distances(elements, path)
    displacement = _audit_piecewise_displacements(elements, path)
    hard_fail = (
        orientation["continuous_orientation_status"] == "fail"
        or pairs["continuous_collision_status"] == "fail"
    )
    review = (
        pairs["continuous_collision_status"] == "review"
        or displacement["status"] == "review"
    )
    certified = not hard_fail and not review
    return {
        "schema_version": INTERNAL_INTERPOLATION_VERSION,
        "path_model": INTERNAL_PATH_MODEL,
        "formula": (
            "all-pair IDPP targets d_ij(t)=(1-t)d_ij(R)+t*d_ij(P_final); "
            "piecewise Cartesian materialization between generated images"
        ),
        "atom_order": "reactant rows at every image",
        "coordinate_units": "angstrom",
        "image_count": len(path),
        "endpoint_inclusive": True,
        "generation": None if generation is None else dict(generation),
        "checks": {
            "orientation": orientation,
            "pair_distance": pairs,
            "displacement": displacement,
        },
        "overall_status": (
            "fail" if hard_fail else "review" if review else "pass"
        ),
        "initial_internal_coordinate_path_certified": certified,
        "optimized_neb_path_certified": False,
        "interpretation": (
            "This audits the IDPP internal-coordinate initial guess and every "
            "linear segment joining its images. It is not an optimized NEB "
            "trajectory or a physical minimum-energy path."
        ),
    }


def audit_linear_interpolation(
    elements: Sequence[str],
    coords_R: np.ndarray,
    coords_P_final: np.ndarray,
    orientation_frames: Sequence[Mapping],
    *,
    image_count: int = 21,
    margin_samples: int = 1001,
) -> dict:
    """Audit continuous row-wise interpolation and return a JSON record."""
    reactant = np.asarray(coords_R, dtype=float)
    product = np.asarray(coords_P_final, dtype=float)
    if reactant.shape != product.shape or reactant.shape != (
            len(elements), 3):
        raise InterpolationError("endpoint coordinates/elements differ")
    if not np.all(np.isfinite(reactant)) or not np.all(np.isfinite(product)):
        raise InterpolationError("endpoint coordinates contain nonfinite data")
    images = linear_images(reactant, product, image_count)
    if not np.array_equal(images[0], reactant) or not np.array_equal(
            images[-1], product):
        raise InterpolationError("interpolation did not preserve endpoints")

    orientation = _audit_orientation_frames(
        reactant,
        product,
        orientation_frames,
        margin_samples=margin_samples,
        image_count=image_count,
    )
    pairs = _audit_pair_distances(elements, reactant, product)
    displacement = _audit_displacements(
        elements, reactant, product, image_count)
    hard_fail = (
        orientation["continuous_orientation_status"] == "fail"
        or pairs["continuous_collision_status"] == "fail"
    )
    review = (
        pairs["continuous_collision_status"] == "review"
        or displacement["status"] == "review"
    )
    return {
        "schema_version": INTERPOLATION_VERSION,
        "path_model": "linear_cartesian_endpoint_preview",
        "formula": "X_i(t) = R_i + t * (P_final_i - R_i)",
        "atom_order": "reactant rows at every image",
        "coordinate_units": "angstrom",
        "image_count": int(image_count),
        "endpoint_inclusive": True,
        "checks": {
            "orientation": orientation,
            "pair_distance": pairs,
            "displacement": displacement,
        },
        "overall_status": (
            "fail" if hard_fail else "review" if review else "pass"
        ),
        "initial_straight_line_certified": not hard_fail and not review,
        "optimized_neb_path_certified": False,
        "interpretation": (
            "This audits the straight-line initial guess only. It is not "
            "an optimized NEB trajectory or a physical minimum-energy path."
        ),
    }

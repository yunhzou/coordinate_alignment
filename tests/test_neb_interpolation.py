from __future__ import annotations

import inspect

import numpy as np
import pytest

from tools.neb_support import neb_interpolation
from tools.neb_support.neb_interpolation import (
    InterpolationError,
    audit_internal_coordinate_interpolation,
    audit_linear_interpolation,
    internal_coordinate_images,
    linear_images,
    write_interpolation_xyz,
    write_linear_xyz,
)


def _frame():
    return {
        "center": 0,
        "neighbors_R_index_order": [1, 2, 3, 4],
        "orientation_model": "affine_four_neighbor_tetrahedron",
        "reactant_normalized_orientation": 1.0,
        "shuffle_blocks": ["test"],
    }


def _legacy_frame():
    return {
        "center": 0,
        "neighbors_R_index_order": [1, 2, 3],
        "orientation_model": "center_to_three_neighbor_vectors",
        "reactant_normalized_orientation": 1.0,
        "shuffle_blocks": ["legacy"],
    }


def test_interpolation_orientation_audit_has_no_volume_threshold():
    parameters = inspect.signature(audit_linear_interpolation).parameters
    internal_parameters = inspect.signature(
        audit_internal_coordinate_interpolation).parameters

    assert "volume_tolerance" not in parameters
    assert "margin_tolerance" not in parameters
    assert "volume_tolerance" not in internal_parameters
    assert "margin_tolerance" not in internal_parameters


def test_linear_images_preserve_exact_endpoints():
    reactant = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    product = np.array([[0.5, 0.2, 0.1], [4.0, 5.0, 6.0]])

    images = linear_images(reactant, product, 5)

    assert images.shape == (5, 2, 3)
    assert np.array_equal(images[0], reactant)
    assert np.array_equal(images[-1], product)
    assert np.allclose(images[2], 0.5 * (reactant + product))


def test_continuous_audit_detects_flip_and_return_between_endpoints():
    # det[u1,u2,u3] = (t - 0.25) * (t - 0.75), so both endpoint
    # signs agree even though the local orientation flips twice inside.
    reactant = np.array([
        [-2.0, -2.0, -2.0],   # center; not the affine determinant origin
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, -0.25, 0.0],
        [0.0, 0.0, -0.75],
    ])
    product = np.array([
        [-2.0, -2.0, -2.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.75, 0.0],
        [0.0, 0.0, 0.25],
    ])

    report = audit_linear_interpolation(
        ["C", "F", "H", "H", "H"],
        reactant,
        product,
        [_frame()],
        image_count=5,
        margin_samples=101,
    )
    orientation = report["checks"]["orientation"]

    assert orientation["continuous_orientation_status"] == "fail"
    assert orientation["frames_with_interior_zero_count"] == 1
    assert orientation["interior_zero_count"] == 2
    assert np.allclose(
        orientation["zero_crossing_frames"][0]["interior_zero_t"],
        [0.25, 0.75],
    )
    assert report["initial_straight_line_certified"] is False


def test_continuous_audit_detects_atom_path_collision():
    reactant = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    product = reactant[::-1].copy()

    report = audit_linear_interpolation(
        ["C", "C"],
        reactant,
        product,
        [],
        image_count=3,
    )
    pairs = report["checks"]["pair_distance"]

    assert pairs["continuous_collision_status"] == "fail"
    assert pairs["hard_collision_count"] == 1
    assert pairs["closest_pair"]["t"] == 0.5
    assert pairs["closest_pair"]["minimum_distance_angstrom"] == 0.0


def test_rigid_translation_path_passes(tmp_path):
    reactant = np.array([
        [-2.0, -2.0, -2.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    product = reactant + np.array([0.2, -0.1, 0.3])
    elements = ["C", "F", "H", "H", "H"]

    report = audit_linear_interpolation(
        elements,
        reactant,
        product,
        [_frame()],
        image_count=5,
        margin_samples=101,
    )
    assert report["overall_status"] == "pass"
    assert report["initial_straight_line_certified"] is True

    path = tmp_path / "linear.xyz"
    write_linear_xyz(
        path,
        elements,
        linear_images(reactant, product, 5),
        step_id="TS_test",
        mechanism_id=1,
    )
    text = path.read_text()
    assert text.count("linear endpoint preview") == 5
    assert text.count("not an optimized NEB path") == 5


def test_legacy_three_neighbor_frame_remains_readable():
    reactant = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    product = reactant + np.array([0.1, 0.2, -0.1])

    report = audit_linear_interpolation(
        ["C", "H", "H", "H"],
        reactant,
        product,
        [_legacy_frame()],
        image_count=3,
        margin_samples=11,
    )

    orientation = report["checks"]["orientation"]
    assert orientation["continuous_orientation_status"] == "pass"
    assert orientation["worst_sampled_frame"]["orientation_model"] == (
        "center_to_three_neighbor_vectors")


def _tetrahedral_endpoints():
    reactant = np.array([
        [-2.0, -2.0, -2.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    product = reactant + np.array([0.2, -0.1, 0.3])
    return reactant, product


def _tetrahedron_distances(coords):
    points = np.asarray(coords)[[1, 2, 3, 4]]
    matrix = np.linalg.norm(
        points[:, None, :] - points[None, :, :], axis=2)
    return matrix[np.triu_indices(4, k=1)]


def test_internal_images_repair_local_mirror_branch(monkeypatch):
    reactant, product = _tetrahedral_endpoints()
    fractions = np.linspace(0.0, 1.0, 5)
    raw = (
        reactant[None]
        + fractions[:, None, None] * (product - reactant)[None]
    )
    # Reflect one ligand through the plane of the other three in image 2.
    raw[2, 4, 2] = -0.7
    wrong_local_distances = _tetrahedron_distances(raw[2])

    def fake_idpp(_elements, _reactant, _product, _image_count):
        return raw.copy(), "test"

    monkeypatch.setattr(
        neb_interpolation, "_ase_idpp_images", fake_idpp)
    images, generation = internal_coordinate_images(
        ["C", "F", "H", "H", "H"],
        reactant,
        product,
        [_frame()],
        image_count=5,
    )

    assert np.array_equal(images[0], reactant)
    assert np.array_equal(images[-1], product)
    assert generation["local_reflection_count"] == 1
    assert generation["repaired_image_count"] == 1
    assert generation["final_hard_frame_violation_count"] == 0
    assert generation["empirical_volume_tolerance_used"] is False
    assert np.allclose(
        _tetrahedron_distances(images[2]), wrong_local_distances)
    assert np.linalg.det(np.stack([
        images[2, 2] - images[2, 1],
        images[2, 3] - images[2, 1],
        images[2, 4] - images[2, 1],
    ])) > 0.0


def test_two_image_internal_path_needs_no_optional_dependency():
    reactant, product = _tetrahedral_endpoints()

    images, generation = internal_coordinate_images(
        ["C", "F", "H", "H", "H"],
        reactant,
        product,
        [_frame()],
        image_count=2,
    )

    assert np.array_equal(images, np.stack([reactant, product]))
    assert generation["ase_version"] is None
    assert generation["endpoint_exact"] is True


def test_internal_generator_rejects_opposite_endpoint_signs():
    reactant, product = _tetrahedral_endpoints()
    product = product.copy()
    product[4, 2] = -0.7

    with pytest.raises(InterpolationError, match="opposite endpoint signs"):
        internal_coordinate_images(
            ["C", "F", "H", "H", "H"],
            reactant,
            product,
            [_frame()],
            image_count=2,
        )


def test_piecewise_internal_audit_detects_flip_and_return():
    reactant = np.array([
        [-2.0, -2.0, -2.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, -0.25, 0.0],
        [0.0, 0.0, -0.75],
    ])
    product = np.array([
        [-2.0, -2.0, -2.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.75, 0.0],
        [0.0, 0.0, 0.25],
    ])

    report = audit_internal_coordinate_interpolation(
        ["C", "F", "H", "H", "H"],
        np.stack([reactant, product]),
        [_frame()],
        margin_samples_per_segment=101,
    )

    orientation = report["checks"]["orientation"]
    assert orientation["interior_zero_count"] == 2
    assert orientation["continuous_orientation_status"] == "fail"
    assert report["initial_internal_coordinate_path_certified"] is False
    assert report["optimized_neb_path_certified"] is False


def test_piecewise_pair_audit_counts_collision_at_internal_knot():
    images = np.array([
        [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    ])

    report = audit_internal_coordinate_interpolation(
        ["C", "C"], images, [])
    pairs = report["checks"]["pair_distance"]

    assert pairs["hard_collision_count"] == 1
    assert pairs["closest_pair"]["t"] == 0.5
    assert pairs["closest_pair"]["minimum_distance_angstrom"] == 0.0


def test_internal_xyz_writer_labels_generated_path(tmp_path):
    reactant, product = _tetrahedral_endpoints()
    images = np.stack([reactant, product])
    path = tmp_path / "internal.xyz"

    write_interpolation_xyz(
        path,
        ["C", "F", "H", "H", "H"],
        images,
        step_id="TS_test",
        mechanism_id=2,
    )

    text = path.read_text()
    assert text.count("IDPP internal-coordinate preview") == 2
    assert text.count("not an optimized NEB path") == 2

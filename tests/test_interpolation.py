import numpy as np

from rxn_core.alignment.interpolation import (
    _dihedral,
    _rotation_fraction,
    _wrap_angle,
    internal_coordinate_interpolation,
    internal_distance_interpolation,
    proper_align_coordinates,
)


def test_half_turn_follows_proper_rotation_without_crossing_center_plane():
    half_turn = np.diag([-1.0, -1.0, 1.0])
    vector = np.array([1.0, 0.0, 0.0])
    path = np.stack([
        vector @ _rotation_fraction(half_turn, t)
        for t in np.linspace(0.0, 1.0, 101)
    ])

    assert np.allclose(np.linalg.norm(path, axis=1), 1.0)
    assert np.linalg.norm(path[50]) == 1.0
    assert np.max(np.linalg.norm(np.diff(path, axis=0), axis=1)) < 0.04


def test_internal_distance_interpolation_preserves_endpoints_and_frame_count():
    coords_r = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    coords_p = np.array([
        [2.0, 3.0, 4.0],
        [2.0, 4.2, 4.0],
        [0.8, 3.0, 4.0],
        [2.0, 3.0, 5.2],
    ])

    result = internal_distance_interpolation(
        coords_r, coords_p, ["C", "H", "H", "H"],
        bonded_pairs=[(0, 1), (0, 2), (0, 3)], n_frames=9)

    assert result["n_frames"] == 9
    assert len(result["frames"]) == 9
    assert np.allclose(result["frames"][0]["coords"], coords_r)
    assert np.allclose(
        sorted(np.linalg.norm(
            np.asarray(result["frames"][-1]["coords"])[0]
            - np.asarray(result["frames"][-1]["coords"])[1:], axis=1)),
        [1.2, 1.2, 1.2])


def test_internal_distance_interpolation_reports_nonbonded_clash():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
    ])
    result = internal_distance_interpolation(
        coords, coords, ["C", "H", "H"], bonded_pairs=[(0, 1)],
        n_frames=2)

    report = result["frames"][0]["clashes"]
    assert report["count"] >= 1
    assert report["pairs"][0]["atoms"] == [0, 2]


def test_midpoint_preserves_redundant_bonds_angles_and_signed_torsion():
    coords_r = np.array([
        [0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.1, 1.1, 0.0],
        [3.1, 1.4, 0.8], [3.8, 2.2, 1.5],
    ])
    coords_p = np.array([
        [0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.2, 1.0, 0.0],
        [3.0, 1.7, -0.9], [3.5, 2.7, -1.7],
    ])
    bonds = [(0, 1), (1, 2), (2, 3), (3, 4)]
    result = internal_coordinate_interpolation(
        coords_r, coords_p, ["C"] * 5, bonded_pairs=bonds,
        persistent_bonded_pairs=bonds, n_frames=3)
    midpoint = np.asarray(result["frames"][1]["coords"])
    residuals = result["frames"][1]["constraint_residuals"]
    assert residuals["max_bond_relative_error"] < 1e-3
    assert residuals["max_angle_error_degrees"] < 0.1
    assert residuals["max_torsion_error_degrees"] < 0.2
    for atoms in ((0, 1, 2, 3), (1, 2, 3, 4)):
        phi_r = _dihedral(coords_r, *atoms)
        phi_p = _dihedral(coords_p, *atoms)
        expected = phi_r + 0.5 * _wrap_angle(phi_p - phi_r)
        assert abs(_wrap_angle(_dihedral(midpoint, *atoms) - expected)) < (
            np.radians(0.2))


def test_high_coordinate_donor_stays_in_anchored_orientation_hemisphere():
    coords_r = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0], [-1.0, 0.0, 0.2],
        [0.0, -1.0, -0.5], [0.5, 0.2, -1.0],
        [1.0, 1.0, 0.0],
    ])
    coords_p = coords_r.copy()
    coords_p[3] = [0.7, -0.1, 0.7]
    bonds = [(0, atom) for atom in range(1, 7)] + [(1, 7), (2, 7)]

    result = internal_coordinate_interpolation(
        coords_r, coords_p, ["Sc"] + ["O"] * 6 + ["C"],
        bonded_pairs=bonds, persistent_bonded_pairs=bonds, n_frames=21)

    assert result["schema_version"] == (
        "rxn_core.fragment_kinematic_interpolation/v8")
    assert result["primitive_counts"]["angles"] >= 15
    determinants = []
    for frame in result["frames"]:
        xyz = np.asarray(frame["coords"])
        determinants.append(np.linalg.det(np.stack([
            xyz[1] - xyz[0], xyz[2] - xyz[0], xyz[3] - xyz[0],
        ])))
    assert min(determinants) > 0
    assert max(
        frame["constraint_residuals"]["max_angle_error_degrees"]
        for frame in result["frames"]) < 2.0
    assert np.allclose(result["frames"][0]["coords"], coords_r)
    assert np.allclose(
        result["frames"][-1]["coords"],
        proper_align_coordinates(coords_p, coords_r))


def test_two_rigid_fragments_follow_one_smooth_bridge_rotation():
    coords_r = np.array([
        [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0],
    ])
    coords_p = coords_r.copy()
    angle = np.radians(140.0)
    rotation = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(angle), np.sin(angle)],
        [0.0, -np.sin(angle), np.cos(angle)],
    ])
    coords_p[3:] = ((coords_r[3:] - coords_r[3]) @ rotation
                    + coords_r[3])
    bonds = [
        (0, 1), (1, 2), (0, 2),
        (3, 4), (4, 5), (3, 5),
        (2, 3),
    ]

    result = internal_coordinate_interpolation(
        coords_r, coords_p, ["C"] * 6, bonded_pairs=bonds,
        persistent_bonded_pairs=bonds, n_frames=101)

    assert result["primitive_counts"]["persistent_fragments"] == 2
    assert result["primitive_counts"]["interfragment_joints"] == 1
    path = np.asarray([frame["coords"] for frame in result["frames"]])
    steps = np.sqrt(np.mean(np.diff(path, axis=0) ** 2, axis=(1, 2)))
    assert np.max(steps) / np.median(steps) < 1.6
    for frame in path:
        for atoms in ((0, 1, 2), (3, 4, 5)):
            endpoint = coords_r[list(atoms)]
            current = frame[list(atoms)]
            assert np.allclose(
                np.sort(np.linalg.norm(
                    current[:, None, :] - current[None, :, :], axis=2),
                        axis=None),
                np.sort(np.linalg.norm(
                    endpoint[:, None, :] - endpoint[None, :, :], axis=2),
                        axis=None), atol=1e-6)

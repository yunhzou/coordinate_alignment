import numpy as np

from rxn_core.alignment.interpolation import (
    _local_internal,
    _wrap_angle,
    _zmatrix_plan,
    internal_coordinate_interpolation,
    internal_distance_interpolation,
    proper_align_coordinates,
)


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


def test_midpoint_interpolates_local_distance_angle_and_signed_torsion():
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
    anchors, entries = _zmatrix_plan(5, bonds, bonds)
    assert len(anchors) == 3
    midpoint = np.asarray(result["frames"][1]["coords"])
    for atom, refs in entries:
        q_r = _local_internal(coords_r, atom, refs)
        q_p = _local_internal(coords_p, atom, refs)
        q_mid = _local_internal(midpoint, atom, refs)
        assert np.isclose(q_mid[0], 0.5 * (q_r[0] + q_p[0]))
        assert np.isclose(q_mid[1], 0.5 * (q_r[1] + q_p[1]))
        expected_phi = q_r[2] + 0.5 * _wrap_angle(q_p[2] - q_r[2])
        assert abs(_wrap_angle(q_mid[2] - expected_phi)) < 1e-8


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

    assert result["coordination_constraints"] == [{
        "center": 0,
        "anchors": [1, 2],
        "persistent_donors": [1, 2, 3, 4, 5, 6],
        "anchor_rule": (
            "persistent_chelating_pair_then_maximum_angular_condition"),
    }]
    determinants = []
    for frame in result["frames"]:
        xyz = np.asarray(frame["coords"])
        determinants.append(np.linalg.det(np.stack([
            xyz[1] - xyz[0], xyz[2] - xyz[0], xyz[3] - xyz[0],
        ])))
    assert min(determinants) > 0
    assert np.allclose(result["frames"][0]["coords"], coords_r)
    assert np.allclose(
        result["frames"][-1]["coords"],
        proper_align_coordinates(coords_p, coords_r))

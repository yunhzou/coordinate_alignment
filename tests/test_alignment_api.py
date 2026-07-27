from pathlib import Path

import numpy as np
import pytest

import rxn_core.alignment.api as api
import rxn_core.alignment.branch as branch
from rxn_core.alignment.index_chirality import IndexChiralityError


def test_analyze_alignment_accepts_multiplicity_not_uhf(tmp_path, monkeypatch):
    calls = []

    def fake_run_xtb(xyz_path, workdir, charge=0, uhf=0):
        calls.append((Path(workdir).name, charge, uhf))
        return ["H"], np.zeros((1, 3)), np.zeros((1, 1))

    monkeypatch.setattr(api, "run_xtb", fake_run_xtb)
    monkeypatch.setattr(api, "align_from_arrays", lambda *a, **kw: {"ok": True})

    out = api.analyze_alignment(
        tmp_path / "r.xyz", tmp_path / "p.xyz", tmp_path / "work",
        charge=-1, multiplicity=4)

    assert out == {"ok": True}
    assert calls == [("R", -1, 3), ("P", -1, 3)]


def _tetrahedral_wbo():
    wbo = np.zeros((5, 5), dtype=float)
    wbo[0, 1:] = 1.0
    wbo[1:, 0] = 1.0
    return wbo


def test_public_chirality_score_uses_wbo_topology_not_distance():
    # Every bonded neighbor is farther than the deleted 1.9-Angstrom cutoff.
    # WBO topology must still define the tetrahedral frame.
    coords = 4.0 * np.array([
        [0.0, 0.0, 0.0],
        [-0.9, -0.9, -0.9],
        [0.9, 0.9, -0.9],
        [0.9, -0.9, 0.9],
        [-0.9, 0.9, 0.9],
    ])
    odd = {index: index for index in range(5)}
    odd[3], odd[4] = odd[4], odd[3]

    assert api._index_chirality_violation_count(
        odd,
        coords,
        coords.copy(),
        _tetrahedral_wbo(),
        _tetrahedral_wbo(),
        graph_floor=0.2,
    ) > 0

    # Conversely, close coordinates do not invent bonds when WBO says there
    # is no graph edge.
    assert api._index_chirality_violation_count(
        odd,
        coords / 10.0,
        coords.copy() / 10.0,
        np.zeros((5, 5)),
        np.zeros((5, 5)),
        graph_floor=0.2,
    ) == 0


def test_public_chirality_score_rejects_nonfinite_coordinates():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [-0.9, -0.9, -0.9],
        [0.9, 0.9, -0.9],
        [0.9, -0.9, 0.9],
        [-0.9, 0.9, 0.9],
    ])
    coords[2, 1] = np.nan
    identity = {index: index for index in range(5)}

    with pytest.raises(
        IndexChiralityError,
        match="only finite values",
    ):
        api._index_chirality_violation_count(
            identity,
            coords,
            coords.copy(),
            _tetrahedral_wbo(),
            _tetrahedral_wbo(),
            graph_floor=0.2,
        )


def test_legacy_distance_chirality_helper_is_removed():
    assert not hasattr(branch, "_chirality_violations")

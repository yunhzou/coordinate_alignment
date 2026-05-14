from pathlib import Path

import numpy as np

import rxn_core.alignment.api as api


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


import json

import numpy as np

from rxn_core.cli import main


def _write_endpoint(path, elements, coordinates, edges, *, modes=False):
    wbo = np.zeros((len(elements), len(elements)))
    for left, right, value in edges:
        wbo[left, right] = wbo[right, left] = value
    payload = {
        "elements": np.asarray(elements),
        "coordinates": np.asarray(coordinates),
        "wbo": wbo,
    }
    if modes:
        payload["frequencies"] = np.asarray([-500.0])
        payload["modes"] = np.ones((1, len(elements), 3))
    np.savez(path, **payload)


def test_typed_cli_writes_rp_ts_and_self_contained_view(tmp_path):
    elements = ("C", "N", "O", "H")
    xyz = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0], [1.0, 1.0, 0.2],
    ])
    reactant, product, target = (
        tmp_path / "R.npz", tmp_path / "P.npz", tmp_path / "TS.npz")
    _write_endpoint(
        reactant, elements, xyz, ((0, 1, 1.0), (2, 3, 1.0)))
    _write_endpoint(
        product, elements, xyz + 0.1, ((0, 2, 1.0), (1, 3, 1.0)))
    _write_endpoint(
        target, elements, xyz + 0.05,
        ((0, 1, 0.5), (2, 3, 0.5), (0, 2, 0.5), (1, 3, 0.5)),
        modes=True)
    output = tmp_path / "result"

    main([
        "--reactant-npz", str(reactant),
        "--product-npz", str(product),
        "--target-npz", str(target),
        "--output", str(output),
        "--seed-count", "1",
        "--workers", "1",
    ])

    rp = json.loads((output / "rp.json").read_text())
    ts = json.loads((output / "ts_001.json").read_text())
    assert rp["schema"] == "rxn_core.rp/v2"
    assert rp["mechanisms"][0]["mapping"] == {
        "0": 0, "1": 1, "2": 2, "3": 3}
    assert ts["schema"] == "rxn_core.ts/v2"
    assert ts["mechanisms"][0]["selected"]["sources"] == [
        "product", "reactant"]
    assert (output / "mechanism_001" / "R.xyz").exists()
    assert (output / "mechanism_001" / "P_aligned.xyz").exists()
    view = (output / "view.html").read_text()
    assert "$3Dmol.createViewer" in view
    assert "<script src=" not in view

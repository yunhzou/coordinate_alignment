import importlib.util
from pathlib import Path
import sys


TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))
SPEC = importlib.util.spec_from_file_location(
    "build_retro_db_viewer", TOOLS / "build_retro_db_viewer.py")
VIEWER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VIEWER)


def test_repeated_precursor_uses_union_of_compressed_symmetry_domains():
    common = {
        "precursor_id": "ligand",
        "smiles": "C",
        "retained_atoms": [0],
        "symmetry_retained_atoms": [0, 1, 2],
        "boundary_bonds": [],
        "leftover_fragments": [],
        "complete": True,
    }
    copies = [
        dict(common, covered_target_atoms=[1], mapping=[[0, 1]],
             target_domains=[[0, [1, 2, 3]]]),
        dict(common, covered_target_atoms=[2], mapping=[[0, 2]],
             target_domains=[[0, [1, 2, 3]]]),
    ]

    group, = VIEWER._group_precursors(copies)

    assert group["multiplicity"] == 2
    assert group["covered_target_atoms"] == [1, 2]
    assert group["symmetry_target_atoms"] == [1, 2, 3]
    assert group["symmetry_retained_atoms"] == [0, 1, 2]


def test_shared_target_support_is_not_assigned_to_first_precursor(monkeypatch):
    monkeypatch.setattr(VIEWER, "mol_3d", lambda *args, **kwargs:
        ("mock sdf", [[0., 0., 0.]] * 3, ["C"] * 3))
    def source(name, atoms):
        return {"precursor_id": name, "smiles": "CCC", "retained_atoms": [0, 1],
                "covered_target_atoms": atoms, "mapping": list(enumerate(atoms)),
                "symmetry_retained_atoms": [0, 1], "target_domains": list(enumerate([[a] for a in atoms])),
                "boundary_bonds": [], "leftover_fragments": [], "complete": True}
    report = {"target_smiles": "CCC", "assemblies": [{
        "precursors": [source("A", [0, 1]), source("B", [1, 2])],
        "formed_bonds": [], "score": {}}], "construction_patterns": [],
        "scan_counts": dict(rows=2, searched=2, matched_precursors=2, fragment_candidates=2, capped=0),
        "recommendation_search_truncated": False}
    payload = VIEWER._payload(report, 20, "Overlap test")
    target = payload["assemblies"][0]["models"][-1]
    assert target["styles"][-1] == {"indices": [1], "color": "#9ca3af"}
    assert target["labels"][1]["text"] == "P1 shared: R1/R2"
    html = VIEWER._html(payload)
    assert "not validated reaction edits" in html

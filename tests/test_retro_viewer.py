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

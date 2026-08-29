from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))

from merge_retro_catalog import _assembly, _assembly_rank  # noqa: E402


def _item(precursor_id, structure_key, retained, total,
          retained_atoms=None, total_atoms=None):
    retained_atoms = retained if retained_atoms is None else retained_atoms
    total_atoms = total if total_atoms is None else total_atoms
    return {
        "precursor_id": precursor_id,
        "structure_key": structure_key,
        "retained_heavy_atoms": retained,
        "total_heavy_atoms": total,
        "retained_atom_count": retained_atoms,
        "total_atom_count": total_atoms,
        "complete": True,
        "boundary_bonds": [],
        "leftover_atom_count": 0,
    }


def test_repeated_copies_do_not_increase_unique_precursor_count():
    one = _assembly([_item("pyrazole", "pyrazole", 5, 5)], [])
    three = _assembly([
        _item("pyrazole", "pyrazole", 5, 5),
        _item("pyrazole", "pyrazole", 5, 5),
        _item("pyrazole", "pyrazole", 5, 5),
    ], [])

    assert one["score"]["unique_precursor_structures"] == 1
    assert three["score"]["unique_precursor_structures"] == 1
    assert one["score"]["set_atom_retention"] == 1.0
    assert three["score"]["set_atom_retention"] == 1.0


def test_repeated_ligands_contribute_to_set_retention():
    reaction_two = _assembly([
        _item("chloroform", "chloroform", 1, 4, 2, 5),
        _item("pyrazole", "pyrazole", 5, 5, 8, 9),
        _item("pyrazole", "pyrazole", 5, 5, 8, 9),
        _item("pyrazole", "pyrazole", 5, 5, 8, 9),
    ], [])

    assert reaction_two["score"]["unique_precursor_structures"] == 2
    assert reaction_two["score"]["set_atom_retention"] == 26 / 32


def test_retention_ranks_chloroform_above_a_large_one_carbon_source():
    chloroform = _assembly([_item("chloroform", "chloroform", 1, 4)], [])
    large_chain = _assembly([_item("chain", "chain", 1, 100)], [])

    assert _assembly_rank(chloroform) < _assembly_rank(large_chain)


def test_unique_structure_count_precedes_retention():
    one_structure = _assembly([_item("one", "one", 1, 4)], [])
    two_structures = _assembly([
        _item("left", "left", 5, 5),
        _item("right", "right", 5, 5),
    ], [])

    assert _assembly_rank(one_structure) < _assembly_rank(two_structures)

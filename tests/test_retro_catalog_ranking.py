from rxn_core.retrosynthesis.ranking import (
    assembly_rank,
    build_ranked_assembly,
    validate_atom_ownership,
)
from rxn_core.retrosynthesis.catalog_index import exact_source_copy_capacity


def _item(precursor_id, structure_key, retained, total,
          retained_atoms=None, total_atoms=None, symmetry_retained_atoms=None,
          symmetry_retained_heavy_atoms=None):
    retained_atoms = retained if retained_atoms is None else retained_atoms
    total_atoms = total if total_atoms is None else total_atoms
    symmetry_retained_atoms = (
        retained_atoms if symmetry_retained_atoms is None
        else symmetry_retained_atoms)
    symmetry_retained_heavy_atoms = (
        retained if symmetry_retained_heavy_atoms is None
        else symmetry_retained_heavy_atoms)
    return {
        "precursor_id": precursor_id,
        "structure_key": structure_key,
        "retained_heavy_atoms": retained,
        "total_heavy_atoms": total,
        "retained_atom_count": retained_atoms,
        "retained_fragments": (tuple(range(retained_atoms)),),
        "total_atom_count": total_atoms,
        "symmetry_retained_atom_count": symmetry_retained_atoms,
        "symmetry_retained_heavy_atoms": symmetry_retained_heavy_atoms,
        "complete": True,
        "chirality_violations": 0,
        "boundary_bonds": [],
        "leftover_atom_count": 0,
    }


def test_repeated_copies_do_not_increase_unique_precursor_count():
    one = build_ranked_assembly([_item("pyrazole", "pyrazole", 5, 5)], [])
    three = build_ranked_assembly([
        _item("pyrazole", "pyrazole", 5, 5),
        _item("pyrazole", "pyrazole", 5, 5),
        _item("pyrazole", "pyrazole", 5, 5),
    ], [])

    assert one["score"]["unique_precursor_structures"] == 1
    assert three["score"]["unique_precursor_structures"] == 1
    assert one["score"]["matched_fragment_count"] == 1
    assert three["score"]["matched_fragment_count"] == 3
    assert one["score"]["set_atom_retention"] == 1.0
    assert three["score"]["set_atom_retention"] == 1.0


def test_repeated_ligands_contribute_to_set_retention():
    reaction_two = build_ranked_assembly([
        _item("chloroform", "chloroform", 1, 4, 2, 5),
        _item("pyrazole", "pyrazole", 5, 5, 8, 9),
        _item("pyrazole", "pyrazole", 5, 5, 8, 9),
        _item("pyrazole", "pyrazole", 5, 5, 8, 9),
    ], [])

    assert reaction_two["score"]["unique_precursor_structures"] == 2
    assert reaction_two["score"]["set_atom_retention"] == 26 / 32


def test_retention_ranks_chloroform_above_a_large_one_carbon_source():
    chloroform = build_ranked_assembly([
        _item("chloroform", "chloroform", 1, 4)], [])
    large_chain = build_ranked_assembly([
        _item("chain", "chain", 1, 100)], [])

    assert assembly_rank(chloroform) < assembly_rank(large_chain)


def test_unique_structure_count_precedes_retention():
    one = _item("one", "one", 2, 4)
    one["retained_fragments"] = ((0,), (1,))
    one_structure = build_ranked_assembly([one], [])
    two_structures = build_ranked_assembly([
        _item("left", "left", 5, 5),
        _item("right", "right", 5, 5),
    ], [])

    assert one_structure["score"]["matched_fragment_count"] == two_structures["score"]["matched_fragment_count"]
    assert assembly_rank(one_structure) < assembly_rank(two_structures)


def test_symmetric_disassembly_precedes_equal_direct_retention():
    symmetric = build_ranked_assembly([
        _item(
            "symmetric-dimer", "symmetric-dimer", 2, 8,
            symmetry_retained_atoms=8,
            symmetry_retained_heavy_atoms=8),
    ], [])
    asymmetric = build_ranked_assembly([
        _item("asymmetric", "asymmetric", 2, 8),
    ], [])

    assert symmetric["score"]["set_atom_retention"] == 0.25
    assert symmetric["score"]["set_symmetry_atom_retention"] == 1.0
    assert assembly_rank(symmetric) < assembly_rank(asymmetric)


def test_symmetry_copy_capacity_uses_correlated_whole_fragment_action():
    swap = dict(enumerate((1, 0, 4, 5, 2, 3, 6)))
    assert exact_source_copy_capacity((0, 2, 3), (swap,))[0] == 2
    assert exact_source_copy_capacity((0, 2, 3, 6), (swap,))[0] == 1


def test_shared_target_claim_is_not_rejected():
    precursors = (
        {"covered_target_atoms": [0, 1],
         "preserved_target_bonds": ((0, 1),),
         "attachment_atoms_target": [0, 1]},
        {"covered_target_atoms": [1, 2],
         "preserved_target_bonds": ((1, 2),),
         "attachment_atoms_target": [1, 2]},
    )

    assert validate_atom_ownership(
        precursors, ((0, 1), (1, 2)), require_attachment_bonds=False,
    ) == []

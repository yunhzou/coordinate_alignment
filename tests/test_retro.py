from collections import Counter

import numpy as np
import pytest

from rxn_core import (
    WeightedGraph,
)
from rxn_core.fragment_matching import (
    FragmentCandidate,
    FragmentDetectionConfig,
    detect_fragments,
)
from rxn_core.retrosynthesis import assemble_fragment_cover


def _matrix(size, edges):
    matrix = np.zeros((size, size), dtype=float)
    for left, right, weight in edges:
        matrix[left, right] = matrix[right, left] = weight
    return matrix


def test_partial_fragment_uses_augmented_copy_when_target_cannot_match():
    precursor = WeightedGraph(
        ["C", "Br"],
        _matrix(2, [(0, 1, 1.0)]),
    )
    target = WeightedGraph(
        ["C", "C"],
        _matrix(2, [(0, 1, 1.0)]),
    )

    result = detect_fragments(
        precursor,
        target,
        source_id="carbon-bromide",
        config=FragmentDetectionConfig(
            minimum_fragment_size=1,
            branch_limit=8,
            maximum_boundary_bonds=1,
            maximum_leftover_fragments=1,
        ),
    )

    assert result.status == "matched"
    assert result.complete
    assert result.best_fragment_size == 1
    assert result.candidates
    candidate = result.candidates[0]
    assert candidate.leftover_fragments == ((1,),)
    assert candidate.boundary_bonds == ((0, 1),)
    assert candidate.copied_residual_placements == ((1, 2),)
    assert candidate.augmented_target_atom_count == 3


def test_augmented_leftover_competes_for_unused_target_atoms():
    precursor = WeightedGraph(
        ["O", "C", "O"],
        _matrix(3, [(0, 1, 2.0), (1, 2, 2.0)]),
    )
    target = WeightedGraph(
        ["O", "C", "O"],
        _matrix(3, [(0, 1, 2.0), (1, 2, 1.0)]),
    )

    result = detect_fragments(
        precursor,
        target,
        source_id="carbon-dioxide",
        config=FragmentDetectionConfig(
            minimum_fragment_size=1,
            iso_tolerance=0.5,
            branch_limit=100,
            candidate_limit=100,
        ),
    )

    assert result.status == "matched"
    assert result.complete
    assert result.best_fragment_size == 3
    assert all(candidate.covered_target_atoms == (0, 1, 2)
               for candidate in result.candidates)
    full = [candidate for candidate in result.candidates
            if candidate.covered_target_atoms == (0, 1, 2)]
    assert full
    assert full[0].retained_atoms == (0, 1, 2)
    assert full[0].retained_fragments in (
        ((0,), (1, 2)),
        ((0, 1), (2,)),
    )
    assert full[0].leftover_fragments == ()
    assert len(full[0].boundary_bonds) == 1
    assert full[0].copied_residual_placements == ()


def test_branch_cap_is_reported_as_incomplete():
    precursor = WeightedGraph(
        ["C", "C"],
        _matrix(2, [(0, 1, 1.0)]),
    )
    target = WeightedGraph(
        ["C"] * 6,
        _matrix(6, [(index, index + 1, 1.0) for index in range(5)]),
    )

    result = detect_fragments(
        precursor,
        target,
        source_id="ambiguous-CC",
        config=FragmentDetectionConfig(branch_limit=2),
    )

    assert result.status == "capped"
    assert not result.complete
    assert result.capped_seed_count > 0
    assert result.maximum_branch_count > result.branch_limit


def test_balanced_williamson_reaction_is_recovered_with_hidden_side_product():
    pytest.importorskip("rdkit")
    from rdkit import Chem
    from rxn_core.smiles import smiles_to_weighted_graph

    reactant_smiles = ("CCBr", "C[O-].[Na+]")
    product_smiles = ("CCOC", "[Na+].[Br-]")

    def composition(smiles_group):
        return Counter(
            atom.GetSymbol()
            for smiles in smiles_group
            for atom in Chem.MolFromSmiles(smiles).GetAtoms()
        )

    assert composition(reactant_smiles) == composition(product_smiles)

    target = smiles_to_weighted_graph("CCOC", expand_hydrogens=False)
    config = FragmentDetectionConfig(
        minimum_fragment_size=2,
        branch_limit=100,
        candidate_limit=100,
        maximum_boundary_bonds=1,
        maximum_leftover_fragments=1,
    )
    bromoethane = detect_fragments(
        smiles_to_weighted_graph("CCBr", expand_hydrogens=False),
        target,
        source_id="bromoethane",
        config=config,
    )
    sodium_methoxide = detect_fragments(
        smiles_to_weighted_graph("C[O-].[Na+]", expand_hydrogens=False),
        target,
        source_id="sodium-methoxide",
        config=config,
    )

    assert bromoethane.complete and bromoethane.best_fragment_size == 2
    assert sodium_methoxide.complete and sodium_methoxide.best_fragment_size == 2
    assert all(candidate.leftover_fragments == ((2,),)
               for candidate in bromoethane.candidates)
    assert all(candidate.leftover_fragments == ((2,),)
               for candidate in sodium_methoxide.candidates)

    assembly_result = assemble_fragment_cover(
        target,
        bromoethane.candidates + sodium_methoxide.candidates,
        maximum_precursors=2,
        require_attachment_bonds=True,
    )

    assert assembly_result.status == "matched"
    assert assembly_result.complete
    assert assembly_result.assemblies
    assembly = assembly_result.assemblies[0]
    assert set(assembly.precursor_ids) == {"bromoethane", "sodium-methoxide"}
    assert len(assembly.formed_bonds) == 1
    assert len(assembly.broken_bonds) == 1


def test_symmetric_cover_can_reuse_one_precursor_three_times():
    target = WeightedGraph(
        ["C", "N", "N", "N"],
        _matrix(4, [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0)]),
    )

    def candidate(precursor_id, source, image):
        return FragmentCandidate(
            source_id=precursor_id,
            mapping=((source, image),),
            retained_atoms=(source,),
            covered_target_atoms=(image,),
            leftover_fragments=(),
            boundary_bonds=(),
            attachment_atoms_source=(source,),
            attachment_atoms_target=(image,),
            copied_residual_placements=(),
            augmented_target_atom_count=4,
        )

    candidates = (
        candidate("carbon-core", 0, 0),
        candidate("same-ligand", 0, 1),
        candidate("same-ligand", 0, 2),
        candidate("same-ligand", 0, 3),
    )
    repeated = assemble_fragment_cover(
        target, candidates, maximum_precursors=4,
        allow_repeated_precursors=True)
    distinct_only = assemble_fragment_cover(
        target, candidates, maximum_precursors=4,
        allow_repeated_precursors=False)

    assert repeated.status == "matched"
    assert repeated.assemblies[0].precursor_ids.count("same-ligand") == 3
    assert len(repeated.assemblies[0].formed_bonds) == 3
    assert distinct_only.status == "no_cover"


def test_chloroform_plus_three_repeated_pyrazoles_covers_target():
    pytest.importorskip("rdkit")
    from rxn_core.smiles import smiles_to_weighted_graph

    target = smiles_to_weighted_graph(
        "N1(C(N2C=CC=N2)N3N=CC=C3)N=CC=C1",
        expand_hydrogens=True)
    config = FragmentDetectionConfig(
        minimum_fragment_size=1, branch_limit=100, candidate_limit=1000)
    core = detect_fragments(
        smiles_to_weighted_graph("ClC(Cl)Cl", expand_hydrogens=True),
        target, source_id="CHCl3", config=config)
    ligand = detect_fragments(
        smiles_to_weighted_graph("c1cn[nH]c1", expand_hydrogens=True),
        target, source_id="pyrazole", config=config)

    result = assemble_fragment_cover(
        target, core.candidates + ligand.candidates,
        maximum_precursors=4, allow_repeated_precursors=True,
        require_attachment_bonds=False)
    strict_first_pass = assemble_fragment_cover(
        target, core.candidates + ligand.candidates,
        maximum_precursors=4, allow_repeated_precursors=True,
        require_attachment_bonds=True)

    assert result.status == "matched"
    assembly = result.assemblies[0]
    assert Counter(assembly.precursor_ids) == {
        "CHCl3": 1, "pyrazole": 3}
    assert len(assembly.formed_bonds) == 3
    assert len(assembly.broken_bonds) == 6
    assert strict_first_pass.status == "no_cover"


def test_known_mcule_suzuki_precursors_cover_ortho_chlorobiphenyl():
    pytest.importorskip("rdkit")
    from rxn_core.smiles import smiles_to_weighted_graph

    # Both exact structures occur in the downloaded Mcule building-block file:
    #   BrC1=CC=CC=C1                 MCULE-5539191636
    #   C1(=CC=CC=C1Cl)B(O)O         MCULE-6011753091
    target = smiles_to_weighted_graph(
        "Clc1ccccc1-c1ccccc1", expand_hydrogens=False)
    bromobenzene = smiles_to_weighted_graph(
        "BrC1=CC=CC=C1", expand_hydrogens=False)
    chlorophenylboronic_acid = smiles_to_weighted_graph(
        "C1(=CC=CC=C1Cl)B(O)O", expand_hydrogens=False)
    config = FragmentDetectionConfig(
        branch_limit=64,
        maximum_boundary_bonds=1,
        maximum_leftover_fragments=1,
    )

    bromine_result = detect_fragments(
        bromobenzene,
        target,
        source_id="MCULE-5539191636",
        config=config,
    )
    boron_result = detect_fragments(
        chlorophenylboronic_acid,
        target,
        source_id="MCULE-6011753091",
        config=config,
    )

    assert bromine_result.complete and bromine_result.best_fragment_size == 6
    assert boron_result.complete and boron_result.best_fragment_size == 7
    assert any(candidate.leftover_fragments == ((0,),)
               for candidate in bromine_result.candidates)
    assert any(len(candidate.leftover_fragments[0]) == 3
               for candidate in boron_result.candidates)

    assembly_result = assemble_fragment_cover(
        target,
        bromine_result.candidates + boron_result.candidates,
        maximum_precursors=2,
        require_attachment_bonds=True,
    )

    assert assembly_result.status == "matched"
    assert assembly_result.complete
    assert assembly_result.assemblies
    assembly = assembly_result.assemblies[0]
    assert set(assembly.precursor_ids) == {
        "MCULE-5539191636",
        "MCULE-6011753091",
    }
    assert len(assembly.formed_bonds) == 1
    assert len(assembly.broken_bonds) == 2

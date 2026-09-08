import itertools
import numpy as np
from rxn_core import AAMProblem,MolecularEndpoint
from golden_evaluation import rank_key
from view_golden_remaining import ranker,INDICES


def test_viewer_rank_matches_existing_ranker_for_partial_and_full_mappings():
    a=MolecularEndpoint(('C','O','H'),np.zeros((3,3)),np.array([[0,2,1],[2,0,0],[1,0,0]]))
    b=MolecularEndpoint(('C','O','H'),np.zeros((3,3)),np.array([[0,1,0],[1,0,1],[0,1,0]]))
    problem=AAMProblem(a,b);score=ranker(problem)
    for size in range(4):
        for left in itertools.combinations(range(3),size):
            for right in itertools.permutations(range(3),size):
                mapping=dict(zip(left,right));expected,events=rank_key(mapping,problem)
                actual,display=score(mapping)
                assert actual==expected
                assert display==dict(broken=events['broken'],formed=events['formed'],order_changed=events['bond_order_changed'])


def test_explicit_remaining_case_scope():
    assert len(INDICES)==len(set(INDICES))==21
    assert 1739 not in INDICES and 1285 in INDICES and 1793 in INDICES


def test_packing_preserves_atom_identity_bonds_and_component_geometry():
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from view_golden_remaining import packed
    mol=Chem.AddHs(Chem.MolFromSmiles('CCO.O.N'))
    rdDepictor.Compute2DCoords(mol)
    for atom in mol.GetAtoms():atom.SetIntProp('original_index',atom.GetIdx())
    moved=packed(mol,True)
    assert Chem.MolToSmiles(moved)==Chem.MolToSmiles(mol)
    assert [a.GetIntProp('original_index') for a in moved.GetAtoms()]==list(range(mol.GetNumAtoms()))
    for atoms in Chem.GetMolFrags(mol):
        for a,b in itertools.combinations(atoms,2):
            before=mol.GetConformer().GetAtomPosition(a)-mol.GetConformer().GetAtomPosition(b)
            after=moved.GetConformer().GetAtomPosition(a)-moved.GetConformer().GetAtomPosition(b)
            assert abs(before.Length()-after.Length())<1e-9

"""A display must retain the archived atom identities after layout/H hiding."""
import sys
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from golden_evaluation import prepare
from view_golden_mapping import molecules, svg


def test_display_uses_exact_benchmark_order_with_disconnected_sources():
    reaction = '[CH3:7][Br:8].[OH2:9].[K+]>>[CH3:7][OH:9]'
    problem, _, _ = prepare(reaction)
    endpoints = molecules(reaction, problem)
    for mol, endpoint in zip(endpoints, (problem.reactant, problem.product)):
        assert tuple(a.GetSymbol() for a in mol.GetAtoms()) == endpoint.elements
        for part in Chem.GetMolFrags(mol, asMols=True):
            for atom in part.GetAtoms():
                i = atom.GetIntProp('original_index')
                assert mol.GetAtomWithIdx(i).GetSymbol() == atom.GetSymbol()


def test_svg_hiding_h_preserves_original_indices_and_has_no_external_assets():
    reaction = '[CH3:7][OH:9]>>[CH2:7]=[O:9]'
    problem, _, _ = prepare(reaction)
    reactant, _ = molecules(reaction, problem)
    for hydrogens in (False, True):
        drawing = svg(reactant, {0:'#009e73'}, hydrogens, 400, 300)
        assert '<svg' in drawing and '</svg>' in drawing
        assert '<image' not in drawing and '<script' not in drawing
    assert reactant.GetNumAtoms() == problem.source_atom_count

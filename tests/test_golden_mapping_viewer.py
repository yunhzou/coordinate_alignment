"""A display must retain the archived atom identities after layout/H hiding."""
import sys
from pathlib import Path

from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from golden_evaluation import prepare
from view_golden_mapping import molecules, svg, mapping_fragments


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


def test_product_hit_targets_keep_true_source_identity_not_same_number():
    reaction = '[CH3:7][OH:9]>>[CH2:7]=[O:9]'
    problem, _, _ = prepare(reaction)
    _, product = molecules(reaction, problem)
    drawing = svg(product, {}, False, 400, 300, mapping={19:0, 27:1},
                  owner={19:4,27:5}, side='P')
    assert 'data-r="19" data-p="0"' in drawing
    assert 'R5:r19 → P:p0' in drawing
    assert 'data-r="27" data-p="1"' in drawing


def test_separate_fragments_in_one_source_are_not_merged_by_molecule():
    problem, _, _ = prepare('[CH3:7][OH:9]>>[CH2:7]=[O:9]')
    record = dict(mapping={0:0,1:1}, context={}, steps=[dict(fragment=[0]), dict(fragment=[1])])
    fragments = mapping_fragments(record, problem)
    assert [f['source'] for f in fragments] == [[0],[1]]
    assert [f['label'] for f in fragments] == ['F1','F2']
    reference = mapping_fragments(dict(mapping={0:0,1:1}, context=None), problem)
    assert all(f['kind'] == 'Reference conserved region' for f in reference)

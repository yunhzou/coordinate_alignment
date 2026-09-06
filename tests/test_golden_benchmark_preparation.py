import importlib.util
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdChemReactions

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('golden_prep', ROOT/'bench/prepare_golden_benchmark.py')
PREP = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PREP)


def audit(smiles):
    reaction = rdChemReactions.ReactionFromSmarts(smiles, useSmiles=True)
    return PREP.audit_block(rdChemReactions.ReactionToRxnBlock(reaction), 7)


def test_reference_labels_never_enter_algorithm_input():
    row = audit('[CH3:1][Br:2].[OH2:3]>>[CH3:1][OH:3].[BrH:2]')
    assert row['status'] == 'parsed'
    assert row['balanced_explicit_H']
    assert row['shared_map_labels'] == 3
    for side in row['input_reaction'].split('>>'):
        assert not any(atom.GetAtomMapNum() for atom in Chem.MolFromSmiles(side).GetAtoms())
    assert ':1]' in row['mapped_reaction']


def test_hydrogen_imbalance_is_reported_without_repair_or_exclusion():
    row = audit('[CH3:1][OH:2]>>[CH2:1]=[O:2]')
    assert row['status'] == 'parsed'
    assert row['balanced_heavy']
    assert not row['balanced_explicit_H']
    assert row['reference_mapped_H'] == [0, 0]


def test_parse_failure_is_explicit_not_a_silent_skip():
    row = PREP.audit_block('invalid', 9)
    assert row['index'] == 9
    assert row['status'] == 'parse_or_sanitize_error'

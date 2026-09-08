import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from golden_competitors import signatures


def test_label_numbers_and_atom_order_do_not_matter():
    a = '[CH3:1][OH:2]>>[CH2:1]=[O:2]'
    b = '[OH:42][CH3:71]>>[O:42]=[CH2:71]'
    assert signatures(a) == signatures(b)


def test_endpoint_symmetry_is_not_a_mapping_error():
    a = '[CH3:1][CH2:2][CH3:3]>>[CH3:1][CH2:2][CH3:3]'
    b = '[CH3:1][CH2:2][CH3:3]>>[CH3:3][CH2:2][CH3:1]'
    assert signatures(a) == signatures(b)


def test_wrong_mapping_detected_with_same_endpoint_chemistry():
    a = '[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:3]'
    b = '[CH3:1][CH2:2][OH:3]>>[CH3:2][CH:1]=[O:3]'
    assert signatures(a)[0] == signatures(b)[0]
    assert signatures(a)[1] != signatures(b)[1]


def test_unmatched_atoms_are_retained():
    a = '[CH3:1][OH:2].O>>[CH2:1]=[O:2].O'
    b = '[CH3:1][OH:2].[OH2:3]>>[CH2:1]=[O:2].[OH2:3]'
    assert signatures(a)[0] == signatures(b)[0]
    assert signatures(a)[1] != signatures(b)[1]


def test_duplicate_heavy_labels_are_not_silently_overwritten():
    with pytest.raises(ValueError, match='Duplicate'):
        signatures('[CH3:1][OH:1]>>[CH2:1]=[O:2]')


def test_changed_chemistry_is_detected():
    assert signatures('CO>>C=O')[0] != signatures('CO>>CO')[0]


def test_agent_field_preserves_all_original_input_components():
    assert signatures('[CH3:1][OH:2].O>>[CH2:1]=[O:2]') == signatures(
        '[CH3:1][OH:2]>O>[CH2:1]=[O:2]')

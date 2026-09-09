import importlib.util
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'bench'))
from elementary_feasibility import valid_mapping


def test_full_element_mapping():
    assert valid_mapping(['C','H','O'], ['O','C','H'], [(0,1),(1,2),(2,0)])


def test_partial_duplicate_and_element_changes_rejected():
    assert not valid_mapping(['C','H'], ['C','H'], [(0,0)])
    assert not valid_mapping(['C','H'], ['C','H'], [(0,0),(1,0)])
    assert not valid_mapping(['C','H'], ['C','H'], [(0,1),(1,0)])

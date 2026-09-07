import sys
from pathlib import Path

import pynauty

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from golden_direction_probe import original_mapping
from golden_evaluation import prepare, project, colored_graph


def test_reverse_mapping_preserves_original_indices_and_unmatched_atoms():
    pairs = [(0,57),(1,76),(2,13)]
    assert original_mapping(pairs,'P_to_R') == {57:0,76:1,13:2}
    assert original_mapping([(57,0),(76,1),(13,2)],'R_to_P') == {57:0,76:1,13:2}


def test_swapped_reference_returns_identical_original_chemical_certificate():
    _, features, reference = prepare('[CH3:1][Br:2].[OH2:3]>>[CH3:1][OH:3]')
    backwards = [(p,r) for r,p in reference.items()]
    restored = original_mapping(backwards,'P_to_R')
    assert pynauty.certificate(colored_graph(features,project(restored,features))) == pynauty.certificate(colored_graph(features,project(reference,features)))

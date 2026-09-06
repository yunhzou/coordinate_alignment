"""Reference-data regressions; no example-specific rules enter the matcher."""
import json
from pathlib import Path

from rdkit import Chem


def test_t05_corrected_nitrile_regioisomer_preserves_ring_attachment():
    case = json.loads((Path(__file__).parents[1] /
                       'docs/example_runs/t05_ground_truth.json').read_text())
    target = Chem.MolFromSmiles(case['target_smiles'])
    source = next(r for r in case['reactants'] if r['id'] == 'INVENTORY-001198')
    assert target.HasSubstructMatch(Chem.MolFromSmiles(source['smiles']))
    # The previously misread nitrile position cannot keep the whole R3 skeleton.
    assert not target.HasSubstructMatch(Chem.MolFromSmiles('N#Cc1ccc2sc(N)nc2c1'))

"""The experimental proposal may reorder seeds, never constrain their pairs."""
import sys
from pathlib import Path
import numpy as np
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from slap_guided_pilot import guided_order


def test_guided_order_retains_all_explicit_atoms_and_is_reproducible():
    r=np.zeros((4,4));r[0,1]=r[1,0]=1;r[2,3]=r[3,2]=1
    p=r.copy();p[0,1]=p[1,0]=0
    problem=AAMProblem(*(MolecularEndpoint(('C','H','C','H'),np.zeros((4,3)),w) for w in (r,p)))
    identity=dict(enumerate(range(4)))
    order=guided_order(problem,identity)
    assert sorted(order)==list(range(4))
    assert order==guided_order(problem,identity)
    assert set(order[:2])=={2,3}
    assert identity==dict(enumerate(range(4)))

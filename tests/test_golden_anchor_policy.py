import numpy as np
from rxn_core import AAMProblem,MolecularEndpoint
from golden_anchor_policy import proposals


def problem(metadata=None):
    elements=('O','C','H','H')
    wbo=np.array([[0,1,1,0],[1,0,0,1],[1,0,0,0],[0,1,0,0]],dtype=float)
    endpoint=MolecularEndpoint(elements,np.zeros((4,3)),wbo,metadata=metadata)
    return AAMProblem(endpoint,endpoint)


def test_proposals_are_reference_blind_reproducible_and_element_compatible():
    p=problem()
    selected=proposals(p,3)
    assert selected==proposals(p,3)
    assert selected==proposals(problem({'reference_mapping':[(0,3)]}),3)
    for policy,pairs in selected.items():
        assert len(pairs)<=3 and len(set(pairs))==len(pairs)
        for a,b in pairs:
            if policy=='random_single_cut':assert p.reactant.wbo[a,b]>.2
            else:assert p.reactant.elements[a]==p.product.elements[b]


def test_hydrogens_are_eligible_and_budget_does_not_duplicate_pairs():
    p=problem();selected=proposals(p,100)
    expected={(a,b) for a,e in enumerate(p.reactant.elements)
              for b,f in enumerate(p.product.elements) if e==f}
    for policy in ('random_anchor','stable_anchor','change_anchor'):
        assert set(selected[policy])==expected
    assert (2,3) in expected

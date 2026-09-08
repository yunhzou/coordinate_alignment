from itertools import permutations
import numpy as np

from rxn_core import AAMProblem,AAMSearchConfig,search_aam
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import minimize_events,bond_events
from golden_evaluation import prepare
from view_golden_remaining import ranker


def endpoint(elements,edges):
    w=np.zeros((len(elements),len(elements)))
    for a,b,v in edges:w[a,b]=w[b,a]=v
    return MolecularEndpoint(tuple(elements),np.zeros((len(elements),3)),w)


def test_tolerant_symmetry_minimum_matches_tiny_exhaustive_oracle():
    mol=endpoint('CCCHHH',[(0,1,2),(1,2,1),(0,2,1),(0,3,1),(1,4,1),(2,5,1)])
    problem=AAMProblem(mol,mol)
    aam=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=1))
    path=next(aam.graph.paths())
    scores=[bond_events(problem,dict(enumerate((*p,*(i+3 for i in p)))))['total']
            for p in permutations(range(3))]
    assert set(scores)=={0,2}
    result=minimize_events(path,problem,seconds=5)
    assert result['method']=='bounded_symbolic'
    assert result['optimal'] and result['upper_bound']==min(scores)


def test_partial_mapping_both_orientations_match_exhaustive_oracle():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOOO',[(0,1,1),(2,3,1)]))
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1)).graph.paths())
    for reverse in (False,True):
        oracle=min(bond_events(problem,dict(enumerate(p)),reverse=reverse)['total']
                   for p in permutations(range(4),2))
        result=minimize_events(path,problem,reverse=reverse,seconds=5)
        assert result['optimal'] and result['upper_bound']==oracle


def test_explicit_h_symmetry_certified_without_solver():
    problem,_,_=prepare('CC(C)(C)C>>CC(C)(C)C')
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10)).graph.paths())
    result=minimize_events(path,problem)
    assert result['method']=='invariance_certificate' and result['solver_queries']==0
    assert result['upper_bound']==0 and len(result['mapping'])==17


def test_event_score_matches_existing_ranker_including_unmatched_atoms():
    problem,_,_=prepare('CCO.O>>CC=O.O')
    mapping={i:i for i in range(min(problem.source_atom_count,problem.target_atom_count)-1)}
    for reverse in (False,True):
        actual=problem if not reverse else AAMProblem(problem.product,problem.reactant)
        oriented=mapping if not reverse else {v:k for k,v in mapping.items()}
        expected=ranker(actual)(oriented)[0][2]
        assert bond_events(problem,mapping,reverse=reverse)['total']==expected

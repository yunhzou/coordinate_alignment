"""Forensic regression: an assignment domain need not equal an automorphism orbit."""
import sys
from pathlib import Path
import numpy as np
import pynauty

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from golden_evaluation import colored_graph,symbolic_path_query,evaluate
from rxn_core import AAMProblem,AAMSearchConfig,search_aam
from rxn_core.domain import MolecularEndpoint


def test_saved_singleton_domain_can_cross_conditioned_automorphism_orbits():
    source=MolecularEndpoint(('O','O'),np.zeros((2,3)),np.zeros((2,2)))
    bonds=np.zeros((4,4));bonds[0,1]=bonds[1,0]=bonds[2,3]=bonds[3,2]=1
    target=MolecularEndpoint(('O',)*4,np.zeros((4,3)),bonds)
    result=search_aam(AAMProblem(source,target),AAMSearchConfig())
    path=next(result.graph.paths());edge=result.graph.transitions[path.transitions[-1]]
    seed=edge.match['fragment'][0];old=path.mapping[seed]
    domain=edge.match['symmetry']['blocks'][0]['p_atoms']
    assert set(domain)=={1,2,3} and old==1
    assert all(g[old]==old for g in edge.match['symmetry']['automorph_generators'])
    assert 2 in domain and not result.graph.capped
    trial={**path.mapping,seed:2}
    assert len(set(trial.values()))==2
    assert bonds[tuple(path.mapping.values())]==1
    assert bonds[tuple(trial.values())]==0
    def feature(n,edges):
        return dict(heavy=list(range(n)),colors=[('O',0,0,0,'')]*n,
                    bonds=[(a,b,(1.0,'STEREONONE')) for a,b in edges])
    score=evaluate(result,[feature(2,[]),feature(4,[(0,1),(2,3)])],trial,
                   reference_side='source',seconds=10)
    assert not score['representative_recovery']
    assert score['reference_recovery']=='recovered'
    assert score['witness_actions']['scope']=='full_explicit'


def test_explicit_singleton_domain_query_can_certify_alternative():
    def feature(n,edges):
        return dict(heavy=list(range(n)),colors=[('O',0,0,0,'')]*n,
                    bonds=[(a,b,(1.0,'STEREONONE')) for a,b in edges])
    features=[feature(2,[]),feature(4,[(0,1),(2,3)])]
    reference={0:0,1:2}
    status,witness=symbolic_path_query({0:0,1:1},[],features,reference,5000,
        independent_singleton_domains={1:(1,2,3)})
    assert status=='recovered'
    assert pynauty.certificate(colored_graph(features,dict(witness['heavy_mapping'])))==pynauty.certificate(colored_graph(features,reference))

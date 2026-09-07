import importlib.util
from itertools import product
from pathlib import Path
import pynauty

spec=importlib.util.spec_from_file_location('golden_evaluation',Path(__file__).parents[1]/'bench/golden_evaluation.py')
E=importlib.util.module_from_spec(spec);spec.loader.exec_module(E)


def feature(n, distinct=False):
    return dict(heavy=list(range(n)),colors=[('C',i if distinct else 0) for i in range(n)],bonds=[])


def test_joint_chemical_equivalence_is_not_independent_atom_membership():
    f=feature(4);f['bonds']=[(0,1,(1.,'')),(1,2,(1.,'')),(2,3,(1.,''))]
    original=pynauty.certificate(E.colored_graph([f,f],dict(enumerate(range(4)))))
    reversed_=pynauty.certificate(E.colored_graph([f,f],dict(enumerate([3,2,1,0]))))
    wrong=pynauty.certificate(E.colored_graph([f,f],dict(enumerate([0,2,1,3]))))
    assert original==reversed_ and original!=wrong


def test_symbolic_conditioned_group_product_matches_tiny_exhaustive_oracle():
    f=feature(3,True)
    groups=[((1,0,2),),((0,2,1),)]
    possible=set()
    for first,second in product([(0,1,2),(1,0,2)],[(0,1,2),(0,2,1)]):
        possible.add(tuple(first[second[i]] for i in range(3)))
    from itertools import permutations
    for values in permutations(range(3)):
        status,_=E.symbolic_path_query(dict(enumerate(range(3))),groups,[f,f],dict(enumerate(values)),5000)
        assert (status=='recovered') == (values in possible)


def test_symbolic_query_preserves_unmatched_domain_and_free_hydrogen_projection():
    f=feature(3,True)
    status,_=E.symbolic_path_query({0:0,1:1},[((1,0,2),)],[f,f],{0:1,1:0},5000)
    assert status=='recovered'
    status,_=E.symbolic_path_query({0:0,1:1},[((1,0,2),)],[f,f],{0:1,2:0},5000)
    assert status=='not_recovered'


def test_inputs_are_unmapped_explicit_H_and_reference_indices_are_separate():
    problem,features,reference=E.prepare('[CH3:7][Br:2].[OH2:5]>>[CH3:7][OH:5].[BrH:2]')
    assert problem.balanced
    assert problem.source_atom_count==8
    assert len(E.project(reference,features))==3
    assert all(':' not in e.metadata['unmapped_smiles'] for e in (problem.reactant,problem.product))


def test_reference_recovery_queries_real_compressed_aam_not_only_representative():
    import numpy as np
    from rxn_core import AAMProblem,AAMSearchConfig,search_aam
    from rxn_core.domain import MolecularEndpoint
    bonds=np.array([[0.,1.,0.],[1.,0.,2.],[0.,2.,0.]])
    endpoint=MolecularEndpoint(('C','C','C'),np.zeros((3,3)),bonds)
    result=search_aam(AAMProblem(endpoint,endpoint),AAMSearchConfig(seed_count=1,cut_floor=10))
    mapping=dict(result.graph.states[result.graph.terminals[0]].mapping)
    reference={a:(2,1,0)[b] for a,b in mapping.items()}
    f=feature(3,True)
    evaluated=E.evaluate(result,[f,f],reference,seconds=10)
    assert not evaluated['representative_recovery']
    assert evaluated['reference_recovery']=='recovered'
    assert evaluated['symbolic_queries']>0


def test_H_only_witness_differences_do_not_repeat_heavy_queries(monkeypatch):
    from types import SimpleNamespace as S
    import numpy as np
    from rxn_core import AAMProblem
    from rxn_core.domain import MolecularEndpoint
    endpoint=MolecularEndpoint(('C','C','H','H'),np.zeros((4,3)),np.zeros((4,4)))
    mappings=[{0:0,1:1,2:2,3:3},{0:0,1:1,2:3,3:2}]
    graph=S(terminals=[0,1],states=[S(mapping=m) for m in mappings],capped=False,
        paths=lambda:iter([S(mapping=m,terminal=i,transitions=[i]) for i,m in enumerate(mappings)]),
        fragment_placement=lambda edge:S(target_generators=[S(images=(1,0,2,3)),S(images=(0,1,3,2))]))
    calls=[]
    def query(*args):calls.append(args);return 'unknown',None
    monkeypatch.setattr(E,'symbolic_path_query',query)
    f=feature(2,True)
    result=E.evaluate(S(graph=graph,problem=AAMProblem(endpoint,endpoint)),[f,f],{0:1,1:0})
    assert result['reference_recovery']=='unknown'
    assert len(calls)==result['symbolic_queries']==1

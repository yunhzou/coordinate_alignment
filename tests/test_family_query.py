from itertools import permutations
import numpy as np

from rxn_core import AAMProblem,AAMSearchConfig,search_aam,query_path
from rxn_core.domain import MolecularEndpoint
from rxn_core.search_graph import AAMSearchGraph,SearchContext,SearchState,FragmentTransition,SearchStop


def endpoint(elements,edges):
    w=np.zeros((len(elements),len(elements)))
    for a,b in edges:w[a,b]=w[b,a]=1
    return MolecularEndpoint(tuple(elements),np.zeros((len(elements),3)),w)


def synthetic_path(blocks,automorph_blocks=()):
    graph=AAMSearchGraph((SearchContext((0,1),(0,1,2),(0,1)),),(0,),
        (SearchState(0,0,(),(),()),SearchState(1,0,((0,0),(1,1)),((0,1),(1,1)),())),
        (FragmentTransition(0,0,1,0,(1,0),dict(fragment=[0,1],symmetry=dict(
            witness={0:0,1:1},blocks=blocks,exact_fixed=[],automorph_blocks=automorph_blocks,
            automorph_generators=[]))),),(SearchStop(1,'objective_met'),))
    return next(graph.paths())


def test_open_pools_keep_same_pair_and_different_pair_assignments():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOOO',[(0,1),(2,3)]))
    result=search_aam(problem,AAMSearchConfig(seed_count=1))
    path=next(result.graph.paths())
    # Exhaustive oracle only for this tiny unit test, never production search.
    for image in permutations(range(4),2):
        status,witness=query_path(path,problem,dict(enumerate(image)),source_atoms=(0,1),timeout_ms=2000)
        assert status=='recovered'
        assert dict(witness['mapping'])==dict(enumerate(image))


def test_correlated_fragment_is_not_independent_atom_membership():
    molecule=endpoint('CCCC',[(0,1),(1,2),(2,3)])
    problem=AAMProblem(molecule,molecule)
    result=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=.5))
    paths=list(result.graph.paths())
    for image in [(0,1,2,3),(3,2,1,0)]:
        assert any(query_path(path,problem,dict(enumerate(image)),source_atoms=tuple(range(4)),
            target_generators=((3,2,1,0),),timeout_ms=2000)[0]=='recovered' for path in paths)
    assert all(query_path(path,problem,{0:0,1:2,2:1,3:3},source_atoms=tuple(range(4)),
        target_generators=((3,2,1,0),),timeout_ms=2000)[0]=='not_recovered' for path in paths)


def test_full_hydrogen_mapping_is_returned_and_one_to_one():
    molecule=endpoint('OHOH',[(0,1),(2,3)])
    problem=AAMProblem(molecule,molecule)
    result=search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10))
    status,witness=query_path(next(result.graph.paths()),problem,{0:2,2:0},source_atoms=(0,2),timeout_ms=3000)
    assert status=='recovered'
    mapping=dict(witness['mapping'])
    assert len(mapping)==len(set(mapping.values()))==4
    assert molecule.wbo[mapping[0],mapping[1]]==molecule.wbo[mapping[2],mapping[3]]==1


def test_one_target_cannot_supply_two_reference_atoms():
    molecule=endpoint('OO',[]);problem=AAMProblem(molecule,molecule)
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1)).graph.paths())
    assert query_path(path,problem,{0:0,1:0},source_atoms=(0,1),timeout_ms=2000)[0]=='not_recovered'


def test_pool_choices_do_not_move_fixed_fragment_slots():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOO',[]))
    path=synthetic_path([dict(r_atoms=[1],p_atoms=[0,1,2])])
    assert query_path(path,problem,{0:0,1:2},source_atoms=(0,1))[0]=='recovered'
    assert query_path(path,problem,{0:2,1:0},source_atoms=(0,1))[0]=='not_recovered'


def test_automorphism_summary_is_not_an_independent_domain():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOO',[]))
    path=synthetic_path([],automorph_blocks=[dict(r_atoms=[0,1],p_atoms=[0,1,2])])
    assert query_path(path,problem,{0:0,1:2},source_atoms=(0,1))[0]=='not_recovered'


def test_heavy_projection_positive_cannot_replace_full_H_check():
    problem=AAMProblem(endpoint('OH',[(0,1)]),endpoint('OHO',[(0,1)]))
    path=synthetic_path([dict(r_atoms=[0],p_atoms=[0,2])])
    assert query_path(path,problem,{0:2},source_atoms=(0,),projected_atoms=(0,))[0]=='recovered'
    assert query_path(path,problem,{0:2},source_atoms=(0,))[0]=='not_recovered'

from itertools import permutations

from rxn_core import AAMProblem,AAMSearchConfig,search_aam
from rxn_core.pattern_collection import PatternEquivalence,extract_path_patterns
from test_family_scoring import endpoint
from golden_evaluation import prepare


def test_exact_classes_match_small_exhaustive_oracle():
    mol=endpoint('CCCHHH',[(0,1,2),(1,2,1),(0,2,1),(0,3,1),(1,4,1),(2,5,1)])
    problem=AAMProblem(mol,mol)
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=1)).graph.paths())
    equivalence=PatternEquivalence(problem)
    expected={equivalence.key(dict(enumerate((*p,*(i+3 for i in p))))) for p in permutations(range(3))}
    result=extract_path_patterns(path,problem,equivalence,seconds=5)
    assert result['complete'],result
    assert {p['key'] for p in result['patterns']}==expected
    assert len(expected)==2
    resumed=extract_path_patterns(path,problem,equivalence,previous=result['patterns'],seconds=5)
    assert resumed['complete'] and len(resumed['patterns'])==2


def test_hydrogen_symmetry_is_certified_without_enumeration():
    problem,_,_=prepare('CC(C)(C)C>>CC(C)(C)C')
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1,cut_floor=10)).graph.paths())
    result=extract_path_patterns(path,problem,PatternEquivalence(problem))
    assert result['complete'] and result['solver_queries']==0
    assert len(result['patterns'])==1


def test_partial_mapping_classes():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OOOO',[(0,1,1),(2,3,1)]))
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1)).graph.paths())
    eq=PatternEquivalence(problem)
    expected={eq.key(dict(enumerate(p))) for p in permutations(range(4),2)}
    result=extract_path_patterns(path,problem,eq,seconds=5)
    assert result['complete']
    assert {p['key'] for p in result['patterns']}==expected


def test_same_score_patterns_and_zero_budget_resume():
    problem=AAMProblem(endpoint('OO',[]),endpoint('OO',[]))
    path=next(search_aam(problem,AAMSearchConfig(seed_count=1)).graph.paths())
    # Supplied endpoint annotations distinguish the two physical assignments,
    # despite identical bond event counts. No charge/resonance normalization.
    eq=PatternEquivalence(problem,atom_tags=({0:'a',1:'b'},{0:'a',1:'b'}))
    partial=extract_path_patterns(path,problem,eq,seconds=0)
    assert not partial['complete'] and len(partial['patterns'])==1
    result=extract_path_patterns(path,problem,eq,previous=partial['patterns'],seconds=5)
    assert result['complete'] and len(result['patterns'])==2
    assert {p['events']['total'] for p in result['patterns']}=={0}


def test_twin_quotient_matches_full_graph_certificate_for_all_small_maps():
    import z3
    mol=endpoint('CCHHHH',[(0,1,1),(0,2,1),(0,3,1),(1,4,1),(1,5,1)])
    eq=PatternEquivalence(AAMProblem(mol,mol))
    assert sorted(map(len,eq.twins[0]))==[1,1,2,2]
    identity=dict(enumerate(range(6)))
    for carbons in permutations(range(2)):
        for hydrogens in permutations(range(2,6)):
            mapping=dict(enumerate((*carbons,*hydrogens)))
            values=[(p,frozenset((p,))) for p in mapping.values()]
            solver=z3.Solver();solver.add(eq.orbit_membership(values,identity))
            assert (solver.check()==z3.sat)==(eq.key(mapping)==eq.key(identity))


def test_twin_quotient_preserves_unmapped_atom_patterns():
    import z3
    from itertools import combinations
    mol=endpoint('OOOO',[(0,1,1),(2,3,1)])
    eq=PatternEquivalence(AAMProblem(mol,mol));known={0:0,1:2}
    for source in combinations(range(4),2):
        for target in permutations(range(4),2):
            mapping=dict(zip(source,target))
            values=[(mapping.get(r,4),frozenset((mapping.get(r,4),))) for r in range(4)]
            solver=z3.Solver();solver.add(eq.orbit_membership(values,known))
            assert (solver.check()==z3.sat)==(eq.key(mapping)==eq.key(known))

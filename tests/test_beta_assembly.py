"""Full correlated supplier assembly and pattern diversity, not first-hit coverage."""
from dataclasses import replace
from itertools import product

from test_beta_retro import FakeBank, block
from rxn_core.retrosynthesis.beta import recommend_big_blocks
from rxn_core.retrosynthesis.beta_assembly import (
    assemble_supplier_copies, completed_assembly_rank, rank_complete_assemblies,
    assembly_metrics,
)


def test_supplier_join_equals_exhaustive_correlated_covers_and_fragment_order():
    first=(block('a',(0,1),refined=True),block('a',(2,3),refined=True))
    second=(block('b',(0,2),refined=True),block('b',(1,2,3),refined=True),
            block('b',(0,1),refined=True))
    second=(*second,replace(second[1],candidate=replace(second[1].candidate,
        retained_fragments=((0,),(1,2)))))
    expected=[pair for pair in product(first,second)
              if set().union(*(p.covered_atoms for p in pair))==set(range(4))]
    actual=list(assemble_supplier_copies((first,second),range(4)))
    assert {tuple(p.key for p in pair) for pair in actual}=={
        tuple(p.key for p in pair) for pair in expected}
    costs=[sum(p.fragment_count for p in pair) for pair in actual]
    assert costs==sorted(costs)


def test_possible_atom_union_does_not_substitute_for_compatible_copy_assignment():
    pool=(block('r',(0,1),refined=True),block('r',(2,3),refined=True))
    assert list(assemble_supplier_copies((pool,),range(4)))==[]
    assert list(assemble_supplier_copies((pool,pool),range(4)))


def test_beta_collects_distinct_patterns_after_first_complete_cover():
    class Bank(FakeBank):
        def query(self,region):
            return (block('a',range(4)),block('b',range(4)))
        def refine(self,selected):
            self.events.append(selected.source_id)
            candidate=selected.candidate
            if selected.source_id=='b':
                candidate=replace(candidate,retained_fragments=((0,1),(2,3)))
            return (replace(selected,candidate=candidate,refined=True),)
    bank=Bank()
    result=recommend_big_blocks(bank,recommendations=2,pattern_limit=2)
    assert len(result.recommendations)==2
    assert bank.events==['a','b']
    assert [sum(p.fragment_count for p in r.placements) for r in result.recommendations]==[1,2]


def test_final_rank_prefers_retention_over_fewer_species():
    target = FakeBank().target
    efficient = (block('a',(0,1),refined=True), block('b',(2,3),refined=True))
    wasteful = tuple(replace(p,candidate=replace(p.candidate,source_id='large',
        leftover_fragments=((10,11,12,13),))) for p in efficient)
    assert completed_assembly_rank(efficient,target) < completed_assembly_rank(wasteful,target)


def test_final_rank_uses_structural_operations_and_rejects_partial():
    import pytest
    target = FakeBank().target
    clean = (block('a',range(4),refined=True),)
    cut = (replace(clean[0],candidate=replace(clean[0].candidate,boundary_bonds=((0,5),))),)
    assert completed_assembly_rank(clean,target) < completed_assembly_rank(cut,target)
    with pytest.raises(ValueError,match='complete'):
        completed_assembly_rank((block('a',(0,1)),),target)


def test_repeated_overlapping_copies_do_not_inflate_retention():
    target = FakeBank().target
    p = block('a',range(4),refined=True)
    one = assembly_metrics((p,),target)
    two = assembly_metrics((p,p),target)
    assert two['retention'] == one['retention']/2
    assert two['species'] == one['species'] == 1


def test_display_order_is_final_rank_even_with_pattern_diversity():
    from rxn_core.retrosynthesis.beta import BetaRecommendation
    target = FakeBank().target
    a = block('a',range(4),refined=True)
    b = replace(a,candidate=replace(a.candidate,source_id='b'))
    c = replace(a,candidate=replace(a.candidate,source_id='c',retained_fragments=((0,1),(2,3))))
    answers = [BetaRecommendation((p,),()) for p in (c,b,a)]
    result = rank_complete_assemblies(answers,target,3,2)
    assert len(result)==3
    ranks = [completed_assembly_rank(r.placements,target) for r in result]
    assert ranks == sorted(ranks)

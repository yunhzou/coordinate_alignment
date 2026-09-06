"""Full correlated supplier assembly and pattern diversity, not first-hit coverage."""
from dataclasses import replace
from itertools import product

from test_beta_retro import FakeBank, block
from rxn_core.retrosynthesis.beta import recommend_big_blocks
from rxn_core.retrosynthesis.beta_assembly import assemble_supplier_copies


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

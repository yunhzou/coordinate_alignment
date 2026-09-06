"""Full correlated supplier assembly and pattern diversity, not first-hit coverage."""
from dataclasses import replace
from itertools import product

from test_beta_retro import FakeBank, block
from rxn_core.retrosynthesis.beta import recommend_big_blocks
from rxn_core.retrosynthesis.beta_assembly import (
    assemble_supplier_copies, rank_complete_assemblies,
    assembly_metrics, dominates, pareto_assembly_ranks, assembly_key,
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
    assert dominates(assembly_metrics(efficient,target), assembly_metrics(wasteful,target))


def test_final_rank_uses_structural_operations_and_rejects_partial():
    import pytest
    target = FakeBank().target
    clean = (block('a',range(4),refined=True),)
    cut = (replace(clean[0],candidate=replace(clean[0].candidate,boundary_bonds=((0,5),))),)
    assert dominates(assembly_metrics(clean,target), assembly_metrics(cut,target))
    with pytest.raises(ValueError,match='complete'):
        from rxn_core.retrosynthesis.beta import BetaRecommendation
        pareto_assembly_ranks((BetaRecommendation((block('a',(0,1)),),(2,3)),),target)


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
    ranks = pareto_assembly_ranks(answers,target)
    assert [ranks[assembly_key(r)] for r in result] == sorted(ranks.values())
    assert {r[0] for r in ranks.values()} == {1}  # partitions do not affect rank


def test_pareto_layers_equal_exhaustive_dominance_and_preserve_tradeoffs():
    from rxn_core.retrosynthesis.beta import BetaRecommendation
    target = FakeBank().target
    answers = []
    for waste in range(5):
        for cuts in range(5):
            for duplicate in range(2):
                p = block(f'r{waste}-{cuts}-{duplicate}',range(4),refined=True)
                c = replace(p.candidate,leftover_fragments=(tuple(range(10,10+waste)),),
                            boundary_bonds=tuple((0,10+i) for i in range(cuts)))
                answers.append(BetaRecommendation((replace(p,candidate=c),),()))
    metrics = {assembly_key(a):assembly_metrics(a.placements,target) for a in answers}
    expected = {}
    remaining = set(metrics)
    layer = 1
    while remaining:
        front = {k for k in remaining if not any(dominates(metrics[j],metrics[k])
                                                 for j in remaining)}
        expected.update((k,layer) for k in front)
        remaining -= front
        layer += 1
    actual = pareto_assembly_ranks(reversed(answers),target)
    assert {k:r[0] for k,r in actual.items()} == expected
    # High retention/more cuts and low retention/fewer cuts are incomparable.
    a,b = metrics[assembly_key(answers[8])],metrics[assembly_key(answers[40])]
    assert not dominates(a,b) and not dominates(b,a)


def test_species_only_breaks_equal_objective_ties_not_pareto_layers():
    from rxn_core.retrosynthesis.beta import BetaRecommendation
    a,b = block('a',(0,1),refined=True),block('b',(2,3),refined=True)
    distinct = BetaRecommendation((a,b),())
    repeated = BetaRecommendation((a,replace(b,candidate=replace(b.candidate,source_id='a'))),())
    ranks = pareto_assembly_ranks((distinct,repeated),FakeBank().target)
    left,right = ranks[assembly_key(repeated)],ranks[assembly_key(distinct)]
    assert left[:3] == right[:3]
    assert left[0] == right[0] == 1
    assert left < right

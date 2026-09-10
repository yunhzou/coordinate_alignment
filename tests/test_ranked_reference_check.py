"""Saved ranking certificates implement the unchanged representative check."""
from types import SimpleNamespace as S

import pytest

from golden_evaluation import prepare, evaluate_planned, project
from publication_analysis import rank_archive
from ranked_reference_check import check_ranked_reference
from rxn_core import AAMProblem, AAMSearchConfig
from rxn_core.search_orientation import AAMSearchPlan


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('partial', [False, True])
@pytest.mark.parametrize('unbalanced', [False, True])
def test_saved_certificate_equals_original_representative_evaluation(reverse, partial, unbalanced):
    reaction = '[CH3:1][OH:2]' + ('.O' if unbalanced else '') + '>>[CH3:1][OH:2]'
    problem, features, reference = prepare(reaction)
    search_problem = AAMProblem(problem.product, problem.reactant) if reverse else problem
    plan = AAMSearchPlan(problem, search_problem, AAMSearchConfig(), reverse)
    mapping = dict(reference)
    if partial:
        mapping.pop(next(iter(mapping)))
    else:
        rh = [i for i,e in enumerate(problem.reactant.elements) if e=='H']
        ph = [i for i,e in enumerate(problem.product.elements) if e=='H']
        mapping.update(zip(rh,ph))
    reference = {r:p for r,p in mapping.items() if problem.reactant.elements[r]!='H'}
    # Different full-H witnesses share a heavy-atom class; selection must retain
    # the same score, coverage and exact relation (including unmatched atoms).
    alternatives = [mapping, dict(reversed(list(mapping.items())))]
    graph = S(terminals=(0,1), states=[S(mapping=tuple(plan.to_search_mapping(m).items()))
              for m in alternatives], capped=False, paths=lambda:iter(()))
    aam = S(graph=graph, problem=search_problem)
    classes = rank_archive(aam,plan)
    actual = check_ranked_reference(classes,plan,dict(features=features,mapping=list(reference.items())),
                                   dict(terminals=2,capped=False))
    expected = evaluate_planned(aam,plan,features,reference,symbolic=False)
    for key in ('reference_pairs','reference_annotation_complete','top1_correct',
                'representative_recovery','reference_recovery','top_terminal','witness_terminal',
                'candidate_terminals','unique_representative_chemistries','best_target_heavy_coverage',
                'best_target_all_atom_coverage','capped','top_events','input_orientation_witness'):
        assert actual[key] == expected[key], key
    if not partial:
        incomplete = dict(reference)
        incomplete.pop(next(iter(incomplete)))
        assert check_ranked_reference(classes,plan,dict(features=features,mapping=incomplete),
                                      dict(terminals=2,capped=False)) is None

from dataclasses import replace
import sys
from pathlib import Path

import numpy as np
import pytest

from rxn_core import AAMProblem,AAMSearchConfig,plan_aam_search,search_aam
from rxn_core.aam import cut_seed,checkpoint_manifest
from rxn_core.domain import MolecularEndpoint
from rxn_core.alignment.branch import _generate_seed_orders
from rxn_core.frag import build_graph


def endpoint(n):
    return MolecularEndpoint(('C',)*n,np.zeros((n,3)),np.zeros((n,n)))


def test_cut_streams_are_distinct_stable_and_order_independent():
    assert cut_seed(())==42
    assert cut_seed(((1,2),(3,4)))==cut_seed(((4,3),(2,1)))
    assert len({cut_seed(((i,i+1),)) for i in range(100)})==100
    matrix=np.ones((5,5))-np.eye(5)
    graph=build_graph(('C',)*5,matrix)
    first=_generate_seed_orders(graph,3,rng_seed=cut_seed(((0,1),)))
    assert first==_generate_seed_orders(graph,3,rng_seed=cut_seed(((1,0),)))
    assert first!=_generate_seed_orders(graph,3,rng_seed=cut_seed(((0,2),)))


def test_smaller_first_inverts_anchors_but_not_the_search_graph_contract():
    original=AAMProblem(endpoint(4),endpoint(2),'unequal')
    config=AAMSearchConfig(anchors=((3,1),))
    plan=plan_aam_search(original,config)
    assert plan.reversed and plan.problem.reactant is original.product
    assert plan.config.anchors==((1,3),) and config.anchors==((3,1),)
    assert plan.to_input_mapping({0:2,1:3})=={2:0,3:1}
    assert plan.to_search_mapping({2:0,3:1})=={0:2,1:3}
    assert not plan_aam_search(AAMProblem(endpoint(2),endpoint(2))).reversed
    assert not plan_aam_search(AAMProblem(endpoint(1),endpoint(2))).reversed


def test_old_checkpoint_seed_policy_is_not_reused(tmp_path):
    import json
    problem=AAMProblem(endpoint(1),endpoint(1))
    config=AAMSearchConfig(cut_floor=10)
    old=checkpoint_manifest(problem,config)
    old.pop('seed_policy');old['schema']='rxn_core.aam_checkpoints/v1'
    (tmp_path/'manifest.json').write_text(json.dumps(old))
    with pytest.raises(ValueError,match='configuration differs'):
        search_aam(problem,config,intermediate_dir=tmp_path,resume=True)


def test_reverse_evaluation_keeps_original_reference_completeness_and_ranking():
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
    from golden_evaluation import prepare,evaluate_planned,rank_key
    problem,features,reference=prepare('[CH3:1][OH:2].[Na+]>>[CH3:1][OH:2]')
    config=AAMSearchConfig(seed_count=3,branch_limit=100,cut_floor=10)
    plan=plan_aam_search(problem,config)
    assert plan.reversed
    result=search_aam(plan.problem,plan.config)
    score=evaluate_planned(result,plan,features,reference)
    assert score['reference_annotation_complete'] and score['top1_correct']
    assert score['reference_recovery']=='recovered'
    assert score['best_target_heavy_coverage']==score['best_target_all_atom_coverage']==1
    mapping=plan.to_input_mapping(result.graph.states[score['top_terminal']].mapping)
    assert score['top_events']==rank_key(mapping,problem)[1]


def test_partial_campaign_scores_saved_cuts_without_claiming_full_search(tmp_path):
    import json
    from types import SimpleNamespace
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
    from golden_policy_campaign import initialize,search,partial
    audit=tmp_path/'audit.jsonl'
    audit.write_text(json.dumps(dict(index=0,mapped_reaction='[CH3:1][OH:2].[Na+]>>[CH3:1][OH:2]'))+'\n')
    args=SimpleNamespace(run=tmp_path/'run',audit=audit,indices=None,cpu_budget=1,index=0)
    initialize(args)
    search(args)
    partial(args)
    score=json.loads((args.run/'0/evaluation.json').read_text())
    assert score['reference_recovery']=='recovered'
    assert score['search_incomplete'] and score['top1_correct'] is None

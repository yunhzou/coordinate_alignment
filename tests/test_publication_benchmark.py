import time
from dataclasses import replace

from golden_evaluation import prepare
from publication_timing import SearchProfiler
from publication_analysis import rank_archive, merge_classes, certificate_id, representative_metrics, union_outcome
from rxn_core import AAMSearchConfig, plan_aam_search, search_aam


def test_profiling_does_not_change_graph_or_parallel_seed_results(tmp_path, monkeypatch):
    from rxn_core import artifacts
    p,features,reference=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    c=AAMSearchConfig(seed_count=10,branch_limit=100)
    expected=search_aam(p,c,workers=1)
    original=artifacts.write_raw_cut
    def slower(*args,**kwargs):
        time.sleep(.01)
        return original(*args,**kwargs)
    monkeypatch.setattr(artifacts,'write_raw_cut',slower)
    with SearchProfiler(tmp_path/'events') as profiler:
        actual=search_aam(p,c,workers=2,intermediate_dir=tmp_path/'cuts',archive_format='checkpoint')
    assert actual.graph.to_record()==expected.graph.to_record()
    metrics=profiler.summary()
    assert metrics['phases']['persistence/workers']['summed_wall_seconds']>=.01
    assert 0<=metrics['compute_cpu_excluding_persistence_and_loading_seconds']<=metrics['total_cpu_including_io_seconds']
    assert artifacts.write_raw_cut is slower
    assert AAMSearchConfig().seed_count==3  # publication settings do not alter defaults


def test_rank_and_mode_merge_keep_reference_out_of_selection():
    p,features,reference=prepare('[CH3:1][OH:2].[Na+]>>[CH3:1][OH:2]')
    plan=plan_aam_search(p,AAMSearchConfig(seed_count=10,cut_floor=10))
    aam=search_aam(plan.problem,plan.config)
    ranked=rank_archive(aam,plan)
    assert ranked
    merged=merge_classes([ranked,ranked])
    assert len(merged)==len(ranked)
    assert len(merged[0]['origins'])==2
    expected=certificate_id(features,reference)
    metrics=representative_metrics(merged,expected)
    assert metrics['topk']['1']


def test_union_never_turns_unknown_into_a_negative():
    assert union_outcome(['not_recovered','unknown'])=='unknown'
    assert union_outcome(['recovered','unknown'])=='recovered'
    assert union_outcome(['not_recovered','not_recovered'])=='not_recovered'


def test_finished_incomplete_ranking_is_unknown_not_pending(tmp_path):
    import json
    from types import SimpleNamespace
    from golden_policy_campaign import save
    from golden_publication import report
    save(tmp_path/'manifest.json',dict(indices=[0,1]))
    (tmp_path/'results/0').mkdir(parents=True)
    save(tmp_path/'results/0/result.json',dict(index=0,modes={
        name:dict(reference_recovery='unknown',ranking_complete=False)
        for name in ('single','bidirectional')}))
    report(SimpleNamespace(run=tmp_path))
    result=json.loads((tmp_path/'summary.json').read_text())
    for mode in result['modes'].values():
        assert mode['top5_family_outcomes']==dict(unknown=1,pending=1)

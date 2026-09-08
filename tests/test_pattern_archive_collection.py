from golden_evaluation import prepare
from view_golden_mapping import molecules
from collect_golden_patterns import collect,first_paths
from rxn_core import search_aam,AAMSearchConfig
from rxn_core.search_orientation import plan_aam_search


def test_archive_completion_resume_and_provenance(tmp_path):
    reaction='CC>>CC'
    problem,_,_=prepare(reaction)
    plan=plan_aam_search(problem,AAMSearchConfig(seed_count=1,cut_floor=10))
    aam=search_aam(plan.problem,plan.config)
    archive=tmp_path/'archive';archive.write_text('immutable archive identity fixture')
    out=tmp_path/'checkpoint.json'
    args=(archive,aam,plan,molecules(reaction,problem),list(aam.graph.terminals),out)
    state,patterns,eq=collect(*args,seconds=5,per_path_seconds=1)
    assert state['complete'] and patterns
    assert all(p['origins'] for p in patterns)
    repeated,again,_=collect(*args,seconds=5,per_path_seconds=1)
    assert repeated['complete']
    assert {p['key'] for p in patterns}=={p['key'] for p in again}
    first=first_paths(aam.graph)
    for terminal in aam.graph.terminals:
        assert first(terminal).transitions==next(aam.graph.paths(terminal)).transitions

"""The public reuse backend preserves the entire conditional search graph."""
import numpy as np
import pytest

from rxn_core import AAMProblem, AAMSearchConfig, search_aam
from rxn_core.domain import MolecularEndpoint


def problem(n=6, m=6):
    def endpoint(count, shift):
        w=np.zeros((count,count))
        for i in range(count-1):w[i,i+1]=w[i+1,i]=1.
        if shift:w[0,1]=w[1,0]=.35
        return MolecularEndpoint(('C',)*(count-2)+('H','H'),np.zeros((count,3)),w)
    return AAMProblem(endpoint(n,False),endpoint(m,True))


@pytest.mark.parametrize('workers',[1,2])
@pytest.mark.parametrize('tol',[.5,1.])
def test_public_backends_have_identical_graphs_and_growth_counts(workers,tol):
    case=problem()
    config=AAMSearchConfig(seed_count=3,branch_limit=3,iso_tolerance=tol)
    reference=search_aam(case,config,workers=workers)
    reused=search_aam(case,config,workers=workers,execution='reused_native')
    assert reference.graph==reused.graph
    for name in ('max_live_branches','max_growth_candidates','subtree_branch_cap_count','cut_count'):
        assert getattr(reference.metrics,name)==getattr(reused.metrics,name)


@pytest.mark.parametrize('sizes',[(5,7),(7,5)])
def test_partial_composition_anchors_and_distinct_cut_floor(sizes):
    case=problem(*sizes)
    config=AAMSearchConfig(seed_count=2,branch_limit=100,anchors=((0,0),),cut_floor=.1,graph_floor=.5)
    assert search_aam(case,config).graph==search_aam(case,config,execution='reused_native').graph


@pytest.mark.parametrize('workers',[1,2])
def test_cross_backend_checkpoint_resume_requires_no_mapping_rerun(tmp_path,monkeypatch,workers):
    import rxn_core.aam as module
    case=problem();config=AAMSearchConfig(seed_count=2)
    original=search_aam(case,config,workers=workers,execution='reused_native',
        intermediate_dir=tmp_path,archive_format='checkpoint')
    def forbidden(*a,**kw):raise AssertionError('completed search reran')
    monkeypatch.setattr(module,'_search_cut',forbidden)
    resumed=search_aam(case,config,workers=workers,intermediate_dir=tmp_path,resume=True,
                       archive_format='checkpoint')
    assert original.graph==resumed.graph


def test_native_requirement_is_checked_before_starting_worker_pool(monkeypatch):
    monkeypatch.setenv('RXN_CORE_NATIVE','0')
    with pytest.raises(ValueError,match='requires the built native engine'):
        search_aam(problem(),workers=2,execution='reused_native')
    with pytest.raises(ValueError,match='unknown AAM execution'):
        search_aam(problem(),execution='unknown')

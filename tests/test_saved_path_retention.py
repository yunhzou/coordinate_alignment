from dataclasses import replace

from golden_evaluation import prepare
from rxn_core import AAMSearchConfig, search_aam
from verify_saved_path_retention import identical_path


def test_retention_requires_a_connected_identical_compressed_path():
    problem,_,_=prepare('[CH3:1][OH:2]>>[CH3:1][OH:2]')
    config=AAMSearchConfig(seed_count=3,cut_floor=10,iso_tolerance=1.)
    original=search_aam(problem,config,execution='reused_native')
    candidate=search_aam(problem,config,execution='shared_policies')
    for path in original.graph.paths():
        found=identical_path(path,candidate.graph)
        assert found is not None
        assert found.mapping==path.mapping
        assert found.fragments==path.fragments
    path=next(original.graph.paths())
    # Same endpoint states/coverage do not justify a certificate when the
    # correlated fragment constraints differ or the connecting edges are absent.
    changed=replace(candidate.graph,transitions=tuple(replace(edge,
        match={**edge.match,'deferred_edges':((0,1),)}) if edge.match else edge
        for edge in candidate.graph.transitions))
    assert identical_path(path,changed) is None
    assert identical_path(path,replace(candidate.graph,transitions=())) is None

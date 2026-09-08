import json
import pytest

from rxn_core import aam
from rxn_core.search_graph import AAMSearchGraph, SearchContext, SearchGraphBuilder
from rxn_core.artifacts import raw_cut_paths,read_raw_cut,write_raw_cut


def test_cut_worker_persists_graph_and_returns_only_reference(tmp_path, monkeypatch):
    graph = SearchGraphBuilder(SearchContext((), (), ())).finish()
    monkeypatch.setattr(aam, '_search_cut', lambda cut: (graph, {
        'search_seconds': 1., 'max_live_branches': 2, 'max_growth_candidates': 3}))
    path = tmp_path/'cut_00007.json'
    index, payload, counts = aam._search_cut_task((7, (), str(path)))
    assert index == 7 and payload == str(path)
    assert AAMSearchGraph.from_record(json.loads(path.read_text())) == graph
    assert counts['max_live_branches'] == 2
    assert counts['checkpoint_seconds'] >= 0
    assert not path.with_suffix('.json.tmp').exists()
    _,owned,_=aam._search_cut_task((7,(),str(path)),in_process=True)
    assert owned is graph  # serial calls do not serialize then reload the graph


def test_in_memory_worker_and_online_profile_reducer(monkeypatch):
    graph = SearchGraphBuilder(SearchContext((), (), ())).finish()
    monkeypatch.setattr(aam, '_search_cut', lambda cut: (graph, {'search_seconds': 0.}))
    assert aam._search_cut_task((4, (), None)) == (4, graph, {
        'search_seconds': 0., 'checkpoint_seconds': 0.})
    reducer = aam._GrowthCounts()
    rows = [{}, {'max_cands_before': 3}, {'result': 'subtree_branch_cap'}, {'max_cands_before': 2}]
    for row in rows: reducer.append(row)
    assert reducer.maximum == max(row.get('max_cands_before', 0) for row in rows)


def test_raw_cut_formats_and_order_are_explicit(tmp_path):
    graph = SearchGraphBuilder(SearchContext((), (), ())).finish()
    binary=tmp_path/'cut_00002.raw.pkl.gz';text=tmp_path/'cut_00001.json'
    write_raw_cut(graph,binary);write_raw_cut(graph,text)
    assert raw_cut_paths(tmp_path)==[text,binary]
    assert read_raw_cut(binary)==read_raw_cut(text)==graph
    with pytest.raises(ValueError,match='Unsupported raw cut format'):
        write_raw_cut(graph,tmp_path/'unknown.bin')
    write_raw_cut(graph,tmp_path/'cut_00002.json')
    with pytest.raises(ValueError,match='Duplicate raw cut'):
        raw_cut_paths(tmp_path)


def test_binary_worker_sends_only_saved_reference(tmp_path,monkeypatch):
    graph = SearchGraphBuilder(SearchContext((), (), ())).finish()
    monkeypatch.setattr(aam,'_search_cut',lambda cut:(graph,{'search_seconds':0.}))
    path=tmp_path/'cut_00000.raw.pkl.gz'
    _,payload,_=aam._search_cut_task((0,(),str(path)))
    assert payload==str(path) and read_raw_cut(path)==graph

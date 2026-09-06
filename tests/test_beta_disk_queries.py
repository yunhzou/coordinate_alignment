"""The distributed disk index preserves the beta's proposal order and evidence."""
import gzip
import importlib.util
import json
from pathlib import Path
import pickle
import sys

from rxn_core.retrosynthesis.beta import FragmentQueryBank, proposal_rank, recommend_big_blocks
from rxn_core.smiles import smiles_to_weighted_graph


def module(name):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).parents[1]/'bench'/f'{name}.py')
    value=importlib.util.module_from_spec(spec)
    sys.modules[name]=value
    spec.loader.exec_module(value)
    return value


def test_disk_query_matches_in_memory_and_returns_a_full_cover(tmp_path,monkeypatch):
    index=module('index_beta_query')
    runner=module('run_beta_distributed')
    query=tmp_path/'query_full'
    for part in ('sources','indexes','parts'):
        (query/part).mkdir(parents=True)
    catalog=tmp_path/'bank.csv.gz'
    with gzip.open(catalog,'wt') as f:
        f.write('Bank ID,SMILES\nmethane,C\nwater,O\n')
    manifest=dict(target_smiles='CO',catalog=str(catalog),id_column='Bank ID',
                  region=list(range(6)),shards=1,bank_rows=2,
                  config={'iso_tolerance':1.0,'branch_limit':100})
    (query/'manifest.json').write_text(json.dumps(manifest))
    graph=lambda s:smiles_to_weighted_graph(s,expand_hydrogens=True)
    memory=FragmentQueryBank([('methane',graph('C')),('water',graph('O'))],graph('CO'))
    blocks=memory.query(range(6))
    for row,source in enumerate(('methane','water')):
        with gzip.open(query/'sources'/f'{row}.blocks.pkl.gz','wb') as f:
            pickle.dump(tuple(b for b in blocks if b.source_id==source),f)
    (query/'parts/0.json').write_text(json.dumps(dict(rows=2,blocks=len(blocks),capped=0)))
    monkeypatch.setattr(sys,'argv',['index','--query',str(query),'--shard','0'])
    index.main()
    disk=runner.DistributedBank(tmp_path,'')
    actual=list(disk.ordered_query(range(6),()))
    expected=list(memory.ordered_query(range(6),()))
    assert {b.key for b in actual}=={b.key for b in expected}
    assert [proposal_rank((b,),6) for b in actual]==[proposal_rank((b,),6) for b in expected]
    # Exercise the lazy alternative frontier while reusing local gap queries;
    # Slurm submission is outside this unit test.
    original=disk.ordered_query
    disk.ordered_query=lambda region,placements: original(region,placements) if len(region)==6 else memory.ordered_query(region,placements)
    result=recommend_big_blocks(disk, recommendations=1)
    assert result.recommendations[0].uncovered_target_atoms==()
    assert all(p.refined for p in result.recommendations[0].placements)
    cached=runner.DistributedBank(tmp_path,'')
    def forbidden(*args,**kwargs):
        raise AssertionError('saved selected-source AAM must not be repeated')
    monkeypatch.setattr(FragmentQueryBank,'detect_selected',forbidden)
    assert cached.detect_selected('methane').candidates

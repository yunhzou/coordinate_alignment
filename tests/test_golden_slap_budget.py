import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bench'))
from golden_slap_budget import prepare
from golden_competitors import signatures


def test_order_and_direction_plan_is_reproducible_and_preserves_endpoints(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    source=tmp_path/'data/aam_benchmarks/golden_original_20260906'
    source.mkdir(parents=True)
    reaction='C[C@H](O)CC.N>>CC(=O)CC.N'
    (source/'audit.jsonl').write_text(json.dumps(dict(index=0,input_reaction=reaction))+'\n')
    for name in ('a','b'):
        prepare(argparse.Namespace(run=tmp_path/name,orders=3,seed=42,symmetry='heavy'))
    assert (tmp_path/'a/inputs.jsonl').read_bytes()==(tmp_path/'b/inputs.jsonl').read_bytes()
    inputs=[json.loads(x) for x in (tmp_path/'a/inputs.jsonl').read_text().splitlines()]
    plans=json.loads((tmp_path/'a/plans.json').read_text())
    assert len(inputs)==6
    for row,plan in zip(inputs,plans):
        text=row['input_reaction']
        if plan['reverse']:text='>>'.join(text.split('>>')[::-1])
        assert signatures(text)[0]==signatures(reaction)[0]
    assert [x['reverse'] for x in plans]==[False,True]*3
    manifest=json.loads((tmp_path/'a/manifest.json').read_text())
    assert manifest['modes']==['slap_binary','slap_weighted']
    assert manifest['output_limit'] is None

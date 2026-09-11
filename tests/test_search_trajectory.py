"""Scientific replay integrity and a manifest pipeline independent of PR7."""
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from rxn_core import AAMProblem,AAMSearchConfig,MolecularEndpoint,search_aam
from rxn_core.artifacts import write_aam_checkpoint, read_aam_checkpoint, aam_json
from rxn_core.search_trajectory import build_trajectory,growth,extension
from rxn_core.viewers import growth_trace_html


def toy_archive(folder, suffix='', translation=0):
    elements=['C','O','O','H'];coords=np.array([[0.,0.,0.],[1.2,0,0],[-1,0,0],[-1,1,0]])
    w=np.array([[0,1.8,.8,0],[1.8,0,0,0],[.8,0,0,.9],[0,0,.9,0.]])
    order=[2,0,3,1]
    r=MolecularEndpoint(elements=elements,coordinates=coords,wbo=w)
    p=MolecularEndpoint(elements=[elements[i] for i in order],coordinates=coords[order]+translation,wbo=w[np.ix_(order,order)])
    result=search_aam(AAMProblem(r,p,name='four_atom_example'),
        AAMSearchConfig(iso_tolerance=.1,graph_floor=.4,branch_limit=32,seed_count=1),workers=1)
    archive=folder/f'toy{suffix}.pkl.gz';write_aam_checkpoint(result,archive)
    return archive


def test_replay_uses_archive_config_and_preserves_inputs(tmp_path):
    archive=toy_archive(tmp_path);before=archive.read_bytes()
    original=growth._extend_sym_cands,extension._dedupe_children
    doc=build_trajectory([dict(archive=str(archive),context=0)])
    assert archive.read_bytes()==before
    assert (growth._extend_sym_cands,extension._dedupe_children)==original
    assert doc['runs'][0]['config']['iso_tolerance']==.1
    assert doc['runs'][0]['config']['graph_floor']==.4
    for path in doc['runs'][0]['paths']:
        assert path['complete'] and len(path['mapping'])==4
        assert path['events']==[]
        assert any(f['kind']=='commit' for f in path['frames'])
        assert all(check['archived_fragment_matches']=='exact' for check in doc['checks'])
    html=growth_trace_html(doc)
    assert 'PR7' not in html and 'data-viewer-layout="growth_trace"' in html


def test_manifest_cli_and_render_only(tmp_path,monkeypatch):
    archive=toy_archive(tmp_path)
    interchange=tmp_path/'toy.json';interchange.write_text(aam_json(read_aam_checkpoint(archive)))
    manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps(dict(selections=[dict(archive=interchange.name,context=0)])))
    root=Path(__file__).resolve().parents[1];output=tmp_path/'output'
    subprocess.run([sys.executable,str(root/'tools/build_search_trajectory.py'),str(manifest),str(output)],
                   check=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},capture_output=True)
    doc=json.loads((output/'trace.json').read_text())
    def fail(*a,**k):raise AssertionError('Rendering must not replay or search')
    monkeypatch.setattr(growth,'grow_island',fail)
    assert growth_trace_html(doc)==(output/'algorithm_trajectory.html').read_text()
    assert json.loads((output/'validation.json').read_text())['status']=='passed'


def test_mismatched_endpoints_are_rejected(tmp_path):
    a,b=toy_archive(tmp_path),toy_archive(tmp_path,'_other',translation=2)
    with pytest.raises(ValueError,match='identical endpoint'):
        build_trajectory([dict(archive=str(a)),dict(archive=str(b))])

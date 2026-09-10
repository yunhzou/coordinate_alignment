"""Check manuscript evidence, mapping invariants, and rendered artifacts."""
import hashlib
import json
import re
from pathlib import Path
import numpy as np
import imageio.v2 as imageio
import pymupdf
from PIL import Image

MAN=Path(__file__).resolve().parents[1]
def read(name):return json.loads((MAN/'evidence'/name).read_text())
sources=read('sources.json')
for name,record in sources.items():
    assert hashlib.sha256((MAN/'evidence'/name).read_bytes()).hexdigest()==record['sha256'],name

seed=read('seed_comparison.json')
assert len(seed['common_case_indices'])==1807
for key,n in [('seeds1',1833),('seeds2',1834),('seeds3',1836),('seeds10',1836)]:
    m=seed['methods'][key]
    assert m['golden_cases']==1851 and m['golden_outcomes']['recovered']==n
    assert sum(m['golden_outcomes'].values())==1851
assert read('slap_sweep.json')['sweep_union_recovered']==1795
a=read('adaptive_no_sweep.json')['totals']['golden_bidirectional']['adaptive']
assert a['reference_recovery']=={'recovered':1723,'unknown':95,'not_recovered':33}
assert a['searches_complete']==a['expected_searches']==3702

traces=read('animation_data.json');checks=0
for trace in traces:
    r,p=(trace['input'][k] for k in ['reactant','product'])
    for frame in trace['frames']:
        for raw in [frame['mapping'],*frame['samples']]:
            m={int(a):int(b) for a,b in raw.items()}
            assert len(m)==len(set(m.values()))
            assert set(frame['active'])<=m.keys()
            assert all(r['elements'][a]==p['elements'][b] for a,b in m.items())
            checks+=1
        assert len(frame['samples'])<=5
        if frame['event']=='terminal':assert len(frame['mapping'])==len(r['elements'])
    assert not any(s['reason']=='branch_cap' for s in trace['graph']['stops'])
assert max(f['candidates'] for f in traces[0]['frames'])==42
assert sum(f['event']=='terminal' for f in traces[0]['frames'])==2
assert sum(f['event']=='consumed' for f in traces[1]['frames'])==8

figures=['fig1_algorithm','fig2_golden','fig3_holdout','fig4_event_windows','fig5_growth']
for name in figures:
    for suffix in ['pdf','png','svg']:assert (MAN/'figs'/f'{name}.{suffix}').stat().st_size>1000
    assert len(pymupdf.open(MAN/'figs'/f'{name}.pdf'))==1

movies=[]
for m in json.loads((MAN/'animations/movies.json').read_text()):
    reader=imageio.get_reader(MAN/'animations'/m['file'])
    meta=reader.get_meta_data();count=reader.count_frames()
    assert meta['size']==(1280,720) and meta['fps']==10
    assert abs(count/meta['fps']-m['duration_seconds'])<.11
    for index in [0,count//2,count-1]:
        frame=reader.get_data(index);assert np.std(frame)>10
    reader.close()
    gif=Image.open(MAN/'animations'/m['file'].replace('.mp4','.gif'))
    assert gif.n_frames==m['events']
    movies.append(dict(file=m['file'],frames=count,duration_seconds=m['duration_seconds']))

doc=pymupdf.open(MAN/'manuscript.pdf');text='\n'.join(p.get_text() for p in doc)
assert len(doc)>=9
assert '??' not in text
for phrase in ['99.03','99.19','96.97','93.08','1,851','1,723','Supporting Information']:
    assert phrase in text,phrase
log=(MAN/'build/preprint.log').read_text()
assert not re.search(r'Overfull|undefined|Missing character|^!',log,re.M)
out=MAN/'build';out.mkdir(exist_ok=True)
(out/'manuscript-text.txt').write_text(text)
result=dict(status='passed',pdf_pages=len(doc),evidence_hashes=len(sources),
            mapping_witnesses_checked=checks,figures=len(figures),movies=movies,
            note='Validates artifact consistency and rendering, not scientific completeness or chemical accuracy.')
(out/'artifact-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

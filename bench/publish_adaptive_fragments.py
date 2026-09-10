"""Collect the bounded experiment's saved evidence; never runs AAM."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from adaptive_fragment_pilot import save


RUNS = {
    'initial': 'adaptive_fragment_20260910_lRLOXa',
    'deeper': 'adaptive_fragment_deeper_20260910_mzYH3n',
    'feedback': 'adaptive_fragment_feedback_20260910_I8Kjaa',
    'fair': 'adaptive_fragment_fair_20260910_fKvIsb',
}


def collect(args):
    root=Path(__file__).resolve().parents[1]
    reference=json.loads((root/'reports/dependency_fragment_repair_20260909/mechanism_check_114.json').read_text())
    records=[];provenance={}
    for name,directory in RUNS.items():
        run=args.bank/directory
        implementation=hashlib.sha256()
        sources=list((run/'engine/src').rglob('*.py'))+list((run/'engine/native/src').rglob('*.h'))+list((run/'engine/native/src').rglob('*.cpp'))
        sources.append(run/'engine/bench/adaptive_fragment_pilot.py')
        for source in sorted(sources):
            implementation.update(str(source.relative_to(run/'engine')).encode()+b'\0')
            implementation.update(hashlib.sha256(source.read_bytes()).digest())
        provenance[name]=dict(manifest=json.loads((run/'manifest.json').read_text()),
            implementation_sha256=implementation.hexdigest(),source_files=len(sources))
        for path in sorted((run/'results').glob('*/summary.json')):
            folder=path.parent;summary=json.loads(path.read_text())
            record=dict(run=name,slot=int(folder.name),path=str(folder),**summary,
                validation=json.loads((folder/'validation.json').read_text()))
            diagnostic=folder/'family_diagnostic.json'
            if diagnostic.exists():
                d=json.loads(diagnostic.read_text())
                record['family_diagnostic']=dict(best=d['best'],paths=len(d['units']),seconds=d['seconds'],scope=d['scope'])
            if summary['index']==114:
                raw=json.loads((run/'inputs/114.json').read_text())
                r,p=[np.asarray(raw[n]['wbo']) for n in ('reactant','product')]
                hits={str(row['id']):None for row in reference['references']}
                label=summary['rows'][-1]['label']
                witnesses=json.loads((folder/f'{label}_witnesses.json').read_text())
                for ordinal,vector in enumerate(witnesses['mappings']):
                    aligned=p[np.ix_(vector,vector)];rb,pb=r>.2,aligned>.2
                    masks=(rb & ~pb,~rb & pb,rb & pb & (np.abs(r-aligned)>.5))
                    signature=[np.argwhere(np.triu(mask,1)).tolist() for mask in masks]
                    for ref in reference['references']:
                        if signature==ref['signature'] and hits[str(ref['id'])] is None:
                            hits[str(ref['id'])]=ordinal
                record['historical_114_patterns']=dict(representative_ordinals=hits,
                    scope='Exact source-index event signatures in saved representatives, not chemical ground truth or exhaustive family exclusion.')
            records.append(record)
    final=args.bank/RUNS['fair']
    tests=ET.parse(final/'tests_final.xml').getroot()
    save(args.output,dict(runs={name:str(args.bank/directory) for name,directory in RUNS.items()},
        provenance=provenance,records=records,tests=[dict(s.attrib) for s in tests.iter('testsuite')],
        default_regression=json.loads((final/'default_regression/comparison.json').read_text()),
        saved_full_sweep_representative_best={'25':5,'76':3,'77':8,'114':4},
        scope='Small fixed-direction pilot. Search times exclude persistence and diagnostic scoring. No proof of optimality, exhaustive coverage or production speedup.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bank',type=Path,default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks'))
    p.add_argument('--output',type=Path,required=True)
    collect(p.parse_args())

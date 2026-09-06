"""Show a failed blind run honestly, optionally beside independently saved validation."""
import argparse
import gzip
import json
from pathlib import Path
import pickle
import subprocess
import sys

from rxn_core.retrosynthesis.beta import BetaResult, BetaRecommendation


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',type=Path,required=True)
    p.add_argument('--reason',required=True)
    p.add_argument('--validation-stem')
    args=p.parse_args()
    case=json.loads(args.case.read_text());run=Path(case['run'])
    manifest=json.loads((run/'query_full/manifest.json').read_text())
    parts=[json.loads(path.read_text()) for path in (run/'query_full/parts').glob('*.json')]
    events=[]
    for path in (run/'refinements').glob('*.pkl.gz'):
        with gzip.open(path,'rb') as stream:
            event,_=pickle.load(stream);events.append(event)
    result=BetaResult((),BetaRecommendation((),tuple(manifest['region'])),
        sum(p['capped'] for p in parts)+sum(bool(e['capped']) for e in events),0,tuple(events))
    with gzip.open(run/'failed_blind.pkl.gz','wb') as stream:
        pickle.dump(result,stream,protocol=pickle.HIGHEST_PROTOCOL)
    (run/'failed_blind.json').write_text(json.dumps(dict(recommendations=[],failure=args.reason))+'\n')
    summary=args.case.with_name('scan_summary.json')
    summary.write_text(json.dumps(dict(scan_counts=dict(matched_precursors=sum(p['matched_sources'] for p in parts))))+'\n')
    command=[sys.executable,'bench/view_beta_result.py','--run',str(run),
        '--result-stem','failed_blind','--case',str(args.case),'--scan-summary',str(summary),
        '--html-output',str(args.case.with_name('viewer.html')),
        '--title',case['id']+' · BLIND RUN FAILED']
    if args.validation_stem:
        command.extend(['--validation-stem',args.validation_stem])
    subprocess.run(command,check=True,timeout=580)


if __name__=='__main__':
    main()

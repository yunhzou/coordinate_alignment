"""Seed and branch-budget ablations on fixed-direction Golden misses."""
import argparse
from collections import Counter
from dataclasses import asdict,replace
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from golden_policy_campaign import load_case,save,guarded
from golden_domain_rescore import score as verify_saved
from rxn_core import AAMSearchConfig,search_aam


def changed_config(config,seeds,cap=None):
    changes=dict(seed_count=seeds)
    if cap is not None:changes['branch_limit']=cap
    return replace(config,**changes)


def init(args):
    args.run.mkdir(parents=True,exist_ok=False)
    baseline=json.loads(args.summary.read_text())
    original=json.loads((args.source/'manifest.json').read_text())
    config=changed_config(AAMSearchConfig(**original['config']),args.seeds,args.cap)
    indices=baseline['not_recovered']
    if args.selection:
        selection=json.loads(args.selection.read_text())
        indices=sorted(selection['not_recovered']+selection['unknown'])
        baseline['all_records']['recovered']=selection['union_recovered']
    assert len(indices)==len(set(indices))
    for index in indices:
        out=args.run/str(index);out.mkdir()
        for name in ('input.json','reference.json','orientation.json'):
            shutil.copy2(args.source/str(index)/name,out/name)
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    save(args.run/'manifest.json',dict(source=str(args.run.resolve()),baseline_source=str(args.source.resolve()),
        indices=indices,config=asdict(config),baseline_config=original['config'],workers=args.workers,
        search_watchdog=300,score_watchdog=240,baseline_summary=baseline,
        scope='Same default orientation and full single-edge sweep; only seed count and optional branch cap change. Selected prior misses: diagnostic, not a full-dataset new-policy benchmark.',
        selection_source=str(args.selection) if args.selection else None,
        source_hashes={str(p.relative_to(args.run/'engine')):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (args.run/'engine').rglob('*.py')}))
    for index in indices:
        _,plan=load_case(args.run,index)
        old=json.loads((args.source/str(index)/'orientation.json').read_text())
        assert plan.direction==old['direction']
        save(args.run/str(index)/'orientation.json',dict(old,config=asdict(config)))
    print(json.dumps(dict(cases=len(indices),seeds=config.seed_count,workers=args.workers)))


def search(args):
    m=json.loads((args.run/'manifest.json').read_text());index=m['indices'][args.slot]
    out,plan=load_case(args.run,index);start=time.perf_counter()
    result=search_aam(plan.problem,plan.config,workers=m['workers'],intermediate_dir=out/'cuts',archive_format='checkpoint')
    save(out/'search.json',dict(seconds=time.perf_counter()-start,metrics=asdict(result.metrics),direction=plan.direction))


def worker(args):
    m=json.loads((args.run/'manifest.json').read_text());out=args.run/str(m['indices'][args.slot])
    start=time.time()
    for phase,limit in [('search',m['search_watchdog']),('score',m['score_watchdog'])]:
        save(out/'status.json',dict(stage=phase,started=start,updated=time.time(),job=os.environ.get('SLURM_JOB_ID')))
        code=guarded([sys.executable,__file__,phase,'--run',str(args.run),'--slot',str(args.slot)],out/f'{phase}.log',limit)
        save(out/f'{phase}_status.json',dict(exit_code=code,finished=time.time()))
    if not (out/'evaluation.json').exists():save(out/'evaluation.json',dict(reference_recovery='unknown',reason='verification interrupted',search_incomplete=not (out/'cuts/aam.pkl.gz').exists()))
    save(out/'status.json',dict(stage='complete',started=start,finished=time.time(),job=os.environ.get('SLURM_JOB_ID')))


def report(args):
    m=json.loads((args.run/'manifest.json').read_text());rows=[];witnesses=[]
    for index in m['indices']:
        out=args.run/str(index);path=out/'evaluation.json'
        e=json.loads(path.read_text()) if path.exists() else dict(reference_recovery='pending')
        search_path=out/'search.json';s=json.loads(search_path.read_text()) if search_path.exists() else {}
        rows.append(dict(index=index,outcome=e['reference_recovery'],search_incomplete=e.get('search_incomplete'),
            top1=e.get('top1_correct'),search_seconds=s.get('seconds'),capped=s.get('metrics',{}).get('subtree_branch_cap_count'),
            evaluation_seconds=e.get('total_seconds'),reason=e.get('reason'),
            phase_status={phase:json.loads((out/f'{phase}_status.json').read_text())
                for phase in ('search','score') if (out/f'{phase}_status.json').exists()},
            evaluation=str(path)))
        if e['reference_recovery']=='recovered':witnesses.append(dict(index=index,**e))
    recovered=[r['index'] for r in rows if r['outcome']=='recovered']
    summary=dict(tested=len(rows),seeds=m['config']['seed_count'],cap=m['config']['branch_limit'],outcomes=dict(Counter(r['outcome'] for r in rows)),
        recovered=recovered,not_recovered=[r['index'] for r in rows if r['outcome']=='not_recovered'],
        unknown=[r['index'] for r in rows if r['outcome']=='unknown'],
        baseline_recovered=m['baseline_summary']['all_records']['recovered'],
        union_recovered=m['baseline_summary']['all_records']['recovered']+len(recovered),
        union_percent=100*(m['baseline_summary']['all_records']['recovered']+len(recovered))/m['baseline_summary']['all_records']['total'],
        note='Diagnostic union with previous saved results; not a full-dataset new-policy benchmark. Original unresolved case 1793 was not part of the requested 32 misses.')
    args.output.mkdir(parents=True,exist_ok=True);save(args.output/'summary.json',summary);save(args.output/'cases.json',rows)
    save(args.output/'witnesses.json',witnesses);save(args.output/'manifest.json',m)
    if args.jobs:(args.output/'slurm_accounting.psv').write_text(subprocess.check_output(['sacct','-j',args.jobs,'-P','--format=JobID,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,Start,End'],text=True))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['init','worker','search','score','report'])
    p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path);p.add_argument('--summary',type=Path)
    p.add_argument('--seeds',type=int,default=10);p.add_argument('--workers',type=int,default=4);p.add_argument('--slot',type=int)
    p.add_argument('--cap',type=int);p.add_argument('--selection',type=Path)
    p.add_argument('--output',type=Path);p.add_argument('--jobs');a=p.parse_args()
    dict(init=init,worker=worker,search=search,score=verify_saved,report=report)[a.mode](a)

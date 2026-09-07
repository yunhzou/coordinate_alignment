"""Reference-blind class ranking and bounded top-five family verification."""
import argparse
from collections import Counter
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pynauty

from golden_evaluation import colored_graph,project,rank_key,evaluate_planned
from golden_policy_campaign import load_case,save,guarded
from rxn_core.artifacts import read_aam_checkpoint


def init(args):
    args.run.mkdir(parents=True,exist_ok=False)
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    original=json.loads((args.source/'manifest.json').read_text())
    save(args.run/'manifest.json',dict(source=str(args.source.resolve()),
        indices=[r['index'] for r in original['records']],watchdog_seconds=300,
        semantics='Rank chemical classes by best saved explicit-H representative, no reference guidance. Exact top-five family recovery is a separate metric.',
        source_hashes={str(p.relative_to(args.run/'engine')):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (args.run/'engine').rglob('*.py')}))


def summarize_classes(classes):
    """Event windows retain the best heavy/all-atom coverage tier."""
    best=classes[0]['score'][:3] if classes else None
    reference=next((c for c in classes if c['reference_equivalent']),None)
    windows={}
    for tolerance in (0,1,2,3,5):
        allowed=[c for c in classes if c['score'][:2]==best[:2]
                 and c['score'][2]<=best[2]+tolerance] if best else []
        windows[str(tolerance)]=dict(count=len(allowed),recovered=any(c['reference_equivalent'] for c in allowed))
    boundary=classes[min(5,len(classes))-1]['score'][:3] if classes else None
    tied=[c for c in classes if c['score'][:3]<=boundary] if boundary else []
    return dict(reference_class_rank=None if reference is None else reference['rank'],
        reference_event_gap=None if reference is None or reference['score'][:2]!=best[:2] else reference['score'][2]-best[2],
        top5_tie_inclusive=dict(count=len(tied),recovered=any(c['reference_equivalent'] for c in tied)),
        event_windows=windows)


def score(args):
    manifest=json.loads((args.run/'manifest.json').read_text());index=manifest['indices'][args.slot]
    directory,plan=load_case(Path(manifest['source']),index);out=args.run/str(index);out.mkdir(exist_ok=True)
    archive=directory/'cuts/aam.pkl.gz';start=time.perf_counter()
    if not archive.exists():
        save(out/'result.json',dict(index=index,status='incomplete_search',top5_family='unknown'))
        return
    aam=read_aam_checkpoint(archive);ref=json.loads((directory/'reference.json').read_text())
    features=ref['features'];expected=pynauty.certificate(colored_graph(features,project(ref['mapping'],features)))
    cache={};groups={}
    for terminal in aam.graph.terminals:
        mapping=plan.to_input_mapping(aam.graph.states[terminal].mapping)
        heavy=tuple(sorted(project(mapping,features).items()))
        if heavy not in cache:cache[heavy]=pynauty.certificate(colored_graph(features,dict(heavy)))
        certificate=cache[heavy];key,events=rank_key(mapping,plan.input_problem)
        item=groups.setdefault(certificate,dict(terminals=[],best=None))
        item['terminals'].append(terminal)
        if item['best'] is None or key<item['best'][0]:item['best']=(key,terminal,events)
    ordered=sorted(groups.items(),key=lambda pair:pair[1]['best'][0]);classes=[]
    for rank,(certificate,item) in enumerate(ordered,1):
        key,terminal,events=item['best']
        classes.append(dict(rank=rank,reference_equivalent=certificate==expected,score=list(key[:3]),
            terminal=terminal,saved_terminals=len(item['terminals']),events=events))
    save(out/'classes.json',classes)
    result=dict(index=index,status='complete',class_count=len(classes),**summarize_classes(classes))
    # Save cheap metrics before the more expensive compressed-family query.
    result.update(top5_family='unknown',ranking_seconds=time.perf_counter()-start)
    save(out/'result.json',result)
    if any(c['reference_equivalent'] for c in classes[:5]):
        result['top5_family']='recovered'
    else:
        eligible={t for _,item in ordered[:5] for t in item['terminals']}
        view=replace(aam,graph=replace(aam.graph,stops=tuple(s for s in aam.graph.stops
            if s.reason not in {'objective_met','stalled'} or s.state in eligible)))
        evaluation=evaluate_planned(view,plan,features,ref['mapping'],seconds=120,query_timeout_ms=5000)
        save(out/'family_verification.json',evaluation)
        result['top5_family']=evaluation['reference_recovery']
    result['total_seconds']=time.perf_counter()-start;save(out/'result.json',result)


def worker(args):
    m=json.loads((args.run/'manifest.json').read_text());index=m['indices'][args.slot]
    out=args.run/str(index);out.mkdir(exist_ok=True)
    save(out/'status.json',dict(stage='running',job=os.environ.get('SLURM_JOB_ID'),started=time.time()))
    code=guarded([sys.executable,__file__,'score','--run',str(args.run),'--slot',str(args.slot)],out/'score.log',300)
    if not (out/'result.json').exists():save(out/'result.json',dict(index=index,status='worker_failure',top5_family='unknown'))
    save(out/'status.json',dict(stage='complete',exit_code=code,finished=time.time()))


def report(args):
    manifest=json.loads((args.run/'manifest.json').read_text());rows=[]
    for index in manifest['indices']:
        path=args.run/str(index)/'result.json'
        rows.append(json.loads(path.read_text()) if path.exists() else dict(index=index,status='pending',top5_family='unknown'))
    total=len(rows);completed=[r for r in rows if r['status']=='complete']
    def metric(count):return dict(recovered=count,total=total,percent=100*count/total)
    summary=dict(records=total,statuses=dict(Counter(r['status'] for r in rows)),
        representative_topk={str(k):metric(sum(r.get('reference_class_rank') is not None and r['reference_class_rank']<=k for r in rows)) for k in (1,3,5,10)},
        top5_family=dict(metric(sum(r['top5_family']=='recovered' for r in rows)),outcomes=dict(Counter(r['top5_family'] for r in rows))),
        top5_tie_inclusive=dict(metric(sum(r['top5_tie_inclusive']['recovered'] for r in completed)),
            mean_candidates=sum(r['top5_tie_inclusive']['count'] for r in completed)/max(1,len(completed))),
        event_windows={str(t):dict(metric(sum(r['event_windows'][str(t)]['recovered'] for r in completed)),
            mean_candidates=sum(r['event_windows'][str(t)]['count'] for r in completed)/max(1,len(completed))) for t in (0,1,2,3,5)},
        note='All 1851 remain in denominators. Incomplete search rankings unknown. Event-window/tie metrics use displayed representatives; top5-family includes full compressed alternatives behind the selected classes, not a globally optimized family score.')
    args.output.mkdir(parents=True,exist_ok=True);save(args.output/'summary.json',summary)
    save(args.output/'cases.json',rows);save(args.output/'manifest.json',manifest)
    if args.jobs:
        (args.output/'slurm_accounting.psv').write_text(subprocess.check_output(['sacct','-j',args.jobs,'-P',
            '--format=JobID,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,Start,End'],text=True))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['init','score','worker','report'])
    p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path)
    p.add_argument('--slot',type=int);p.add_argument('--output',type=Path);p.add_argument('--jobs')
    a=p.parse_args();dict(init=init,score=score,worker=worker,report=report)[a.mode](a)

"""Rescore saved Golden archives with full explicit domain/group verification."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from golden_policy_campaign import load_case,save,guarded
from golden_evaluation import evaluate_planned
from rxn_core.artifacts import read_aam_checkpoint,raw_cut_paths,read_raw_cut


def init(args):
    args.run.mkdir(parents=True,exist_ok=False)
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    old=json.loads((args.source/'progress.json').read_text())
    jobs=[]
    for p in sorted(args.source.glob('*/evaluation.json'),key=lambda p:int(p.parent.name)):
        e=json.loads(p.read_text())
        selected = (not e.get('reference_annotation_complete') if args.one_sided
                    else e.get('reference_annotation_complete') and e['reference_recovery']!='recovered')
        if selected:jobs.append(int(p.parent.name))
    save(args.run/'manifest.json',dict(source=str(args.source.resolve()),indices=jobs,
        cached_complete_reference_positives=old['complete_reference_recovered'],
        complete_reference_total=old['complete_reference_total'],watchdog_seconds=300,
        interpretation=('Score every one-sided reference with exact unmatched-atom semantics' if args.one_sided
                        else 'Reuse independently certified positives; rescore all remaining complete-reference records, no new AAM searches'),
        source_hashes={str(p.relative_to(args.run/'engine')):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (args.run/'engine').rglob('*.py')}))
    print(len(jobs))


def score(args):
    m=json.loads((args.run/'manifest.json').read_text());index=m['indices'][args.slot]
    source,plan=load_case(Path(m['source']),index);out=args.run/str(index);out.mkdir(exist_ok=True)
    ref=json.loads((source/'reference.json').read_text());archive=source/'cuts/aam.pkl.gz'
    start=time.perf_counter()
    if archive.exists():
        aam=read_aam_checkpoint(archive)
        e=evaluate_planned(aam,plan,ref['features'],ref['mapping'],seconds=220,query_timeout_ms=5000)
        e.update(search_incomplete=False,archive=str(archive),total_seconds=time.perf_counter()-start)
        save(out/'evaluation.json',e)
        return
    from rxn_core.search_symmetry import finalize_graph_symmetry
    from rxn_core.frag import build_graph
    from rxn_core.domain import AAMResult,AAMSearchMetrics
    target=build_graph(plan.problem.product.elements,plan.problem.product.wbo,plan.config.graph_floor)
    checks=[]
    for cut in raw_cut_paths(source/'cuts'):
        remaining=230-(time.perf_counter()-start)
        if remaining<=0:break
        graph=read_raw_cut(cut)
        graph,_=finalize_graph_symmetry(graph,target,iso_tolerance=plan.config.iso_tolerance)
        aam=AAMResult(plan.problem,plan.config,graph,AAMSearchMetrics.from_record({},0))
        e=evaluate_planned(aam,plan,ref['features'],ref['mapping'],seconds=min(30,remaining),query_timeout_ms=3000)
        e.update(search_incomplete=True,archive=str(cut),top1_correct=None)
        checks.append(e);save(out/'cut_checks.json',checks)
        if e['reference_recovery']=='recovered':save(out/'evaluation.json',e);return
    save(out/'evaluation.json',dict(reference_recovery='unknown',search_incomplete=True,
        checked_cuts=len(checks),top1_correct=None,total_seconds=time.perf_counter()-start))


def worker(args):
    m=json.loads((args.run/'manifest.json').read_text());index=m['indices'][args.slot]
    out=args.run/str(index);out.mkdir(exist_ok=True)
    save(out/'status.json',dict(stage='scoring',job=os.environ.get('SLURM_JOB_ID'),updated=time.time()))
    code=guarded([sys.executable,__file__,'score','--run',str(args.run),'--slot',str(args.slot)],out/'score.log',300)
    if not (out/'evaluation.json').exists():save(out/'evaluation.json',dict(reference_recovery='unknown',reason=code))
    save(out/'status.json',dict(stage='complete',exit_code=code,updated=time.time()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['init','worker','score'])
    p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path);p.add_argument('--slot',type=int)
    p.add_argument('--one-sided',action='store_true')
    a=p.parse_args();dict(init=init,worker=worker,score=score)[a.mode](a)

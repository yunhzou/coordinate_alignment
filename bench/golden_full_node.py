"""Resume fixed-budget Golden searches and verify independent cuts in parallel."""
import argparse
from dataclasses import asdict
import gc
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
import time

from golden_policy_campaign import load_case, save, guarded
from golden_evaluation import evaluate_planned
from rxn_core import search_aam
from rxn_core.aam import _initialize_finalization, _restore_finalized_cut
from rxn_core.alignment.sweep import cut_sweep_items
from rxn_core.domain import AAMResult, AAMSearchMetrics
from rxn_core.artifacts import raw_cut_paths


def aggregate(checks, expected):
    if any(e['reference_recovery']=='recovered' for e in checks):return 'recovered'
    if len(checks)==expected and all(e['reference_recovery']=='not_recovered' for e in checks):return 'not_recovered'
    return 'unknown'


def init(args):
    original=json.loads((args.source/'manifest.json').read_text())
    summary=json.loads(args.summary.read_text())
    args.run.mkdir(parents=True,exist_ok=False)
    for index in summary['unknown']:
        src=args.source/str(index);out=args.run/str(index);out.mkdir()
        for name in ('input.json','reference.json','orientation.json'):shutil.copy2(src/name,out/name)
        # Immutable checkpoints are shared; atomic replacement keeps source files intact.
        shutil.copytree(src/'cuts',out/'cuts',copy_function=os.link,
                        ignore=shutil.ignore_patterns('*.tmp'))
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    save(args.run/'manifest.json',dict(source=str(args.run.resolve()),baseline_source=str(args.source.resolve()),
        indices=summary['unknown'],config=original['config'],workers=48,
        baseline_summary={'all_records':{'recovered':summary['union_recovered'],
            'total':original['baseline_summary']['all_records']['total']}},
        search_watchdog=300,score_watchdog=240,
        scope='Same search configuration; resume saved cuts with 48 workers and parallel independent-cut verification. No global top-one ranking.',
        source_hashes={str(p.relative_to(args.run/'engine')):hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (args.run/'engine').rglob('*.py')}))


def search(args):
    m=json.loads((args.run/'manifest.json').read_text());index=m['indices'][args.slot]
    out,plan=load_case(args.run,index)
    if (out/'cuts/aam.pkl.gz').exists():
        save(out/'search.json',dict(reused_complete_archive=True,seconds=0,direction=plan.direction));return
    start=time.perf_counter()
    result=search_aam(plan.problem,plan.config,workers=48,intermediate_dir=out/'cuts',
                      resume=True,archive_format='checkpoint')
    save(out/'search.json',dict(seconds=time.perf_counter()-start,metrics=asdict(result.metrics),
                                direction=plan.direction,resumed=True))


def setup(run,index):
    global _PLAN,_REF,_OUT
    gc.disable()
    _OUT,_PLAN=load_case(Path(run),index)
    _REF=json.loads((_OUT/'reference.json').read_text())
    _initialize_finalization(_PLAN.problem,_PLAN.config)


def check_cut(raw):
    raw=Path(raw);number=int(raw.name.split('_')[1].split('.')[0])
    name=f'cut_{number:05d}'
    finalized=raw.with_name(name+'.finalized.pkl.gz')
    _,graph,_=_restore_finalized_cut((number,str(raw),str(finalized)))
    aam=AAMResult(_PLAN.problem,_PLAN.config,graph,AAMSearchMetrics.from_record({},0))
    e=evaluate_planned(aam,_PLAN,_REF['features'],_REF['mapping'],seconds=200,query_timeout_ms=5000)
    e.update(archive=str(finalized),cut_index=number,top1_correct=None,
             witness_scope='terminal and path IDs are local to this cut archive')
    save(_OUT/'checks'/f'{name}.json',e)
    return e


def score(args):
    m=json.loads((args.run/'manifest.json').read_text());index=m['indices'][args.slot]
    out,plan=load_case(args.run,index);(out/'checks').mkdir(exist_ok=True)
    raw=raw_cut_paths(out/'cuts')
    if not raw:return
    with mp.get_context('fork').Pool(min(48,len(raw)),initializer=setup,
                                    initargs=(str(args.run),index)) as pool:
        for e in pool.imap_unordered(check_cut,map(str,raw),chunksize=1):
            if e['reference_recovery']=='recovered':break


def finish(out,plan):
    expected=len(cut_sweep_items(plan.problem.reactant.wbo,plan.config.cut_floor))
    checks=[json.loads(p.read_text()) for p in sorted((out/'checks').glob('cut_*.json'))]
    outcome=aggregate(checks,expected)
    result=next((dict(e) for e in checks if e['reference_recovery']=='recovered'),{})
    result.update(reference_recovery=outcome,top1_correct=None,expected_cuts=expected,
        checked_cuts=len(checks),search_incomplete=len(raw_cut_paths(out/'cuts'))!=expected,
        verification_incomplete=len(checks)!=expected or any(e['reference_recovery']=='unknown' for e in checks))
    save(out/'evaluation.json',result)


def worker(args):
    m=json.loads((args.run/'manifest.json').read_text());out,plan=load_case(args.run,m['indices'][args.slot])
    start=time.time()
    for phase in ('search','score'):
        save(out/'status.json',dict(stage=phase,started=start,job=os.environ.get('SLURM_JOB_ID')))
        code=guarded([sys.executable,__file__,phase,'--run',str(args.run),'--slot',str(args.slot)],
                     out/f'{phase}.log',m[f'{phase}_watchdog'])
        save(out/f'{phase}_status.json',dict(exit_code=code,finished=time.time()))
    finish(out,plan)
    save(out/'status.json',dict(stage='complete',started=start,finished=time.time()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['init','search','score','worker'])
    p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path)
    p.add_argument('--summary',type=Path);p.add_argument('--slot',type=int)
    a=p.parse_args();dict(init=init,search=search,score=score,worker=worker)[a.mode](a)

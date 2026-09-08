"""Matched-budget random versus distance-biased seed experiment."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import statistics
import sys
from types import SimpleNamespace

from golden_policy_campaign import save
from golden_more_seeds import init as initialize_policy


def init(args):
    args.run.mkdir(parents=True)
    baseline=json.loads(args.summary.read_text())
    misses=sorted(set(baseline['not_recovered']+baseline['unknown']))
    # Fixed index-spaced controls; never inspect their reference mappings.
    controls=[i for i in range(0,1851,150) if i not in misses]
    selection=dict(baseline,not_recovered=sorted(misses+controls))
    save(args.run/'selection.json',selection)
    save(args.run/'design.json',dict(misses=misses,controls=controls,
        seeds=10,cap=100,workers_per_policy=8,
        note='Paired policies run concurrently on one 16-CPU allocation. Same orientation, tolerance and sweep; references are used only for evaluation. Diagnostic subset, not a full-benchmark policy score.'))
    for policy in ('random','distance'):
        initialize_policy(SimpleNamespace(run=args.run/policy,source=args.source,
            summary=args.run/'selection.json',seeds=10,cap=100,seed_selection=policy,
            selection=None,workers=8))


def worker(args):
    children=[]
    for policy in ('random','distance'):
        run=(args.run/policy).resolve()
        env=dict(os.environ,PYTHONPATH=f'{run}/engine/src:{run}/engine/bench',
            RXN_CORE_NATIVE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        children.append(subprocess.Popen([sys.executable,str(run/'engine/bench/golden_more_seeds.py'),
            'worker','--run',str(run),'--slot',str(args.slot)],env=env))
    codes=[p.wait() for p in children]
    if any(codes):raise SystemExit(f'Worker exits: {codes}')


def report(args):
    design=json.loads((args.run/'design.json').read_text());rows=[]
    for index in sorted(design['misses']+design['controls']):
        row=dict(index=index,group='miss' if index in design['misses'] else 'control')
        for policy in ('random','distance'):
            root=args.run/policy/str(index)
            def read(name):
                path=root/name
                return json.loads(path.read_text()) if path.exists() else {}
            row[policy]=dict(evaluation=read('evaluation.json'),search=read('search.json'),
                status=read('status.json'))
        rows.append(row)
    args.output.mkdir(parents=True,exist_ok=True)
    save(args.output/'cases.json',rows);save(args.output/'design.json',design)
    outcomes={policy:{group:{outcome:sum(
        r['group']==group and r[policy]['evaluation'].get('reference_recovery','pending')==outcome for r in rows)
        for outcome in ('recovered','not_recovered','unknown','pending')} for group in ('miss','control')}
        for policy in ('random','distance')}
    paired=[r for r in rows if all(r[p]['search'] for p in ('random','distance'))]
    summary=dict(outcomes=outcomes,
        distance_only=[r['index'] for r in rows if r['distance']['evaluation'].get('reference_recovery')=='recovered'
            and r['random']['evaluation'].get('reference_recovery')=='not_recovered'],
        random_only=[r['index'] for r in rows if r['random']['evaluation'].get('reference_recovery')=='recovered'
            and r['distance']['evaluation'].get('reference_recovery')=='not_recovered'],
        paired_completed_searches=len(paired),
        paired_search_seconds={p:sum(r[p]['search']['seconds'] for r in paired) for p in ('random','distance')},
        paired_median_search_ratio=statistics.median(r['distance']['search']['seconds']/r['random']['search']['seconds']
            for r in paired) if paired else None)
    save(args.output/'summary.json',summary)
    for policy in ('random','distance'):
        save(args.output/f'{policy}_manifest.json',json.loads((args.run/policy/'manifest.json').read_text()))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['init','worker','report']);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--source',type=Path);p.add_argument('--summary',type=Path)
    p.add_argument('--slot',type=int);p.add_argument('--output',type=Path)
    args=p.parse_args();globals()[args.mode](args)

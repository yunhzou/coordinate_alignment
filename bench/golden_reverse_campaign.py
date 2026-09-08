"""Opposite-direction campaign over the fixed 21 unresolved Golden records."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

from golden_policy_campaign import load_case, save, guarded
from view_golden_remaining import INDICES


def initialize(args):
    args.run.mkdir(parents=True, exist_ok=False)
    (args.run/'cases').mkdir()
    (args.run/'status').mkdir()
    for folder in ('src','bench'):
        shutil.copytree(folder, args.run/'engine'/folder,
                        ignore=shutil.ignore_patterns('__pycache__'))
    rows = []
    for index in INDICES:
        _, plan = load_case(args.source, index)
        rows.append(dict(index=index, original_direction=plan.direction,
            direction='R_to_P' if plan.reversed else 'P_to_R',
            input_R_atoms=plan.input_problem.source_atom_count,
            input_P_atoms=plan.input_problem.target_atom_count, reused=index == 590))
    previous = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden590_reverse_20260908')
    design = json.loads((previous/'design.json').read_text())
    assert design['direction'] == 'P_to_R' and design['config']['seed_count'] == 10 and design['config']['branch_limit'] == 100
    assert (previous/'evaluation.json').exists()
    (args.run/'cases/590').symlink_to(previous, target_is_directory=True)
    manifest = dict(source=str(args.source.resolve()), cases=rows, seeds=10, cap=100,
        workers=args.workers, search_watchdog=300, score_watchdog=260,
        baseline_recovered=1830, benchmark_total=1851,
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
        note='Fixed unresolved cohort; opposite search direction only. Reuse case 590. '
             'Diagnostic recovery union, not a new full-dataset policy score.')
    save(args.run/'manifest.json',manifest)


def worker(args):
    m = json.loads((args.run/'manifest.json').read_text())
    row = [r for r in m['cases'] if not r['reused']][args.slot]
    index = row['index']
    output = args.run/'cases'/str(index)
    status_file = args.run/'status'/f'{index}.json'
    start = time.time()
    status = dict(index=index, job=os.environ.get('SLURM_JOB_ID'), started=start,
                  direction=row['direction'], phase='search')
    save(status_file,status)
    common = ['--source',m['source'],'--output',str(output),'--index',str(index),
              '--direction',row['direction'],'--seeds',str(m['seeds']),
              '--cap',str(m['cap']),'--workers',str(m['workers'])]
    script = str(args.run/'engine/bench/golden_forced_direction.py')
    status['search_exit'] = guarded([sys.executable,script,'search',*common],
        args.run/'status'/f'{index}.search.log',m['search_watchdog'])
    status['phase'] = 'score'
    save(status_file,status)
    if (output/'design.json').exists():
        status['score_exit'] = guarded([sys.executable,script,'score',*common],
            args.run/'status'/f'{index}.score.log',m['score_watchdog'])
        if not (output/'evaluation.json').exists():
            save(output/'evaluation.json',dict(reference_recovery='unknown',
                reason=f"Score process exit: {status['score_exit']}",
                search_incomplete=not (output/'cuts/aam.pkl.gz').exists()))
    status.update(phase='finished',finished=time.time(),wall_seconds=time.time()-start)
    save(status_file,status)


def submit(args):
    m = json.loads((args.run/'manifest.json').read_text())
    count = sum(not r['reused'] for r in m['cases'])
    engine = args.run.resolve()/'engine'
    command = shlex.join(['env',f'PYTHONPATH={engine}/src:{engine}/bench',
        'RXN_CORE_NATIVE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1',
        sys.executable,str(engine/'bench/golden_reverse_campaign.py'),'worker',
        '--run',str(args.run.resolve()),'--slot'])+' "$SLURM_ARRAY_TASK_ID"'
    job = subprocess.check_output(['sbatch','--parsable','--partition=cpunodes',
        f'--cpus-per-task={m["workers"]}','--mem=48G','--time=00:10:00',
        '--job-name=gold_reverse21',f'--array=0-{count-1}%10',
        f'--output={args.run.resolve()}/status/slurm_%A_%a.log','--wrap',command],text=True).strip()
    save(args.run/'job.json',dict(job=job,submitted=time.time()))
    print(job,flush=True)


def report(args):
    m = json.loads((args.run/'manifest.json').read_text())
    rows = []
    for row in m['cases']:
        output = args.run/'cases'/str(row['index'])
        def read(p):
            return json.loads(p.read_text()) if p.exists() else {}
        evaluation = read(output/'evaluation.json')
        status = read(args.run/'status'/f'{row["index"]}.json')
        rows.append(dict(**row,evaluation=evaluation,search=read(output/'search.json'),status=status,
            outcome=evaluation.get('reference_recovery','unknown' if status.get('phase') == 'finished' else 'pending'),
            output=str(output.resolve())))
    recovered = [r['index'] for r in rows if r['outcome'] == 'recovered']
    summary = dict(outcomes=dict(Counter(r['outcome'] for r in rows)),newly_recovered=recovered,
        not_recovered=[r['index'] for r in rows if r['outcome'] == 'not_recovered'],
        unknown=[r['index'] for r in rows if r['outcome'] == 'unknown'],
        pending=[r['index'] for r in rows if r['outcome'] == 'pending'],
        baseline_recovered=m['baseline_recovered'],diagnostic_union=m['baseline_recovered']+len(recovered),
        benchmark_total=m['benchmark_total'],
        diagnostic_union_percent=100*(m['baseline_recovered']+len(recovered))/m['benchmark_total'],
        capped_cases=[r['index'] for r in rows if r['search'].get('capped')],
        search_watchdog_cases=[r['index'] for r in rows if r['status'].get('search_exit') == 'timeout'],
        score_watchdog_cases=[r['index'] for r in rows if r['status'].get('score_exit') == 'timeout'],
        note='Only certified positives add to the earlier recovery union. Capped negatives are not exhaustive.')
    save(args.run/'cases.json',rows)
    save(args.run/'summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('initialize','submit','worker','report'))
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--source',type=Path)
    p.add_argument('--workers',type=int,default=16)
    p.add_argument('--slot',type=int)
    a = p.parse_args()
    globals()[a.mode](a)

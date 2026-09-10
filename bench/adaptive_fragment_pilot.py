"""Bounded adaptive closure ablations; all graph snapshots survive evaluation."""
import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import sys
import time


def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2)+'\n')


def prepare(args):
    root=Path(__file__).resolve().parents[1]
    for folder in ('src','native','bench','tests','tools','benchmarks','docs/example_runs'):
        shutil.copytree(root/folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    source=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/conditioned_reuse_full_20260910_0o09px')
    (args.run/'inputs').mkdir();(args.run/'status').mkdir()
    hashes={}
    for i in (25,76,77,114):
        path=args.run/f'inputs/{i}.json';shutil.copy2(source/f'inputs/{i}.json',path)
        hashes[i]=hashlib.sha256(path.read_bytes()).hexdigest()
    tasks=[dict(index=i,direction=d,policy=p) for i,d in
        ((25,'R_to_P'),(76,'P_to_R'),(77,'P_to_R'),(114,'R_to_P'))
        for p in args.policies]
    save(args.run/'tasks.json',tasks)
    save(args.run/'manifest.json',dict(parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        source=str(source),input_sha256=hashes,branch_cap=100,seeds=1,iso_tolerance=1.,explicit_H=True,
        root_seed=42,work_budgets=args.work_budgets,soft_search_seconds=args.search_seconds,policy='No event-score pruning',
        watchdog=dict(worker_seconds=300,kill_after_seconds=5,slurm_minutes=10,requeue=False),
        native_sha256=hashlib.sha256(next((args.run/'engine/src/rxn_core').glob('_engine*.so')).read_bytes()).hexdigest()))


def submit(args):
    n=len(json.loads((args.run/'tasks.json').read_text()))
    command=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','PYTHONHASHSEED=0',
        f'PYTHONPATH={args.run}/engine/src:{args.run}/engine/bench','timeout','--kill-after=5s','300',
        sys.executable,str(args.run/'engine/bench/adaptive_fragment_pilot.py'),'worker','--run',str(args.run),'--slot']
    options=['sbatch','--parsable','--partition=cpunodes','--nodelist=bosque6','--nodes=1','--cpus-per-task=1',
        '--mem=8G','--time=00:10:00','--no-requeue',f'--array=0-{n-1}', '--job-name=adaptive_closure',
        f'--output={args.run}/status/%A_%a.out','--wrap',shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip();save(args.run/'submission.json',dict(job=job,command=options));print(job)


def worker(args):
    import numpy as np
    from rxn_core import AAMProblem,AAMSearchConfig
    from rxn_core.domain import MolecularEndpoint
    from rxn_core.adaptive_search import AdaptiveFragmentSearch
    from rxn_core.artifacts import write_aam_checkpoint
    from compare_elementary_outputs import event_counts
    spec=json.loads((args.run/'tasks.json').read_text())[args.slot]
    manifest=json.loads((args.run/'manifest.json').read_text())
    raw=json.loads((args.run/f"inputs/{spec['index']}.json").read_text())
    endpoints=[MolecularEndpoint(**{k:v for k,v in raw[n].items() if k in ('elements','coordinates','wbo')})
               for n in ('reactant','product')]
    if spec['direction']=='P_to_R':endpoints.reverse()
    problem=AAMProblem(*endpoints,name=raw['name'])
    folder=args.run/f'results/{args.slot}';folder.mkdir(parents=True,exist_ok=False)
    cpu=time.process_time();wall=time.perf_counter()
    session=AdaptiveFragmentSearch(problem,AAMSearchConfig(seed_count=1,branch_limit=100,iso_tolerance=1.),policy=spec['policy'])
    compute_cpu=time.process_time()-cpu;compute_wall=time.perf_counter()-wall
    rows=[]
    def snapshot(label):
        nonlocal compute_cpu,compute_wall
        cpu=time.process_time();wall=time.perf_counter();result=session.snapshot()
        compute_cpu+=time.process_time()-cpu;compute_wall+=time.perf_counter()-wall
        t=time.process_time();write_aam_checkpoint(result.aam,folder/f'{label}.pkl.gz');io=time.process_time()-t
        save(folder/f'{label}_pending.json',result.pending)
        vectors=sorted({tuple(dict(result.aam.graph.states[t].mapping)[i] for i in range(problem.atom_count))
            for t in result.aam.graph.terminals if len(result.aam.graph.states[t].mapping)==problem.atom_count})
        t=time.process_time();scores=event_counts(problem.reactant.wbo,problem.product.wbo,vectors) if vectors else []
        scoring=time.process_time()-t
        row=dict(label=label,work=result.work,growth=result.growth_calls,closures=result.closure_calls,
            reused_states=result.reused_states,pending=len(result.pending),exhausted=result.exhausted,
            states=len(result.aam.graph.states),terminals=len(result.aam.graph.terminals),capped=result.aam.graph.capped,
            best_representative=min((int(sum(v)) for v in scores),default=None),witnesses=len(vectors),
            compute_cpu=compute_cpu,compute_wall=compute_wall,checkpoint_cpu=io,scoring_cpu=scoring,
            max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        save(folder/f'{label}_witnesses.json',dict(mappings=vectors,events=[list(map(int,s)) for s in scores]))
        rows.append(row);save(folder/'summary.json',dict(**spec,rows=rows,host=os.uname().nodename,
            scope='Search+setup+symmetry CPU/wall exclude persistence/evaluation. Scores are representatives, not family minima.'))
        print(json.dumps(row),flush=True)
    for budget in manifest['work_budgets']:
        while session.work<budget and session.agenda and compute_wall<manifest['soft_search_seconds']:
            cpu=time.process_time();wall=time.perf_counter();session.advance()
            compute_cpu+=time.process_time()-cpu;compute_wall+=time.perf_counter()-wall
            if session.work%25==0:save(folder/'progress.json',dict(work=session.work,cpu=compute_cpu,wall=compute_wall,pending=len(session.agenda)))
        snapshot(f'work_{session.work:05d}')
        if not session.agenda or compute_wall>=manifest['soft_search_seconds']:break


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','submit','worker'))
    p.add_argument('--run',type=Path,required=True);p.add_argument('--slot',type=int)
    p.add_argument('--work-budgets',type=int,nargs='+',default=[100,400,1600])
    p.add_argument('--search-seconds',type=float,default=30.)
    p.add_argument('--policies',nargs='+',default=['largest_first','smallest_first'],
                   choices=['largest_first','smallest_first','event_guided','fair_depth'])
    args=p.parse_args();globals()[args.command](args)

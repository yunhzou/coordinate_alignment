"""Public-API integration check; every raw cut and final AAM result is saved."""
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
    args.run.mkdir(exist_ok=True,parents=True)
    for name in ('src','bench','tests','tools','benchmarks','docs/example_runs'):
        shutil.copytree(root/name,args.run/'engine'/name,ignore=shutil.ignore_patterns('__pycache__'))
    (args.run/'inputs').mkdir()
    source=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_tol1_20260909/inputs')
    for i in (77,114):shutil.copy2(source/f'{i}/input.json',args.run/f'inputs/{i}.json')
    save(args.run/'tasks.json',[dict(index=77,direction='P_to_R'),dict(index=114,direction='R_to_P')])
    save(args.run/'manifest.json',dict(parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        native_sha256=hashlib.sha256(next((args.run/'engine/src/rxn_core').glob('_engine*.so')).read_bytes()).hexdigest(),
        seeds=10,branch_cap=100,iso_tolerance=1.,workers=10,explicit_H=True))
    (args.run/'status').mkdir()


def submit(args):
    command=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','PYTHONHASHSEED=0',
        f'PYTHONPATH={args.run}/engine/src:{args.run}/engine/bench',sys.executable,
        str(args.run/'engine/bench/public_reuse_pilot.py'),'paired','--run',str(args.run),'--slot']
    options=['sbatch','--parsable','--partition=cpunodes','--nodes=1','--cpus-per-task=10','--mem=96G',
        '--time=00:10:00','--array=0-1','--job-name=public_reuse',f'--output={args.run}/status/%A_%a.out',
        '--wrap',shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=options));print(job)


def paired(args):
    modes=['reference','reused_native']
    if args.slot%2:modes.reverse()
    for mode in modes:
        subprocess.run(['timeout','--kill-after=5s','300',sys.executable,__file__,'worker',
            '--run',str(args.run),'--slot',str(args.slot),'--mode',mode],check=True)


def submit_tests(args):
    command=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',
        f'PYTHONPATH={args.run}/engine/src:{args.run}/engine/bench',
        'timeout','--kill-after=5s','300',sys.executable,'-m','pytest','-q',
        str(args.run/'engine/tests'),f'--junitxml={args.run}/{args.test_label}.xml']
    options=['sbatch','--parsable','--partition=cpunodes','--nodes=1','--cpus-per-task=1',
        '--mem=8G','--time=00:10:00','--job-name=reuse_tests',f'--output={args.run}/{args.test_label}.out',
        '--wrap',shlex.join(command)]
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/f'{args.test_label}_submission.json',dict(job=job,command=options));print(job)


def worker(args):
    from rxn_core import AAMProblem,AAMSearchConfig,search_aam
    from rxn_core.domain import MolecularEndpoint
    import rxn_core.artifacts as artifacts
    spec=json.loads((args.run/'tasks.json').read_text())[args.slot]
    raw=json.loads((args.run/f"inputs/{spec['index']}.json").read_text())
    endpoints=[MolecularEndpoint(**{k:v for k,v in raw[name].items() if k in ('elements','coordinates','wbo')})
               for name in ('reactant','product')]
    if spec['direction']=='P_to_R':endpoints.reverse()
    problem=AAMProblem(*endpoints,name=raw['name'])
    folder=args.run/f'results/{args.slot}/{args.mode}';folder.mkdir(parents=True,exist_ok=False)
    # Disjoint persistence entry points, timed inside each producing process.
    for name in ('write_raw_cut','write_graph_checkpoint','write_aam_checkpoint'):
        original=getattr(artifacts,name)
        def timed(*a,_fn=original,_name=name,**kw):
            cpu=time.process_time();wall=time.perf_counter();result=_fn(*a,**kw)
            row=dict(function=_name,cpu=time.process_time()-cpu,wall=time.perf_counter()-wall)
            with (folder/f'io_{os.getpid()}.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
            return result
        setattr(artifacts,name,timed)
    def cpu():
        return sum(getattr(resource.getrusage(who),axis) for who in (resource.RUSAGE_SELF,resource.RUSAGE_CHILDREN)
                   for axis in ('ru_utime','ru_stime'))
    start_cpu=cpu();started=time.perf_counter()
    result=search_aam(problem,AAMSearchConfig(seed_count=10,branch_limit=100,iso_tolerance=1.),
                      workers=10,execution=args.mode,intermediate_dir=folder/'aam',archive_format='checkpoint')
    elapsed=time.perf_counter()-started;total_cpu=cpu()-start_cpu
    io=[json.loads(line) for p in folder.glob('io_*.jsonl') for line in p.read_text().splitlines()]
    save(folder/'timing.json',dict(**spec,execution=args.mode,name=raw['name'],host=os.uname().nodename,
        cpu_including_persistence=total_cpu,checkpoint_cpu=sum(r['cpu'] for r in io),
        compute_cpu_excluding_checkpoint_writes=total_cpu-sum(r['cpu'] for r in io),
        elapsed_including_persistence=elapsed,metrics=asdict(result.metrics),
        scope='Parent+all worker CPU. Excludes measured checkpoint writes; includes restoration, IPC and graph merging. Actual elapsed includes persistence.'))
    digest=hashlib.sha256()
    for chunk in json.JSONEncoder(sort_keys=True,separators=(',',':')).iterencode(result.graph.to_record(copy=False)):
        digest.update(chunk.encode())
    save(folder/'validation.json',dict(digest=digest.hexdigest(),states=len(result.graph.states),
        terminals=len(result.graph.terminals),capped=result.graph.capped))
    print(args.mode,'DONE',total_cpu,elapsed,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','submit','paired','worker','submit_tests'))
    p.add_argument('--run',type=Path,required=True);p.add_argument('--slot',type=int)
    p.add_argument('--test-label',default='tests')
    p.add_argument('--mode',choices=('reference','reused_native'))
    args=p.parse_args();globals()[args.command](args)

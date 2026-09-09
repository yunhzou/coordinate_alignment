"""Six-vertex graph-only oracle: identical binary inputs, exact optimum known.

The 720 permutations are used solely as a tiny test oracle, never in production
mapping or holdout family extraction. These are graphs, not chemical reactions.
"""
import argparse
from itertools import permutations
import json
from pathlib import Path
import time
import shutil
import shlex
import subprocess
import sys
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor
import numpy as np


def save(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def generate(args):
    from slapmapper.core import SlapMapper,LabeledGraph
    args.run.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(20260909);n=6;perms=np.array(list(permutations(range(n))))
    a,b=np.triu_indices(n,1)
    def graph():
        while True:
            w=np.zeros((n,n),dtype=int)
            for i in range(1,n):j=int(rng.integers(i));w[i,j]=w[j,i]=1
            for i,j in zip(a,b):
                if w[i].sum()<4 and w[j].sum()<4 and rng.random()<.18:w[i,j]=w[j,i]=1
            if max(w.sum(axis=1))<=4:return w
    for index in range(100):
        r,p=graph(),graph();costs=np.sum(r[a,b]!=p[perms[:,a],perms[:,b]],axis=1)
        optimal=int(costs.min());oracle=perms[int(np.argmin(costs))].tolist()
        pair=[LabeledGraph({i:{j:int(w[i,j]) for j in range(n) if w[i,j]} for i in range(n)},[6]*n) for w in (r,p)]
        mapper=SlapMapper(binary=True);start=time.perf_counter();mapper.get_maps(pair,break_sym_targets=list(range(n)))
        elapsed=time.perf_counter()-start;results=[]
        for result in mapper.results:
            l,rlabels=result['lgp'];mapping={i:rlabels.label2idxs[label][0] for i,label in enumerate(l.labels)}
            assert all(len(v)==1 for v in rlabels.label2idxs.values())
            m=np.array([mapping[i] for i in range(n)])
            results.append(dict(mapping=m.tolist(),edits=int(np.sum(r[a,b]!=p[m[a],m[b]])),native_cost=float(result['val'])))
        save(args.run/f'{index}.json',dict(index=index,R=r.tolist(),P=p.tolist(),exact_min=optimal,
            oracle_witness=oracle,slap=results,slap_seconds=elapsed))
    print('100 fixed-seed six-vertex graph pairs saved')


def aam(args):
    from rxn_core import AAMProblem,AAMSearchConfig,search_aam
    from rxn_core.domain import MolecularEndpoint
    row=json.loads((args.run/f'{args.index}.json').read_text());n=len(row['R'])
    problem=AAMProblem(*(MolecularEndpoint(('C',)*n,np.zeros((n,3)),row[k]) for k in ('R','P')))
    out=args.run/f'aam_{args.index}';out.mkdir(exist_ok=False)
    start=time.perf_counter();result=search_aam(problem,AAMSearchConfig(seed_count=10,branch_limit=100,iso_tolerance=.5),
        workers=1,intermediate_dir=out/'cuts',archive_format='checkpoint')
    r,p=np.array(row['R']),np.array(row['P']);a,b=np.triu_indices(n,1);best=None;witness=None
    for terminal in result.graph.terminals:
        m=dict(result.graph.states[terminal].mapping)
        if len(m)!=n:continue
        vector=np.array([m[i] for i in range(n)]);cost=int(np.sum(r[a,b]!=p[vector[a],vector[b]]))
        if best is None or cost<best:best=cost;witness=vector.tolist()
    save(out/'summary.json',dict(best_saved_events=best,mapping=witness,capped=result.graph.capped,
        terminals=len(result.graph.terminals),elapsed=time.perf_counter()-start,exact_min=row['exact_min']))


def aam_all(args):
    tasks=[SimpleNamespace(run=args.run,index=i) for i in range(100)]
    with ProcessPoolExecutor(max_workers=16) as pool:list(pool.map(aam,tasks))


def expand_slap(args):
    from slapmapper.core import SlapMapper,LabeledGraph
    rng=np.random.default_rng(20260909);results=[]
    for index in range(100):
        row=json.loads((args.run/f'{index}.json').read_text());r,p=np.array(row['R']),np.array(row['P']);n=len(r)
        a,b=np.triu_indices(n,1);best=10**9;attempts=[]
        for order in range(10):
            ro=np.arange(n) if order==0 else rng.permutation(n)
            po=np.arange(n) if order==0 else rng.permutation(n)
            for reverse in (False,True):
                pair=[LabeledGraph({i:{j:int(w[i,j]) for j in range(n) if w[i,j]} for i in range(n)},[6]*n)
                      for w in (r[np.ix_(ro,ro)],p[np.ix_(po,po)])]
                mapper=SlapMapper(binary=True);mapper.get_maps(pair[::-1] if reverse else pair,break_sym_targets=list(range(n)))
                found=[]
                for result in mapper.results:
                    left,right=result['lgp'];m={i:right.label2idxs[label][0] for i,label in enumerate(left.labels)}
                    if reverse:m={v:k for k,v in m.items()}
                    mapped={int(ro[i]):int(po[j]) for i,j in m.items()};vector=np.array([mapped[i] for i in range(n)])
                    cost=int(np.sum(r[a,b]!=p[vector[a],vector[b]]));best=min(best,cost)
                    found.append(dict(mapping=vector.tolist(),events=cost))
                attempts.append(dict(order=order,reverse=reverse,results=found))
        results.append(dict(index=index,exact_min=row['exact_min'],best=best,attempts=attempts))
    save(args.run/'slap_expanded.json',results)
    print('Expanded SLAP recovered exact optimum:',sum(r['best']==r['exact_min'] for r in results),'/100')


def submit(args):
    shutil.copy2(__file__,args.run/'oracle_driver.py')
    command=['env','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','PYTHONHASHSEED=0','RXN_CORE_NATIVE=1',
        'timeout','--kill-after=5s','300',sys.executable,str(args.run/'oracle_driver.py'),'aam_all','--run',str(args.run)]
    options=['sbatch','--parsable','--partition=cpunodes','--exclude=bosque5,bosque6,bosque8',
        '--cpus-per-task=16','--mem=32G','--time=00:10:00','--job-name=aam_oracle',
        f'--output={args.run}/aam_%j.out','--wrap',shlex.join(command)]
    job=subprocess.check_output(options,text=True).strip();save(args.run/'submission.json',dict(job=job,command=options));print(job)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=('generate','aam','aam_all','submit','expand_slap'))
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--index',type=int)
    args=parser.parse_args();globals()[args.command](args)

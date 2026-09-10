"""Same-node default single-seed regression check across two frozen engines."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

from adaptive_fragment_pilot import save


def worker(args):
    from rxn_core import AAMProblem, AAMSearchConfig
    from rxn_core.domain import MolecularEndpoint
    import rxn_core.aam as aam
    from rxn_core.artifacts import write_graph_checkpoint
    from rxn_core.conditioned_symmetry import ConditionedSymmetryWorkspace
    from rxn_core.search_symmetry import finalize_graph_symmetry

    rows=[]
    for repeat in range(3):
        for index,direction in ((25,'R_to_P'),(76,'P_to_R'),(77,'P_to_R'),(114,'R_to_P')):
            raw=json.loads((args.inputs/f'{index}.json').read_text())
            endpoints=[MolecularEndpoint(**{k:v for k,v in raw[name].items()
                if k in ('elements','coordinates','wbo')}) for name in ('reactant','product')]
            if direction=='P_to_R':endpoints.reverse()
            problem=AAMProblem(*endpoints,name=raw['name'])
            config=AAMSearchConfig(seed_count=1,branch_limit=100,iso_tolerance=1.)
            cpu=time.process_time();wall=time.perf_counter()
            aam._initialize_search(problem,config,'reused_native')
            workspace=ConditionedSymmetryWorkspace(aam._SEARCH_CONTEXT[2],1.)
            graph,counts=aam._search_cut(())
            graph,_=finalize_graph_symmetry(graph,aam._SEARCH_CONTEXT[2],iso_tolerance=1.,workspace=workspace)
            elapsed=time.perf_counter()-wall;compute=time.process_time()-cpu
            payload=json.dumps(graph.to_record(copy=False),sort_keys=True,separators=(',',':')).encode()
            write_graph_checkpoint(graph,args.output/f'{index}_{repeat}.pkl.gz')
            rows.append(dict(index=index,repeat=repeat,cpu=compute,wall=elapsed,
                digest=hashlib.sha256(payload).hexdigest(),states=len(graph.states),
                terminals=len(graph.terminals),capped=graph.capped))
            save(args.output/'measurements.json',rows)
    save(args.output/'environment.json',dict(host=os.uname().nodename,affinity=list(os.sched_getaffinity(0)),
        engine_file=aam.__file__,native_sha256=hashlib.sha256(
            next(Path(aam.__file__).parent.glob('_engine*.so')).read_bytes()).hexdigest()))


def compare(args):
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    for name,source in (('baseline',args.baseline),('current',args.current)):
        output=args.output/name;output.mkdir(parents=True,exist_ok=False)
        env=dict(os.environ,PYTHONPATH=str(source),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        subprocess.run(['timeout','--kill-after=5s','300',sys.executable,__file__,'worker',
            '--inputs',str(args.inputs),'--output',str(output)],env=env,check=True)
    data={name:json.loads((args.output/name/'measurements.json').read_text()) for name in ('baseline','current')}
    rows=[]
    for index in (25,76,77,114):
        selected={name:[r for r in records if r['index']==index] for name,records in data.items()}
        baseline=statistics.median(r['cpu'] for r in selected['baseline'])
        current=statistics.median(r['cpu'] for r in selected['current'])
        rows.append(dict(index=index,baseline_median_cpu=baseline,current_median_cpu=current,
            ratio=current/baseline,exact_graph_equality=all(a['digest']==b['digest']
                for a,b in zip(selected['baseline'],selected['current']))))
    save(args.output/'comparison.json',dict(rows=rows,
        scope='Three fresh-cache trials per case/engine, one pinned CPU; setup+uncut search+symmetry. Imports, input loading, digest and persistence excluded. Small smoke test, not a full-benchmark performance claim.'))
    assert all(r['exact_graph_equality'] for r in rows)
    print(json.dumps(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('worker','compare'))
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--baseline',type=Path)
    p.add_argument('--current',type=Path)
    args=p.parse_args();globals()[args.command](args)

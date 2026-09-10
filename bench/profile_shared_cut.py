"""Profile one reference-blind cut, including cyclic-GC cost, without scoring labels."""
import argparse
import cProfile
import gc
import json
from pathlib import Path
import pstats
import time
from types import SimpleNamespace

from adaptive_full_benchmark import problem_plan, read, save
from rxn_core import aam


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--index',type=int,required=True)
    p.add_argument('--direction',choices=['R_to_P','P_to_R'],required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    spec=dict(dataset='golden',index=args.index,direction=args.direction)
    _,plan=problem_plan(SimpleNamespace(run=args.run,method='adaptive'),spec)
    aam._initialize_search(plan.problem,plan.config,'shared_policies')
    times=[]
    def observe(phase, info):
        times.append((phase,time.process_time()))
    gc.callbacks.append(observe)
    profiler=cProfile.Profile()
    start=time.process_time()
    profiler.enable()
    graph,counts=aam._search_cut(())
    profiler.disable()
    cpu=time.process_time()-start
    gc.callbacks.remove(observe)
    args.output.mkdir(parents=True,exist_ok=True)
    profiler.dump_stats(str(args.output/'profile.pstats'))
    gc_cpu=sum(b[1]-a[1] for a,b in zip(times,times[1:]) if a[0]=='start' and b[0]=='stop')
    save(args.output/'summary.json',dict(cpu=cpu,gc_cpu=gc_cpu,gc_callbacks=len(times),
        states=len(graph.states),transitions=len(graph.transitions),counts=counts))
    print(json.dumps(dict(cpu=cpu,gc_cpu=gc_cpu,states=len(graph.states))),flush=True)
    pstats.Stats(profiler).strip_dirs().sort_stats('cumulative').print_stats(25)

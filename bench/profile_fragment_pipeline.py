"""Profile unchanged seed/cut search; save every compressed graph separately.

The seed index refers to the original ten-order schedule, not a new seed policy.
cProfile instruments search only. These are diagnostic, not speed-claim timings.
"""
import argparse
import cProfile
import hashlib
import json
import os
from pathlib import Path
import pstats
import resource
import time
import numpy as np

from rxn_core.aam import cut_seed
from rxn_core.alignment.branch import _generate_seed_orders, find_islands
from rxn_core.alignment.sweep import cut_sweep_items
from rxn_core.artifacts import write_graph_checkpoint
from rxn_core.frag import build_graph
from rxn_core.matcher import _nauty_orbits
from rxn_core.search_symmetry import finalize_graph_symmetry


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def run(args):
    prefix='profiles' if args.backend=='python' else 'profiles_native'
    root=args.run/f'{prefix}/{args.index}/{args.direction}/seed_{args.seed_index:02d}'
    root.mkdir(parents=True,exist_ok=False)
    raw=json.loads((args.source/f'inputs/{args.index}/input.json').read_text())
    reactant,product=(raw[k] for k in ('reactant','product'))
    if args.direction=='P_to_R':reactant,product=product,reactant
    source=build_graph(reactant['elements'],np.asarray(reactant['wbo']),bond_cut=.2)
    target=build_graph(product['elements'],np.asarray(product['wbo']),bond_cut=.2)
    orbits=_nauty_orbits(target,wbo_tol=1.)
    cuts=cut_sweep_items(np.asarray(reactant['wbo']),.2)
    profiler=cProfile.Profile()
    matcher=find_islands
    if args.backend=='native':
        from rxn_core.native_search import find_islands_native
        matcher=find_islands_native
    records=[]
    for ordinal,cut in enumerate(cuts):
        r=source.copy();r.remove_edges_from(cut)
        r_orbits=_nauty_orbits(r,wbo_tol=1.)
        order=_generate_seed_orders(r,10,rng_seed=cut_seed(cut))[args.seed_index]
        growth=[]
        cpu=time.process_time();wall=time.perf_counter()
        profiler.enable()
        try:
            graph=matcher(r,target,order,graph_floor=.2,iso_tol=1.,
                max_branches=100,p_orbits=orbits,r_orbits=r_orbits,cuts=cut,profile=growth)
        finally:
            profiler.disable()
        timings=dict(search_cpu=time.process_time()-cpu,search_wall=time.perf_counter()-wall)
        cpu=time.process_time()
        graph,_=finalize_graph_symmetry(graph,target,iso_tolerance=1.)
        timings['symmetry_cpu']=time.process_time()-cpu
        cpu=time.process_time()
        digest=hashlib.sha256(json.dumps(graph.to_record(copy=False),sort_keys=True,separators=(',',':')).encode()).hexdigest()
        write_graph_checkpoint(graph,root/f'cut_{ordinal:04d}.pkl.gz')
        timings['audit_and_persistence_cpu']=time.process_time()-cpu
        records.append(dict(cut=cut,digest=digest,timings=timings,capped=graph.capped,
            states=len(graph.states),transitions=len(graph.transitions),terminals=len(graph.terminals),
            growth_calls=sum(p.get('growth_calls',1) for p in growth),
            native_including_conversion_wall=sum(p.get('extend_elapsed_sec',0) for p in growth),
            native_batch_phases={k:sum(p.get(k,0) for p in growth) for k in
                ('native_schedule_and_export_cpu','python_graph_cpu')}))
        save(root/'progress.json',dict(completed=len(records),total=len(cuts),last=records[-1]))
    profiler.dump_stats(root/'search.pstats')
    stats=pstats.Stats(profiler)
    functions=[dict(file=f,line=l,function=n,primitive_calls=cc,calls=nc,self_seconds=tt,cumulative_seconds=ct)
               for (f,l,n),(cc,nc,tt,ct,callers) in stats.stats.items()]
    functions.sort(key=lambda r:-r['self_seconds'])
    save(root/'summary.json',dict(index=args.index,direction=args.direction,seed_index=args.seed_index,backend=args.backend,
        name=raw['name'],host=os.uname().nodename,settings=dict(seeds=10,cap=100,match_tolerance=1.,
        score_tolerance=.5,explicit_hydrogen=True,all_cuts=True),cuts=records,functions=functions,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='cProfile-instrumented search, all cuts for one seed index. Diagnostic timing, not speedup evidence.'))
    print(json.dumps(functions[:20],indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path,default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_tol1_20260909'))
    parser.add_argument('--index',type=int,required=True)
    parser.add_argument('--direction',choices=('R_to_P','P_to_R'),required=True)
    parser.add_argument('--seed-index',type=int,default=0)
    parser.add_argument('--backend',choices=('python','native'),default='python')
    run(parser.parse_args())

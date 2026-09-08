"""Bounded single-process core-AAM memory probe; search rules are unchanged."""
import argparse
import gc
import hashlib
import json
from pathlib import Path
import pickle
import resource
import signal
import sys
import time
import tracemalloc


def prepare(args):
    from golden_policy_campaign import load_case
    from rxn_core.frag import build_graph
    from rxn_core.aam import cut_seed
    from rxn_core.alignment.branch import _generate_seed_orders
    _,plan=load_case(args.source,833)
    r=build_graph(plan.problem.reactant.elements,plan.problem.reactant.wbo,plan.config.graph_floor)
    p=build_graph(plan.problem.product.elements,plan.problem.product.wbo,plan.config.graph_floor)
    orders=list(_generate_seed_orders(r,n_trials=10,rng_seed=cut_seed(())))
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'input.pkl').write_bytes(pickle.dumps((r,p,orders),protocol=5))


def probe(args):
    sys.path.insert(0,str(args.package/'src'))
    from rxn_core.alignment.branch import find_islands
    from rxn_core.growth.native import built
    assert built(), 'Native engine must be present for a fair comparison'
    data=args.input.read_bytes();r,p,orders=pickle.loads(data)
    builders=[];results=[];profile=[];start=time.monotonic()
    if not args.legacy:
        from rxn_core.search_graph import SearchGraphBuilder
        original=SearchGraphBuilder.__init__
        def capture(self,*a,**kw):
            original(self,*a,**kw);builders.append(self)
        SearchGraphBuilder.__init__=capture
    if args.trace:tracemalloc.start(1)
    def emit(event,**extra):
        rss=int(next(l for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('VmRSS:')).split()[1])/1024
        row=dict(event=event,seconds=time.monotonic()-start,rss_mib=rss,
            peak_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            seeds_finished=len(results),profile_rows=len(profile),
            states=sum(len(b.states) for b in builders),transitions=sum(len(b.transitions) for b in builders),
            stops=sum(len(b.stops) for b in builders),partitions=sum(len(b._partitions) for b in builders),**extra)
        print(json.dumps(row),flush=True)
        return row
    class Limit(Exception):pass
    def tick(*_):
        row=emit('sample')
        if row['seconds']>=args.seconds or row['rss_mib']>=8192:raise Limit()
    signal.signal(signal.SIGALRM,tick);signal.setitimer(signal.ITIMER_REAL,5,5)
    status='complete'
    try:
        for order in orders[:args.seeds]:
            result=find_islands(r,p,order,graph_floor=.2,iso_tol=1.,max_branches=args.cap,profile=profile)
            results.append(result);emit('seed_complete')
    except Limit:status='diagnostic_limit'
    finally:signal.setitimer(signal.ITIMER_REAL,0)
    row=emit('search_end',status=status,cap=args.cap,requested_seeds=args.seeds,
             input_sha256=hashlib.sha256(data).hexdigest(),legacy=args.legacy,traced=args.trace)
    if args.fingerprint:
        fingerprints=[]
        for graph in results:
            digest=hashlib.sha256()
            for chunk in json.JSONEncoder(sort_keys=True).iterencode(graph.to_record(copy=False)):
                digest.update(chunk.encode())
            fingerprints.append(digest.hexdigest())
        row['graph_sha256']=fingerprints
    if args.trace:
        row['python_traced_bytes']=tracemalloc.get_traced_memory()
        row['allocation_sites']=[dict(site=str(s.traceback),bytes=s.size,count=s.count)
                                 for s in tracemalloc.take_snapshot().statistics('lineno')[:25]]
        tracemalloc.stop()
    if not args.legacy:
        audits=[]
        for b in builders:
            g=b.finish();live=set(g.ancestor_transitions(g.terminals))
            audits.append(dict(states=len(g.states),transitions=len(g.transitions),terminals=len(g.terminals),
                capped_stops=sum(s.reason=='capped' for s in g.stops),
                terminal_ancestor_transitions=len(live),
                fragment_transitions=sum(e.match is not None for e in g.transitions),
                join_transitions=sum(e.match is None for e in g.transitions)))
        row['graphs']=audits
    args.output.write_text(json.dumps(row,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','probe'])
    p.add_argument('--source',type=Path);p.add_argument('--input',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--package',type=Path);p.add_argument('--legacy',action='store_true');p.add_argument('--trace',action='store_true')
    p.add_argument('--fingerprint',action='store_true')
    p.add_argument('--cap',type=int,default=100);p.add_argument('--seeds',type=int,default=10);p.add_argument('--seconds',type=int,default=90)
    a=p.parse_args();dict(prepare=prepare,probe=probe)[a.mode](a)

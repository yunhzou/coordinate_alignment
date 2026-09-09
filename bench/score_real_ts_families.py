"""Common all-H event scoring of saved TS families, without mapper reruns.

AAM uses its correlated path constraints. SLAP uses final labels AND returned
LAP fingerprint counts, not arbitrary independent H-label permutations.
"""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import pickle
import time
from types import SimpleNamespace

import numpy as np
import z3

from compare_real_ts_mappings import load_slap, save
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.artifacts import read_graph_checkpoint
from rxn_core.family_query import compile_path
from rxn_core.family_scoring import bond_events, event_objective, invariant_score
from rxn_core.search_graph import frozen_value


def compile_slap(result, ro, po, problem):
    left,right=result['lgp']; n=problem.atom_count
    solver=z3.Solver(); values=[None]*n; targets=defaultdict(list); variables=[]
    for b,label in enumerate(right.labels):targets[label].append(int(po[b]))
    for a,label in enumerate(left.labels):
        source=int(ro[a]); pool=tuple(targets[label])
        if len(pool)==1:value=pool[0]
        else:
            assert problem.reactant.elements[source]=='H'
            value=z3.Int(f'h{source}'); solver.add(z3.Or(*(value==p for p in pool)))
        values[source]=(value,frozenset(pool))
        variables.append(z3.IntVal(value) if isinstance(value,int) else value)
    solver.add(z3.Distinct(*variables))
    # Retain all returned assignment correlations, including counts coupling
    # multiple H-neighborhood classes inside the same final label.
    fingerprints=0
    for solutions,infos in result['lap_sols'].values():
        source_groups=[[int(ro[a]) for a in info['idxs']] for info in infos[0].values()]
        target_groups=[[int(po[b]) for b in info['idxs']] for info in infos[1].values()]
        alternatives=[]
        for sol in solutions:
            if sol['fingerprint'] is None:
                alternatives.append(z3.BoolVal(True));continue
            fp=sol['fingerprint']; terms=[]
            for i,aa in enumerate(source_groups):
                for j,bb in enumerate(target_groups):
                    terms.append(z3.Sum([z3.If(z3.Or(*(values[a][0]==b for b in bb)),1,0)
                                         for a in aa])==int(fp[i,j]))
            alternatives.append(z3.And(*terms));fingerprints+=1
        solver.add(z3.Or(*alternatives))
    return solver,values,fingerprints


def minimize(solver,values,problem,realize,seconds):
    start=time.perf_counter(); solver.set(timeout=max(1,int(seconds*1000)))
    status=solver.check()
    if status!=z3.sat:
        return dict(status='inconsistent_native_result' if status==z3.unsat else 'unresolved',
                    lower=0,upper=None,seconds=time.perf_counter()-start)
    mapping=realize(solver.model())
    assert len(mapping)==len(set(mapping.values()))==problem.atom_count
    assert all(problem.reactant.elements[a]==problem.product.elements[b] for a,b in mapping.items())
    upper=bond_events(problem,mapping)['total']
    proxy=SimpleNamespace(problem=problem,values=values,representative=mapping)
    objective,lower,_=event_objective(proxy)
    assert lower<=upper
    while lower<upper and time.perf_counter()-start<seconds:
        mid=(lower+upper)//2;solver.push();solver.add(objective<=mid)
        solver.set(timeout=max(1,int(1000*(seconds-(time.perf_counter()-start)))))
        status=solver.check()
        if status==z3.sat:
            mapping=realize(solver.model()); upper=bond_events(problem,mapping)['total']
            assert upper<=mid
        elif status==z3.unsat:lower=mid+1
        else:solver.pop();break
        solver.pop()
    return dict(status='optimal_in_family' if lower==upper else 'bounded',lower=int(lower),upper=int(upper),
                mapping=sorted(mapping.items()),events=bond_events(problem,mapping),
                seconds=time.perf_counter()-start)


def task(args):
    spec=json.loads((args.run/'tasks.json').read_text())[args.slot]
    raw=json.loads((args.run/f"inputs/{spec['index']}.json").read_text())
    problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    method=spec['method']; folder=args.run/f"results/{spec['index']}/{method}"
    start=time.perf_counter(); cpu=time.process_time(); units=[]; seen=set(); complete=True
    if method.startswith('aam'):
        data=[dict(raw[k]) for k in ('reactant','product')]
        if method=='aam_binary':
            for d in data:d['wbo']=(np.asarray(d['wbo'])>=.2).astype(float)
        search_problem=AAMProblem(*(MolecularEndpoint(**d) for d in data))
        for filename in sorted(folder.glob('cut_*.pkl.gz')):
            if time.perf_counter()-start>120:complete=False;break
            graph=read_graph_checkpoint(filename)
            for path in graph.paths():
                if time.perf_counter()-start>120:complete=False;break
                if len(path.mapping)!=problem.atom_count:continue
                key=(tuple(sorted(path.mapping.items())),path.context.cuts,
                     tuple(frozen_value(f) for f in path.fragments))
                if key in seen:continue
                seen.add(key)
                invariant,_=invariant_score(path,problem)
                if invariant:
                    events=bond_events(problem,path.mapping)
                    result=dict(status='optimal_in_family',lower=events['total'],upper=events['total'],
                                mapping=sorted(path.mapping.items()),events=events,method='invariance')
                else:
                    compiled=compile_path(path,search_problem,{},source_atoms=(),complete_reference=False)
                    result=minimize(compiled.solver,compiled.values,problem,
                        lambda model:dict(compiled.realize(model)['mapping']),.5)
                units.append(dict(**result,graph=filename.name,terminal=path.terminal,transitions=path.transitions))
    else:
        load_slap(args.run,method)
        for filename in sorted(folder.glob('order_*.pkl.gz')):
            with gzip.open(filename,'rb') as stream:native=pickle.load(stream)
            for index,result in enumerate(native['results']):
                if time.perf_counter()-start>120:complete=False;break
                solver,values,count=compile_slap(result,native['source_order'],native['target_order'],problem)
                def realize(model):
                    return {i:(v if isinstance(v,int) else model.eval(v).as_long()) for i,(v,_) in enumerate(values)}
                scored=minimize(solver,values,problem,realize,.5)
                units.append(dict(**scored,order=filename.name,candidate=index,lap_fingerprints=count))
    feasible=[u for u in units if u.get('upper') is not None]
    best=min(feasible,key=lambda u:u['upper'],default=None)
    lower=min((u['lower'] for u in units),default=0) if complete else 0
    result=dict(**spec,complete_traversal=complete,units=len(units),
        unresolved_units=sum(u['status'] in ('bounded','unresolved') for u in units),
        inconsistent_native_units=sum(u['status']=='inconsistent_native_result' for u in units),
        lower=lower,upper=best['upper'] if best else None,
        minimum_proven_over_saved_families=best is not None and lower==best['upper'] and complete,
        best=best,postprocess_cpu=time.process_time()-cpu,
        wall_including_loading=time.perf_counter()-start,
        scope='All explicit H; event threshold .5, graph floor .2. No new mapping searches. Native LAP fingerprints retained.')
    save(folder/'full_atom_scores.json',dict(summary=result,families=units))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--slot',type=int,required=True)
    task(p.parse_args())

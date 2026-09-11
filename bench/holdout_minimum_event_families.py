"""Bidirectional membership of observed minimum-event patterns in saved families.

An UNSAT conclusion covers all relevant saved families only after a complete
scan; time-limited queries remain unresolved. No mapping search is repeated.
"""
import argparse
from collections import defaultdict
from functools import lru_cache
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import time
from types import SimpleNamespace

import pynauty
import z3
from holdout_minimum_events import (AAM,SWEEP,PYTHON,KINDS,EventPatterns,read,save,sha,
    compare_sets,AAMProblem,MolecularEndpoint,read_aam_checkpoint)
from golden_evaluation import colored_graph,exact_action
from rxn_core.family_query import compile_path,SymbolicActions
from rxn_core.family_scoring import event_objective
from rxn_core.search_graph import frozen_value


def slap_model(raw,candidate):
    problem = AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    solver = z3.Solver()
    targets,variables = defaultdict(list),defaultdict(list)
    for atom,label in enumerate(candidate['graphs'][1]['labels']):targets[label].append(atom)
    values = []
    for atom,label in enumerate(candidate['graphs'][0]['labels']):
        support = targets[label]
        if len(support)==1:value = support[0]
        else:
            assert raw['reactant']['elements'][atom]=='H'
            value = z3.FreshInt('native_h')
            solver.add(z3.Or(*(value==p for p in support)))
            variables[label].append(value)
        values.append((value,frozenset(support)))
    for group in variables.values():solver.add(z3.Distinct(*group))
    representative = {}
    left = defaultdict(list)
    for atom,label in enumerate(candidate['graphs'][0]['labels']):left[label].append(atom)
    for label,atoms in left.items():representative.update(zip(atoms,targets[label]))
    return SimpleNamespace(problem=problem,solver=solver,values=values,representative=representative)


def original_values(compiled,reverse):
    if not reverse:return compiled.values
    values = []
    for atom in range(compiled.problem.target_atom_count):
        possible = [(a,v[0]) for a,v in enumerate(compiled.values) if atom in v[1]]
        fixed = next((a for a,v in possible if isinstance(v,int) and v==atom),None)
        if fixed is not None:values.append((fixed,frozenset((fixed,))));continue
        assert possible
        expr = possible[-1][0]
        for a,v in reversed(possible[:-1]):expr = z3.If(v==atom,a,expr)
        values.append((expr,frozenset(a for a,v in possible)))
    return values


def condition_for_event(kind,va,vb,weight,product_wbo):
    allowed = defaultdict(list)
    for x in va[1]:
        for y in vb[1]:
            if x==y:continue
            p = product_wbo[x,y]
            valid = (p<=.2 if kind=='broken' else p>.2 if kind=='formed' else
                     p>.2 and p-weight>.5 if kind=='strengthened' else p>.2 and weight-p>.5)
            if valid:allowed[x].append(y)
    return z3.Or(*(z3.And(va[0]==x,z3.Or(*(vb[0]==y for y in ys))) for x,ys in allowed.items()))


def query(compiled,canonical,pattern,generators,*,reverse=False,timeout_ms=1500):
    started = time.perf_counter()
    solver = compiled.solver
    solver.push()
    encoder = SymbolicActions(solver,fresh=True)
    event_atoms = sorted({a for edges in pattern['events'].values() for edge in edges for a in edge})
    transformed = encoder.act([encoder.constant(a) for a in event_atoms],generators,'event_equivalence')
    values = dict(enumerate(original_values(compiled,reverse)))
    mapped = {a:encoder.lookup(v,values) for a,v in zip(event_atoms,transformed,strict=True)}
    for kind,edges in pattern['events'].items():
        for a,b in edges:solver.add(condition_for_event(kind,mapped[a],mapped[b],canonical.r[a,b],canonical.p))
    remaining = timeout_ms-int(1000*(time.perf_counter()-started))
    if remaining<=0:
        solver.pop()
        return 'unknown',None
    solver.set(timeout=remaining)
    status = solver.check()
    witness = None
    if status==z3.sat:
        model = solver.model()
        if hasattr(compiled,'realize'):
            mapping = dict(compiled.realize(model)['mapping'])
            if reverse:mapping = {p:r for r,p in mapping.items()}
        else:
            mapping = {a:(v[0] if isinstance(v[0],int) else model.eval(v[0]).as_long()) for a,v in enumerate(compiled.values)}
        witness = canonical.describe([mapping[a] for a in range(canonical.n)])
        assert witness['total']==pattern['total'] and witness['id']==pattern['id']
    solver.pop()
    return ('represented' if status==z3.sat else 'absent' if status==z3.unsat else 'unknown'),witness


def native_membership(raw,canonical,patterns,method,snapshot,generators,deadline):
    root = AAM if method=='native_slap' else SWEEP
    results = {k:dict(status='represented',method='saved_representative',witness=v)
               for k,v in snapshot['methods'][method]['patterns'].items()}
    minimum = snapshot['methods'][method]['minimum']
    pending = {k:v for k,v in patterns.items() if v['total']==minimum and k not in results}
    for k,v in patterns.items():
        if v['total']!=minimum:results[k] = dict(status='different_score',method='outside_method_minimum')
    if not pending:return results,dict(queries=0,encoded=0,exhausted=True)
    native = read(root/f"slap_xyz/{snapshot['index']}.json")['candidates']
    scored = read(root/f"slap_scoring/{snapshot['index']}.refined.json")['slap']
    unknown = set()
    queries,encoded = 0,0
    exhausted = True
    for ordinal,(candidate,score) in enumerate(zip(native,scored,strict=True)):
        assert score['hydrogen_score_optimization']['optimal']
        if score['events']['total']!=minimum:continue
        if not pending:break
        if time.perf_counter()>deadline:exhausted=False;break
        compiled = slap_model(raw,candidate)
        objective,lower,_ = event_objective(compiled)
        compiled.solver.add(objective==minimum)
        encoded += 1
        for key,pattern in list(pending.items()):
            remaining = int(1000*(deadline-time.perf_counter()))
            if remaining<=0:exhausted=False;break
            status,witness = query(compiled,canonical,pattern,generators,timeout_ms=min(2000,remaining))
            queries += 1
            if status=='represented':
                results[key] = dict(status=status,method='native_label_family_query',candidate=ordinal,witness=witness)
                del pending[key]
            elif status=='unknown':unknown.add(key)
        if not exhausted:break
    for key in pending:
        results[key] = dict(status='excluded_from_saved_families' if exhausted and key not in unknown else 'unresolved',
                            method='native_label_family_query')
    return results,dict(queries=queries,encoded=encoded,exhausted=exhausted)


def aam_membership(raw,canonical,patterns,snapshot,generators,deadline,
                   directions=('R_to_P','P_to_R'),known_results=None,archive_root=None):
    archive_root=AAM if archive_root is None else Path(archive_root)
    results = {k:dict(status='represented',method='saved_terminal',witness=v)
               for k,v in snapshot['methods']['aam']['patterns'].items()}
    if known_results:
        results.update(known_results)
    minimum = snapshot['methods']['aam']['minimum']
    pending = {k:v for k,v in patterns.items() if v['total']==minimum and k not in results}
    for k,v in patterns.items():
        if v['total']!=minimum:results[k] = dict(status='different_score',method='outside_method_minimum')
    if not pending:return results,dict(queries=0,encoded=0,exhausted=True)
    problem = AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    unknown = set()
    queries,encoded,paths_checked,invariant_terminals = 0,0,0,0
    exhausted = True
    for direction in directions:
        if not pending:break
        if time.perf_counter()>deadline:exhausted=False;break
        reverse = direction=='P_to_R'
        search_problem = AAMProblem(problem.product,problem.reactant) if reverse else problem
        feature = canonical.features[0 if reverse else 1]
        graph = read_aam_checkpoint(archive_root/f"results/holdout/{snapshot['index']}/{direction}/original/cuts/aam.pkl.gz").graph
        @lru_cache(None)
        def invariant_action(images):return exact_action(images,feature)
        @lru_cache(None)
        def invariant_edge(edge):
            placement = graph.fragment_placement(edge)
            if placement is None:return True
            assert placement.target_generators is not None
            if any(not invariant_action(tuple(g.images)) for g in placement.target_generators):return False
            for domain in placement.symmetry_domains:
                if domain.source=='exact_automorph_group':continue
                for a,b in zip(domain.p_atoms,domain.p_atoms[1:]):
                    images = list(range(canonical.n))
                    images[a],images[b] = b,a
                    if not invariant_action(tuple(images)):return False
            return True
        incoming = [[] for _ in graph.states]
        for edge in graph.transitions:incoming[edge.target].append(edge)
        @lru_cache(None)
        def can_change(state):
            return any(can_change(e.source) or not invariant_edge(e.id) for e in incoming[state])
        seen = set()
        for terminal in graph.terminals:
            if len(graph.states[terminal].mapping)!=canonical.n:continue
            if time.perf_counter()>deadline:exhausted=False;break
            if not can_change(terminal):invariant_terminals+=1;continue
            for path in graph.paths(terminal):
                if time.perf_counter()>deadline:exhausted=False;break
                paths_checked += 1
                if all(invariant_edge(e) for e in path.transitions):continue
                key = (tuple(sorted(path.mapping.items())),tuple(path.context.cuts),
                       path.context.graph_floor,path.context.iso_tolerance,frozen_value(path.fragments))
                if key in seen:continue
                seen.add(key)
                compiled = compile_path(path,search_problem,{},source_atoms=(),complete_reference=False)
                objective,lower,_ = event_objective(compiled,reverse=reverse)
                encoded += 1
                if lower>minimum:continue
                compiled.solver.add(objective==minimum)
                for key,pattern in list(pending.items()):
                    remaining = int(1000*(deadline-time.perf_counter()))
                    if remaining<=0:exhausted=False;break
                    status,witness = query(compiled,canonical,pattern,generators,reverse=reverse,timeout_ms=min(2000,remaining))
                    queries += 1
                    if status=='represented':
                        results[key] = dict(status=status,method='compressed_path_query',direction=direction,
                            terminal=terminal,transitions=path.transitions,witness=witness)
                        del pending[key]
                    elif status=='unknown':unknown.add(key)
                if not pending or not exhausted:break
            if not pending or not exhausted:break
        if not exhausted:break
    for key in pending:
        results[key] = dict(status='excluded_from_saved_families' if exhausted and key not in unknown else 'unresolved',
                            method='compressed_path_query')
    return results,dict(queries=queries,encoded=encoded,paths_checked=paths_checked,
                       invariant_terminals=invariant_terminals,exhausted=exhausted)


def families(args):
    started = time.perf_counter()
    status_path = args.run/f'status/families_{args.index}.json'
    assert not status_path.exists()
    status = dict(index=args.index,started=time.time(),host=socket.gethostname())
    save(status_path,status)
    snapshot = read(args.run/f'snapshots/{args.index}.json')
    raw = read(AAM/f'inputs/{args.index}/input.json')
    canonical = EventPatterns(raw)
    generators = tuple(tuple(g[:canonical.n]) for g in pynauty.autgrp(colored_graph([canonical.features[0]]))[0])
    patterns = {k:v for method in snapshot['methods'].values() for k,v in method['patterns'].items()}
    result = dict(index=args.index,name=raw['name'],patterns=patterns,methods={},query_metrics={},
        scope='Membership of the union of observed minimum-event patterns in saved compressed outputs. '
              'Does not enumerate unmaterialized patterns absent from all three representative sets.')
    for method in ('native_slap','slap_sweep','aam'):
        deadline = time.perf_counter()+args.seconds if method=='aam' else time.perf_counter()+180
        if method=='aam':found,metrics = aam_membership(raw,canonical,patterns,snapshot,generators,deadline)
        else:found,metrics = native_membership(raw,canonical,patterns,method,snapshot,generators,deadline)
        result['methods'][method] = found
        result['query_metrics'][method] = metrics
        result['seconds'] = time.perf_counter()-started
        save(args.run/f'families/{args.index}.json',result)
    save(status_path,dict(status,finished=time.time(),exit=0))
    print(json.dumps(dict(index=args.index,seconds=result['seconds'],
        counts={m:{s:sum(v['status']==s for v in values.values()) for s in
            ('represented','excluded_from_saved_families','unresolved','different_score')}
            for m,values in result['methods'].items()})),flush=True)


def submit(args):
    assert not (args.run/'family_submission.json').exists()
    assert len(list((args.run/'snapshots').glob('*.json')))==140
    shutil.copy2(__file__,args.run/'family_driver.py')
    shutil.copy2(Path(__file__).with_name('holdout_minimum_events.py'),args.run/'holdout_minimum_events.py')
    command = ['sbatch','--parsable','--partition=cpunodes_nia','--exclude=bosque49,bosque50,bosque56',
        '--nodes=1','--cpus-per-task=1','--mem=8G','--time=00:25:00','--no-requeue',
        '--array=0-139%32','--job-name=min_event_families',f'--output={args.run}/status/families_%A_%a.out',
        '--wrap',shlex.join(['env',f'PYTHONPATH={AAM}/original/src:{AAM}/engine/bench',
            'PYTHONHASHSEED=0','PYTHONDONTWRITEBYTECODE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1',
            'timeout','--kill-after=5s','1300',PYTHON,str(args.run/'family_driver.py'),'families',
            '--run',str(args.run),'--seconds',str(args.seconds),'--index'])+' "$SLURM_ARRAY_TASK_ID"']
    job = subprocess.check_output(command,text=True).strip()
    save(args.run/'family_submission.json',dict(job=job,command=command,
        drivers_sha256={n:sha(args.run/n) for n in ('family_driver.py','holdout_minimum_events.py')}))
    print(job,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('families','submit'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int)
    parser.add_argument('--seconds',type=int,default=600)
    args = parser.parse_args()
    globals()[args.command](args)

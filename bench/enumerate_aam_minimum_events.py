"""Bounded event-pattern enumeration across every saved full AAM path family.

Invariant subgraphs are certified without expansion. Noninvariant families are
queried at the recorded benchmark minimum, blocking complete literal event sets.
Incomplete scans remain lower bounds; they never certify missing patterns.
"""
import argparse
from functools import lru_cache
from pathlib import Path
import os
import shlex
import shutil
import socket
import subprocess
import time
from types import SimpleNamespace

import z3
from holdout_minimum_event_families import *
from enumerate_slap_minimum_events import pattern_constraint


def case(args):
    started=time.perf_counter();deadline=started+args.seconds
    status_path=args.run/f'status/aam_enumeration_{args.index}.json'
    assert not status_path.exists()
    status=dict(index=args.index,started=time.time(),host=socket.gethostname())
    save(status_path,status)
    raw=read(AAM/f'inputs/{args.index}/input.json');canonical=EventPatterns(raw)
    problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    snapshot=read(args.run/f'snapshots/{args.index}.json')
    enriched=args.run/f'families_enriched/{args.index}.json'
    membership=read(enriched if enriched.exists() else args.run/f'families/{args.index}.json')
    minimum=snapshot['methods']['aam']['minimum']
    patterns={k:v['witness'] for k,v in membership['methods']['aam'].items() if v['status']=='represented'}
    metrics=dict(encoded_paths=0,scanned_paths=0,invariant_terminals=0,invariant_paths=0,
                 additional_models=0,unknown_families=0,completed_families=0,duplicate_families=0)
    complete=True;reasons=[]
    directions=getattr(args,'directions',('R_to_P','P_to_R'))
    def checkpoint():
        save(args.run/f'aam_enumeration/{args.index}.json',dict(index=args.index,minimum=minimum,patterns=patterns,
            complete=complete,metrics=metrics,reasons=reasons,seconds=time.perf_counter()-started,directions=directions,
            scope='All saved path families at the recorded benchmark minimum if complete; otherwise verified event-pattern lower bound.'))
    for direction in directions:
        if time.perf_counter()>deadline:complete=False;reasons.append('case budget');break
        reverse=direction=='P_to_R';search_problem=AAMProblem(problem.product,problem.reactant) if reverse else problem
        feature=canonical.features[0 if reverse else 1]
        graph=read_aam_checkpoint(AAM/f'results/holdout/{args.index}/{direction}/original/cuts/aam.pkl.gz').graph
        @lru_cache(None)
        def invariant_action(images):return exact_action(images,feature)
        @lru_cache(None)
        def invariant_edge(edge):
            placement=graph.fragment_placement(edge)
            if placement is None:return True
            assert placement.target_generators is not None
            if any(not invariant_action(tuple(g.images)) for g in placement.target_generators):return False
            for domain in placement.symmetry_domains:
                if domain.source=='exact_automorph_group':continue
                for a,b in zip(domain.p_atoms,domain.p_atoms[1:]):
                    images=list(range(canonical.n));images[a],images[b]=b,a
                    if not invariant_action(tuple(images)):return False
            return True
        incoming=[[] for _ in graph.states]
        for edge in graph.transitions:incoming[edge.target].append(edge)
        @lru_cache(None)
        def can_change(state):return any(can_change(e.source) or not invariant_edge(e.id) for e in incoming[state])
        seen=set()
        for terminal in graph.terminals:
            if len(graph.states[terminal].mapping)!=canonical.n:continue
            if time.perf_counter()>deadline:complete=False;reasons.append('case budget');break
            if not can_change(terminal):metrics['invariant_terminals']+=1;continue
            for path in graph.paths(terminal):
                if time.perf_counter()>deadline:complete=False;reasons.append('case budget');break
                metrics['scanned_paths']+=1
                if all(invariant_edge(e) for e in path.transitions):metrics['invariant_paths']+=1;continue
                family_key=(tuple(sorted(path.mapping.items())),tuple(path.context.cuts),
                    path.context.graph_floor,path.context.iso_tolerance,frozen_value(path.fragments))
                if family_key in seen:metrics['duplicate_families']+=1;continue
                seen.add(family_key)
                compiled=compile_path(path,search_problem,{},source_atoms=(),complete_reference=False)
                objective,lower,_=event_objective(compiled,reverse=reverse)
                metrics['encoded_paths']+=1
                if lower>minimum:metrics['completed_families']+=1;continue
                compiled.solver.add(objective==minimum)
                view=SimpleNamespace(values=original_values(compiled,reverse))
                physical=dict(path.mapping)
                if reverse:physical={b:a for a,b in physical.items()}
                seed=canonical.describe([physical[a] for a in range(canonical.n)])
                if seed['total']==minimum:
                    patterns.setdefault(seed['id'],dict(seed,direction=direction,terminal=terminal))
                    compiled.solver.add(z3.Not(pattern_constraint(view,canonical,seed)))
                family_complete=False;models=0
                while time.perf_counter()<deadline:
                    compiled.solver.set(timeout=max(1,min(3000,int(1000*(deadline-time.perf_counter())))))
                    result=compiled.solver.check()
                    if result==z3.unsat:family_complete=True;break
                    if result==z3.unknown:
                        reasons.append(dict(direction=direction,terminal=terminal,reason=compiled.solver.reason_unknown()))
                        break
                    physical=dict(compiled.realize(compiled.solver.model())['mapping'])
                    if reverse:physical={b:a for a,b in physical.items()}
                    pattern=canonical.describe([physical[a] for a in range(canonical.n)])
                    assert pattern['total']==minimum
                    patterns.setdefault(pattern['id'],dict(pattern,direction=direction,terminal=terminal,
                        transitions=path.transitions,source='compressed_family_enumeration'))
                    compiled.solver.add(z3.Not(pattern_constraint(view,canonical,pattern)))
                    metrics['additional_models']+=1;models+=1
                    if models>=10000:
                        reasons.append(dict(direction=direction,terminal=terminal,reason='literal-pattern limit'))
                        break
                if family_complete:metrics['completed_families']+=1
                else:complete=False;metrics['unknown_families']+=1
                if metrics['encoded_paths']%50==0:checkpoint()
            if time.perf_counter()>deadline:break
        checkpoint()
    if time.perf_counter()>deadline:complete=False
    checkpoint()
    save(status_path,dict(status,finished=time.time(),exit=0))
    print(dict(index=args.index,patterns=len(patterns),complete=complete,**metrics),flush=True)


def submit(args):
    assert not (args.run/'aam_enumeration_submission.json').exists()
    shutil.copy2(__file__,args.run/'enumerate_aam.py')
    shutil.copy2(Path(__file__).with_name('enumerate_slap_minimum_events.py'),args.run/'enumerate_slap_minimum_events.py')
    command=['sbatch','--parsable','--partition=cpunodes_nia','--exclude=bosque49,bosque50,bosque56,bosque59,bosque70',
        '--nodes=1','--cpus-per-task=1','--mem=8G','--time=00:10:00','--no-requeue',
        '--array=0-139%32','--job-name=aam_min_events',f'--output={args.run}/status/aam_enumeration_%A_%a.out',
        '--wrap',shlex.join(['env',f'PYTHONPATH={AAM}/original/src:{AAM}/engine/bench','PYTHONHASHSEED=0',
            'PYTHONDONTWRITEBYTECODE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1',
            'timeout','--kill-after=5s','500',PYTHON,str(args.run/'enumerate_aam.py'),'case',
            '--run',str(args.run),'--seconds',str(args.seconds),'--index'])+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(command,text=True).strip()
    save(args.run/'aam_enumeration_submission.json',dict(job=job,command=command,
        drivers_sha256={n:sha(args.run/n) for n in ('enumerate_aam.py','enumerate_slap_minimum_events.py')}))
    print(job,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('case','submit'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int)
    parser.add_argument('--seconds',type=int,default=240)
    args=parser.parse_args();globals()[args.command](args)

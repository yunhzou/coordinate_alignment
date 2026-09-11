"""Enumerate minimum-score event patterns in native SLAP label families.

Block event patterns, not atom permutations. Safe twin-atom symmetry breaking
preserves every event class. An unfinished family is explicitly unresolved.
"""
import argparse
from collections import defaultdict
from functools import lru_cache
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import time

import z3
from holdout_minimum_event_families import (slap_model,condition_for_event,original_values,exact_action,
    event_objective,EventPatterns,read,save,sha,AAM,SWEEP,PYTHON)


def family_key(candidate):
    groups = [defaultdict(list),defaultdict(list)]
    for side,g in enumerate(candidate['graphs']):
        for atom,label in enumerate(g['labels']):groups[side][label].append(atom)
    return tuple(sorted((tuple(left),tuple(groups[1][label])) for label,left in groups[0].items()))


def pattern_constraint(compiled,canonical,pattern):
    return z3.And(*(condition_for_event(kind,compiled.values[a],compiled.values[b],canonical.r[a,b],canonical.p)
        for kind,edges in pattern['events'].items() for a,b in edges))


def add_twin_constraints(compiled,candidate,canonical,twin):
    groups = [defaultdict(list),defaultdict(list)]
    for side,graph in enumerate(candidate['graphs']):
        for atom,label in enumerate(graph['labels']):groups[side][label].append(atom)
    inverse = original_values(compiled,True)
    restrictions = 0
    for side in range(2):
        values = compiled.values if side==0 else inverse
        for group in groups[side].values():
            clusters = []
            for atom in group:
                for cluster in clusters:
                    if twin(side,cluster[0],atom):cluster.append(atom);break
                else:clusters.append([atom])
            for cluster in clusters:
                for a,b in zip(cluster,cluster[1:]):
                    compiled.solver.add(values[a][0]<values[b][0])
                    restrictions += 1
    return restrictions


def enumerate_family(raw,canonical,candidate,score,twin,deadline,check_ms=3000):
    compiled = slap_model(raw,candidate)
    minimum = score['events']['total']
    objective,_,_ = event_objective(compiled)
    compiled.solver.add(objective==minimum)
    restrictions = add_twin_constraints(compiled,candidate,canonical,twin)
    mapping = dict(score['mapping'])
    first = canonical.describe([mapping[a] for a in range(canonical.n)])
    patterns = {first['id']:first}
    compiled.solver.add(z3.Not(pattern_constraint(compiled,canonical,first)))
    models = 0
    complete = False
    reason = None
    while time.perf_counter()<deadline:
        compiled.solver.set(timeout=max(1,min(check_ms,int(1000*(deadline-time.perf_counter())))))
        status = compiled.solver.check()
        if status==z3.unsat:complete=True;break
        if status==z3.unknown:reason=compiled.solver.reason_unknown();break
        model = compiled.solver.model()
        vector = [v[0] if isinstance(v[0],int) else model.eval(v[0]).as_long() for v in compiled.values]
        pattern = canonical.describe(vector)
        assert pattern['total']==minimum
        for a,b in enumerate(vector):
            assert candidate['graphs'][0]['labels'][a]==candidate['graphs'][1]['labels'][b]
        patterns.setdefault(pattern['id'],pattern)
        compiled.solver.add(z3.Not(pattern_constraint(compiled,canonical,pattern)))
        models += 1
        if models>=10000:reason='literal-pattern limit';break
    if not complete and reason is None:reason='case time budget'
    return dict(patterns=patterns,complete=complete,reason=reason,additional_models=models,twin_constraints=restrictions)


def case(args):
    started = time.perf_counter()
    status_path = args.run/f'status/enumeration_{args.index}.json'
    assert not status_path.exists()
    status = dict(index=args.index,started=time.time(),host=socket.gethostname())
    save(status_path,status)
    raw = read(AAM/f'inputs/{args.index}/input.json')
    canonical = EventPatterns(raw)
    @lru_cache(None)
    def twin(side,a,b):
        images = list(range(canonical.n))
        images[a],images[b] = b,a
        return exact_action(images,canonical.features[side])
    cache,methods = {},{}
    deadline = started+args.seconds
    for method,root in [('native_slap',AAM),('slap_sweep',SWEEP)]:
        native = read(root/f'slap_xyz/{args.index}.json')['candidates']
        scored = read(root/f'slap_scoring/{args.index}.refined.json')
        minimum = scored['slap_min_representative_events']
        patterns,families = {},[]
        for ordinal,(candidate,score) in enumerate(zip(native,scored['slap'],strict=True)):
            assert score['hydrogen_score_optimization']['optimal']
            if score['events']['total']!=minimum:continue
            key = family_key(candidate),minimum
            if key not in cache:
                cache[key] = enumerate_family(raw,canonical,candidate,score,twin,deadline)
            result = cache[key]
            for k,v in result['patterns'].items():patterns.setdefault(k,dict(v,candidate=ordinal))
            families.append(dict(candidate=ordinal,**{k:v for k,v in result.items() if k!='patterns'},
                                 pattern_ids=sorted(result['patterns'])))
        methods[method] = dict(minimum=minimum,patterns=patterns,families=families,
            complete=all(f['complete'] for f in families))
        save(args.run/f'enumeration/{args.index}.json',dict(index=args.index,methods=methods,seconds=time.perf_counter()-started,
            scope='All event classes in saved SLAP minimum-score label families when complete; otherwise a certified lower bound.'))
    save(status_path,dict(status,finished=time.time(),exit=0))
    print({m:dict(patterns=len(v['patterns']),complete=v['complete']) for m,v in methods.items()},flush=True)


def submit(args):
    assert not (args.run/'enumeration_submission.json').exists()
    shutil.copy2(__file__,args.run/'enumerate_slap.py')
    shutil.copy2(Path(__file__).with_name('holdout_minimum_event_families.py'),args.run/'holdout_minimum_event_families.py')
    command = ['sbatch','--parsable','--partition=cpunodes_nia','--exclude=bosque49,bosque50,bosque56',
        '--nodes=1','--cpus-per-task=1','--mem=8G','--time=00:10:00','--no-requeue',
        '--array=0-139%32','--job-name=min_event_enumerate',f'--output={args.run}/status/enumeration_%A_%a.out',
        '--wrap',shlex.join(['env',f'PYTHONPATH={AAM}/original/src:{AAM}/engine/bench','PYTHONHASHSEED=0',
            'PYTHONDONTWRITEBYTECODE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1',
            'timeout','--kill-after=5s','500',PYTHON,str(args.run/'enumerate_slap.py'),'case','--run',str(args.run),
            '--seconds',str(args.seconds),'--index'])+' "$SLURM_ARRAY_TASK_ID"']
    job = subprocess.check_output(command,text=True).strip()
    save(args.run/'enumeration_submission.json',dict(job=job,command=command,
        drivers_sha256={n:sha(args.run/n) for n in ('enumerate_slap.py','holdout_minimum_event_families.py')}))
    print(job,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('case','submit'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int)
    parser.add_argument('--seconds',type=int,default=240)
    args = parser.parse_args()
    globals()[args.command](args)

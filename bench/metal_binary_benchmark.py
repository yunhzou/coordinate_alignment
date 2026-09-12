"""Frozen-engine binary-metal ablation; raw WBO remains the event authority."""
import argparse
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import time
from types import SimpleNamespace

import numpy as np
from metal_binary_events import binary_metal_input, DeltaPatterns, scalar_events, membership, recanonicalize_patterns

DATA=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
FROZEN=DATA/'holdout_cap1000_seed1_20260910'
GOLDEN=DATA/'aam_one_seed_bidirectional_20260910'
PYTHON='/project/yunhengzou/coordinate_alignment/.venv/bin/python'
VARIANTS=('original','metal_binary')


def read(p):return json.loads(Path(p).read_text())
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp')
    t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def variant_root(run,dataset,variant):return run/'runs'/dataset/variant
def folder(run,task,variant):return variant_root(run,task['dataset'],variant)/f"results/{task['dataset']}/{task['index']}/{task['direction']}/original"
def environment(run):
    return dict(os.environ,PYTHONPATH=f'{run}/original/src:{run}/engine/bench',PYTHONHASHSEED='0',
                PYTHONDONTWRITEBYTECODE='1',RXN_CORE_NATIVE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')


def prepare(args):
    args.run.mkdir(parents=True,exist_ok=False)
    source=read(FROZEN/'manifest.json');assert source['original_commit'].startswith('98b01b1')
    for name in ('original','engine'):
        shutil.copytree(FROZEN/name,args.run/name,ignore=shutil.ignore_patterns('__pycache__'))
    for relative,expected in source['frozen_sha256'].items():assert sha(args.run/relative)==expected,relative
    for name in ('metal_binary_benchmark.py','metal_binary_events.py'):
        shutil.copy2(Path(__file__).with_name(name),args.run/name)
    changed={};hashes={};tasks=[]
    for dataset,count,root in [('holdout',140,FROZEN),('golden',1851,GOLDEN)]:
        changed[dataset]=[];hashes[dataset]={}
        inputs={}
        for index in range(count):
            src=root/f'inputs/{dataset}/{index}/input.json';raw=read(src);binary=binary_metal_input(raw)
            is_changed=raw!=binary
            if is_changed:changed[dataset].append(index)
            hashes[dataset][str(index)]=dict(original_input_sha256=sha(src),search_input_changes=is_changed)
            if dataset=='golden' and not is_changed:continue
            inputs[index]=(raw,binary)
            dest=args.run/f'raw/{dataset}/{index}';dest.mkdir(parents=True)
            shutil.copy2(src,dest/'input.json')
            if dataset=='golden':shutil.copy2(src.with_name('reference.json'),dest/'reference.json')
        order=([101,11,64,104]+[i for i in inputs if i not in (101,11,64,104)]) if dataset=='holdout' else list(inputs)
        for index in order:
            for direction in (('R_to_P',) if dataset=='holdout' else ('R_to_P','P_to_R')):
                tasks.append(dict(dataset=dataset,index=index,direction=direction))
        for variant in VARIANTS:
            dest=variant_root(args.run,dataset,variant);dest.mkdir(parents=True)
            config=dict(source['original_config'],branch_limit=1000 if dataset=='holdout' else 100)
            assert config['iso_tolerance']==1.0 and config['seed_count']==1
            save(dest/'manifest.json',dict(original_config=config,original_workers=8,
                original_execution=source['original_execution'],original_commit=source['original_commit'],root_seed=42))
            for index,(raw,binary) in inputs.items():
                d=dest/f'inputs/{dataset}/{index}';d.mkdir(parents=True)
                if variant=='original':(d/'input.json').symlink_to(args.run/f'raw/{dataset}/{index}/input.json')
                else:save(d/'input.json',binary)
    assert changed['golden']==[691,692,904,937,959,1093,1312,1418,1447,1474,1475]
    for dataset in ('holdout','golden'):
        for variant in VARIANTS:save(variant_root(args.run,dataset,variant)/'tasks.json',tasks)
    for phase in ('search','analyze','query'):(args.run/f'status/{phase}').mkdir(parents=True)
    save(args.run/'tasks.json',tasks)
    shutil.copy2(DATA/'holdout_missing_pattern_seeds_20260910/targets.json',args.run/'historical_targets.json')
    save(args.run/'input_audit.json',dict(changed=changed,hashes=hashes))
    save(args.run/'manifest.json',dict(schema='binary_metal_ablation/v1',original_commit=source['original_commit'],
        frozen_source=str(FROZEN),golden_baseline=str(GOLDEN),tasks=len(tasks),variants=VARIANTS,
        holdout=dict(cases=140,directions=['R_to_P'],cap=1000),golden=dict(cases=changed['golden'],directions=['R_to_P','P_to_R'],cap=100),
        seeds=1,iso_tolerance=1.0,metal_connectivity_threshold=.2,graph_floor=.2,workers=8,
        transformation='Every metal-containing pair is 1 if original WBO >= 0.2, otherwise 0; all other pairs unchanged.',
        scoring='Original raw-WBO classify_bonds: decrease >= pair threshold is broken/weakened; increase >= threshold is formed/strengthened. Ordinary 0.5, metal-containing 0.3, inclusive. No separate event bond-floor gate.',
        historical_score_caveat='Older benchmark helper counted graph-floor changes plus order shifts >0.5, without a metal override. Both arms are freshly rescored here; old event totals are not directly comparable.',
        timing='Fresh paired searches on same CPU allocation; alternating order by slot. Parent plus child CPU, excluding measured checkpoint persistence/loading. Input preprocessing and offline scoring reported separately.',
        changed_inputs=changed,search_watchdog=300,analysis_watchdog=300,query_watchdog=240,
        original_config=source['original_config'],frozen_sha256=source['frozen_sha256'],
        drivers_sha256={name:sha(args.run/name) for name in ('metal_binary_benchmark.py','metal_binary_events.py')},
        targets='Three historical missing witnesses are reclassified using raw 0.5/0.3 deltas; targets never enter search.'))
    print(dict(run=str(args.run),tasks=len(tasks),changed_inputs=changed),flush=True)


def submit(args):
    assert not (args.run/'submissions.json').exists()
    tasks=read(args.run/'tasks.json');jobs=[]
    for phase in ('search','analyze'):
        cpus=8 if phase=='search' else 1
        command=['sbatch','--parsable','--partition=cpunodes_nia','--nodelist=bosque83,bosque84','--nodes=1',
            f'--cpus-per-task={cpus}','--mem=16G','--time=00:15:00','--no-requeue',
            f'--array=0-{len(tasks)-1}%{16 if phase=="search" else 32}',f'--job-name=metal_binary_{phase}',
            f'--output={args.run}/status/{phase}_%A_%a.out']
        if jobs:command.append('--dependency=afterany:'+jobs[0]['job'])
        command+=['--wrap',shlex.join(['env',f'PYTHONPATH={args.run}/original/src:{args.run}/engine/bench',
            'PYTHONDONTWRITEBYTECODE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1',PYTHON,
            str(args.run/'metal_binary_benchmark.py'),'worker','--run',str(args.run),'--phase',phase,'--slot'])+' "$SLURM_ARRAY_TASK_ID"']
        job=subprocess.check_output(command,text=True).strip().split(';')[0]
        jobs.append(dict(phase=phase,job=job,command=command));save(args.run/'submissions.json',jobs);print(phase,job,flush=True)


def worker(args):
    task=read(args.run/'tasks.json')[args.slot];status=args.run/f'status/{args.phase}/{args.slot}.json';assert not status.exists()
    row=dict(**task,slot=args.slot,host=socket.gethostname(),affinity=sorted(os.sched_getaffinity(0)),
             job=os.environ.get('SLURM_JOB_ID'),started=time.time(),variants={});save(status,row)
    for variant in VARIANTS[::(-1 if args.slot%2 else 1)]:
        if args.phase=='search':
            command=[PYTHON,str(args.run/'engine/bench/adaptive_full_benchmark.py'),'search','--run',str(variant_root(args.run,task['dataset'],variant)),
                     '--slot',str(args.slot),'--method','original']
        else:
            if not (folder(args.run,task,variant)/'search.json').exists():
                row['variants'][variant]=dict(status='search_incomplete');save(status,row);continue
            command=[PYTHON,str(args.run/'metal_binary_benchmark.py'),'analyze','--run',str(args.run),'--slot',str(args.slot),'--variant',variant]
        start=time.monotonic()
        with status.with_name(f'{args.slot}_{variant}.log').open('w') as log:
            code=subprocess.run(['timeout','--kill-after=5s','300',*command],env=environment(args.run),stdout=log,stderr=subprocess.STDOUT).returncode
        row['variants'][variant]=dict(exit=code,wall_including_startup_io=time.monotonic()-start);save(status,row)
    row['finished']=time.time();save(status,row)


def audit(aam):
    r,p=aam.problem.reactant,aam.problem.product;g=aam.graph;checks=0;generators=set()
    active=g.ancestor_transitions(g.terminals);unfinalized_dead=0
    for state in g.states:
        m=dict(state.mapping)
        assert len(m)==len(set(m.values()))
        assert all(r.elements[a]==p.elements[b] for a,b in m.items())
    for edge in g.transitions:
        if not edge.match:continue
        placement=g.fragment_placement(edge.id);m=dict(placement.representative_assignments)
        for a,b in edge.preserved_bonds:
            assert p.wbo[m[a],m[b]]>=aam.config.graph_floor
            assert abs(r.wbo[a,b]-p.wbo[m[a],m[b]])<=aam.config.iso_tolerance+1e-9
            checks+=1
        if placement.target_generators is None:
            assert edge.id not in active, 'A returned path lacks finalized symmetry'
            unfinalized_dead+=1
        else:
            generators.update(tuple(gen.images) for gen in placement.target_generators)
    adjacency=p.wbo>=aam.config.graph_floor
    for gen in generators:
        assert sorted(gen)==list(range(p.atom_count))
        assert all(p.elements[a]==p.elements[b] for a,b in enumerate(gen))
        assert np.array_equal(adjacency,adjacency[np.ix_(gen,gen)])
    return dict(states=len(g.states),preserved_representative_bonds=checks,unique_generators=len(generators),
        unfinalized_discarded_transitions=unfinalized_dead,violations=0,
        scope='All state elements/injectivity, stored transition representative preserved bonds, and full target adjacency/element invariance of stored generators. Does not enumerate every family member.')


def analyze(args):
    from rxn_core.artifacts import read_aam_checkpoint
    task=read(args.run/'tasks.json')[args.slot];f=folder(args.run,task,args.variant)
    raw=read(args.run/f"raw/{task['dataset']}/{task['index']}/input.json")
    start=time.monotonic();aam=read_aam_checkpoint(f/'cuts/aam.pkl.gz');loading=time.monotonic()-start
    start=time.monotonic();cpu=time.process_time()
    if task['dataset']=='holdout':
        canonical=DeltaPatterns(raw);vectors=sorted({tuple(dict(aam.graph.states[t].mapping)[a] for a in range(canonical.n))
            for t in aam.graph.terminals if len(aam.graph.states[t].mapping)==canonical.n})
        best=None;patterns={};histogram={}
        target=None
        if task['index'] in (11,64,101):
            old=read(args.run/'historical_targets.json')[str(task['index'])]['pattern']
            target=canonical.describe(old['mapping'])
        target_witnesses={}
        for offset in range(0,len(vectors),128):
            batch=vectors[offset:offset+128];scores=canonical.counts(batch).sum(axis=1)
            for vector,score in zip(batch,scores):
                score=int(score);histogram[score]=histogram.get(score,0)+1
                if best is None or score<best:best=score;patterns={}
                if score==best or target is not None and score==target['total']:
                    pattern=canonical.describe(vector)
                    if score==best:patterns.setdefault(pattern['id'],pattern)
                    if target is not None and pattern['id']==target['id']:target_witnesses[pattern['id']]=pattern
        result=dict(full_mapping=bool(vectors),full_representatives=len(vectors),best_representative_events=best,
            minimum_patterns=patterns,event_histogram=histogram,target=target,target_representatives=target_witnesses,
            scoring='raw_delta_0p5_ordinary_0p3_metal_inclusive')
    else:
        from golden_evaluation import evaluate,project,colored_graph
        import pynauty
        ref=read(args.run/f"raw/golden/{task['index']}/reference.json");reverse=task['direction']=='P_to_R'
        reference=dict(ref['mapping']);reference={v:k for k,v in reference.items()} if reverse else reference
        features=list(reversed(ref['features'])) if reverse else ref['features']
        def rank(mapping):
            m=dict(mapping);m={v:k for k,v in m.items()} if reverse else m
            events=scalar_events(raw,m)
            heavy=sum(raw['product']['elements'][p]!='H' for p in m.values())
            return (-heavy,-len(m),events['total'],tuple(sorted(m.items()))),dict(broken=len(events['broken']),formed=len(events['formed']))
        result=evaluate(aam,features,reference,seconds=60,query_timeout_ms=1500,
                        reference_side='source' if reverse else 'target',ranker=rank)
        if result['reference_recovery']=='recovered':
            m=dict(aam.graph.states[result['witness_terminal']].mapping) if result['representative_recovery'] else dict(result['witness_actions']['mapping'])
            m={v:k for k,v in m.items()} if reverse else m
            assert len(m)==len(set(m.values()))
            assert all(raw['reactant']['elements'][r]==raw['product']['elements'][p] for r,p in m.items())
            assert pynauty.certificate(colored_graph(ref['features'],project(m,ref['features'])))==pynauty.certificate(colored_graph(ref['features'],project(ref['mapping'],ref['features'])))
            result['input_witness']=sorted(m.items())
    result.update(analysis_cpu=time.process_time()-cpu,analysis_wall=time.monotonic()-start,loading_wall=loading)
    result['audit']=audit(aam)
    save(f/'raw_evaluation.json',result)
    print(dict(task=task,variant=args.variant,best=result.get('best_representative_events'),recovery=result.get('reference_recovery')),flush=True)


def query(args):
    from rxn_core.artifacts import read_aam_checkpoint
    task=read(args.run/'tasks.json')[args.slot];assert task['dataset']=='holdout'
    raw=read(args.run/f"raw/holdout/{task['index']}/input.json");canonical=DeltaPatterns(raw)
    analyses={v:read(folder(args.run,task,v)/'raw_evaluation.json') for v in VARIANTS}
    # Stored evaluations can predate the sparse, versioned ID schema.
    for value in analyses.values():
        for field in ('minimum_patterns', 'target_representatives'):
            value[field] = recanonicalize_patterns(canonical, value[field])
        if value['target'] is not None:
            value['target'] = canonical.describe(value['target']['mapping'])
    patterns={}
    for v,value in analyses.items():
        patterns.update(value['minimum_patterns'])
        if value['target'] is not None:patterns[value['target']['id']]=value['target']
    result=dict(**task,patterns=patterns,variants={},complete=False)
    output=args.run/f"family_queries/{task['index']}.json"
    save(output,result)
    for variant in VARIANTS:
        f=folder(args.run,task,variant);aam=read_aam_checkpoint(f/'cuts/aam.pkl.gz')
        # A counterpart's minimum may be a non-minimum representative here.
        # Scan every requested score before using score-invariance shortcuts.
        saved=dict(analyses[variant]['minimum_patterns'],**analyses[variant]['target_representatives'])
        wanted={p['total'] for p in patterns.values()}
        vectors=sorted({tuple(dict(aam.graph.states[t].mapping)[a] for a in range(canonical.n))
            for t in aam.graph.terminals if len(aam.graph.states[t].mapping)==canonical.n})
        for offset in range(0,len(vectors),128):
            batch=vectors[offset:offset+128]
            for vector,score in zip(batch,canonical.counts(batch).sum(axis=1)):
                if int(score) in wanted:
                    pattern=canonical.describe(vector)
                    if pattern['id'] in patterns:saved[pattern['id']]=pattern
        result['variants'][variant]=membership(aam,canonical,patterns,saved,args.seconds)
        save(output,result)
    result['complete']=True;save(output,result)
    print({v:{k:r['status'] for k,r in x['results'].items()} for v,x in result['variants'].items()},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','submit','worker','analyze','query'))
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--slot',type=int)
    parser.add_argument('--phase',choices=('search','analyze'));parser.add_argument('--variant',choices=VARIANTS)
    parser.add_argument('--seconds',type=float,default=120)
    args=parser.parse_args();globals()[args.command](args)

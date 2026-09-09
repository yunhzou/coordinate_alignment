"""Fresh cut replay ablations; exact graph equality, not just minimum score.

One Slurm task handles one reaction/direction/policy. Every cut graph is saved.
Matching, symmetry, scoring, persistence and audit encoding are timed separately.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
from types import SimpleNamespace

import numpy as np
from rxn_core.domain import AAMProblem, MolecularEndpoint
from rxn_core.aam import cut_seed
from rxn_core.alignment.branch import _generate_seed_orders, find_islands
from rxn_core.alignment.sweep import cut_sweep_items
from rxn_core.artifacts import write_graph_checkpoint
from rxn_core.cut_replay import CutReplay
from rxn_core.frag import build_graph
from rxn_core.matcher import _nauty_orbits
from rxn_core.search_graph import AAMSearchGraph
from rxn_core.search_symmetry import finalize_graph_symmetry, SymmetryWorkspace
from compare_elementary_outputs import event_counts

SOURCE=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_tol1_20260909')
POLICIES=('independent','shared_fresh','shared_whole','shared_checkpoint',
          'shared_all','independent_all','shared_rolling','independent_rolling')


def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2)+'\n')


def task(args):
    spec=json.loads((args.run/'tasks.json').read_text())[args.slot]
    index,direction,policy,seeds=(spec[k] for k in ('index','direction','policy','seeds'))
    folder=args.run/f"results/{index}/{direction}/{policy}"
    if seeds>1 and args.seed_index is None:
        folder.mkdir(parents=True,exist_ok=False)
        started=time.perf_counter()
        children=[SimpleNamespace(run=args.run,slot=args.slot,seed_index=i) for i in range(seeds)]
        with ProcessPoolExecutor(max_workers=seeds) as pool:list(pool.map(task,children))
        rows=[json.loads((folder/f'seed_{i:02d}/summary.json').read_text()) for i in range(seeds)]
        result=dict(rows[0]);result.update(spec)
        result['timings']={k:{axis:sum(r['timings'][k][axis] for r in rows) for axis in ('cpu','wall')}
                           for k in rows[0]['timings']}
        result['cuts']=[]
        for cuts in zip(*(r['cuts'] for r in rows),strict=True):
            result['cuts'].append(dict(cut=cuts[0]['cut'],
                digest=hashlib.sha256(''.join(c['digest'] for c in cuts).encode()).hexdigest(),
                best=min((c['best'] for c in cuts if c['best'] is not None),default=None),
                terminals=sum(c['terminals'] for c in cuts),states=sum(c['states'] for c in cuts),
                capped=any(c['capped'] for c in cuts)))
        result['best']=min((r['best'] for r in rows if r['best'] is not None),default=None)
        result['peak_rss_kib']=max(r['peak_rss_kib'] for r in rows)
        result['peak_rss_scope']='maximum single-worker high water mark, not whole job RSS'
        result['witnesses']=sum(r['witnesses'] for r in rows)
        result['witnesses_scope']='sum across seeds; not interseed deduplicated'
        result['replay']=({k:sum(r['replay'][k] for r in rows) for k in rows[0]['replay']}
                          if rows[0]['replay'] is not None else None)
        result['workers']=seeds
        result['job_elapsed_including_io']=time.perf_counter()-started
        result['timing_scope']='CPU phases sum across seed workers; phase wall is summed worker duration, not job latency. Every worker has a private replay cache. Imports, input load and queue excluded.'
        save(folder/'summary.json',result)
        return
    if args.seed_index is not None:folder=folder/f'seed_{args.seed_index:02d}'
    folder.mkdir(parents=True,exist_ok=False)
    raw=json.loads((SOURCE/f'inputs/{index}/input.json').read_text())
    original=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    problem=AAMProblem(original.product,original.reactant) if direction=='P_to_R' else original
    timings={}
    def measure(label,fn):
        cpu=time.process_time();wall=time.perf_counter()
        value=fn()
        row=timings.setdefault(label,dict(cpu=0.,wall=0.))
        row['cpu']+=time.process_time()-cpu;row['wall']+=time.perf_counter()-wall
        return value
    def prepare():
        r=build_graph(problem.reactant.elements,problem.reactant.wbo,bond_cut=.2)
        p=build_graph(problem.product.elements,problem.product.wbo,bond_cut=.2)
        po=_nauty_orbits(p,wbo_tol=1.)
        replay=(CutReplay(r,p,po,checkpoints=policy!='shared_whole',reference_only=not policy.endswith('rolling'))
                if policy not in ('independent','shared_fresh') else None)
        return r,p,po,replay,_generate_seed_orders(r,seeds,rng_seed=42)
    source,target,p_orbits,replay,shared_orders=measure('setup',prepare)
    workspace=SymmetryWorkspace(target,1.) if policy.endswith(('all','rolling')) else None
    cuts=cut_sweep_items(problem.reactant.wbo,.2)
    records=[];witnesses=set()
    for ordinal,cut in enumerate(cuts):
        def search():
            view=replay.for_cut(cut) if replay else None
            r=view.source if view else source.copy()
            if view is None:r.remove_edges_from(cut)
            ro=_nauty_orbits(r,wbo_tol=1.)
            orders=(_generate_seed_orders(r,seeds,rng_seed=cut_seed(cut))
                    if policy.startswith('independent') else shared_orders)
            if args.seed_index is not None:orders=[orders[args.seed_index]]
            graphs=[find_islands(r,target,order,graph_floor=.2,iso_tol=1.,max_branches=100,
                        r_orbits=ro,p_orbits=p_orbits,cuts=cut,growth_replay=view) for order in orders]
            return AAMSearchGraph.combine(graphs)
        graph=measure('search',search)
        graph,groups=measure('symmetry',lambda:finalize_graph_symmetry(graph,target,iso_tolerance=1.,workspace=workspace))
        def score():
            vectors=[]
            for terminal in graph.terminals:
                m=dict(graph.states[terminal].mapping)
                if len(m)!=original.atom_count:continue
                if direction=='P_to_R':m={b:a for a,b in m.items()}
                vectors.append(tuple(m[i] for i in range(original.atom_count)))
            witnesses.update(vectors)
            counts=event_counts(original.reactant.wbo,original.product.wbo,vectors) if vectors else []
            return min((int(sum(row)) for row in counts),default=None)
        best=measure('score',score)
        digest=measure('audit_encoding',lambda:hashlib.sha256(
            json.dumps(graph.to_record(copy=False),sort_keys=True,separators=(',',':')).encode()).hexdigest())
        measure('persistence',lambda:write_graph_checkpoint(graph,folder/f'cut_{ordinal:04d}.pkl.gz'))
        records.append(dict(cut=cut,digest=digest,best=best,terminals=len(graph.terminals),
                            capped=graph.capped,states=len(graph.states),groups=groups))
        save(folder/'progress.json',dict(completed=len(records),total=len(cuts),timings=timings))
    counts=event_counts(original.reactant.wbo,original.product.wbo,sorted(witnesses)) if witnesses else []
    save(folder/'witnesses.json',dict(mappings=sorted(witnesses),events=[list(map(int,row)) for row in counts]))
    save(folder/'summary.json',dict(**spec,name=raw['name'],cuts=records,timings=timings,
        best=min((r['best'] for r in records if r['best'] is not None),default=None),
        replay=None if replay is None else replay.stats(),witnesses=len(witnesses),
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        host=os.uname().nodename,pid=os.getpid(),seed_index=args.seed_index,settings=dict(branch_cap=100,match_tolerance=1.,
        event_tolerance=.5,bond_floor=.2,explicit_hydrogen=True,cache_bytes=64*1024*1024,
        reference_only=not policy.endswith('rolling')),
        timing_scope='Fresh process per reaction/direction/policy. Queue/import/input loading excluded. Persistence and audit encoding separated; no subtraction from summed worker time.'))


def prepare(args):
    tasks=[dict(index=i,direction=d,policy=p,seeds=args.seeds)
           for i in args.indices for d in ('R_to_P','P_to_R') for p in args.policies]
    save(args.run/'tasks.json',tasks)
    (args.run/'status').mkdir(exist_ok=True)
    root=Path(__file__).resolve().parents[1]
    shutil.copytree(root/'src',args.run/'engine/src',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(root/'bench',args.run/'engine/bench',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(root/'native',args.run/'engine/native')
    save(args.run/'manifest.json',dict(parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        source=str(SOURCE),tasks=len(tasks),seeds=args.seeds,policies=args.policies,
        native_sha256=hashlib.sha256(next((args.run/'engine/src/rxn_core').glob('_engine*.so')).read_bytes()).hexdigest()))
    print(len(tasks),'tasks prepared')


def submit(args):
    tasks=json.loads((args.run/'tasks.json').read_text())
    workers=max(t['seeds'] for t in tasks)
    cmd=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',
         'PYTHONHASHSEED=0','RXN_CORE_NATIVE=1',
         f'PYTHONPATH={args.run}/engine/src:{args.run}/engine/bench',
         'timeout','--kill-after=5s','300',sys.executable,str(args.run/'engine/bench/cut_replay_pilot.py'),
         'task','--run',str(args.run),'--slot']
    options=['sbatch','--parsable','--partition=cpunodes,cpunodes_nia',
             '--exclude=bosque5,bosque6,bosque8,bosque10',f'--cpus-per-task={workers}',f'--mem={max(4,workers*2)}G',
             '--time=00:10:00',f'--array=0-{len(tasks)-1}%{max(1,192//workers)}','--job-name=cut_replay',
             f'--output={args.run}/status/%A_%a.out','--wrap',shlex.join(cmd)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=options));print(job)


def report(args):
    tasks=json.loads((args.run/'tasks.json').read_text());rows=[];missing=[]
    scheduled = (set(slot for pair in json.loads((args.run/'pairs.json').read_text()) for slot in pair['slots'])
                 if (args.run/'pairs.json').exists() else set(range(len(tasks))))
    for slot,t in enumerate(tasks):
        if slot not in scheduled:continue
        path=args.run/f"results/{t['index']}/{t['direction']}/{t['policy']}/summary.json"
        if path.exists():rows.append(json.loads(path.read_text()))
        else:missing.append(slot)
    totals={}
    for p in POLICIES:
        chosen=[r for r in rows if r['policy']==p]
        if not chosen:continue
        totals[p]=dict(tasks=len(chosen),timings={k:sum(r['timings'][k]['cpu'] for r in chosen)
            for k in chosen[0]['timings']},best=[(r['index'],r['direction'],r['best']) for r in chosen],
            peak_rss_mib=max(r['peak_rss_kib'] for r in chosen)/1024,
            replay={k:sum(r['replay'][k] for r in chosen) for k in
                    ('calls','full_hits','prefix_hits','logical_extensions','reused_extensions',
                     'logical_certificates','reused_certificates','evictions')} if chosen[0]['replay'] else None)
    reference_rows=[]
    if args.reference is not None:
        reference_rows=[json.loads(f.read_text()) for f in args.reference.glob('results/*/*/*/summary.json')]
    base={(r['index'],r['direction'],r['policy']):r for r in reference_rows+rows}
    mismatches=[];checked=0
    for row in rows:
        if row['policy'] in ('shared_fresh','independent'):continue
        reference=base.get((row['index'],row['direction'],
                           'independent' if row['policy'].startswith('independent') else 'shared_fresh'))
        if reference is None:continue
        for a,b in zip(reference['cuts'],row['cuts'],strict=True):
            checked+=1
            if a['digest']!=b['digest']:mismatches.append((row['index'],row['direction'],row['policy'],a['cut']))
    comparisons=[]
    for policy,baseline in (('shared_all','shared_fresh'),('independent_all','independent'),
                            ('shared_rolling','shared_fresh'),('independent_rolling','independent')):
        pairs=[(base.get((r['index'],r['direction'],baseline)),r) for r in rows if r['policy']==policy]
        pairs=[(a,b) for a,b in pairs if a is not None]
        if not pairs:continue
        phases={p:{'before':sum(a['timings'][p]['cpu'] for a,b in pairs),
                   'after':sum(b['timings'][p]['cpu'] for a,b in pairs)}
                for p in ('setup','search','symmetry','score')}
        before=sum(v['before'] for v in phases.values());after=sum(v['after'] for v in phases.values())
        comparisons.append(dict(policy=policy,baseline=baseline,pairs=len(pairs),phases=phases,
            compute_cpu_before=before,compute_cpu_after=after,speedup=before/after,
            changed_best=[(a['index'],a['direction'],a['best'],b['best']) for a,b in pairs if a['best']!=b['best']]))
    summary=dict(completed=len(rows),scheduled=len(scheduled),missing=missing,totals=totals,comparisons=comparisons,
                 reference=None if args.reference is None else str(args.reference),
                 exact_cut_graph_comparisons=checked,mismatches=mismatches)
    save(args.run/'summary.json',summary);print(json.dumps(summary,indent=2))


def relocate(args):
    submission=json.loads((args.run/'submission.json').read_text())
    jobs=subprocess.check_output(['squeue','-h','-j',submission['job'],'-t','CONFIGURING','-o','%i'],text=True).split()
    tasks=json.loads((args.run/'tasks.json').read_text());slots=[];hosts=set()
    for job in jobs:
        slot=int(job.split('_')[1]);t=tasks[slot]
        output=args.run/f"results/{t['index']}/{t['direction']}/{t['policy']}"
        if output.exists():continue
        hosts.update(subprocess.check_output(['squeue','-h','-j',job,'-o','%N'],text=True).split())
        subprocess.run(['scancel',job],check=True);slots.append(slot)
    if not slots:print('No unstarted jobs to relocate');return
    options=[('--array='+','.join(map(str,slots))) if x.startswith('--array=') else
             x+','+','.join(sorted(hosts)) if x.startswith('--exclude=') else x
             for x in submission['command']]
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/f'relocation_{job}.json',dict(job=job,command=options,unstarted_slots=slots))
    print(job,slots)


def submit_paired(args):
    tasks=json.loads((args.run/'tasks.json').read_text())
    stress=[(25,'R_to_P'),(76,'P_to_R'),(77,'P_to_R'),(114,'R_to_P')]
    pairs=[]
    for i,d in stress:
        for policies in (('independent','independent_all','independent_rolling'),
                         ('shared_fresh','shared_all','shared_rolling')):
            slots=[next(j for j,t in enumerate(tasks) if (t['index'],t['direction'],t['policy'])==(i,d,p)) for p in policies]
            # Counterbalance which method runs first, without changing search.
            shift=len(pairs)%3;slots=slots[shift:]+slots[:shift]
            pairs.append(dict(index=i,direction=d,slots=slots))
    save(args.run/'pairs.json',pairs)
    cmd=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','PYTHONHASHSEED=0','RXN_CORE_NATIVE=1',
         f'PYTHONPATH={args.run}/engine/src:{args.run}/engine/bench',sys.executable,
         str(args.run/'engine/bench/cut_replay_pilot.py'),'paired','--run',str(args.run),'--slot']
    options=['sbatch','--parsable','--partition=cpunodes','--exclude=bosque5,bosque6,bosque8,bosque10',
             '--nodes=1','--cpus-per-task=10','--mem=24G','--time=00:10:00',
             f'--array=0-{len(pairs)-1}%8','--job-name=replay_paired',
             f'--output={args.run}/status/%A_%a.out','--wrap',shlex.join(cmd)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=options));print(job)


def paired(args):
    spec=json.loads((args.run/'pairs.json').read_text())[args.slot]
    started=time.perf_counter()
    for slot in spec['slots']:
        print('starting task',slot,flush=True)
        subprocess.run(['timeout','--kill-after=5s','300',sys.executable,__file__,
                        'task','--run',str(args.run),'--slot',str(slot)],check=True)
    save(args.run/f'pair_{args.slot}.json',dict(**spec,elapsed=time.perf_counter()-started,
        host=os.uname().nodename,cpu_info=Path('/proc/cpuinfo').read_text()))


def mechanisms(args):
    """Check saved witness event patterns; never use references during search."""
    input_path=SOURCE/f'inputs/{args.index}/input.json'
    raw=json.loads(input_path.read_text())
    r,p=(np.asarray(raw[k]['wbo']) for k in ('reactant','product'))
    a,b=np.triu_indices(len(r),1)
    def signature(mapping):
        mapped=p[np.ix_(mapping,mapping)]
        masks=((r[a,b]>.2)&(mapped[a,b]<=.2),
               (r[a,b]<=.2)&(mapped[a,b]>.2),
               (r[a,b]>.2)&(mapped[a,b]>.2)&(abs(r[a,b]-mapped[a,b])>.5))
        return tuple(tuple(zip(a[mask].tolist(),b[mask].tolist())) for mask in masks)
    references=[]
    for record in json.loads(args.mechanism_reference.read_text())['mechanisms']:
        mapping=[record['mapping_RP'][str(i)] for i in range(len(r))]
        references.append(dict(id=record['id'],signature=signature(mapping)))
    rows=[]
    for folder in sorted((args.run/f'results/{args.index}').glob('*/*')):
        matches={};count=0
        paths=sorted(folder.glob('seed_*/witnesses.json')) or sorted(folder.glob('witnesses.json'))
        for path in paths:
            for mapping in json.loads(path.read_text())['mappings']:
                count+=1
                matches.setdefault(signature(mapping),dict(mapping=mapping,artifact=str(path)))
        rows.append(dict(direction=folder.parent.name,policy=folder.name,witnesses=count,
            event_patterns=len(matches),reference_matches=[dict(id=ref['id'],
                recovered=ref['signature'] in matches,witness=matches.get(ref['signature'])) for ref in references]))
    result=dict(index=args.index,input=str(input_path),reference=str(args.mechanism_reference),
        input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
        reference_sha256=hashlib.sha256(args.mechanism_reference.read_bytes()).hexdigest(),
        bond_floor=.2,event_tolerance=.5,references=references,rows=rows,
        scope='Exact reactant-index event signatures of saved terminal witnesses, not chemical ground truth or an exhaustive expansion of compressed families.')
    save(args.run/f'mechanism_check_{args.index}.json',result)
    for row in rows:print(row['direction'],row['policy'],[m['recovered'] for m in row['reference_matches']])


def archive(args):
    """Export small evidence and frozen sources, leaving full graphs in place."""
    selected=set();indices=set()
    for run in args.runs:
        indices.update(t['index'] for t in json.loads((run/'tasks.json').read_text()))
        for pattern in ('*.json','*.log','*.txt','status/*','results/*/*/*/summary.json',
                        'results/*/*/*/seed_*/summary.json','results/114/*/*/seed_*/witnesses.json',
                        'engine/bench/cut_replay*.py','engine/native/src/*.h','engine/native/src/*.cpp',
                        'engine/src/rxn_core/cut_replay.py','engine/src/rxn_core/search_symmetry.py'):
            selected.update((path,run.name+'/'+str(path.relative_to(run))) for path in run.glob(pattern) if path.is_file())
    for index in indices:
        path=SOURCE/f'inputs/{index}/input.json'
        selected.add((path,f'inputs/{index}/input.json'))
    for run in args.runs:
        for path in run.glob('mechanism_check_*.json'):
            reference=Path(json.loads(path.read_text())['reference'])
            selected.add((reference,'references/'+reference.parent.name+'/'+reference.name))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with tarfile.open(args.output,'x:gz') as bundle:
        for path,name in sorted(selected,key=lambda item:item[1]):bundle.add(path,arcname=name)
    print(args.output,len(selected),'files',args.output.stat().st_size,'bytes')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','submit','task','report','relocate','submit_paired','paired','mechanisms','archive'))
    p.add_argument('--run',type=Path);p.add_argument('--slot',type=int)
    p.add_argument('--runs',type=Path,nargs='+');p.add_argument('--output',type=Path)
    p.add_argument('--seed-index',type=int)
    p.add_argument('--reference',type=Path)
    p.add_argument('--index',type=int);p.add_argument('--mechanism-reference',type=Path)
    p.add_argument('--indices',type=int,nargs='+',default=[0,4,59,64,114,123,135,136])
    p.add_argument('--seeds',type=int,default=1)
    p.add_argument('--policies',nargs='+',choices=POLICIES,default=list(POLICIES))
    a=p.parse_args();globals()[a.command](a)

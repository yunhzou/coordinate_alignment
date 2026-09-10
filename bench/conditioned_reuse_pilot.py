"""Paired exact-output / compute-CPU ablations, with frozen engines and outputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np

SOURCE = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_tol1_20260909')
OLD_ENGINE = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/real_ts_compare_20260909_EwGzgT/run/engine/src')
CASES = ((25, 'R_to_P'), (76, 'P_to_R'), (77, 'P_to_R'), (114, 'R_to_P'))
MODES = ('baseline', 'current', 'repair', 'shared_symmetry', 'symmetry', 'extensions', 'combined')


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2)+'\n')


def prepare(args):
    root = Path(__file__).resolve().parents[1]
    args.run.mkdir(parents=True, exist_ok=True)
    shutil.copytree(OLD_ENGINE, args.run/'baseline/src', ignore=shutil.ignore_patterns('__pycache__'))
    for name in ('src', 'native', 'bench'):
        shutil.copytree(root/name, args.run/'engine'/name, ignore=shutil.ignore_patterns('__pycache__'))
    (args.run/'inputs').mkdir()
    input_hashes={}
    for index,_ in CASES:
        dest=args.run/f'inputs/{index}.json'
        shutil.copy2(SOURCE/f'inputs/{index}/input.json', dest)
        input_hashes[index]=hashlib.sha256(dest.read_bytes()).hexdigest()
    tasks = [dict(index=i, direction=d, seed=s, repeat=r) for r in range(args.repeats)
             for i,d in CASES for s in range(args.seeds)]
    save(args.run/'tasks.json', tasks)
    save(args.run/'manifest.json', dict(parent_commit=subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
        source=str(SOURCE), old_engine=str(OLD_ENGINE), modes=args.modes, tasks=len(tasks),
        input_sha256=input_hashes,
        original_seed_count=10, tested_seed_indices=list(range(args.seeds)), repeats=args.repeats,
        iso_tolerance=1., event_tolerance=.5, branch_cap=100, explicit_H=True,
        seed_policy='unchanged cut_seed(cut), original ten-order generation; select specified seed index',
        cuts='uncut plus every source edge', whole_cache_bytes=64*2**20, extension_cache_bytes=64*2**20,
        baseline_native_sha256=hashlib.sha256(next((args.run/'baseline/src/rxn_core').glob('_engine*.so')).read_bytes()).hexdigest(),
        optimized_native_sha256=hashlib.sha256(next((args.run/'engine/src/rxn_core').glob('_engine*.so')).read_bytes()).hexdigest()))
    (args.run/'status').mkdir()
    print(args.run, len(tasks), 'paired tasks')


def submit(args):
    tasks=json.loads((args.run/'tasks.json').read_text())
    command=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','PYTHONHASHSEED=0',
        sys.executable,str(args.run/'engine/bench/conditioned_reuse_pilot.py'),'paired','--run',str(args.run),'--slot']
    options=['sbatch','--parsable','--partition=cpunodes','--nodes=1','--cpus-per-task=1','--mem=6G',
        '--time=00:10:00',f'--array=0-{len(tasks)-1}','--job-name=conditioned_reuse',
        f'--output={args.run}/status/%A_%a.out','--wrap',shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=options));print(job)


def paired(args):
    modes=json.loads((args.run/'manifest.json').read_text())['modes']
    shift=args.slot%len(modes);modes=modes[shift:]+modes[:shift]
    started=time.perf_counter()
    for mode in modes:
        if mode=='shared_symmetry':
            recover_one(args.run,args.slot,os.uname().nodename)
        if (args.run/f'results/{args.slot}/{mode}/summary.json').exists():
            print('SAVED',mode,flush=True);continue
        folder='baseline' if mode=='baseline' else 'engine'
        env=dict(os.environ,PYTHONPATH=f'{args.run}/{folder}/src:{args.run}/engine/bench',RXN_CORE_NATIVE='1')
        print('START',mode,flush=True)
        subprocess.run(['timeout','--kill-after=5s','300',sys.executable,__file__,'worker',
            '--run',str(args.run),'--slot',str(args.slot),'--mode',mode],env=env,check=True)
    save(args.run/f'pair_{args.slot}.json',dict(slot=args.slot,order=modes,host=os.uname().nodename,
        elapsed_including_io=time.perf_counter()-started,cpu_info=Path('/proc/cpuinfo').read_text()))


def worker(args):
    from rxn_core.aam import cut_seed
    from rxn_core.alignment.branch import _generate_seed_orders, find_islands
    from rxn_core.alignment.sweep import cut_sweep_items
    from rxn_core.artifacts import write_graph_checkpoint
    from rxn_core.cut_replay import FragmentRepair
    from rxn_core.frag import build_graph
    from rxn_core.matcher import _nauty_orbits
    from rxn_core.native_search import find_islands_native
    from rxn_core.search_symmetry import finalize_graph_symmetry, SymmetryWorkspace

    spec=json.loads((args.run/'tasks.json').read_text())[args.slot]
    folder=args.run/f'results/{args.slot}/{args.mode}';folder.mkdir(parents=True,exist_ok=False)
    raw=json.loads((args.run/f"inputs/{spec['index']}.json").read_text())
    left,right=(raw[k] for k in ('reactant','product'))
    if spec['direction']=='P_to_R':left,right=right,left
    phases={}
    def measure(label,fn):
        cpu=time.process_time();wall=time.perf_counter();result=fn()
        value=phases.setdefault(label,dict(cpu=0.,wall=0.))
        value['cpu']+=time.process_time()-cpu;value['wall']+=time.perf_counter()-wall
        return result
    def setup():
        r=build_graph(left['elements'],np.asarray(left['wbo']),bond_cut=.2)
        p=build_graph(right['elements'],np.asarray(right['wbo']),bond_cut=.2)
        po=_nauty_orbits(p,wbo_tol=1.)
        repair=None;workspace=None
        if args.mode not in ('baseline','current'):
            repair=FragmentRepair(r,p,po,extension_cache_bytes=64*2**20 if args.mode in ('extensions','combined') else 0)
        if args.mode in ('symmetry','combined'):
            from rxn_core.conditioned_symmetry import ConditionedSymmetryWorkspace
            workspace=ConditionedSymmetryWorkspace(p,1.)
        elif args.mode=='shared_symmetry':
            workspace=SymmetryWorkspace(p,1.)
        return r,p,po,repair,workspace
    source,target,orbits,repair,workspace=measure('setup',setup)
    cuts=cut_sweep_items(np.asarray(left['wbo']),.2)
    rows=[]
    for ordinal,cut in enumerate(cuts):
        def search():
            view=repair.for_cut(cut) if repair else None
            r=view.source if view else source.copy()
            if view is None:r.remove_edges_from(cut)
            ro=_nauty_orbits(r,wbo_tol=1.)
            order=_generate_seed_orders(r,10,rng_seed=cut_seed(cut))[spec['seed']]
            matcher=find_islands if args.mode in ('baseline','current') else find_islands_native
            return matcher(r,target,order,graph_floor=.2,iso_tol=1.,max_branches=100,
                r_orbits=ro,p_orbits=orbits,cuts=cut,growth_replay=view)
        graph=measure('search',search)
        graph,groups=measure('symmetry',lambda:finalize_graph_symmetry(graph,target,iso_tolerance=1.,workspace=workspace))
        digest=measure('audit_encoding',lambda:hashlib.sha256(json.dumps(graph.to_record(copy=False),
            sort_keys=True,separators=(',',':')).encode()).hexdigest())
        measure('persistence',lambda:write_graph_checkpoint(graph,folder/f'cut_{ordinal:04d}.pkl.gz'))
        rows.append(dict(cut=cut,digest=digest,states=len(graph.states),terminals=len(graph.terminals),
                         capped=graph.capped,groups=groups))
        save(folder/'progress.json',dict(completed=len(rows),total=len(cuts),phases=phases))
    save(folder/'summary.json',dict(**spec,mode=args.mode,name=raw['name'],phases=phases,cuts=rows,
        repair=None if repair is None else repair.stats(),symmetry=(None if workspace is None else
            dict(group_solves=workspace.coloring_count(),cached_partitions=len(workspace.coloring_cache))
            if args.mode=='shared_symmetry' else workspace.stats()),
        host=os.uname().nodename,peak_worker_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        compute_cpu=sum(phases[k]['cpu'] for k in ('setup','search','symmetry')),
        scope='one CPU; no queue/import/input-loading/persistence/audit encoding in compute; full graphs saved'))
    print(args.mode, 'DONE',sum(phases[k]['cpu'] for k in ('setup','search','symmetry')),flush=True)


def report(args):
    tasks=json.loads((args.run/'tasks.json').read_text())
    modes=json.loads((args.run/'manifest.json').read_text())['modes']
    totals={m:dict(compute_cpu=0.,search_cpu=0.,symmetry_cpu=0.,peak_worker_rss_kib=None,
                   rss_measurements=0) for m in modes}
    rows=[];mismatches=[];missing=[];checks=0
    for slot,spec in enumerate(tasks):
        found={}
        for mode in modes:
            path=args.run/f'results/{slot}/{mode}/summary.json'
            if not path.exists():missing.append((slot,mode));continue
            data=json.loads(path.read_text());found[mode]=data
            totals[mode]['compute_cpu']+=data['compute_cpu']
            for phase in ('search','symmetry'):totals[mode][phase+'_cpu']+=data['phases'][phase]['cpu']
            if data['peak_worker_rss_kib'] is not None:
                totals[mode]['peak_worker_rss_kib']=max(totals[mode]['peak_worker_rss_kib'] or 0,data['peak_worker_rss_kib'])
                totals[mode]['rss_measurements']+=1
        if 'baseline' in found:
            for mode,data in found.items():
                if mode=='baseline':continue
                for ordinal,(a,b) in enumerate(zip(found['baseline']['cuts'],data['cuts'],strict=True)):
                    checks+=1
                    if a['digest']!=b['digest']:mismatches.append((slot,mode,ordinal))
        rows.append(dict(**spec,slot=slot,methods={m:dict(compute_cpu=d['compute_cpu'],phases=d['phases'],
            peak_worker_rss_kib=d['peak_worker_rss_kib'],repair=d['repair'],symmetry=d['symmetry']) for m,d in found.items()}))
    result=dict(totals=totals,missing=missing,mismatches=mismatches,exact_cut_checks=checks,rows=rows)
    save(args.run/'analysis.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


def recover_one(run, slot, host):
    """Recover only the known post-persistence statistics-export failure."""
    from rxn_core.artifacts import read_graph_checkpoint
    folder=run/f'results/{slot}/shared_symmetry'
    if (folder/'summary.json').exists() or not (folder/'progress.json').exists():return False
    progress=json.loads((folder/'progress.json').read_text())
    if progress['completed']!=progress['total']:return False
    spec=json.loads((run/'tasks.json').read_text())[slot]
    rows=[]
    for ordinal in range(progress['total']):
        graph=read_graph_checkpoint(folder/f'cut_{ordinal:04d}.pkl.gz')
        rows.append(dict(cut=graph.contexts[0].cuts,
            digest=hashlib.sha256(json.dumps(graph.to_record(copy=False),sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            states=len(graph.states),terminals=len(graph.terminals),capped=graph.capped,groups=None))
    raw=json.loads((run/f"inputs/{spec['index']}.json").read_text())
    phases=progress['phases']
    save(folder/'summary.json',dict(**spec,mode='shared_symmetry',name=raw['name'],phases=phases,cuts=rows,
        repair=None,symmetry=None,host=host,peak_worker_rss_kib=None,
        compute_cpu=sum(phases[k]['cpu'] for k in ('setup','search','symmetry')),
        recovered_from_complete_checkpoints=True,
        scope='Original completed search phase timers; statistics exporter failed after persistence. No mapper rerun. RSS/cache counters unavailable.'))
    return True


def recover(args):
    """Serial recovery utility; large campaigns recover in parallel in paired()."""
    tasks=json.loads((args.run/'tasks.json').read_text())
    job=json.loads((args.run/'submission.json').read_text())['job']
    accounting=subprocess.check_output(['sacct','-j',job,'-nP','--format=JobID,State,NodeList'],text=True)
    hosts={int(row[0].split('_')[1]):row[2] for line in accounting.splitlines()
           if len(row:=line.split('|'))>=3 and '_' in row[0] and '.' not in row[0]}
    recovered=[slot for slot in range(len(tasks)) if recover_one(args.run,slot,hosts[slot])]
    save(args.run/'report_recovery.json',dict(slots=recovered,accounting=accounting))
    print('recovered',len(recovered),flush=True)


def resume(args):
    tasks=json.loads((args.run/'tasks.json').read_text())
    modes=json.loads((args.run/'manifest.json').read_text())['modes']
    original=json.loads((args.run/'submission.json').read_text())['job']
    accounting=subprocess.check_output(['sacct','-j',original,'-nP','--format=JobID,State,NodeList'],text=True)
    by_host={}
    for line in accounting.splitlines():
        row=line.split('|')
        if len(row)<3 or '_' not in row[0] or '.' in row[0]:continue
        slot=int(row[0].split('_')[1])
        if row[1] not in ('FAILED','COMPLETED'):continue
        if any(not (args.run/f'results/{slot}/{mode}/summary.json').exists() for mode in modes):
            by_host.setdefault(row[2],[]).append(slot)
    submissions=[]
    for host,slots in by_host.items():
        command=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','PYTHONHASHSEED=0',
            sys.executable,str(Path(__file__).resolve()),'paired','--run',str(args.run),'--slot']
        options=['sbatch','--parsable','--partition=cpunodes','--nodes=1','--nodelist',host,
            '--cpus-per-task=1','--mem=6G','--time=00:10:00','--array',','.join(map(str,slots)),
            '--job-name=conditioned_resume',f'--output={args.run}/status/resume_%A_%a.out',
            '--wrap',shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
        job=subprocess.check_output(options,text=True).strip()
        submissions.append(dict(job=job,host=host,slots=slots,command=options))
    save(args.run/'resume_submission.json',dict(submissions=submissions,driver=str(Path(__file__).resolve()),
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print(submissions)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','submit','paired','worker','report','recover','resume'))
    p.add_argument('--run',type=Path,required=True);p.add_argument('--slot',type=int)
    p.add_argument('--mode',choices=MODES);p.add_argument('--modes',choices=MODES,nargs='+',default=list(MODES))
    p.add_argument('--seeds',type=int,default=1);p.add_argument('--repeats',type=int,default=1)
    args=p.parse_args();globals()[args.command](args)

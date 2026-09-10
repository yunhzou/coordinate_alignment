"""Upstream SLAP plus blind single-edge sweeps, scored on original endpoints.

No AAM search or reference labels enter a mapper worker. All native alternatives
are checkpointed before the next cut; partial runs remain inspectable.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
from queue import SimpleQueue
import shlex
import shutil
import signal
import socket
import subprocess
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
ENV=BASE/'competitor_env_20260908'
PACKAGE=ENV/'lib/python3.10/site-packages/slapmapper'
DATASET=Path('/h/399/yunhengzou/coordinate_alignment/data/aam_benchmarks/golden_original_20260906')


def read(path):
    return json.loads(path.read_text())


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(path)


def records(path):
    if not path.exists():return []
    # An interrupted final append is not a committed cut; expected-cut counts
    # expose it as unfinished. Earlier native frames remain readable by offset.
    return [json.loads(line) for line in path.read_bytes().splitlines(keepends=True)
            if line.endswith(b'\n')]


def folder(run,task):
    return run/f"outputs/{task['case']}/{task['direction']}/{task['mode']}"


def cut_graphs(base,edge):
    from slapmapper.core import LabeledGraph
    graph={a:dict(nbrs) for a,nbrs in base[0].graph.items()}
    if edge is not None:
        a,b=edge
        del graph[a][b]
        del graph[b][a]
    # Rebuild initial WL labels. Mutating a copied LabeledGraph's adjacency
    # alone would retain stale, uncut initialization caches.
    return [LabeledGraph(graph,list(base[0].labels)),base[1].copy()]


def prepare(args):
    from rdkit import Chem
    rows=[json.loads(line) for line in (args.dataset/'audit.jsonl').read_text().splitlines()]
    args.run.mkdir(parents=True,exist_ok=False)
    (args.run/'status').mkdir()
    inputs=[dict(index=r['index'],input_reaction=r['input_reaction']) for r in rows]
    assert [r['index'] for r in inputs]==list(range(len(inputs)))
    save(args.run/'inputs.json',inputs)
    tasks=[]
    for row in inputs:
        mols=[Chem.AddHs(Chem.MolFromSmiles(s)) for s in row['input_reaction'].split('>>')]
        assert all(a.GetAtomMapNum()==0 for mol in mols for a in mol.GetAtoms())
        for reverse in (False,True):
            mol=mols[int(reverse)]
            for mode in ('binary','weighted'):
                tasks.append(dict(slot=len(tasks),case=row['index'],
                    direction='P_to_R' if reverse else 'R_to_P',mode=mode,
                    expected_cuts=1+mol.GetNumBonds(),
                    schedule_weight=(1+mol.GetNumBonds())*max(m.GetNumAtoms() for m in mols)**2))
    save(args.run/'tasks.json',tasks)
    ordered=sorted(tasks,key=lambda t:(-t['schedule_weight'],t['slot']))
    save(args.run/'batches.json',[[t['slot'] for t in ordered[i::args.batches]] for i in range(args.batches)])
    shutil.copytree(PACKAGE,args.run/'vendor/slapmapper',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(__file__,args.run/'slap_edge_sweep.py')
    for name in ('golden_competitors.py','golden_evaluation.py'):
        shutil.copy2(ROOT/'bench'/name,args.run/name)
    shutil.copy2(ROOT/'reports/real_ts_comparison_20260909/SLAP_LICENSE',args.run/'vendor/SLAP_LICENSE')
    save(args.run/'manifest.json',dict(cases=len(rows),tasks=len(tasks),
        expected_calls=sum(t['expected_cuts'] for t in tasks),dataset=str(args.dataset.resolve()),
        audit_sha256=hashlib.sha256((args.dataset/'audit.jsonl').read_bytes()).hexdigest(),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        upstream_commit='ea248fd9494f52f4865193e87a98cc92c62b5f9e',
        vendor_sha256={str(p.relative_to(args.run)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (args.run/'vendor').rglob('*.py')},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        python=str(ENV/'bin/python'),workers_per_node=args.workers,batches=args.batches,
        seed=20260910,hash_seed=0,watchdog_seconds=300,add_Hs=True,break_sym='heavy',
        cut_scope='Uncut plus every single source edge, including H bonds; no atoms removed',
        objective_scope='SLAP optimizes each perturbed graph; outputs restored to original endpoints for evaluation',
        balancing='Unchanged upstream smiles2lgp element-balancing placeholders',
        native_archive='native.bin: independent gzip/pickle frames; records.jsonl: committed offsets and lengths',
        evaluation='Original strict heavy-atom relation modulo endpoint chemical symmetry, including unmatched atoms',
        scheduling='Largest label-free structural work estimates distributed across nodes; no search-policy effect'))
    print(json.dumps(dict(cases=len(rows),tasks=len(tasks),calls=sum(t['expected_cuts'] for t in tasks))),flush=True)


def case(args):
    import random
    import numpy as np
    from slapmapper.core import SlapMapper
    from slapmapper.aam._smiles import smiles2lgp,get_numbered_rxn_smiles
    random.seed(20260910);np.random.seed(20260910)
    task=read(args.run/'tasks.json')[args.slot]
    original=read(args.run/'inputs.json')[task['case']]['input_reaction']
    reverse=task['direction']=='P_to_R'
    reaction='>>'.join(original.split('>>')[::-1]) if reverse else original
    destination=folder(args.run,task)
    destination.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter()
    base=smiles2lgp(reaction,add_Hs=True)
    edges=sorted((a,b) for a,nbrs in base[0].graph.items() for b in nbrs if a<b)
    assert len(edges)+1==task['expected_cuts']
    targets=[i for i,z in enumerate(base[0].props['atomic numbers']) if z>1]
    save(destination/'input.json',dict(**task,reaction=reaction,original_reaction=original,
        input_graphs=[dict(labels=g.labels,edges=[(a,b,w) for a,n in g.graph.items() for b,w in n.items() if a<b],
                          atom_counts=g.props['natoms slices']) for g in base],
        preparation_wall=time.perf_counter()-start,host=socket.gethostname(),
        affinity=sorted(os.sched_getaffinity(0)),slurm_job=os.environ.get('SLURM_JOB_ID')))
    mapper=SlapMapper(binary=task['mode']=='binary')
    with (destination/'records.jsonl').open('wb') as journal,(destination/'native.bin').open('wb') as archive:
        for ordinal,edge in enumerate([None,*edges]):
            cpu,wall=time.process_time(),time.perf_counter()
            graphs=cut_graphs(base,edge)
            graph_cpu,graph_wall=time.process_time()-cpu,time.perf_counter()-wall
            cpu,wall=time.process_time(),time.perf_counter()
            try:
                mapper.get_maps(graphs,break_sym_targets=targets)
                status,error='mapped',None
            except Exception:
                status,error='mapping_error',traceback.format_exc()
            mapping_cpu,mapping_wall=time.process_time()-cpu,time.perf_counter()-wall
            cpu,wall=time.process_time(),time.perf_counter()
            candidates=[];native=[]
            if status=='mapped':
                for result in mapper.results:
                    labels=[list(g.labels) for g in result['lgp']]
                    mapped=get_numbered_rxn_smiles(reaction,labels)
                    if reverse:mapped='>>'.join(mapped.split('>>')[::-1])
                    candidates.append(dict(mapped_rxn=mapped,native_cost=float(result['val'])))
                    native.append(dict(labels=labels,val=result['val'],lap_sols=result['lap_sols']))
            export_cpu,export_wall=time.process_time()-cpu,time.perf_counter()-wall
            cpu,wall=time.process_time(),time.perf_counter()
            blob=gzip.compress(pickle.dumps(native,protocol=5),compresslevel=1,mtime=0)
            encoding_cpu,encoding_wall=time.process_time()-cpu,time.perf_counter()-wall
            offset=archive.tell();archive.write(blob);archive.flush()
            row=dict(ordinal=ordinal,cut=edge,status=status,error=error,candidates=candidates,
                graph_cpu=graph_cpu,graph_wall=graph_wall,mapping_cpu=mapping_cpu,mapping_wall=mapping_wall,
                export_cpu=export_cpu,export_wall=export_wall,encoding_cpu=encoding_cpu,encoding_wall=encoding_wall,
                native=dict(offset=offset,length=len(blob),sha256=hashlib.sha256(blob).hexdigest()))
            journal.write((json.dumps(row)+'\n').encode());journal.flush()
    save(destination/'complete.json',dict(**task,complete=True,cuts=len(edges)+1))


def batch(args):
    manifest=read(args.run/'manifest.json');tasks=read(args.run/'tasks.json')
    cores=SimpleQueue()
    affinity=sorted(os.sched_getaffinity(0))
    assert len(affinity)>=manifest['workers_per_node']
    for core in affinity[:manifest['workers_per_node']]:cores.put(core)
    env=dict(os.environ,PYTHONPATH=str(args.run/'vendor'),PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',
             OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def one(slot):
        path=args.run/f'status/{slot}.json'
        if path.exists():return
        core=cores.get();started=time.perf_counter()
        save(path,dict(slot=slot,started=time.time(),host=socket.gethostname(),core=core))
        try:
            with (args.run/f'status/{slot}.log').open('w') as log:
                p=subprocess.Popen(['taskset','-c',str(core),manifest['python'],str(args.run/'slap_edge_sweep.py'),
                    'case','--run',str(args.run),'--slot',str(slot)],env=env,stdout=log,stderr=subprocess.STDOUT,
                    start_new_session=True)
                try:code=p.wait(timeout=manifest['watchdog_seconds'])
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid,signal.SIGKILL);p.wait();code=124
            save(path,dict(slot=slot,exit=code,finished=time.time(),host=socket.gethostname(),
                elapsed_including_startup_io=time.perf_counter()-started,
                completed_cuts=len(records(folder(args.run,tasks[slot])/'records.jsonl')),
                expected_cuts=tasks[slot]['expected_cuts']))
        finally:cores.put(core)
    with ThreadPoolExecutor(max_workers=manifest['workers_per_node']) as pool:
        list(pool.map(one,read(args.run/'batches.json')[args.slot]))


def submit(args):
    m=read(args.run/'manifest.json')
    command=['sbatch','--parsable','--partition=cpunodes_nia','--exclude=bosque49,bosque56',
        '--nodes=1',f"--cpus-per-task={m['workers_per_node']}",'--mem=96G','--time=00:45:00','--no-requeue',
        f"--array=0-{m['batches']-1}%{m['batches']}",'--job-name=slap_edge_sweep',
        f'--output={args.run}/status/%A_%a.out','--wrap',
        shlex.join([m['python'],str(args.run/'slap_edge_sweep.py'),'batch','--run',str(args.run),'--slot'])+
        ' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(command,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=command));print(job,flush=True)


def evaluate(args):
    from golden_competitors import signatures
    manifest=read(args.run/'manifest.json')
    audit_path=Path(manifest['dataset'])/'audit.jsonl'
    assert hashlib.sha256(audit_path.read_bytes()).hexdigest()==manifest['audit_sha256']
    audit=[json.loads(x) for x in audit_path.read_text().splitlines()]
    tasks=read(args.run/'tasks.json')
    for index in range(args.slot,len(audit),args.shards):
        expected_endpoints,expected=signatures(audit[index]['mapped_reaction'])
        attempts=[];classes=set();uncut=set();cut_only=set();cpu=0.;workflow=0.;calls=0;errors=0;invalid=0
        for task in tasks[4*index:4*index+4]:
            rows=records(folder(args.run,task)/'records.jsonl');hits=[];direction_classes=set()
            for row in rows:
                calls+=1;cpu+=row['mapping_cpu'];workflow+=row['mapping_cpu']+row['graph_cpu']+row['export_cpu']
                errors+=row['status']!='mapped'
                for number,candidate in enumerate(row['candidates']):
                    try:
                        endpoints,actual=signatures(candidate['mapped_rxn'])
                        if endpoints!=expected_endpoints:raise ValueError('Endpoint chemistry changed')
                    except Exception:
                        invalid+=1;continue
                    classes.add(actual);direction_classes.add(actual)
                    (uncut if row['cut'] is None else cut_only).add(actual)
                    if actual==expected:hits.append(dict(cut=row['cut'],candidate=number,ordinal=row['ordinal']))
            attempts.append(dict(**task,completed_cuts=len(rows),hits=hits,classes=len(direction_classes)))
        save(args.run/f'evaluations/{index}.json',dict(index=index,recovered=expected in classes,
            uncut_recovered=expected in uncut,cut_recovered=expected in cut_only,
            unique_classes=len(classes),uncut_classes=len(uncut),attempts=attempts,
            completed_calls=calls,mapping_errors=errors,invalid_predictions=invalid,
            mapping_cpu=cpu,workflow_cpu_excluding_io=workflow))


def summary(args):
    m=read(args.run/'manifest.json');rows=[read(p) for p in (args.run/'evaluations').glob('*.json')]
    baselines=[read(ROOT/f'reports/golden_slap_budget_20260908/{name}_summary.json')
               for name in ('expanded','all_atoms')]
    missing_sets=[set(v['unresolved_cases']) for summary in baselines for v in summary.values()
                  if isinstance(v,dict) and 'unresolved_cases' in v]
    old_missing=set.intersection(*missing_sets)
    recovered={r['index'] for r in rows if r['recovered']}
    result=dict(cases=m['cases'],evaluated=len(rows),uncut_recovered=sum(r['uncut_recovered'] for r in rows),
        sweep_union_recovered=len(recovered),expanded_baseline_recovered=m['cases']-len(old_missing),
        combined_with_expanded_baseline=m['cases']-len(old_missing-recovered),
        new_cases_beyond_expanded_baseline=sorted(old_missing & recovered),
        completed_calls=sum(r['completed_calls'] for r in rows),expected_calls=m['expected_calls'],
        mapping_cpu=sum(r['mapping_cpu'] for r in rows),workflow_cpu_excluding_io=sum(r['workflow_cpu_excluding_io'] for r in rows),
        mapping_errors=sum(r['mapping_errors'] for r in rows),invalid_predictions=sum(r['invalid_predictions'] for r in rows),
        unrecovered=sorted(set(range(m['cases']))-recovered),
        scope='Fixed full denominator; alternatives scored on original endpoints; incomplete cuts remain untested')
    save(args.run/'summary.json',result);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','case','batch','submit','evaluate','summary'))
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--dataset',type=Path,default=DATASET)
    p.add_argument('--slot',type=int,default=0)
    p.add_argument('--shards',type=int,default=1)
    p.add_argument('--workers',type=int,default=32)
    p.add_argument('--batches',type=int,default=32)
    args=p.parse_args();globals()[args.command](args)

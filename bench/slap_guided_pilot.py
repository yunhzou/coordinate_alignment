"""Exploratory cheap-proposal / no-cut AAM ablation. Never changes production search."""
import argparse
from collections import Counter
from dataclasses import replace
import gzip
import importlib.util
import json
from pathlib import Path
import pickle
import random
import shutil
import subprocess
import shlex
import sys
import time
from types import SimpleNamespace

import numpy as np
from rxn_core import AAMProblem, AAMSearchConfig
from rxn_core.domain import MolecularEndpoint
from rxn_core.aam import _initialize_search, _search_cut
from rxn_core.alignment.branch import find_islands
from rxn_core.frag import build_graph
from rxn_core.matcher import _nauty_orbits
from rxn_core.search_symmetry import finalize_graph_symmetry
from rxn_core.artifacts import write_graph_checkpoint
from compare_elementary_outputs import event_counts, events_row, features, certificate, refine

SOURCE=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_tol1_20260909')
SLAP=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/competitor_env_20260908/lib/python3.10/site-packages/slapmapper/core.py')


def save(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2)+'\n')


def propose(problem,raw,run,index):
    from rdkit import Chem
    spec=importlib.util.spec_from_file_location('slap_core_pilot',run/'slap_core.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    graphs=[]
    for endpoint in (problem.reactant,problem.product):
        w=endpoint.wbo;n=len(w)
        graphs.append(mod.LabeledGraph({i:{j:1 for j in range(n) if i!=j and w[i,j]>.2} for i in range(n)},
            [Chem.GetPeriodicTable().GetAtomicNumber(e) for e in endpoint.elements]))
    start=time.process_time();wall=time.perf_counter()
    mapper=mod.SlapMapper(binary=True)
    mapper.get_maps(graphs,break_sym_targets=[i for i,e in enumerate(problem.reactant.elements) if e!='H'])
    cpu=time.process_time()-start;elapsed=time.perf_counter()-wall
    records=[];native=[]
    for i,result in enumerate(mapper.results):
        left,right=result['lgp'];mapping={a:b for label,atoms in left.label2idxs.items() for a,b in zip(atoms,right.label2idxs[label])}
        assert all(len(atoms)==1 or all(problem.reactant.elements[a]=='H' for a in atoms) for atoms in left.label2idxs.values())
        vector=[mapping[a] for a in range(problem.atom_count)]
        records.append(dict(candidate=i,mapping=sorted(mapping.items()),
            events=events_row(event_counts(problem.reactant.wbo,problem.product.wbo,[vector])[0]),
            all_h_label_permutations_score_invariant=False))
        native.append(dict(graphs=[dict(labels=list(g.labels)) for g in (left,right)],native_cost=float(result['val'])))
    # Preserve internal LAP alternatives as well as label pools.
    with gzip.open(run/f'proposal/{index}.native.pkl.gz','wb') as f:
        pickle.dump([dict(val=r['val'],lap_sols=r.get('lap_sols'),labels=[list(g.labels) for g in r['lgp']]) for r in mapper.results],f)
    save(run/f'slap_xyz/{index}.json',dict(candidates=native))
    save(run/f'proposal/{index}.json',dict(slap=records,scope='Direct SLAP on WBO>0.2 binary connectivity; no XYZ perception or SMILES.'))
    start=time.process_time();wall=time.perf_counter()
    refine(SimpleNamespace(run=run/'proposal',source=run,index=index))
    refinement_cpu=time.process_time()-start;refinement_wall=time.perf_counter()-wall
    refined=json.loads((run/f'proposal/{index}.refined.json').read_text())
    best=min(refined['slap'],key=lambda r:r['events']['total'])
    return dict(best=best,all=refined['slap'],search_cpu=cpu,search_wall=elapsed,
                refinement_cpu=refinement_cpu,refinement_wall=refinement_wall)


def guided_order(problem,mapping):
    """Soft proposal only: grow first away from its putative edit endpoints."""
    r,p=problem.reactant.wbo,problem.product.wbo;n=len(r);penalty=np.zeros(n)
    for a in range(n):
        for b in range(a+1,n):
            x,y=r[a,b],p[mapping[a],mapping[b]]
            if (x>.2)!=(y>.2) or x>.2 and y>.2 and abs(x-y)>.5:
                penalty[a]+=1;penalty[b]+=1
    order=list(range(n));random.Random(42).shuffle(order)
    return sorted(order,key=lambda a:(penalty[a],-int(np.count_nonzero(r[a]>.2))))


def case(args):
    raw=json.loads((args.run/f'inputs/{args.index}/input.json').read_text())
    original=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    feat=features(raw);heavy=[i for i,e in enumerate(original.reactant.elements) if e!='H']
    baseline=json.loads((args.source/f'comparison/{args.index}.refined.json').read_text())
    if args.index==123:
        baseline=json.loads((args.source/'case123_cap200/comparison/123.refined.json').read_text())
    extra=json.loads((args.source/f'overlap/{args.index}.json').read_text())['extra_aam_pattern_absent_from_saved_slap']
    reference=[baseline['aam_best_example']]+([extra] if extra else [])
    refcert=[certificate(feat,[dict(r['mapping'])[i] for i in range(original.atom_count)],heavy) for r in reference]
    proposal=propose(original,raw,args.run,args.index)
    config=AAMSearchConfig(seed_count=10,branch_limit=100,iso_tolerance=1.)
    results=[]
    for direction in ('R_to_P','P_to_R'):
        reverse=direction=='P_to_R';problem=AAMProblem(original.product,original.reactant) if reverse else original
        for policy,seeds in [('random1',1),('random10',10),('guided1',1)]:
            start=time.process_time();wall=time.perf_counter()
            cfg=replace(config,seed_count=seeds)
            if policy=='guided1':
                mapping=dict(proposal['best']['mapping'])
                if reverse:mapping={b:a for a,b in mapping.items()}
                order=guided_order(problem,mapping)
                r=build_graph(problem.reactant.elements,problem.reactant.wbo,bond_cut=.2)
                p=build_graph(problem.product.elements,problem.product.wbo,bond_cut=.2)
                graph=find_islands(r,p,order,graph_floor=.2,iso_tol=1.,max_branches=100,
                    p_orbits=_nauty_orbits(p,wbo_tol=1.),r_orbits=_nauty_orbits(r,wbo_tol=1.),cuts=())
            else:
                _initialize_search(problem,cfg);graph,_=_search_cut(())
                p=build_graph(problem.product.elements,problem.product.wbo,bond_cut=.2)
            search_cpu=time.process_time()-start;search_wall=time.perf_counter()-wall
            start=time.process_time();wall=time.perf_counter()
            graph,_=finalize_graph_symmetry(graph,p,iso_tolerance=1.)
            symmetry_cpu=time.process_time()-start;symmetry_wall=time.perf_counter()-wall
            folder=args.run/f'cases/{args.index}/{policy}/{direction}';folder.mkdir(parents=True)
            write_graph_checkpoint(graph,folder/'graph.pkl.gz')
            start=time.process_time();vectors=[]
            for t in graph.terminals:
                m=dict(graph.states[t].mapping)
                if len(m)!=original.atom_count:continue
                if reverse:m={b:a for a,b in m.items()}
                assert sorted(m.values())==list(range(original.atom_count))
                vectors.append([m[i] for i in range(original.atom_count)])
            scores=event_counts(original.reactant.wbo,original.product.wbo,vectors) if vectors else []
            certs={};best=None;bestmap=None
            for v,score in zip(vectors,scores):
                key=tuple(v[i] for i in heavy)
                if key not in certs:certs[key]=certificate(feat,v,heavy)
                total=int(sum(score))
                if best is None or total<best:best=total;bestmap=v
            row=dict(policy=policy,direction=direction,best=best,mapping=bestmap,
                full_witnesses=len(vectors),heavy_witness_classes=len(set(certs.values())),
                baseline_best_pattern=refcert[0] in certs.values(),
                baseline_extra_pattern=refcert[1] in certs.values() if len(refcert)>1 else None,
                capped=graph.capped,search_cpu=search_cpu,search_wall=search_wall,
                symmetry_cpu=symmetry_cpu,symmetry_wall=symmetry_wall,
                evaluation_cpu=time.process_time()-start)
            save(folder/'summary.json',row);results.append(row)
    save(args.run/f'{args.index}.json',dict(index=args.index,name=raw['name'],baseline_best=baseline['aam_min_saved_events'],
        proposal=proposal,results=results,scope='No-cut proposal-guided search pilot, not full-family equivalence. No hard anchors, no production algorithm edits.'))


def submit(args):
    args.run.mkdir(parents=True,exist_ok=False)
    (args.run/'proposal').mkdir();(args.run/'status').mkdir()
    shutil.copytree(args.source/'inputs',args.run/'inputs')
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(SLAP,args.run/'slap_core.py')
    env=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','PYTHONHASHSEED=0',
         f'PYTHONPATH={args.run}/dependencies:{args.run}/engine/src:{args.run}/engine/bench']
    cmd=[*env,'timeout','--kill-after=5s','300',sys.executable,str(args.run/'engine/bench/slap_guided_pilot.py'),
         'case','--run',str(args.run),'--source',str(args.source),'--index']
    options=['sbatch','--parsable','--partition=cpunodes,cpunodes_nia','--exclude=bosque5,bosque6,bosque8',
        '--cpus-per-task=1','--mem=8G','--time=00:10:00','--array=0-139%32','--job-name=slap_pilot',
        f'--output={args.run}/status/%A_%a.out','--wrap',shlex.join(cmd)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=options,source=str(args.source),
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()))
    print(job)


def resume(args):
    previous=json.loads((args.run/'submission.json').read_text())
    slots=[i for i in range(140) if not (args.run/f'{i}.json').exists()]
    cmd=[('--array='+','.join(map(str,slots))+'%32' if s.startswith('--array=') else
          f'PYTHONPATH={args.run}/dependencies:{args.run}/engine/src:{args.run}/engine/bench' if s.startswith('PYTHONPATH=') else s)
         for s in previous['command']]
    # --wrap contains a shell command, so update its explicit environment assignment too.
    cmd[-1]=cmd[-1].replace(f'PYTHONPATH={args.run}/engine/src:',f'PYTHONPATH={args.run}/dependencies:{args.run}/engine/src:')
    job=subprocess.check_output(cmd,text=True).strip()
    save(args.run/'resubmission.json',dict(job=job,command=cmd,slots=slots,
         reason='Initial tasks failed importing absent scipy before matching; isolated SciPy installed, smoke case passed.'))
    print(job)


def report(args):
    rows=[json.loads((args.run/f'{i}.json').read_text()) for i in range(140)]
    outputs={}
    for policy in ('random1','random10','guided1'):
        for hybrid in (False,True):
            name=('slap_plus_' if hybrid else '')+policy;details=[]
            for row in rows:
                found=[r for r in row['results'] if r['policy']==policy]
                costs=[r['best'] for r in found if r['best'] is not None]
                proposal=row['proposal']
                if hybrid:costs.append(proposal['best']['events']['total'])
                best=min(costs) if costs else None
                cpu=sum(r['search_cpu']+r['symmetry_cpu'] for r in found)
                # Guidance must pay proposal search and hydrogen-score refinement.
                if hybrid or policy=='guided1':cpu+=proposal['search_cpu']+proposal['refinement_cpu']
                details.append(dict(index=row['index'],best=best,baseline=row['baseline_best'],cpu=cpu,
                    best_pattern_seen=any(r['baseline_best_pattern'] for r in found),
                    extra_pattern_seen=any(r['baseline_extra_pattern'] for r in found)))
            # Pattern retention below concerns AAM witnesses only, even in hybrid mode.
            outputs[name]=dict(cases=140,full=sum(r['best'] is not None for r in details),
                equal_or_better=sum(r['best'] is not None and r['best']<=r['baseline'] for r in details),
                within_plus1=sum(r['best'] is not None and r['best']<=r['baseline']+1 for r in details),
                worse_cases=[r['index'] for r in details if r['best'] is None or r['best']>r['baseline']],
                total_cpu=sum(r['cpu'] for r in details),
                aam_best_witness_pattern_retained=sum(r['best_pattern_seen'] for r in details),
                aam_extra_witness_pattern_retained=sum(r['extra_pattern_seen'] for r in details),
                per_case=details)
    outputs['slap_direct']=dict(equal_or_better=sum(r['proposal']['best']['events']['total']<=r['baseline_best'] for r in rows),
        total_cpu=sum(r['proposal']['search_cpu']+r['proposal']['refinement_cpu'] for r in rows),
        unresolved_hydrogen_pools=sum(not c['hydrogen_score_optimization']['optimal'] for r in rows for c in r['proposal']['all']))
    phases=Counter()
    for file in args.source.glob('directions/*/*/search.json'):
        for name,value in json.loads(file.read_text())['phases'].items():phases[name]+=value['cpu_seconds']
    save(args.run/'summary.json',dict(policies=outputs,full_search_phases_cpu=dict(phases),
        scope='Exploratory workload reductions, not equally complete acceleration. AAM pattern retention checks one saved baseline-best witness and one saved extra per case, not exhaustive family containment. Hybrid score union does not merge compressed data structures. Initial dependency failures excluded; complete retry outputs retained.'))
    print(json.dumps({k:{a:b for a,b in v.items() if a!='per_case'} for k,v in outputs.items()},indent=2))


def relocate(args):
    old=json.loads((args.run/('refinement_submission.json' if args.refinement else 'resubmission.json')).read_text())
    jobs=subprocess.check_output(['squeue','-h','-j',old['job'],'-t','CONFIGURING','-o','%i'],text=True).split()
    slots=[]
    for job in jobs:
        i=int(job.split('_')[1])
        if args.refinement:
            assert not (args.run/f'local/{i}.json').exists()
            assert not (args.run/f'local/{i}').exists()
        else:
            assert not (args.run/f'proposal/{i}.refined.json').exists()
            assert not (args.run/f'cases/{i}').exists()
        subprocess.run(['scancel',job],check=True);slots.append(i)
    assert slots
    command=[('--array='+','.join(map(str,slots)) if s.startswith('--array=') else
              '--exclude=bosque5,bosque6,bosque7,bosque8' if s.startswith('--exclude=') else s) for s in old['command']]
    job=subprocess.check_output(command,text=True).strip()
    save(args.run/('refinement_relocation.json' if args.refinement else 'relocation.json'),dict(job=job,command=command,unstarted_slots=slots))
    print(job)


def refine_local(args):
    pilot=json.loads((args.run/f'{args.index}.json').read_text())
    raw=json.loads((args.run/f'inputs/{args.index}/input.json').read_text())
    original=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    proposal_maps=[dict(pilot['proposal']['best']['mapping'])]
    proposal_maps += [dict(enumerate(r['mapping'])) for r in pilot['results'] if r['policy']=='random1' and r['mapping'] is not None]
    cfg=AAMSearchConfig(seed_count=1,branch_limit=100,iso_tolerance=1.)
    rows=[]
    for direction in ('R_to_P','P_to_R'):
        reverse=direction=='P_to_R';problem=AAMProblem(original.product,original.reactant) if reverse else original
        direct=set();centers=set();r,p=problem.reactant.wbo,problem.product.wbo
        for physical in proposal_maps:
            m={b:a for a,b in physical.items()} if reverse else physical
            for a in range(len(r)):
                for b in range(a+1,len(r)):
                    x,y=r[a,b],p[m[a],m[b]]
                    if (x>.2)!=(y>.2) or x>.2 and y>.2 and abs(x-y)>.5:
                        centers.update((a,b))
                        if x>.2:direct.add((a,b))
        incident={(a,b) for a in range(len(r)) for b in range(a+1,len(r)) if r[a,b]>.2 and (a in centers or b in centers)}
        target=build_graph(problem.product.elements,p,bond_cut=.2)
        _initialize_search(problem,cfg)
        for a,b in sorted(incident):
            start=time.process_time();wall=time.perf_counter()
            graph,_=_search_cut(((a,b),))
            search_cpu=time.process_time()-start
            start=time.process_time();graph,_=finalize_graph_symmetry(graph,target,iso_tolerance=1.)
            symmetry_cpu=time.process_time()-start;elapsed=time.perf_counter()-wall
            vectors=[]
            for t in graph.terminals:
                m=dict(graph.states[t].mapping)
                if len(m)!=original.atom_count:continue
                if reverse:m={b:a for a,b in m.items()}
                vectors.append([m[i] for i in range(original.atom_count)])
            scores=event_counts(original.reactant.wbo,original.product.wbo,vectors) if vectors else []
            index=min(range(len(vectors)),key=lambda i:sum(scores[i])) if vectors else None
            folder=args.run/f'local/{args.index}/{direction}';folder.mkdir(parents=True,exist_ok=True)
            write_graph_checkpoint(graph,folder/f'{a}_{b}.pkl.gz')
            rows.append(dict(direction=direction,cut=[a,b],direct=(a,b) in direct,
                best=None if index is None else int(sum(scores[index])),mapping=None if index is None else vectors[index],
                capped=graph.capped,search_cpu=search_cpu,symmetry_cpu=symmetry_cpu,elapsed=elapsed))
    save(args.run/f'local/{args.index}.json',dict(index=args.index,cuts=rows,
        scope='One seed per selected cut, both directions. Cuts from changed bonds or source bonds incident to proposed event atoms; no reference used for selection.'))


def submit_refinement(args):
    (args.run/'local').mkdir()
    shutil.copy2(__file__,args.run/'refinement_driver.py')
    cmd=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','PYTHONHASHSEED=0',
         f'PYTHONPATH={args.run}/engine/src:{args.run}/engine/bench','timeout','--kill-after=5s','300',
         sys.executable,str(args.run/'refinement_driver.py'),'refine_local','--run',str(args.run),'--index']
    options=['sbatch','--parsable','--partition=cpunodes,cpunodes_nia','--exclude=bosque5,bosque6,bosque8',
        '--cpus-per-task=1','--mem=8G','--time=00:10:00','--array=0-139%32','--job-name=slap_local',
        f'--output={args.run}/local/slurm_%A_%a.out','--wrap',shlex.join(cmd)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip();save(args.run/'refinement_submission.json',dict(job=job,command=options));print(job)


def report_local(args):
    from rxn_core.family_scoring import bond_events
    summary=json.loads((args.run/'summary.json').read_text());policies={}
    audit=0
    for policy in ('changed_edges','incident_edges'):
        total=summary['policies']['slap_plus_random1']['total_cpu'];hits=0;cuts=0;details=[]
        for index in range(140):
            row=json.loads((args.run/f'{index}.json').read_text())
            local=json.loads((args.run/f'local/{index}.json').read_text())
            chosen=[r for r in local['cuts'] if policy=='incident_edges' or r['direct']]
            base=[r['best'] for r in row['results'] if r['policy']=='random1' and r['best'] is not None]
            scores=base+[row['proposal']['best']['events']['total']]+[r['best'] for r in chosen if r['best'] is not None]
            best=min(scores);hits+=best<=row['baseline_best'];cuts+=len(chosen)
            cpu=sum(r['search_cpu']+r['symmetry_cpu'] for r in chosen);total+=cpu
            details.append(dict(index=index,best=best,baseline=row['baseline_best'],additional_cpu=cpu,cuts=len(chosen)))
        policies[policy]=dict(equal_or_better=hits,total_search_group_proposal_cpu=total,selected_directional_cuts=cuts,
                              missing=[r['index'] for r in details if r['best']>r['baseline']],per_case=details)
    for index in range(140):
        raw=json.loads((args.run/f'inputs/{index}/input.json').read_text())
        problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
        row=json.loads((args.run/f'{index}.json').read_text());local=json.loads((args.run/f'local/{index}.json').read_text())
        records=[(r['mapping'],r['best']) for r in row['results'] if r['mapping'] is not None]
        records += [(r['mapping'],r['best']) for r in local['cuts'] if r['mapping'] is not None]
        records += [([dict(r['mapping'])[i] for i in range(problem.atom_count)],r['events']['total']) for r in row['proposal']['all']]
        for vector,expected in records:
            m=dict(enumerate(vector));assert sorted(m.values())==list(range(problem.atom_count))
            assert all(problem.reactant.elements[a]==problem.product.elements[b] for a,b in m.items())
            assert bond_events(problem,m)['total']==expected;audit+=1
    save(args.run/'local_summary.json',dict(policies=policies,independently_rescored_witnesses=audit,
        scope='All 140 cases; proposal-directed cuts chosen without baseline answers. Kernel CPU includes matching, group finalization and SLAP search/H-refinement; excludes artifact IO, benchmark certificate evaluation, import/setup and local cut-list construction. Not equal-completeness or end-to-end latency claims.'))
    print(json.dumps({k:{a:b for a,b in v.items() if a!='per_case'} for k,v in policies.items()},indent=2));print('audited',audit)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('case','submit','resume','report','relocate','refine_local','submit_refinement','report_local'))
    p.add_argument('--run',required=True,type=Path);p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--index',type=int)
    p.add_argument('--refinement',action='store_true')
    args=p.parse_args();globals()[args.command](args)

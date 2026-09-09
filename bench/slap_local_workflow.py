"""Fresh end-to-end experimental proposal + one-seed local-cut workflow."""
import argparse
import json
from pathlib import Path
import shutil
import shlex
import subprocess
import sys
import time

from slap_guided_pilot import propose, save, SOURCE
from rxn_core import AAMProblem, AAMSearchConfig
from rxn_core.domain import MolecularEndpoint
from rxn_core.aam import _initialize_search, _search_cut
from rxn_core.frag import build_graph
from rxn_core.search_symmetry import finalize_graph_symmetry
from rxn_core.artifacts import write_graph_checkpoint
from compare_elementary_outputs import event_counts


def case(args):
    raw=json.loads((args.run/f'inputs/{args.index}/input.json').read_text())
    problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    start_cpu=time.process_time();start_wall=time.perf_counter()
    proposal=propose(problem,raw,args.run,args.index)
    cfg=AAMSearchConfig(seed_count=1,branch_limit=100,iso_tolerance=1.)
    results=[];proposals=[dict(proposal['best']['mapping'])]
    def run_cut(direction,cut):
        reverse=direction=='P_to_R';search_problem=AAMProblem(problem.product,problem.reactant) if reverse else problem
        graph,_=_search_cut(cut)
        target=build_graph(search_problem.product.elements,search_problem.product.wbo,bond_cut=.2)
        graph,_=finalize_graph_symmetry(graph,target,iso_tolerance=1.)
        mappings=[]
        for t in graph.terminals:
            m=dict(graph.states[t].mapping)
            if len(m)!=problem.atom_count:continue
            if reverse:m={b:a for a,b in m.items()}
            mappings.append([m[i] for i in range(problem.atom_count)])
        counts=event_counts(problem.reactant.wbo,problem.product.wbo,mappings) if mappings else []
        best=min(range(len(mappings)),key=lambda i:sum(counts[i])) if mappings else None
        folder=args.run/f'graphs/{args.index}/{direction}';folder.mkdir(parents=True,exist_ok=True)
        name='none' if not cut else '_'.join(map(str,cut[0]))
        write_graph_checkpoint(graph,folder/f'{name}.pkl.gz')
        row=dict(direction=direction,cut=cut,capped=graph.capped,
                 best=None if best is None else int(sum(counts[best])),mapping=None if best is None else mappings[best])
        results.append(row)
        return row
    for direction in ('R_to_P','P_to_R'):
        search_problem=AAMProblem(problem.product,problem.reactant) if direction=='P_to_R' else problem
        _initialize_search(search_problem,cfg)
        row=run_cut(direction,())
        if row['mapping'] is not None:proposals.append(dict(enumerate(row['mapping'])))
    for direction in ('R_to_P','P_to_R'):
        reverse=direction=='P_to_R';search_problem=AAMProblem(problem.product,problem.reactant) if reverse else problem
        r,p=search_problem.reactant.wbo,search_problem.product.wbo;cuts=set()
        for physical in proposals:
            m={b:a for a,b in physical.items()} if reverse else physical
            for a in range(len(r)):
                for b in range(a+1,len(r)):
                    if r[a,b]>.2 and (p[m[a],m[b]]<=.2 or abs(r[a,b]-p[m[a],m[b]])>.5):cuts.add((a,b))
        _initialize_search(search_problem,cfg)
        for edge in sorted(cuts):run_cut(direction,(edge,))
    elapsed=time.perf_counter()-start_wall;cpu=time.process_time()-start_cpu
    best=min([proposal['best']['events']['total']]+[r['best'] for r in results if r['best'] is not None])
    # Full-search answer is read only after the complete experimental workflow.
    baseline=json.loads((args.source/f'comparison/{args.index}.refined.json').read_text())
    if args.index==123:baseline=json.loads((args.source/'case123_cap200/comparison/123.refined.json').read_text())
    save(args.run/f'{args.index}.json',dict(index=args.index,best=best,baseline=baseline['aam_min_saved_events'],
         cpu_including_internal_artifact_io=cpu,wall_including_internal_artifact_io=elapsed,
         proposal=proposal,results=results,scope='Fresh one-CPU workflow, including proposals, H-refinement, AAM search, symmetry, WBO scoring and internal graph/proposal persistence; excludes initial imports/input load and final summary write.'))


def submit(args):
    args.run.mkdir(parents=True,exist_ok=False);(args.run/'proposal').mkdir();(args.run/'status').mkdir()
    (args.run/'inputs').symlink_to(args.source/'inputs',target_is_directory=True)
    shutil.copy2(args.pilot/'slap_core.py',args.run/'slap_core.py')
    for f in ('slap_local_workflow.py','slap_guided_pilot.py'):shutil.copy2(Path(__file__).with_name(f),args.run/f)
    cmd=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','PYTHONHASHSEED=0',
         f'PYTHONPATH={args.pilot}/dependencies:{args.pilot}/engine/src:{args.pilot}/engine/bench',
         'timeout','--kill-after=5s','300',sys.executable,str(args.run/'slap_local_workflow.py'),'case','--run',str(args.run),'--index']
    options=['sbatch','--parsable','--partition=cpunodes,cpunodes_nia','--exclude=bosque5,bosque6,bosque7,bosque8',
        '--cpus-per-task=1','--mem=8G','--time=00:10:00','--array=0-139%32','--job-name=slap_fresh',
        f'--output={args.run}/status/%A_%a.out','--wrap',shlex.join(cmd)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip();save(args.run/'submission.json',dict(job=job,command=options));print(job)


def relocate(args):
    previous=json.loads((args.run/'submission.json').read_text())
    tasks=subprocess.check_output(['squeue','-h','-j',previous['job'],'-t','CONFIGURING','-o','%i'],text=True).split()
    slots=[]
    for job in tasks:
        i=int(job.split('_')[1]);assert not (args.run/f'graphs/{i}').exists()
        assert not (args.run/f'proposal/{i}.refined.json').exists()
        subprocess.run(['scancel',job],check=True);slots.append(i)
    assert slots
    options=[('--array='+','.join(map(str,slots)) if s.startswith('--array=') else
              '--exclude=bosque5,bosque6,bosque7,bosque8,bosque10' if s.startswith('--exclude=') else s) for s in previous['command']]
    job=subprocess.check_output(options,text=True).strip();save(args.run/'relocation.json',dict(job=job,command=options,unstarted_slots=slots));print(job)


def report(args):
    from rxn_core.family_scoring import bond_events
    rows=[json.loads((args.run/f'{i}.json').read_text()) for i in range(140)]
    checked=0
    for row in rows:
        raw=json.loads((args.run/f"inputs/{row['index']}/input.json").read_text())
        problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
        for r in row['results']:
            if r['mapping'] is None:continue
            m=dict(enumerate(r['mapping']))
            assert sorted(m.values())==list(range(problem.atom_count))
            assert all(problem.reactant.elements[a]==problem.product.elements[b] for a,b in m.items())
            assert bond_events(problem,m)['total']==r['best'];checked+=1
        for r in row['proposal']['all']:
            assert bond_events(problem,dict(r['mapping']))['total']==r['events']['total'];checked+=1
    cpu=sum(r['cpu_including_internal_artifact_io'] for r in rows)
    elapsed=sum(r['wall_including_internal_artifact_io'] for r in rows)
    summary=dict(cases=140,equal=sum(r['best']==r['baseline'] for r in rows),
                 better=sum(r['best']<r['baseline'] for r in rows),worse=sum(r['best']>r['baseline'] for r in rows),
                 cpu_seconds=cpu,mean_cpu_seconds=cpu/140,mean_elapsed_seconds=elapsed/140,
                 max_elapsed_seconds=max(r['wall_including_internal_artifact_io'] for r in rows),
                 audited_witnesses=checked,workers_per_case=1,scope=rows[0]['scope'])
    save(args.run/'summary.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('case','submit','relocate','report'))
    p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path,default=SOURCE)
    p.add_argument('--pilot',type=Path,default=SOURCE.parent/'slap_guided_pilot_20260909');p.add_argument('--index',type=int)
    args=p.parse_args();globals()[args.command](args)

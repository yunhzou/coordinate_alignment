"""Saved-family diagnosis and controlled cap/direction ablations; no reference-guided search."""
import argparse
from collections import Counter
from dataclasses import replace,asdict
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pynauty
from golden_policy_campaign import load_case,save,guarded
from golden_evaluation import project,colored_graph,endpoint_generators,evaluate_planned
from rxn_core import AAMProblem,AAMSearchPlan,search_aam
from rxn_core.artifacts import read_aam_checkpoint


def initialize(args):
    args.run.mkdir(exist_ok=False,parents=True)
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    summary=json.loads(args.summary.read_text());jobs=[]
    for index in summary['not_recovered']+summary['unknown']:
        source=args.source/str(index)
        score=json.loads((source/'evaluation.json').read_text())
        jobs.append(dict(index=index,variant='archive'))
        if index in summary['not_recovered']:
            jobs.append(dict(index=index,variant='opposite'))
            if score['capped']:jobs.append(dict(index=index,variant='cap2000'))
    save(args.run/'manifest.json',dict(source=str(args.source.resolve()),jobs=jobs,
        revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        purpose='Separate diagnosis, not a replacement benchmark; 3 seeds, tolerance 1.0 unchanged',
        search_watchdog=300,score_watchdog=180,workers=2))
    print(len(jobs))


def context(args):
    manifest=json.loads((args.run/'manifest.json').read_text());job=manifest['jobs'][args.slot]
    source=Path(manifest['source']);directory,base=load_case(source,job['index'])
    plan=base
    if job['variant']=='cap2000':plan=replace(base,config=replace(base.config,branch_limit=2000))
    elif job['variant']=='opposite':
        problem=AAMProblem(base.problem.product,base.problem.reactant,base.problem.name)
        plan=AAMSearchPlan(base.input_problem,problem,
            replace(base.config,anchors=tuple((p,r) for r,p in base.config.anchors)),not base.reversed)
    output=args.run/f"{job['index']}_{job['variant']}"
    output.mkdir(exist_ok=True)
    return manifest,job,directory,plan,output


def search(args):
    manifest,job,source,plan,out=context(args)
    save(out/'orientation.json',dict(reversed=plan.reversed,direction=plan.direction,config=asdict(plan.config)))
    start=time.perf_counter()
    result=search_aam(plan.problem,plan.config,workers=manifest['workers'],
        intermediate_dir=out/'cuts',archive_format='checkpoint')
    save(out/'search.json',dict(seconds=time.perf_counter()-start,metrics=asdict(result.metrics)))


def diagnose(args):
    manifest,job,source,plan,out=context(args)
    archive=(source if job['variant']=='archive' else out)/'cuts/aam.pkl.gz'
    if not archive.exists():
        save(out/'diagnosis.json',dict(reference_recovery='unknown',reason='incomplete search; no full archive',
            saved_cuts=len(list(archive.parent.glob('cut_*.json')))))
        return
    start=time.perf_counter();aam=read_aam_checkpoint(archive)
    ref=json.loads((source/'reference.json').read_text())
    features=list(reversed(ref['features'])) if plan.reversed else ref['features']
    expected=project(plan.to_search_mapping(ref['mapping']),features)
    orbits=pynauty.autgrp(colored_graph([features[0]]))[3][:len(features[0]['heavy'])]
    def source_key(mapping):return tuple(sorted(Counter(orbits[r] for r in mapping).items()))
    required=source_key(expected)
    eligible={t for t in aam.graph.terminals if source_key(project(aam.graph.states[t].mapping,features))==required}
    source_eligible=len(eligible)
    # A union of target action orbits is an overapproximation: reject only,
    # never use it to assert a valid correlated mapping.
    degree=len(features[1]['heavy']);parents=list(range(degree))
    def root(a):
        while parents[a]!=a:parents[a]=parents[parents[a]];a=parents[a]
        return a
    def join(g):
        for a,b in enumerate(g):parents[root(a)]=root(b)
    for g in endpoint_generators(features[1]):join(g)
    target_index={a:i for i,a in enumerate(features[1]['heavy'])}
    ancestors=aam.graph.ancestor_transitions(eligible)
    for edge in aam.graph.transitions:
        if edge.id not in ancestors or edge.match is None:continue
        placement=aam.graph.fragment_placement(edge.id)
        for g in placement.target_generators:join(tuple(target_index[g.images[a]] for a in features[1]['heavy']))
        # Assignment domains can cross conditioned automorphism orbits.
        # Include them in this rejection-only overapproximation as well.
        for domain in placement.symmetry_domains+placement.automorph_domains:
            atoms=[target_index[p] for p in domain.p_atoms if p in target_index]
            for atom in atoms[1:]:parents[root(atom)]=root(atoms[0])
    target_orbits=[root(i) for i in range(degree)]
    def pairs(mapping):return Counter((orbits[r],target_orbits[p]) for r,p in mapping.items())
    wanted=pairs(expected)
    eligible={t for t in eligible if pairs(project(aam.graph.states[t].mapping,features))==wanted}
    info=dict(index=job['index'],variant=job['variant'],direction=plan.direction,archive=str(archive),
        source_heavy=len(features[0]['heavy']),target_heavy=len(features[1]['heavy']),
        capped=aam.graph.capped,terminals=len(aam.graph.terminals),source_eligible=source_eligible,
        joint_orbit_eligible=len(eligible),precheck_seconds=time.perf_counter()-start,
        scope='Saved compressed families; reference used only for evaluation; top-1 not recomputed')
    save(out/'filter.json',info)
    if not eligible:
        info.update(reference_recovery='not_recovered',reason='source occupation impossible' if not source_eligible else 'joint endpoint orbit occupation impossible')
    else:
        view=replace(aam,graph=replace(aam.graph,stops=tuple(s for s in aam.graph.stops
            if s.reason not in {'objective_met','stalled'} or s.state in eligible)))
        score=evaluate_planned(view,plan,ref['features'],ref['mapping'],seconds=130,
            query_timeout_ms=15000)
        save(out/'verification.json',score)
        info.update(reference_recovery='recovered' if score['reference_recovery']=='recovered' else 'unknown',
            generator_only_outcome=score['reference_recovery'],symbolic_queries=score['symbolic_queries'],
            unknown_queries=score['unknown_queries'],reason='compressed-family verification after necessary filters')
    info['seconds']=time.perf_counter()-start
    save(out/'diagnosis.json',info)


def worker(args):
    manifest,job,source,plan,out=context(args)
    def phase(name,limit):
        save(out/'status.json',dict(stage=name,job=os.environ.get('SLURM_JOB_ID'),updated=time.time()))
        code=guarded([sys.executable,__file__,name,'--run',str(args.run),'--slot',str(args.slot)],out/f'{name}.log',limit)
        save(out/f'{name}_status.json',dict(exit_code=code,ended=time.time()))
        return code
    if job['variant']!='archive':phase('search',manifest['search_watchdog'])
    code=phase('diagnose',manifest['score_watchdog'])
    save(out/'status.json',dict(stage='complete' if code==0 else 'evaluation_interrupted',updated=time.time()))


def check_partial(args):
    from rxn_core.search_graph import AAMSearchGraph
    manifest,job,source,plan,out=context(args)
    ref=json.loads((source/'reference.json').read_text())
    features=list(reversed(ref['features'])) if plan.reversed else ref['features']
    expected=pynauty.certificate(colored_graph(features,project(plan.to_search_mapping(ref['mapping']),features)))
    start=time.perf_counter();checked=0
    for cut in sorted((out/'cuts').glob('cut_*.json')):
        if time.perf_counter()-start>90:break
        graph=AAMSearchGraph.from_record(json.loads(cut.read_bytes()),copy=False);seen=set()
        for terminal in graph.terminals:
            if time.perf_counter()-start>90:break
            mapping=project(graph.states[terminal].mapping,features);key=tuple(sorted(mapping.items()))
            if key in seen:continue
            seen.add(key);checked+=1
            if pynauty.certificate(colored_graph(features,mapping))==expected:
                save(out/'partial_witness.json',dict(reference_recovery='recovered',cut=str(cut),terminal=terminal,
                    input_orientation_witness=sorted(plan.to_input_mapping(graph.states[terminal].mapping).items()),
                    scope='Positive certificate from a completed cut; full sweep incomplete',seconds=time.perf_counter()-start))
                return
    save(out/'partial_witness.json',dict(reference_recovery='unknown',checked=checked,
        scope='Partial representative-only check; not evidence of absence',seconds=time.perf_counter()-start))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['init','worker','search','diagnose','partial'])
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--source',type=Path)
    p.add_argument('--summary',type=Path)
    p.add_argument('--slot',type=int)
    args=p.parse_args()
    dict(init=initialize,worker=worker,search=search,diagnose=diagnose,partial=check_partial)[args.mode](args)

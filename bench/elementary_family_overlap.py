"""Certified positive family overlap and bounded exclusion checks on saved AAM."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np
import pynauty
from compare_elementary_outputs import features,certificate,event_counts,events_row
from golden_evaluation import colored_graph
from golden_competitors import save
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.family_query import query_path


def analyze(args):
    start=time.perf_counter();deadline=start+180
    comparison=json.loads((args.comparison/f'{args.index}.refined.json').read_text())
    root=Path(comparison['aam_source'])
    raw=json.loads((root/'inputs'/str(args.index)/'input.json').read_text())
    problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    feat=features(raw);heavy=[i for i,e in enumerate(raw['reactant']['elements']) if e!='H']
    groups=[pynauty.autgrp(colored_graph([f])) for f in feat]
    orbits=[g[3][:len(f['colors'])] for g,f in zip(groups,feat)]
    generators=[tuple(tuple(g[:len(f['colors'])]) for g in group[0]) for group,f in zip(groups,feat)]
    def orbit_key(mapping):
        return tuple(sorted(Counter((orbits[0][i],orbits[1][mapping[i]]) for i in heavy).items()))
    targets={};target_orbits=set()
    for candidate in comparison['slap']:
        mapping=dict(candidate['mapping']);vector=[mapping[i] for i in range(len(mapping))]
        key=certificate(feat,vector,heavy)
        targets.setdefault(key,dict(mapping={i:mapping[i] for i in heavy},status='pending',candidates=[]))['candidates'].append(candidate['candidate'])
        target_orbits.add(orbit_key(vector))
    cache={};extra=None;scanned=0
    for direction in ('R_to_P','P_to_R'):
        with gzip.open(root/'directions'/str(args.index)/direction/'terminal_mappings.jsonl.gz','rt') as stream:
            while True:
                rows=[];vectors=[]
                for _ in range(2048):
                    line=stream.readline()
                    if not line:break
                    row=json.loads(line)
                    if not row['full_element_bijection']:continue
                    m=dict(row['mapping']);vectors.append([m[i] for i in range(len(m))]);rows.append(row)
                if not vectors:break
                scores=event_counts(problem.reactant.wbo,problem.product.wbo,vectors)
                for row,vector,score in zip(rows,vectors,scores):
                    scanned+=1;h=tuple(vector[i] for i in heavy)
                    if h not in cache:
                        cache[h]=(certificate(feat,vector,heavy) if orbit_key(vector) in target_orbits else None)
                    key=cache[h]
                    if key in targets:
                        if targets[key]['status']=='pending':targets[key].update(status='represented',
                            method='saved_terminal_certificate',direction=direction,terminal=row['terminal'],mapping_witness=row['mapping'])
                    elif extra is None or int(sum(score))<extra['events']['total']:
                        extra=dict(direction=direction,terminal=row['terminal'],mapping=row['mapping'],events=events_row(score))
    invariant_cache={};symbolic_queries=0;paths_checked=0;complete=True;unknown=set()
    for direction in ('R_to_P','P_to_R'):
        if all(v['status']=='represented' for v in targets.values()):break
        if time.perf_counter()>deadline:complete=False;break
        reverse=direction=='P_to_R';ordered=feat[::-1] if reverse else feat
        search_problem=AAMProblem(problem.product,problem.reactant) if reverse else problem
        hs=[i for i,e in enumerate(search_problem.reactant.elements) if e!='H']
        ht=[i for i,e in enumerate(search_problem.product.elements) if e!='H']
        source_g,target_g=generators[::-1] if reverse else generators
        edge_labels={tuple(sorted((a,b))):label for a,b,label in ordered[1]['bonds']}
        colors=ordered[1]['colors']
        def preserves(g):
            k=(direction,tuple(g))
            if k not in invariant_cache:
                invariant_cache[k]=(all(g[i]==i for i in ht) or
                    all(colors[i]==colors[g[i]] for i in range(len(colors))) and
                    all(edge_labels.get(tuple(sorted((g[a],g[b]))))==label for (a,b),label in edge_labels.items()))
            return invariant_cache[k]
        def invariant(path):
            for edge in path.transitions:
                placement=path.graph.fragment_placement(edge)
                if placement is None:continue
                assert placement.target_generators is not None
                if any(not preserves(g.images) for g in placement.target_generators):return False
                for domain in placement.symmetry_domains:
                    if domain.source=='exact_automorph_group':continue
                    for a,b in zip(domain.p_atoms,domain.p_atoms[1:]):
                        g=list(range(len(colors)));g[a],g[b]=b,a
                        if not preserves(g):return False
            return True
        graph=read_aam_checkpoint(root/'directions'/str(args.index)/direction/'cuts/aam.pkl.gz').graph
        for path in graph.paths():
            if time.perf_counter()>deadline:complete=False;break
            paths_checked+=1
            if invariant(path):continue  # Entire path has its already-checked representative class.
            for key,target in targets.items():
                if target['status']=='represented':continue
                reference={b:a for a,b in target['mapping'].items()} if reverse else target['mapping']
                status,witness=query_path(path,search_problem,reference,source_atoms=hs,
                    source_generators=source_g,target_generators=target_g,
                    timeout_ms=max(1,min(2000,int((deadline-time.perf_counter())*1000))))
                symbolic_queries+=1
                if status=='recovered':
                    physical=dict(witness['mapping'])
                    if reverse:physical={b:a for a,b in physical.items()}
                    assert certificate(feat,[physical[i] for i in range(len(physical))],heavy)==key
                    target.update(status='represented',method='compressed_family_query',direction=direction,
                        terminal=path.terminal,mapping_witness=sorted(physical.items()),query_witness=witness)
                elif status=='unknown':unknown.add(key)
            if all(v['status']=='represented' for v in targets.values()):break
        if not complete:break
    for key,target in targets.items():
        if target['status']=='pending':target['status']='excluded_from_saved_families' if complete and key not in unknown else 'unresolved'
    if extra:extra['above_aam_saved_min']=extra['events']['total']-comparison['aam_min_saved_events']
    result=dict(index=args.index,name=raw['name'],targets=list(targets.values()),
        targets_unique_classes=len(targets),statuses=dict(Counter(v['status'] for v in targets.values())),
        extra_aam_pattern_absent_from_saved_slap=extra,scanned_witnesses=scanned,
        symbolic_queries=symbolic_queries,paths_checked=paths_checked,seconds=time.perf_counter()-start,
        scope='Heavy mapping families modulo shared score-preserving endpoint graph symmetry. Explicit H retained in feasibility. '
              'SLAP heavy labels are singleton in these saved outputs. Exclusion concerns saved families only, never all possible searches. '
              f'No chemical correctness labels; search source: {root}.')
    args.run.mkdir(exist_ok=True,parents=True);save(args.run/f'{args.index}.json',result)


def submit(args):
    args.run.mkdir(exist_ok=False,parents=True)
    shutil.copy2(__file__,args.run/'driver.py')
    command=['env','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1',f'PYTHONPATH={Path(__file__).resolve().parent}',
        'timeout','--kill-after=5s','300',sys.executable,str(args.run/'driver.py'),'analyze',
        '--comparison',str(args.comparison),'--run',str(args.run),'--index']
    options=['sbatch','--parsable','--partition=cpunodes','--exclude=bosque8','--cpus-per-task=1',
        '--mem=8G','--time=00:10:00','--array=0-139%32','--job-name=elem_overlap',
        f'--output={args.run}/slurm_%A_%a.out','--wrap',shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=options,sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print(job)


def report(args):
    rows=[json.loads((args.run/f'{i}.json').read_text()) for i in range(140)]
    totals=Counter(s for r in rows for s in [t['status'] for t in r['targets']])
    save(args.run/'summary.json',dict(cases=140,target_class_statuses=dict(totals),
        all_slap_classes_represented_cases=sum(all(t['status']=='represented' for t in r['targets']) for r in rows),
        cases_with_excluded_slap_class=[r['index'] for r in rows if any(t['status']=='excluded_from_saved_families' for t in r['targets'])],
        cases_with_unresolved_slap_class=[r['index'] for r in rows if any(t['status']=='unresolved' for t in r['targets'])],
        aam_extra_cases_by_event_window={str(w):sum(r['extra_aam_pattern_absent_from_saved_slap'] is not None and
            r['extra_aam_pattern_absent_from_saved_slap']['above_aam_saved_min']<=w for r in rows) for w in (0,1,2,5)},
        symbolic_queries=sum(r['symbolic_queries'] for r in rows),seconds_sum=sum(r['seconds'] for r in rows)))
    print((args.run/'summary.json').read_text())


def relocate(args):
    previous=json.loads((args.run/'submission.json').read_text())
    tasks=subprocess.check_output(['squeue','-h','-j',previous['job'],'-t','CONFIGURING','-o','%i'],text=True).split()
    slots=[]
    for job in tasks:
        index=int(job.split('_')[1]);assert not (args.run/f'{index}.json').exists()
        subprocess.run(['scancel',job],check=True);slots.append(index)
    if not slots:return
    command=[('--exclude=bosque5,bosque6,bosque8' if x.startswith('--exclude=') else
              '--array='+','.join(map(str,slots))+'%32' if x.startswith('--array=') else x) for x in previous['command']]
    job=subprocess.check_output(command,text=True).strip()
    save(args.run/'relocation.json',dict(original_job=previous['job'],unstarted_slots=slots,job=job,command=command));print(job)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('analyze','submit','report','relocate'))
    p.add_argument('--run',type=Path,required=True);p.add_argument('--index',type=int)
    p.add_argument('--comparison',type=Path,default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_output_comparison_20260908'))
    args=p.parse_args();globals()[args.command](args)

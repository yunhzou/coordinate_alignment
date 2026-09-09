"""Blind mapping of cached real endpoint/TS pairs, with separate saved evaluation.

No SMILES, synthetic oracle, quantum calculation, or production-core changes.
Every AAM cut graph and every native SLAP result (including LAP internals) is saved.
"""
import argparse
from collections import defaultdict
from dataclasses import asdict
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np

METHODS = ('aam_weighted', 'aam_binary', 'slap_binary', 'slap_weighted', 'slap_float_cost')


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2)+'\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_slap(run, method):
    name = 'slap_core_experiment'
    spec = importlib.util.spec_from_file_location(name, run/'slap_core.py')
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    if method == 'slap_float_cost':
        # Explicit numerical diagnostic, not an unmodified upstream result.
        # Keep upstream search/perturbation/tie rules; replace only the matrix
        # returned to the LAP solver with the pre-cast floating neighborhood cost.
        def floating_cost(self, graphs, label):
            _, masks, infos = original(self, graphs, label)
            left, right = [list(info.values()) for info in infos]
            matrix = np.block([
                [np.full((len(a['idxs']), len(b['idxs'])),
                         self._diff_nbrs(a['nbrs'], b['nbrs'])) for b in right]
                for a in left
            ])
            return matrix, masks, infos
        original = module.SlapMapper._get_cost_matrix
        module.SlapMapper._get_cost_matrix = floating_cost
    return module


def prepare(args):
    from rxn_core.chemistry_computations.xtb import load_cached_xtb
    args.run.mkdir(parents=True, exist_ok=False)
    cases = ('pr1.tempo_ts1', 'pr16.carbocation_ts5', 'pr7.V.dodh_ts910')
    inputs, references = [], []
    for name in cases:
        work = args.source/'work'/name
        stage_path = args.source/'stages'/name/'ts_stage.json'
        stage = json.loads(stage_path.read_text())
        for side in ('R', 'P'):
            endpoints, files = [], {}
            for label, folder in ((side, work/'endpoints'/side), ('TS', work/'targets/GT_sp')):
                elements, xyz, wbo, path = load_cached_xtb(folder)
                endpoints.append(dict(elements=list(elements), coordinates=xyz.tolist(),
                                      wbo=wbo.tolist(), label=label))
                files[str(path)] = sha(path)
                files[str(folder/'wbo')] = sha(folder/'wbo')
            index = len(inputs)
            data = dict(index=index, name=name, direction=side+'_to_TS',
                        reactant=endpoints[0], product=endpoints[1], input_sha256=files)
            save(args.run/f'inputs/{index}.json', data)
            refs = []
            for m in stage['mechanisms']:
                if m['gt'] is None:
                    refs.append(dict(mechanism=m['id'], status='no_historical_TS_assignment'))
                    continue
                mapping = {int(a): int(b) for a,b in m['gt']['core_map'].items()}
                if side == 'P':
                    mapping = {int(m['mapping_RP'][str(a)]): b for a,b in mapping.items()}
                refs.append(dict(mechanism=m['id'], status='available',
                    heavy_core={a:b for a,b in mapping.items() if endpoints[0]['elements'][a]!='H'},
                    full_core=mapping))
            inputs.append(dict(index=index, name=name, direction=side+'_to_TS'))
            references.append(dict(index=index, source=str(stage_path), sha256=sha(stage_path),
                                   references=refs))
    save(args.run/'references.json', references)
    tasks = [dict(slot=i*len(METHODS)+j, index=i, method=m) for i in range(len(inputs)) for j,m in enumerate(METHODS)]
    save(args.run/'tasks.json', tasks)
    root = Path(__file__).resolve().parents[1]
    for folder in ('src', 'bench'):
        shutil.copytree(root/folder, args.run/'engine'/folder,
                        ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(args.slap_core, args.run/'slap_core.py')
    save(args.run/'manifest.json', dict(inputs=inputs, methods=METHODS, tasks=len(tasks),
        seed_count=10, root_seed=42, input_order_seed=20260909, branch_cap=100,
        iso_tolerance=1.0, graph_floor=.2, scoring_event_tolerance=.5,
        explicit_hydrogen=True, aam_cuts='uncut plus every source edge, both AAM weight modes',
        slap_orders=10, slap_break_sym='heavy',
        float_variant='Only LAP input matrix dtype corrected; all other upstream heuristics unchanged.',
        reference_role='Historical algorithm-generated TS core assignments; evaluation only, not curated truth.',
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        driver_sha256=sha(Path(__file__)), slap_sha256=sha(args.slap_core),
        native_sha256=sha(next((root/'src/rxn_core').glob('_engine*.so')))))
    print(args.run)


def submit(args):
    (args.run/'status').mkdir()
    driver=args.run/'engine/bench/compare_real_ts_mappings.py'
    env=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',
         'PYTHONHASHSEED=0','RXN_CORE_NATIVE=1',
         f'PYTHONPATH={args.dependencies}:{args.run}/engine/src:{args.run}/engine/bench']
    cmd=env+['timeout','--kill-after=5s','300',sys.executable,str(driver),'task','--run',str(args.run),'--slot']
    options=['sbatch','--parsable','--partition=cpunodes','--nodelist=bosque7',
        '--cpus-per-task=1','--mem=6G','--time=00:10:00','--array=0-29%30',
        '--job-name=ts_compare',f'--output={args.run}/status/%A_%a.out',
        '--wrap',shlex.join(cmd)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json', dict(job=job, command=options))
    print(job)


def task(args):
    from rxn_core import AAMProblem, AAMSearchConfig
    from rxn_core.domain import MolecularEndpoint
    from rxn_core.frag import build_graph
    from rxn_core.aam import _initialize_search, _search_cut
    from rxn_core.alignment.sweep import cut_sweep_items
    from rxn_core.search_symmetry import finalize_graph_symmetry
    from rxn_core.artifacts import write_graph_checkpoint
    from rdkit import Chem
    spec=json.loads((args.run/'tasks.json').read_text())[args.slot]
    raw=json.loads((args.run/f"inputs/{spec['index']}.json").read_text())
    folder=args.run/f"results/{spec['index']}/{spec['method']}"
    folder.mkdir(parents=True, exist_ok=False)
    method=spec['method']; n=len(raw['reactant']['elements'])
    problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
    phases=defaultdict(lambda: dict(cpu=0.,wall=0.))
    def timed(label, fn):
        cpu=time.process_time(); start=time.perf_counter()
        value=fn()
        phases[label]['cpu']+=time.process_time()-cpu
        phases[label]['wall']+=time.perf_counter()-start
        return value
    start=time.perf_counter(); vectors={}; progress=[]
    if method.startswith('aam'):
        data=[dict(raw[k]) for k in ('reactant','product')]
        if method=='aam_binary':
            for d in data:d['wbo']=(np.asarray(d['wbo'])>=.2).astype(float)
        search_problem=AAMProblem(*(MolecularEndpoint(**d) for d in data))
        cfg=AAMSearchConfig(seed_count=10,branch_limit=100,iso_tolerance=1.)
        timed('prepare',lambda:_initialize_search(search_problem,cfg))
        target=build_graph(search_problem.product.elements, search_problem.product.wbo,bond_cut=.2)
        # Freeze identical source edge sets and cut identifiers across weight modes.
        cuts=cut_sweep_items(problem.reactant.wbo,.2)
        for ordinal,cut in enumerate(cuts):
            graph,_=timed('search',lambda:_search_cut(cut))
            graph,_=timed('symmetry',lambda:finalize_graph_symmetry(graph,target,iso_tolerance=1.))
            timed('persistence',lambda:write_graph_checkpoint(graph,folder/f'cut_{ordinal:04d}.pkl.gz'))
            def collect():
                for terminal in graph.terminals:
                    mapping=dict(graph.states[terminal].mapping)
                    if len(mapping)!=n:continue
                    key=tuple(mapping[i] for i in range(n))
                    vectors.setdefault(key,dict(cut=ordinal,terminal=terminal))
            timed('collect',collect)
            progress.append(dict(cut=cut,terminals=len(graph.terminals),capped=graph.capped,
                                 states=len(graph.states)))
            save(folder/'progress.json',dict(completed=len(progress),total=len(cuts),phases=dict(phases)))
    else:
        slap=load_slap(args.run,method)
        rng=np.random.default_rng(20260909+spec['index'])
        for order in range(10):
            ro=np.arange(n) if order==0 else rng.permutation(n)
            po=np.arange(n) if order==0 else rng.permutation(n)
            graphs=[]
            for endpoint,permutation in zip((problem.reactant,problem.product),(ro,po)):
                w=endpoint.wbo[np.ix_(permutation,permutation)]
                graphs.append(slap.LabeledGraph(
                    {i:{j:(1.0 if method=='slap_binary' else float(w[i,j]))
                        for j in range(n) if i!=j and w[i,j]>=.2} for i in range(n)},
                    [Chem.GetPeriodicTable().GetAtomicNumber(endpoint.elements[i]) for i in permutation]))
            mapper=slap.SlapMapper(binary=method=='slap_binary')
            timed('search',lambda:mapper.get_maps(graphs,
                break_sym_targets=[i for i,a in enumerate(ro) if problem.reactant.elements[a]!='H']))
            def persist():
                with gzip.open(folder/f'order_{order:02d}.pkl.gz','wb',compresslevel=1) as stream:
                    pickle.dump(dict(results=mapper.results,source_order=ro,target_order=po),stream)
            timed('persistence',persist)
            def collect():
                for ordinal,result in enumerate(mapper.results):
                    left,right=result['lgp']
                    assert all(len(atoms)==1 or all(problem.reactant.elements[ro[a]]=='H' for a in atoms)
                               for atoms in left.label2idxs.values())
                    mapping={int(ro[a]):int(po[b]) for label,atoms in left.label2idxs.items()
                             for a,b in zip(atoms,right.label2idxs[label])}
                    assert len(mapping)==len(set(mapping.values()))==n
                    assert all(problem.reactant.elements[a]==problem.product.elements[b] for a,b in mapping.items())
                    vectors.setdefault(tuple(mapping[i] for i in range(n)),dict(order=order,candidate=ordinal,
                                                                                native_cost=float(result['val'])))
            timed('collect',collect)
            progress.append(dict(order=order,results=len(mapper.results)))
            save(folder/'progress.json',dict(completed=len(progress),total=10,phases=dict(phases)))
    timed('persistence',lambda:save(folder/'witnesses.json',[
        dict(mapping=list(k),origin=v) for k,v in sorted(vectors.items())]))
    save(folder/'summary.json',dict(**spec,name=raw['name'],direction=raw['direction'],
        phases=dict(phases),units=progress,full_witnesses=len(vectors),host=os.uname().nodename,
        elapsed_including_output=time.perf_counter()-start,
        timing_scope='One CPU. Search and symmetry separate from collection/persistence; import/input loading excluded.',
        witness_scope='AAM terminal representatives, not exhaustive family expansion; SLAP heavy labels singleton.'))


def analyze(args):
    from compare_elementary_outputs import features, certificate, event_counts, events_row
    references=json.loads((args.run/'references.json').read_text())
    rows=[]
    for item in json.loads((args.run/'manifest.json').read_text())['inputs']:
        index=item['index']; raw=json.loads((args.run/f'inputs/{index}.json').read_text())
        r,p=[np.asarray(raw[k]['wbo']) for k in ('reactant','product')]
        heavy=[i for i,e in enumerate(raw['reactant']['elements']) if e!='H']
        feat=features(raw); a,b=np.triu_indices(len(r),1)
        ha,hb=np.array([(x,y) for x,y in zip(a,b) if x in heavy and y in heavy]).T
        entry=dict(**item, methods={})
        for method in METHODS:
            folder=args.run/f'results/{index}/{method}'
            if not (folder/'summary.json').exists():
                entry['methods'][method]=dict(status='incomplete');continue
            summary=json.loads((folder/'summary.json').read_text())
            candidates=json.loads((folder/'witnesses.json').read_text())
            classes={}; matches=defaultdict(list); best=None
            for c in candidates:
                m=np.array(c['mapping']); key=hashlib.sha256(certificate(feat,m,heavy)).hexdigest()
                events=events_row(event_counts(r,p,[m])[0])
                rw,pw=r[ha,hb],p[m[ha],m[hb]]
                he=int(np.sum(((rw>.2)!=(pw>.2))|((rw>.2)&(pw>.2)&(np.abs(rw-pw)>.5))))
                l1=float(np.sum(np.abs(rw-pw)))
                record=dict(**c,events=events,heavy_events=he,heavy_wbo_l1=l1,certificate=key)
                classes.setdefault(key,record)
                if best is None or (he,l1)<(best['heavy_events'],best['heavy_wbo_l1']):best=record
                for ref in references[index]['references']:
                    if ref['status']!='available':continue
                    core={int(x):int(y) for x,y in ref['heavy_core'].items()}
                    if all(m[x]==y for x,y in core.items()):matches[str(ref['mechanism'])].append(record)
            phase=summary['phases']
            result=dict(status='complete',full_witnesses=len(candidates),heavy_classes=len(classes),
                class_scope='Score-preserving endpoint-symmetry classes; AAM witness lower bound, not full-family count.',
                compute_cpu=sum(v['cpu'] for k,v in phase.items() if k!='persistence'),
                phases=phase, capped_units=sum(x.get('capped',False) for x in summary['units']),
                total_units=len(summary['units']),best=best,
                historical_core={str(ref['mechanism']):dict(
                    status=ref['status'] if ref['status']!='available' else
                    ('witness_recovered' if matches[str(ref['mechanism'])] else 'not_in_witnesses'),
                    example=matches[str(ref['mechanism'])][0] if matches[str(ref['mechanism'])] else None)
                    for ref in references[index]['references']})
            save(folder/'analysis.json',result)
            save(folder/'heavy_classes.json',classes)
            entry['methods'][method]=result
        if all(v['status']=='complete' for v in entry['methods'].values()):
            keys={m:set(json.loads((args.run/f'results/{index}/{m}/heavy_classes.json').read_text())) for m in METHODS}
            entry['overlap']={a:{b:dict(shared=len(keys[a]&keys[b]),only_left=len(keys[a]-keys[b]))
                                   for b in METHODS} for a in METHODS}
        rows.append(entry)
    save(args.run/'analysis.json',dict(cases=rows,
        scope='Actual R/TS and P/TS mappings. Common heavy-edge score and saved historical core recovery, not chemical accuracy.'))
    for row in rows:
        print(row['name'],row['direction'])
        for m,v in row['methods'].items():
            print(m, v['status'], 'classes',v.get('heavy_classes'), 'best_heavy_events',
                  v.get('best',{}).get('heavy_events') if v.get('best') else None,
                  'cpu',round(v.get('compute_cpu',0),3),'refs',
                  {k:w['status'] for k,w in v.get('historical_core',{}).items()})


def relocate(args):
    previous=json.loads((args.run/'submission.json').read_text())
    tasks=json.loads((args.run/'tasks.json').read_text())
    jobs=subprocess.check_output(['squeue','-h','-j',previous['job'],'-t','CONFIGURING','-o','%i'],text=True).split()
    slots=[]
    for job in jobs:
        slot=int(job.split('_')[1]); task_spec=tasks[slot]
        assert not (args.run/f"results/{task_spec['index']}/{task_spec['method']}").exists()
        subprocess.run(['scancel',job],check=True); slots.append(slot)
    if not slots:return
    command=[('--nodelist=bosque11' if v.startswith('--nodelist=') else
              '--array='+','.join(map(str,slots)) if v.startswith('--array=') else v)
             for v in previous['command']]
    job=subprocess.check_output(command,text=True).strip()
    save(args.run/'relocation.json',dict(original_job=previous['job'],unstarted_slots=slots,
                                       job=job,command=command))
    print(job)


def membership(args):
    from rxn_core import AAMProblem
    from rxn_core.domain import MolecularEndpoint
    from rxn_core.artifacts import read_graph_checkpoint
    from rxn_core.family_query import query_path
    from rxn_core.search_graph import frozen_value
    refs=json.loads((args.run/'references.json').read_text())
    reports=json.loads((args.run/'analysis.json').read_text())
    results=[]
    for row in reports['cases']:
        raw=json.loads((args.run/f"inputs/{row['index']}.json").read_text())
        for method in ('aam_weighted','aam_binary'):
            record=row['methods'][method]
            if record['status']!='complete':continue
            data=[dict(raw[k]) for k in ('reactant','product')]
            if method=='aam_binary':
                for d in data:d['wbo']=(np.asarray(d['wbo'])>=.2).astype(float)
            problem=AAMProblem(*(MolecularEndpoint(**d) for d in data))
            folder=args.run/f"results/{row['index']}/{method}"
            for reference in refs[row['index']]['references']:
                mid=str(reference['mechanism'])
                if reference['status']!='available':continue
                if record['historical_core'][mid]['status']=='witness_recovered':
                    results.append(dict(index=row['index'],method=method,mechanism=mid,
                                        status='recovered',via='terminal_witness'))
                    continue
                target={int(a):int(b) for a,b in reference['heavy_core'].items()}
                start=time.perf_counter(); seen=set(); unknown=False; found=None; scanned=0
                for filename in sorted(folder.glob('cut_*.pkl.gz')):
                    if time.perf_counter()-start>60:unknown=True;break
                    graph=read_graph_checkpoint(filename)
                    for path in graph.paths():
                        if time.perf_counter()-start>60:unknown=True;break
                        if len(path.mapping)!=problem.source_atom_count:continue
                        if not set(target)<=set(path.mapping):continue
                        key=(tuple(sorted(path.mapping.items())),path.context.cuts,
                             tuple(frozen_value(f) for f in path.fragments))
                        if key in seen:continue
                        seen.add(key);scanned+=1
                        status,witness=query_path(path,problem,target,source_atoms=tuple(target),
                            complete_reference=False,timeout_ms=500)
                        if status=='recovered':
                            found=dict(graph=str(filename),terminal=path.terminal,
                                       transitions=path.transitions,witness=witness);break
                        unknown|=status=='unknown'
                    if found:break
                result=dict(index=row['index'],method=method,mechanism=mid,
                    status='recovered' if found else 'unresolved' if unknown else 'excluded_from_saved_families',
                    via='compressed_family_query',paths_queried=scanned,seconds=time.perf_counter()-start,
                    example=found,reference=target)
                results.append(result)
                save(args.run/'membership.json',results)
    save(args.run/'membership.json',results)
    print(json.dumps([{k:v for k,v in r.items() if k not in ('example','reference')} for r in results],indent=2))


def core_equivalence(args):
    """Check historical cores modulo common score-preserving endpoint symmetry."""
    import pynauty
    import z3
    from compare_elementary_outputs import features
    from golden_evaluation import colored_graph
    from rxn_core.family_query import SymbolicActions, query_path
    from rxn_core import AAMProblem
    from rxn_core.domain import MolecularEndpoint
    from rxn_core.artifacts import read_graph_checkpoint
    from rxn_core.search_graph import frozen_value
    reports=json.loads((args.run/'analysis.json').read_text())
    refs=json.loads((args.run/'references.json').read_text()); results=[]
    for row in reports['cases']:
        raw=json.loads((args.run/f"inputs/{row['index']}.json").read_text())
        n=len(raw['reactant']['elements'])
        generators=[tuple(tuple(g[:n]) for g in pynauty.autgrp(colored_graph([f]))[0])
                    for f in features(raw)]
        for method,analysis in row['methods'].items():
            folder=args.run/f"results/{row['index']}/{method}"
            candidates=json.loads((folder/'witnesses.json').read_text())
            for reference in refs[row['index']]['references']:
                if reference['status']!='available':continue
                mid=str(reference['mechanism']); target={int(a):int(b) for a,b in reference['heavy_core'].items()}
                result=dict(index=row['index'],method=method,mechanism=mid,reference=target)
                if analysis['historical_core'][mid]['status']=='witness_recovered':
                    result.update(status='recovered',via='literal_terminal_witness')
                    results.append(result);continue
                start=time.perf_counter(); unknown=False; found=None; seen=set()
                for c in candidates:
                    mapping=c['mapping']; key=tuple(mapping[i] for i,e in enumerate(raw['reactant']['elements']) if e!='H')
                    if key in seen:continue
                    seen.add(key)
                    solver=z3.Solver(); encoder=SymbolicActions(solver)
                    values=encoder.act([encoder.constant(a) for a in target],generators[0],'source')
                    values=[encoder.lookup(v,{i:encoder.constant(p) for i,p in enumerate(mapping)}) for v in values]
                    values=encoder.act(values,generators[1],'target')
                    solver.add(*(v[0]==b for v,b in zip(values,target.values())))
                    solver.set(timeout=500); status=solver.check()
                    if status==z3.sat:found=c;break
                    unknown|=status==z3.unknown
                if not found and method.startswith('aam'):
                    data=[dict(raw[k]) for k in ('reactant','product')]
                    if method=='aam_binary':
                        for d in data:d['wbo']=(np.asarray(d['wbo'])>=.2).astype(float)
                    problem=AAMProblem(*(MolecularEndpoint(**d) for d in data)); seen=set()
                    for filename in sorted(folder.glob('cut_*.pkl.gz')):
                        if time.perf_counter()-start>60:unknown=True;break
                        graph=read_graph_checkpoint(filename)
                        for path in graph.paths():
                            if time.perf_counter()-start>60:unknown=True;break
                            if len(path.mapping)!=n:continue
                            key=(tuple(sorted(path.mapping.items())),path.context.cuts,
                                 tuple(frozen_value(f) for f in path.fragments))
                            if key in seen:continue
                            seen.add(key)
                            status,witness=query_path(path,problem,target,source_atoms=tuple(target),
                                source_generators=generators[0],target_generators=generators[1],
                                complete_reference=False,timeout_ms=500)
                            if status=='recovered':
                                found=dict(graph=str(filename),terminal=path.terminal,witness=witness);break
                            unknown|=status=='unknown'
                        if found:break
                result.update(status='recovered' if found else 'unresolved' if unknown else 'excluded_from_saved_output',
                              via='symmetry_normalized_query',seconds=time.perf_counter()-start,example=found)
                results.append(result)
                save(args.run/'core_equivalence.json',results)
    save(args.run/'core_equivalence.json',results)
    print(json.dumps([{k:v for k,v in r.items() if k not in ('example','reference')} for r in results],indent=2))


def submit_scores(args):
    driver=args.run/'score_real_ts_families.py'
    shutil.copy2(Path(__file__).with_name(driver.name),driver)
    cmd=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',
         'PYTHONHASHSEED=0',f'PYTHONPATH={args.dependencies}:{args.run}/engine/src:{args.run}/engine/bench',
         'timeout','--kill-after=5s','180',sys.executable,str(driver),'--run',str(args.run),'--slot']
    options=['sbatch','--parsable','--partition=cpunodes','--nodelist=bosque11',
             '--cpus-per-task=1','--mem=6G','--time=00:05:00','--array=0-29%30',
             '--job-name=ts_score',f'--output={args.run}/status/score_%A_%a.out',
             '--wrap',shlex.join(cmd)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'scoring_submission.json',dict(job=job,command=options,driver_sha256=sha(driver)))
    print(job)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','submit','task','analyze','relocate','membership','core_equivalence','submit_scores'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--slot',type=int)
    parser.add_argument('--source',type=Path,default=Path('/h/399/yunhengzou/appendix_final/bgcp_rerank_latest'))
    parser.add_argument('--slap-core',type=Path,default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/competitor_env_20260908/lib/python3.10/site-packages/slapmapper/core.py'))
    parser.add_argument('--dependencies',type=Path,default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/slap_guided_pilot_20260909/dependencies'))
    args=parser.parse_args();globals()[args.command](args)

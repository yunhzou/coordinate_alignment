"""Reference-free 140-step feasibility benchmark. No accuracy labels are inferred."""
import argparse
from collections import Counter
from dataclasses import asdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import statistics
import tarfile
import subprocess
import sys
import time

from golden_competitors import save

SOURCE = Path('/h/399/yunhengzou/appendix_final')
BENCH_ROOT = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def valid_mapping(elements_r, elements_p, mapping):
    mapping = dict(mapping)
    return (set(mapping) == set(range(len(elements_r)))
            and set(mapping.values()) == set(range(len(elements_p)))
            and len(mapping) == len(elements_p)
            and all(elements_r[a] == elements_p[b] for a, b in mapping.items()))


def prepare(args):
    import numpy as np
    from rxn_core import AAMSearchConfig
    from rxn_core.chemistry_computations.xyz import parse_xyz
    from rxn_core.chemistry_computations.xtb import load_cached_xtb
    args.run.mkdir(parents=True, exist_ok=False)
    (args.run/'inputs').mkdir(); (args.run/'status').mkdir()
    selection = SOURCE/'aam_neb_preservation/sweep_selections_140.json'
    selected = json.loads(selection.read_text())['cases']
    records = []
    for index, case in enumerate(selected):
        directory = args.run/'inputs'/str(index); directory.mkdir()
        endpoints, audit = [], []
        for side, key in [('R', 'reactant_xyz'), ('P', 'product_xyz')]:
            cache = SOURCE/'bgcp_rerank_latest/work'/case['step_id']/'endpoints'/side
            elements, coords, wbo, xyz = load_cached_xtb(cache)
            expected_elements, expected_coords = parse_xyz(case[key])
            assert list(elements) == list(expected_elements) and np.array_equal(coords, expected_coords)
            for source, name in [(xyz, f'{side}.xyz'), (cache/'wbo', f'{side}.wbo')]:
                shutil.copy2(source, directory/name)
            endpoints.append(dict(elements=list(elements), coordinates=coords.tolist(),
                                  wbo=wbo.tolist(), label=side))
            audit.append(dict(side=side, source_cache=str(cache), atoms=len(elements)))
        assert Counter(endpoints[0]['elements']) == Counter(endpoints[1]['elements'])
        save(directory/'input.json', dict(name=case['step_id'], reactant=endpoints[0], product=endpoints[1]))
        # Original component files are supplied to SLAP's native XYZ interface.
        components = {}
        for side, role in [('R', 'reactants'), ('P', 'products')]:
            paths = []
            for ordinal, source in enumerate(sorted((SOURCE/'benchmark'/case['step_id']/role).glob('*.xyz'))):
                dest = directory/f'{side}_component_{ordinal}.xyz'; shutil.copy2(source, dest)
                paths.append(str(dest.resolve()))
            assembled = [e for p in paths for e in parse_xyz(p)[0]]
            assert assembled == endpoints[0 if side == 'R' else 1]['elements']
            components[side] = paths
        save(directory/'components.json', components)
        records.append(dict(index=index, step_id=case['step_id'], endpoints=audit,
            full_composition=True, atoms=len(endpoints[0]['elements']),
            sha256={p.name:digest(p) for p in directory.iterdir() if p.is_file()}))
    config=asdict(AAMSearchConfig(seed_count=10, branch_limit=100, iso_tolerance=0.5))
    save(args.run/'inputs/manifest.json', dict(config=config, records=records))
    for folder in ('src', 'bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    save(args.run/'manifest.json',dict(records=records, config=config, root_seed=42,
        source_selection=str(selection), source_selection_sha256=digest(selection),
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        engine_sha256={str(p.relative_to(args.run/'engine')):digest(p)
            for p in (args.run/'engine').rglob('*') if p.is_file()},
        reference_available=False, accuracy=None, watchdog_seconds=300,
        protocols=dict(aam='Cached continuous WBO; explicit H; two directions; seed 10; cap 100; tol 0.5; sweep cut.',
            slap_xyz='Unmodified native map_3d; original components; binary adjacency; bond_scale=1.2; heavy symmetry.'),
        previous_predictions_used=False, groundtruth_ts_used=False))
    print(json.dumps(dict(cases=len(records))))


def aam(args):
    from golden_publication import search
    search(args)  # Existing search, profiler and checkpoint contract, no reference reads.
    from golden_policy_campaign import load_case
    from rxn_core.artifacts import read_aam_checkpoint
    _, plan = load_case(args.run/'inputs', args.index)
    problem = plan.input_problem
    out=args.run/'directions'/str(args.index)/args.direction
    result=read_aam_checkpoint(out/'cuts/aam.pkl.gz')
    started=time.perf_counter(); cpu=time.process_time()
    count=0
    with gzip.open(out/'terminal_mappings.jsonl.gz','wt') as stream:
        for terminal in result.graph.terminals:
            mapping=dict(result.graph.states[terminal].mapping)
            if args.direction=='P_to_R': mapping={b:a for a,b in mapping.items()}
            valid=valid_mapping(problem.reactant.elements,problem.product.elements,mapping)
            count+=valid
            stream.write(json.dumps(dict(terminal=terminal,mapping=sorted(mapping.items()),
                                         full_element_bijection=valid))+'\n')
    save(out/'feasibility.json',dict(valid_full_terminal_mappings=count,
        terminals=len(result.graph.terminals), capped=result.graph.capped, accuracy=None,
        note='Terminal witnesses are not an exhaustive list of symmetry alternatives; full archives retained.',
        witness_export_including_io_seconds=time.perf_counter()-started,
        witness_export_including_io_cpu_seconds=time.process_time()-cpu))


def slap_xyz(args):
    from slapmapper.aam import SlapAAM
    files=json.loads((args.run/'inputs'/str(args.index)/'components.json').read_text())
    mapper=SlapAAM(binary=True)
    out=args.run/'slap_xyz';out.mkdir(exist_ok=True)
    wall,cpu=time.perf_counter(),time.process_time()
    try:
        mapper.map_3d(files['R'],files['P'],break_sym='heavy',base=0)
        seconds,cpu_seconds=time.perf_counter()-wall,time.process_time()-cpu
        candidates=[]
        for result in mapper.results:
            graphs=[dict(labels=[int(x) for x in g.labels],
                elements=[int(x) for x in g.props['atomic numbers']],
                edges=[(int(a),int(b),float(w)) for a,nb in g.graph.items() for b,w in nb.items() if a<b])
                for g in result['lgp']]
            # Equal labeled element multisets establish existence of a full
            # element-preserving assignment within the native compressed groups.
            valid=Counter(zip(graphs[0]['labels'],graphs[0]['elements']))==Counter(zip(graphs[1]['labels'],graphs[1]['elements']))
            candidates.append(dict(mapping=result['mapping'],val=float(result['val']),
                cd=float(result['cd']),graphs=graphs,full_element_assignment_exists=valid))
        row=dict(status='mapped',candidates=candidates,mapping_seconds=seconds,mapping_cpu_seconds=cpu_seconds)
    except Exception as error:
        import traceback
        row=dict(status='mapping_error',error=str(error),traceback=traceback.format_exc(),
            mapping_seconds=time.perf_counter()-wall,mapping_cpu_seconds=time.process_time()-cpu)
    save(out/f'{args.index}.json',dict(row,index=args.index,accuracy=None))


def worker(args):
    if args.method=='aam': args.index,args.direction=divmod(args.slot,2);args.direction=('R_to_P','P_to_R')[args.direction]
    else: args.index=args.slot
    label=f'{args.method}_{args.index}_{args.direction}'
    status=args.run/'status'/f'{label}.json'
    if status.exists() and json.loads(status.read_text())['exit']==0:
        return
    command=[sys.executable,str(Path(__file__).resolve()),args.method,'--run',str(args.run),
             '--index',str(args.index),'--direction',args.direction]
    with (args.run/'status'/f'{label}.log').open('a') as log:
        child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            code=child.wait(timeout=300)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid,signal.SIGKILL);child.wait();code='timeout'
    save(status,dict(exit=code,index=args.index,direction=args.direction,
        method=args.method,job=os.environ.get('SLURM_JOB_ID'),finished=time.time()))


def submit(args):
    manifest=json.loads((args.run/'manifest.json').read_text());n=len(manifest['records'])
    engine=args.run/'engine';jobs=[]
    for method in ('aam','slap_xyz'):
        cpus=16 if method=='aam' else 1
        python=(Path(sys.executable) if method=='aam' else BENCH_ROOT/'competitor_env_20260908/bin/python')
        env=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',
             'PYTHONHASHSEED=0','RXN_CORE_NATIVE=1',f'PYTHONPATH={args.run}/dependencies:{engine}/src:{engine}/bench',
             'CUDA_VISIBLE_DEVICES=']
        command=[str(python),str(engine/'bench/elementary_feasibility.py'),'worker','--run',str(args.run),
                 '--method',method,'--slot']
        count=n*2 if method=='aam' else n
        options=['sbatch','--parsable','--partition=cpunodes_nia,cpunodes','--nodes=1',
            f'--cpus-per-task={cpus}','--mem=32G' if method=='aam' else '--mem=8G',
            '--time=00:10:00',f'--array=0-{count-1}%32',f'--job-name=elem_{method}',
            f'--output={args.run}/status/{method}_%A_%a.out',
            '--wrap',shlex.join([*env,*command])+' "$SLURM_ARRAY_TASK_ID"']
        job=subprocess.check_output(options,text=True).strip()
        jobs.append(dict(method=method,job=job,command=options));save(args.run/'jobs.json',jobs)
    print(json.dumps(jobs,indent=2))


def freeze(args):
    """Refresh driver snapshot only before any cluster submissions."""
    assert not (args.run/'jobs.json').exists()
    for name in ('elementary_feasibility.py','golden_competitors.py'):
        shutil.copy2(Path(__file__).with_name(name),args.run/'engine/bench'/name)
    manifest=json.loads((args.run/'manifest.json').read_text())
    manifest['engine_sha256']={str(p.relative_to(args.run/'engine')):digest(p)
        for p in (args.run/'engine').rglob('*') if p.is_file()}
    manifest['driver_freeze_git_commit']=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    save(args.run/'manifest.json',manifest)


def report(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    records=[]
    for case in manifest['records']:
        index=case['index'];row=dict(index=index,step_id=case['step_id'],atoms=case['atoms'],methods={})
        for direction in ('R_to_P','P_to_R'):
            out=args.run/'directions'/str(index)/direction
            status=args.run/'status'/f'aam_{index}_{direction}.json'
            result=dict(status='pending')
            if status.exists():
                result['status']=str(json.loads(status.read_text())['exit'])
            if (out/'feasibility.json').exists():
                result.update(json.loads((out/'feasibility.json').read_text()),status='finished')
                timing=json.loads((out/'search.json').read_text())
                result['compute_cpu_seconds']=timing['compute_cpu_excluding_persistence_and_loading_seconds']
                result['elapsed_including_io_seconds']=timing['elapsed_wall_including_io_seconds']
            row['methods']['aam_'+direction]=result
        for method in ('slap_xyz',):
            path=args.run/method/f'{index}.json'
            if not path.exists():
                row['methods'][method]=dict(status='pending');continue
            raw=json.loads(path.read_text());result={k:v for k,v in raw.items() if k not in ('candidates','traceback')}
            candidates=raw.get('candidates',[]);result['returned_candidates']=len(candidates)
            result['valid_full_candidates']=sum(c['full_element_assignment_exists'] for c in candidates)
            row['methods'][method]=result
        records.append(row)
    summaries={}
    for method in records[0]['methods']:
        values=[r['methods'][method] for r in records]
        summaries[method]=dict(statuses=dict(Counter(v['status'] for v in values)),
            cases_with_full_assignment=sum(v.get('valid_full_terminal_mappings',v.get('valid_full_candidates',0))>0 for v in values),
            returned_call_cpu_seconds=sum(v.get('compute_cpu_seconds',v.get('mapping_cpu_seconds',0)) for v in values),
            accuracy=None)
        times=[v.get('compute_cpu_seconds',v.get('mapping_cpu_seconds')) for v in values]
        times=[x for x in times if x is not None]
        summaries[method]['cpu_seconds_per_completed_call']=dict(count=len(times),
            mean=statistics.mean(times) if times else None,
            median=statistics.median(times) if times else None,max=max(times) if times else None)
        summaries[method]['capped_cases']=(sum(v.get('capped',False) for v in values)
                                          if method.startswith('aam_') else None)
    save(args.run/'feasibility_per_case.json',records)
    save(args.run/'feasibility_summary.json',dict(cases=len(records),methods=summaries,accuracy=None,
        note='Full assignment is structural feasibility, not correctness. Only native XYZ/WBO and native XYZ runs are included. '
             'Unrequested MOL-derived SMILES runs were cancelled and excluded. '
             'CPU excludes killed workers without final readings; inspect Slurm accounting.'))
    print(json.dumps(summaries,indent=2))


def publish(args):
    report(args)
    destination=Path('reports/elementary140_feasibility_20260908');destination.mkdir(exist_ok=True)
    for name in ('feasibility_per_case.json','feasibility_summary.json','manifest.json',
                 'jobs.json','relocation.json','WITHDRAWN_SMILES_RUNS.md'):
        shutil.copy2(args.run/name,destination/name)
    jobs=[j['job'] for filename in ('jobs.json','relocation.json')
          for j in json.loads((args.run/filename).read_text()) if j['method'] in ('aam','slap_xyz')]
    accounting=subprocess.check_output(['sacct','-j',','.join(jobs),'-P',
        '--format=JobID,State,ExitCode,ElapsedRaw,TotalCPU,CPUTimeRAW,AllocCPUS,MaxRSS,NodeList'],text=True)
    (args.run/'slurm_native_accounting.psv').write_text(accounting)
    shutil.copy2(args.run/'slurm_native_accounting.psv',destination/'slurm_native_accounting.psv')
    shutil.copy2(__file__,args.run/'report_driver.py')
    with tarfile.open(destination/'native_slap_outputs.tar.gz','w:gz') as archive:
        for path in sorted((args.run/'slap_xyz').glob('*.json')):
            archive.add(path,arcname='slap_xyz/'+path.name)
    save(destination/'report_provenance.json',dict(script_sha256=digest(Path(__file__)),
        ase_version='3.26.0',ase_location=str(args.run/'dependencies'),
        source_run=str(args.run),accuracy=None))


def relocate(args):
    """Move only unstarted native allocations; never repeat completed searches."""
    jobs=[j for j in json.loads((args.run/'jobs.json').read_text()) if j['method'] in ('aam','slap_xyz')]
    recovered=json.loads((args.run/'relocation.json').read_text()) if (args.run/'relocation.json').exists() else []
    jobs=jobs+list(recovered)
    for job in jobs:
        pending=subprocess.check_output(['squeue','-h','-j',job['job'],'-t','CONFIGURING','-o','%i'],text=True).split()
        slots=[]
        for identifier in pending:
            slot=int(identifier.split('_')[1]);index=slot//2 if job['method']=='aam' else slot
            direction=('R_to_P','P_to_R')[slot%2] if job['method']=='aam' else 'R_to_P'
            if job['method']=='aam':
                assert not (args.run/'directions'/str(index)/direction/'environment.json').exists()
            else:
                assert not (args.run/'slap_xyz'/f'{index}.json').exists()
            subprocess.run(['scancel',identifier],check=True);slots.append(slot)
        if not slots: continue
        command=[('--partition=cpunodes' if x.startswith('--partition=') else
                  '--array='+','.join(map(str,slots))+'%12' if x.startswith('--array=') else x)
                 for x in job['command']]
        if args.exclude: command.insert(2,'--exclude='+args.exclude)
        replacement=subprocess.check_output(command,text=True).strip()
        recovered.append(dict(method=job['method'],replaced_unstarted_slots=slots,
            original_job=job['job'],job=replacement,command=command,reason='Node setup stalled; no search outputs created.'))
        save(args.run/'relocation.json',recovered)
    print(json.dumps(recovered,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','aam','slap_xyz','worker','submit','freeze','report','relocate','publish'))
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--index',type=int);p.add_argument('--slot',type=int)
    p.add_argument('--exclude')
    p.add_argument('--direction',default='R_to_P');p.add_argument('--method')
    args=p.parse_args();globals()[args.command](args)

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
import subprocess
import sys
import time

from golden_competitors import save

SOURCE = Path('/h/399/yunhengzou/appendix_final')
BENCH_ROOT = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
SMILES_METHODS = ('rxnmapper', 'localmapper', 'slap_binary', 'slap_weighted', 'indigo', 'chython', 'rdt')


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
    from rdkit import Chem, RDLogger
    from rxn_core import AAMSearchConfig
    from rxn_core.chemistry_computations.xyz import parse_xyz
    from rxn_core.chemistry_computations.xtb import load_cached_xtb
    RDLogger.DisableLog('rdApp.*')
    args.run.mkdir(parents=True, exist_ok=False)
    (args.run/'inputs').mkdir(); (args.run/'status').mkdir()
    selection = SOURCE/'aam_neb_preservation/sweep_selections_140.json'
    selected = json.loads(selection.read_text())['cases']
    records, smiles_rows = [], []
    for index, case in enumerate(selected):
        directory = args.run/'inputs'/str(index); directory.mkdir()
        endpoints, smiles, audit = [], [], []
        for side, key in [('R', 'reactant_xyz'), ('P', 'product_xyz')]:
            cache = SOURCE/'bgcp_rerank_latest/work'/case['step_id']/'endpoints'/side
            elements, coords, wbo, xyz = load_cached_xtb(cache)
            expected_elements, expected_coords = parse_xyz(case[key])
            assert list(elements) == list(expected_elements) and np.array_equal(coords, expected_coords)
            for source, name in [(xyz, f'{side}.xyz'), (cache/'wbo', f'{side}.wbo'),
                                 (cache/'xtbtopo.mol', f'{side}.mol')]:
                shutil.copy2(source, directory/name)
            endpoints.append(dict(elements=list(elements), coordinates=coords.tolist(),
                                  wbo=wbo.tolist(), label=side))
            mol = Chem.MolFromMolFile(str(directory/f'{side}.mol'), sanitize=False, removeHs=False)
            assert mol is not None and [a.GetSymbol() for a in mol.GetAtoms()] == list(elements)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0); atom.SetNoImplicit(True)
            mol.UpdatePropertyCache(strict=False)
            text = Chem.MolToSmiles(mol, canonical=False, allHsExplicit=True)
            check = Chem.Mol(mol)
            try:
                Chem.SanitizeMol(check); sanitization_error = None
            except Exception as error:
                sanitization_error = str(error)
            smiles.append(text)
            audit.append(dict(side=side, mol_sanitization_error=sanitization_error,
                mol_bond_types=dict(Counter(str(b.GetBondType()) for b in mol.GetBonds())),
                source_cache=str(cache), atoms=len(elements)))
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
        smiles_rows.append(dict(index=index, input_reaction='>>'.join(smiles)))
    config=asdict(AAMSearchConfig(seed_count=10, branch_limit=100, iso_tolerance=0.5))
    save(args.run/'inputs/manifest.json', dict(config=config, records=records))
    with (args.run/'inputs.jsonl').open('w') as stream:
        for row in smiles_rows: stream.write(json.dumps(row)+'\n')
    for folder in ('src', 'bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    save(args.run/'manifest.json',dict(records=records, config=config, root_seed=42,
        source_selection=str(selection), source_selection_sha256=digest(selection),
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        engine_sha256={str(p.relative_to(args.run/'engine')):digest(p)
            for p in (args.run/'engine').rglob('*') if p.is_file()},
        reference_available=False, accuracy=None, watchdog_seconds=300,
        protocols=dict(aam='Cached continuous WBO; explicit H; two directions; seed 10; cap 100; tol 0.5; sweep cut.',
            slap_xyz='Unmodified native map_3d; original components; binary adjacency; bond_scale=1.2; heavy symmetry.',
            smiles='Cached xTB MOL exported without valence repairs; all explicit atoms, no implicit H. '
                   'Input compatibility only; not equivalent to native WBO graphs or curated chemical SMILES.'),
        previous_predictions_used=False, groundtruth_ts_used=False))
    print(json.dumps(dict(cases=len(records), endpoints_with_valence_error=sum(
        bool(e['mol_sanitization_error']) for r in records for e in r['endpoints']))))


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
    old=json.loads((BENCH_ROOT/'golden_competitors_full_20260908/rdt_submission.json').read_text())
    for method in ('aam','slap_xyz',*SMILES_METHODS):
        cpus=16 if method=='aam' else 8 if method in ('rdt','chython') else 1
        python=(Path(sys.executable) if method=='aam' else BENCH_ROOT/(
            'neural_competitor_env_20260908' if method in ('rxnmapper','localmapper')
            else 'competitor_env_20260908')/'bin/python')
        env=['env','OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',
             'PYTHONHASHSEED=0','RXN_CORE_NATIVE=1',f'PYTHONPATH={args.run}/dependencies:{engine}/src:{engine}/bench',
             'CUDA_VISIBLE_DEVICES=']
        if method=='rdt': env += [f'RDT_JAVA={old["java"]}',f'RDT_CLASSPATH={old["classpath"]}']
        if method in ('aam','slap_xyz'):
            command=[str(python),str(engine/'bench/elementary_feasibility.py'),'worker','--run',str(args.run),
                     '--method',method,'--slot']
        else:
            command=[str(python),str(engine/'bench/golden_competitors.py'),'worker','--run',str(args.run),
                     '--method',method,'--shards',str(n),'--shard']
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


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','aam','slap_xyz','worker','submit','freeze'))
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--index',type=int);p.add_argument('--slot',type=int)
    p.add_argument('--direction',default='R_to_P');p.add_argument('--method')
    args=p.parse_args();globals()[args.command](args)

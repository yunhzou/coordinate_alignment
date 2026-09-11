"""Single-edge SLAP ablation on the unchanged native XYZ graphs of 140 cases.

Both directions, uncut controls, and every individual edge including H bonds.
Mapping sees only XYZ inputs. All scores use the original full-H WBO matrices.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import shlex
import shutil
import socket
import subprocess
import time
import traceback

DATA = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
PAIRED = DATA/'holdout_cap1000_seed1_20260910'
PYTHON = '/project/yunhengzou/coordinate_alignment/.venv/bin/python'
SLAP_PYTHON = str(DATA/'competitor_env_20260908/bin/python')
PACKAGE = DATA/'competitor_env_20260908/lib/python3.10/site-packages/slapmapper'
DEPENDENCIES = DATA/'elementary140_feasibility_20260908/dependencies'
DIRECTIONS = ('R_to_P','P_to_R')


def read(path):
    return json.loads(path.read_text())


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(path)


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def records(path):
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_bytes().splitlines(keepends=True) if line.endswith(b'\n')]


def prepare(args):
    args.run.mkdir(parents=True,exist_ok=False)
    for name in ('status','slap_xyz','slap_scoring','evaluations'):
        (args.run/name).mkdir()
    source = Path(read(PAIRED/'manifest.json')['source'])
    (args.run/'inputs').symlink_to(source/'inputs',target_is_directory=True)
    shutil.copytree(PACKAGE,args.run/'vendor/slapmapper',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(Path(__file__).with_name('slap_edge_sweep.py'),args.run/'slap_edge_sweep.py')
    shutil.copy2(PAIRED/'driver.py',args.run/'score_driver.py')
    shutil.copy2(__file__,args.run/'driver.py')
    save(args.run/'manifest.json',dict(schema='native_xyz_slap_sweep/v1',cases=140,
        paired_baseline=str(PAIRED),source=str(source),reference_available=False,
        input_hashes=read(PAIRED/'input_hashes.json'),
        graph_protocol='Unchanged upstream geoms2lgp, original XYZ components, bond_scale=1.2, base=0; '
                       'binary native XYZ adjacency. No WBO graph substitution.',
        cut_protocol='Both directions; uncut plus every individual source edge, including H bonds. '
                     'Rebuild initial WL labels after each cut. Restore original endpoints for scoring.',
        mapping_protocol='Unmodified SlapAAM(binary=True).get_maps with heavy symmetry breaking.',
        seed=20260910,python_hash_seed=0,search_watchdog_per_direction=300,scoring_watchdog=600,
        timing='Single CPU. Sum initial XYZ graph preparation, per-cut graph construction, mapping and '
               'native label export CPU; exclude compressed artifact encoding/writing and offline scoring. '
               'Each case runs on the same host as its paired AAM/native-SLAP baseline.',
        vendor_sha256={str(p.relative_to(args.run)):sha(p) for p in (args.run/'vendor').rglob('*.py')},
        drivers_sha256={name:sha(args.run/name) for name in ('driver.py','score_driver.py','slap_edge_sweep.py')}))
    print(args.run,flush=True)


def direction(args):
    import random
    import numpy as np
    from slapmapper.aam import SlapAAM
    from slapmapper.aam._geom import geoms2lgp
    from slap_edge_sweep import cut_graphs
    random.seed(20260910)
    np.random.seed(20260910)
    folder = args.run/f'outputs/{args.index}/{args.direction}'
    folder.mkdir(parents=True,exist_ok=False)
    files = read(args.run/f'inputs/{args.index}/components.json')
    raw = read(args.run/f'inputs/{args.index}/input.json')
    reverse = args.direction == 'P_to_R'
    cpu,wall = time.process_time(),time.perf_counter()
    base = geoms2lgp(files['P'] if reverse else files['R'],files['R'] if reverse else files['P'],base=0,bond_scale=1.2)
    preparation_cpu,preparation_wall = time.process_time()-cpu,time.perf_counter()-wall
    edges = sorted((a,b) for a,ns in base[0].graph.items() for b in ns if a<b)
    elements = [list(g.props['atomic numbers']) for g in base]
    assert all(w == 1 for g in base for ns in g.graph.values() for w in ns.values())
    save(folder/'input.json',dict(index=args.index,direction=args.direction,expected_calls=len(edges)+1,
        edges=edges,preparation_cpu=preparation_cpu,preparation_wall=preparation_wall,
        original_graphs=[dict(labels=list(g.labels),edges=[(a,b,w) for a,ns in g.graph.items()
                                                        for b,w in ns.items() if a<b]) for g in base]))
    mapper = SlapAAM(binary=True)
    targets = [i for i,z in enumerate(elements[0]) if z>1]
    with (folder/'records.jsonl').open('wb') as journal,(folder/'native.bin').open('wb') as archive:
        for ordinal,edge in enumerate([None,*edges]):
            cpu,wall = time.process_time(),time.perf_counter()
            graphs = cut_graphs(base,edge)
            graph_cpu,graph_wall = time.process_time()-cpu,time.perf_counter()-wall
            cpu,wall = time.process_time(),time.perf_counter()
            try:
                mapper.get_maps(graphs,break_sym_targets=targets,base=0)
                status,error = 'mapped',None
            except Exception:
                status,error = 'mapping_error',traceback.format_exc()
            mapping_cpu,mapping_wall = time.process_time()-cpu,time.perf_counter()-wall
            cpu,wall = time.process_time(),time.perf_counter()
            native = []
            if status == 'mapped':
                for result in mapper.results:
                    labels = [list(map(int,g.labels)) for g in result['lgp']]
                    assert Counter(zip(labels[0],elements[0])) == Counter(zip(labels[1],elements[1]))
                    if reverse:
                        labels.reverse()
                    assert len(labels[0]) == len(raw['reactant']['elements'])
                    assert len(labels[1]) == len(raw['product']['elements'])
                    native.append(dict(graphs=[dict(labels=values) for values in labels],
                        val=float(result['val']),lap_sols=result['lap_sols']))
            export_cpu,export_wall = time.process_time()-cpu,time.perf_counter()-wall
            cpu,wall = time.process_time(),time.perf_counter()
            blob = gzip.compress(pickle.dumps(native,protocol=5),compresslevel=1,mtime=0)
            encoding_cpu,encoding_wall = time.process_time()-cpu,time.perf_counter()-wall
            offset = archive.tell()
            archive.write(blob)
            archive.flush()
            row = dict(ordinal=ordinal,cut=edge,status=status,error=error,candidate_count=len(native),
                graph_cpu=graph_cpu,graph_wall=graph_wall,mapping_cpu=mapping_cpu,mapping_wall=mapping_wall,
                export_cpu=export_cpu,export_wall=export_wall,encoding_cpu=encoding_cpu,encoding_wall=encoding_wall,
                native=dict(offset=offset,length=len(blob),sha256=hashlib.sha256(blob).hexdigest()))
            journal.write((json.dumps(row)+'\n').encode())
            journal.flush()
    save(folder/'complete.json',dict(complete=True,calls=len(edges)+1))


def label_key(candidate):
    groups = [defaultdict(list),defaultdict(list)]
    for side,graph in enumerate(candidate['graphs']):
        for atom,label in enumerate(graph['labels']):
            groups[side][label].append(atom)
    assert set(groups[0]) == set(groups[1])
    return tuple(sorted((tuple(left),tuple(groups[1][label])) for label,left in groups[0].items()))


def aggregate(args):
    candidates,sources,known = [],[],{}
    cpu,wall = time.process_time(),time.perf_counter()
    for direction_name in DIRECTIONS:
        folder = args.run/f'outputs/{args.index}/{direction_name}'
        with (folder/'native.bin').open('rb') as stream:
            for row in records(folder/'records.jsonl'):
                stream.seek(row['native']['offset'])
                blob = stream.read(row['native']['length'])
                assert hashlib.sha256(blob).hexdigest() == row['native']['sha256']
                values = pickle.loads(gzip.decompress(blob))
                assert len(values) == row['candidate_count']
                for ordinal,candidate in enumerate(values):
                    key = label_key(candidate)
                    if key not in known:
                        known[key] = len(candidates)
                        candidates.append(dict(graphs=candidate['graphs']))
                        sources.append([])
                    sources[known[key]].append(dict(direction=direction_name,cut=row['cut'],
                        call_ordinal=row['ordinal'],native_candidate=ordinal))
    save(args.run/f'slap_xyz/{args.index}.json',dict(status='mapped',candidates=candidates))
    save(args.run/f'evaluations/{args.index}_sources.json',dict(sources=sources,
        aggregation_cpu_including_read_decode=time.process_time()-cpu,
        aggregation_wall_including_read_decode=time.perf_counter()-wall))


def case(args):
    status = args.run/f'status/{args.index}.json'
    assert not status.exists(), 'Never overwrite an attempted case'
    paired = read(PAIRED/f'status/{args.index}.json')
    assert socket.gethostname() == paired['host'], 'Match the paired baseline host'
    record = dict(index=args.index,host=socket.gethostname(),affinity=sorted(os.sched_getaffinity(0)),
        started=time.time(),job=os.environ.get('SLURM_JOB_ID'),phases={})
    save(status,record)
    env = dict(os.environ,PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',
        OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
        PYTHONPATH=f'{args.run}/vendor:{DEPENDENCIES}')
    score_env = dict(env,PYTHONPATH=f'{PAIRED}/original/src:{PAIRED}/engine/bench')
    commands = [(d,300,[SLAP_PYTHON,str(args.run/'driver.py'),'direction','--run',str(args.run),
                       '--index',str(args.index),'--direction',d],env) for d in DIRECTIONS]
    commands.extend([
        ('aggregate',120,[PYTHON,str(args.run/'driver.py'),'aggregate','--run',str(args.run),'--index',str(args.index)],env),
        ('score',600,[PYTHON,str(args.run/'score_driver.py'),'score_slap','--run',str(args.run),'--index',str(args.index)],score_env)])
    for label,watchdog,command,phase_env in commands:
        start = time.perf_counter()
        with (args.run/f'status/{args.index}_{label}.log').open('w') as log:
            code = subprocess.run(['timeout','--kill-after=5s',str(watchdog),*command],env=phase_env,
                                  stdout=log,stderr=subprocess.STDOUT).returncode
        record['phases'][label] = dict(exit=code,elapsed_including_startup_io=time.perf_counter()-start)
        save(status,record)
    record['finished'] = time.time()
    save(status,record)


def submit(args):
    assert not (args.run/'submissions.json').exists()
    groups = defaultdict(list)
    for index in range(140):
        path = args.run/f'status/{index}.json'
        if path.exists():
            value = read(path)
            assert 'finished' in value and all(p['exit'] == 0 for p in value['phases'].values())
            continue
        groups[read(PAIRED/f'status/{index}.json')['host']].append(index)
    submissions = []
    for host,indices in sorted(groups.items()):
        command = ['sbatch','--parsable','--partition=cpunodes_nia',f'--nodelist={host}',
            '--nodes=1','--cpus-per-task=1','--mem=8G','--time=00:25:00','--no-requeue',
            '--array='+','.join(map(str,indices))+'%8','--job-name=holdout_slap_sweep',
            f'--output={args.run}/status/slurm_%A_%a.out','--wrap',
            shlex.join([PYTHON,str(args.run/'driver.py'),'case','--run',str(args.run),'--index'])+' "$SLURM_ARRAY_TASK_ID"']
        job = subprocess.check_output(command,text=True).strip()
        submissions.append(dict(host=host,indices=indices,job=job,command=command))
        save(args.run/'submissions.json',submissions)
        print(host,job,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','direction','aggregate','case','submit'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int)
    parser.add_argument('--direction',choices=DIRECTIONS)
    args = parser.parse_args()
    globals()[args.command](args)

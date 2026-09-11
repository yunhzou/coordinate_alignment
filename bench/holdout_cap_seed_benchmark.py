"""Fresh paired 140-case cap/seed benchmark against native XYZ SLAP.

Reuses the frozen stable search and scoring adapters. Search never reads prior
scores. Every case runs both methods on the same Slurm allocation; CPU, search
elapsed time, and offline scoring are recorded separately.
"""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import socket
import statistics
import subprocess
import sys
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
DATA = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
BASELINE = DATA/'adaptive_full_20260910'
SOURCE = DATA/'elementary140_tol1_20260909'
PYTHON = '/project/yunhengzou/coordinate_alignment/.venv/bin/python'
SLAP_PYTHON = str(DATA/'competitor_env_20260908/bin/python')


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def prepare(args):
    old = read(BASELINE/'manifest.json')
    assert old['original_commit'].startswith('98b01b1')
    args.run.mkdir(parents=True, exist_ok=False)
    shutil.copytree(BASELINE/'original', args.run/'original',
                    ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(BASELINE/'engine/bench', args.run/'engine/bench',
                    ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(__file__, args.run/'driver.py')
    (args.run/'inputs').mkdir()
    (args.run/'inputs/holdout').symlink_to(BASELINE/'inputs/holdout', target_is_directory=True)
    input_hashes = {}
    for index in range(140):
        assert sha(BASELINE/f'inputs/holdout/{index}/input.json') == sha(SOURCE/f'inputs/{index}/input.json')
        (args.run/f'inputs/{index}').symlink_to(SOURCE/f'inputs/{index}', target_is_directory=True)
        for path in (SOURCE/f'inputs/{index}').iterdir():
            if path.is_file():
                input_hashes[f'{index}/{path.name}'] = sha(path)
        for paths in read(SOURCE/f'inputs/{index}/components.json').values():
            for name in paths:
                input_hashes[name] = sha(Path(name))
    tasks = [dict(dataset='holdout', index=i, direction=d)
             for i in range(140) for d in ('R_to_P', 'P_to_R')]
    save(args.run/'tasks.json', tasks)
    save(args.run/'input_hashes.json', input_hashes)
    config = dict(old['original_config'], seed_count=args.seeds, branch_limit=args.cap)
    for folder in ('status', 'slap_xyz', 'slap_scoring'):
        (args.run/folder).mkdir()
    manifest = dict(schema='holdout_cap_seed/v1', cases=140, original_commit=old['original_commit'],
        original_config=config, baseline_config=old['original_config'],
        changed_search_fields=['seed_count', 'branch_limit'],
        original_execution=old['original_execution'], original_workers=old['original_workers'],
        root_seed=old['root_seed'], python_hash_seed=0, source=str(SOURCE), baseline=str(BASELINE),
        reference_available=False, chemical_accuracy=None, search_watchdog_seconds=300,
        analysis_watchdog_seconds=240, fresh_slap=True,
        slap_protocol='Unmodified map_3d; original XYZ components; binary=True; break_sym=heavy; base=0; default bond_scale=1.2.',
        score='Shared full-H WBO events: graph threshold 0.2, change threshold 0.5. '
              'AAM best saved full terminal representative over both directions. '
              'SLAP native H-label assignments optimized with the frozen bounded scorer.',
        timing='AAM summed parent plus child compute CPU, excluding measured archive persistence/loading, '
               'over both directions. SLAP native map_3d CPU includes XYZ input processing. '
               'Scoring measured separately. Imports/dispatch excluded from primary compute timers. '
               'Same host per paired case; AAM eight workers, SLAP one. '
               'Wall times are observed search durations, not equal-core latency.',
        frozen_sha256={str(p.relative_to(args.run)):sha(p)
            for directory in ('original', 'engine/bench') for p in (args.run/directory).rglob('*')
            if p.is_file()}, driver_sha256=sha(args.run/'driver.py'))
    save(args.run/'manifest.json', manifest)
    print(json.dumps(dict(run=str(args.run), config=config)), flush=True)


def submit(args):
    assert not (args.run/'submission.json').exists()
    slots = []
    for index in range(140):
        status = args.run/f'status/{index}.json'
        if status.exists():
            row = read(status)
            assert 'finished' in row and all(p['exit'] == 0 for p in row['phases'].values())
        else:
            slots.append(index)
    assert slots
    shutil.copy2(__file__, args.run/'driver.py')
    manifest = read(args.run/'manifest.json')
    manifest['driver_sha256'] = sha(args.run/'driver.py')
    manifest['pilot_cases'] = sorted(set(range(140))-set(slots))
    package = DATA/'competitor_env_20260908/lib/python3.10/site-packages/slapmapper'
    manifest['slap_source_sha256'] = {str(p.relative_to(package)):sha(p) for p in package.rglob('*.py')}
    manifest['slap_python'] = SLAP_PYTHON
    save(args.run/'manifest.json', manifest)
    command = ['sbatch', '--parsable', '--partition=cpunodes_nia', '--exclude=bosque49,bosque56',
        '--nodes=1', '--cpus-per-task=8', '--mem=32G', '--time=00:25:00', '--no-requeue',
        '--array='+','.join(map(str,slots))+'%16', '--job-name=holdout_cap1000',
        f'--output={args.run}/status/slurm_%A_%a.out', '--wrap',
        shlex.join([PYTHON, str(args.run/'driver.py'), 'case', '--run', str(args.run), '--index'])
        +' "$SLURM_ARRAY_TASK_ID"']
    job = subprocess.check_output(command, text=True).strip()
    save(args.run/'submission.json', dict(job=job, command=command))
    print(job, flush=True)


def case(args):
    status = args.run/f'status/{args.index}.json'
    assert not status.exists(), 'Do not overwrite an attempted case'
    record = dict(index=args.index, host=socket.gethostname(), started=time.time(),
        job=os.environ.get('SLURM_JOB_ID'), affinity=sorted(os.sched_getaffinity(0)),
        cpu_models=sorted({line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                           if line.startswith('model name')}), phases={})
    save(status, record)
    bench = args.run/'engine/bench'
    env = dict(os.environ, PYTHONPATH=f'{args.run}/original/src:{bench}', RXN_CORE_NATIVE='1',
        PYTHONHASHSEED='0', PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1',
        OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', CUDA_VISIBLE_DEVICES='')
    commands = [('slap', 300, [SLAP_PYTHON, str(bench/'elementary_feasibility.py'), 'slap_xyz',
                             '--run', str(args.run), '--index', str(args.index)])]
    for direction in range(2):
        for phase, watchdog in [('search', 300), ('analyze', 240)]:
            commands.append((f'{phase}_{direction}', watchdog,
                [PYTHON, str(bench/'adaptive_full_benchmark.py'), phase, '--run', str(args.run),
                 '--slot', str(2*args.index+direction), '--method', 'original']))
    commands.append(('slap_score', 240, [PYTHON, str(args.run/'driver.py'), 'score_slap',
                                        '--run', str(args.run), '--index', str(args.index)]))
    for label, watchdog, command in commands:
        phase_env = dict(env)
        if label == 'slap':
            phase_env['PYTHONPATH'] = str(DATA/'elementary140_feasibility_20260908/dependencies')+':'+env['PYTHONPATH']
        start = time.perf_counter()
        with (args.run/f'status/{args.index}_{label}.log').open('w') as log:
            code = subprocess.run(['timeout', '--kill-after=5s', str(watchdog), *command],
                env=phase_env, stdout=log, stderr=subprocess.STDOUT).returncode
        record['phases'][label] = dict(exit=code, elapsed_including_startup_io=time.perf_counter()-start)
        save(status, record)
    record['finished'] = time.time()
    record['child_cpu_including_all_phases_and_io'] = sum((
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime,
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime))
    save(status, record)


def score_slap(args):
    # The initial label representative and invariance check reproduce the
    # frozen compare_elementary_outputs.compare adapter, without AAM I/O.
    import numpy as np
    from compare_elementary_outputs import event_counts, events_row, features, refine
    raw = read(args.run/f'inputs/{args.index}/input.json')
    native = read(args.run/f'slap_xyz/{args.index}.json')
    assert native['status'] == 'mapped'
    r, p = [np.asarray(raw[k]['wbo']) for k in ('reactant', 'product')]
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    feat = features(raw)
    results = []
    for ordinal, candidate in enumerate(native['candidates']):
        groups = []
        for graph in candidate['graphs']:
            partition = defaultdict(list)
            for i, label in enumerate(graph['labels']):
                partition[label].append(i)
            groups.append(partition)
        mapping = dict(pair for label, left in groups[0].items() for pair in zip(left, groups[1][label]))
        assert sorted(mapping) == list(range(len(r))) and sorted(mapping.values()) == list(range(len(p)))
        assert all(raw['reactant']['elements'][a] == raw['product']['elements'][b] for a,b in mapping.items())
        assert all(sum(raw['reactant']['elements'][i] != 'H' for i in group) <= 1 for group in groups[0].values())
        labels = {(a,b):c for a,b,c in feat[1]['bonds']}
        def edge(a,b):
            return labels.get(tuple(sorted((a,b))))
        invariant = True
        for group in groups[1].values():
            for a,b in zip(group, group[1:]):
                perm = list(range(len(p)))
                perm[a],perm[b] = b,a
                invariant &= all(edge(i,j) == edge(perm[i],perm[j])
                                 for i in (a,b) for j in range(len(p)) if i != j)
        results.append(dict(candidate=ordinal, mapping=sorted(mapping.items()),
            events=events_row(event_counts(r,p,[[mapping[i] for i in range(len(r))]])[0]),
            all_h_label_permutations_score_invariant=bool(invariant)))
    folder = args.run/'slap_scoring'
    save(folder/f'{args.index}.json', dict(index=args.index, slap=results, scope='Fresh native SLAP outputs.'))
    refine(SimpleNamespace(run=folder, source=args.run, index=args.index))
    row = read(folder/f'{args.index}.refined.json')
    row.update(scoring_cpu_seconds=time.process_time()-start_cpu,
               scoring_wall_including_io_seconds=time.perf_counter()-start_wall)
    save(folder/f'{args.index}.refined.json', row)


def collect(args):
    rows = []
    for index in range(140):
        status = read(args.run/f'status/{index}.json')
        assert 'finished' in status
        directions = {}
        for direction in ('R_to_P', 'P_to_R'):
            folder = args.run/f'results/holdout/{index}/{direction}/original'
            search = read(folder/'search.json') if (folder/'search.json').exists() else {}
            score = read(folder/'full_sweep_evaluation.json') if (folder/'full_sweep_evaluation.json').exists() else {}
            timing = search['rows'][-1] if search.get('complete') else {}
            directions[direction] = dict(search_complete=bool(search.get('complete')), evaluation_complete=bool(score),
                best_events=score.get('best_events'), full_representatives=score.get('valid_full_representatives'),
                cpu_seconds=timing.get('compute_cpu_excluding_persistence_and_loading_seconds'),
                search_wall_including_io=timing.get('elapsed_wall_including_io_seconds'),
                scoring_cpu_seconds=score.get('analysis_cpu'), capped=timing.get('capped'))
        native = read(args.run/f'slap_xyz/{index}.json') if (args.run/f'slap_xyz/{index}.json').exists() else {}
        scored_path = args.run/f'slap_scoring/{index}.refined.json'
        scored = read(scored_path) if scored_path.exists() else {}
        old = read(SOURCE/f'comparison/{index}.refined.json')
        complete = all(d['search_complete'] and d['evaluation_complete'] for d in directions.values())
        def total(key):
            return sum(d[key] for d in directions.values()) if all(d[key] is not None for d in directions.values()) else None
        aam = dict(directions=directions, complete=complete,
            best_events=min((d['best_events'] for d in directions.values() if d['best_events'] is not None), default=None),
            cpu_seconds=total('cpu_seconds'), search_wall_including_io=total('search_wall_including_io'),
            scoring_cpu_seconds=total('scoring_cpu_seconds'))
        slap = dict(complete=native.get('status') == 'mapped' and bool(scored),
            best_events=scored.get('slap_min_representative_events'), cpu_seconds=native.get('mapping_cpu_seconds'),
            search_wall_including_io=native.get('mapping_seconds'), scoring_cpu_seconds=scored.get('scoring_cpu_seconds'),
            unresolved_h_candidates=sum(not s['hydrogen_score_optimization']['optimal'] for s in scored.get('slap',[])))
        winner = 'unresolved'
        if complete and slap['complete'] and aam['best_events'] is not None:
            winner = 'aam' if aam['best_events'] < slap['best_events'] else 'slap' if aam['best_events'] > slap['best_events'] else 'tie'
        rows.append(dict(index=index, name=read(args.run/f'inputs/{index}/input.json')['name'],
            aam=aam, slap=slap, winner=winner, previous_aam_cap100_seed10=old['aam_min_saved_events'],
            previous_slap=old['slap_min_representative_events'], environment=status))
    complete = [r for r in rows if r['aam']['complete'] and r['slap']['complete']]
    timing = {}
    for method in ('aam', 'slap'):
        values = [r[method]['cpu_seconds'] for r in complete]
        timing[method] = dict(total_cpu_seconds=sum(values), mean_cpu_seconds=statistics.mean(values),
            median_cpu_seconds=statistics.median(values), max_cpu_seconds=max(values),
            scoring_cpu_seconds=sum(r[method]['scoring_cpu_seconds'] for r in complete),
            summed_search_wall_including_io=sum(r[method]['search_wall_including_io'] for r in complete))
    summary = dict(cases=140, event_comparison=dict(Counter(r['winner'] for r in rows)),
        paired_complete_cases=len(complete), timing=timing,
        aam_cpu_over_slap=timing['aam']['total_cpu_seconds']/timing['slap']['total_cpu_seconds'],
        full_mapping_cases={m:sum(r[m]['best_events'] is not None for r in rows) for m in ('aam','slap')},
        differs_from_previous_slap=[r['index'] for r in rows if r['slap']['best_events'] != r['previous_slap']],
        unresolved_slap_h_candidates=sum(r['slap']['unresolved_h_candidates'] for r in rows),
        phase_exits=dict(Counter(str(p['exit']) for r in rows for p in r['environment']['phases'].values())),
        chemical_accuracy=None, scope=read(args.run/'manifest.json')['timing'])
    save(args.run/'case_metrics.json', rows)
    save(args.run/'summary.json', summary)
    with (args.run/'per_case.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index','name','aam_events','slap_events','winner','aam_cpu_seconds','slap_cpu_seconds',
                         'aam_scoring_cpu_seconds','slap_scoring_cpu_seconds','host'])
        for r in rows:
            writer.writerow([r['index'],r['name'],r['aam']['best_events'],r['slap']['best_events'],r['winner'],
                r['aam']['cpu_seconds'],r['slap']['cpu_seconds'],r['aam']['scoring_cpu_seconds'],
                r['slap']['scoring_cpu_seconds'],r['environment']['host']])
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare','submit','case','score_slap','collect'))
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--index', type=int)
    parser.add_argument('--seeds', type=int, default=1)
    parser.add_argument('--cap', type=int, default=1000)
    args = parser.parse_args()
    globals()[args.command](args)

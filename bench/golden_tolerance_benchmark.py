"""Paired Golden tolerance ablation using the unchanged frozen search/evaluator.

Searches see structures only. Each directional pair runs on the same CPU set,
alternating tolerance order. Reference queries run in a separate SLURM phase.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import sys
import time

from aam_seed_ablation import claim, read, save, sha, result_folder

BASE = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
BASELINE = BASE / 'aam_one_seed_bidirectional_20260910'
PYTHON = '/project/yunhengzou/coordinate_alignment/.venv/bin/python'
TOLERANCES = {'tol_1p0': 1.0, 'tol_1p5': 1.5}


def environment(run):
    return dict(os.environ, PYTHONPATH=f'{run}/original/src:{run}/engine/bench',
                RXN_CORE_NATIVE='1', PYTHONHASHSEED='0', PYTHONDONTWRITEBYTECODE='1',
                OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')


def prepare(args):
    baseline = read(BASELINE / 'manifest.json')
    assert baseline['original_commit'].startswith('98b01b1')
    config = baseline['original_config']
    assert (config['seed_count'], config['branch_limit'], config['iso_tolerance']) == (1, 100, 1.0)
    tasks = [t for t in read(BASELINE / 'tasks.json') if t['dataset'] == 'golden']
    assert len(tasks) == 3702
    args.run.mkdir(parents=True, exist_ok=False)
    for name in ('original', 'engine'):
        shutil.copytree(BASELINE / name, args.run / name,
                        ignore=shutil.ignore_patterns('__pycache__'))
    # Validate the copied engine and adapters against the archived baseline.
    for category, directory in [('original_sha256', 'original'), ('adapter_sha256', 'engine/bench')]:
        for relative, expected in baseline[category].items():
            assert sha(args.run / directory / relative) == expected, relative
    shutil.copy2(__file__, args.run / 'driver.py')
    shutil.copy2(Path(__file__).with_name('aam_seed_ablation.py'), args.run / 'aam_seed_ablation.py')
    save(args.run / 'tasks.json', tasks)
    weights = []
    hashes = {}
    for slot, task in enumerate(tasks):
        folder = BASELINE / f"inputs/golden/{task['index']}"
        raw = read(folder / 'input.json')
        weights.append((max(len(raw[side]['elements']) for side in ('reactant', 'product')), slot))
        if task['index'] not in hashes:
            hashes[task['index']] = {name: sha(folder / name) for name in ('input.json', 'reference.json')}
    save(args.run / 'input_hashes.json', hashes)
    save(args.run / 'work_order.json', [slot for _, slot in sorted(weights, reverse=True)])
    for phase in ('search', 'analyze'):
        (args.run / f'status/{phase}').mkdir(parents=True)
        save(args.run / f'{phase}_queue.json', {'next': 0})
    for key, tolerance in TOLERANCES.items():
        folder = args.run / 'runs' / key
        folder.mkdir(parents=True)
        (folder / 'inputs').symlink_to(BASELINE / 'inputs', target_is_directory=True)
        save(folder / 'tasks.json', tasks)
        variant = dict(config, iso_tolerance=tolerance)
        assert [k for k in config if variant[k] != config[k]] == ([] if tolerance == 1.0 else ['iso_tolerance'])
        save(folder / 'manifest.json', dict(original_config=variant,
             original_workers=baseline['original_workers'], original_execution=baseline['original_execution'],
             root_seed=baseline['root_seed'], original_commit=baseline['original_commit']))
    save(args.run / 'manifest.json', dict(schema='golden_tolerance/v1', baseline=str(BASELINE),
         cases=1851, tasks=len(tasks), tolerances=TOLERANCES, original_config=config,
         original_commit=baseline['original_commit'], original_execution=baseline['original_execution'],
         workers=8, root_seed=42, python_hash_seed=0, explicit_H=True,
         search='Original uncut plus single-edge sweep; both independently searched directions',
         changed_search_fields=['iso_tolerance'],
         tolerance_scope='Fragment WBO compatibility and symmetry grouping; graph floor and scoring unchanged',
         search_watchdog=300, analysis_watchdog=240, cpu_budget=args.cpu_budget,
         node_cpus=32, analysis_node_cpus=8,
         timing='Fresh paired runs on the same host and CPU set; order alternates by directional slot. '
                'CPU includes parent plus children and excludes measured archive persistence/loading. '
                'Each tolerance is measured once per direction; analysis is separate.',
         evaluation='Unchanged frozen Golden heavy-reference chemical-symmetry and partial-annotation evaluator; '
                    '60-second query budget, 1500 ms per symbolic query. Unknowns stay unknown.',
         driver_sha256=sha(args.run / 'driver.py'), baseline_manifest_sha256=sha(BASELINE / 'manifest.json'),
         original_sha256=baseline['original_sha256'], adapter_sha256=baseline['adapter_sha256']))
    print(args.run, flush=True)


def submit(args):
    manifest = read(args.run / 'manifest.json')
    assert not (args.run / 'submissions.json').exists()
    allocations = manifest['cpu_budget'] // manifest['node_cpus']
    assert allocations > 0
    jobs = []
    for phase in ('search', 'analyze'):
        cpus = manifest['node_cpus'] if phase == 'search' else manifest['analysis_node_cpus']
        command = ['sbatch', '--parsable', '--partition=cpunodes_nia', '--exclude=bosque49,bosque56',
                   '--nodes=1', f'--cpus-per-task={cpus}', '--mem=64G', '--time=00:30:00', '--no-requeue',
                   f'--array=0-{allocations-1}%{allocations}', f'--job-name=golden_tol15_{phase}',
                   f'--output={args.run}/status/{phase}_%A_%a.out']
        if jobs:
            command.append('--dependency=afterany:' + jobs[0]['job'].split(';')[0])
        command += ['--wrap', shlex.join([PYTHON, str(args.run / 'driver.py'), 'batch',
                                         '--run', str(args.run), '--phase', phase])]
        job = subprocess.check_output(command, text=True).strip()
        jobs.append(dict(phase=phase, job=job, command=command))
        save(args.run / 'submissions.json', jobs)
        print(phase, job, flush=True)


def batch(args):
    manifest = read(args.run / 'manifest.json')
    tasks = read(args.run / 'tasks.json')
    cores = sorted(os.sched_getaffinity(0))
    width = manifest['workers'] if args.phase == 'search' else 1
    count = manifest['node_cpus'] if args.phase == 'search' else manifest['analysis_node_cpus']
    assert len(cores) >= count and count % width == 0

    def work(affinity):
        while (slot := claim(args.run, args.phase)) is not None:
            task = tasks[slot]
            status = args.run / f'status/{args.phase}/{slot}.json'
            assert not status.exists(), 'Refuse duplicate task execution'
            row = dict(slot=slot, **task, phase=args.phase, host=socket.gethostname(),
                       affinity=affinity, job=os.environ.get('SLURM_JOB_ID'), started=time.time(), variants={})
            save(status, row)
            keys = list(TOLERANCES)
            if slot % 2:
                keys.reverse()
            for key in keys:
                variant = args.run / 'runs' / key
                if args.phase == 'analyze' and not (result_folder(variant, task) / 'search.json').exists():
                    row['variants'][key] = dict(status='search_incomplete')
                    save(status, row)
                    continue
                command = ['taskset', '-c', ','.join(map(str, affinity)), 'timeout', '--kill-after=5s',
                           str(manifest['search_watchdog' if args.phase == 'search' else 'analysis_watchdog']),
                           PYTHON, str(args.run / 'engine/bench/adaptive_full_benchmark.py'), args.phase,
                           '--run', str(variant), '--slot', str(slot), '--method', 'original']
                start = time.monotonic()
                with (status.parent / f'{slot}_{key}.log').open('w') as log:
                    code = subprocess.run(command, env=environment(args.run), stdout=log, stderr=subprocess.STDOUT).returncode
                row['variants'][key] = dict(exit=code, wall_including_startup_io=time.monotonic()-start)
                save(status, row)
            row['finished'] = time.time()
            save(status, row)

    with ThreadPoolExecutor(max_workers=count // width) as pool:
        list(pool.map(work, [cores[i:i+width] for i in range(0, count, width)]))


def union(outcomes):
    return ('recovered' if 'recovered' in outcomes else
            'not_recovered' if outcomes == ['not_recovered', 'not_recovered'] else 'unknown')


def collect(args):
    rows = []
    roots = {key: args.run / 'runs' / key for key in TOLERANCES}
    roots['archived_1p0'] = BASELINE
    for index in range(1851):
        variants = {}
        for key, root in roots.items():
            directions = {}
            for direction in ('R_to_P', 'P_to_R'):
                folder = result_folder(root, dict(dataset='golden', index=index, direction=direction))
                search = read(folder / 'search.json') if (folder / 'search.json').exists() else {}
                evaluation = read(folder / 'full_sweep_evaluation.json') if (folder / 'full_sweep_evaluation.json').exists() else {}
                timing = search.get('rows', [{}])[-1]
                directions[direction] = dict(search_complete=search.get('complete', False),
                    evaluation_complete=bool(evaluation), recovery=evaluation.get('reference_recovery', 'unknown'),
                    mapping_cpu=timing.get('compute_cpu_excluding_persistence_and_loading_seconds'),
                    mapping_wall=timing.get('elapsed_wall_including_io_seconds'),
                    capped=timing.get('capped'), states=timing.get('states'), terminals=timing.get('terminals'),
                    representative_recovery=evaluation.get('representative_recovery'),
                    top1_correct=evaluation.get('top1_correct'), unknown_queries=evaluation.get('unknown_queries'),
                    analysis_cpu=evaluation.get('analysis_cpu'))
            completed = all(d['search_complete'] for d in directions.values())
            variants[key] = dict(directions=directions,
                recovery=union([d['recovery'] for d in directions.values()]),
                searches_complete=completed,
                mapping_cpu=sum(d['mapping_cpu'] for d in directions.values()) if completed else None,
                capped=any(d['capped'] for d in directions.values()))
        rows.append(dict(index=index, variants=variants))
    paired = [r for r in rows if all(r['variants'][k]['searches_complete'] for k in TOLERANCES)]
    metrics = {}
    for key in roots:
        values = [r['variants'][key] for r in rows]
        metrics[key] = dict(recovery=dict(Counter(v['recovery'] for v in values)),
            completed_search_cases=sum(v['searches_complete'] for v in values),
            completed_directional_evaluations=sum(d['evaluation_complete'] for v in values for d in v['directions'].values()),
            capped_cases=sum(v['capped'] for v in values),
            completed_mapping_cpu=sum(v['mapping_cpu'] for v in values if v['searches_complete']),
            unknown_indices=[r['index'] for r in rows if r['variants'][key]['recovery'] == 'unknown'],
            not_recovered_indices=[r['index'] for r in rows if r['variants'][key]['recovery'] == 'not_recovered'])
    paired_cpu = {k: sum(r['variants'][k]['mapping_cpu'] for r in paired) for k in TOLERANCES}
    changes = [dict(index=r['index'], baseline=r['variants']['tol_1p0']['recovery'],
                    trial=r['variants']['tol_1p5']['recovery']) for r in rows
               if r['variants']['tol_1p0']['recovery'] != r['variants']['tol_1p5']['recovery']]
    summary = dict(cases=1851, metrics=metrics, paired_complete_cases=len(paired),
                   paired_mapping_cpu=paired_cpu,
                   cpu_ratio_1p5_over_1p0=paired_cpu['tol_1p5']/paired_cpu['tol_1p0'] if paired_cpu['tol_1p0'] else None,
                   recovery_changes=changes,
                   baseline_reproduction_changes=[r['index'] for r in rows
                       if r['variants']['tol_1p0']['recovery'] != r['variants']['archived_1p0']['recovery']])
    save(args.run / 'comparison/case_metrics.json', rows)
    save(args.run / 'comparison/summary.json', summary)
    with (args.run / 'comparison/per_case.csv').open('w') as stream:
        fields = ['index', 'variant', 'recovery', 'searches_complete', 'mapping_cpu', 'capped']
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for r in rows:
            for key, value in r['variants'].items():
                writer.writerow(dict(index=r['index'], variant=key, **value))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'submit', 'batch', 'collect'))
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--phase', choices=('search', 'analyze'))
    parser.add_argument('--cpu-budget', type=int, default=512)
    args = parser.parse_args()
    globals()[args.command](args)

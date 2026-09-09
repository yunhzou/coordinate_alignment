"""One-factor holdout rerun: frozen original AAM engine, matching tolerance 1.0."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from types import SimpleNamespace


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def prepare(args):
    source = args.source
    args.run.mkdir(parents=True, exist_ok=False)
    for folder in ('engine', 'inputs', 'slap_xyz'):
        shutil.copytree(source / folder, args.run / folder,
                        ignore=shutil.ignore_patterns('__pycache__'))
    (args.run / 'status').mkdir()
    for name in ('manifest.json', 'inputs/manifest.json'):
        original = json.loads((source / name).read_text())
        assert original['config']['iso_tolerance'] == .5
        original['config']['iso_tolerance'] = 1.0
        if name == 'manifest.json':
            original['protocols']['aam'] = 'Frozen original engine; explicit H; both directions; seed 10; cap 100; matching tolerance 1.0; event threshold 0.5; sweep cut.'
            original['tolerance_ab'] = dict(source=str(source), changed_search_settings={'iso_tolerance': [.5, 1.0]},
                native_slap='Byte-identical saved outputs reused; its search does not use our matching tolerance.',
                cap123='No substitution with old cap-200 run; all new cases start at cap 100.')
        save(args.run / name, original)
    for path in (source / 'inputs').glob('*/input.json'):
        assert path.read_bytes() == (args.run / 'inputs' / path.parent.name / path.name).read_bytes()
    # Preserve exact original executable sources, not a newly modified core.
    for path in (source / 'engine/src').rglob('*'):
        if path.is_file() and '__pycache__' not in path.parts:
            assert path.read_bytes() == (args.run / 'engine/src' / path.relative_to(source / 'engine/src')).read_bytes()
    shutil.copy2(__file__, args.run / 'campaign.py')
    for folder in ('comparison', 'overlap'):
        (args.run / folder).mkdir()
    print(args.run)


def submit(args):
    assert not (args.run / 'jobs.json').exists()
    engine = args.run / 'engine'
    env = ['env', 'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1', 'MKL_NUM_THREADS=1',
           'PYTHONHASHSEED=0', 'RXN_CORE_NATIVE=1', f'PYTHONPATH={engine}/src:{engine}/bench', 'CUDA_VISIBLE_DEVICES=']
    command = [*env, sys.executable, str(engine / 'bench/elementary_feasibility.py'),
               'worker', '--run', str(args.run), '--method', 'aam', '--slot']
    options = ['sbatch', '--parsable', '--partition=cpunodes', '--nodes=1', '--cpus-per-task=16',
               '--mem=32G', '--time=00:10:00', '--array=0-279%32', '--job-name=elem_tol1',
               f'--output={args.run}/status/aam_%A_%a.out', '--wrap', shlex.join(command) + ' "$SLURM_ARRAY_TASK_ID"']
    if args.exclude: options.insert(2, '--exclude=' + args.exclude)
    job = subprocess.check_output(options, text=True).strip()
    save(args.run / 'jobs.json', [dict(method='aam', job=job, command=options)])
    print(job)


def analyze(args):
    from compare_elementary_outputs import compare, refine
    from elementary_family_overlap import analyze as overlap
    full = 0
    for direction in ('R_to_P', 'P_to_R'):
        result = json.loads((args.run / 'directions' / str(args.index) / direction / 'feasibility.json').read_text())
        full += result['valid_full_terminal_mappings']
    if not full:
        raise RuntimeError(f'Case {args.index}: neither direction returned a full mapping; investigate before comparing')
    job = SimpleNamespace(index=args.index, source=args.run, cap200=args.run, run=args.run / 'comparison')
    compare(job)
    refine(job)
    overlap(SimpleNamespace(index=args.index, comparison=args.run / 'comparison', run=args.run / 'overlap'))


def status(args):
    rows = [json.loads(p.read_text()) for p in (args.run / 'status').glob('aam_*.json')]
    print(json.dumps(dict(completed=len(rows), exits=dict(Counter(str(r['exit']) for r in rows)),
                         failed=[r for r in rows if r['exit'] != 0],
                         analyzed=len(list((args.run / 'comparison').glob('*.refined.json')))), indent=2))


def submit_analysis(args):
    dest = args.run / 'analysis_engine'
    dest.mkdir()
    shutil.copytree(Path(__file__).resolve().parent, dest / 'bench', ignore=shutil.ignore_patterns('__pycache__'))
    search = json.loads((args.run / 'jobs.json').read_text())[0]['job']
    command = ['env', 'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1', 'MKL_NUM_THREADS=1',
               f'PYTHONPATH={args.run}/engine/src:{dest}/bench', 'timeout', '--kill-after=5s', '300',
               sys.executable, str(dest / 'bench/elementary_tolerance_ab.py'), 'analyze',
               '--run', str(args.run), '--index']
    options = ['sbatch', '--parsable', '--partition=cpunodes', '--cpus-per-task=1', '--mem=8G',
               '--time=00:10:00', '--array=0-139%32', '--job-name=tol1_compare',
               '--dependency=afterany:' + search, '--exclude=bosque5,bosque6,bosque8',
               f'--output={args.run}/comparison/slurm_%A_%a.out', '--wrap', shlex.join(command) + ' "$SLURM_ARRAY_TASK_ID"']
    job = subprocess.check_output(options, text=True).strip()
    save(args.run / 'analysis_submission.json', dict(job=job, command=options))
    print(job)


def summarize(args):
    from compare_elementary_outputs import report as comparison_report
    from elementary_family_overlap import report as overlap_report
    comparison_report(SimpleNamespace(run=args.run / 'comparison'))
    overlap_report(SimpleNamespace(run=args.run / 'overlap'))
    rows, timings = [], []
    previous = args.source.parent / 'elementary140_output_comparison_20260908'
    for index in range(140):
        old = json.loads((previous / f'{index}.refined.json').read_text())
        new = json.loads((args.run / 'comparison' / f'{index}.refined.json').read_text())
        rows.append(dict(index=index, name=new['name'], old_events=old['aam_min_saved_events'],
                         new_events=new['aam_min_saved_events'], slap_events=new['slap_min_representative_events'],
                         shared_best_heavy_pattern=new['best_aam_heavy_pattern_matches_any_slap_modulo_score_preserving_symmetry'],
                         old_cap=200 if index == 123 else 100, new_cap=100))
        for direction in ('R_to_P', 'P_to_R'):
            item = dict(index=index, direction=direction)
            for label, root in [('old', args.source), ('new', args.run)]:
                folder = root / 'directions' / str(index) / direction
                time = json.loads((folder / 'search.json').read_text())
                env = json.loads((folder / 'environment.json').read_text())
                feasibility = json.loads((folder / 'feasibility.json').read_text())
                item[label] = dict(cpu_seconds=time['compute_cpu_excluding_persistence_and_loading_seconds'],
                                   wall_including_io=time['elapsed_wall_including_io_seconds'],
                                   cpu_models=env['cpu_models'], host=env['host'], capped=feasibility['capped'],
                                   full_witnesses=feasibility['valid_full_terminal_mappings'])
            timings.append(item)
    save(args.run / 'tolerance_per_case.json', rows)
    save(args.run / 'tolerance_timings.json', timings)
    summary = dict(cases=140, improved=sum(r['new_events'] < r['old_events'] for r in rows),
                   equal=sum(r['new_events'] == r['old_events'] for r in rows),
                   worse=sum(r['new_events'] > r['old_events'] for r in rows),
                   total_compute_cpu_seconds={k:sum(r[k]['cpu_seconds'] for r in timings) for k in ('old','new')},
                   mean_direction_wall_including_io={k:sum(r[k]['wall_including_io'] for r in timings)/280 for k in ('old','new')},
                   max_direction_wall_including_io={k:max(r[k]['wall_including_io'] for r in timings) for k in ('old','new')},
                   full_directional_calls={k:sum(r[k]['full_witnesses']>0 for r in timings) for k in ('old','new')},
                   matched_cpu_model_calls=sum(r['old']['cpu_models']==r['new']['cpu_models'] for r in timings),
                   timing_scope='280 directional search calls; compute CPU excludes measured persistence/loading; wall includes IO but not scheduling. Hardware recorded separately. Old timing is cap100, including case123; old event score for case123 uses the separately labeled cap200 follow-up.',
                   chemical_accuracy=None)
    save(args.run / 'tolerance_summary.json', summary)
    followup = json.loads((args.run / 'case123_cap200/comparison/123.refined.json').read_text())
    matched = [dict(row) for row in rows]
    matched[123].update(new_events=followup['aam_min_saved_events'], new_cap=200,
                        shared_best_heavy_pattern=followup['best_aam_heavy_pattern_matches_any_slap_modulo_score_preserving_symmetry'])
    save(args.run / 'matched_cap_policy_per_case.json', matched)
    save(args.run / 'matched_cap_policy_summary.json', dict(
        scope='Same cap policy as preceding comparison: cap100 except separately rerun cap200 for case123.',
        improved=sum(r['new_events']<r['old_events'] for r in matched),
        equal=sum(r['new_events']==r['old_events'] for r in matched),
        worse=sum(r['new_events']>r['old_events'] for r in matched),
        vs_slap=dict(equal=sum(r['new_events']==r['slap_events'] for r in matched),
                     aam_fewer=sum(r['new_events']<r['slap_events'] for r in matched),
                     slap_fewer=sum(r['new_events']>r['slap_events'] for r in matched))))
    print(json.dumps(summary, indent=2))


def cap_followup(args):
    out = args.run / 'case123_cap200'
    out.mkdir(); (out / 'status').mkdir()
    shutil.copytree(args.run / 'inputs', out / 'inputs')
    (out / 'engine').symlink_to(args.run / 'engine', target_is_directory=True)
    for name in ('manifest.json', 'inputs/manifest.json'):
        manifest = json.loads((args.run / name).read_text())
        manifest['config']['branch_limit'] = 200
        if name == 'manifest.json':
            manifest['followup_scope'] = 'Only case123, both directions; cap200, tolerance1; never mixed into fixed cap100 results.'
        save(out / name, manifest)
    env = ['env', 'OMP_NUM_THREADS=1', 'OPENBLAS_NUM_THREADS=1', 'MKL_NUM_THREADS=1', 'PYTHONHASHSEED=0',
           'RXN_CORE_NATIVE=1', f'PYTHONPATH={out}/engine/src:{out}/engine/bench']
    cmd = [*env, sys.executable, str(out / 'engine/bench/elementary_feasibility.py'), 'worker',
           '--run', str(out), '--method', 'aam', '--slot']
    options = ['sbatch', '--parsable', '--partition=cpunodes', '--exclude=bosque5,bosque6,bosque8',
               '--cpus-per-task=16', '--mem=32G', '--time=00:10:00', '--array=246-247',
               '--job-name=tol1_cap200', f'--output={out}/status/slurm_%A_%a.out',
               '--wrap', shlex.join(cmd) + ' "$SLURM_ARRAY_TASK_ID"']
    job = subprocess.check_output(options, text=True).strip()
    save(out / 'jobs.json', [dict(method='aam', job=job, command=options)])
    print(job)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=('prepare', 'submit', 'analyze', 'status', 'submit_analysis', 'summarize', 'cap_followup'))
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--source', type=Path, default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_feasibility_20260908'))
    p.add_argument('--exclude'); p.add_argument('--index', type=int)
    args = p.parse_args(); globals()[args.command](args)

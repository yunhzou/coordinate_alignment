"""Publish a completed Golden tolerance run and check its recovered witnesses."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shutil
import statistics
import subprocess

from golden_tolerance_benchmark import TOLERANCES, PYTHON, read, save, sha, result_folder


def publish(args):
    manifest = read(args.run / 'manifest.json')
    phase_statuses = {}
    for phase in ('search', 'analyze'):
        statuses = [read(args.run / f'status/{phase}/{slot}.json') for slot in range(manifest['tasks'])]
        assert all('finished' in row for row in statuses), f'{phase} still running'
        phase_statuses[phase] = dict(
            tasks=len(statuses), outcomes=dict(Counter(str(v.get('exit', v.get('status')))
                for row in statuses for v in row['variants'].values())),
            execution_span=max(s['finished'] for s in statuses)-min(s['started'] for s in statuses),
            exceptions=[dict(slot=s['slot'], index=s['index'], direction=s['direction'], variant=k, **v)
                        for s in statuses for k, v in s['variants'].items() if v.get('exit') != 0])
    subprocess.run([PYTHON, str(args.run / 'driver.py'), 'collect', '--run', str(args.run)],
                   check=True, stdout=subprocess.DEVNULL)
    rows = read(args.run / 'comparison/case_metrics.json')
    summary = read(args.run / 'comparison/summary.json')
    summary['phase_statuses'] = phase_statuses
    summary['directional_recovery'] = {key: {direction: dict(Counter(
        r['variants'][key]['directions'][direction]['recovery'] for r in rows))
        for direction in ('R_to_P', 'P_to_R')} for key in TOLERANCES}
    paired = [r for r in rows if all(r['variants'][key]['searches_complete'] for key in TOLERANCES)]
    summary['paired_case_cpu_median'] = {k: statistics.median(r['variants'][k]['mapping_cpu'] for r in paired)
                                        for k in TOLERANCES}
    summary['faster_cases_at_1p5'] = sum(r['variants']['tol_1p5']['mapping_cpu'] < r['variants']['tol_1p0']['mapping_cpu'] for r in paired)
    summary['representative_and_graph_counts'] = {key: {
        metric: sum(r['variants'][key]['directions'][direction][metric] or 0
                    for r in rows for direction in ('R_to_P', 'P_to_R'))
        for metric in ('states', 'terminals')} for key in TOLERANCES}
    # Recheck explicit saved witnesses against the fixed chemical reference
    # graph, independently of the evaluator's reported recovery flag.
    import pynauty
    from golden_evaluation import colored_graph, project
    witness_checks = Counter()
    weight_case_counts = Counter()
    added_weight_pairs = []
    changes = {r['index'] for r in summary['recovery_changes']}
    changed_evaluations = []
    for index in range(manifest['cases']):
        input_dir = args.run / f'runs/tol_1p0/inputs/golden/{index}'
        raw = read(input_dir / 'input.json')
        weights = [sorted({w for row in raw[side]['wbo'] for w in row if w >= .2})
                   for side in ('reactant', 'product')]
        weight_case_counts.update(str(w) for w in set(weights[0] + weights[1]))
        added = [(a, b) for a in weights[0] for b in weights[1] if 1.0 < abs(a-b) <= 1.5]
        if added:
            added_weight_pairs.append(dict(index=index, added_weight_pairs=added))
        reference = read(input_dir / 'reference.json')
        features = reference['features']
        expected = pynauty.certificate(colored_graph(features, project(reference['mapping'], features)))
        for key in TOLERANCES:
            for direction in ('R_to_P', 'P_to_R'):
                folder = result_folder(args.run / 'runs' / key,
                                       dict(dataset='golden', index=index, direction=direction))
                path = folder / 'full_sweep_evaluation.json'
                if not path.exists():
                    continue
                result = read(path)
                if index in changes:
                    changed_evaluations.append(dict(index=index, variant=key, direction=direction, evaluation=result))
                if result['reference_recovery'] != 'recovered':
                    continue
                mapping = dict(result['input_orientation_witness'])
                assert len(set(mapping.values())) == len(mapping), (index, key, direction, 'injectivity')
                assert all(raw['reactant']['elements'][r] == raw['product']['elements'][p]
                           for r, p in mapping.items()), (index, key, direction, 'elements')
                assert pynauty.certificate(colored_graph(features, project(mapping, features))) == expected, (index, key, direction, 'reference')
                witness_checks[key] += 1
    summary['witness_validation'] = dict(checked_directional_witnesses=dict(witness_checks), violations=0,
        scope='Every recovered directional output: saved input-oriented witness is injective, element preserving, '
              'and matches the fixed heavy chemical reference modulo endpoint chemical symmetry. '
              'This does not exhaustively audit every non-witness compressed family member.')
    extended = None
    if (args.run / 'extended_1599.json').exists():
        extended = read(args.run / 'extended_1599.json')
        assert extended['index'] == 1599 and extended['direction'] == 'P_to_R'
        assert extended['variant'] == 'tol_1p5' and not extended['search_rerun']
        assert sha(Path(extended['archive'])) == extended['archive_sha256']
        outcomes = []
        for row in rows:
            variant = row['variants']['tol_1p5']
            if row['index'] == 1599:
                statuses = [variant['directions']['R_to_P']['recovery'], extended['evaluation']['reference_recovery']]
                outcome = ('recovered' if 'recovered' in statuses else
                           'not_recovered' if all(s == 'not_recovered' for s in statuses) else 'unknown')
            else:
                outcome = variant['recovery']
            outcomes.append(outcome)
        summary['targeted_extended_verification'] = dict(index=1599, direction='P_to_R',
            seconds=extended['verification_seconds'], search_rerun=False,
            outcome=extended['evaluation']['reference_recovery'],
            updated_trial_recovery=dict(Counter(outcomes)),
            scope='Additional verification only; primary fixed-budget metrics above are unchanged.')
    # Detect source or input changes during the campaign.
    for category, directory in [('original_sha256', 'original'), ('adapter_sha256', 'engine/bench')]:
        for relative, expected in manifest[category].items():
            assert sha(args.run / directory / relative) == expected, relative
    for index, files in read(args.run / 'input_hashes.json').items():
        for name, expected in files.items():
            assert sha(args.run / f'runs/tol_1p0/inputs/golden/{index}/{name}') == expected, (index, name)
    summary['frozen_sources_and_inputs_unchanged'] = True
    save(args.run / 'comparison/summary.json', summary)
    args.destination.mkdir(parents=True, exist_ok=False)
    for name in ('manifest.json', 'input_hashes.json', 'submissions.json'):
        shutil.copy2(args.run / name, args.destination / name)
    for name in ('summary.json', 'case_metrics.json', 'per_case.csv'):
        if name.endswith('.csv'):
            (args.destination / name).write_text((args.run / 'comparison' / name).read_text())
        else:
            shutil.copy2(args.run / 'comparison' / name, args.destination / name)
    save(args.destination / 'changed_case_evaluations.json', changed_evaluations)
    if extended:
        for name in ('extended_1599.json', 'extended_1599_submission.json'):
            shutil.copy2(args.run / name, args.destination / name)
    save(args.destination / 'input_tolerance_effects.json', dict(
        endpoint_positive_weight_case_counts=dict(weight_case_counts),
        cases_with_new_cross_endpoint_weight_pairs=added_weight_pairs,
        scope='Input bond-weight gate only; ignores atom labels and subgraph context. '
              'This does not predict whether an added edge assignment survives full search.'))
    save(args.destination / 'artifacts.json', dict(raw_run=str(args.run), frozen_driver=str(args.run / 'driver.py'),
        publisher_sha256=sha(Path(__file__)), baseline=manifest['baseline']))
    jobs = ','.join(r['job'].split(';')[0] for r in read(args.run / 'submissions.json'))
    accounting = subprocess.check_output(['sacct', '-j', jobs, '--array', '--noheader', '--parsable2',
        '--format=JobID,State,ElapsedRaw,TotalCPU,AllocCPUS,Submit,Start,End,MaxRSS'], text=True)
    (args.destination / 'slurm_accounting.psv').write_text(accounting)
    first, second = [summary['metrics'][k] for k in TOLERANCES]
    def count(metric, row):
        return row['recovery'].get(metric, 0)
    lines = ['# Golden: matching tolerance 1.0 versus 1.5', '',
        'Fresh paired measurements on all 1,851 Golden cases. Frozen engine `98b01b1`, one seed, '
        'branch cap 100, explicit H, eight workers per directional call, original uncut/single-edge sweep, '
        'both directions. Only `iso_tolerance` changes; this affects fragment compatibility and symmetry grouping. '
        'Reference scoring, event thresholds, graph floor, seed orders, and branch cap stay fixed.', '',
        '| Metric | Tolerance 1.0 | Tolerance 1.5 |', '|---|---:|---:|',
        f"| Recovered reference | {count('recovered',first)}/1851 ({100*count('recovered',first)/1851:.2f}%) | {count('recovered',second)}/1851 ({100*count('recovered',second)/1851:.2f}%) |",
        f"| Not recovered from saved output | {count('not_recovered',first)} | {count('not_recovered',second)} |",
        f"| Unknown | {count('unknown',first)} | {count('unknown',second)} |",
        f"| Completed bidirectional searches | {first['completed_search_cases']} | {second['completed_search_cases']} |",
        f"| Mapping CPU on paired completed cases, seconds | {summary['paired_mapping_cpu']['tol_1p0']:.3f} | {summary['paired_mapping_cpu']['tol_1p5']:.3f} |",
        f"| Median paired case mapping CPU, seconds | {summary['paired_case_cpu_median']['tol_1p0']:.3f} | {summary['paired_case_cpu_median']['tol_1p5']:.3f} |",
        f"| Cases with branch-cap flags | {first['capped_cases']} | {second['capped_cases']} |", '',
        f"Paired completed cases: {len(paired)}. CPU ratio (1.5/1.0): {summary['cpu_ratio_1p5_over_1p0']:.4f}. "
        'CPU is parent plus child compute excluding measured archive persistence/loading; it already includes all eight workers. '
        'These are single paired measurements with alternating run order on the same host and CPU set. '
        'Offline evaluation and initial scheduler queue time are separate.', '',
        'Golden input weights are discrete bond orders (1, 1.5, 2, 3). The added edge-weight '
        'compatibility is aromatic/triple (1.5 versus 3); both occur across endpoints in '
        f'{len(added_weight_pairs)} cases. This input-level count ignores element and subgraph constraints. '
        'It differs from the continuous WBO setting of the TS holdout.', '',
        '## Case changes', '', '| Golden index (zero-based) | At 1.0 | At 1.5 |', '|---|---|---|']
    lines += [f"| {c['index']} | {c['baseline']} | {c['trial']} |" for c in summary['recovery_changes']]
    if not summary['recovery_changes']:
        lines += ['| None | | |']
    if extended:
        lines += ['', '### Targeted verification of case 1599', '',
            f"Only the saved tolerance-1.5 P-to-R archive was rechecked with a {extended['verification_seconds']:g}-second "
            f"verification budget. It returned `{extended['evaluation']['reference_recovery']}` after "
            f"{extended['evaluation']['evaluation_seconds']:.2f} seconds. Mapping was not rerun. "
            f"After this extra check, trial recovery counts are {summary['targeted_extended_verification']['updated_trial_recovery']}. "
            'The primary table retains the original equal-budget results.']
    lines += ['', 'A not-recovered result describes the saved search families under this cap; it is not a proof '
        'that the algorithm cannot find the reference. Evaluator timeouts remain unknown. '
        'Coverage is recovery anywhere in the saved output, not top-ranked mapping accuracy.', '',
        f"The fresh 1.0 control differs from the archived 1.0 recovery status on indices: {summary['baseline_reproduction_changes']}. "
        'Both the historical control and the fresh paired control are retained in `case_metrics.json`.', '',
        f"All reported recovered witnesses passed element, injectivity, and fixed-reference certificate checks: "
        f"{dict(witness_checks)} directional witnesses. Frozen search/adapters and input/reference hashes were unchanged. "
        'This audit covers recovered witnesses; it does not enumerate every compressed family member.', '',
        'Directional counts, caps, phase exceptions, and all changed-case witnesses are in the accompanying JSON files. '
        'The existing paper baseline and engine defaults are unchanged.', '',
        '## Reproduction', '', 'From the paper checkout:', '', '```bash',
        f'{PYTHON} bench/golden_tolerance_benchmark.py prepare --run /path/to/new_run',
        f'{PYTHON} bench/golden_tolerance_benchmark.py submit --run /path/to/new_run',
        '# After both phases finish:',
        f'PYTHONPATH=/path/to/new_run/original/src:/path/to/new_run/engine/bench {PYTHON} bench/publish_golden_tolerance.py --run /path/to/new_run --destination /path/to/new_report',
        '```', '', f'Raw archives and logs: `{args.run}`.', '']
    (args.destination / 'README.md').write_text('\n'.join(lines))
    print(json.dumps(summary, indent=2))
    print(args.destination)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    publish(parser.parse_args())

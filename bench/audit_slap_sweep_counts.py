"""Re-score saved SLAP sweep journals without running or changing any mapper.

Uses the campaign's frozen strict evaluator so the metric stays comparable.
Independently reconstructs every case/mode union and validates journal coverage,
task status, saved hit lists, and one native witness frame per recovered case.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import time


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def initialize(run):
    global RUN, TASKS, AUDIT, SIGNATURES
    RUN = Path(run)
    sys.path.insert(0, str(RUN))
    from golden_competitors import signatures
    SIGNATURES = signatures
    TASKS = read(RUN / 'tasks.json')
    manifest = read(RUN / 'manifest.json')
    AUDIT = [json.loads(line) for line in
             (Path(manifest['dataset']) / 'audit.jsonl').read_text().splitlines()]


def case(index):
    expected_endpoints, expected = SIGNATURES(AUDIT[index]['mapped_reaction'])
    previous = read(RUN / f'evaluations/{index}.json')
    modes, directions, witnesses = set(), set(), []
    calls = errors = invalid = candidates = 0
    uncut = cut = False
    mapping_cpu = workflow_cpu = 0.0
    incomplete, variants, journal_hashes = [], [], []
    first_frame_checked = False
    tasks = TASKS[4*index:4*index+4]
    assert {(t['direction'], t['mode']) for t in tasks} == {
        (direction, mode) for direction in ('R_to_P', 'P_to_R') for mode in ('binary', 'weighted')}
    for task, old in zip(tasks, previous['attempts']):
        assert task['case'] == index and task['slot'] == old['slot']
        folder = RUN / f"outputs/{index}/{task['direction']}/{task['mode']}"
        payload = (folder / 'records.jsonl').read_bytes()
        lines = payload.splitlines(keepends=True)
        rows = [json.loads(line) for line in lines if line.endswith(b'\n')]
        assert len(rows) <= task['expected_cuts']
        original = read(folder / 'input.json')
        expected_cuts = [None] + sorted([sorted(edge[:2]) for edge in original['input_graphs'][0]['edges']])
        assert len(expected_cuts) == task['expected_cuts']
        hits, next_offset = [], 0
        archive_size = (folder / 'native.bin').stat().st_size
        for ordinal, row in enumerate(rows):
            assert row['ordinal'] == ordinal and row['cut'] == expected_cuts[ordinal]
            frame = row['native']
            assert frame['offset'] == next_offset
            next_offset += frame['length']
            assert next_offset <= archive_size
            calls += 1
            errors += row['status'] != 'mapped'
            mapping_cpu += row['mapping_cpu']
            workflow_cpu += row['mapping_cpu'] + row['graph_cpu'] + row['export_cpu']
            for number, prediction in enumerate(row['candidates']):
                candidates += 1
                try:
                    endpoints, actual = SIGNATURES(prediction['mapped_rxn'])
                    if endpoints != expected_endpoints:
                        raise ValueError('Endpoint chemistry differs from reference')
                except Exception:
                    invalid += 1
                    continue
                if actual != expected:
                    continue
                hits.append(dict(cut=row['cut'], candidate=number, ordinal=ordinal))
                modes.add(task['mode'])
                directions.add(task['direction'])
                uncut |= row['cut'] is None
                cut |= row['cut'] is not None
                if not first_frame_checked:
                    with (folder / 'native.bin').open('rb') as f:
                        f.seek(frame['offset'])
                        blob = f.read(frame['length'])
                    assert hashlib.sha256(blob).hexdigest() == frame['sha256']
                    first_frame_checked = True
                    witnesses.append(dict(slot=task['slot'], direction=task['direction'], mode=task['mode'],
                                          ordinal=ordinal, candidate=number, cut=row['cut'], **frame))
        assert hits == old['hits'], f'Hit-list disagreement: case {index}, slot {task["slot"]}'
        assert len(rows) == old['completed_cuts']
        status = read(RUN / f"status/{task['slot']}.json")
        assert status['completed_cuts'] == len(rows)
        assert status['expected_cuts'] == task['expected_cuts']
        if status['exit'] == 0:
            assert read(folder / 'complete.json')['cuts'] == len(rows) == task['expected_cuts']
        detail = dict(slot=task['slot'], case=index, direction=task['direction'], mode=task['mode'],
                      completed=len(rows), expected=task['expected_cuts'], exit=status['exit'],
                      elapsed_seconds=status['elapsed_including_startup_io'],
                      next_ordinal=len(rows) if len(rows) < len(expected_cuts) else None,
                      next_cut=expected_cuts[len(rows)] if len(rows) < len(expected_cuts) else None,
                      hits=len(hits), uncommitted_tail_bytes=sum(len(x) for x in lines if not x.endswith(b'\n')))
        if len(rows) != task['expected_cuts']:
            incomplete.append(detail)
        variants.append(detail)
        journal_hashes.append(dict(slot=task['slot'], sha256=hashlib.sha256(payload).hexdigest()))
    recovered = bool(modes)
    for key, value in dict(recovered=recovered, uncut_recovered=uncut, cut_recovered=cut,
                           completed_calls=calls, mapping_errors=errors, invalid_predictions=invalid).items():
        assert previous[key] == value, (index, key, previous[key], value)
    assert math.isclose(previous['mapping_cpu'], mapping_cpu, abs_tol=1e-6)
    assert math.isclose(previous['workflow_cpu_excluding_io'], workflow_cpu, abs_tol=1e-6)
    return dict(index=index, recovered=recovered, binary='binary' in modes, weighted='weighted' in modes,
                r_to_p='R_to_P' in directions, p_to_r='P_to_R' in directions, uncut=uncut, cut_only=cut,
                completed_calls=calls, candidate_outputs=candidates, mapping_errors=errors,
                invalid_predictions=invalid, mapping_cpu=mapping_cpu, workflow_cpu=workflow_cpu,
                fully_swept_valid=not incomplete and not errors and not invalid,
                incomplete=incomplete, variants=variants, witnesses=witnesses, journals=journal_hashes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = read(args.run / 'manifest.json')
    audit_file = Path(manifest['dataset']) / 'audit.jsonl'
    assert digest(audit_file) == manifest['audit_sha256']
    tasks = read(args.run / 'tasks.json')
    assert [t['slot'] for t in tasks] == list(range(manifest['tasks']))
    assert len(tasks) == 4 * manifest['cases']
    assert sum(t['expected_cuts'] for t in tasks) == manifest['expected_calls']
    assert sorted(int(p.stem) for p in (args.run / 'evaluations').glob('*.json')) == list(range(manifest['cases']))
    started = time.monotonic()
    results = []
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize, initargs=(str(args.run),)) as pool:
        for result in pool.map(case, range(manifest['cases'])):
            results.append(result)
            if len(results) % 100 == 0:
                print(f"Validated {len(results)}/{manifest['cases']} cases in {time.monotonic()-started:.1f}s", flush=True)
    successes = {r['index'] for r in results if r['recovered']}
    missing = set(range(manifest['cases'])) - successes
    incomplete = [v for r in results for v in r['incomplete']]
    summary = dict(status='passed', run=str(args.run), baseline_upstream_commit=manifest['upstream_commit'],
        cases=len(results), recovered=len(successes), uncut_recovered=sum(r['uncut'] for r in results),
        binary_recovered=sum(r['binary'] for r in results), weighted_recovered=sum(r['weighted'] for r in results),
        both_modes_recovered=sum(r['binary'] and r['weighted'] for r in results),
        cut_only_recovered=sum(r['cut_only'] for r in results),
        completed_misses=[r['index'] for r in results if not r['recovered'] and r['fully_swept_valid']],
        unresolved_misses=[r['index'] for r in results if not r['recovered'] and not r['fully_swept_valid']],
        incomplete_variants=incomplete, task_statuses=dict(Counter(str(v['exit']) for r in results for v in r['variants'])),
        fully_swept_valid_cases=sum(r['fully_swept_valid'] for r in results),
        completed_calls=sum(r['completed_calls'] for r in results), expected_calls=manifest['expected_calls'],
        candidate_outputs_rechecked=sum(r['candidate_outputs'] for r in results),
        mapping_errors=sum(r['mapping_errors'] for r in results), invalid_predictions=sum(r['invalid_predictions'] for r in results),
        mapping_cpu=sum(r['mapping_cpu'] for r in results), workflow_cpu=sum(r['workflow_cpu'] for r in results),
        native_witness_frames_verified=sum(len(r['witnesses']) for r in results),
        recovered_case_indices=sorted(successes), elapsed_audit_wall_seconds=time.monotonic()-started,
        original_summary_sha256=digest(args.run / 'summary.json'),
        dataset_sha256=digest(audit_file),
        evaluator_sha256={name:digest(args.run / name) for name in ['golden_competitors.py','golden_evaluation.py']},
        audit_script_sha256=digest(Path(__file__)),
        metric='Frozen strict heavy-atom relation modulo endpoint chemical symmetry; all native exported alternatives, not top-1.',
        scope='Fresh re-scoring of every committed candidate with the original evaluator; no mapping rerun or watchdog-budget change. One native frame per recovered case hash-checked, not every native frame decoded. Recorded CPU omits interrupted in-flight work.')
    original = read(args.run / 'summary.json')
    assert summary['recovered'] == original['sweep_union_recovered']
    assert summary['binary_recovered'] == original['sweep_recovery_by_mode']['binary']
    assert summary['weighted_recovered'] == original['sweep_recovery_by_mode']['weighted']
    assert missing == set(original['unrecovered'])
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    (args.output / 'case_checks.json').write_text(json.dumps(results, separators=(',',':'))+'\n')
    fields = ['index','recovered','binary','weighted','r_to_p','p_to_r','uncut','cut_only','completed_calls',
              'candidate_outputs','mapping_errors','invalid_predictions','fully_swept_valid']
    with (args.output / 'cases.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader();writer.writerows(results)
    print(json.dumps({k:v for k,v in summary.items() if not isinstance(v,(list,dict))}, indent=2), flush=True)


if __name__ == '__main__':
    main()

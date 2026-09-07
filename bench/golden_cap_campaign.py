"""Isolate branch-cap effects on previously investigated Golden cases."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from investigate_golden_mapping import save


def initialize(run, retry_from=None, slot=None):
    run.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    base = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
    cases = [11, 602, 745, 853, 986, 1285, 1405, 1654, 1665, 1740]
    jobs = [dict(index=i, seeds=3) for i in cases]
    jobs += [dict(index=i, seeds=100) for i in (1665, 1740)]
    for job in jobs:
        job.update(cap=2000, tolerance=1.0, workers=8,
                   source=str(base/'golden_full_20260906'/str(job['index'])))
    if retry_from is not None:
        jobs = [json.loads((retry_from/'manifest.json').read_text())['jobs'][slot]]
    save(run/'manifest.json', dict(jobs=jobs, watchdog_seconds=570,
         purpose='Cap-only comparison: unchanged seeds, tolerance 1.0, ordinary sweep cuts; no reference-directed search',
         selection='Five previously empty cases; four capped cap1000 misses; unresolved blind case1740. Also repeat 1665/1740 at their previous 100 seeds.',
         git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()))
    if retry_from is not None:
        save(run/'retry.json', dict(original_run=str(retry_from), slot=slot,
             reason='32 GiB cgroup killed workers; retry with 128 GiB, same search configuration'))
    for directory in ('src', 'bench'):
        shutil.copytree(root/directory, run/'engine'/directory,
                        ignore=shutil.ignore_patterns('__pycache__'))


def worker(run, slot):
    manifest = json.loads((run/'manifest.json').read_text())
    job = manifest['jobs'][slot]
    output = run/f"case{job['index']}_seeds{job['seeds']}_cap{job['cap']}"
    command = [sys.executable, str(Path(__file__).with_name('golden_mapping_probe.py')),
               '--source', job['source'], '--output', str(output),
               '--seeds', str(job['seeds']), '--cap', str(job['cap']),
               '--tolerance', str(job['tolerance']), '--workers', str(job['workers'])]
    started = time.perf_counter()
    cgroup = Path('/sys/fs/cgroup') / Path('/proc/self/cgroup').read_text().strip().split('::')[-1].lstrip('/')
    events = next((p/'memory.events' for p in (cgroup, *cgroup.parents)
                   if p.name.startswith('job_')), None)
    def oom_kills():
        return int(dict(line.split() for line in events.read_text().splitlines()).get('oom_kill', 0)) if events else 0
    initial_oom = oom_kills()
    with (run/f'slot{slot}.log').open('w') as log:
        child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        code = None
        while code is None:
            try:
                code = child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                if oom_kills() > initial_oom:
                    code = 'oom_guard'
                elif time.perf_counter()-started >= manifest['watchdog_seconds']:
                    code = 'watchdog'
        if code in ('oom_guard', 'watchdog'):
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
    save(run/f'slot{slot}_status.json', dict(job=job, exit_code=code,
         seconds=time.perf_counter()-started, slurm_job=os.environ.get('SLURM_JOB_ID'),
         hostname=os.uname().nodename))


def report(run, output):
    manifest = json.loads((run/'manifest.json').read_text())
    rows = []
    for slot, job in enumerate(manifest['jobs']):
        directory = run/f"case{job['index']}_seeds{job['seeds']}_cap{job['cap']}"
        evaluation_path = directory/'evaluation.json'
        evaluation = json.loads(evaluation_path.read_text()) if evaluation_path.exists() else {}
        search_path = directory/'search.json'
        search = json.loads(search_path.read_text()) if search_path.exists() else {}
        status_path = run/f'slot{slot}_status.json'
        status = json.loads(status_path.read_text()) if status_path.exists() else {}
        partial_path = directory/'partial_reference_check.json'
        partial = json.loads(partial_path.read_text()) if partial_path.exists() else {}
        rows.append(dict(index=job['index'], seeds=job['seeds'], cap=job['cap'],
            reference_recovery=evaluation.get('reference_recovery',
                'recovered' if partial.get('recovered') else
                'not_evaluated' if status.get('exit_code') else 'unfinished'),
            scope='full_sweep' if evaluation else 'partial_completed_cuts',
            capped=evaluation.get('capped'), terminals=evaluation.get('candidate_terminals'),
            cap_stops=search.get('metrics', {}).get('subtree_branch_cap_count'),
            search_seconds=search.get('seconds'), evaluation_seconds=evaluation.get('evaluation_seconds'),
            total_seconds=status.get('seconds'), exit_code=status.get('exit_code'),
            archive=str(directory/'cuts/aam.pkl.gz') if (directory/'cuts/aam.pkl.gz').exists() else str(directory/'cuts'),
            evaluation=str(evaluation_path) if evaluation_path.exists() else str(partial_path) if partial_path.exists() else ''))
    output.mkdir(parents=True, exist_ok=True)
    save(output/'results.json', dict(manifest=manifest, results=rows))
    with (output/'results.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(rows, indent=2))


def check_partial(run, slot):
    """Positive-only certificate check of completed cuts; not a full evaluation."""
    import gc
    import pynauty
    from golden_evaluation import colored_graph, project
    from rxn_core.search_graph import AAMSearchGraph
    job = json.loads((run/'manifest.json').read_text())['jobs'][slot]
    directory = run/f"case{job['index']}_seeds{job['seeds']}_cap{job['cap']}"
    reference = json.loads((Path(job['source'])/'reference.json').read_text())
    features = reference['features']
    expected = pynauty.certificate(colored_graph(features, project(reference['mapping'], features)))
    rows = []
    gc.disable()
    for path in sorted((directory/'cuts').glob('cut_*.json')):
        graph = AAMSearchGraph.from_record(json.loads(path.read_bytes()), copy=False)
        seen = set()
        hit = None
        for terminal in graph.terminals:
            mapping = project(graph.states[terminal].mapping, features)
            key = tuple(sorted(mapping.items()))
            if key in seen:
                continue
            seen.add(key)
            if pynauty.certificate(colored_graph(features, mapping)) == expected:
                hit = dict(terminal=terminal, mapping=graph.states[terminal].mapping)
                break
        rows.append(dict(cut=str(path), terminals=len(graph.terminals),
                         representative_witness=hit))
        save(directory/'partial_reference_check.json', dict(
            scope='Completed cut representatives only. Absence is inconclusive; no full symmetry-family verification.',
            recovered=any(r['representative_witness'] is not None for r in rows), cuts=rows))
        del graph


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['init', 'worker', 'report', 'partial'])
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--slot', type=int)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--retry-from', type=Path)
    args = parser.parse_args()
    if args.mode == 'init':
        initialize(args.run, args.retry_from, args.slot)
    elif args.mode == 'worker':
        worker(args.run, args.slot)
    elif args.mode == 'report':
        report(args.run, args.output)
    else:
        check_partial(args.run, args.slot)

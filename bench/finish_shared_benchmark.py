"""Finish interrupted benchmark calls from their verified cut checkpoints.

No alternate search policy is used. Completed cuts are never rematched.
Continuation timing is explicitly distinct from a fresh-run measurement.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

from adaptive_full_benchmark import read, save, spec_and_folder, problem_plan, analyze
from publication_timing import SearchProfiler


def prepare(args):
    tasks = read(args.run/'tasks.json')
    path = args.run/'continuation_tasks.json'
    pending = read(path) if path.exists() else []
    for slot, spec in enumerate(tasks):
        path = args.run/f'status/adaptive_{slot}.json'
        status = read(path) if path.exists() else {}
        if not status.get('complete'):
            continue  # active/queued work belongs to the original allocation
        if slot not in pending and (status.get('search',{}).get('exit') != 0 or status.get('analyze',{}).get('exit') != 0):
            pending.append(slot)
    save(args.run/'continuation_tasks.json', pending)
    print(json.dumps(dict(pending=pending, count=len(pending))))


def finish(args):
    from rxn_core import search_aam
    from rxn_core.artifacts import read_aam_checkpoint
    from dataclasses import asdict
    args.slot = read(args.run/'continuation_tasks.json')[args.ordinal]
    args.method = 'adaptive'
    spec, folder = spec_and_folder(args)
    status_path = args.run/f'status/adaptive_{args.slot}.json'
    previous = read(status_path)
    history = args.run/f'continuations/{args.slot}'
    history.mkdir(parents=True,exist_ok=False)
    save(history/'previous_status.json',previous)
    manifest = read(args.run/'manifest.json')
    raw, plan = problem_plan(args,spec)
    collector_commit = subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    record = dict(previous, continuation=dict(collector_commit=collector_commit,started=time.time(),
        source='Existing raw/finalized cut checkpoints; same search config and execution backend'))
    started = time.perf_counter()
    if previous.get('search',{}).get('exit') != 0:
        if (folder/'search.json').exists():
            save(history/'previous_search.json',read(folder/'search.json'))
        with SearchProfiler(history/'timing_events') as profiler:
            if (folder/'cuts/aam.pkl.gz').exists():
                result = read_aam_checkpoint(folder/'cuts/aam.pkl.gz')
            else:
                result = search_aam(plan.problem,plan.config,workers=manifest['adaptive_workers'],
                    execution='shared_policies',intermediate_dir=folder/'cuts',archive_format='checkpoint',resume=True)
        completed = []
        for path in (folder/'timing_events').glob('*.jsonl'):
            completed.extend(json.loads(line) for line in path.read_text().splitlines())
        previous_completed_cpu = sum(r['cpu_seconds'] for r in completed
            if not r['parent'] and r['phase'] in ('cut_search_inclusive','symmetry'))
        timing = profiler.summary()
        timing['compute_cpu_excluding_persistence_and_loading_seconds'] += previous_completed_cpu
        row = dict(label='full_sweep',archive='cuts/aam.pkl.gz',**timing,
            metrics=asdict(result.metrics),states=len(result.graph.states),terminals=len(result.graph.terminals),
            capped=result.graph.capped,timing_source='checkpoint_continuation',
            previous_completed_worker_phase_cpu=previous_completed_cpu,
            timing_scope='Completed earlier cut/group phases plus continuation CPU. Earlier setup/IPC and '
                'abandoned parent phases are not reconstructible; exclude this row from fresh paired timing.')
        save(folder/'search.json',dict(**spec,method='adaptive',rows=[row],complete=True,
                                       collector_commit=collector_commit))
        record['search'] = dict(exit=0,elapsed_including_io=time.perf_counter()-started,continued=True)
        save(status_path,record)
        del result
    started = time.perf_counter()
    analyze(args)
    record['analyze'] = dict(exit=0,elapsed_including_io=time.perf_counter()-started,continued=True)
    record['continuation']['finished'] = time.time()
    record['complete'] = True
    save(status_path,record)
    save(history/'finished.json',record)
    print(json.dumps(dict(slot=args.slot,**spec,complete=True)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['prepare','finish'])
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--ordinal',type=int)
    args=p.parse_args()
    globals()[args.command](args)

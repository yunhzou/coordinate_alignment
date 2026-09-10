"""Full Golden + XYZ/WBO holdout confirmation against a frozen original engine.

Search never reads reference labels. Both directions, full checkpoints, fixed
policies, and failed attempts are retained. Analysis runs in a separate process
after persistence; a failed evaluation cannot destroy a successful search.
"""
import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import errno
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
DIRECTIONS = ('R_to_P', 'P_to_R')
SOURCES = dict(golden=DATA/'golden_publication_20260908',
               holdout=DATA/'elementary140_tol1_20260909')
COUNTS = dict(golden=1851, holdout=140)


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n')
    temporary.replace(path)


def read(path):
    # Live status files are atomically replaced on NFS. A reader can briefly
    # receive ESTALE for the replaced inode; reopen that same pathname only.
    for attempt in range(3):
        try:
            return json.loads(path.read_text())
        except OSError as error:
            if error.errno != errno.ESTALE or attempt == 2:
                raise
            time.sleep(.01)


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def environment():
    return dict(host=os.uname().nodename, python=sys.version,
        affinity=sorted(os.sched_getaffinity(0)),
        cpu_models=sorted({line.split(':', 1)[1].strip() for line in
            Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')}),
        slurm={k:os.environ.get(k) for k in ('SLURM_JOB_ID', 'SLURM_CPUS_PER_TASK',
                                         'SLURM_ARRAY_TASK_ID', 'SLURM_JOB_PARTITION')})


def prepare(args):
    from rxn_core import AAMSearchConfig
    args.run.mkdir(parents=True, exist_ok=False)
    for folder in ('src', 'native', 'bench', 'tests', 'tools', 'benchmarks', 'docs/example_runs'):
        shutil.copytree(ROOT/folder, args.run/'engine'/folder,
                        ignore=shutil.ignore_patterns('__pycache__'))
    original = args.run/'original'
    original.mkdir()
    commit = subprocess.check_output(['git', 'rev-parse', args.original_commit], text=True).strip()
    archive = subprocess.check_output(['git', 'archive', commit, 'src', 'native'])
    subprocess.run(['tar', '-x', '-C', str(original)], input=archive, check=True)
    # Build the actual frozen baseline, not a current binary loaded into old
    # Python or an approximately matching historical build.
    subprocess.run(['timeout', '--kill-after=5s', '300', sys.executable,
                    str(original/'native/build_engine.py')], check=True)
    tasks, hashes = [], []
    for dataset, count in COUNTS.items():
        for index in range(count):
            destination = args.run/f'inputs/{dataset}/{index}'
            destination.mkdir(parents=True)
            for name in ('input.json', 'reference.json') if dataset == 'golden' else ('input.json',):
                shutil.copy2(SOURCES[dataset]/f'inputs/{index}/{name}', destination/name)
                hashes.append(dict(dataset=dataset, index=index, file=name, sha256=sha(destination/name)))
            raw = read(destination/'input.json')
            small = 'P_to_R' if len(raw['reactant']['elements']) > len(raw['product']['elements']) else 'R_to_P'
            tasks.extend(dict(dataset=dataset, index=index, direction=d, smaller_first=small) for d in DIRECTIONS)
    # Interleave datasets so the entire holdout is not stuck behind Golden.
    tasks.sort(key=lambda t:(t['index'], t['dataset'], t['direction']))
    save(args.run/'tasks.json', tasks)
    save(args.run/'input_hashes.json', hashes)
    save(args.run/'manifest.json', dict(schema='adaptive_full_benchmark/v1', counts=COUNTS,
        sources={k:str(v) for k,v in SOURCES.items()}, tasks=len(tasks), original_commit=commit,
        experimental_commit=subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
        original_config=asdict(AAMSearchConfig(seed_count=10, branch_limit=100, iso_tolerance=1.)),
        adaptive_config=asdict(AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1.)),
        original_execution='reused_native', original_workers=args.workers,
        adaptive_work_budgets=args.work_budgets, adaptive_cpu_seconds=args.cpu_seconds,
        root_seed=42, python_hash_seed=0, explicit_H=True,
        watchdog=dict(search_seconds=300, analysis_seconds=240, kill_after_seconds=5, slurm_minutes=10),
        branch_cap_scope=dict(original='Native growth and synchronized live frontier',
                              adaptive='Native growth only; agenda has declared work/CPU budget'),
        modes=dict(single='smaller explicit endpoint to larger; R_to_P on ties',
                   bidirectional='union of directional outputs, no splicing or reference selection'),
        metrics='Unchanged Golden heavy-reference chemical symmetry and partial-annotation evaluator; '
                'full explicit-H event scoring at 0.5. Holdout has NO annotated ground truth: '
                'compare full mapping feasibility, event scores and original score-response heavy classes. '
                'Absent representative classes are NOT claimed absent from compressed families.',
        timing='CPU excludes measured checkpoint persistence/loading. Search and analysis separate. '
               'Parallel elapsed includes IO, excludes queue; never subtract summed worker IO from elapsed. '
               'Original uses multiple CPUs, adaptive one: compare CPU costs, not equal-resource latency.',
        engine_sha256={str(p.relative_to(args.run)):sha(p) for directory in ('engine','original')
                       for p in (args.run/directory).rglob('*') if p.is_file()
                       and 'build' not in p.relative_to(args.run).parts}))
    (args.run/'status').mkdir()
    print(json.dumps(dict(run=str(args.run), directional_pairs=len(tasks), search_calls=2*len(tasks))), flush=True)


def prepare_revision(args):
    """Freeze a fresh candidate; retain the already measured original baseline."""
    reference = args.reference_run.resolve()
    previous = read(reference/'manifest.json')
    assert previous['counts'] == COUNTS
    args.run.mkdir(parents=True, exist_ok=False)
    for name in ('src', 'native', 'bench', 'tests', 'tools', 'benchmarks', 'docs/example_runs'):
        shutil.copytree(ROOT/name, args.run/'engine'/name,
                        ignore=shutil.ignore_patterns('__pycache__'))
    for name in ('inputs', 'original'):
        (args.run/name).symlink_to(reference/name, target_is_directory=True)
    for name in ('tasks.json', 'input_hashes.json'):
        shutil.copy2(reference/name, args.run/name)
    tasks = read(args.run/'tasks.json')
    for task in tasks:
        folder = Path(f"results/{task['dataset']}/{task['index']}/{task['direction']}")
        (args.run/folder).mkdir(parents=True)
        (args.run/folder/'original').symlink_to(reference/folder/'original', target_is_directory=True)
    manifest = dict(previous, reference_run=str(reference), methods=['adaptive'],
        adaptive_policy=args.adaptive_policy, adaptive_work_budgets=args.work_budgets,
        adaptive_cpu_seconds=args.cpu_seconds,
        experimental_commit=subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
        engine_sha256={str(p.relative_to(args.run/'engine')):sha(p)
            for p in (args.run/'engine').rglob('*') if p.is_file()},
        baseline_reuse='Original artifacts and their original recorded timings; no mapper rerun or retiming.',
        batch_workers=args.batch_workers, batch_tasks=args.batch_tasks)
    if args.adaptive_policy == 'shared_policies':
        manifest['adaptive_config'] = dict(previous['original_config'])
        manifest['branch_cap_scope'] = dict(original='Native growth and synchronized live frontier',
            adaptive='Identical per-policy native growth and synchronized live frontier')
    save(args.run/'manifest.json', manifest)
    (args.run/'status').mkdir()
    print(json.dumps(dict(run=str(args.run), candidate_calls=len(tasks), original_artifacts=str(reference))), flush=True)


def spec_and_folder(args):
    spec = read(args.run/'tasks.json')[args.slot]
    return spec, args.run/f"results/{spec['dataset']}/{spec['index']}/{spec['direction']}/{args.method}"


def problem_plan(args, spec):
    from rxn_core import AAMProblem, AAMSearchConfig
    from rxn_core.domain import MolecularEndpoint
    from rxn_core.search_orientation import AAMSearchPlan
    raw = read(args.run/f"inputs/{spec['dataset']}/{spec['index']}/input.json")
    endpoints = [MolecularEndpoint(**{k:v for k,v in raw[n].items()
                    if k in ('elements', 'coordinates', 'wbo', 'label', 'metadata')})
                 for n in ('reactant', 'product')]
    original = AAMProblem(*endpoints, name=raw['name'])
    reverse = spec['direction'] == 'P_to_R'
    problem = AAMProblem(*reversed(endpoints), name=raw['name']) if reverse else original
    config = AAMSearchConfig(**read(args.run/'manifest.json')[args.method+'_config'])
    return raw, AAMSearchPlan(original, problem, config, reverse)


def search(args):
    from rxn_core import search_aam
    from rxn_core.artifacts import write_aam_checkpoint
    from publication_timing import SearchProfiler
    spec, folder = spec_and_folder(args)
    _, plan = problem_plan(args, spec)
    manifest = read(args.run/'manifest.json')
    folder.mkdir(parents=True, exist_ok=False)
    save(folder/'environment.json', environment())
    rows = []
    if args.method == 'original':
        with SearchProfiler(folder/'timing_events') as profiler:
            result = search_aam(plan.problem, plan.config, workers=manifest['original_workers'],
                execution=manifest['original_execution'], intermediate_dir=folder/'cuts', archive_format='checkpoint')
        rows.append(dict(label='full_sweep', archive='cuts/aam.pkl.gz', **profiler.summary(),
                         metrics=asdict(result.metrics), states=len(result.graph.states),
                         terminals=len(result.graph.terminals), capped=result.graph.capped))
    else:
        from rxn_core.adaptive_seed_search import AdaptiveSeedSearch
        from rxn_core.adaptive_cut_search import AdaptiveCutSearch
        cpu, wall = time.process_time(), time.perf_counter()
        policy = manifest.get('adaptive_policy', 'seed_frontier')
        session = (AdaptiveCutSearch(plan.problem, plan.config, policy='shared')
                   if policy == 'shared_policies' else
                   {'seed_frontier':AdaptiveSeedSearch,
                    'cut_interleaved':AdaptiveCutSearch}[policy](plan.problem, plan.config))
        cpu_used, wall_used = time.process_time()-cpu, time.perf_counter()-wall
        for budget in manifest['adaptive_work_budgets']:
            while session.agenda and (budget == 0 or session.work < budget) and cpu_used < manifest['adaptive_cpu_seconds']:
                cpu, wall = time.process_time(), time.perf_counter()
                session.advance()
                cpu_used += time.process_time()-cpu
                wall_used += time.perf_counter()-wall
                if session.work % 100 == 0:
                    save(folder/'progress.json', dict(work=session.work, compute_cpu=cpu_used,
                                                      compute_wall=wall_used, pending=len(session.agenda)))
            cpu, wall = time.process_time(), time.perf_counter()
            result = session.snapshot()
            cpu_used += time.process_time()-cpu
            wall_used += time.perf_counter()-wall
            label = f'work_{session.work:05d}'
            cpu, wall = time.process_time(), time.perf_counter()
            write_aam_checkpoint(result.aam, folder/f'{label}.pkl.gz')
            save(folder/f'{label}_pending.json', result.pending)
            row = dict(label=label, archive=f'{label}.pkl.gz', work=result.work,
                compute_cpu_excluding_persistence_and_loading_seconds=cpu_used,
                compute_wall_excluding_io_seconds=wall_used,
                checkpoint_cpu_seconds=time.process_time()-cpu, checkpoint_wall_seconds=time.perf_counter()-wall,
                growth_calls=result.growth_calls, reused_states=result.reused_states,
                pending=len(result.pending), agenda_exhausted=result.exhausted,
                states=len(result.aam.graph.states), terminals=len(result.aam.graph.terminals),
                capped=result.aam.graph.capped,
                stop='agenda_exhausted' if result.exhausted else
                     'cpu_budget' if cpu_used >= manifest['adaptive_cpu_seconds'] else 'work_budget')
            if policy in ('cut_interleaved', 'shared_policies'):
                row.update(visited_cuts=len(session.sessions), total_cuts=len(session.cuts),
                           native_reuse=session.repair.stats())
            rows.append(row)
            save(folder/'search.json', dict(**spec, method=args.method, rows=rows, complete=False))
            if result.exhausted or cpu_used >= manifest['adaptive_cpu_seconds']:
                break
    save(folder/'search.json', dict(**spec, method=args.method, rows=rows, complete=True,
                                   max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


def holdout_classes(aam, plan, raw):
    """One score per saved terminal, canonicalizing only its heavy relation."""
    from compare_elementary_outputs import features, certificate, event_counts
    feat = features(raw)
    heavy = [i for i,e in enumerate(raw['reactant']['elements']) if e != 'H']
    count = plan.input_problem.source_atom_count
    vectors = sorted({tuple(mapping[i] for i in range(count)) for terminal in aam.graph.terminals
        if len(mapping := plan.to_input_mapping(aam.graph.states[terminal].mapping)) == count})
    groups, cache = {}, {}
    for offset in range(0, len(vectors), 256):
        batch = vectors[offset:offset+256]
        scores = event_counts(plan.input_problem.reactant.wbo, plan.input_problem.product.wbo, batch)
        for vector, counts in zip(batch, scores, strict=True):
            assert len(set(vector)) == count
            assert all(raw['reactant']['elements'][i] == raw['product']['elements'][p] for i,p in enumerate(vector))
            relation = tuple(vector[a] for a in heavy)
            if relation not in cache:
                cache[relation] = hashlib.sha256(certificate(feat, vector, heavy)).hexdigest()
            cid = cache[relation]
            score = int(sum(counts))
            if cid not in groups or score < groups[cid]['events']:
                groups[cid] = dict(id=cid, events=score, counts=list(map(int, counts)), mapping=vector)
    return dict(valid_full_representatives=len(vectors), classes=list(groups.values()),
                best_events=min((r['events'] for r in groups.values()), default=None))


def analyze(args):
    from rxn_core.artifacts import read_aam_checkpoint
    from golden_evaluation import evaluate_planned
    from publication_analysis import rank_archive
    spec, folder = spec_and_folder(args)
    raw, plan = problem_plan(args, spec)
    search_record = read(folder/'search.json')
    # Final output first; smaller saved budgets are separate anytime results.
    for row in reversed(search_record['rows']):
        label = row['label']
        cpu, wall = time.process_time(), time.perf_counter()
        aam = read_aam_checkpoint(folder/row['archive'])
        loading = dict(cpu=time.process_time()-cpu, wall=time.perf_counter()-wall)
        cpu, wall = time.process_time(), time.perf_counter()
        if spec['dataset'] == 'golden':
            classes = rank_archive(aam, plan)
            save(folder/f'{label}_classes.json', classes)
            reference = read(args.run/f"inputs/golden/{spec['index']}/reference.json")
            result = evaluate_planned(aam, plan, reference['features'], reference['mapping'],
                                      seconds=60, query_timeout_ms=1500)
        else:
            result = holdout_classes(aam, plan, raw)
            save(folder/f'{label}_classes.json', result.pop('classes'))
        result.update(label=label, loading=loading, analysis_cpu=time.process_time()-cpu,
                      analysis_wall=time.perf_counter()-wall)
        save(folder/f'{label}_evaluation.json', result)


def worker(args):
    args.slot += args.offset
    spec, folder = spec_and_folder(args)
    status = args.run/f'status/{args.method}_{args.slot}.json'
    record = dict(**spec, method=args.method, environment=environment(), started=time.time())
    save(status, record)
    source = args.run/('original' if args.method == 'original' else 'engine')/'src'
    env = dict(os.environ, PYTHONPATH=f'{source}:{args.run}/engine/bench', RXN_CORE_NATIVE='1',
               PYTHONHASHSEED='0', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    common = ['--run', str(args.run), '--slot', str(args.slot), '--method', args.method]
    for phase, seconds in (('search', 300), ('analyze', 240)):
        if phase == 'analyze' and not (folder/'search.json').exists():
            break
        wall = time.perf_counter()
        with (args.run/f'status/{args.method}_{args.slot}_{phase}.log').open('w') as log:
            completed = subprocess.run(['timeout', '--kill-after=5s', str(seconds), sys.executable,
                str(args.run/'engine/bench/adaptive_full_benchmark.py'), phase, *common],
                env=env, stdout=log, stderr=subprocess.STDOUT)
        record[phase] = dict(exit=completed.returncode, elapsed_including_io=time.perf_counter()-wall)
        save(status, record)
    record.update(complete=True, finished=time.time())
    save(status, record)


def batch_worker(args):
    """Reuse a Slurm allocation, not reaction results, across independent calls.

    Each worker is pinned to one CPU and retains its independent search and
    evaluation watchdogs. Existing attempts are never silently rerun.
    """
    from concurrent.futures import ThreadPoolExecutor
    from queue import SimpleQueue
    manifest = read(args.run/'manifest.json')
    cores = sorted(os.sched_getaffinity(0))
    count = manifest['batch_workers']
    assert len(cores) >= count
    available = SimpleQueue()
    for core in cores[:count]:
        available.put(core)
    start = args.slot*manifest['batch_tasks']
    end = min(start+manifest['batch_tasks'], len(read(args.run/'tasks.json')))
    def run_one(slot):
        status = args.run/f'status/adaptive_{slot}.json'
        if status.exists():
            return dict(slot=slot, existing_attempt=True)
        core = available.get()
        try:
            command = ['taskset', '-c', str(core), 'timeout', '--kill-after=5s', '550',
                sys.executable, str(args.run/'engine/bench/adaptive_full_benchmark.py'),
                'worker', '--run', str(args.run), '--method', 'adaptive', '--slot', str(slot)]
            code = subprocess.run(command).returncode
            row = dict(slot=slot, exit=code, core=core)
            print(json.dumps(row), flush=True)
            return row
        finally:
            available.put(core)
    with ThreadPoolExecutor(max_workers=count) as pool:
        results = list(pool.map(run_one, range(start, end)))
    save(args.run/f'batches/{args.slot}.json', dict(tasks=results, environment=environment()))


def submit_batches(args):
    manifest = read(args.run/'manifest.json')
    tasks = len(read(args.run/'tasks.json'))
    batches = (tasks+manifest['batch_tasks']-1)//manifest['batch_tasks']
    workers = manifest['batch_workers']
    command = [sys.executable, str(args.run/'engine/bench/adaptive_full_benchmark.py'),
               'batch_worker', '--run', str(args.run), '--slot']
    options = ['sbatch', '--parsable', '--partition='+args.partition,
        '--nodes=1', f'--cpus-per-task={workers}', f'--mem={4*workers}G',
        '--time=01:00:00', '--no-requeue',
        f'--array=0-{batches-1}%{max(1,args.cpu_budget//workers)}',
        '--job-name=adaptive_cut_batch', f'--output={args.run}/status/%A_%a.out',
        '--wrap', shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
    if args.exclude:
        options.insert(2, '--exclude='+args.exclude)
    job = subprocess.check_output(options, text=True).strip()
    save(args.run/'batch_submission.json', dict(job=job, command=options,
        allocation_scope='Batch lifetime only; every individual AAM retains a 300-second hard watchdog.'))
    print(job, flush=True)


def submit(args):
    tasks = read(args.run/'tasks.json')
    manifest = read(args.run/'manifest.json')
    jobs = []
    for method in ('original', 'adaptive'):
        cpus = manifest['original_workers'] if method == 'original' else 1
        for start in range(0, len(tasks), 1000):
            stop = min(start+1000, len(tasks))-1
            command = [sys.executable, str(args.run/'engine/bench/adaptive_full_benchmark.py'),
                       'worker', '--run', str(args.run), '--method', method, '--offset', str(start), '--slot']
            options = ['sbatch', '--parsable', '--partition='+args.partition, '--nodes=1',
                f'--cpus-per-task={cpus}', '--mem='+('24G' if method == 'original' else '6G'),
                '--time=00:10:00', '--no-requeue', f'--array=0-{stop-start}%{max(1,args.cpu_budget//8//cpus)}',
                '--job-name=full_'+method, f'--output={args.run}/status/%A_%a.out',
                '--wrap', shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
            if args.exclude:
                options.insert(2, '--exclude='+args.exclude)
            job = subprocess.check_output(options, text=True).strip()
            jobs.append(dict(method=method, job=job, command=options))
            save(args.run/'jobs.json', jobs)
            print(method, start, stop, job, flush=True)


def status(args):
    rows = [read(p) for p in (args.run/'status').glob('*.json')]
    print(json.dumps(dict(started=len(rows), completed=sum(r.get('complete', False) for r in rows),
        searches=dict(Counter((r['method']+':'+str(r['search']['exit'])) for r in rows if 'search' in r)),
        analysis=dict(Counter((r['method']+':'+str(r['analyze']['exit'])) for r in rows if 'analyze' in r))), indent=2))


def compare_case(payload):
    """One saved case; independent of other cases and safe to analyze in parallel."""
    from publication_analysis import union_outcome, merge_classes, certificate_id
    run, specification = payload
    dataset, index = specification['dataset'], specification['index']
    methods = {}
    method_classes = {}
    for method in ('original', 'adaptive'):
        directed = {}
        for direction in DIRECTIONS:
            folder = run/f'results/{dataset}/{index}/{direction}/{method}'
            search_path = folder/'search.json'
            record = read(search_path) if search_path.exists() else {}
            last = record['rows'][-1] if record else {}
            label = last.get('label', 'missing')
            evaluation_path = folder/f'{label}_evaluation.json'
            classes_path = folder/f'{label}_classes.json'
            evaluation = read(evaluation_path) if evaluation_path.exists() else {}
            classes = read(classes_path) if classes_path.exists() else None
            directed[direction] = dict(search=record, evaluation=evaluation, classes=classes)
        method_classes[method] = {direction:None if value['classes'] is None else
            {r['id']:tuple(r['key'][:3]) if dataset=='golden' else r['events'] for r in value['classes']}
            for direction,value in directed.items()}
        modes = {}
        for mode, directions in (('single', [specification['smaller_first']]),
                                 ('bidirectional', list(DIRECTIONS))):
            values = [directed[d] for d in directions]
            available = [v['classes'] for v in values if v['classes'] is not None]
            data = dict(searches_complete=sum(v['search'].get('complete', False) for v in values),
                expected_searches=len(directions), classes_complete=len(available) == len(directions),
                capped_directions=sum(v['search']['rows'][-1]['capped'] for v in values if v['search']),
                compute_cpu=sum(v['search']['rows'][-1]['compute_cpu_excluding_persistence_and_loading_seconds']
                                for v in values if v['search'].get('complete')))
            if dataset == 'golden':
                reference = read(run/f'inputs/golden/{index}/reference.json')
                expected = certificate_id(reference['features'], reference['mapping'])
                classes = merge_classes(available)
                data['reference_recovery'] = union_outcome([
                    v['evaluation'].get('reference_recovery', 'unknown') for v in values])
                best = classes[0]['key'][:3] if classes else None
                hit = next((r for r in classes if r['id'] == expected), None)
                data['representative_event_windows'] = {str(delta):bool(hit and
                    hit['key'][:2] == best[:2] and hit['key'][2] <= best[2]+delta)
                    for delta in range(11)}
                data['class_count'] = len(classes)
            else:
                merged = {}
                for group in available:
                    for row in group:
                        if row['id'] not in merged or row['events'] < merged[row['id']]['events']:
                            merged[row['id']] = row
                data['class_events'] = {k:v['events'] for k,v in merged.items()}
                data['best_events'] = min(data['class_events'].values(), default=None)
                data['full_mapping_found'] = bool(merged)
            modes[mode] = data
        methods[method] = modes
    row = dict(dataset=dataset, index=index, methods=methods)
    row['directional_class_comparison'] = {}
    for direction in DIRECTIONS:
        before, after = (method_classes[method][direction] for method in ('original','adaptive'))
        row['directional_class_comparison'][direction] = (
            dict(status='unavailable') if before is None or after is None else
            dict(status='compared', original=len(before), candidate=len(after),
                missing=sorted(before.keys()-after.keys()),
                worse_scores=sorted(k for k in before.keys() & after.keys() if after[k]>before[k]),
                new=len(after.keys()-before.keys())))
    if dataset == 'holdout':
        row['comparison'] = {}
        for mode in ('single', 'bidirectional'):
            before, after = (methods[m][mode] for m in ('original', 'adaptive'))
            limit = before['best_events']
            windows = []
            if limit is not None:
                for delta in range(3):
                    expected = {k for k,v in before['class_events'].items() if v <= limit+delta}
                    observed = {k for k,v in after['class_events'].items() if v <= limit+delta}
                    windows.append(dict(delta=delta, original=len(expected), shared=len(expected & observed),
                        unresolved=sorted(expected-observed), new=len(observed-expected)))
            row['comparison'][mode] = dict(event_windows=windows,
                scope='Saved representative heavy classes under unchanged original equivalence. '
                      'Unresolved classes require compressed-family queries, not assumed absent.')
    return row


def compare(args):
    """Read saved outputs in parallel; retain the complete fixed denominators."""
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    tasks = read(args.run/'tasks.json')
    specifications = {(t['dataset'],t['index']):t for t in tasks}
    payloads = [(args.run,specifications[(dataset,index)])
                for dataset,count in COUNTS.items() for index in range(count)]
    with ProcessPoolExecutor(max_workers=args.report_workers, mp_context=mp.get_context('spawn')) as pool:
        rows = list(pool.map(compare_case, payloads, chunksize=4))
    totals = {}
    for dataset, count in COUNTS.items():
        selected = [r for r in rows if r['dataset'] == dataset]
        for mode in ('single', 'bidirectional'):
            summary = {}
            for method in ('original', 'adaptive'):
                values = [r['methods'][method][mode] for r in selected]
                summary[method] = dict(cases=count, searches_complete=sum(v['searches_complete'] for v in values),
                    expected_searches=sum(v['expected_searches'] for v in values),
                    compute_cpu_completed=sum(v['compute_cpu'] for v in values),
                    capped_directions=sum(v['capped_directions'] for v in values),
                    classes_complete=sum(v['classes_complete'] for v in values))
                if dataset == 'golden':
                    summary[method]['reference_recovery'] = dict(Counter(v['reference_recovery'] for v in values))
                    summary[method]['event_windows'] = {str(d):sum(v['representative_event_windows'][str(d)]
                                                                 for v in values) for d in range(11)}
                else:
                    summary[method]['full_mapping_found'] = sum(v['full_mapping_found'] for v in values)
            if dataset == 'golden':
                summary['original_recovered_not_yet_recovered_by_adaptive'] = [r['index'] for r in selected
                    if r['methods']['original'][mode]['reference_recovery'] == 'recovered'
                    and r['methods']['adaptive'][mode]['reference_recovery'] != 'recovered']
            else:
                summary['event_comparison'] = dict(Counter('unknown' if None in (
                    a := r['methods']['original'][mode]['best_events'],
                    b := r['methods']['adaptive'][mode]['best_events']) else
                    'adaptive_worse' if b > a else 'adaptive_better' if b < a else 'equal'
                    for r in selected))
            totals[dataset+'_'+mode] = summary
    save(args.run/'comparison/per_case.json', rows)
    save(args.run/'comparison/summary.json', dict(totals=totals,
        scope='Full fixed denominators. Live report; unavailable results remain unknown. '
              'No change to original equivalence, no reference-guided search, no top-k truncation. '
              'CPU is completed measured work only; failed/censored cost is in worker status and Slurm accounting.'))
    print(json.dumps(totals, indent=2))


def elapsed_seconds(value):
    days, clock = value.split('-', 1) if '-' in value else ('0', value)
    seconds = 0
    for component in clock.split(':'):
        seconds = 60*seconds+int(component)
    return 86400*int(days)+seconds


def watch(args):
    """Bounded campaign monitor; only replace allocations that never started.

    Every replacement uses the same frozen worker and limits. A started or
    failed search is never silently retried. Scheduler cleanup state is not a
    substitute for actual saved worker completion records.
    """
    jobs = read(args.run/'jobs.json')
    expected = 2*len(read(args.run/'tasks.json'))
    deadline = time.monotonic()+args.monitor_seconds
    recovered = set()
    journal = []
    last_report = -1
    while time.monotonic() < deadline:
        rows = [read(p) for p in (args.run/'status').glob('*.json')]
        complete = sum(r.get('complete', False) for r in rows)
        save(args.run/'monitor.json', dict(started=len(rows), complete=complete, expected=expected,
             updated=time.time(), replacements=journal, monitor_deadline_reached=False))
        print(json.dumps(dict(started=len(rows), complete=complete, expected=expected)), flush=True)
        if complete == expected or complete//1000 > last_report:
            with (args.run/'monitor_report.log').open('w') as log:
                subprocess.run(['timeout', '--kill-after=5s', '300', sys.executable, __file__,
                    'compare', '--run', str(args.run)], stdout=log, stderr=subprocess.STDOUT)
            last_report = complete//1000
        if complete == expected:
            return
        try:
            listing = subprocess.run(['squeue', '-h', '-j', ','.join(j['job'] for j in jobs),
                '-t', 'CONFIGURING', '-o', '%i|%M|%N'], text=True, capture_output=True, timeout=30)
        except subprocess.TimeoutExpired:
            print('Scheduler query timed out; worker watchdogs remain independent', flush=True)
            time.sleep(30)
            continue
        if listing.returncode == 0:
            for line in listing.stdout.splitlines():
                job_id, elapsed, host = line.split('|')
                if '_' not in job_id or job_id in recovered:
                    continue
                if elapsed_seconds(elapsed) < 600:
                    continue
                parent, local_slot = job_id.split('_')
                job = next(j for j in jobs if j['job'] == parent)
                command = shlex.split(job['command'][-1])
                offset = int(command[command.index('--offset')+1])
                slot = offset+int(local_slot)
                method = job['method']
                worker_status = args.run/f'status/{method}_{slot}.json'
                if worker_status.exists():
                    continue
                subprocess.run(['scancel', job_id], check=True, timeout=30)
                recovered.add(job_id)
                entry = dict(job=job_id, slot=slot, method=method, host=host, elapsed=elapsed,
                             reason='CONFIGURING >= 10 minutes, worker never started', time=time.time())
                # Recheck after cancellation to avoid racing a late worker start.
                if worker_status.exists():
                    entry['replacement'] = None
                    entry['reason'] += '; late-start race, no automatic search retry'
                else:
                    replacement = list(job['command'])
                    replacement = [a for a in replacement if not a.startswith('--array=')]
                    replacement[replacement.index('--partition=cpunodes_nia')] = '--partition=cpunodes'
                    replacement[-1] = shlex.join(command[:-1]+[str(local_slot)])
                    entry['replacement'] = subprocess.check_output(replacement, text=True, timeout=30).strip()
                journal.append(entry)
                save(args.run/'scheduler_recovery.json', journal)
        time.sleep(30)
    save(args.run/'monitor.json', dict(started=len(rows), complete=complete, expected=expected,
         updated=time.time(), replacements=journal, monitor_deadline_reached=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('prepare', 'prepare_revision', 'submit', 'submit_batches',
        'worker', 'batch_worker', 'search', 'analyze', 'status', 'compare', 'watch'))
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--original-commit', default='98b01b1')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--work-budgets', type=int, nargs='+', default=[400, 1600, 6400],
                        help='Saved decision budgets; 0 means exhaust the declared policies, subject to watchdog/CPU budget')
    parser.add_argument('--cpu-seconds', type=float, default=30.)
    parser.add_argument('--cpu-budget', type=int, default=1024)
    parser.add_argument('--partition', default='cpunodes_nia')
    parser.add_argument('--exclude')
    parser.add_argument('--slot', type=int)
    parser.add_argument('--offset', type=int, default=0)
    parser.add_argument('--method', choices=('original', 'adaptive'))
    parser.add_argument('--monitor-seconds', type=int, default=7200)
    parser.add_argument('--reference-run', type=Path)
    parser.add_argument('--adaptive-policy', choices=('seed_frontier','cut_interleaved','shared_policies'), default='cut_interleaved')
    parser.add_argument('--batch-workers', type=int, default=16)
    parser.add_argument('--batch-tasks', type=int, default=64)
    parser.add_argument('--report-workers', type=int, default=16)
    arguments = parser.parse_args()
    globals()[arguments.command](arguments)

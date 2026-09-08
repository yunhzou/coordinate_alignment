"""Reference-blind competitor runs; persisted predictions and separate evaluation.

The mapper process never receives reference atom labels. Each reaction has a
hard process watchdog. Startup, mapping, serialization, and evaluation are
separate measurements. No repairs, balancing, or candidate reranking are used.
"""
import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import random
import socket
import subprocess
import time
import traceback
from functools import lru_cache
import signal
import sys


METHODS = ('rxnmapper', 'localmapper', 'slap_binary', 'slap_weighted', 'indigo', 'chython', 'rdt')


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def mapper_process(connection, method):
    os.setsid()  # The watchdog also owns subprocesses such as RDT's JVM.
    start = time.perf_counter()
    random.seed(20260908)
    import numpy as np
    np.random.seed(20260908)
    try:
        if method in ('rxnmapper', 'localmapper'):
            import torch
            torch.manual_seed(20260908)
            torch.set_num_threads(1)
            if method == 'rxnmapper':
                from rxnmapper import RXNMapper
                mapper = RXNMapper()
                invoke = lambda reaction: mapper.get_attention_guided_atom_maps([reaction])[0]
            else:
                from localmapper import localmapper
                mapper = localmapper(device='cpu')
                invoke = lambda reaction: mapper.get_atom_map(reaction, return_dict=True)
        elif method.startswith('slap_'):
            from slapmapper.aam import SlapAAM
            mapper = SlapAAM(binary=method == 'slap_binary')

            def invoke(reaction):
                mapper.reset()
                mapper.map_smiles(reaction)
                return mapper.results
        elif method == 'indigo':
            from indigo import Indigo
            mapper = Indigo()

            def invoke(reaction):
                value = mapper.loadReaction(reaction)
                value.automap('discard')
                return dict(mapped_rxn=value.smiles())
        elif method == 'chython':
            from chython import smiles

            def invoke(reaction):
                value = smiles(reaction)
                score = value.attention_mapping()
                return dict(mapped_rxn=format(value, 'm'), mapping_return=str(score))
        elif method == 'rdt':
            import base64
            java = subprocess.Popen([os.environ['RDT_JAVA'], '-Xmx4g', '-XX:ActiveProcessorCount=8',
                '-cp', os.environ['RDT_CLASSPATH'], 'RDTGolden'],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
            for line in java.stdout:
                if line.strip() == 'RDT_READY':
                    break
            else:
                raise RuntimeError('RDT JVM exited before readiness')

            def invoke(reaction):
                java.stdin.write(reaction+'\n')
                java.stdin.flush()
                for line in java.stdout:
                    if line.startswith('RDT_RESULT\t'):
                        _, status, wall, cpu, encoded = line.rstrip().split('\t')
                        value = base64.b64decode(encoded).decode()
                        if status != 'OK':
                            raise RuntimeError(value)
                        return dict(mapped_rxn=value, engine_mapping_seconds=float(wall),
                                    engine_cpu_seconds=float(cpu))
                raise RuntimeError('RDT JVM exited during mapping')
        versions = {d.metadata['Name']: d.version for d in importlib.metadata.distributions()}
        connection.send(dict(status='ready', startup_seconds=time.perf_counter()-start,
                             versions=versions, python=platform.python_version()))
    except Exception:
        connection.send(dict(status='setup_error', error=traceback.format_exc()))
        return
    while True:
        reaction = connection.recv()
        if reaction is None:
            return
        wall, cpu = time.perf_counter(), time.process_time()
        try:
            raw = invoke(reaction)
            # Stop mapping clocks before normalizing or serializing output.
            elapsed, cpu_elapsed = time.perf_counter()-wall, time.process_time()-cpu
            if method == 'rdt':
                cpu_elapsed += raw['engine_cpu_seconds']
            if method.startswith('slap_'):
                raw = [dict(mapped_rxn=r['smiles'], cd=float(r['cd']), val=float(r['val']),
                    compressed_graphs=[dict(labels=[int(x) for x in g.labels],
                        edges=[(int(a),int(b),float(w)) for a,neighbors in g.graph.items()
                               for b,w in neighbors.items() if a < b]) for g in r['lgp']])
                        for r in raw]
            candidates = raw if isinstance(raw, list) else [raw]
            connection.send(dict(status='mapped', candidates=candidates,
                                 mapping_seconds=elapsed, mapping_cpu_seconds=cpu_elapsed))
        except Exception:
            connection.send(dict(status='mapping_error', error=traceback.format_exc(),
                mapping_seconds=time.perf_counter()-wall, mapping_cpu_seconds=time.process_time()-cpu))


def initialize(args):
    args.run.mkdir(parents=True, exist_ok=False)
    audit = args.dataset / 'audit.jsonl'
    rows = [json.loads(line) for line in audit.read_text().splitlines()]
    # Deliberately only write label-free inputs into the worker data file.
    (args.run/'inputs.jsonl').write_text(''.join(json.dumps(dict(index=r['index'],
        input_reaction=r['input_reaction']))+'\n' for r in rows))
    save(args.run/'manifest.json', dict(dataset=str(args.dataset.resolve()),
        dataset_manifest=json.loads((args.dataset/'manifest.json').read_text()),
        audit_sha256=hashlib.sha256(audit.read_bytes()).hexdigest(),
        input_sha256=hashlib.sha256((args.run/'inputs.jsonl').read_bytes()).hexdigest(),
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seed=20260908, methods=METHODS,
        settings='upstream defaults; CPU; 1 thread except Chython ONNX default 8 threads',
        slap_variants='binary=True default and binary=False separately; add_Hs=True; break_sym=heavy',
        watchdog_seconds=args.timeout, model_startup_timeout_seconds=300,
        ranking='upstream output order, unmodified; SLAP first is not asserted to be a ranked top-one',
        evaluation='exact heavy-atom relation modulo endpoint chemical symmetry; unmatched atoms retained',
        note='Local reproduction on pinned Golden, not an assertion of identical historical paper protocol.'))


def worker(args):
    output = args.run/args.method
    output.mkdir(exist_ok=True)
    rows = [json.loads(line) for line in (args.run/'inputs.jsonl').read_text().splitlines()]
    selected = set(map(int,args.indices.split(','))) if args.indices else None
    rows = [r for r in rows if (selected is None and r['index'] % args.shards == args.shard)
            or (selected is not None and r['index'] in selected)]
    context = mp.get_context('spawn')
    process = connection = None
    startup_serial = 0
    try:
        for row in rows:
            path = output/f"{row['index']}.json"
            if path.exists():
                continue
            if process is None:
                connection, child = context.Pipe()
                process = context.Process(target=mapper_process, args=(child,args.method))
                process.start()
                child.close()
                if not connection.poll(300):
                    raise TimeoutError('Model startup exceeded 300 seconds')
                startup = connection.recv()
                save(output/f'startup_{args.shard}_{time.time_ns()}_{startup_serial}.json', startup)
                startup_serial += 1
                if startup['status'] != 'ready':
                    raise RuntimeError(startup['error'])
            wall = time.perf_counter()
            connection.send(row['input_reaction'])
            if connection.poll(args.timeout):
                try:
                    result = connection.recv()
                except EOFError:
                    result = dict(status='worker_crash')
                    process.join()
                    connection.close()
                    process = None
            else:
                os.killpg(process.pid, signal.SIGKILL)
                process.join()
                connection.close()
                process = None
                result = dict(status='timeout', watchdog_seconds=args.timeout)
            result.update(index=row['index'], method=args.method,
                request_seconds=time.perf_counter()-wall, hostname=socket.gethostname(),
                slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
                slurm_array_task_id=os.environ.get('SLURM_ARRAY_TASK_ID'))
            save(path,result)
            print(json.dumps({k:v for k,v in result.items() if k not in ('candidates','error')}), flush=True)
            if result['status'] == 'mapping_error' and process is not None:
                # A library/JVM error must not poison later reactions in this shard.
                os.killpg(process.pid, signal.SIGKILL)
                process.join()
                connection.close()
                process = None
    finally:
        if process is not None:
            if process.is_alive():
                os.killpg(process.pid, signal.SIGKILL)
            process.join()


@lru_cache(maxsize=8192)
def signatures(reaction):
    """Use the same endpoint colors and exact mapping relation as our evaluator."""
    import pynauty
    from rdkit import Chem
    from golden_evaluation import prepare, project, colored_graph
    from collections import Counter
    reactants, agents, products = reaction.split('>')
    # Golden supplies all input components on the left. RDT can move spectator
    # components into the SMILES agent field; retain them in the input endpoint.
    reaction = '.'.join(x for x in (reactants,agents) if x)+'>>'+products
    for side in reaction.split('>>'):
        mol = Chem.MolFromSmiles(side)
        if mol is None:
            raise ValueError('Unparseable mapped endpoint')
        labels = [a.GetAtomMapNum() for a in mol.GetAtoms()
                  if a.GetAtomicNum() != 1 and a.GetAtomMapNum()]
        if len(labels) != len(set(labels)):
            raise ValueError('Duplicate heavy-atom map labels on one endpoint')
    _, features, mapping = prepare(reaction)
    # Nauty's certificate encodes the colored partition, not semantic color
    # names. Unlike our fixed-endpoint AAM evaluation, competitors can change
    # endpoints, so the chemical palette must accompany the certificate.
    palette = tuple((tuple(sorted(Counter(map(tuple,f['colors'])).items(),key=repr)),
                     tuple(sorted(Counter(tuple(c) for _,_,c in f['bonds']).items(),key=repr)))
                    for f in features)
    return ((palette,pynauty.certificate(colored_graph(features))),
            (palette,pynauty.certificate(colored_graph(features, project(mapping,features)))))


def report(args):
    # Evaluation runs in our existing environment, not in competitors' envs.
    from collections import Counter
    import statistics
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.warning')
    manifest = json.loads((args.run/'manifest.json').read_text())
    audit = Path(manifest['dataset'])/'audit.jsonl'
    assert hashlib.sha256(audit.read_bytes()).hexdigest() == manifest['audit_sha256']
    rows = [json.loads(line) for line in audit.read_text().splitlines()]
    summaries = {}
    for method in METHODS:
        if args.method and method != args.method:
            continue
        records = []
        wall_times, cpu_times, request_times = [], [], []
        for row in rows:
            path = args.run/method/f"{row['index']}.json"
            if not path.exists():
                continue
            result = json.loads(path.read_text())
            request_times.append(result['request_seconds'])
            if result['status'] == 'mapped':
                wall_times.append(result['mapping_seconds'])
                cpu_times.append(result['mapping_cpu_seconds'])
            evaluation = dict(index=row['index'], status=result['status'])
            if result['status'] == 'mapped':
                endpoints, expected = signatures(row['mapped_reaction'])
                hits, invalid = [], []
                for index, candidate in enumerate(result['candidates']):
                    try:
                        actual_endpoints, actual = signatures(candidate['mapped_rxn'])
                        if actual_endpoints != endpoints:
                            invalid.append(dict(candidate=index, error='Endpoint chemistry changed'))
                            continue
                        if actual == expected:
                            hits.append(index)
                    except Exception as error:
                        invalid.append(dict(candidate=index,error=str(error)))
                evaluation.update(first_correct=0 in hits, any_correct=bool(hits),
                    matching_candidate_indices=hits, candidates=len(result['candidates']), invalid=invalid)
            records.append(evaluation)
        if not records:
            continue
        save(args.run/f'{method}_evaluation.json', records)
        summaries[method] = dict(total=len(rows), finished=len(records),
            statuses=dict(Counter(r['status'] for r in records)),
            first_correct=sum(r.get('first_correct',False) for r in records),
            any_correct=sum(r.get('any_correct',False) for r in records),
            invalid_candidates=sum(len(r.get('invalid',[])) for r in records),
            successful_mapping_wall_sum_seconds=sum(wall_times),
            successful_mapping_cpu_hours=sum(cpu_times)/3600,
            successful_mapping_median_seconds=statistics.median(wall_times) if wall_times else None,
            successful_mapping_max_seconds=max(wall_times,default=None),
            all_request_wall_sum_seconds=sum(request_times),
            timing_note='Per-call compute only; startup, queue, evaluation and artifact saving excluded. Sum is not parallel elapsed.')
    save(args.run/(f'{args.method}_summary.json' if args.method else 'summary.json'), summaries)
    print(json.dumps(summaries, indent=2))


def environment(args):
    """Freeze package versions, VCS pins, model weights and template data."""
    packages, artifacts = {}, {}
    for distribution in importlib.metadata.distributions():
        name = distribution.metadata['Name']
        packages[name] = dict(version=distribution.version,
                             direct_url=distribution.read_text('direct_url.json'))
        if name.lower().replace('_','-') not in ('rxnmapper','localmapper','chython-rxnmap'):
            continue
        for entry in distribution.files or ():
            if Path(entry).suffix not in ('.bin','.pth','.pt','.onnx','.json','.pkl','.safetensors'):
                continue
            path = Path(distribution.locate_file(entry))
            digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda:stream.read(1024*1024),b''):
                    digest.update(block)
            artifacts[str(path)] = dict(sha256=digest.hexdigest(),bytes=path.stat().st_size)
    save(args.run/(Path(sys.prefix).name+'_provenance.json'),
         dict(python=sys.executable,version=sys.version,packages=packages,artifacts=artifacts))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('init','worker','report','environment'))
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--dataset', type=Path, default=Path('data/aam_benchmarks/golden_original_20260906'))
    parser.add_argument('--method', choices=METHODS)
    parser.add_argument('--shard', type=int, default=0)
    parser.add_argument('--shards', type=int, default=1)
    parser.add_argument('--indices')
    parser.add_argument('--timeout', type=int, default=300)
    args = parser.parse_args()
    {'init':initialize, 'worker':worker, 'report':report, 'environment':environment}[args.command](args)


if __name__ == '__main__':
    main()

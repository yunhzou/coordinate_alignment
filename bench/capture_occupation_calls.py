#!/usr/bin/env python3
"""Diagnostic capture at the native occupation boundary; no search modifications."""
import argparse
import gzip
import json
import os
from pathlib import Path
import pickle
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
import search_mcule_retro as scanner
from rxn_core.fragment_matching import symmetry


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--prior-run', type=Path, required=True)
    parser.add_argument('--source-id', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, required=True)
    args = parser.parse_args()
    source = next(row for row in json.loads((args.prior_run / 'combined_progress.json').read_text())
                  ['unfinished_sources'] if row['source_id'] == args.source_id)
    config = json.loads(next((args.prior_run / 'parts').glob('*.summary.json')).read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'manifest.json').write_text(json.dumps(dict(source=source,
        configuration=config, scope='unchanged detector; native call inputs saved before execution'),
        indent=2) + '\n')
    scanner._worker_init(config['target_smiles'], config['config'], None, True, None)
    real_orbit = symmetry.occupation_orbit
    real_materialize = symmetry.materialize_target_coverage_orbit
    context = {}
    counter = 0

    def materialize(candidate, target, **kwargs):
        nonlocal context
        context = dict(target_count=len(scanner._TARGET.graph),
                       target_elements=[target.nodes[i]['element'] for i in range(len(target))],
                       retained_fragments=candidate.retained_fragments,
                       fragment_classes=candidate.fragment_classes)
        return real_materialize(candidate, target, **kwargs)

    def orbit(*parameters):
        nonlocal counter
        index = counter
        counter += 1
        name = f'{os.getpid()}_{index}'
        path = args.output / f'{name}.pkl.gz'
        with gzip.open(path, 'xb', compresslevel=1) as stream:
            pickle.dump(dict(parameters=parameters, context=context), stream, protocol=5)
        started = time.perf_counter()
        print(json.dumps(dict(event='orbit_start', input=str(path), pid=os.getpid(),
                              unix_time=time.time())), flush=True)
        result = real_orbit(*parameters)
        measured = dict(input=str(path), seconds=time.perf_counter()-started,
                        states=len(result))
        (args.output / f'{name}.done.json').write_text(json.dumps(measured) + '\n')
        print(json.dumps(dict(event='orbit_done', **measured)), flush=True)
        return result

    # The normal forked worker pool inherits instrumentation. Neither function
    # changes the arguments, result, seed policy, or caps of the real call.
    symmetry.materialize_target_coverage_orbit = materialize
    symmetry.occupation_orbit = orbit
    result = scanner._detect_one((source['row_index'], source['smiles'], source['source_id']),
                                seed_workers=args.workers - 1)
    with gzip.open(args.output / 'detection.pkl.gz', 'xb', compresslevel=1) as stream:
        pickle.dump(result, stream, protocol=5)


if __name__ == '__main__':
    main()

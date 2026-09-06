#!/usr/bin/env python3
"""Opt-in big-block/gap-first beta; leaves full-bank detection unchanged."""
import argparse
import csv
from dataclasses import asdict
import gzip
import json
from pathlib import Path
import pickle

from rxn_core.fragment_matching import FragmentDetectionConfig
from rxn_core.retrosynthesis.beta import FragmentQueryBank, recommend_big_blocks
from rxn_core.retrosynthesis.config import DEFAULT_ISO_TOLERANCE
from rxn_core.smiles import smiles_to_weighted_graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target-smiles', required=True)
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--id-column', default='Bank ID')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--iso-tolerance', type=float, default=DEFAULT_ISO_TOLERANCE)
    parser.add_argument('--branch-limit', type=int, default=100)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--recommendations', type=int, default=20)
    parser.add_argument('--patterns', type=int, help='Construction patterns (default: min(4, recommendations))')
    args = parser.parse_args()
    config = FragmentDetectionConfig(iso_tolerance=args.iso_tolerance,
                                    branch_limit=args.branch_limit)
    opener = gzip.open if args.catalog.suffix == '.gz' else open
    with opener(args.catalog, 'rt') as stream:
        rows = list(csv.DictReader(stream))
    args.output_dir.mkdir(parents=True, exist_ok=False)
    (args.output_dir / 'checkpoints').mkdir()
    manifest = dict(workflow='big_blocks_beta/v1', target_smiles=args.target_smiles,
        catalog=str(args.catalog.resolve()), sources=[
            {'source_id': r[args.id_column], 'smiles': r['SMILES']} for r in rows],
        explicit_hydrogens=True, config=asdict(config), workers=args.workers,
        recommendations=args.recommendations, exhaustive=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    sources = [(r[args.id_column], smiles_to_weighted_graph(r['SMILES'], expand_hydrogens=True))
               for r in rows]
    target = smiles_to_weighted_graph(args.target_smiles, expand_hydrogens=True)
    with (args.output_dir / 'events.jsonl').open('w') as log:
        def checkpoint(event, evidence):
            path = args.output_dir / 'checkpoints' / f'{checkpoint.count:08d}.pkl.gz'
            with gzip.open(path, 'wb') as stream:
                pickle.dump((event, evidence), stream, protocol=pickle.HIGHEST_PROTOCOL)
            checkpoint.count += 1
            record = dict(event, checkpoint=str(path.resolve()))
            log.write(json.dumps(record) + '\n')
            log.flush()
            print(json.dumps(record), flush=True)
        checkpoint.count = 0
        bank = FragmentQueryBank(sources, target, config=config,
                                 checkpoint=checkpoint, workers=args.workers)
        result = recommend_big_blocks(bank, recommendations=args.recommendations,pattern_limit=args.patterns)
    with gzip.open(args.output_dir / 'result.pkl.gz', 'wb') as stream:
        pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)
    def record(recommendation):
        return dict(uncovered_target_atoms=recommendation.uncovered_target_atoms,
            placements=[dict(source_id=p.candidate.source_id, mapping=p.mapping,
                retained_fragments=p.candidate.retained_fragments,
                target_fragments=[sorted(dict(p.mapping)[a] for a in fragment)
                                  for fragment in p.candidate.retained_fragments],
                target_atoms=p.target_atoms, refined=p.refined)
                for p in recommendation.placements])
    report = dict(workflow='big_blocks_beta/v1', elapsed_seconds=result.elapsed_seconds,
        exhaustive=False, capped_searches=result.capped_searches,
        recommendations=[record(r) for r in result.recommendations],
        best_partial=record(result.best_partial))
    (args.output_dir / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()

"""Export the three audited omissions for the shared reaction viewer.

Reads published witnesses and endpoint geometry only; never runs a mapping search.
"""
import argparse
import gzip
import json
from pathlib import Path

from holdout_minimum_events import AAM, EventPatterns, read, save, sha
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events


ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / 'reports/holdout_forward_cap1000_seed1_20260910'
DIAGNOSTIC = ROOT / 'reports/holdout_missing_pattern_seeds_20260910'


def build(args):
    baseline = {r['index']: r for r in read(args.baseline / 'per_case.json')}
    targets = read(args.diagnostic / 'targets.json')
    diagnostic = read(args.diagnostic / 'per_case.json')
    with gzip.open(args.baseline / 'enumeration.json.gz', 'rt') as stream:
        sweep = json.load(stream)
    sources = [args.baseline / 'per_case.json', args.baseline / 'enumeration.json.gz',
               args.diagnostic / 'targets.json', args.diagnostic / 'per_case.json']
    cases, checks = [], []
    for index in (11, 64, 101):
        source = args.inputs / str(index) / 'input.json'
        sources.append(source)
        raw = read(source)
        endpoints = [raw[k] for k in ('reactant', 'product')]
        canonical = EventPatterns(raw)
        problem = AAMProblem(*(MolecularEndpoint(**e) for e in endpoints))
        base = baseline[index]
        assert base['catalogue_complete']['aam']
        assert base['minima']['aam'] == base['minima']['slap_sweep']
        target = targets[str(index)]['pattern']['id']
        assert target in base['comparisons']['slap_sweep']['slap_only_proven']
        rows = sorted((r for r in diagnostic if r['fresh'] and r['index'] == index),
                      key=lambda r: r['seeds'])
        assert [r['seeds'] for r in rows] == [1, 3, 10, 30]
        assert all(r['cap'] == 1000 and r['direction'] == 'R_to_P' for r in rows)
        records = []

        def record(witness, name, role, provenance, membership):
            actual = canonical.describe(witness['mapping'])
            assert actual['id'] == witness['id']
            assert actual['total'] == base['minima']['aam']
            counts = actual['counts']
            assert bond_events(problem, dict(enumerate(actual['mapping']))) == dict(
                broken=counts['broken'], formed=counts['formed'],
                order_changed=counts['strengthened'] + counts['weakened'], total=actual['total'])
            events = []
            for kind, pairs in actual['events'].items():
                for a, b in pairs:
                    pa, pb = actual['mapping'][a], actual['mapping'][b]
                    events.append(dict(kind=kind, r=[a, b], p=[pa, pb],
                                       wbo=[endpoints[0]['wbo'][a][b], endpoints[1]['wbo'][pa][pb]]))
            checks.append(dict(index=index, name=name, pattern=actual['id'], total=actual['total']))
            records.append(dict(name=name, role=role, pattern=actual['id'],
                                mapping=actual['mapping'], events=events, counts=counts,
                                total=actual['total'], provenance=provenance, membership=membership))

        for ordinal, key in enumerate(base['pattern_ids']['aam'], 1):
            record(base['patterns'][key], f'AAM · 1 seed · alternative {ordinal}', 'aam',
                   dict(source=str(args.baseline / 'per_case.json'), direction='R_to_P',
                        seeds=1, cap=1000), base['membership'][key])
        target_record = sweep[str(index)]['methods']['slap_sweep']['patterns'][target]
        assert all(s['direction'] == 'R_to_P' for s in target_record['sources'])
        record(target_record, 'SLAP + sweep · missing at 1 seed', 'missing',
               dict(source=str(args.baseline / 'enumeration.json.gz'),
                    candidate=target_record['candidate'], sources=target_record['sources']),
               base['membership'][target])
        for row in rows:
            if row['target_result']['status'] == 'represented':
                witness = row['target_result']['witness']
                assert witness['id'] == target
                record(witness, f"AAM · {row['seeds']} seeds · recovered alternative", 'recovered',
                       dict(source=str(args.diagnostic / 'per_case.json'), seeds=row['seeds'],
                            cap=1000, direction='R_to_P', archive=row['archive'],
                            archive_sha256=row['archive_sha256'], terminal=witness.get('terminal')),
                       base['membership'][target])
        cases.append(dict(index=index, name=raw['name'], endpoints=endpoints, records=records,
                          baseline_aam_patterns=len(base['pattern_ids']['aam']),
                          target=target, minimum=base['minima']['aam'],
                          recovery=[dict(seeds=r['seeds'], status=r['target_result']['status'],
                                         search_cpu=r['search_cpu']) for r in rows]))
    args.output.mkdir(parents=True, exist_ok=True)
    save(args.output / 'viewer_data.json', cases)
    print(f'Exported {len(checks)} independently rescored display witnesses.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, default=BASELINE)
    parser.add_argument('--diagnostic', type=Path, default=DIAGNOSTIC)
    parser.add_argument('--inputs', type=Path, default=AAM / 'inputs')
    parser.add_argument('--output', type=Path, default=DIAGNOSTIC)
    build(parser.parse_args())

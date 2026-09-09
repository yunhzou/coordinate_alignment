"""Independently rescore saved comparison witnesses; never rerun either mapper."""
import argparse
from itertools import permutations
import json
from pathlib import Path
import numpy as np
from compare_elementary_outputs import features, certificate
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events


def audit(base, overlap=None, comparison=None):
    overlap = overlap or base / 'elementary_family_overlap_20260909'
    comparison = comparison or base / 'elementary140_output_comparison_20260908'
    checked = 0
    for index in range(140):
        row = json.loads((overlap / f'{index}.json').read_text())
        original = json.loads((comparison / f'{index}.refined.json').read_text())
        raw = json.loads((Path(original['aam_source']) / 'inputs' / str(index) / 'input.json').read_text())
        problem = AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant', 'product')))
        feat = features(raw)
        heavy = [i for i, e in enumerate(problem.reactant.elements) if e != 'H']
        n = len(problem.reactant.elements)
        def cert(pairs):
            m = {int(a): int(b) for a, b in pairs}
            assert set(m) == set(range(n)) == set(m.values())
            assert all(problem.reactant.elements[i] == problem.product.elements[j] for i, j in m.items())
            return certificate(feat, [m[i] for i in range(n)], heavy)
        slap_keys = {c['candidate']: cert(c['mapping']) for c in original['slap']}
        for target in row['targets']:
            if target['status'] == 'represented':
                assert cert(target['mapping_witness']) == slap_keys[target['candidates'][0]]
                checked += 1
        extra = row['extra_aam_pattern_absent_from_saved_slap']
        if extra:
            assert cert(extra['mapping']) not in set(slap_keys.values())
            events = bond_events(problem, dict(extra['mapping']))
            assert all(events[k] == extra['events'][k] for k in ('broken', 'formed', 'order_changed'))
            checked += 1
    oracle = base / 'aam_slap_small_graph_oracle_20260909'
    expanded = json.loads((oracle / 'slap_expanded.json').read_text())
    counts = dict(aam=0, slap_default=0, slap_expanded=0)
    perms = np.array(list(permutations(range(6))))
    a, b = np.triu_indices(6, 1)
    for index in range(100):
        row = json.loads((oracle / f'{index}.json').read_text())
        r, p = np.array(row['R']), np.array(row['P'])
        exact = int(np.sum(r[a,b] != p[perms[:,a], perms[:,b]], axis=1).min())
        assert exact == row['exact_min']
        def cost(vector):
            m = np.array(vector)
            assert sorted(m.tolist()) == list(range(6))
            return int(np.sum(r[a,b] != p[m[a],m[b]]))
        own = json.loads((oracle / f'aam_{index}' / 'summary.json').read_text())
        assert cost(own['mapping']) == own['best_saved_events']
        for c in row['slap']: assert cost(c['mapping']) == c['edits']
        exp = expanded[index]
        assert exp['index'] == index
        values = []
        for attempt in exp['attempts']:
            for c in attempt['results']:
                assert cost(c['mapping']) == c['events']
                values.append(c['events'])
        assert min(values) == exp['best']
        counts['aam'] += own['best_saved_events'] == exact
        counts['slap_default'] += min(c['edits'] for c in row['slap']) == exact
        counts['slap_expanded'] += exp['best'] == exact
    return dict(positive_overlap_and_extra_witnesses_checked=checked, oracle_cases=100,
                exact_minimum_recovery=counts,
                scope='Independent witness/certificate and score audit, not an independent proof of negative family queries.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', type=Path, required=True)
    p.add_argument('--overlap', type=Path)
    p.add_argument('--comparison', type=Path)
    args = p.parse_args()
    print(json.dumps(audit(args.base, args.overlap, args.comparison), indent=2))

"""Offline native-geometry viewer of saved mappings; no search or SMILES conversion."""
import argparse
import json
from pathlib import Path

import numpy as np
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events


def build(base, output):
    cases = []
    for index in (135, 59, 64):
        row = json.loads((base / 'elementary140_output_comparison_20260908' / f'{index}.refined.json').read_text())
        raw = json.loads((Path(row['aam_source']) / 'inputs' / str(index) / 'input.json').read_text())
        endpoints = [raw[k] for k in ('reactant', 'product')]
        problem = AAMProblem(*(MolecularEndpoint(**e) for e in endpoints))
        records = []
        for name, record in [('Our AAM', row['aam_best_example']),
                             ('SLAP', min(row['slap'], key=lambda c: c['events']['total']))]:
            mapping = dict(record['mapping'])
            assert bond_events(problem, mapping) == record['events']
            assert set(mapping) == set(mapping.values()) == set(range(row['atoms']))
            events = []
            r, p = [np.asarray(e['wbo']) for e in endpoints]
            for a in range(len(r)):
                for b in range(a + 1, len(r)):
                    w, v = float(r[a,b]), float(p[mapping[a],mapping[b]])
                    kind = ('broken' if w > .2 and v <= .2 else
                            'formed' if w <= .2 and v > .2 else
                            'order_changed' if w > .2 and v > .2 and abs(w-v) > .5 else None)
                    if kind: events.append(dict(kind=kind, r=[a,b], p=[mapping[a],mapping[b]], wbo=[w,v]))
            assert len(events) == record['events']['total']
            records.append(dict(name=name, mapping=record['mapping'], events=events,
                                counts=record['events'], provenance=record))
        cases.append(dict(index=index, name=row['name'], endpoints=endpoints, records=records,
                          source=row['aam_source']))
    output.parent.mkdir(parents=True, exist_ok=True)
    library = (Path(__file__).resolve().parents[1] / 'src/rxn_core/static/3Dmol-min.js').read_text()
    output.write_text(HTML.replace('__LIBRARY__', library).replace('__DATA__', json.dumps(cases)))
    print(output.resolve())


HTML = (Path(__file__).with_name('elementary_comparison_viewer.html')).read_text()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(args.base, args.output)

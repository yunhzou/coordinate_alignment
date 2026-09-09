"""Read-only numerical audit of real cached TS graphs; no mapping benchmark.

Calls the unmodified pinned SLAP initial-cost builder on continuous WBO input.
Saves both its integer matrix and its own pre-cast neighborhood costs. Neither
mapper is modified. Archived TS assignments illustrate information retention,
not independently verified mapping accuracy. No quantum calculations are run.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from rdkit import Chem

from rxn_core.chemistry_computations.xtb import load_cached_xtb
from rxn_core.matcher.primitives import _growth_edge_supported


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cost_audit(mapper, graphs):
    records = []
    for label in sorted(graphs[0].label2idxs):
        actual, _, infos = mapper._get_cost_matrix(graphs, label)
        left, right = [list(info.values()) for info in infos]
        floating = np.block([
            [np.full((len(a['idxs']), len(b['idxs'])),
                     mapper._diff_nbrs(a['nbrs'], b['nbrs'])) for b in right]
            for a in left
        ])
        assert actual.dtype.kind in 'iu'
        assert np.array_equal(actual, floating.astype(actual.dtype))
        fractional = np.abs(floating - actual) > 1e-10
        lost = (floating > 1e-10) & (actual == 0)
        positions = np.argwhere(lost)
        row_atoms = [i for a in left for i in a['idxs']]
        col_atoms = [i for b in right for i in b['idxs']]
        records.append(dict(
            element=Chem.GetPeriodicTable().GetElementSymbol(label),
            row_atoms=row_atoms, col_atoms=col_atoms,
            returned_dtype=str(actual.dtype), actual_cost=actual.tolist(),
            pre_cast_cost=floating.tolist(), cells=int(actual.size),
            fractional_cells=int(fractional.sum()), positive_to_zero=int(lost.sum()),
            examples=[dict(source_atom=row_atoms[i], target_atom=col_atoms[j],
                           pre_cast=float(floating[i, j]), returned=int(actual[i, j]))
                      for i, j in positions[:3]],
        ))
    return records


def main(args):
    args.output.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('audited_slap', args.slap_core)
    slap = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(slap)
    started, cpu = time.perf_counter(), time.process_time()
    cases = []
    for name in args.cases:
        work = args.source / 'work' / name
        endpoints, graphs, hashes = {}, [], {}
        for side, relative in [('R', 'endpoints/R'), ('TS', 'targets/GT_sp')]:
            folder = work / relative
            elements, coordinates, wbo, xyz = load_cached_xtb(folder)
            endpoints[side] = dict(elements=list(elements), coordinates=coordinates.tolist(),
                                   wbo=wbo.tolist(), xyz=str(xyz))
            hashes[side] = {str(xyz): digest(xyz), str(folder/'wbo'): digest(folder/'wbo')}
            n = len(elements)
            graphs.append(slap.LabeledGraph(
                {i: {j: float(wbo[i, j]) for j in range(n)
                     if i != j and wbo[i, j] >= args.graph_floor} for i in range(n)},
                [Chem.GetPeriodicTable().GetAtomicNumber(e) for e in elements],
            ))
        costs = cost_audit(slap.SlapMapper(binary=False), graphs)
        stage_path = args.source / 'stages' / name / 'ts_stage.json'
        stage = json.loads(stage_path.read_text())
        terms, unscored = [], []
        for mechanism in stage['mechanisms']:
            if mechanism['gt'] is None:
                unscored.append(mechanism['id'])
                continue
            for event in mechanism['gt']['event_terms']:
                wr, wp, wt = [float(event[k]) for k in ('wbo_R', 'wbo_P', 'wbo_T')]
                tr = event['T_pair']
                # The stage's selected mapping supplies these labels. Verify the
                # quoted weights against the original cached arrays, not our prose.
                rp = event['R_pair']
                assert np.isclose(wr, endpoints['R']['wbo'][rp[0]][rp[1]])
                assert np.isclose(wt, endpoints['TS']['wbo'][tr[0]][tr[1]])
                terms.append(dict(
                    mechanism=mechanism['id'], kind=event['kind'],
                    R_pair=rp, P_pair=event['P_pair'], TS_pair=tr,
                    weights_R_TS_P=[wr, wt, wp],
                    binary_R_TS_P=[int(w >= args.graph_floor) for w in (wr, wt, wp)],
                    active_R_edge_supported_tol1=(
                        _growth_edge_supported(wr, wt, 1.0, args.graph_floor)
                        if wr >= args.graph_floor else None),
                    archived_ts_progress=event['ts_progress'],
                ))
        record = dict(name=name, atoms=len(endpoints['R']['elements']),
                      inputs=endpoints, sha256=hashes, costs=costs,
                      archived_event_source=str(stage_path),
                      archived_event_sha256=digest(stage_path), event_examples=terms,
                      archived_mechanisms_without_ts_assignment=unscored)
        (args.output/f'{name}.json').write_text(json.dumps(record, indent=2)+'\n')
        cases.append(dict(name=name, atoms=record['atoms'],
                          **{key: sum(c[key] for c in costs)
                             for key in ('cells', 'fractional_cells', 'positive_to_zero')}))
    result = dict(
        scope='Initial-cost numerical audit of cached reference TS geometries, not mapping accuracy.',
        graph_floor=args.graph_floor, slap_core=str(args.slap_core),
        slap_sha256=digest(args.slap_core), source=str(args.source), cases=cases,
        parent_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        audit_cpu_seconds=time.process_time()-cpu,
        audit_wall_including_case_writes=time.perf_counter()-started,
        notes=[
            'No AAM/SLAP full search, permutation oracle, TS optimization, or chemistry accuracy evaluation.',
            'Initial element-labelled neighborhood matrices only; later refinements can retain other distinctions.',
            'Weighted SMILES uses scaled discrete orders; this audit supplies raw continuous WBO directly.',
            'Integer storage truncates each neighborhood cost, not each individual input edge.',
            'Archived TS assignments are prior algorithm outputs, not independently curated atom identities.',
        ],
    )
    (args.output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--slap-core', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', nargs='+', required=True)
    parser.add_argument('--graph-floor', type=float, default=0.2)
    main(parser.parse_args())

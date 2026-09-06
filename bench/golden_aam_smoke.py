"""Small full-search API smoke test, explicitly not an accuracy benchmark."""
import argparse
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem
from rxn_core import AAMProblem, AAMSearchConfig, search_aam
from rxn_core.artifacts import aam_record
from rxn_core.domain import MolecularEndpoint


def endpoint(smiles, label):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if any(a.GetAtomMapNum() for a in molecule.GetAtoms()):
        raise ValueError('Reference labels must not enter AAM input')
    n = molecule.GetNumAtoms()
    bonds = np.zeros((n, n))
    for bond in molecule.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bonds[a, b] = bonds[b, a] = bond.GetBondTypeAsDouble()
    # search_aam is graph-only. No coordinates or hydrogen identities are
    # inferred from ground truth; stereochemical postprocessing is not run.
    return MolecularEndpoint(tuple(a.GetSymbol() for a in molecule.GetAtoms()),
                             np.zeros((n, 3)), bonds, label=label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--indices', type=int, nargs='+', default=[1001, 979, 267])
    parser.add_argument('--worker', type=int)
    args = parser.parse_args()
    if args.worker is not None:
        row = next(json.loads(s) for s in args.audit.read_text().splitlines()
                   if json.loads(s)['index'] == args.worker)
        left, right = row['input_reaction'].split('>>')
        problem = AAMProblem(endpoint(left, 'R'), endpoint(right, 'P'), f'golden_{args.worker}')
        start = time.perf_counter()
        result = search_aam(problem, AAMSearchConfig(), workers=1,
                            intermediate_dir=args.output/f'{args.worker}_cuts')
        seconds = time.perf_counter()-start
        (args.output/f'{args.worker}.aam.json').write_text(json.dumps(aam_record(result))+'\n')
        report = dict(index=args.worker, explicit_atoms=problem.atom_count,
            wall_search_and_intermediate_seconds=seconds, core_metrics=vars(result.metrics),
            peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            hostname=platform.node(), status='search_completed',
            accuracy_evaluated=False)
        (args.output/f'{args.worker}.json').write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(report))
        return
    args.output.mkdir(parents=True, exist_ok=False)
    records = []
    for index in args.indices:
        command = [sys.executable, __file__, '--audit', str(args.audit),
                   '--output', str(args.output), '--worker', str(index)]
        with (args.output/f'{index}.log').open('w') as stream:
            try:
                completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                    timeout=300, env=dict(os.environ, RXN_CORE_NATIVE='1',
                        OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
                record = (json.loads((args.output/f'{index}.json').read_text())
                          if completed.returncode == 0 else
                          dict(index=index, status='error', exit_code=completed.returncode))
            except subprocess.TimeoutExpired:
                record = dict(index=index, status='timeout', seconds=300)
        records.append(record)
        print(json.dumps(record), flush=True)
    (args.output/'summary.json').write_text(json.dumps(records, indent=2)+'\n')


if __name__ == '__main__':
    main()

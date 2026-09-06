"""Audit saved selected mappings and publish a portable eight-case HTML bundle."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import pickle
import shutil

from rdkit import Chem
from beta_suite_status import update


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert update(args.root), 'Suite is still running'
    progress = json.loads((args.root/'progress.json').read_text())
    summary = []
    bank = None
    for row in progress:
        directory = args.root/row['case']
        case = json.loads((directory/'case.json').read_text())
        run = Path(case['run'])
        manifest = json.loads((run/'query_full/manifest.json').read_text())
        if bank is None:
            with gzip.open(manifest['catalog'], 'rt') as stream:
                bank = {r[manifest['id_column']]: r['SMILES'] for r in csv.DictReader(stream)}
        target = Chem.AddHs(Chem.MolFromSmiles(case['target_smiles']))
        item = dict(case=row['case'], status=row['stage'], elapsed_including_restarts=row['seconds'])
        blind_stem = 'failed_blind' if row['stage'] == 'failed' else 'assemblies_eight'
        for label, stem in [('blind', blind_stem), ('reference', 'reference_eight')]:
            with gzip.open(run/f'{stem}.pkl.gz', 'rb') as stream:
                result = pickle.load(stream)
            for recommendation in result.recommendations:
                covered = set()
                for placement in recommendation.placements:
                    candidate = placement.candidate
                    source = Chem.AddHs(Chem.MolFromSmiles(bank[candidate.source_id]))
                    mapping = dict(candidate.mapping)
                    assert len(mapping) == len(candidate.mapping) == len(set(mapping.values()))
                    assert set(mapping) == set(candidate.retained_atoms)
                    assert set(mapping.values()) == set(candidate.covered_target_atoms)
                    parts = [a for fragment in candidate.retained_fragments for a in fragment]
                    assert len(parts) == len(set(parts)) and set(parts) == set(mapping)
                    for a, b in mapping.items():
                        assert source.GetAtomWithIdx(a).GetAtomicNum() == target.GetAtomWithIdx(b).GetAtomicNum()
                    for a, b in candidate.preserved_source_bonds:
                        sb = source.GetBondBetweenAtoms(a, b)
                        tb = target.GetBondBetweenAtoms(mapping[a], mapping[b])
                        assert sb is not None and tb is not None
                        assert abs(sb.GetBondTypeAsDouble()-tb.GetBondTypeAsDouble()) <= 1.0
                    covered.update(mapping.values())
                assert covered == set(range(target.GetNumAtoms()))
                assert not recommendation.uncovered_target_atoms
            item[label] = dict(complete_assemblies=len(result.recommendations),
                               capped_searches=result.capped_searches)
        item['explicit_target_atoms'] = target.GetNumAtoms()
        summary.append(item)
        destination = args.output/row['case']
        destination.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(directory/'viewer.html', destination/'viewer.html')
    shutil.copyfile(args.root/'index.html', args.output/'index.html')
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

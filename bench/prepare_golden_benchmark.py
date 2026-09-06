"""Pin and audit original Golden data without repairing or excluding records."""
import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import zipfile

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

REPOSITORY = 'Laboratoire-de-Chemoinformatique/Reaction_Data_Cleaning'


def audit_block(block, index):
    row = dict(index=index)
    try:
        _, marker, rxn = block.partition('$RXN')
        if not marker:
            raise ValueError('Missing RDF reaction block')
        reaction = rdChemReactions.ReactionFromRxnBlock(
            '$RXN' + rxn.split('$DTYPE', 1)[0],
            sanitize=False, removeHs=False, strictParsing=True)
        sides = [list(reaction.GetReactants()), list(reaction.GetProducts())]
        for side in sides:
            for mol in side:
                Chem.SanitizeMol(mol)
        row['mapped_reaction'] = '>>'.join('.'.join(
            Chem.MolToSmiles(mol) for mol in side) for side in sides)
        row['has_both_sides'] = all(sides)
        compositions, heavy, labels = [], [], []
        row['reference_mapped_H'] = []
        for side in sides:
            expanded = [Chem.AddHs(mol) for mol in side]
            compositions.append(Counter(a.GetSymbol() for mol in expanded for a in mol.GetAtoms()))
            heavy.append(Counter(a.GetSymbol() for mol in side for a in mol.GetAtoms()
                                 if a.GetAtomicNum() != 1))
            labels.append([a.GetAtomMapNum() for mol in side for a in mol.GetAtoms()
                           if a.GetAtomMapNum()])
            row['reference_mapped_H'].append(sum(a.GetAtomicNum() == 1 and a.GetAtomMapNum() > 0
                                                for mol in side for a in mol.GetAtoms()))
        row.update(explicit_compositions=[dict(x) for x in compositions],
                   balanced_explicit_H=compositions[0] == compositions[1],
                   balanced_heavy=heavy[0] == heavy[1],
                   duplicate_map_labels=any(len(x) != len(set(x)) for x in labels),
                   shared_map_labels=len(set(labels[0]) & set(labels[1])))
        # Erase all reference labels before constructing any algorithm input.
        for side in sides:
            for mol in side:
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)
        row['input_reaction'] = '>>'.join('.'.join(sorted(
            Chem.MolToSmiles(mol) for mol in side)) for side in sides)
        row['status'] = 'parsed'
    except (ValueError, RuntimeError) as error:
        row.update(status='parse_or_sanitize_error', error=str(error))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--commit', help='Use a previously recorded upstream commit')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    commit = args.commit
    if commit is None:
        with urllib.request.urlopen(f'https://api.github.com/repos/{REPOSITORY}/commits/master') as stream:
            commit = json.load(stream)['sha']
    url = f'https://raw.githubusercontent.com/{REPOSITORY}/{commit}/data/golden_dataset.zip'
    with urllib.request.urlopen(url) as stream:
        archive = stream.read()
    (args.output/'golden_dataset.zip').write_bytes(archive)
    with zipfile.ZipFile(io.BytesIO(archive)) as zipped:
        rdf = zipped.read('golden_dataset.rdf')
    (args.output/'golden_dataset.rdf').write_bytes(rdf)
    rows = [audit_block(block, i) for i, block in enumerate(rdf.decode().split('$RFMT')[1:])]
    with (args.output/'audit.jsonl').open('w') as stream:
        for row in rows:
            stream.write(json.dumps(row)+'\n')
    parsed = [r for r in rows if r['status'] == 'parsed']
    summary = dict(repository=REPOSITORY, upstream_commit=commit, download_url=url,
        zip_sha256=hashlib.sha256(archive).hexdigest(), rdf_sha256=hashlib.sha256(rdf).hexdigest(),
        rdkit_version=rdBase.rdkitVersion, total_records=len(rows), parsed_records=len(parsed),
        parse_or_sanitize_errors=len(rows)-len(parsed),
        both_sides=sum(r['has_both_sides'] for r in parsed),
        balanced_heavy=sum(r['has_both_sides'] and r['balanced_heavy'] for r in parsed),
        balanced_explicit_H=sum(r['has_both_sides'] and r['balanced_explicit_H'] for r in parsed),
        duplicate_map_label_records=sum(r['duplicate_map_labels'] for r in parsed),
        duplicate_input_records=len(parsed)-len({r['input_reaction'] for r in parsed}),
        reference_mapped_H_records=sum(any(r['reference_mapped_H']) for r in parsed),
        note='Input audit only, not AAM accuracy. No repairs, filtering, or byproduct completion applied.')
    (args.output/'manifest.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

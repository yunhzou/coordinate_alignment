"""Re-rank saved beta covers without rerunning detection or injecting validation sets."""
import argparse
from dataclasses import replace
import gzip
import json
from pathlib import Path
import pickle
import time

from rxn_core.retrosynthesis.beta_assembly import (
    assembly_metrics, completed_assembly_rank, placement_pattern, rank_complete_assemblies,
)
from rxn_core.smiles import smiles_to_weighted_graph
from rxn_core.subgraph import _coerce_graph


def read(path):
    with gzip.open(path, 'rb') as stream:
        return pickle.load(stream)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--recommendations', type=int, default=20)
    parser.add_argument('--patterns', type=int, default=4)
    parser.add_argument('--validation-stem')
    args = parser.parse_args()
    started = time.perf_counter()
    root = args.run
    manifest = json.loads((root/'query_full/manifest.json').read_text())
    target = _coerce_graph(smiles_to_weighted_graph(manifest['target_smiles'], expand_hydrogens=True), .2)
    previous = read(root/'assemblies.pkl.gz')
    answers = [read(path) for path in sorted((root/'assembly_candidates').glob('*.pkl.gz'))]
    selected = rank_complete_assemblies(answers, target, args.recommendations, args.patterns)
    result = replace(previous, recommendations=selected, best_partial=selected[0],
                     elapsed_seconds=time.perf_counter()-started)
    with gzip.open(root/'assemblies_ranked.pkl.gz', 'wb') as stream:
        pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)
    patterns = {}
    records = []
    for answer in selected:
        pattern = placement_pattern(answer.placements, target)
        patterns.setdefault(pattern, len(patterns)+1)
        metrics = assembly_metrics(answer.placements, target)
        metrics['retention'] = float(metrics['retention'])
        records.append(dict(construction_pattern=patterns[pattern], metrics=metrics,
            sources=[p.source_id for p in answer.placements]))
    report = dict(workflow='big_blocks_beta/ranked', seconds=result.elapsed_seconds,
        pool_size=len(answers), recommendations=records,
        objective=['fragments', '-explicit_atom_retention', 'cuts + connections', 'distinct_species'],
        scope='Saved blind complete assemblies only; validation sets are NOT admitted to this pool')
    if args.validation_stem:
        report['validation_comparison'] = []
        for answer in read(root/f'{args.validation_stem}.pkl.gz').recommendations:
            key = completed_assembly_rank(answer.placements, target)
            metrics = assembly_metrics(answer.placements, target)
            metrics['retention'] = float(metrics['retention'])
            report['validation_comparison'].append(dict(metrics=metrics,
                outranks_saved_blind=sum(key < completed_assembly_rank(a.placements,target) for a in answers)))
    (root/'assemblies_ranked.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()

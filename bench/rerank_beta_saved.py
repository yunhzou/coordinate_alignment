"""Re-rank saved beta covers without rerunning detection or injecting validation sets."""
import argparse
from dataclasses import replace
import gzip
import json
from pathlib import Path
import pickle
import time

from rxn_core.retrosynthesis.beta_assembly import (
    assembly_metrics, assembly_key, dominates, pareto_assembly_ranks, placement_pattern, rank_complete_assemblies,
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
    parser.add_argument('--result-stem', default='assemblies_pareto')
    args = parser.parse_args()
    started = time.perf_counter()
    root = args.run
    manifest = json.loads((root/'query_full/manifest.json').read_text())
    target = _coerce_graph(smiles_to_weighted_graph(manifest['target_smiles'], expand_hydrogens=True), .2)
    previous = read(root/'assemblies.pkl.gz')
    answers = [read(path) for path in sorted((root/'assembly_candidates').glob('*.pkl.gz'))]
    selected = rank_complete_assemblies(answers, target, args.recommendations, args.patterns)
    ranks = pareto_assembly_ranks(answers, target)
    result = replace(previous, recommendations=selected, best_partial=selected[0],
                     elapsed_seconds=time.perf_counter()-started)
    with gzip.open(root/f'{args.result_stem}.pkl.gz', 'wb') as stream:
        pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)
    patterns = {}
    records = []
    for answer in selected:
        pattern = placement_pattern(answer.placements, target)
        patterns.setdefault(pattern, len(patterns)+1)
        metrics = assembly_metrics(answer.placements, target)
        metrics['retention'] = float(metrics['retention'])
        records.append(dict(construction_pattern=patterns[pattern], metrics=metrics,
            pareto_layer=ranks[assembly_key(answer)][0],
            sources=[p.source_id for p in answer.placements]))
    report = dict(workflow='big_blocks_beta/ranked', seconds=result.elapsed_seconds,
        pool_size=len(answers), recommendations=records,
        objective=['maximize explicit-H retention', 'minimize cuts + connections'],
        ordering='Pareto layers; within-layer order is display only; fewer fragments then species break equal-objective ties',
        nondominated_in_pool=sum(rank[0]==1 for rank in ranks.values()),
        nondominated_displayed=sum(ranks[assembly_key(a)][0]==1 for a in selected),
        scope='Saved blind complete assemblies only; validation sets are NOT admitted to this pool')
    if args.validation_stem:
        report['validation_comparison'] = []
        for answer in read(root/f'{args.validation_stem}.pkl.gz').recommendations:
            metrics = assembly_metrics(answer.placements, target)
            blind_metrics = [assembly_metrics(a.placements,target) for a in answers]
            counts = dict(dominates_saved_blind=sum(dominates(metrics,m) for m in blind_metrics),
                dominated_by_saved_blind=sum(dominates(m,metrics) for m in blind_metrics),
                incomparable_saved_blind=sum(not dominates(metrics,m) and not dominates(m,metrics)
                    and (metrics['retention'],metrics['cuts']+metrics['connections']) !=
                        (m['retention'],m['cuts']+m['connections']) for m in blind_metrics))
            metrics['retention'] = float(metrics['retention'])
            report['validation_comparison'].append(dict(metrics=metrics, **counts))
    (root/f'{args.result_stem}.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()

"""Assemble bank-recorded occupations of a specified validation supplier set."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import pickle
import time

from run_beta_distributed import DistributedBank
from rxn_core.retrosynthesis.beta import BetaRecommendation, BetaResult
from rxn_core.retrosynthesis.beta_assembly import assemble_supplier_copies, placement_pattern


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--patterns',type=int,default=4)
    parser.add_argument('--result-stem',default='ground_truth_assemblies')
    args=parser.parse_args()
    started=time.perf_counter()
    case=json.loads(args.case.read_text())
    bank=DistributedBank(args.run,'')
    assert bank.manifest['target_smiles']==case['target_smiles']
    wanted=[r['id'] for r in case['reactants']]
    with gzip.open(bank.manifest['catalog'],'rt') as stream:
        rows={row[bank.manifest['id_column']]:(i,row['SMILES'])
              for i,row in enumerate(csv.DictReader(stream)) if row[bank.manifest['id_column']] in wanted}
    pools=[]
    sources=[]
    for source_id in wanted:
        index,smiles=rows[source_id]
        path=args.run/'query_full/sources'/f'{index}.blocks.pkl.gz'
        with gzip.open(path,'rb') as stream:
            initial=pickle.load(stream)
        options={}
        cached=(args.run/'refinements'/f'{source_id}.pkl.gz').exists()
        for anchor in initial:
            for placement in bank.refine(anchor):
                options.setdefault(placement.key,placement)
        pools.append(tuple(options.values()))
        sources.append(dict(source_id=source_id,smiles=smiles,bank_row=index,
            connected_checkpoint=str(path),connected_occupations=len(initial),
            augmented_occupations=len(options),reused_augmentation=cached))
        print(json.dumps(sources[-1]),flush=True)
    chosen={}
    for placements in assemble_supplier_copies(pools,range(len(bank.target))):
        key=placement_pattern(placements,bank.target)
        if key in chosen:
            continue
        chosen[key]=BetaRecommendation(placements,())
        if len(chosen)>=args.patterns:
            break
    recommendations=tuple(chosen.values())
    result=BetaResult(recommendations,
        recommendations[0] if recommendations else BetaRecommendation((),tuple(range(len(bank.target)))),
        bank.capped_searches,time.perf_counter()-started,tuple(bank.events))
    with gzip.open(args.run/f'{args.result_stem}.pkl.gz','wb',compresslevel=1) as stream:
        pickle.dump(result,stream,protocol=pickle.HIGHEST_PROTOCOL)
    report=dict(scope='Specified ground-truth supplier-set assembly, separate from blind recommendation rank',
        sources=sources,target_atoms=len(bank.target),covered=bool(recommendations),
        capped_searches=bank.capped_searches,seconds=result.elapsed_seconds,
        recommendations=[dict(construction_pattern=i,
            fragments=sum(p.fragment_count for p in r.placements),
            uncovered_target_atoms=r.uncovered_target_atoms,
            placements=[dict(source_id=p.source_id,mapping=p.mapping,
                retained_fragments=p.candidate.retained_fragments) for p in r.placements])
            for i,r in enumerate(recommendations,1)])
    (args.run/f'{args.result_stem}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()

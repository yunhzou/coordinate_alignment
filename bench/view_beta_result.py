#!/usr/bin/env python3
"""Render the beta's saved selected occupations in the existing assembled viewer."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import pickle
import sys

from rdkit import Chem
from rxn_core.retrosynthesis.catalog_index import candidate_entry
from rxn_core.retrosynthesis.compressed_coverage import candidate_target_domains
from rxn_core.retrosynthesis.ranking import build_ranked_assembly, validate_atom_ownership
from rxn_core.smiles import smiles_to_weighted_graph
from rxn_core.subgraph import _coerce_graph


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--result-stem',default='assemblies')
    parser.add_argument('--validation-stem')
    parser.add_argument('--case',type=Path)
    parser.add_argument('--html-output',type=Path)
    args=parser.parse_args()
    root=args.run
    manifest=json.loads((root/'query_full/manifest.json').read_text())
    with gzip.open(manifest['catalog'],'rt') as stream:
        bank={r[manifest['id_column']]:(i,r['SMILES']) for i,r in enumerate(csv.DictReader(stream))}
    with gzip.open(root/f'{args.result_stem}.pkl.gz','rb') as stream:
        result=pickle.load(stream)
    metadata=json.loads((root/f'{args.result_stem}.json').read_text())
    validations=[]
    if args.validation_stem:
        with gzip.open(root/f'{args.validation_stem}.pkl.gz','rb') as stream:
            validations=list(pickle.load(stream).recommendations)
    target=_coerce_graph(smiles_to_weighted_graph(manifest['target_smiles'],expand_hydrogens=True),.2)
    assemblies=[]
    statuses={}
    all_recommendations=list(result.recommendations)+validations
    for recommendation in all_recommendations:
        entries=[]
        for placed in recommendation.placements:
            c=placed.candidate
            row,smiles=bank[c.source_id]
            molecule=Chem.AddHs(Chem.MolFromSmiles(smiles))
            if c.source_id not in statuses:
                with gzip.open(root/'refinements'/f'{c.source_id}.pkl.gz','rb') as stream:
                    _,detection=pickle.load(stream)
                statuses[c.source_id]=(detection.complete,detection.status,detection.best_fragment_size)
            complete,status,best_size=statuses[c.source_id]
            raw={name:getattr(c,name) for name in ('mapping','retained_atoms','covered_target_atoms',
                'leftover_fragments','boundary_bonds','attachment_atoms_source',
                'attachment_atoms_target','preserved_source_bonds','retained_fragments')}
            raw['attachment_trimmed_target_atoms']=c.covered_target_atoms
            record=dict(source_id=c.source_id,representation=smiles,row_index=row,
                complete=complete,status=status,best_fragment_size=best_size)
            item=candidate_entry(record,raw,molecule,Chem.MolToSmiles(Chem.RemoveHs(molecule)),
                candidate_target_domains(c),(),(1,c.retained_atoms))
            mapping=dict(c.mapping)
            item['preserved_target_bonds']=tuple(tuple(sorted((mapping[a],mapping[b])))
                                                 for a,b in c.preserved_source_bonds)
            entries.append(item)
        formed=validate_atom_ownership(entries,target.edges())
        assert formed is not None
        assembly=build_ranked_assembly(entries,formed)
        # The beta does not compute a symmetry-retention rank. Do not display
        # an invented estimate; the viewer shows actual selected occupations.
        assembly['score'].pop('set_symmetry_atom_retention')
        assembly['score'].pop('set_symmetry_heavy_atom_retention')
        position=len(assemblies)
        if position<len(result.recommendations):
            assembly['construction_pattern']=metadata['recommendations'][position]['construction_pattern']
        else:
            assembly['construction_pattern']=f'GT-{position-len(result.recommendations)+1}'
        assemblies.append(assembly)
    parts=[json.loads(p.read_text()) for p in (root/'query_full/parts').glob('*.json')]
    patterns={}
    for rank,assembly in enumerate(assemblies[:len(result.recommendations)],1):
        key=assembly['construction_pattern']
        item=patterns.setdefault(key,dict(pattern=key,recommendation_ranks=[],
            fragment_sizes=[len(r['covered_target_atoms']) for r in assembly['precursors']],
            matched_fragment_pattern=[[sorted(dict(r['mapping'])[a] for a in fragment)
                for fragment in r['retained_fragments']] for r in assembly['precursors']]))
        item['recommendation_ranks'].append(rank)
    report=dict(target_smiles=manifest['target_smiles'],assemblies=assemblies[:len(result.recommendations)],
        validation_assemblies=assemblies[len(result.recommendations):],
        search_scope='Beta final order: fewer fragments, higher explicit-H retention, fewer cuts + connections, fewer distinct species. Ranked among discovered full covers; ground-truth validation is separate.',
        scan_counts=dict(rows=sum(p['rows'] for p in parts),searched=sum(p['rows'] for p in parts),
            matched_precursors=sum(1 for _ in bank),fragment_candidates=sum(p['blocks'] for p in parts),
            capped=result.capped_searches),recommendation_search_truncated=True,
        construction_patterns=list(patterns.values()),uncovered_target_atoms=[])
    # Matched-source count is not recoverable from shard block totals alone.
    report['scan_counts']['matched_precursors']=sum(
        json.loads(p.read_text())['blocks']>0 for p in (root/'query_full/sources').glob('*.json'))
    if args.case:
        expected=[r['id'] for r in json.loads(args.case.read_text())['reactants']]
        report.update(expected_ids=expected,expected_ids_found={r:r in bank for r in expected})
    (root/f'{args.result_stem}_viewer.json').write_text(json.dumps(report,indent=2)+'\n')
    sys.path.insert(0,str(Path(__file__).parents[1]/'tools'))
    from build_retro_db_viewer import _payload,_html
    payload=_payload(report,len(result.recommendations),'Example 5 · assembled construction patterns',
        'covered' if validations else 'not-evaluated',
        'Ground-truth panels validate actual compatible supplier occupations and full target coverage. '
        'They are separate from blind recommendation ranks. Overlapping suppliers are labelled explicitly.')
    output = args.html_output or root/f'{args.result_stem}.html'
    output.write_text(_html(payload))
    print(output,flush=True)


if __name__=='__main__':
    main()

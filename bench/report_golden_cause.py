"""Summarize isolated causal probes without changing reference or search."""
import argparse
import json
from pathlib import Path
import pynauty

from golden_policy_campaign import load_case,save
from golden_evaluation import colored_graph,project
from prepare_golden_benchmark import audit_block
from rxn_core.artifacts import read_graph_checkpoint


def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    original=args.rdf.read_text().split('$RFMT')[1:]
    for index in (986,1285):
        root=args.run/f'golden_cause{index}_audit_20260908'
        details=json.loads((root/'results.json').read_text())
        directory,plan=load_case(args.source,index)
        raw=json.loads((directory/'reference.json').read_text())
        reference=plan.to_search_mapping(dict(raw['mapping']))
        graphs={p:read_graph_checkpoint(root/f'{p}.pkl.gz')
                for p in ('ordinary','ordinary_native','reference_constrained')}
        def certificate(mapping):
            return pynauty.certificate(colored_graph(raw['features'],
                project(plan.to_input_mapping(mapping),raw['features'])))
        def bond_changes(mapping):
            a,b=plan.problem.reactant.wbo,plan.problem.product.wbo
            return int(sum(a[r,s]!=b[mapping[r],mapping[s]]
                           for r in reference for s in reference if r<s))
        witnesses=[dict(graphs['ordinary'].states[t].mapping) for t in graphs['ordinary'].terminals]
        row=dict(details,original_rdf_verified=audit_block(original[index],index)['mapped_reaction']==
            json.loads((root/'input_audit.json').read_text())['mapped_reaction'],
            native_python_graph_equal=graphs['ordinary'].to_record()==graphs['ordinary_native'].to_record(),
            reference_heavy_bond_changes=bond_changes(reference),
            ordinary_heavy_bond_changes=[bond_changes(m) for m in witnesses],
            ordinary_representative_reference_equivalent=any(certificate(m)==certificate(reference) for m in witnesses),
            archive=str(root))
        if index==986:
            swapped=dict(reference);swapped[11],swapped[12]=reference[12],reference[11]
            row['diagnostic_swap_only']=dict(source_atoms=[11,12],
                ordinary_equivalent=any(certificate(m)==certificate(swapped) for m in witnesses),
                heavy_bond_changes=bond_changes(swapped),original_reference_unchanged=True)
        # Necessary reference violations in the observed ordinary trace.
        violations=[];seed=None
        for e in json.loads((root/'ordinary.events.json').read_text()):
            if e['type']=='seed_start':seed=e['seed']
            if e['type']!='commit':continue
            r,s=e['edge']['frag_atom'],e['edge']['ext_atom']
            if r not in reference or s not in reference:continue
            source=float(plan.problem.reactant.wbo[r,s]);target=float(plan.problem.product.wbo[reference[r],reference[s]])
            if abs(source-target)>plan.config.iso_tolerance:
                violations.append(dict(seed=seed,source_edge=[r,s],target_edge=[reference[r],reference[s]],
                    source_wbo=source,target_wbo=target,fragment=e['fragment']))
        row['ordinary_trace_reference_edge_conflicts']=violations
        save(args.output/f'{index}.json',row)
        print(json.dumps(row,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--rdf',type=Path,default=Path('data/aam_benchmarks/golden_original_20260906/golden_dataset.rdf'))
    main(p.parse_args())

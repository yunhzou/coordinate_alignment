"""Same-archive, same-annotations comparison against the previous collector."""
import argparse
import json
from pathlib import Path
import subprocess
import time

from collect_golden_patterns import first_paths
from golden_publication import plans
from golden_policy_campaign import save
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.pattern_collection import PatternEquivalence,extract_path_patterns
from view_golden_mapping import molecules


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    old={'__package__':'rxn_core'}
    exec(subprocess.check_output(['git','show','7eec906:src/rxn_core/pattern_collection.py'],text=True),old)
    audit={r['index']:r for r in map(json.loads,Path('data/aam_benchmarks/golden_original_20260906/audit.jsonl').read_text().splitlines())}
    results=[]
    for index,direction in [(44,'P_to_R'),(865,'P_to_R'),(1033,'P_to_R')]:
        plan=plans(args.source,index)[0][direction];base=args.source/'directions'/str(index)/direction
        aam=read_aam_checkpoint(base/'cuts/aam.pkl.gz');path=first_paths(aam.graph)(json.loads((base/'classes.json').read_text())[0]['terminal'])
        mols=molecules(audit[index]['mapped_reaction'],plan.input_problem)
        atoms=tuple({a.GetIdx():(a.GetFormalCharge(),a.GetIsotope(),a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else '') for a in m.GetAtoms()} for m in mols)
        bonds=tuple({tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))):str(b.GetStereo()) for b in m.GetBonds()} for m in mols)
        if plan.reversed:atoms=atoms[::-1];bonds=bonds[::-1]
        for label,cls in [('previous',old['PatternEquivalence']),('twin_quotient',PatternEquivalence)]:
            eq=cls(plan.problem,atom_tags=atoms,bond_tags=bonds)
            result=extract_path_patterns(path,plan.problem,eq,seconds=3,reverse=plan.reversed)
            record=dict(index=index,direction=direction,method=label,terminal=path.terminal,
                **{k:v for k,v in result.items() if k!='patterns'},pattern_keys=[p['key'] for p in result['patterns']])
            results.append(record);save(args.output,results);print(json.dumps(record),flush=True)

"""Package saved adaptive-seed evidence without rerunning matching."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from adaptive_fragment_pilot import save


RUNS={
    'depth_first':'adaptive_seed_frontier_20260910_WPQlDv',
    'fair_depth':'adaptive_seed_fair_20260910_IF7GXy',
    'novel_regions':'adaptive_seed_novel_20260910_yb0TW1',
    'one_route':'adaptive_seed_route_20260910_4bbcer',
    'balanced_regions':'adaptive_seed_balanced_20260910_Nogs0S',
}


def main(args):
    records=[];provenance={}
    root=Path(__file__).resolve().parents[1]
    for name,directory in RUNS.items():
        run=args.bank/directory
        digest=hashlib.sha256()
        for path in sorted((run/'engine/src').rglob('*.py')):
            digest.update(str(path.relative_to(run/'engine')).encode()+b'\0')
            digest.update(hashlib.sha256(path.read_bytes()).digest())
        provenance[name]=dict(path=str(run),manifest=json.loads((run/'manifest.json').read_text()),
            python_source_sha256=digest.hexdigest())
        for summary in sorted((run/'results').glob('*/summary.json')):
            folder=summary.parent
            row=dict(run=name,slot=int(folder.name),summary=json.loads(summary.read_text()),evidence={})
            for filename in ('validation.json','pattern_comparison.json','pattern_family_queries.json',
                'pair_pattern_comparison.json','pair_pattern_family_queries.json',
                'pair_work_00400_pattern_comparison.json','pair_work_00400_pattern_family_queries.json'):
                path=folder/filename
                if path.exists():row['evidence'][filename]=json.loads(path.read_text())
            records.append(row)
    final=args.bank/RUNS['balanced_regions']
    reference=json.loads((root/'reports/dependency_fragment_repair_20260909/mechanism_check_114.json').read_text())
    raw=json.loads((final/'inputs/114.json').read_text());r,p=[np.array(raw[k]['wbo']) for k in ('reactant','product')]
    witnesses=json.loads((final/'results/3/work_00400_witnesses.json').read_text())
    found={str(ref['id']):None for ref in reference['references']}
    for ordinal,vector in enumerate(witnesses['mappings']):
        aligned=p[np.ix_(vector,vector)];rb,pb=r>.2,aligned>.2
        signature=[np.argwhere(np.triu(mask,1)).tolist() for mask in
            (rb & ~pb,~rb & pb,rb & pb & (np.abs(r-aligned)>.5))]
        for ref in reference['references']:
            if signature==ref['signature'] and found[str(ref['id'])] is None:found[str(ref['id'])]=ordinal
    save(args.output,dict(provenance=provenance,records=records,
        full_tests=[s.attrib for s in ET.parse(final/'tests.xml').getroot().iter('testsuite')],
        vanadium_literal_event_signatures_at_400=found,
        scope='Frozen four-case pilot; original and element-pair score-equivalence results both retained. Not equality of full labelled result sets, all H-transfer patterns, higher-event spectrum or chemical ground-truth coverage.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank',type=Path,default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks'))
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args())

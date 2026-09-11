"""Independent finite checks for the event-family postprocessor on real outputs."""
import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import permutations,product
from pathlib import Path
import time

from enumerate_slap_minimum_events import *
from holdout_minimum_events import KINDS,event_counts,colored_graph,pynauty


def main(args):
    checks = []
    symmetry_checks = 0
    for index,ordinals in [(0,[0]),(25,[0,2,9,10,12]),(101,[0,16]),(132,[0,3])]:
        raw = read(AAM/f'inputs/{index}/input.json')
        canonical = EventPatterns(raw)
        native = read(SWEEP/f'slap_xyz/{index}.json')['candidates']
        scored = read(SWEEP/f'slap_scoring/{index}.refined.json')
        generators = pynauty.autgrp(colored_graph([canonical.features[0]]))[0]
        @lru_cache(None)
        def twin(side,a,b):
            images = list(range(canonical.n))
            images[a],images[b] = b,a
            return exact_action(images,canonical.features[side])
        for ordinal in ordinals:
            candidate = native[ordinal]
            groups = [defaultdict(list),defaultdict(list)]
            for side,graph in enumerate(candidate['graphs']):
                for atom,label in enumerate(graph['labels']):groups[side][label].append(atom)
            labels = sorted(groups[0])
            vectors = []
            for choices in product(*(permutations(groups[1][label]) for label in labels)):
                vector = [None]*canonical.n
                for label,chosen in zip(labels,choices,strict=True):
                    for atom,target in zip(groups[0][label],chosen,strict=True):vector[atom] = target
                vectors.append(vector)
            scores = event_counts(canonical.r,canonical.p,vectors).sum(axis=1)
            minimum = int(min(scores))
            assert minimum==scored['slap'][ordinal]['events']['total']
            patterns = {}
            for vector,score in zip(vectors,scores,strict=True):
                if score==minimum:
                    pattern = canonical.describe(vector)
                    patterns[pattern['id']] = pattern
            symbolic = enumerate_family(raw,canonical,candidate,scored['slap'][ordinal],twin,
                                        time.perf_counter()+60,check_ms=30000)
            assert symbolic['complete'] and set(symbolic['patterns'])==set(patterns)
            for pattern in patterns.values():
                events = tuple(tuple(tuple(pair) for pair in pattern['events'][kind]) for kind in KINDS)
                for g in generators:
                    transported = tuple(tuple(sorted(tuple(sorted((g[a],g[b]))) for a,b in edges)) for edges in events)
                    assert canonical.certificate(transported)==pattern['id']
                    symmetry_checks += 1
                for a,b in [(0,1),(2,3)]:
                    if len(events[a])==len(events[b]):continue
                    swapped = list(events)
                    swapped[a],swapped[b] = swapped[b],swapped[a]
                    assert canonical.certificate(tuple(swapped))!=pattern['id']
                    symmetry_checks += 1
            checks.append(dict(index=index,candidate=ordinal,all_assignments=len(vectors),minimum=minimum,
                minimum_assignments=int(sum(scores==minimum)),event_classes=len(patterns),
                symbolic_complete=True,exact_pattern_set_agreement=True))
    result = dict(families=checks,explicit_assignments_checked=sum(r['all_assignments'] for r in checks),
        symmetry_and_event_type_checks=symmetry_checks,failures=0,
        scope='Exhaustive within-label assignments for ten small real SLAP families validate event blocking and twin symmetry reduction. '
              'This is a postprocessing verification, not a mapper accuracy benchmark.')
    save(args.run/'validation/finite_checks.json',result)
    print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    main(parser.parse_args())

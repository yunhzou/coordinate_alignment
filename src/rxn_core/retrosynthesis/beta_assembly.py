"""Assembly of whole beta occupations, independent of fragment discovery."""
from heapq import heappop, heappush
from itertools import count

from .assembly import construction_pattern
from ..search_graph import frozen_value


def placement_pattern(placements, target):
    items=[]
    for placement in placements:
        mapping=dict(placement.mapping)
        items.append(dict(target_fragment_atoms=tuple(
            tuple(sorted(mapping[a] for a in part))
            for part in placement.candidate.retained_fragments),
            preserved_target_bonds=tuple(tuple(sorted((mapping[a],mapping[b])))
                for a,b in placement.candidate.preserved_source_bonds)))
    labels=tuple((data['element'],frozen_value(data.get('features',{})))
                 for _,data in target.nodes(data=True))
    bonds=tuple((a,b,float(target.graph['wbo_matrix'][a,b])) for a,b in target.edges())
    return construction_pattern(items,labels,bonds)


def assemble_supplier_copies(option_pools, target_atoms):
    """Lazily rank full covers for specified source copies by fragment count.

    Each pool is the allowed *whole* occupations of one copy. Repeating a pool
    allows another copy. Overlap is allowed. No independent atom choices,
    sampled witnesses, or coverage-only deduplication replace the saved tuples.
    All copies have fixed input atom cost; the fragment-count suffix bound is
    therefore sufficient for ranking these assemblies. This is an explicit
    supplier-set query, not a claim about its blind recommendation rank.
    """
    pools=tuple(tuple(sorted({p.key:p for p in options}.values(),
                            key=lambda p:(p.fragment_count,p.key))) for options in option_pools)
    if not pools or any(not pool for pool in pools):
        return
    full=frozenset(target_atoms)
    unions=[frozenset()]*(len(pools)+1)
    lower=[0]*(len(pools)+1)
    for i in range(len(pools)-1,-1,-1):
        unions[i]=unions[i+1] | frozenset(a for p in pools[i] for a in p.covered_atoms)
        lower[i]=lower[i+1]+pools[i][0].fragment_count
    queue=[]
    serial=count()
    def push(prefix,covered,cost,index):
        level=len(prefix)
        if index == len(pools[level]):
            return
        block=pools[level][index]
        bound=cost+block.fragment_count+lower[level+1]
        heappush(queue,(bound,next(serial),prefix,covered,cost,index))
    if full <= unions[0]:
        push((),frozenset(),0,0)
    while queue:
        _,_,prefix,covered,cost,index=heappop(queue)
        level=len(prefix)
        push(prefix,covered,cost,index+1)
        block=pools[level][index]
        selection=prefix+(block,)
        covered=covered | block.covered_atoms
        if not full <= covered | unions[level+1]:
            continue
        if len(selection)==len(pools):
            if full <= covered:
                yield selection
        else:
            push(selection,covered,cost+block.fragment_count,0)

"""Assembly of whole beta occupations, independent of fragment discovery."""
from heapq import heappop, heappush
from itertools import count
from fractions import Fraction

from .assembly import construction_pattern
from ..search_graph import frozen_value


def assembly_metrics(placements, target):
    """Whole-set metrics; overlapping claims never multiply product yield.

    Cuts and connections describe geometric support, not balanced reaction edits.
    Every physical source copy contributes its explicit-H input atom count.
    """
    covered = frozenset(a for p in placements for a in p.covered_atoms)
    carried = set()
    for p in placements:
        mapping = dict(p.mapping)
        carried.update(tuple(sorted((mapping[a], mapping[b])))
                       for a, b in p.candidate.preserved_source_bonds)
    edges = {tuple(sorted(edge)) for edge in target.edges()}
    return dict(
        uncovered=len(target) - len(covered),
        fragments=sum(p.fragment_count for p in placements),
        retention=Fraction(len(covered), sum(p.input_atom_count for p in placements)),
        cuts=sum(len(p.candidate.boundary_bonds) for p in placements),
        connections=len(edges - carried),
        species=len({p.source_id for p in placements}))


def completed_assembly_rank(placements, target):
    """Final beta objective, separate from the discovery queue heuristic.

    No fitted weights or source identity bonuses. Exact fractions avoid float
    ties. Source IDs only provide a reproducible tie-break after all objectives.
    """
    m = assembly_metrics(placements, target)
    if m['uncovered']:
        raise ValueError('Final ranking requires complete target coverage')
    return (m['fragments'], -m['retention'], m['cuts'] + m['connections'],
            m['species'], tuple(sorted(p.key for p in placements)))


def rank_complete_assemblies(answers, target, recommendations, pattern_limit):
    """Rank saved complete sets, reserving a representative per chosen pattern."""
    rank = lambda answer: completed_assembly_rank(answer.placements, target)
    unique = {tuple(sorted(p.key for p in a.placements)): a for a in answers}
    buckets = {}
    for answer in sorted(unique.values(), key=rank):
        pattern = placement_pattern(answer.placements, target)
        buckets.setdefault(pattern, []).append(answer)
    chosen = list(buckets.values())[:pattern_limit]
    selected = [bucket[0] for bucket in chosen]
    alternatives = sorted((a for bucket in chosen for a in bucket[1:]), key=rank)
    selected.extend(alternatives[:max(0, recommendations - len(selected))])
    return tuple(sorted(selected, key=rank))


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

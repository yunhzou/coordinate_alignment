"""Symbolic realization of compressed AAM paths: domains plus correlated groups.

Automorph-domain summaries are not independent assignment pools. Only native
symmetry_domains provide local pool permutations; exact group factors retain
their correlations. No group elements or complete bijections are enumerated.
"""
from functools import lru_cache
from dataclasses import dataclass
from itertools import combinations
import time


@lru_cache(maxsize=512)
def group_factors(generators):
    from sympy.combinatorics import Permutation,PermutationGroup
    if not generators:return ()
    group=PermutationGroup([Permutation(list(g)) for g in generators])
    return tuple(tuple(tuple(p(i) for i in range(group.degree)) for p in level.values())
                 for level in reversed(group.basic_transversals) if len(level)>1)


def compile_path(path, problem, reference, *, source_atoms, source_generators=(),
               target_generators=(), complete_reference=True,
               projected_atoms=None):
    """Compile the shared correlated path constraints, without solving.

    Endpoint generators normalize reference identity only. They do not change
    the returned physical mapping. ``source_atoms`` is the reference-scored
    atom set (e.g. heavy atoms); all mapped hydrogens still enter the path model.
    ``projected_atoms`` requests a rejection-only relaxation: its UNSAT result
    is conclusive for the full model, but its SAT result requires a full query.
    """
    import z3
    start=time.perf_counter();solver=z3.Solver();serial=0;program=[]
    nr,np_=problem.source_atom_count,problem.target_atom_count
    representative=path.mapping
    if projected_atoms is not None:
        selected=set(projected_atoms)
        representative={r:p for r,p in representative.items() if r in selected}
    def constant(v):return (int(v),frozenset((int(v),)))
    def symbol(domain):
        nonlocal serial
        domain=tuple(sorted(set(domain)))
        if len(domain)==1:return constant(domain[0])
        value=z3.Int(f'q{serial}');serial+=1
        solver.add(z3.Or(*(value==p for p in domain)))
        return value,frozenset(domain)
    def lookup(value,table):
        expression,support=value
        if isinstance(expression,int):return table.get(expression,value)
        images=[(p,table.get(p,constant(p))) for p in sorted(support)]
        if all(isinstance(v[0],int) and v[0]==p for p,v in images):return value
        out=images[-1][1][0]
        for p,v in reversed(images[:-1]):out=z3.If(expression==p,v[0],out)
        return out,frozenset().union(*(v[1] for _,v in images))
    def act(values,generators,label):
        for factor in group_factors(tuple(tuple(g) for g in generators)):
            active=set().union(*(v[1] for v in values))
            moved={p for p in active if p<len(factor[0]) and any(g[p]!=p for g in factor)}
            if not moved:continue
            choice=symbol(range(len(factor)))
            table={p:lookup(choice,{i:constant(g[p]) for i,g in enumerate(factor)}) for p in moved}
            values=[lookup(v,table) for v in values]
            program.append(('group',label,choice[0],factor))
        return values
    values=[constant(representative.get(r,np_)) for r in range(nr)]
    placements=[]
    for edge in path.transitions:
        placement=path.graph.fragment_placement(edge)
        if placement is not None:
            if placement.target_generators is None:
                raise ValueError('Finalize path symmetry before querying its family')
            placements.append((edge,placement))
    # Later decisions are in the coordinate frame fixed by earlier witnesses.
    for edge,placement in reversed(placements):
        for domain in placement.symmetry_domains:
            if domain.source=='exact_automorph_group':continue
            pool=domain.p_atoms
            if len(pool)<2:continue
            if not set(pool).intersection(set().union(*(v[1] for v in values))):continue
            table={p:symbol(pool) for p in pool}
            solver.add(z3.Distinct(*(v[0] for v in table.values())))
            locked=dict(path.graph.states[path.graph.transitions[edge].source].mapping)
            block_atoms={r for block in placement.symmetry_domains for r in block.r_atoms
                         if block.source!='exact_automorph_group'}
            fixed=set(locked.values())|{dict(placement.representative_assignments)[r]
                for r in placement.exact_fixed}|{p for r,p in placement.representative_assignments
                                                if r not in block_atoms}
            for p in fixed.intersection(pool):solver.add(table[p][0]==p)
            values=[lookup(v,table) for v in values]
            program.append(('domain',edge,{p:v[0] for p,v in table.items()},pool))
        values=act(values,[g.images for g in placement.target_generators],edge)
    mapped=tuple(representative)
    if len(mapped)>1:
        solver.add(z3.Distinct(*(z3.IntVal(values[r][0]) if isinstance(values[r][0],int)
                                else values[r][0] for r in mapped)))
    for r in mapped:
        solver.add(z3.Or(*(values[r][0]==p for p in values[r][1]
                          if p<np_ and problem.reactant.elements[r]==problem.product.elements[p])))
    cuts={tuple(sorted(e)) for e in path.context.cuts}
    required={}
    for _,fragment in placements:
        deferred={tuple(sorted(e)) for e in fragment.deferred_edges}
        for a,b in combinations(fragment.r_atoms,2):
            if a not in representative or b not in representative:continue
            pair=tuple(sorted((a,b)));w=problem.reactant.wbo[a,b]
            if pair in cuts or pair in deferred or w<path.context.graph_floor:continue
            required[pair]=float(w)
    for (a,b),w in required.items():
        va,vb=values[a],values[b]
        allowed=[(x,y) for x in va[1] for y in vb[1]
            if x!=y and x<np_ and y<np_ and problem.product.wbo[x,y]>=path.context.graph_floor
            and abs(w-problem.product.wbo[x,y])<=path.context.iso_tolerance+1e-9]
        solver.add(z3.Or(*(z3.And(va[0]==x,vb[0]==y) for x,y in allowed)))
    raw_values=values
    # Normalize endpoint equivalence without touching the physical witness.
    source_choices=act([constant(r) for r in source_atoms],source_generators,'source_equivalence')
    normalized=[lookup(v,dict(enumerate(raw_values))) for v in source_choices]
    normalized=act(normalized,target_generators,'target_equivalence')
    for atom,value in zip(source_atoms,normalized):
        if complete_reference or atom in reference:solver.add(value[0]==reference.get(atom,np_))
    return CompiledFamily(solver,raw_values,representative,program,required,path,problem,
                          projected_atoms,time.perf_counter()-start)


@dataclass
class CompiledFamily:
    solver: object
    values: list
    representative: dict
    program: list
    required: dict
    path: object
    problem: object
    projected_atoms: object
    encoding_seconds: float

    def realize(self, model):
        return _realize(self,model)


def query_path(path, problem, reference, *, source_atoms, source_generators=(),
               target_generators=(), complete_reference=True, timeout_ms=10000,
               projected_atoms=None):
    """Return a certified full explicit-atom realization, UNSAT, or unknown."""
    import z3
    start=time.perf_counter()
    compiled=compile_path(path,problem,reference,source_atoms=source_atoms,
        source_generators=source_generators,target_generators=target_generators,
        complete_reference=complete_reference,projected_atoms=projected_atoms)
    solver=compiled.solver
    left=timeout_ms-1000*(time.perf_counter()-start)
    if left<=0:return 'unknown',dict(reason='encoding budget',seconds=time.perf_counter()-start)
    solver.set(timeout=max(1,int(left)))
    status=solver.check()
    if status!=z3.sat:return ('not_recovered' if status==z3.unsat else 'unknown'),dict(
        reason='unsat' if status==z3.unsat else solver.reason_unknown(),seconds=time.perf_counter()-start)
    result=compiled.realize(solver.model())
    result['seconds']=time.perf_counter()-start
    return 'recovered',result


def _realize(compiled,model):
    raw_values=compiled.values;representative=compiled.representative
    program=compiled.program;required=compiled.required
    path=compiled.path;problem=compiled.problem;mapped=tuple(representative)
    def number(v):return v if isinstance(v,int) else model.eval(v,model_completion=True).as_long()
    actual=dict(representative);actions=[]
    for kind,edge,choice,data in program:
        if edge in ('source_equivalence','target_equivalence'):continue
        permutation={p:number(v) for p,v in choice.items()} if kind=='domain' else dict(enumerate(data[number(choice)]))
        actual={r:permutation.get(p,p) for r,p in actual.items()}
        actions.append(dict(kind=kind,transition=edge,permutation=sorted(permutation.items())))
    assert actual=={r:number(raw_values[r][0]) for r in mapped}
    assert len(set(actual.values()))==len(actual)
    assert all(problem.reactant.elements[r]==problem.product.elements[p] for r,p in actual.items())
    assert all(problem.product.wbo[actual[a],actual[b]]>=path.context.graph_floor and
        abs(w-problem.product.wbo[actual[a],actual[b]])<=path.context.iso_tolerance+1e-9
        for (a,b),w in required.items())
    return dict(mapping=sorted(actual.items()),actions=actions,
        checked_fragment_edges=len(required),
        scope='full_explicit' if compiled.projected_atoms is None else 'projection_only')

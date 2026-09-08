"""Exact bounded event scoring of one compressed path, without orbit expansion.

This is an optional postprocessor. Matching and the stored search graph are
unchanged. Score invariance is certified before constructing a solver model.
"""
from collections import defaultdict
from itertools import combinations
import time

from .family_query import compile_path


def bond_events(problem,mapping,*,reverse=False,bond_floor=.2,event_tolerance=.5):
    r,p=problem.reactant,problem.product
    mapping=dict(mapping)
    if reverse:r,p=p,r;mapping={v:k for k,v in mapping.items()}
    inverse={v:k for k,v in mapping.items()};broken=formed=changed=0
    for a,b in combinations(range(r.atom_count),2):
        w=r.wbo[a,b]
        if w<=bond_floor:continue
        if a not in mapping or b not in mapping:broken+=int(a in mapping or b in mapping)
        elif p.wbo[mapping[a],mapping[b]]<=bond_floor:broken+=1
        elif abs(w-p.wbo[mapping[a],mapping[b]])>event_tolerance:changed+=1
    for a,b in combinations(range(p.atom_count),2):
        if p.wbo[a,b]<=bond_floor:continue
        if a not in inverse or b not in inverse or r.wbo[inverse[a],inverse[b]]<=bond_floor:formed+=1
    return dict(broken=broken,formed=formed,order_changed=changed,total=broken+formed+changed)


def invariant_score(path,problem,*,bond_floor=.2,event_tolerance=.5):
    """Certify each generating action preserves all relevant bond-score labels."""
    r,p=problem.reactant,problem.product;n=p.atom_count
    weights=sorted({float(r.wbo[a,b]) for a,b in combinations(range(r.atom_count),2)
                    if r.wbo[a,b]>bond_floor})
    def label(w):
        return (True,*(abs(w-v)>event_tolerance for v in weights)) if w>bond_floor else (False,)
    labels=[[label(float(w)) for w in row] for row in p.wbo]
    tested=0
    def preserves(g):
        nonlocal tested
        tested+=1
        return all(labels[a][b]==labels[g[a]][g[b]]
                   for a in range(n) if g[a]!=a for b in range(n))
    for edge in path.transitions:
        placement=path.graph.fragment_placement(edge)
        if placement is None:continue
        if placement.target_generators is None:raise ValueError('Path symmetry must be finalized')
        for generator in placement.target_generators:
            if not preserves(generator.images):return False,tested
        for domain in placement.symmetry_domains:
            if domain.source=='exact_automorph_group':continue
            for a,b in zip(domain.p_atoms,domain.p_atoms[1:]):
                g=list(range(n));g[a],g[b]=g[b],g[a]
                if not preserves(g):return False,tested
    return True,tested


def validate_representative(path,problem):
    mapping=path.mapping
    assert len(set(mapping.values()))==len(mapping)
    assert all(problem.reactant.elements[r]==problem.product.elements[p] for r,p in mapping.items())
    cuts={tuple(sorted(e)) for e in path.context.cuts}
    for edge in path.transitions:
        fragment=path.graph.fragment_placement(edge)
        if fragment is None:continue
        deferred={tuple(sorted(e)) for e in fragment.deferred_edges}
        for a,b in combinations(fragment.r_atoms,2):
            if a not in mapping or b not in mapping:continue
            w=problem.reactant.wbo[a,b]
            if tuple(sorted((a,b))) in cuts|deferred or w<path.context.graph_floor:continue
            assert problem.product.wbo[mapping[a],mapping[b]]>=path.context.graph_floor
            assert abs(w-problem.product.wbo[mapping[a],mapping[b]])<=path.context.iso_tolerance+1e-9


def event_objective(compiled,*,reverse=False,bond_floor=.2,event_tolerance=.5):
    """Local bond tables only; constant bond terms are eliminated exactly."""
    import z3
    r,p=compiled.problem.reactant,compiled.problem.product
    values=compiled.values
    if reverse:
        inverse=[];sentinel=r.atom_count
        full=len(compiled.representative)==p.atom_count
        for target in range(p.atom_count):
            candidates=[(a,v[0]) for a,v in enumerate(values) if target in v[1]]
            fixed=next((a for a,v in candidates if isinstance(v,int) and v==target),None)
            if fixed is not None:inverse.append((fixed,frozenset((fixed,))));continue
            expr=sentinel
            for a,v in reversed(candidates):expr=z3.If(v==target,a,expr)
            support={a for a,v in candidates}
            if not full:support.add(sentinel)
            inverse.append((expr,frozenset(support)))
        values=inverse;r,p=p,r
    sentinel=p.atom_count
    constant=sum(p.wbo[a,b]>bond_floor for a,b in combinations(range(p.atom_count),2))
    terms=[];lower=int(constant);removed=0
    for a,b in combinations(range(r.atom_count),2):
        w=float(r.wbo[a,b])
        if w<=bond_floor:continue
        va,vb=values[a],values[b];groups=defaultdict(list)
        for x in va[1]:
            for y in vb[1]:
                if x==y and x!=sentinel:continue # injectivity is already constrained
                if x==sentinel and y==sentinel:cost=0
                elif x==sentinel or y==sentinel:cost=1
                elif p.wbo[x,y]<=bond_floor:cost=1
                else:cost=0 if abs(w-p.wbo[x,y])>event_tolerance else -1
                groups[cost].append((x,y))
        assert groups
        lower+=min(groups)
        if len(groups)==1:
            constant+=next(iter(groups));removed+=1;continue
        default=max(groups,key=lambda k:len(groups[k]));term=default
        for cost,pairs in groups.items():
            if cost==default:continue
            rows=defaultdict(list)
            for x,y in pairs:rows[x].append(y)
            condition=z3.Or(*(z3.And(va[0]==x,z3.Or(*(vb[0]==y for y in ys)))
                              for x,ys in rows.items()))
            term=z3.If(condition,cost,term)
        terms.append(term)
    objective=z3.IntVal(int(constant))+z3.Sum(terms)
    return objective,max(0,lower),dict(variable_bond_terms=len(terms),constant_bond_terms=removed)


def tighten_supports(compiled):
    """Arc consistency and forced-image propagation; never split group factors.

    Every deletion is entailed by the existing fragment-edge, element, or
    injectivity constraints. Coupled group choices remain in the solver.
    """
    import z3
    problem=compiled.problem;sentinel=problem.target_atom_count
    domains=[{p for p in support if p==sentinel or
              problem.reactant.elements[r]==problem.product.elements[p]}
             for r,(_,support) in enumerate(compiled.values)]
    before=sum(map(len,domains));changed=True
    while changed:
        changed=False
        fixed={next(iter(d)) for d in domains if len(d)==1 and sentinel not in d}
        for r,domain in enumerate(domains):
            if len(domain)>1:
                new=domain-fixed
                if new!=domain:domains[r]=new;changed=True
        for (a,b),w in compiled.required.items():
            for source,target in ((a,b),(b,a)):
                new={p for p in domains[source] if p<sentinel and any(
                    q<sentinel and p!=q and
                    problem.product.wbo[p,q]>=compiled.path.context.graph_floor and
                    abs(w-problem.product.wbo[p,q])<=compiled.path.context.iso_tolerance+1e-9
                    for q in domains[target])}
                if new!=domains[source]:domains[source]=new;changed=True
        assert all(domains), 'Stored representative must satisfy the compiled path'
    narrowed=[]
    for (expression,support),domain in zip(compiled.values,domains):
        if domain!=set(support):compiled.solver.add(z3.Or(*(expression==p for p in sorted(domain))))
        narrowed.append((next(iter(domain)) if len(domain)==1 else expression,frozenset(domain)))
    compiled.values=narrowed
    return dict(atom_support_before=before,atom_support_after=sum(map(len,domains)))


def minimize_events(path,problem,*,reverse=False,seconds=10.,bond_floor=.2,event_tolerance=.5):
    """Return a proven minimum or explicit bounds and an achievable witness."""
    import z3
    start=time.perf_counter();deadline=start+seconds
    validate_representative(path,problem)
    witness=path.mapping;events=bond_events(problem,witness,reverse=reverse,
        bond_floor=bond_floor,event_tolerance=event_tolerance);upper=events['total']
    invariant,tested=invariant_score(path,problem,bond_floor=bond_floor,event_tolerance=event_tolerance)
    invariant_seconds=time.perf_counter()-start
    result=dict(representative_score=upper,generators_tested=tested,invariant_seconds=invariant_seconds,
                method='invariance_certificate' if invariant else 'bounded_symbolic',solver_queries=0)
    if invariant:
        result.update(lower_bound=upper,upper_bound=upper,optimal=True,events=events,
                      mapping=sorted(witness.items()),seconds=time.perf_counter()-start)
        return result
    compiled=compile_path(path,problem,{},source_atoms=(),complete_reference=False)
    result['constraint_encoding_seconds']=compiled.encoding_seconds
    stamp=time.perf_counter();result.update(tighten_supports(compiled))
    result['support_propagation_seconds']=time.perf_counter()-stamp
    stamp=time.perf_counter()
    objective,lower,metrics=event_objective(compiled,reverse=reverse,
        bond_floor=bond_floor,event_tolerance=event_tolerance)
    result.update(metrics,objective_encoding_seconds=time.perf_counter()-stamp)
    assert lower<=upper
    while lower<upper and time.perf_counter()<deadline:
        mid=(lower+upper)//2;solver=compiled.solver;solver.push();solver.add(objective<=mid)
        solver.set(timeout=max(1,int(1000*(deadline-time.perf_counter()))))
        status=solver.check();result['solver_queries']+=1
        if status==z3.sat:
            model=solver.model();realized=compiled.realize(model)
            witness=dict(realized['mapping'])
            events=bond_events(problem,witness,reverse=reverse,bond_floor=bond_floor,event_tolerance=event_tolerance)
            upper=events['total'];assert upper==model.eval(objective).as_long() and upper<=mid
        elif status==z3.unsat:lower=mid+1
        else:result['unknown_reason']=solver.reason_unknown()
        solver.pop()
        if status==z3.unknown:break
    result.update(lower_bound=lower,upper_bound=upper,optimal=lower==upper,events=events,
                  mapping=sorted(witness.items()),seconds=time.perf_counter()-start)
    if lower<upper and 'unknown_reason' not in result:result['unknown_reason']='time budget'
    return result

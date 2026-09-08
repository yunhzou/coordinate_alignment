"""Experimental independent fragment detection and symbolic compatible-cover assembly.

Sequential AAM is unchanged. Families are indivisible: no implicit trimming or
singleton completion. A negative result concerns only the supplied family pool.
"""
from dataclasses import dataclass
import time

from .fragment import FragmentMatchConfig, FragmentMatchContext, match_fragment
from .frag import build_graph
from .matcher import _nauty_orbits
from .alignment.branch import _Branch
from .search_graph import SearchContext, SearchGraphBuilder, frozen_value
from .search_symmetry import finalize_graph_symmetry
from .family_query import group_factors, query_path


@dataclass(frozen=True)
class IndependentDetection:
    graph: object
    attempts: tuple


def detect_independent(problem, *, seeds, cuts=(), branch_limit=100,
                       iso_tolerance=1.0, graph_floor=.2):
    """One cut, every requested seed from empty context; dedupe during admission."""
    source = build_graph(problem.reactant.elements, problem.reactant.wbo, graph_floor)
    target = build_graph(problem.product.elements, problem.product.wbo, graph_floor)
    source.remove_edges_from(cuts)
    context = FragmentMatchContext(source_orbits=_nauty_orbits(source, wbo_tol=iso_tolerance),
                                   target_orbits=_nauty_orbits(target, wbo_tol=iso_tolerance))
    config = FragmentMatchConfig(graph_floor=graph_floor, iso_tolerance=iso_tolerance,
                                 branch_limit=branch_limit)
    builder = SearchGraphBuilder(SearchContext(tuple(source), tuple(target), tuple(seeds),
        cuts=tuple(cuts), graph_floor=graph_floor, iso_tolerance=iso_tolerance,
        branch_limit=branch_limit))
    seen = set()
    attempts = []
    for position, seed in enumerate(seeds):
        result = match_fragment(source, target, seed=seed, context=context, config=config)
        attempts.append(dict(seed=seed, capped=result.capped, seconds=result.elapsed_seconds,
                             returned=len(result.matches)))
        builder.seed = seed
        builder.step = (0, position)
        if result.capped:
            root = _Branch(graph=builder)
            builder.stop(root, 'capped', stage='independent_fragment_growth',
                         count=result.branch_count, limit=result.branch_limit)
        for placement in result.matches:
            key = (tuple(sorted(placement.mapping.items())), tuple(sorted(placement.fragment)),
                   tuple(sorted(placement.deferred_edges)), frozen_value(placement.symmetry),
                   placement.preserved_bonds)
            if key in seen:
                continue
            seen.add(key)
            branch = _Branch(graph=builder)
            branch.commit(placement)
            builder.stop(branch, 'stalled')  # terminal of detection, not a complete AAM
    graph, _ = finalize_graph_symmetry(builder.finish(), target, iso_tolerance=iso_tolerance)
    return IndependentDetection(graph, tuple(attempts))


def independent_paths(graphs):
    """Deduplicate identical finalized relations, retaining one saved provenance.

    Different witnesses are not assumed equivalent. No symmetry expansion.
    """
    seen = set()
    for graph in graphs:
        for terminal in graph.terminals:
            path = next(graph.paths(terminal))
            placements = [graph.fragment_placement(e) for e in path.transitions
                          if graph.transitions[e].match is not None]
            assert len(placements) == 1
            edge = next(graph.transitions[e] for e in path.transitions if graph.transitions[e].match is not None)
            effective_cuts = tuple(cut for cut in path.context.cuts
                                   if all(r in path.mapping for r in cut))
            key = (frozen_value(edge.match), effective_cuts,
                   path.context.iso_tolerance, path.context.graph_floor)
            if key not in seen:
                seen.add(key)
                yield path


def _encode_family(path, prefix):
    """Compile ONE independent fragment's correlated target action to SMT.

    Uses the same exact group-factor representation as the existing family query.
    No endpoint-equivalence normalization and no independent per-atom orbit product.
    """
    import z3
    edge = next(e for e in path.transitions if path.graph.transitions[e].match is not None)
    fragment = path.graph.fragment_placement(edge)
    constraints = []
    serial = 0
    def symbol(domain):
        nonlocal serial
        domain = tuple(sorted(set(domain)))
        if len(domain) == 1:
            return (domain[0], frozenset(domain))
        v = z3.Int(f'{prefix}_{serial}'); serial += 1
        constraints.append(z3.Or(*(v == p for p in domain)))
        return v, frozenset(domain)
    def lookup(value, table):
        expression, support = value
        if isinstance(expression, int):
            return table.get(expression, value)
        images = [(p, table.get(p, (p, frozenset((p,))))) for p in support]
        out = images[-1][1][0]
        for p, v in reversed(images[:-1]):
            out = z3.If(expression == p, v[0], out)
        return out, frozenset().union(*(v[1] for _,v in images))
    values = {r:(p, frozenset((p,))) for r,p in path.mapping.items()}
    for domain in fragment.symmetry_domains:
        if domain.source == 'exact_automorph_group' or len(domain.p_atoms) < 2:
            continue
        pool = domain.p_atoms
        table = {p:symbol(pool) for p in pool}
        constraints.append(z3.Distinct(*(v[0] for v in table.values())))
        block_atoms = {r for d in fragment.symmetry_domains if d.source != 'exact_automorph_group' for r in d.r_atoms}
        fixed = {p for r,p in fragment.representative_assignments if r in fragment.exact_fixed or r not in block_atoms}
        constraints.extend(table[p][0] == p for p in fixed.intersection(pool))
        values = {r:lookup(v, table) for r,v in values.items()}
    for factor in group_factors(tuple(g.images for g in fragment.target_generators)):
        active = set().union(*(v[1] for v in values.values()))
        moved = {p for p in active if any(g[p] != p for g in factor)}
        if not moved:
            continue
        choice = symbol(range(len(factor)))
        table = {p:lookup(choice, {i:(g[p], frozenset((g[p],))) for i,g in enumerate(factor)}) for p in moved}
        values = {r:lookup(v, table) for r,v in values.items()}
    return values, constraints


def assemble_independent(problem, paths, *, seconds=60, required_mapping=None):
    """Joint compressed compatible cover. No ranking or solution enumeration.

    All atoms on the smaller side are mapped; all source atoms are required here.
    Caller must orient source no larger than target. A required mapping is solely
    an explicitly labelled diagnostic existence query, never used in detection.
    """
    import z3
    start = time.monotonic()
    nr, np = problem.source_atom_count, problem.target_atom_count
    if nr > np:
        raise ValueError('Orient the smaller endpoint as source for assembly')
    paths = tuple(paths)
    solver = z3.Solver()
    atom = [z3.Int(f'atom_{r}') for r in range(nr)]
    for r in range(nr):
        solver.add(z3.Or(*(atom[r] == p for p in range(np)
                          if problem.reactant.elements[r] == problem.product.elements[p])))
    solver.add(z3.Distinct(*atom))
    incidence = [[] for _ in range(nr)]
    selectors = []
    for i, path in enumerate(paths):
        if time.monotonic()-start >= seconds:
            return dict(status='unknown', reason='encoding watchdog', seconds=time.monotonic()-start)
        selected = z3.Bool(f'family_{i}')
        selectors.append(selected)
        values, constraints = _encode_family(path, f'f{i}')
        for r,(expression, support) in values.items():
            incidence[r].append(selected)
            constraints.append(atom[r] == expression)
        # Reject group/domain combinations that break the fragment's actual bonds.
        edge = next(e for e in path.transitions if path.graph.transitions[e].match is not None)
        placement = path.graph.fragment_placement(edge)
        cuts = {tuple(sorted(e)) for e in path.context.cuts} | set(placement.deferred_edges)
        for r in values:
            for s in values:
                w = problem.reactant.wbo[r,s]
                if r >= s or w < path.context.graph_floor or (r,s) in cuts:
                    continue
                vr, vs = values[r], values[s]
                constraints.append(z3.Or(*(z3.And(vr[0] == p, vs[0] == q)
                    for p in vr[1] for q in vs[1] if p != q
                    and problem.product.wbo[p,q] >= path.context.graph_floor
                    and abs(w-problem.product.wbo[p,q]) <= path.context.iso_tolerance+1e-9)))
        solver.add(z3.Implies(selected, z3.And(*constraints)))
    missing = [r for r,v in enumerate(incidence) if not v]
    if missing:
        return dict(status='not_covered', missing_source_atoms=missing, seconds=time.monotonic()-start)
    for choices in incidence:
        # Overlapping claims are legal if they agree on the same global mapping.
        # Global Distinct still forbids different source atoms sharing a target.
        solver.add(z3.Or(*choices))
    for r,p in (required_mapping or {}).items():
        solver.add(atom[int(r)] == int(p))
    solver.set(timeout=max(1, int(1000*(seconds-(time.monotonic()-start)))))
    status = solver.check()
    if status != z3.sat:
        return dict(status='not_covered' if status == z3.unsat else 'unknown',
                    reason=str(status), seconds=time.monotonic()-start)
    model = solver.model()
    mapping = {r:model.eval(v).as_long() for r,v in enumerate(atom)}
    selected = [i for i,v in enumerate(selectors) if z3.is_true(model.eval(v))]
    # Independent certification through the existing AAM family-query code.
    for i in selected:
        partial = {r:mapping[r] for r in paths[i].mapping}
        left = seconds-(time.monotonic()-start)
        if left <= 0:
            return dict(status='unknown', reason='certificate watchdog', seconds=time.monotonic()-start)
        certificate, _ = query_path(paths[i], problem, partial, source_atoms=tuple(partial),
                                     timeout_ms=max(1, int(left*1000)))
        if certificate != 'recovered':
            if certificate == 'not_recovered':
                raise AssertionError('Assembly encoder disagrees with existing family query')
            return dict(status='unknown', reason='certificate timeout', seconds=time.monotonic()-start)
    return dict(status='covered', mapping=sorted(mapping.items()), selected=selected,
                seconds=time.monotonic()-start, reference_constrained=required_mapping is not None)

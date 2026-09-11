"""Raw-WBO delta events, separate from the graph used for fragment matching."""
from collections import defaultdict
from functools import cached_property, lru_cache
import copy
import hashlib
import time

import numpy as np
import pynauty

from rxn_core.frag import is_metal_element, classify_bonds
from golden_evaluation import colored_graph, exact_action


def binary_metal_input(raw, threshold=.2):
    """Copy input, changing only metal-pair search weights; never edit raw WBO."""
    out = copy.deepcopy(raw)
    for side in ('reactant', 'product'):
        endpoint = out[side]
        w = np.asarray(endpoint['wbo'], dtype=float)
        metal = np.array([is_metal_element(e) for e in endpoint['elements']])
        mask = metal[:, None] | metal[None, :]
        w[mask] = (w[mask] >= threshold).astype(float)
        np.fill_diagonal(w, 0.)
        endpoint['wbo'] = w.tolist()
    return out


def scalar_events(raw, mapping):
    r, p = raw['reactant'], raw['product']
    broken, formed, _, _ = classify_bonds(dict(mapping), np.asarray(r['wbo']), np.asarray(p['wbo']),
        dwbo_threshold=.5, elements_R=r['elements'], elements_P=p['elements'], metal_dwbo_threshold=.3)
    inverse = {p: r for r, p in dict(mapping).items()}
    return dict(broken=broken, formed=formed, total=len(broken)+len(formed),
                formed_on_R=[tuple(sorted((inverse[a], inverse[b]))) for a,b,*_ in formed
                             if a in inverse and b in inverse])


class DeltaPatterns:
    """Signed delta events modulo exact symmetry of the score-response graph.

    Every pair, including sub-floor contacts, enters the response labels.
    A generator can swap atoms only if all possible raw-WBO event responses
    remain identical. No graph-floor test substitutes for the delta threshold.
    """
    def __init__(self, raw):
        self.raw = raw
        self.r, self.p = [np.asarray(raw[s]['wbo']) for s in ('reactant', 'product')]
        self.elements = raw['reactant']['elements']
        self.n = len(self.r)
        assert len(self.p) == self.n
        self.a, self.b = np.triu_indices(self.n, 1)
        metal = np.array([is_metal_element(e) for e in self.elements])
        self.tau = np.where(metal[self.a] | metal[self.b], .3, .5)
        self.cache = {}

    def counts(self, vectors):
        vectors = np.asarray(vectors, dtype=int)
        delta = self.p[vectors[:, self.a], vectors[:, self.b]] - self.r[self.a, self.b]
        return np.column_stack(((delta <= -self.tau).sum(axis=1), (delta >= self.tau).sum(axis=1)))

    @cached_property
    def features(self):
        features = []
        for side, other in (('reactant', 'product'), ('product', 'reactant')):
            ep, op = self.raw[side], self.raw[other]
            w, v = np.asarray(ep['wbo']), np.asarray(op['wbo'])
            values = defaultdict(set)
            for a in range(self.n):
                for b in range(a+1, self.n):
                    values[tuple(sorted((op['elements'][a], op['elements'][b])))].add(float(v[a,b]))
            values = {k: np.array(sorted(v)) for k,v in values.items()}
            bonds = []
            for a in range(self.n):
                for b in range(a+1, self.n):
                    pair = tuple(sorted((ep['elements'][a], ep['elements'][b])))
                    tau = .3 if any(is_metal_element(e) for e in pair) else .5
                    delta = values[pair] - w[a,b]
                    response = np.where(delta <= -tau, -1, np.where(delta >= tau, 1, 0)).astype(np.int8)
                    bonds.append((a,b,(response.tobytes().hex(),)))
            features.append(dict(colors=[(e,) for e in ep['elements']], bonds=bonds))
        return features

    @cached_property
    def base(self):
        return colored_graph([self.features[0]])

    @cached_property
    def generators(self):
        return tuple(tuple(g[:self.n]) for g in pynauty.autgrp(self.base)[0])

    def describe(self, vector):
        vector = list(map(int, vector))
        assert sorted(vector) == list(range(self.n))
        assert all(e == self.raw['product']['elements'][vector[a]] for a,e in enumerate(self.elements))
        delta = self.p[np.ix_(vector, vector)][self.a,self.b] - self.r[self.a,self.b]
        events = tuple(tuple((int(a), int(b)) for a,b in zip(self.a[mask], self.b[mask]))
                       for mask in (delta <= -self.tau, delta >= self.tau))
        if events not in self.cache:
            base = self.base
            adjacency = {a:set(ns) for a,ns in base.adjacency_dict.items()}
            colors = [set(c) for c in base.vertex_coloring]
            n = base.number_of_vertices
            for edges in events:
                marker = n; n += 1
                adjacency[marker] = set(); colors.append({marker}); gadgets = set()
                for a,b in edges:
                    adjacency[n] = {a,b,marker}
                    for v in (a,b,marker): adjacency[v].add(n)
                    gadgets.add(n); n += 1
                if gadgets: colors.append(gadgets)
            graph = pynauty.Graph(n, adjacency_dict={a:sorted(adjacency.get(a,())) for a in range(n)}, vertex_coloring=colors)
            self.cache[events] = hashlib.sha256(pynauty.certificate(graph)).hexdigest()
        result = dict(id=self.cache[events], mapping=vector, events=dict(zip(('broken','formed'), events)),
                      counts=dict(zip(('broken','formed'), map(len,events))), total=sum(map(len,events)))
        scalar = scalar_events(self.raw, dict(enumerate(vector)))
        assert scalar['total'] == result['total']
        assert {(a,b) for a,b,*_ in scalar['broken']} == set(events[0])
        assert set(scalar['formed_on_R']) == set(events[1])
        return result


def delta_objective(compiled, canonical):
    """Score symbolic mappings with raw WBOs; leave search constraints intact."""
    import z3
    constant = 0; lower = 0; terms = []
    for a,b,tau in zip(canonical.a,canonical.b,canonical.tau):
        va,vb = compiled.values[a],compiled.values[b]
        groups = defaultdict(list)
        for x in va[1]:
            for y in vb[1]:
                if x==y: continue
                delta = canonical.p[x,y]-canonical.r[a,b]
                cost = int(delta <= -tau or delta >= tau)
                groups[cost].append((x,y))
        assert groups
        lower += min(groups)
        if len(groups)==1:
            constant += next(iter(groups)); continue
        default = max(groups,key=lambda k:len(groups[k])); term = default
        for cost,pairs in groups.items():
            if cost==default: continue
            rows = defaultdict(list)
            for x,y in pairs: rows[x].append(y)
            condition = z3.Or(*(z3.And(va[0]==x,z3.Or(*(vb[0]==y for y in ys))) for x,ys in rows.items()))
            term = z3.If(condition,cost,term)
        terms.append(term)
    return z3.IntVal(constant)+z3.Sum(terms), lower


def query_pattern(compiled, canonical, pattern, objective, timeout_ms):
    import z3
    from rxn_core.family_query import SymbolicActions
    started = time.perf_counter(); solver = compiled.solver; solver.push()
    solver.add(objective==pattern['total'])
    encoder = SymbolicActions(solver,fresh=True)
    atoms = sorted({a for edges in pattern['events'].values() for edge in edges for a in edge})
    normalized = encoder.act([encoder.constant(a) for a in atoms],canonical.generators,'delta_event_equivalence')
    mapped = {a:encoder.lookup(v,dict(enumerate(compiled.values))) for a,v in zip(atoms,normalized)}
    for kind,edges in pattern['events'].items():
        for a,b in edges:
            tau = .3 if is_metal_element(canonical.elements[a]) or is_metal_element(canonical.elements[b]) else .5
            va,vb = mapped[a],mapped[b]; rows = defaultdict(list)
            for x in va[1]:
                for y in vb[1]:
                    delta = canonical.p[x,y]-canonical.r[a,b]
                    if x!=y and (delta <= -tau if kind=='broken' else delta >= tau): rows[x].append(y)
            solver.add(z3.Or(*(z3.And(va[0]==x,z3.Or(*(vb[0]==y for y in ys))) for x,ys in rows.items())))
    remaining = timeout_ms-int(1000*(time.perf_counter()-started))
    if remaining<=0: solver.pop(); return 'unknown',None
    solver.set(timeout=remaining); status=solver.check(); witness=None
    if status==z3.sat:
        mapping=dict(compiled.realize(solver.model())['mapping'])
        witness=canonical.describe([mapping[a] for a in range(canonical.n)])
        assert witness['id']==pattern['id'] and witness['total']==pattern['total']
    solver.pop()
    return ('represented' if status==z3.sat else 'absent' if status==z3.unsat else 'unknown'),witness


def membership(aam, canonical, patterns, saved_patterns, seconds=120):
    """Query supplied event patterns in saved families, with explicit unknowns."""
    from rxn_core.family_query import compile_path
    from rxn_core.search_graph import frozen_value
    start=time.perf_counter(); deadline=start+seconds
    results={k:dict(status='represented',method='saved_representative',witness=saved_patterns[k])
             for k in patterns if k in saved_patterns}
    pending={k:v for k,v in patterns.items() if k not in results}
    graph=aam.graph; unknown=set(); seen=set(); queries=encoded=scanned=0; exhausted=True
    @lru_cache(None)
    def invariant_action(images): return exact_action(images,canonical.features[1])
    @lru_cache(None)
    def invariant_edge(edge):
        placement=graph.fragment_placement(edge)
        if placement is None:return True
        if any(not invariant_action(tuple(g.images)) for g in placement.target_generators):return False
        for domain in placement.symmetry_domains:
            if domain.source=='exact_automorph_group':continue
            for a,b in zip(domain.p_atoms,domain.p_atoms[1:]):
                images=list(range(canonical.n)); images[a],images[b]=b,a
                if not invariant_action(tuple(images)):return False
        return True
    incoming=[[] for _ in graph.states]
    for edge in graph.transitions:incoming[edge.target].append(edge)
    @lru_cache(None)
    def can_change(state):return any(can_change(e.source) or not invariant_edge(e.id) for e in incoming[state])
    for terminal in graph.terminals:
        if not pending:break
        if time.perf_counter()>deadline:exhausted=False;break
        if len(graph.states[terminal].mapping)!=canonical.n or not can_change(terminal):continue
        for path in graph.paths(terminal):
            if time.perf_counter()>deadline:exhausted=False;break
            scanned+=1
            if all(invariant_edge(e) for e in path.transitions):continue
            key=(tuple(sorted(path.mapping.items())),path.context.cuts,frozen_value(path.fragments))
            if key in seen:continue
            seen.add(key)
            compiled=compile_path(path,aam.problem,{},source_atoms=(),complete_reference=False)
            objective,lower=delta_objective(compiled,canonical); encoded+=1
            for key,pattern in list(pending.items()):
                if lower>pattern['total']:continue
                remaining=int(1000*(deadline-time.perf_counter()))
                if remaining<=0:exhausted=False;break
                status,witness=query_pattern(compiled,canonical,pattern,objective,min(3000,remaining)); queries+=1
                if status=='represented':
                    results[key]=dict(status=status,method='compressed_family',terminal=terminal,
                        transitions=list(path.transitions),witness=witness); del pending[key]
                elif status=='unknown':unknown.add(key)
            if not pending or not exhausted:break
        if not exhausted:break
    for key in pending:
        results[key]=dict(status='excluded_from_saved_families' if exhausted and key not in unknown else 'unresolved')
    return dict(results=results,queries=queries,encoded=encoded,paths_scanned=scanned,exhausted=exhausted,
                elapsed=time.perf_counter()-start,scope='Supplied patterns only; no enumeration of all unseen patterns or global optimum claim.')

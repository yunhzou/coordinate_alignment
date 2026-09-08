"""Resumable extraction of exact-symmetry mapping classes from a path.

The search family stays compressed. Each exclusion removes a whole exact
endpoint-symmetry orbit using quantified, factored actions, not one bijection.
Certificates include the full explicit-atom mapping and unmatched context.
"""
from collections import defaultdict
import hashlib
from itertools import combinations
import time

import pynauty

from .family_query import SymbolicActions,compile_path
from .family_scoring import bond_events,validate_representative
from .search_graph import frozen_value


class PatternEquivalence:
    """Exact supplied endpoint graphs; optional chemical tags refine identity.

    No tolerance or resonance merging is applied. Callers with formal-charge,
    isotope or stereochemical annotations should pass them as atom/bond tags.
    """
    def __init__(self,problem,*,atom_tags=None,bond_tags=None,bond_floor=.2):
        endpoints=(problem.reactant,problem.product)
        atom_tags=atom_tags if atom_tags is not None else ({},{})
        bond_tags=bond_tags if bond_tags is not None else ({},{})
        self.problem=problem;self.atoms=[];self.bonds=[]
        for side,e in enumerate(endpoints):
            self.atoms.append(tuple((e.elements[i],frozen_value(atom_tags[side].get(i,())))
                                    for i in range(e.atom_count)))
            self.bonds.append({(a,b):(float(e.wbo[a,b]),frozen_value(bond_tags[side].get((a,b),())))
                               for a,b in combinations(range(e.atom_count),2) if e.wbo[a,b]>bond_floor})
        self.generators=tuple(tuple(tuple(g[:len(self.atoms[side])])
            for g in pynauty.autgrp(self._graph((side,)))[0]) for side in (0,1))
        self._invariance={}

    def _graph(self,sides=(0,1),mapping=None):
        adjacency=defaultdict(set);colors=defaultdict(set);offsets={};count=0
        for side in sides:
            offsets[side]=count
            for label in self.atoms[side]:colors[('atom',side,label)].add(count);count+=1
        def link(a,b,label):
            nonlocal count
            adjacency[count].update((a,b));adjacency[a].add(count);adjacency[b].add(count)
            colors[label].add(count);count+=1
        for side in sides:
            for (a,b),label in self.bonds[side].items():
                link(offsets[side]+a,offsets[side]+b,('bond',side,label))
        if mapping is not None:
            for a,b in sorted(dict(mapping).items()):link(offsets[0]+a,offsets[1]+b,('mapping',))
        return pynauty.Graph(count,adjacency_dict={i:sorted(adjacency[i]) for i in range(count)},
            vertex_coloring=[colors[k] for k in sorted(colors,key=repr)])

    def key(self,mapping):
        return hashlib.sha256(pynauty.certificate(self._graph(mapping=mapping))).hexdigest()

    def preserves_target(self,g):
        g=tuple(g)
        if g not in self._invariance:
            atoms,bonds=self.atoms[1],self.bonds[1]
            self._invariance[g]=(all(atoms[a]==atoms[g[a]] for a in range(len(atoms))) and
                all(bonds.get(tuple(sorted((g[a],g[b]))))==label for (a,b),label in bonds.items()))
        return self._invariance[g]

    def invariant_path(self,path):
        for edge in path.transitions:
            placement=path.graph.fragment_placement(edge)
            if placement is None:continue
            if placement.target_generators is None:raise ValueError('Finalize path symmetry first')
            if any(not self.preserves_target(g.images) for g in placement.target_generators):return False
            for domain in placement.symmetry_domains:
                if domain.source=='exact_automorph_group':continue
                for a,b in zip(domain.p_atoms,domain.p_atoms[1:]):
                    g=list(range(len(self.atoms[1])));g[a],g[b]=g[b],g[a]
                    if not self.preserves_target(g):return False
        return True

    def orbit_membership(self,values,mapping):
        """Existence of an exact endpoint action, expressed without orbit lists."""
        import z3
        solver=z3.Solver();encoder=SymbolicActions(solver,fresh=True)
        nr=len(self.atoms[0]);sentinel=len(self.atoms[1]);mapping=dict(mapping)
        source=encoder.act([encoder.constant(r) for r in range(nr)],self.generators[0],'source')
        normalized=[encoder.lookup(v,{r:encoder.constant(mapping.get(r,sentinel)) for r in range(nr)})
                    for v in source]
        normalized=encoder.act(normalized,self.generators[1],'target')
        equality=z3.And(*solver.assertions(),*(x[0]==y[0] for x,y in zip(values,normalized)))
        return z3.Exists(encoder.variables,equality) if encoder.variables else equality


def extract_path_patterns(path,problem,equivalence,*,previous=(),seconds=10.,reverse=False,
                          on_pattern=None):
    """Keep one physical witness per exact full-mapping class, with resume state.

    A complete flag is issued only by an invariance certificate or UNSAT after
    class-level exclusions. Timeouts preserve all earlier candidates. Solver
    timeouts are soft; callers can also impose a process watchdog.
    """
    import z3
    start=time.perf_counter();deadline=start+seconds
    validate_representative(path,problem)
    patterns={p['key']:p for p in previous};queries=0
    def add(mapping,actions):
        key=equivalence.key(mapping)
        if key in patterns:return False
        record=dict(key=key,mapping=sorted(dict(mapping).items()),actions=actions,
            events=bond_events(problem,mapping,reverse=reverse),
            terminal=path.terminal,transitions=list(path.transitions))
        patterns[key]=record
        if on_pattern is not None:on_pattern(record)
        return True
    add(path.mapping,[])
    if equivalence.invariant_path(path):
        return dict(patterns=list(patterns.values()),complete=True,reason='exact_pattern_invariance',
                    solver_queries=0,seconds=time.perf_counter()-start)
    compiled=compile_path(path,problem,{},source_atoms=(),complete_reference=False)
    for record in patterns.values():
        compiled.solver.add(z3.Not(equivalence.orbit_membership(compiled.values,record['mapping'])))
    reason='time_budget';complete=False
    while time.perf_counter()<deadline:
        compiled.solver.set(timeout=max(1,int(1000*(deadline-time.perf_counter()))))
        status=compiled.solver.check();queries+=1
        if status==z3.unsat:complete=True;reason='all_exact_pattern_classes_excluded';break
        if status==z3.unknown:reason=compiled.solver.reason_unknown();break
        witness=compiled.realize(compiled.solver.model());mapping=dict(witness['mapping'])
        assert add(mapping,witness['actions']), 'Orbit exclusion must eliminate the whole previous class'
        compiled.solver.add(z3.Not(equivalence.orbit_membership(compiled.values,mapping)))
    return dict(patterns=list(patterns.values()),complete=complete,reason=reason,
                solver_queries=queries,seconds=time.perf_counter()-start,
                constraint_encoding_seconds=compiled.encoding_seconds)

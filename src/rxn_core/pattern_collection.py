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
        self.twins=tuple(self._twins(side) for side in (0,1))
        self.block_of=tuple({a:i for i,block in enumerate(blocks) for a in block}
                            for blocks in self.twins)
        self.quotient_generators=tuple(tuple(sorted({
            tuple(self.block_of[side][g[block[0]]] for block in self.twins[side])
            for g in self.generators[side]}-{tuple(range(len(self.twins[side])))})) for side in (0,1))
        self.block_orbits=tuple(self._orbits(len(self.twins[side]),self.quotient_generators[side])
                                for side in (0,1))
        self._invariance={}

    @staticmethod
    def _orbits(n,generators):
        groups=[{i} for i in range(n)]
        for g in generators:
            for a,b in enumerate(g):
                merged=groups[a]|groups[b]
                for i in merged:groups[i]=merged
        return tuple(frozenset(g) for g in groups)

    def _twins(self,side):
        """Maximal exact interchangeable blocks (including attached H groups).

        Swapping two members fixes every other atom. The quotient therefore
        removes only a full symmetric-group kernel, not correlated ring moves.
        """
        atoms,bonds=self.atoms[side],self.bonds[side];blocks=[]
        def edge(a,b):return bonds.get(tuple(sorted((a,b))))
        for a in range(len(atoms)):
            for block in blocks:
                b=block[0]
                if atoms[a]==atoms[b] and all(edge(a,k)==edge(b,k)
                    for k in range(len(atoms)) if k not in (a,b)):
                    block.append(a);break
            else:blocks.append([a])
        return tuple(tuple(b) for b in blocks)

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
        """Exact orbit membership on a twin-block contingency matrix.

        Counts fully characterize an injective mapping modulo independent
        within-block permutations. Only the remaining correlated block actions
        are quantified; factorial hydrogen permutations never reach the solver.
        """
        import z3
        solver=z3.Solver();encoder=SymbolicActions(solver,fresh=True)
        mapping=dict(mapping);nr,np_=map(len,self.twins)
        source=encoder.act([encoder.constant(r) for r in range(nr)],self.quotient_generators[0],'source')
        target=encoder.act([encoder.constant(p) for p in range(np_)],self.quotient_generators[1],'target')
        counts=defaultdict(int)
        for r,p in mapping.items():counts[self.block_of[0][r],self.block_of[1][p]]+=1
        cells=set()
        for r,(_,support) in enumerate(values):
            cells.update((self.block_of[0][r],self.block_of[1][p])
                         for p in support if p in self.block_of[1])
        for r,p in counts:
            cells.update((a,b) for a in self.block_orbits[0][r] for b in self.block_orbits[1][p])
        conditions=[]
        for r,p in sorted(cells):
            block,pblock=self.twins[0][r],self.twins[1][p]
            terms=[]
            for a in block:
                expression,support=values[a];allowed=support.intersection(pblock)
                if not allowed:continue
                terms.append(1 if support.issubset(pblock) else z3.If(z3.Or(
                    *(expression==b for b in allowed)),1,0))
            rows={i:encoder.lookup(target[p],{j:encoder.constant(counts[i,j])
                  for j in target[p][1]}) for i in source[r][1]}
            expected=encoder.lookup(source[r],rows)[0]
            actual=z3.Sum(terms) if terms else 0
            conditions.append(actual==expected)
        equality=z3.And(*solver.assertions(),*conditions)
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

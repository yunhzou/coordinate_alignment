"""Frozen dense reference for differential tests only; never used by AAM runners.

Derived from metal_binary_events.py at 7cc3be4. Kept to make the replacement
independently falsifiable, including historical certificate compatibility.
"""
from collections import defaultdict
from functools import cached_property
import hashlib
import numpy as np
import pynauty
from rxn_core.frag import is_metal_element
from golden_evaluation import colored_graph
from metal_binary_events import scalar_events

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


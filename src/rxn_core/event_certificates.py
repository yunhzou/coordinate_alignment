"""Exact event-class certificates over the recorded, possibly looser, search family.

Reachable-image sets are rejection-only relaxations. Positive class-invariance
certificates use exact response automorphisms and preserve action order. Neither
operation expands group elements or complete mappings.
"""
from collections import OrderedDict

import numpy as np


class EventFamilyCertificates:
    def __init__(self, index, cache_size=16384):
        self.index = index
        self.cache_size = cache_size
        self._graph = None
        self._edges = OrderedDict()
        self._orbits = OrderedDict()
        self._events = OrderedDict()
        self._source_actions = OrderedDict()

    def _remember(self, cache, key, value):
        cache[key] = value
        cache.move_to_end(key)
        if len(cache) > self.cache_size:
            cache.popitem(last=False)
        return value

    def actions(self, path):
        """Actions in application order, including the recorded pool locks.

        Search graphs are treated as immutable, as in the archive/query API.
        Edge metadata is cached only while using the same graph object.
        """
        if self._graph is not path.graph:
            self._graph = path.graph
            self._edges.clear()
        actions = []
        for edge in reversed(path.transitions):
            if edge not in self._edges:
                placement = path.graph.fragment_placement(edge)
                row = []
                if placement is not None:
                    if placement.target_generators is None:
                        raise ValueError('finalize saved path symmetry first')
                    for domain in placement.symmetry_domains:
                        if domain.source == 'exact_automorph_group' or len(domain.p_atoms) < 2:
                            continue
                        locked = dict(path.graph.states[path.graph.transitions[edge].source].mapping)
                        block_atoms = {r for block in placement.symmetry_domains
                                       if block.source != 'exact_automorph_group' for r in block.r_atoms}
                        assignments = dict(placement.representative_assignments)
                        fixed = (set(locked.values()) | {assignments[r] for r in placement.exact_fixed}
                                 | {p for r, p in placement.representative_assignments if r not in block_atoms})
                        row.append(('pool', tuple(domain.p_atoms), tuple(sorted(fixed.intersection(domain.p_atoms)))))
                    generators = tuple(tuple(g.images) for g in placement.target_generators)
                    if generators:
                        row.append(('group', generators))
                self._remember(self._edges, edge, tuple(row))
            actions.extend(self._edges[edge])
        return tuple(actions)

    def group_orbits(self, generators):
        if generators not in self._orbits:
            parent = list(range(self.index.n))

            def root(a):
                while parent[a] != a:
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                return a

            for generator in generators:
                for a, b in enumerate(generator):
                    ra, rb = root(a), root(b)
                    if ra != rb:
                        parent[ra] = rb
            groups = {}
            for a in range(self.index.n):
                groups.setdefault(root(a), set()).add(a)
            self._remember(self._orbits, generators,
                           tuple(frozenset(groups[root(a)]) for a in range(self.index.n)))
        return self._orbits[generators]

    def lower_bound(self, vector, actions, stop_above):
        """Count only events unavoidable across overapproximated atom images.

        Correlations and matching constraints are relaxed. Impossible diagonal
        image pairs are also allowed, which can only weaken this lower bound.
        Pairs without a representative event can be omitted safely.
        """
        index = self.index
        vector = tuple(vector)
        if vector not in self._events:
            values = np.asarray(vector)
            delta = index.p[values[index.a], values[index.b]] - index.r[index.a, index.b]
            self._remember(self._events, vector,
                           tuple(np.flatnonzero((delta <= -index.tau) | (delta >= index.tau))))
        pairs = self._events[vector]
        if len(pairs) <= stop_above:
            return 0
        supports = [frozenset((p,)) for p in vector]
        for action in actions:
            if action[0] == 'pool':
                free = frozenset(set(action[1]) - set(action[2]))
                supports = [(s - free) | free if s & free else s for s in supports]
            else:
                orbits = self.group_orbits(action[1])
                supports = [frozenset().union(*(orbits[p] for p in s)) for s in supports]
        lower = 0
        for pair in pairs:
            a, b = int(index.a[pair]), int(index.b[pair])
            xs, ys = tuple(supports[a]), tuple(supports[b])
            delta = index.p[np.ix_(xs, ys)] - index.r[a, b]
            tau = index.tau[pair]
            if np.all((delta <= -tau) | (delta >= tau)):
                lower += 1
                if lower > stop_above:
                    break
        return lower

    def invariant_class(self, vector, actions):
        """Certify actions factor as target symmetry · m · source symmetry.

        An inner prefix may be conjugate through m to exact source response
        automorphisms. Remaining outer actions must preserve target responses.
        Arbitrary interleaving is unsafe because the groups need not commute.
        """
        index = self.index
        inverse = [0] * index.n
        for a, p in enumerate(vector):
            inverse[p] = a

        def source_preserves(generator):
            images = tuple(inverse[generator[p]] for p in vector)
            if images not in self._source_actions:
                labels = index.response_data[0]['labels']
                preserves = (all(index.elements[a] == index.elements[b] for a, b in enumerate(images))
                             and np.array_equal(labels, labels[np.ix_(images, images)]))
                self._remember(self._source_actions, images, preserves)
            return self._source_actions[images]

        outer_target = False
        for action in actions:
            if action[0] == 'group':
                generators = action[1]
            else:
                free = sorted(set(action[1]) - set(action[2]))
                generators = []
                for a, b in zip(free, free[1:]):
                    images = list(range(index.n))
                    images[a], images[b] = b, a
                    generators.append(tuple(images))
            if not outer_target and all(source_preserves(g) for g in generators):
                continue
            if all(index.preserves_target_action(g) for g in generators):
                outer_target = True
            else:
                return False
        return True

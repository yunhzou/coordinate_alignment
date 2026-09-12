"""Unordered final fragment pairs with a lossless union of mapping families.

Branch identity excludes growth history. A family retains the allowed action
program and the final preservation constraints, never a chosen history. Only
provably commuting actions are reordered. Intrinsic pair symmetry is rebuilt
from matching-response labels and is shared independently of growth prefixes;
it does not silently enlarge the saved mapping relation.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import combinations
import heapq
import hashlib
import json


def _problem_fingerprint(problem):
    import numpy as np
    digest = hashlib.sha256()
    for endpoint in (problem.reactant, problem.product):
        digest.update(json.dumps(endpoint.elements, separators=(',', ':')).encode())
        digest.update(np.asarray(endpoint.wbo, dtype='<f8').tobytes())
    return digest.hexdigest()


def final_fragment_pairs(state):
    mapping = dict(state.mapping)
    owners = dict(state.islands)
    groups = defaultdict(list)
    for r in sorted(mapping):
        # Unlabelled anchors are singleton fragments.
        groups[('island', owners[r]) if r in owners else ('single', r)].append(r)
    return tuple(sorted((tuple(atoms), tuple(sorted(mapping[r] for r in atoms)))
                        for atoms in groups.values()))


def _support(action):
    kind, data = action
    if kind == 'pool':
        return frozenset(data)
    return frozenset(a for g in data for a, b in enumerate(g) if a != b)


def canonical_actions(actions):
    """Lexical normal form of independent actions; dependent order is retained."""
    actions = tuple(actions)
    support = tuple(_support(a) for a in actions)
    outgoing = [[] for _ in actions]
    incoming = [0] * len(actions)
    for j in range(len(actions)):
        for i in range(j):
            if support[i] & support[j]:
                outgoing[i].append(j)
                incoming[j] += 1
    ready = [(actions[i], i) for i, degree in enumerate(incoming) if not degree]
    heapq.heapify(ready)
    result = []
    while ready:
        action, i = heapq.heappop(ready)
        # A full subgroup followed immediately by itself is the same subgroup.
        if not result or result[-1] != action:
            result.append(action)
        for j in outgoing[i]:
            incoming[j] -= 1
            if not incoming[j]:
                heapq.heappush(ready, (actions[j], j))
    return tuple(result)


@dataclass
class FinalFamily:
    mapping: tuple
    required_edges: tuple
    actions: tuple
    policy: tuple
    provenance: list = field(default_factory=list)

    def validate_representative(self, problem):
        mapping = dict(self.mapping)
        assert set(mapping) == set(range(problem.source_atom_count))
        assert set(mapping.values()) == set(range(problem.target_atom_count))
        assert all(problem.reactant.elements[r] == problem.product.elements[p] for r, p in mapping.items())
        floor, iso = self.policy
        assert all(problem.product.wbo[mapping[a], mapping[b]] >= floor and
                   abs(problem.reactant.wbo[a, b] - problem.product.wbo[mapping[a], mapping[b]]) <= iso + 1e-9
                   for a, b in self.required_edges)

    def as_path(self, problem, source_edges=None):
        """A flat equivalent family for the existing exact decoder/query API."""
        from .search_graph import (AAMSearchGraph, SearchContext, SearchState,
                                   FragmentTransition, SearchStop)
        n = problem.source_atom_count
        floor, iso = self.policy
        required = set(self.required_edges)
        if source_edges is None:
            source_edges = tuple((a, b) for a, b in combinations(range(n), 2)
                                 if problem.reactant.wbo[a, b] >= floor)
        cuts = tuple(pair for pair in source_edges if pair not in required)
        context = SearchContext(tuple(range(n)), tuple(range(n)), (), cuts=cuts,
                                graph_floor=floor, iso_tolerance=iso,
                                objective='final_family_postprocessing')
        def match(atoms, blocks=(), generators=(), witness=()):
            return dict(fragment=atoms, deferred_edges=(), symmetry=dict(
                witness=dict(witness), blocks=blocks, exact_fixed=[],
                automorph_blocks=[], automorph_generators=generators))
        records = [match(tuple(range(n)), witness=self.mapping)]
        for kind, data in reversed(self.actions):
            records.append(match((), blocks=(dict(r_atoms=(), p_atoms=data),))
                           if kind == 'pool' else match((), generators=data))
        states = [SearchState(i, 0, self.mapping if i == len(records) else (),
                              tuple((r, 0) for r, _ in self.mapping) if i == len(records) else (), ())
                  for i in range(len(records) + 1)]
        edges = tuple(FragmentTransition(i, i, i + 1, None, (i, 0), record,
                                         self.required_edges if i == 0 else ())
                      for i, record in enumerate(records))
        graph = AAMSearchGraph((context,), (0,), tuple(states), edges,
                              (SearchStop(len(records), 'objective_met'),))
        return next(graph.paths())

    @classmethod
    def from_record(cls, record):
        actions = tuple((kind, tuple(data) if kind == 'pool' else tuple(tuple(g) for g in data))
                        for kind, data in record['actions'])
        return cls(tuple(tuple(p) for p in record['mapping']),
                   tuple(tuple(p) for p in record['required_edges']), actions,
                   tuple(record['policy']), list(record.get('provenance', ())))

    def to_record(self):
        return dict(mapping=self.mapping, required_edges=self.required_edges,
                    actions=self.actions, policy=self.policy, provenance=self.provenance)


class FinalBranchCatalogue:
    """One problem, many saved searches; all histories remain as references."""
    schema = 'rxn_core.final_fragment_branches/v1'

    def __init__(self, problem):
        if (problem.source_atom_count != problem.target_atom_count or
                Counter(problem.reactant.elements) != Counter(problem.product.elements)):
            raise ValueError('Final-family postprocessing requires balanced endpoints')
        self.problem = problem
        self.problem_sha256 = _problem_fingerprint(problem)
        self.branches = {}
        self.coupled_branches = {}
        self.families = []
        self._family_ids = {}
        self._required = {}
        self._programs = {}
        self._factors = {}
        self.symmetries = {}
        self.path_count = 0
        self.incomplete_paths = 0
        self.mode = 'lossless_saved_family_union'
        self.reconstruction = {}
        self._labels = {}

    def _intern(self, cache, value):
        return cache.setdefault(value, value)

    def add_aam(self, aam, archive='in_memory'):
        if _problem_fingerprint(aam.problem) != self.problem_sha256:
            raise ValueError('Cannot combine different endpoint problems')
        return self.add_graph(aam.graph, archive)

    def add_graph(self, graph, archive='in_memory'):
        n = self.problem.source_atom_count
        constraints, actions, pairs = {}, {}, {}
        for path in graph.paths():
            state = graph.states[path.terminal]
            if len(state.mapping) != n:
                self.incomplete_paths += 1
                continue
            self.path_count += 1
            context = path.context
            policy = context.graph_floor, context.iso_tolerance
            if path.terminal not in pairs:
                pairs[path.terminal] = final_fragment_pairs(state)
            pair_key = policy, pairs[path.terminal]
            branch = self.branches.setdefault(pair_key, set())
            required = set()
            program = []
            for eid in path.transitions:
                edge = graph.transitions[eid]
                if edge.match is None:
                    continue
                if eid not in constraints:
                    raw = edge.match
                    ignored = {tuple(sorted(e)) for e in context.cuts}
                    ignored.update(tuple(sorted(e)) for e in raw.get('deferred_edges', ()))
                    atoms = sorted(raw['fragment'])
                    constraints[eid] = tuple((a, b) for a, b in combinations(atoms, 2)
                                             if (a, b) not in ignored and
                                             self.problem.reactant.wbo[a, b] >= context.graph_floor)
                required.update(constraints[eid])
            for eid in reversed(path.transitions):
                edge = graph.transitions[eid]
                if edge.match is None:
                    continue
                if eid not in actions:
                    raw = edge.match['symmetry']
                    if raw.get('automorph_generators') is None:
                        raise ValueError('Finalize saved symmetry before final-branch deduplication')
                    blocks = [b for b in raw.get('blocks', ())
                              if b.get('source') != 'exact_automorph_group']
                    block_atoms = {r for b in blocks for r in b.get('r_atoms', ())}
                    witness = dict(raw.get('witness', {}))
                    fixed = set(dict(graph.states[edge.source].mapping).values())
                    fixed.update(witness[r] for r in raw.get('exact_fixed', ()))
                    fixed.update(p for r, p in witness.items() if r not in block_atoms)
                    factors = []
                    for b in blocks:
                        free = tuple(sorted(set(b.get('p_atoms', ())) - fixed))
                        if len(free) > 1:
                            factors.append(self._intern(self._factors, ('pool', free)))
                    generators = tuple(sorted({tuple(g) for g in raw['automorph_generators']
                                               if any(i != p for i, p in enumerate(g))}))
                    if generators:
                        factors.append(self._intern(self._factors, ('group', generators)))
                    actions[eid] = tuple(factors)
                program.extend(actions[eid])
            required = self._intern(self._required, tuple(sorted(required)))
            program = self._intern(self._programs, canonical_actions(program))
            mapping = tuple(sorted(state.mapping))
            key = policy, mapping, required, program
            if key not in self._family_ids:
                self._family_ids[key] = len(self.families)
                self.families.append(FinalFamily(mapping, required, program, policy))
            family_id = self._family_ids[key]
            branch.add(family_id)
            self.families[family_id].provenance.append(
                dict(archive=str(archive), terminal=path.terminal, transitions=path.transitions))
        return self

    def rebuild_symmetry(self, policy, pair):
        """Intrinsic target automorphisms preserving every within-pair match test.

        Labels encode target adjacency and compatibility with all reactant bond
        weights of the same element pair. Thus these are exact automorphisms of
        a finite labelled graph even though matching uses a broad tolerance.
        No earlier assignments, seed, growth order or event threshold enters.
        """
        key = policy, pair
        if key in self.symmetries and key in self._labels:
            return self.symmetries[key]
        import pynauty
        floor, iso = policy
        r_atoms, p_atoms = pair
        r, p = self.problem.reactant, self.problem.product
        weights = defaultdict(set)
        for a, b in combinations(r_atoms, 2):
            if r.wbo[a, b] >= floor:
                weights[tuple(sorted((r.elements[a], r.elements[b])))].add(float(r.wbo[a, b]))
        weights = {k: tuple(sorted(v)) for k, v in weights.items()}
        adjacency = {i: set() for i in range(len(p_atoms))}
        colors = defaultdict(set)
        for i, atom in enumerate(p_atoms):
            colors[('atom', p.elements[atom])].add(i)
        import numpy as np
        labels = np.zeros((len(p_atoms), len(p_atoms)), dtype=np.int32)
        label_ids = {}
        next_vertex = len(p_atoms)
        for i, j in combinations(range(len(p_atoms)), 2):
            a, b = p_atoms[i], p_atoms[j]
            w = float(p.wbo[a, b])
            if w < floor:
                continue
            elements = tuple(sorted((p.elements[a], p.elements[b])))
            response = tuple(abs(v - w) <= iso + 1e-9 for v in weights.get(elements, ()))
            label = ('bond', elements, response)
            labels[i, j] = labels[j, i] = label_ids.setdefault(label, len(label_ids) + 1)
            colors[label].add(next_vertex)
            adjacency[next_vertex] = {i, j}
            adjacency[i].add(next_vertex)
            adjacency[j].add(next_vertex)
            next_vertex += 1
        graph = pynauty.Graph(next_vertex,
                             adjacency_dict={i: sorted(v) for i, v in adjacency.items()},
                             vertex_coloring=[colors[k] for k in sorted(colors, key=repr)])
        generators = tuple(sorted({tuple(g[:len(p_atoms)]) for g in pynauty.autgrp(graph)[0]
                                   if any(i != g[i] for i in range(len(p_atoms)))}))
        result = dict(r_atoms=r_atoms, p_atoms=p_atoms, generators=generators,
                      graph_floor=floor, iso_tolerance=iso,
                      policy='intrinsic_matching_response', local_target_indices=True)
        self.symmetries[key] = result
        self._labels[key] = labels
        return result

    def rebuild_all_symmetries(self):
        for policy, pairs in self.branches:
            for pair in pairs:
                self.rebuild_symmetry(policy, pair)
        return self.symmetries

    @property
    def branch_count(self):
        return len(self.branches) + len(self.coupled_branches)

    @classmethod
    def from_record(cls, problem, record):
        if record.get('schema') != cls.schema:
            raise ValueError('Unsupported final-branch schema')
        result = cls(problem)
        if record.get('problem_sha256') != result.problem_sha256:
            raise ValueError('Final branches belong to a different endpoint problem')
        result.mode = record.get('mode', 'lossless_saved_family_union')
        result.path_count = int(record['input_paths'])
        result.incomplete_paths = int(record.get('incomplete_paths', 0))
        result.reconstruction = dict(record.get('reconstruction', {}))
        result.families = [FinalFamily.from_record(f) for f in record['families']]
        for branch in record['branches']:
            pairs = tuple((tuple(rs), tuple(ps)) for rs, ps in branch['fragment_pairs'])
            result.branches[tuple(branch['policy']), pairs] = set(branch['families'])
        for branch in record.get('coupled_branches', ()):
            result.coupled_branches[branch['family']] = {
                tuple((tuple(rs), tuple(ps)) for rs, ps in pairs)
                for pairs in branch['source_fragment_pairs']}
        for ids in (*result.branches.values(), result.coupled_branches):
            if any(i < 0 or i >= len(result.families) for i in ids):
                raise ValueError('Invalid family reference')
        for raw in record.get('intrinsic_symmetries', ()):
            group = {**raw, 'r_atoms': tuple(raw['r_atoms']), 'p_atoms': tuple(raw['p_atoms']),
                     'generators': tuple(tuple(g) for g in raw['generators'])}
            key = ((group['graph_floor'], group['iso_tolerance']),
                   (group['r_atoms'], group['p_atoms']))
            result.symmetries[key] = group
        for i, f in enumerate(result.families):
            result._family_ids[f.policy, f.mapping, f.required_edges, f.actions] = i
        return result

    def to_record(self):
        return dict(schema=self.schema, problem_sha256=self.problem_sha256, mode=self.mode, reconstruction=self.reconstruction, input_paths=self.path_count,
                    incomplete_paths=self.incomplete_paths,
                    branches=[dict(policy=policy, fragment_pairs=pairs, families=sorted(families))
                              for (policy, pairs), families in self.branches.items()],
                    coupled_branches=[dict(family=k, source_fragment_pairs=list(v)) for k, v in self.coupled_branches.items()],
                    families=[f.to_record() for f in self.families],
                    intrinsic_symmetries=list(self.symmetries.values()))


def deduplicate_final_branches(aam):
    return FinalBranchCatalogue(aam.problem).add_aam(aam)

"""Mechanism-independent fragment-decision DAG and deterministic persistence.

Runtime frontier objects reference this graph; they never copy entire paths.
Joins are inserted only at synchronized frontier admission. Independent seed
and cut contexts retain separate nodes when graphs are combined.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, replace


def frozen_value(value):
    if isinstance(value, dict):
        return tuple(sorted((str(key), frozen_value(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(frozen_value(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return tuple(sorted((frozen_value(item) for item in value), key=repr))
    return value


@dataclass(frozen=True)
class SearchContext:
    source_atoms: tuple[int, ...]
    target_atoms: tuple[int, ...]
    seed_order: tuple[int, ...]
    cuts: tuple = ()
    core_atoms: tuple = ()
    anchors: tuple = ()
    graph_floor: float = 0.2
    iso_tolerance: float = 0.5
    branch_limit: int = 100
    objective: str = 'full_source'


@dataclass(frozen=True)
class SearchState:
    id: int
    context: int
    mapping: tuple[tuple[int, int], ...]
    islands: tuple[tuple[int, int], ...]
    deferred_edges: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class FragmentTransition:
    id: int
    source: int
    target: int
    seed: int | None
    step: tuple
    match: dict | None
    preserved_bonds: tuple = ()

    @property
    def placement(self):
        """Typed fragment view; a join has no fragment placement."""
        if self.match is None:
            return None
        from .alignment.post_aam import AAMHierarchy
        return AAMHierarchy.from_record({'fragments': (self.match,)}).fragments[0]


@dataclass(frozen=True)
class SearchStop:
    state: int
    reason: str
    seed: int | None = None
    stage: str = ""
    count: int = 0
    limit: int = 0
    step: tuple[int, ...] = ()


@dataclass(frozen=True)
class SearchPath:
    graph: AAMSearchGraph
    terminal: int
    transitions: tuple[int, ...]

    @property
    def mapping(self):
        return dict(self.graph.states[self.terminal].mapping)

    @property
    def deferred_edges(self):
        return self.graph.states[self.terminal].deferred_edges

    @property
    def fragments(self):
        return tuple(self.graph.transitions[index].match for index in self.transitions
                     if self.graph.transitions[index].match is not None)

    @property
    def symmetry_paths(self):
        return (self.fragments,)

    @property
    def hierarchy(self):
        from .alignment.post_aam import AAMHierarchy
        return AAMHierarchy.from_record({"fragments": self.fragments})

    @property
    def context(self):
        return self.graph.contexts[self.graph.states[self.terminal].context]

    def to_reference(self, graph_id):
        return {'graph': graph_id, 'terminal': self.terminal,
                'transitions': list(self.transitions)}

    @classmethod
    def from_reference(cls, record, graphs):
        return cls(graphs[record['graph']], record['terminal'], tuple(record['transitions']))

    def realize(self, generator_words):
        """Apply recorded correlated choices, transporting later decisions.

        generator_words maps transition IDs to sequences of generator indices.
        For chronological local actions h1,h2,... the full target action is
        h1 composed with h2 composed with ... . Later groups fix the original
        locked prefix; reversing that order would mix incompatible frames.
        """
        atoms = self.context.target_atoms
        action = {atom: atom for atom in atoms}
        if set(generator_words) - set(self.transitions):
            raise ValueError('generator choice is not on this path')
        for edge_id in self.transitions:
            for index in generator_words.get(edge_id, ()):
                placement = self.graph.transitions[edge_id].placement
                if placement is None or placement.target_generators is None:
                    raise ValueError('transition has no finalized exact group')
                generator = placement.target_generators[index].images
                action = {atom: action[generator[atom]] for atom in atoms}
        return PathRealization(self,
            tuple(sorted((atom, action[image]) for atom, image in self.mapping.items())),
            tuple(sorted(action.items())),
            tuple((edge, tuple(word)) for edge, word in generator_words.items()))

    def sample(self, rng, *, steps_per_fragment=1):
        """Nonuniform random walk in conditioned groups; no orbit expansion."""
        choices = {}
        for edge_id in self.transitions:
            placement = self.graph.transitions[edge_id].placement
            if placement is None:
                continue
            if placement.target_generators is None:
                raise ValueError('finalize exact groups before sampling a path')
            if placement.target_generators:
                choices[edge_id] = tuple(rng.randrange(len(placement.target_generators))
                                         for _ in range(steps_per_fragment))
        return self.realize(choices)


@dataclass(frozen=True)
class PathRealization:
    path: SearchPath
    mapping: tuple
    target_action: tuple
    generator_words: tuple

    @property
    def hierarchy(self):
        return self.path.hierarchy.relabel_target(self.target_action)


@dataclass(frozen=True)
class SearchBranch:
    """One literal matched-fragment relation and every path discovering it."""

    paths: tuple[SearchPath, ...]

    @property
    def representative(self):
        from .domain import AtomAssignment
        return AtomAssignment.from_mapping(self.paths[0].mapping)

    @property
    def hierarchy(self):
        return self.paths[0].hierarchy


@dataclass(frozen=True)
class AAMSearchGraph:
    contexts: tuple[SearchContext, ...]
    roots: tuple[int, ...]
    states: tuple[SearchState, ...]
    transitions: tuple[FragmentTransition, ...]
    stops: tuple[SearchStop, ...]

    @property
    def terminals(self):
        return tuple(dict.fromkeys(stop.state for stop in self.stops
                                   if stop.reason in {"objective_met", "stalled"}))

    @property
    def capped(self):
        return any(stop.reason == "capped" for stop in self.stops)

    def paths(self, terminal=None):
        """Lazily unfold recorded decision paths, never atom permutations."""
        incoming = [[] for _state in self.states]
        for edge in self.transitions:
            incoming[edge.target].append(edge)
        roots = set(self.roots)
        for end in self.terminals if terminal is None else (terminal,):
            stack = [(end, ())]
            while stack:
                state, suffix = stack.pop()
                if state in roots:
                    yield SearchPath(self, end, suffix)
                else:
                    for edge in reversed(incoming[state]):
                        stack.append((edge.source, (edge.id,) + suffix))

    def branches(self):
        """Deduplicate identical matched relations, retaining seed/cut paths."""
        relations = {}
        for path in self.paths():
            fragments = []
            for edge_id in path.transitions:
                edge = self.transitions[edge_id]
                if edge.match is None:
                    continue
                symmetry = {key: value for key, value in edge.match['symmetry'].items()
                            if key not in {'multiplicity', 'automorph_group_source'}}
                fragments.append((tuple(edge.match['fragment']),
                                  edge.preserved_bonds, frozen_value(symmetry),
                                  frozen_value(edge.match['deferred_edges'])))
            key = (tuple(sorted(path.mapping.items())), tuple(fragments))
            relations.setdefault(key, []).append(path)
        return tuple(SearchBranch(tuple(paths)) for paths in relations.values())

    def to_record(self):
        from dataclasses import asdict
        return {"schema": "rxn_core.aam_search_graph/v1",
                **asdict(self)}

    @classmethod
    def initial_fragment_search(cls, source, target, seed, result, config):
        """Record an unanchored one-seed call, without inventing continuation."""
        context = SearchContext(tuple(sorted(source)), tuple(sorted(target)),
            (int(seed),), graph_floor=config.graph_floor,
            iso_tolerance=config.iso_tolerance, branch_limit=config.branch_limit,
            objective='seeded_fragment')
        states = [SearchState(0, 0, (), (), ())]
        edges, stops = [], []
        for placement in result.matches:
            node = len(states)
            states.append(SearchState(node, 0, tuple(sorted(placement.items())),
                tuple((atom, 1) for atom in sorted(placement.fragment)),
                tuple(sorted(placement.deferred_edges))))
            record = {'island_idx': 1, 'fragment': sorted(placement.fragment),
                      'symmetry': copy.deepcopy(placement.symmetry),
                      'deferred_edges': sorted(placement.deferred_edges)}
            edges.append(FragmentTransition(len(edges), 0, node, seed, (1, 0),
                                            record, placement.preserved_bonds))
            stops.append(SearchStop(node, 'objective_met', seed, 'seeded_fragment'))
        if result.capped:
            stops.append(SearchStop(0, 'capped', seed, 'fragment_growth',
                                    result.branch_count, result.branch_limit))
        elif not result.matches:
            stops.append(SearchStop(0, 'stalled', seed, 'seeded_fragment'))
        return cls((context,), (0,), tuple(states), tuple(edges), tuple(stops))

    @classmethod
    def from_record(cls, record):
        if record["schema"] != "rxn_core.aam_search_graph/v1":
            raise ValueError("unsupported AAM search graph schema")
        def pairs(items):
            return tuple(tuple(item) for item in items)
        contexts = tuple(SearchContext(
            **{**item, "source_atoms": tuple(item["source_atoms"]),
               "target_atoms": tuple(item["target_atoms"]),
               "seed_order": tuple(item["seed_order"]), "cuts": pairs(item["cuts"]),
               "anchors": pairs(item["anchors"]), "core_atoms": tuple(item["core_atoms"])})
            for item in record["contexts"])
        states = tuple(SearchState(**{**item, "mapping": pairs(item["mapping"]),
            "islands": pairs(item["islands"]), "deferred_edges": pairs(item["deferred_edges"])})
            for item in record["states"])
        transitions = []
        for item in record["transitions"]:
            match = copy.deepcopy(item["match"])
            if match is not None:
                symmetry = match["symmetry"]
                symmetry["witness"] = {int(a): int(b) for a, b in symmetry["witness"].items()}
            transitions.append(FragmentTransition(**{**item, "match": match,
                "step": tuple(item["step"]), "preserved_bonds": pairs(item["preserved_bonds"])}))
        return cls(contexts, tuple(record["roots"]), states, tuple(transitions),
                   tuple(SearchStop(**{**item, 'step': tuple(item.get('step', ()))})
                         for item in record["stops"]))

    @classmethod
    def combine(cls, graphs):
        contexts, roots, states, transitions, stops = [], [], [], [], []
        for graph in graphs:
            c, s, t = len(contexts), len(states), len(transitions)
            contexts.extend(graph.contexts)
            roots.extend(root + s for root in graph.roots)
            states.extend(replace(state, id=state.id+s, context=state.context+c)
                          for state in graph.states)
            transitions.extend(replace(edge, id=edge.id+t, source=edge.source+s,
                                       target=edge.target+s) for edge in graph.transitions)
            stops.extend(replace(stop, state=stop.state+s) for stop in graph.stops)
        return cls(tuple(contexts), tuple(roots), tuple(states), tuple(transitions), tuple(stops))


class SearchGraphBuilder:
    def __init__(self, context):
        self.context = context
        self.states = []
        self.transitions = []
        self.stops = []
        self.roots = []
        self.step = (0, 0)
        self.seed = None

    def state(self, branch):
        node = SearchState(len(self.states), 0, tuple(sorted(branch.mapping.items())),
                           tuple(sorted(branch.islands_R.items())),
                           tuple(sorted(branch.deferred_edges)))
        self.states.append(node)
        return node.id

    def root(self, branch):
        node = self.state(branch)
        self.roots.append(node)
        return node

    def commit(self, parent, branch, match, graph):
        node = self.state(branch)
        atoms = set(match["fragment"])
        deferred = {tuple(sorted(edge)) for edge in match["deferred_edges"]}
        bonds = tuple(sorted(tuple(sorted((a, b))) for a, b in graph.edges()
                             if a in atoms and b in atoms
                             and tuple(sorted((a, b))) not in deferred))
        self.transitions.append(FragmentTransition(len(self.transitions), parent,
            node, self.seed, self.step, copy.deepcopy(match), bonds))
        return node

    def join(self, kept, other):
        if kept.node == other.node:
            return
        node = self.state(kept)
        for parent in (kept.node, other.node):
            self.transitions.append(FragmentTransition(len(self.transitions), parent,
                                                       node, self.seed, self.step, None))
        kept.node = node

    def stop(self, branch, reason, *, stage="", count=0, limit=0):
        self.stops.append(SearchStop(branch.node, reason, self.seed, stage, count, limit,
                                     self.step))

    def finish(self):
        return AAMSearchGraph((self.context,), tuple(self.roots), tuple(self.states),
                              tuple(self.transitions), tuple(self.stops))

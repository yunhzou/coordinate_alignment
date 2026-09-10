"""Shared conditional growth with fair exploration of source-cut conditions.

Cuts belong to the search domain, not a recovery fallback. Every condition is
scheduled without reference labels. This is a finite-budget adaptive policy,
not a guarantee of equivalence to ten independently seeded full cut sweeps.
"""
from collections import deque

from .aam import cut_seed
from .adaptive_search import AdaptiveResult, AdaptiveSearchCondition
from .adaptive_seed_search import AdaptiveSeedSearch
from .alignment.branch import _generate_seed_orders
from .alignment.sweep import cut_sweep_items
from .conditioned_symmetry import ConditionedSymmetryWorkspace
from .cut_replay import FragmentRepair
from .domain import AAMResult, AAMSearchConfig, AAMSearchMetrics
from .frag import build_graph
from .matcher import _nauty_orbits
from .search_graph import AAMSearchGraph


class AdaptiveCutSearch:
    """Advance one conditional fragment growth, sharing exact reusable work.

    Each cut has its own sequential mapping/DAG, because changing topology
    changes future constraints. A dependency-checked native cache shares full
    compressed growth results across cuts only when their inputs are equivalent.
    No witnesses are substituted for correlated placements.
    """
    def __init__(self, problem, config=None, *, policy='adaptive'):
        self.problem = problem
        self.config = config or AAMSearchConfig(seed_count=1)
        if policy not in ('adaptive', 'shared'):
            raise ValueError('Unknown cut search policy')
        if (policy == 'adaptive' and self.config.seed_count != 1) or self.config.anchors:
            raise ValueError('Adaptive cut search uses one guide per cut and no external anchors')
        self.policy = policy
        self.source, self.target = [build_graph(e.elements, e.wbo,
            bond_cut=self.config.graph_floor) for e in (problem.reactant, problem.product)]
        self.orbits = _nauty_orbits(self.target, wbo_tol=self.config.iso_tolerance)
        self.workspace = ConditionedSymmetryWorkspace(self.target, self.config.iso_tolerance)
        self.repair = FragmentRepair(self.source, self.target, self.orbits)
        self.cuts = tuple(tuple(c) for c in cut_sweep_items(problem.reactant.wbo, self.config.cut_floor))
        self.sessions = {}
        self.agenda = deque(range(len(self.cuts)))
        self.work = 0

    def _session(self, index):
        if index not in self.sessions:
            cut = self.cuts[index]
            active_cut = tuple(e for e in cut if self.source.has_edge(*e))
            view = self.repair.for_cut(active_cut)
            order = tuple(_generate_seed_orders(view.source, 1, rng_seed=cut_seed(cut),
                                                 seed_selection=self.config.seed_selection)[0])
            condition = AdaptiveSearchCondition(view.source, self.target, order, self.orbits,
                                                self.workspace, cut, view)
            if self.policy == 'shared':
                from .shared_seed_search import SharedSeedSearch
                self.sessions[index] = SharedSeedSearch(self.problem, self.config, condition=condition)
            else:
                self.sessions[index] = AdaptiveSeedSearch(self.problem, self.config, condition=condition)
        return self.sessions[index]

    def advance(self):
        if not self.agenda:
            return False
        index = self.agenda.popleft()
        session = self._session(index)
        session.advance()
        self.work += 1
        if session.agenda:
            self.agenda.append(index)
        return True

    def snapshot(self):
        results = [(i, self.sessions[i].snapshot()) for i in sorted(self.sessions)]
        graph = AAMSearchGraph.combine([r.aam.graph for _,r in results])
        pending = tuple(dict(cut_index=i, **p) for i,r in results for p in r.pending)
        pending += tuple(dict(kind='unvisited_cut', cut_index=i) for i in self.agenda if i not in self.sessions)
        metrics = AAMSearchMetrics.from_record(dict(cuts=len(self.cuts),
            raw_result_count=len(graph.terminals), retained_branch_count=len(graph.terminals),
            subtree_branch_cap_count=sum(s.reason == 'capped' for s in graph.stops)),
            sum(s.elapsed for s in self.sessions.values()))
        return AdaptiveResult(AAMResult(self.problem, self.config, graph, metrics), self.work,
            sum(r.growth_calls for _,r in results), 0, sum(r.reused_states for _,r in results),
            pending, not self.agenda)

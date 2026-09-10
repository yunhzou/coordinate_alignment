"""Experimental persistent sequential fragment agenda; no external sweep cuts.

Queue priority is a heuristic, not an admissible event bound. No family is
discarded by a sampled mapping's score. Work exhaustion is explicitly reported.
"""
from dataclasses import dataclass, replace
import heapq
import itertools
import time

from .alignment.branch import _Branch, _generate_seed_orders
from .domain import AAMResult, AAMSearchConfig, AAMSearchMetrics
from .frag import build_graph
from .fragment import FragmentMatchConfig, FragmentMatchContext
from .fragment_choices import FragmentChoices
from .matcher import _nauty_orbits
from .search_graph import SearchContext, SearchGraphBuilder, FragmentTransition, SearchStop
from .search_symmetry import finalize_graph_symmetry
from .conditioned_symmetry import ConditionedSymmetryWorkspace


@dataclass(frozen=True)
class AdaptiveResult:
    aam: AAMResult
    work: int
    growth_calls: int
    closure_calls: int
    reused_states: int
    pending: tuple
    exhausted: bool


class AdaptiveFragmentSearch:
    """One resumable search session with compressed, condition-dependent choices.

    advance() consumes one growth/closure operation. Callers own time/work
    budgets and can save snapshots without restarting completed work.
    """
    def __init__(self, problem, config=None, *, policy='largest_first'):
        self.problem = problem
        self.config = config or AAMSearchConfig(seed_count=1)
        if self.config.seed_count != 1:
            raise ValueError('Adaptive search uses one persistent seed order')
        if policy not in ('largest_first', 'smallest_first'):
            raise ValueError('unknown closure priority')
        self.policy = policy
        self.source, self.target = [build_graph(endpoint.elements, endpoint.wbo,
            bond_cut=self.config.graph_floor) for endpoint in (problem.reactant, problem.product)]
        self.orbits = _nauty_orbits(self.target, wbo_tol=self.config.iso_tolerance)
        self.order = tuple(_generate_seed_orders(self.source, 1, rng_seed=42,
                                                 seed_selection=self.config.seed_selection)[0])
        self.builder = SearchGraphBuilder(SearchContext(tuple(sorted(self.source)), tuple(sorted(self.target)),
            self.order, anchors=self.config.anchors, graph_floor=self.config.graph_floor,
            iso_tolerance=self.config.iso_tolerance, branch_limit=self.config.branch_limit))
        self.workspace = ConditionedSymmetryWorkspace(self.target, self.config.iso_tolerance)
        self.agenda = []
        self.serial = itertools.count()
        self.seen = {}
        self.work = self.growth_calls = self.closure_calls = self.reused_states = 0
        self.elapsed = 0.
        root = _Branch(self.builder, self.config.anchors)
        self._schedule(root, 0, 0)

    def _push(self, priority, payload):
        heapq.heappush(self.agenda, (priority, next(self.serial), payload))

    def _schedule(self, branch, position, depth):
        while position < len(self.order) and self.order[position] in branch.mapping:
            position += 1
        key = (branch.state_key(), position)
        previous = self.seen.get(key)
        if previous is not None:
            self.reused_states += 1
            if previous.node != branch.node:
                self.builder.transitions.append(FragmentTransition(len(self.builder.transitions),
                    branch.node, previous.node, self.builder.seed, self.builder.step, None))
            return
        self.seen[key] = branch
        if position == len(self.order):
            self.builder.stop(branch, 'objective_met' if len(branch.mapping) == len(self.source) else 'stalled')
        else:
            self._push((depth, 0, -len(branch.mapping)), ('grow', branch, position, depth, None, None))

    def _commit(self, branch, result, position, depth):
        if result.capped:
            self.builder.stop(branch, 'capped', stage='fragment_growth',
                              count=result.branch_count, limit=result.branch_limit)
            return
        if not result.matches:
            self._schedule(branch, position + 1, depth)
        for match in result.matches:
            child = branch.fork()
            child.commit(match)
            if child.state_key() != branch.state_key():
                self._schedule(child, 0, depth)

    def advance(self):
        if not self.agenda:
            return False
        started = time.perf_counter()
        _, _, (kind, branch, position, depth, session, closure) = heapq.heappop(self.agenda)
        seed = self.order[position]
        self.builder.seed = seed
        self.builder.step = (depth, position, self.work)
        self.work += 1
        if kind == 'grow':
            self.growth_calls += 1
            session = FragmentChoices(self.source, self.target, seed=seed,
                context=FragmentMatchContext(branch.mapping, branch.islands_R,
                    tuple(branch.deferred_edges), target_orbits=self.orbits),
                config=FragmentMatchConfig(graph_floor=self.config.graph_floor,
                    iso_tolerance=self.config.iso_tolerance, branch_limit=self.config.branch_limit))
            for option in session.closures:
                rank = -option.size if self.policy == 'largest_first' else option.size
                self._push((depth + 1, 1, rank), ('close', branch, position, depth + 1, session, option))
            self._commit(branch, session.normal, position, depth)
        else:
            self.closure_calls += 1
            self._commit(branch, session.close(closure), position, depth)
        self.elapsed += time.perf_counter() - started
        return True

    def snapshot(self):
        pending = tuple(dict(kind=task[0], state=task[1].node, position=task[2], depth=task[3],
                             closure=None if task[5] is None else task[5].index)
                        for _, _, task in sorted(self.agenda))
        graph = self.builder.finish()
        if pending:
            graph = replace(graph, stops=graph.stops + tuple(SearchStop(state, 'budget', stage='adaptive_agenda')
                for state in sorted({p['state'] for p in pending})))
        graph, _ = finalize_graph_symmetry(graph, self.target, iso_tolerance=self.config.iso_tolerance,
                                          workspace=self.workspace)
        metrics = AAMSearchMetrics.from_record(dict(cuts=1, raw_result_count=len(graph.terminals),
            retained_branch_count=len(graph.terminals),
            subtree_branch_cap_count=sum(s.reason == 'capped' for s in graph.stops)), self.elapsed)
        return AdaptiveResult(AAMResult(self.problem, self.config, graph, metrics), self.work,
            self.growth_calls, self.closure_calls, self.reused_states, pending, not pending)

"""Experimental persistent sequential fragment agenda; no external sweep cuts.

Queue priority is a heuristic, not an admissible event bound. No family is
discarded by a sampled mapping's score. Work exhaustion is explicitly reported.
"""
from dataclasses import dataclass, replace
import heapq
import itertools
import time
import numpy as np
from collections import Counter

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


@dataclass(frozen=True)
class AdaptiveSearchCondition:
    """Prepared directed graphs and one explicit topology/seed condition."""
    source: object
    target: object
    order: tuple
    orbits: object
    workspace: object
    cuts: tuple = ()
    growth_replay: object = None

    @classmethod
    def uncut(cls, problem, config):
        source, target = [build_graph(endpoint.elements, endpoint.wbo,
            bond_cut=config.graph_floor) for endpoint in (problem.reactant, problem.product)]
        order = tuple(_generate_seed_orders(source, 1, rng_seed=42,
                                             seed_selection=config.seed_selection)[0])
        return cls(source, target, order, _nauty_orbits(target, wbo_tol=config.iso_tolerance),
                   ConditionedSymmetryWorkspace(target, config.iso_tolerance))


class _ChoiceAgenda:
    """Stable queues per discrepancy depth; optionally share work fairly."""
    def __init__(self, fair):
        self.levels = {}
        self.served = Counter()
        self.serial = itertools.count()
        self.count = 0
        self.fair = fair

    def __len__(self):
        return self.count

    def __iter__(self):
        return itertools.chain.from_iterable(self.levels.values())

    def push(self, priority, payload):
        heapq.heappush(self.levels.setdefault(priority[0], []), (priority, next(self.serial), payload))
        self.count += 1

    def pop(self):
        depth = min(self.levels, key=lambda d: (self.served[d], d)) if self.fair else min(self.levels)
        item = heapq.heappop(self.levels[depth])
        self.served[depth] += 1
        self.count -= 1
        if not self.levels[depth]:
            del self.levels[depth]
        return item

    def reprioritize(self, function):
        for depth, heap in self.levels.items():
            self.levels[depth] = [(function(priority, task), serial, task) for priority, serial, task in heap]
            heapq.heapify(self.levels[depth])


class AdaptiveFragmentSearch:
    """One resumable search session with compressed, condition-dependent choices.

    advance() consumes one growth/closure operation. Callers own time/work
    budgets and can save snapshots without restarting completed work.

    This experiment retains the native per-growth branch cap, not the mature
    scheduler's synchronized live-frontier cap. The agenda is budgeted by the
    caller. Snapshots persist results and pending descriptors, not native
    sessions for cross-process resume. No cut sweep or symmetry-repair pass is
    run by this scheduler.
    """
    def __init__(self, problem, config=None, *, policy='largest_first', condition=None):
        self.problem = problem
        self.config = config or AAMSearchConfig(seed_count=1)
        if self.config.seed_count != 1:
            raise ValueError('Adaptive search uses one persistent seed order')
        if self.config.anchors:
            raise ValueError('Adaptive search does not yet support anchored growth')
        if policy not in ('largest_first', 'smallest_first', 'event_guided', 'fair_depth'):
            raise ValueError('unknown closure priority')
        self.policy = policy
        self.condition = condition or AdaptiveSearchCondition.uncut(problem, self.config)
        self.source, self.target = self.condition.source, self.condition.target
        self.orbits, self.order = self.condition.orbits, self.condition.order
        self.builder = SearchGraphBuilder(SearchContext(tuple(sorted(self.source)), tuple(sorted(self.target)),
            self.order, cuts=self.condition.cuts, anchors=self.config.anchors, graph_floor=self.config.graph_floor,
            iso_tolerance=self.config.iso_tolerance, branch_limit=self.config.branch_limit))
        self.workspace = self.condition.workspace
        self.agenda = _ChoiceAgenda(fair=policy == 'fair_depth')
        self.seen = {}
        self.work = self.growth_calls = self.closure_calls = self.reused_states = 0
        self.elapsed = 0.
        self.best_score = None
        self.error_edges = ()
        root = _Branch(self.builder, self.config.anchors)
        self._schedule(root, 0, 0)

    def _push(self, priority, payload):
        self.agenda.push(priority, payload)

    def _closure_priority(self, depth, session, option):
        rank = option.size if self.policy == 'smallest_first' else -option.size
        affected = 0.
        if self.policy in ('event_guided', 'fair_depth'):
            released = session.frontier_atoms - option.atoms
            if released:
                affected = sum(a in released or b in released for a, b in self.error_edges) / len(released)
        return (depth, 1, -affected, rank)

    def _feedback(self, mapping):
        if self.policy not in ('event_guided', 'fair_depth') or not self.problem.balanced:
            return
        r = self.problem.reactant.wbo
        images = [mapping[i] for i in range(len(r))]
        p = self.problem.product.wbo[np.ix_(images, images)]
        rb, pb = r > self.config.graph_floor, p > self.config.graph_floor
        errors = np.triu((rb != pb) | (rb & pb & (np.abs(r-p) > self.config.event_threshold)), 1)
        score = int(errors.sum())
        if self.best_score is None or score < self.best_score:
            self.best_score = score
            self.error_edges = tuple(zip(*np.where(errors)))
            self.agenda.reprioritize(lambda priority, task:
                self._closure_priority(task[3], task[4], task[5]) if task[0] == 'close' else priority)

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
            if len(branch.mapping) == len(self.source):
                self._feedback(branch.mapping)
        else:
            self._push((depth, 0, 0, -len(branch.mapping)), ('grow', branch, position, depth, None, None))

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
        _, _, (kind, branch, position, depth, session, closure) = self.agenda.pop()
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
                self._push(self._closure_priority(depth + 1, session, option),
                           ('close', branch, position, depth + 1, session, option))
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

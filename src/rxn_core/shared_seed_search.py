"""Exact policy sharing for sequential AAM, without rerunning common decisions.

Policies retain the original shuffled order, repeated passes and synchronized
frontier cap. Only equal conditional states/decisions and fragment records are
shared. This is not a smaller seed sample or a reference-guided recovery pass.
"""
from collections import deque
from dataclasses import replace
import time

from .aam import cut_seed
from .adaptive_search import AdaptiveResult
from .alignment.branch import _Branch, _generate_seed_orders
from .domain import AAMResult, AAMSearchMetrics
from .fragment import match_fragment, FragmentMatchConfig, FragmentMatchContext
from .matcher import _nauty_orbits
from .search_graph import SearchContext, SearchGraphBuilder, FragmentTransition, SearchStop, frozen_value
from .search_symmetry import finalize_graph_symmetry


class _SharedGraphBuilder(SearchGraphBuilder):
    """Hash-cons exact states and equal fragment edges, preserving correlations."""
    def __init__(self, context):
        super().__init__(context)
        self.nodes = {}
        self.edges = set()
        self.stop_keys = set()

    def state(self, branch):
        key = branch.state_key()
        if key not in self.nodes:
            self.nodes[key] = super().state(branch)
        return self.nodes[key]

    def commit(self, parent, branch, match, preserved_bonds):
        node = self.state(branch)
        key = (parent, node, frozen_value(match), tuple(preserved_bonds))
        if key not in self.edges:
            self.edges.add(key)
            self.transitions.append(FragmentTransition(len(self.transitions), parent,
                node, self.seed, self.step, match, preserved_bonds))
        return node

    def stop(self, branch, reason, *, stage='', count=0, limit=0):
        key = (branch.node, reason, stage, count, limit)
        if key not in self.stop_keys:
            self.stop_keys.add(key)
            super().stop(branch, reason, stage=stage, count=count, limit=limit)


class SharedSeedSearch:
    """One cut's original seed policies sharing an exact decision DAG.

The future of a seed decision depends on mapping, island partition, deferred
edges and the cut condition, not the history or which policy reached it.
Policy cursors/frontier admission remain separate: their cap decisions are
unchanged. Identical growth output is never expanded into bijections.
"""
    def __init__(self, problem, config, *, condition):
        if config.anchors:
            raise ValueError('Shared policy search currently accepts unanchored AAM only')
        self.problem, self.config, self.condition = problem, config, condition
        self.source, self.target = condition.source, condition.target
        self.source_orbits = _nauty_orbits(self.source, wbo_tol=config.iso_tolerance)
        self.orders = _generate_seed_orders(self.source, config.seed_count,
            rng_seed=cut_seed(condition.cuts), seed_selection=config.seed_selection)
        # Multiple policies use this context; edge.step records policy/pass/position.
        self.builder = _SharedGraphBuilder(SearchContext(tuple(sorted(self.source)),
            tuple(sorted(self.target)), (), cuts=condition.cuts, graph_floor=config.graph_floor,
            iso_tolerance=config.iso_tolerance, branch_limit=config.branch_limit))
        self.root = _Branch(self.builder)
        self.branches = {self.root.state_key():self.root}
        self.decisions = {}
        self.work = self.growth_calls = self.reused_states = 0
        self.elapsed = 0.
        self.frontiers = {}
        self.agenda = deque()
        for index, order in enumerate(self.orders):
            policy = self._policy(index, order)
            request = next(policy, None)
            if request is not None:
                self.agenda.append((policy, request))

    def _expand(self, branch, seed):
        key = (branch.node, seed)
        if key in self.decisions:
            self.reused_states += 1
            return self.decisions[key]
        self.growth_calls += 1
        result = match_fragment(self.source, self.target, seed=seed,
            context=FragmentMatchContext(branch.mapping, branch.islands_R,
                tuple(branch.deferred_edges), self.source_orbits, self.condition.orbits,
                growth_replay=self.condition.growth_replay),
            config=FragmentMatchConfig(graph_floor=self.config.graph_floor,
                iso_tolerance=self.config.iso_tolerance, branch_limit=self.config.branch_limit))
        if result.capped:
            self.builder.stop(branch, 'capped', stage='fragment_growth',
                              count=result.branch_count, limit=result.branch_limit)
            children, changed = (), False
        elif not result.matches:
            children, changed = (branch,), False
        else:
            children, changed = [], False
            for match in result.matches:
                child = branch.fork()
                child.commit(match)
                changed |= child.state_key() != branch.state_key()
                children.append(self.branches.setdefault(child.state_key(), child))
            children = tuple(children)
        self.decisions[key] = (children, changed)
        return children, changed

    def _policy(self, index, order):
        """Original synchronized admission, with growth supplied by the shared DAG."""
        branches, progressed, pass_no = [self.root], True, 0
        cap = self.config.branch_limit
        while progressed:
            progressed = False
            pass_no += 1
            for position, seed in enumerate(order):
                if not any(seed not in branch.mapping for branch in branches) and len(branches) <= cap:
                    continue
                admitted, seen = [], set()
                for branch in branches:
                    self.frontiers[index] = tuple(branches) + tuple(admitted)
                    if seed in branch.mapping:
                        subtree, changed = (branch,), False
                    else:
                        subtree, changed = yield (index, pass_no, position, seed, branch)
                    additions, local = [], set()
                    for child in subtree:
                        if child.node not in seen and child.node not in local:
                            additions.append(child)
                            local.add(child.node)
                    if len(admitted)+len(additions) > cap:
                        for child in additions:
                            self.builder.stop(child, 'capped', stage='combined_live_leaves',
                                count=len(admitted)+len(additions), limit=cap)
                        continue
                    admitted.extend(additions)
                    seen.update(local)
                    progressed |= bool(additions) and changed
                branches = sorted(admitted, key=lambda branch:-len(branch.mapping))
        for branch in branches:
            self.builder.stop(branch, 'objective_met' if len(branch.mapping)==len(self.source) else 'stalled')
        self.frontiers.pop(index, None)

    def advance(self):
        if not self.agenda:
            return False
        started = time.perf_counter()
        policy, (index, pass_no, position, seed, branch) = self.agenda.popleft()
        self.builder.seed, self.builder.step = seed, (index, pass_no, position)
        self.work += 1
        output = self._expand(branch, seed)
        try:
            request = policy.send(output)
        except StopIteration:
            pass
        else:
            self.agenda.append((policy, request))
        self.elapsed += time.perf_counter()-started
        return True

    def snapshot(self):
        pending = tuple(dict(kind='policy_frontier', policy=index, state=state)
            for index, frontier in sorted(self.frontiers.items())
            for state in sorted({branch.node for branch in frontier}))
        graph = self.builder.finish()
        if pending:
            graph = replace(graph, stops=graph.stops + tuple(SearchStop(state, 'budget', stage='shared_policy')
                for state in sorted({p['state'] for p in pending})))
        graph, _ = finalize_graph_symmetry(graph, self.target, iso_tolerance=self.config.iso_tolerance,
                                          workspace=self.condition.workspace)
        metrics = AAMSearchMetrics.from_record(dict(cuts=1, raw_result_count=len(graph.terminals),
            retained_branch_count=len(graph.terminals),
            subtree_branch_cap_count=sum(s.reason=='capped' for s in graph.stops)), self.elapsed)
        return AdaptiveResult(AAMResult(self.problem, self.config, graph, metrics), self.work,
            self.growth_calls, 0, self.reused_states, pending, not self.agenda)

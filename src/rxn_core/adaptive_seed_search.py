"""Experimental shared-state alternatives for the next sequential fragment seed."""
import time
import heapq
import itertools
from collections import Counter

from .adaptive_search import AdaptiveFragmentSearch
from .fragment import match_fragment, FragmentMatchConfig, FragmentMatchContext
from .matcher import _nauty_orbits
from .search_graph import FragmentTransition


class _SeedAgenda:
    """Finish one route, then yield to a different conditional seed choice.

    Normal siblings remain pending. Finishing all of them before revisiting a
    seed would expand a whole search tree before trying another fragment order.
    """
    def __init__(self):
        self.ordinary=[]
        self.levels={}
        self.served=Counter()
        self.serial=itertools.count()
        self.yield_to_alternatives=False

    def __len__(self):
        return len(self.ordinary)+sum(map(len,self.levels.values()))

    def __iter__(self):
        return itertools.chain(self.ordinary,itertools.chain.from_iterable(self.levels.values()))

    def push(self,priority,payload):
        item=(priority,next(self.serial),payload)
        if priority[1]==0:
            self.ordinary.append(item)
        else:
            heapq.heappush(self.levels.setdefault(priority[0],[]),item)

    def pop(self):
        if self.levels and (self.yield_to_alternatives or not self.ordinary):
            depth=min(self.levels,key=lambda d:(self.levels[d][0][0][2],self.served[d],d))
            self.served[depth] += 1
            self.yield_to_alternatives=False
            item=heapq.heappop(self.levels[depth])
            if not self.levels[depth]:del self.levels[depth]
            return item
        return self.ordinary.pop()


class AdaptiveSeedSearch(AdaptiveFragmentSearch):
    """Explore seed choices without restarting completed conditional prefixes.

    Source-region coverage orders untried seeds; it never removes a seed or
    compressed placement. There is one lazy seed task per discovered state.
    Native per-growth cap and caller-owned work/time budgets still apply.
    """
    def __init__(self, problem, config=None):
        self.seed_trials = {}
        super().__init__(problem, config)
        pending=tuple(self.agenda)
        self.agenda=_SeedAgenda()
        for priority,_,payload in pending:self.agenda.push(priority,payload)
        self.source_orbits = _nauty_orbits(self.source, wbo_tol=self.config.iso_tolerance)

    def _next_seed(self, branch, trial):
        candidates=[(atom in trial['covered'], i) for i,atom in enumerate(self.order)
                    if atom not in branch.mapping and atom not in trial['tried']]
        return min(candidates)[1] if candidates else None

    def _schedule(self, branch, position, depth):
        if len(branch.mapping)==len(self.source):
            self.agenda.yield_to_alternatives=True
        key=branch.state_key()
        previous=self.seen.get(key)
        if previous is not None:
            self.reused_states += 1
            if previous.node != branch.node:
                self.builder.transitions.append(FragmentTransition(len(self.builder.transitions),
                    branch.node, previous.node, self.builder.seed, self.builder.step, None))
            return
        self.seen[key]=branch
        if len(branch.mapping)==len(self.source):
            self.builder.stop(branch,'objective_met')
            return
        trial=dict(tried=set(),covered=set(),productive=False,depth=depth)
        self.seed_trials[key]=trial
        position=self._next_seed(branch,trial)
        self._push((depth,0,0,-len(branch.mapping)),('grow',branch,position,depth,None,None))

    def advance(self):
        if not self.agenda:
            return False
        started=time.perf_counter()
        _,_,(_,branch,position,depth,_,_)=self.agenda.pop()
        seed=self.order[position]
        trial=self.seed_trials[branch.state_key()]
        self.builder.seed=seed
        self.builder.step=(depth,position,self.work)
        self.work += 1
        self.growth_calls += 1
        trial['tried'].add(seed)
        result=match_fragment(self.source,self.target,seed=seed,
            context=FragmentMatchContext(branch.mapping,branch.islands_R,
                tuple(branch.deferred_edges),self.source_orbits,self.orbits),
            config=FragmentMatchConfig(graph_floor=self.config.graph_floor,
                iso_tolerance=self.config.iso_tolerance,branch_limit=self.config.branch_limit))
        for match in result.matches:
            trial['covered'].update(match.fragment)
        trial['productive'] |= bool(result.matches)
        next_seed=self._next_seed(branch,trial)
        if next_seed is not None:
            alternative_depth=trial['depth']+1
            covered=self.order[next_seed] in trial['covered']
            self._push((alternative_depth,1,int(covered),len(branch.mapping)),
                ('grow',branch,next_seed,alternative_depth,None,None))
        elif not trial['productive']:
            self.builder.stop(branch,'stalled')
        if result.capped or result.matches:
            self._commit(branch,result,position,depth)
        self.elapsed += time.perf_counter()-started
        return True

"""Frozen unoptimized event encoder, used only as an independent test oracle."""
from collections import defaultdict
import time

def _event_model(compiled, index, deadline):
    """Small local pair tables; never enumerate complete assignments."""
    import z3
    terms, objective, variable_pairs, table_entries, lower = {}, [], 0, 0, 0
    for a, b, tau in zip(index.a, index.b, index.tau):
        if time.perf_counter() >= deadline:
            raise TimeoutError('event encoding budget')
        va, vb = compiled.values[a], compiled.values[b]
        groups = defaultdict(lambda: defaultdict(list))
        for x in sorted(va[1]):
            if time.perf_counter() >= deadline:
                raise TimeoutError('event encoding budget')
            for y in sorted(vb[1]):
                if x == y:
                    continue
                delta = index.p[x, y] - index.r[a, b]
                sign = -1 if delta <= -tau else 1 if delta >= tau else 0
                groups[sign][x].append(y)
                table_entries += 1
        if not groups:
            compiled.solver.add(False)
            term = 0
        elif len(groups) == 1:
            term = next(iter(groups))
        else:
            default = max(groups, key=lambda k: sum(map(len, groups[k].values())))
            term = default
            for sign, rows in groups.items():
                if sign != default:
                    condition = z3.Or(*(z3.And(va[0] == x, z3.Or(*(vb[0] == y for y in ys)))
                                        for x, ys in rows.items()))
                    term = z3.If(condition, sign, term)
            variable_pairs += 1
        lower += min((int(sign != 0) for sign in groups), default=0)
        terms[int(a), int(b)] = term
        objective.append(int(term != 0) if isinstance(term, int) else z3.If(term != 0, 1, 0))
    return terms, z3.Sum(objective) if objective else z3.IntVal(0), dict(
        variable_pairs=variable_pairs, pair_table_entries=table_entries, event_lower_bound=lower)


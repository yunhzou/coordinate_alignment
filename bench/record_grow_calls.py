"""Record every grow_island call (inputs + outputs) of a serial search_aam run.

    .venv/bin/python bench/record_grow_calls.py <case> <out.pkl>

Inputs are stored in plain-Python form (graphs as element lists, WBO matrices,
bond cut; orbit maps as dicts with their bucket tables) so a native engine can
be replayed on each call and compared with bench/compare_grow_calls.py.
"""
import hashlib, pickle, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "bench"))
from cases import CASES
import rxn_core.growth.island as island_mod
import rxn_core.alignment.branch as branch_mod
from rxn_core.aam import search_aam
from rxn_core.domain import AAMProblem, AAMSearchConfig
from rxn_core.growth.result import IslandBranchLimitExceeded


def graph_record(g, cache):
    # keyed by content: object ids are recycled once a cut graph is freed
    nodes = sorted(g.nodes())
    rec = {
        "elements": [g.nodes[v].get("element") for v in nodes],
        "wbo": g.graph["wbo_matrix"].tolist(),
        "bond_cut": float(g.graph.get("bond_cut", 0.2)),
        "edges": sorted((min(a, b), max(a, b)) for a, b in g.edges()),
    }
    key = hashlib.blake2b(repr(rec).encode(), digest_size=16).hexdigest()
    cache.setdefault(key, rec)
    return key


def orbit_record(o, cache):
    if o is None:
        return None
    rec = {
        "orbits": dict(o),
        "wbo_buckets": {k: v for k, v in getattr(o, "wbo_buckets", {}).items()},
        "zero_bucket": getattr(o, "zero_bucket", None),
        "wbo_tol": getattr(o, "wbo_tol", None),
    }
    key = hashlib.blake2b(repr(sorted(rec["orbits"].items())).encode() + repr(sorted(rec["wbo_buckets"].items())).encode()
                          + repr((rec["zero_bucket"], rec["wbo_tol"])).encode(), digest_size=16).hexdigest()
    cache.setdefault(key, rec)
    return key


def iso_record(iso):
    return {
        "mapping": dict(iso),
        "deferred_edges": sorted(tuple(e) for e in iso.deferred_edges),
        "fragment": sorted(iso.fragment),
        "symmetry": iso.symmetry,
    }


def main():
    case, out = sys.argv[1], sys.argv[2]
    R, P = CASES[case]()
    graphs, orbits, calls = {}, {}, []
    orig = island_mod.grow_island

    def recording(g_R, g_P, seed, mapping, graph_floor=0.2, iso_tol=0.5, min_lock_size=1,
                  max_branches=1_000_000, events=None, islands_R=None, p_orbits=None,
                  r_orbits=None, prior_deferred_edges=None, node_policy=None,
                  allow_mapped_seed=False, profile=None, profile_context=None):
        record = {
            "g_R": graph_record(g_R, graphs), "g_P": graph_record(g_P, graphs),
            "seed": int(seed), "mapping": dict(mapping), "graph_floor": graph_floor,
            "iso_tol": iso_tol, "min_lock_size": min_lock_size, "max_branches": max_branches,
            "islands_R": dict(islands_R) if islands_R is not None else None,
            "p_orbits": orbit_record(p_orbits, orbits), "r_orbits": orbit_record(r_orbits, orbits),
            "prior_deferred_edges": sorted(tuple(e) for e in (prior_deferred_edges or ())),
            "allow_mapped_seed": bool(allow_mapped_seed),
        }
        try:
            result = orig(g_R, g_P, seed, mapping, graph_floor=graph_floor, iso_tol=iso_tol,
                          min_lock_size=min_lock_size, max_branches=max_branches, events=events,
                          islands_R=islands_R, p_orbits=p_orbits, r_orbits=r_orbits,
                          prior_deferred_edges=prior_deferred_edges, node_policy=node_policy,
                          allow_mapped_seed=allow_mapped_seed, profile=profile,
                          profile_context=profile_context)
        except IslandBranchLimitExceeded as exc:
            record["raised"] = ("IslandBranchLimitExceeded", int(exc.count), int(exc.limit))
            calls.append(record)
            raise
        record["result"] = [iso_record(iso) for iso in result]
        calls.append(record)
        return result

    island_mod.grow_island = recording
    branch_mod.grow_island = recording
    t0 = time.perf_counter()
    search_aam(AAMProblem(R, P, name=case), AAMSearchConfig(), workers=1)
    elapsed = time.perf_counter() - t0
    with open(out, "wb") as handle:
        pickle.dump({"case": case, "graphs": graphs, "orbits": orbits, "calls": calls}, handle,
                    protocol=pickle.HIGHEST_PROTOCOL)
    print(f"{case}: recorded {len(calls)} grow_island calls in {elapsed:.1f}s -> {out}")


if __name__ == "__main__":
    main()

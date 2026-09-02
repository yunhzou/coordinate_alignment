"""Replay recorded grow_island calls through the native engine and compare.

    .venv/bin/python bench/compare_grow_calls.py <corpus.pkl> [--stop-at-first]

Every recorded call (inputs and the Python engine's outputs, see
bench/record_grow_calls.py) is rebuilt as networkx graphs and orbit maps and
run through rxn_core.growth.native.grow_island.  Mappings, deferred edges,
fragments, symmetry states (order-insensitive JSON) and cap behaviour must
match exactly.  On a mismatch the first differing call is dumped.
"""
import json, pickle, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
from rxn_core.frag import build_graph
from rxn_core.matcher.orbits import _OrbitMap
from rxn_core.growth import native
from rxn_core.growth.result import IslandBranchLimitExceeded


def rebuild_graph(rec):
    g = build_graph(rec["elements"], np.asarray(rec["wbo"], dtype=float), bond_cut=rec["bond_cut"])
    # the recorded graph may be a cut graph: match its edge set exactly
    wanted = set(tuple(e) for e in rec["edges"])
    for a, b in list(g.edges()):
        if (min(a, b), max(a, b)) not in wanted:
            g.remove_edge(a, b)
    return g


def rebuild_orbits(rec):
    if rec is None:
        return None
    return _OrbitMap(rec["orbits"], wbo_buckets=rec["wbo_buckets"],
                     zero_bucket=rec["zero_bucket"], wbo_tol=rec["wbo_tol"])


def iso_key(iso):
    return (
        tuple(sorted((int(k), int(v)) for k, v in iso["mapping"].items())),
        tuple(sorted(tuple(e) for e in iso["deferred_edges"])),
        tuple(sorted(iso["fragment"])),
        json.dumps(iso["symmetry"], sort_keys=True, default=list),
    )


def result_key(isos):
    return [iso_key({"mapping": dict(iso), "deferred_edges": iso.deferred_edges,
                     "fragment": iso.fragment, "symmetry": iso.symmetry}) for iso in isos]


def main():
    corpus = pickle.load(open(sys.argv[1], "rb"))
    stop = "--stop-at-first" in sys.argv
    graphs = {k: rebuild_graph(v) for k, v in corpus["graphs"].items()}
    orbits = {k: rebuild_orbits(v) for k, v in corpus["orbits"].items()}
    assert native.available(), "native engine not importable"
    mismatches = 0
    t_native = 0.0
    for index, call in enumerate(corpus["calls"]):
        g_R, g_P = graphs[call["g_R"]], graphs[call["g_P"]]
        p_orbits = orbits[call["p_orbits"]]
        t0 = time.perf_counter()
        try:
            out = native.grow_island(
                g_R, g_P, call["seed"], call["mapping"], graph_floor=call["graph_floor"],
                iso_tol=call["iso_tol"], min_lock_size=call["min_lock_size"],
                max_branches=call["max_branches"], islands_R=call["islands_R"],
                p_orbits=p_orbits, prior_deferred_edges=call["prior_deferred_edges"],
                allow_mapped_seed=call["allow_mapped_seed"], profile=None,
                profile_context=None)
            got = ("ok", result_key(out) if out is not None else None)
        except IslandBranchLimitExceeded as exc:
            got = ("raised", (exc.count, exc.limit))
        t_native += time.perf_counter() - t0
        if "raised" in call:
            expected = ("raised", (call["raised"][1], call["raised"][2]))
        else:
            expected = ("ok", [iso_key(iso) for iso in call["result"]])
        if got != expected:
            mismatches += 1
            print(f"MISMATCH call {index}: seed {call['seed']} mapped {len(call['mapping'])} "
                  f"islands {len(call['islands_R'] or {})} deferred {len(call['prior_deferred_edges'])}")
            if got[0] != expected[0]:
                print("  kind:", got[0], "vs", expected[0], got[1] if got[0] == "raised" else "", expected[1] if expected[0] == "raised" else "")
            elif got[1] is None:
                print("  native declined the call")
            else:
                print(f"  isos native {len(got[1])} vs python {len(expected[1])}")
                for i, (a, b) in enumerate(zip(got[1], expected[1])):
                    if a != b:
                        for label, x, y in zip(("mapping", "deferred", "fragment", "symmetry"), a, b):
                            if x != y:
                                print(f"  iso {i} {label} differs:\n    native {str(x)[:400]}\n    python {str(y)[:400]}")
                        break
            if stop:
                break
    n = len(corpus["calls"])
    print(f"{corpus['case']}: {n} calls, {mismatches} mismatches, native time {t_native:.2f}s")
    sys.exit(1 if mismatches else 0)


if __name__ == "__main__":
    main()

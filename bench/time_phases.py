"""Time the phases of search_aam (parallel or serial) on a benchmark case."""
import sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "bench"))
from cases import CASES
import rxn_core.aam as aam_mod
from rxn_core.domain import AAMProblem, AAMSearchConfig

def main():
    if "--fork" in sys.argv:
        import multiprocessing
        multiprocessing.set_start_method("fork", force=True)
        sys.argv.remove("--fork")
    case = sys.argv[1] if len(sys.argv) > 1 else "tempo"
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    R, P = CASES[case]()
    problem = AAMProblem(R, P, name=case)
    phases = {}
    def timed(name, fn):
        def wrapper(*a, **k):
            t0 = time.perf_counter(); out = fn(*a, **k); phases[name] = phases.get(name, 0.0) + time.perf_counter() - t0
            return out
        return wrapper
    aam_mod.finalize_graph_symmetry = timed("exact_groups", aam_mod.finalize_graph_symmetry)
    t0 = time.perf_counter()
    result = aam_mod.search_aam(problem, AAMSearchConfig(), workers=workers)
    from rxn_core import group_mechanisms
    timed('optional_grouping', group_mechanisms)(result)
    total = time.perf_counter() - t0
    m = result.metrics
    print(f"case {case} workers {workers}: total {total:.2f}s | " + " | ".join(f"{k} {v:.2f}s" for k, v in phases.items())
          + f" | parent_merge {m.parent_merge_seconds:.3f}s | branches {m.retained_branch_count} | group_calcs {m.completed_group_calculations}")


if __name__ == "__main__":
    main()

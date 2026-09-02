"""cProfile one cut work unit (default: no-cut, all seed orders) of the TEMPO example."""
import cProfile, pstats, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rxn_core.cli import _endpoint_cache
from rxn_core.alignment.sweep import run_cut_sweep_chunk

R = _endpoint_cache(str(ROOT / "docs/example_runs/pr1.tempo_ts3/work/endpoints/R"), "R")
P = _endpoint_cache(str(ROOT / "docs/example_runs/pr1.tempo_ts3/work/endpoints/P"), "P")
cuts = [()] if len(sys.argv) < 2 else [()] + [((int(a), int(b)),) for a, b in (c.split("-") for c in sys.argv[1:])]
prof = cProfile.Profile()
t0 = time.perf_counter()
prof.enable()
pool, metrics = run_cut_sweep_chunk(
    list(R.elements), R.wbo, list(P.elements), P.wbo, cuts,
    n_workers=1, return_metrics=True, symmetry_wbo_tol=1.0)
prof.disable()
print(f"wall {time.perf_counter()-t0:.2f}s  mechanisms {len(pool)}  metrics {metrics}")
stats = pstats.Stats(prof)
print("\n=== top 35 by tottime ===")
stats.sort_stats("tottime").print_stats(35)
print("\n=== top 30 by cumtime ===")
stats.sort_stats("cumtime").print_stats(30)

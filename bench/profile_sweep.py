"""cProfile the whole serial search_aam on the TEMPO example; curated report."""
import cProfile, pstats, sys, time, io
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rxn_core.cli import _endpoint_cache
from rxn_core.domain import AAMProblem, AAMSearchConfig
from rxn_core.aam import search_aam

R = _endpoint_cache(str(ROOT / "docs/example_runs/pr1.tempo_ts3/work/endpoints/R"), "R")
P = _endpoint_cache(str(ROOT / "docs/example_runs/pr1.tempo_ts3/work/endpoints/P"), "P")
problem = AAMProblem(R, P, name="tempo")
prof = cProfile.Profile()
t0 = time.perf_counter()
prof.enable()
result = search_aam(problem, AAMSearchConfig(), workers=1)
prof.disable()
wall = time.perf_counter() - t0
print(f"wall {wall:.2f}s mechanisms {len(result.mechanisms)}")
stats = pstats.Stats(prof)
WATCH = [
    "sweep.py:1166(_run_cut_work)", "branch.py:574(find_islands)", "island.py:33(grow_island)",
    "extend.py:43(_extend_sym_cands)", "_collect_free_target_entries", "_dedupe_children",
    "_dedup_sym_cands", "_dedupe_certificates", "_boundary_signature", "canonical.py", "(certificate)",
    "(role_key)", "_p_relation_signature_from_parts", "_support_witness_for_value",
    "_children_from_context_group", "_children_from_block_join", "state.py:49(__init__)",
    "_score_branch_mapping", "symmetry_repair_mapping", "_mechanism_signature", "classify_bonds",
    "_MechanismEventCanonicalizer", "frag.py:146(build_graph)", "_nauty_orbits", "_nauty_atom_generators",
    "_wbo_tolerance_bucket_lookup", "_generate_seed_orders", "attach_completed_candidate_groups",
    "atom_generators", "_pool_add", "_analytical_branch", "_branch_analytical_derivations",
    "graph_cert", "autgrp", "_result_from_pool", "_branch_from_record", "expand_mapping",
    "_set_unique", "_frontier_boundary_edges", "_push_edges_from", "_progress_key", "commit", "fork",
    "with_automorph_equivalent", "_refine_sym_assignments", "_r_compatible_with_block", "_sym_block_indexes",
    "compatible", "_edge_wbo", "_node_attr",
]
buf = io.StringIO()
stats.stream = buf
stats.sort_stats("cumtime").print_stats(400)
lines = buf.getvalue().splitlines()
print(f"{'cumtime':>8} {'tottime':>8} {'ncalls':>9}  function")
seen = set()
for line in lines:
    for w in WATCH:
        if w in line and w not in seen:
            parts = line.split()
            if len(parts) >= 6:
                ncalls, tottime, _, cumtime = parts[0], parts[1], parts[2], parts[3]
                name = " ".join(parts[5:])
                print(f"{float(cumtime):8.2f} {float(tottime):8.2f} {ncalls:>9}  {name[-95:]}")
                seen.add(w)
            break

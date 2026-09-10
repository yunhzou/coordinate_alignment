"""Snapshot benchmark evidence; record a small, illustrative frozen-engine run.

The illustration is NOT added to benchmark scores. No reference labels are read
by the replay. Execute with the project Python, without plotting dependencies.
"""
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

MAN = Path(__file__).resolve().parents[1]
ROOT = MAN.parents[1]
REPORT = MAN.parent / 'reports'
RUNS = ROOT / 'aam_benchmarks'
OUT = MAN / 'evidence'
OUT.mkdir(exist_ok=True)
case = int(sys.argv[1]) if len(sys.argv) > 1 else 64

sources = {
    'seed_comparison.json': REPORT / 'aam_seed1_2_20260910/comparison.json',
    'slap_sweep.json': REPORT / 'slap_sweep_cut_20260910/summary.json',
    'publication.json': REPORT / 'golden_publication_20260908/analytics.json',
    'collection.json': REPORT / 'golden_pattern_benchmark_20260908/benchmark_report.json',
    'competitors.json': REPORT / 'golden_competitors_20260908/summary.json',
    'holdout_scores.json': RUNS / 'elementary140_tol1_20260909/tolerance_per_case.json',
    'holdout_overlap.json': RUNS / 'elementary140_tol1_20260909/overlap/summary.json',
    'holdout_manifest.json': RUNS / 'elementary140_tol1_20260909/manifest.json',
    'adaptive_no_sweep.json': RUNS / 'adaptive_full_20260910/comparison/summary.json',
    'adaptive_manifest.json': RUNS / 'adaptive_full_20260910/manifest.json',
    f'case{case}_input.json': RUNS / f'elementary140_tol1_20260909/inputs/{case}/input.json',
}
manifest = json.loads((OUT/'sources.json').read_text()) if (OUT/'sources.json').exists() else {}
for name, path in sources.items():
    payload = path.read_bytes()
    (OUT / name).write_bytes(payload)
    manifest[name] = {'source': str(path), 'sha256': hashlib.sha256(payload).hexdigest()}

engine = RUNS / 'aam_one_seed_bidirectional_20260910/original/src'
sys.path.insert(0, str(engine))
import rxn_core
# The frozen source snapshot stores growth binaries separately. Symmetry
# finalization uses the shared, unchanged bookkeeping extension; hash it below.
rxn_core.__path__.append(str(ROOT / 'src/rxn_core'))
from rxn_core.frag import build_graph
from rxn_core.alignment.branch import find_islands, _generate_seed_orders
from rxn_core.search_symmetry import finalize_graph_symmetry
from rxn_core.matcher import _nauty_orbits
import numpy as np

raw = json.loads((OUT / f'case{case}_input.json').read_text())
R, P = [build_graph(raw[k]['elements'], np.asarray(raw[k]['wbo']), bond_cut=.2)
        for k in ['reactant', 'product']]
order = next(iter(_generate_seed_orders(R, n_trials=1, rng_seed=42)))
events = []
g = find_islands(R, P, order, iso_tol=1., max_branches=100, events=events,
                p_orbits=_nauty_orbits(P, wbo_tol=1.), cuts=())
g, metrics = finalize_graph_symmetry(g, P, iso_tolerance=1.)
record = g.to_record()
out = {'case': case, 'name': raw['name'], 'input': raw, 'seed_order': order,
       'events': events, 'graph': record,
       'scope': 'Illustrative frozen-engine uncut replay, root seed 42, cap 100, '
                'matching tolerance 1.0, explicit H. Event trace follows the '
                'scheduler first branch only; graph stores all retained branches. '
                'Not a benchmark rerun, molecular dynamics, or reference mapping.',
       'engine_commit': '98b01b175eeed31f70d13e7cbf178b80bf07c9e0',
       'bookkeeping_binary_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in (ROOT/'src/rxn_core').glob('_group_ops*.so')},
       'engine_files': {str(p.relative_to(engine)): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in (engine / 'rxn_core').rglob('*.py')}}
(OUT / ('growth_trace.json' if case == 64 else f'growth_trace_{case}.json')).write_text(json.dumps(out, indent=2, default=lambda x: int(x))+'\n')
(OUT / 'sources.json').write_text(json.dumps(manifest, indent=2)+'\n')
print(json.dumps({'events': Counter(e['type'] for e in events),
                  'states': len(g.states), 'transitions': len(g.transitions),
                  'terminals': len(g.terminals), 'capped': g.capped}))

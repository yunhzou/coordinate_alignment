"""See what fragment sizes / uniqueness we get for a given step."""
import sys
from pathlib import Path
from collections import Counter
from rxn_core_frag import (
    run_xtb, build_graph, find_fragment_anchors, merge_anchors,
)

step = sys.argv[1] if len(sys.argv) > 1 else "pr13.Cyclobutane_JOC2023_TS-CD_step1"
root = Path("/Users/yunhengz/empty_for_claude/Benchmark") / step / "plain" / "stage0"
work = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_frag") / step

elR, xyzR, wboR = run_xtb(root / "reactant.xyz", work / "R")
elP, xyzP, wboP = run_xtb(root / "product.xyz", work / "P")
gR = build_graph(elR, wboR)
gP = build_graph(elP, wboP)

anchors = find_fragment_anchors(gR, gP, max_radius=4, min_radius=1, wbo_tol=0.15)
print(f"total anchors: {len(anchors)}")

size_counter = Counter(a['size'] for a in anchors)
unique_counter = Counter((a['size'], a['unique']) for a in anchors)
print("\nfragment-size distribution:")
for sz in sorted(size_counter):
    n_total = size_counter[sz]
    n_unique = unique_counter.get((sz, True), 0)
    n_amb = unique_counter.get((sz, False), 0)
    print(f"  size={sz:>2}  n={n_total:>4}  unique={n_unique:>4}  ambiguous={n_amb:>4}")

print(f"\nlargest 5 unique anchors:")
for a in sorted([x for x in anchors if x['unique']], key=lambda x: -x['size'])[:5]:
    print(f"  size={a['size']} radius={a['radius']} root R[{a['root_R']}]({elR[a['root_R']]}) -> P[{a['root_P']}]({elP[a['root_P']]})")
print(f"\nlargest 5 ambiguous anchors:")
for a in sorted([x for x in anchors if not x['unique']], key=lambda x: -x['size'])[:5]:
    print(f"  size={a['size']} radius={a['radius']} root R[{a['root_R']}]({elR[a['root_R']]}) -> P[{a['root_P']}]({elP[a['root_P']]})")

# What does merge produce if we only use UNIQUE anchors?
unique_only = [a for a in anchors if a['unique']]
m_unique = merge_anchors(unique_only, len(elR), len(elP))
m_all = merge_anchors(anchors, len(elR), len(elP))
print(f"\nmapping size with unique-only anchors: {len(m_unique)} / {len(elR)}")
print(f"mapping size with all anchors:         {len(m_all)} / {len(elR)}")

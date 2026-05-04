"""Inspect what the algorithm actually does on a step."""
import sys
from pathlib import Path
from rxn_core_wbo import (
    run_xtb, signatures, find_anchors, propagate, classify_bonds,
)

step = sys.argv[1] if len(sys.argv) > 1 else "pr1.tempo_ts1"
root = Path("/Users/yunhengz/empty_for_claude/Benchmark") / step / "plain" / "stage0"
work = Path("/Users/yunhengz/empty_for_claude/rxn_core/work") / step

elR, xyzR, wboR = run_xtb(root / "reactant.xyz", work / "R")
elP, xyzP, wboP = run_xtb(root / "product.xyz", work / "P")

sigs_R, nb_R = signatures(elR, wboR, max_radius=4)
sigs_P, nb_P = signatures(elP, wboP, max_radius=4)

anchors = find_anchors(sigs_R, sigs_P, elR, elP, 4, min_anchor_radius=2)
print("Anchors:")
for i, j, r in anchors:
    print(f"  R[{i}]({elR[i]}) <-> P[{j}]({elP[j]})  at radius {r}")

from rxn_core_wbo import cleanup_pass
mapping = propagate(anchors, sigs_R, sigs_P, nb_R, nb_P, elR, elP,
                    max_radius=4, min_match_radius=0)
mapping = cleanup_pass(mapping, nb_R, nb_P, elR, elP)

print(f"\nMapped {len(mapping)} / {len(elR)} atoms")
print("Mapping:")
for i in sorted(mapping):
    print(f"  R[{i:>2}]({elR[i]}) -> P[{mapping[i]:>2}]({elP[mapping[i]]})")
unmapped_R = [i for i in range(len(elR)) if i not in mapping]
inv = {v: k for k, v in mapping.items()}
unmapped_P = [j for j in range(len(elP)) if j not in inv]
print(f"Unmapped R: {unmapped_R}  ({[elR[i] for i in unmapped_R]})")
print(f"Unmapped P: {unmapped_P}  ({[elP[j] for j in unmapped_P]})")

print("\nNeighbors of unmapped R atoms (j, wbo, mapped?):")
for i in unmapped_R[:10]:
    nbs = [(j, round(w,2), mapping.get(j, "?")) for (j, w) in nb_R[i]]
    print(f"  R[{i}]({elR[i]}): {nbs}")

print("\nSignature counts at each radius:")
for r in range(5):
    from collections import Counter
    cR = Counter(sigs_R[r].values())
    cP = Counter(sigs_P[r].values())
    n_unique_R = sum(1 for v in cR.values() if v == 1)
    n_unique_P = sum(1 for v in cP.values() if v == 1)
    print(f"  r={r}: R unique sigs = {n_unique_R}/{len(elR)}, P unique sigs = {n_unique_P}/{len(elP)}")

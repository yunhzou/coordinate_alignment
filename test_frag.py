"""Test the fragment matcher on the three reference cases."""
import sys
import time
from pathlib import Path
from rxn_core_frag import analyze

BENCH = Path("/Users/yunhengz/empty_for_claude/Benchmark")
WORK = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_frag")

cases = sys.argv[1:] or [
    "pr1.tempo_ts1",
    "pr11.cycloadditions_tsIa",
    "pr13.Cyclobutane_JOC2023_TS-CD_step1",
]

for c in cases:
    t = time.time()
    r = BENCH / c / "plain" / "stage0" / "reactant.xyz"
    p = BENCH / c / "plain" / "stage0" / "product.xyz"
    try:
        res = analyze(r, p, WORK / c)
        print(f"\n=== {c} ({time.time()-t:.1f}s) ===")
        print(f"  N atoms: R={len(res['elements_R'])} P={len(res['elements_P'])}")
        print(f"  fragment anchors: {res['n_anchors']}")
        print(f"  after merge: {res['n_after_merge']}")
        print(f"  after expand: {len(res['mapping'])}")
        print(f"  broken bonds: {len(res['broken'])}")
        for (i, j, wR, wP) in res['broken'][:15]:
            wPs = '—' if wP is None else f"{wP:.2f}"
            print(f"    R[{i}]-R[{j}] WBO_R={wR:.2f} WBO_P={wPs}  "
                  f"({res['elements_R'][i]}-{res['elements_R'][j]})")
        print(f"  formed bonds: {len(res['formed'])}")
        for (i, j, wR, wP) in res['formed'][:15]:
            wRs = '—' if wR is None else f"{wR:.2f}"
            print(f"    P[{i}]-P[{j}] WBO_R={wRs} WBO_P={wP:.2f}  "
                  f"({res['elements_P'][i]}-{res['elements_P'][j]})")
        print(f"  core_R = {res['core_R']}")
        print(f"  core_P = {res['core_P']}")
    except Exception as e:
        print(f"\n=== {c} FAILED ({time.time()-t:.1f}s) ===")
        print(f"  {type(e).__name__}: {e}")
        import traceback; traceback.print_exc()

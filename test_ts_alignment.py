"""
Test: align a TS candidate to R and P using the same algorithm, then
verify the R→TS→P bond evolution is internally consistent.

For each step, we take the rank-0 candidate TS, run:
   analyze(R, TS) -> mapping R-idx to TS-idx
   analyze(TS, P) -> mapping TS-idx to P-idx
   analyze(R, P)  -> mapping R-idx to P-idx (reference)
and check that the composed R→TS→P mapping matches R→P.

Atom indices in the TS xyz are NOT assumed to match R/P — alignment is
exactly the mapping algorithm's job.
"""
from __future__ import annotations
import re
import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from rxn_core_frag import analyze
from build_tsdisco_viewer import step_inputs, concat_xyz, _parse_xyz_text


TSDISCO_INDEX = Path(
    "/Users/yunhengz/empty_for_claude/tsdisco_benchmark_visualization_plain_portable/index.html"
)
WORK = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_tsalign")
WORK.mkdir(parents=True, exist_ok=True)


def load_data():
    html = TSDISCO_INDEX.read_text()
    m = re.search(r"const TSDISCO_DATA = (\{.*?\});\s*\n", html, re.DOTALL)
    return json.loads(m.group(1))


def scramble_xyz(xyz_text, seed=42):
    """Randomly permute atom order in a single-fragment xyz (worst-case
    atom labeling for the alignment algorithm)."""
    elements, coords = _parse_xyz_text(xyz_text)
    n = len(elements)
    perm = list(range(n))
    random.Random(seed).shuffle(perm)
    body = "\n".join(f"{elements[p]}  {coords[p][0]:.6f}  {coords[p][1]:.6f}  {coords[p][2]:.6f}"
                     for p in perm)
    return f"{n}\nscrambled\n{body}\n", perm


def test_step(step, scramble_ts=True):
    name = f"{step['dataset']}/{step['step_id']}"
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)
    chg = step.get("charge", 0) or 0
    uhf = max(0, (step.get("multiplicity", 1) or 1) - 1)

    # Reactant and product
    rxyz_text, pxyz_text, _, _ = step_inputs(step)
    (wd / "reactant.xyz").write_text(rxyz_text)
    (wd / "product.xyz").write_text(pxyz_text)

    # TS = best candidate (rank 0). Sometimes called "passed" candidates.
    cands = step.get("candidates", [])
    if not cands:
        return f"  no candidates"
    ts = cands[0]
    ts_xyz = ts.get("xyz", "")
    if not ts_xyz:
        return f"  rank-0 cand has no xyz"

    # Optionally scramble the TS atom ordering to test the alignment
    perm = None
    if scramble_ts:
        ts_xyz, perm = scramble_xyz(ts_xyz)

    (wd / "ts.xyz").write_text(ts_xyz)

    # 1. R->P reference
    res_RP = analyze(wd / "reactant.xyz", wd / "product.xyz", wd / "RP",
                     charge=chg, uhf=uhf)
    rp_br = len(res_RP["broken"]); rp_fm = len(res_RP["formed"])

    # 2. R->TS
    try:
        res_RT = analyze(wd / "reactant.xyz", wd / "ts.xyz", wd / "RT",
                         charge=chg, uhf=uhf)
        rt_br = len(res_RT["broken"]); rt_fm = len(res_RT["formed"])
        rt_mapped = len(res_RT["mapping"])
    except Exception as e:
        return f"  R->TS FAILED: {e}"

    # 3. TS->P
    try:
        res_TP = analyze(wd / "ts.xyz", wd / "product.xyz", wd / "TP",
                         charge=chg, uhf=uhf)
        tp_br = len(res_TP["broken"]); tp_fm = len(res_TP["formed"])
        tp_mapped = len(res_TP["mapping"])
    except Exception as e:
        return f"  TS->P FAILED: {e}"

    # 4. Compose R -> TS -> P and check against R -> P
    m_RT = res_RT["mapping"]
    m_TP = res_TP["mapping"]
    m_RP = res_RP["mapping"]
    composed = {r: m_TP[m_RT[r]] for r in m_RT if m_RT[r] in m_TP}
    matches = sum(1 for r, p in composed.items() if m_RP.get(r) == p)
    n_compose = len(composed)
    consistency = f"{matches}/{n_compose}" if n_compose else "0/0"

    scrambled_tag = " (TS scrambled)" if scramble_ts else ""
    return (f"  R->P {rp_br}/{rp_fm}  R->TS {rt_br}/{rt_fm} "
            f"(map={rt_mapped})  TS->P {tp_br}/{tp_fm} (map={tp_mapped})  "
            f"composed-vs-direct: {consistency}{scrambled_tag}")


def main():
    test_step_ids = sys.argv[1:] or [
        "Benchmark/pr1.tempo_ts2",
        "Benchmark/pr11.cycloadditions_tsIa",
        "Benchmark/pr12.Co_Silylation_JACS2015_TS_Dstar-Estar",
        "Benchmark/pr14.Pd_hydroamination_JOC2025_TS14_step4_reductive_elimination",
        "Benchmark/pr9.carbene.rearr_ts47a",
    ]
    data = load_data()
    by_id = {f"{s['dataset']}/{s['step_id']}": s for s in data["steps"]}

    for sid in test_step_ids:
        step = by_id.get(sid)
        if not step:
            print(f"{sid}  NOT FOUND")
            continue
        t = time.time()
        result = test_step(step, scramble_ts=True)
        print(f"{sid}  ({time.time()-t:.1f}s):")
        print(result)


if __name__ == "__main__":
    main()

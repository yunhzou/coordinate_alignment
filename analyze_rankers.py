"""
A/B different ranker formulas on the per-IG records produced by
evaluate_alignment.py.

The fixed-aligner output contains, per step, every IG with its picked
imag mode's (beta, rho, kappa, n_imag) plus its mass-weighted cosine
to GT's picked imag mode (the "ground-truth signal" we want to rank
toward).

For each ranker formula:
  - sort IGs by that formula per step
  - take top-1 (and top-2's max)
  - record top-1's mwc_to_GT
Aggregate over the 155 steps.

Question: does the full formula's complexity (rho, kappa, n_imag^p
factors) actually buy anything over beta alone, now that the alignment
is correct?

Usage:
  python analyze_rankers.py
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np


PROJECT = Path(__file__).resolve().parent
IN_JSON = PROJECT / "out" / "bgcp_alignment_eval.json"


def s_beta(r):              return r["beta"]
def s_beta_div_nimag(r):    return r["beta"] / max(r["n_imag"], 1) ** 0.3
def s_beta_x_rho(r):        return r["beta"] * (1 + r["rho"])
def s_beta_x_kappa(r):      return r["beta"] * (1 + 0.2 * r["kappa"])
def s_beta_x_rho_kappa(r):  return r["beta"] * (1 + r["rho"]) * (1 + 0.2 * r["kappa"])
def s_full(r):
    return r["beta"] * (1 + r["rho"]) * (1 + 0.2 * r["kappa"]) / max(r["n_imag"], 1) ** 0.3


RANKERS = [
    ("beta only",                        s_beta),
    ("beta / n_imag^0.3",                s_beta_div_nimag),
    ("beta * (1 + rho)",                 s_beta_x_rho),
    ("beta * (1 + 0.2 kappa)",           s_beta_x_kappa),
    ("beta * (1+rho)(1+0.2 kappa)",      s_beta_x_rho_kappa),
    ("FULL: beta(1+rho)(1+0.2k)/n^0.3",  s_full),
]


def main():
    rows = json.loads(IN_JSON.read_text())
    ok = [r for r in rows if "error" not in r or not r["error"]]
    n = len(ok)
    print(f"loaded {n} steps\n")

    # Oracle: best possible top-1 mwc for each step
    oracle = []
    for rec in ok:
        valid = [g for g in rec["igs"] if "beta" in g]
        if not valid: continue
        oracle.append(max(g["mwc_to_GT"] for g in valid))
    oracle = np.array(oracle)

    print(f"Oracle (best mwc-to-GT among 20 IGs per step):")
    print(f"  mean={oracle.mean():.3f}  median={np.median(oracle):.3f}  "
          f">0.7: {(oracle > 0.7).sum()}/{len(oracle)}  "
          f">0.5: {(oracle > 0.5).sum()}/{len(oracle)}")
    print()

    print(f"{'ranker':40s}  {'top-1 mwc (mean)':>16s}  "
          f"{'med':>5s}  {'top-1 >0.7':>10s}  {'top-1 >0.5':>10s}  "
          f"{'top-2 best mwc (mean)':>22s}  {'>0.7':>5s}")
    print('-' * 130)
    for name, fn in RANKERS:
        t1_mwc = []; t2_best_mwc = []
        for rec in ok:
            valid = [g for g in rec["igs"] if "beta" in g]
            if not valid: continue
            ranked = sorted(valid, key=lambda r: -fn(r))
            top1 = ranked[0]
            top2 = ranked[1] if len(ranked) > 1 else None
            t1_mwc.append(top1["mwc_to_GT"])
            t2_best_mwc.append(max(top1["mwc_to_GT"],
                                   top2["mwc_to_GT"] if top2 else 0.0))
        a1 = np.array(t1_mwc); a2 = np.array(t2_best_mwc)
        print(f"{name:40s}  {a1.mean():>16.3f}  {np.median(a1):>5.3f}  "
              f"{(a1>0.7).sum():>3d}/{len(a1):<3d}  "
              f"{(a1>0.5).sum():>3d}/{len(a1):<3d}  "
              f"{a2.mean():>22.3f}  {(a2>0.7).sum():>3d}/{len(a2)}")

    # Spearman correlation between each ranker and the GT-similarity
    # signal (per-IG, pooled across all steps)
    print()
    print("Ranker vs mwc-to-GT (Spearman) -- pooled across all (step, IG) pairs")
    pool = [(g, fn_for_n) for rec in ok for g in rec["igs"] if "beta" in g
            for fn_for_n in [None]]  # just collect IGs
    igs = [g for rec in ok for g in rec["igs"] if "beta" in g]
    mwc_arr = np.array([g["mwc_to_GT"] for g in igs])
    from scipy.stats import spearmanr
    print(f"  N(IG, step) pairs = {len(igs)}\n")
    for name, fn in RANKERS:
        scores = np.array([fn(g) for g in igs])
        rho, p = spearmanr(scores, mwc_arr)
        print(f"  {name:40s}  Spearman rho = {rho:+.3f}  p={p:.2e}")


if __name__ == "__main__":
    main()

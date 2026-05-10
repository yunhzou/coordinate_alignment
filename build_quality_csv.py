"""
Build a merged quality-measurement CSV that combines:

  - per-step verifier top-1 / top-2 / oracle mwc-to-GT (picked-imag vs
    picked-imag, both reindexed to R-frame), under the post-_set_unique-fix
    alignment, computed from out/bgcp_alignment_eval.json
  - human strict and lenient annotations (IG#1, IG#2) from the user's
    existing final_quality_measurement-strict-version.csv

GROUND TRUTH = human strict labels (expert computational-chemist judgment
of whether the picked imaginary mode actually represents the reaction).
mwc is INFORMATIONAL ONLY -- it's a loose numerical similarity that can
be moderate (0.4-0.6) even for chemically-correct guesses because of
vibrational mixing with spectator modes, rotational coupling, etc. Do
NOT interpret low mwc as "wrong"; consult the human label.

mwc semantics: for each IG, we compare ONLY its picked imag mode
(max-beta imag) against GT's picked imag mode. One value per IG. The
mwc_>=_threshold columns are kept as a coarse signal but should not
override the human label.

Output: out/bgcp_quality_merged.csv

Usage:
  python build_quality_csv.py [--threshold 0.7]
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path


PROJECT = Path(__file__).resolve().parent
IN_JSON = PROJECT / "out" / "bgcp_alignment_eval.json"
HUMAN_CSV = Path("/Users/yunhengz/Downloads/final_quality_measurement-strict-version.csv")
OUT_CSV = PROJECT / "out" / "bgcp_quality_merged.csv"


W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3


def full_score(r):
    n = max(r.get("n_imag", 0), 1)
    return (r.get("beta", 0.0)
            * (1 + W_RXN * r.get("rho", 0.0))
            * (1 + W_CORE * r.get("kappa", 0.0))
            / n ** IMAG_PEN)


def load_human(path: Path) -> dict:
    """Returns {step: {'IG1_strict':..., 'IG2_strict':...,
                       'IG1_lenient':..., 'IG2_lenient':...,
                       'old_top1_picked_mwc': float}}.
    The reference CSV has two duplicated header rows; column order:
    step, IG#1 (strict), IG#2 (strict), IG#1 (lenient), IG#2 (lenient), ...
    column 47 is verifier_top1_picked = old picked-imag mwc-to-GT."""
    out = {}
    if not path.exists():
        print(f"  (no human CSV at {path}; skipping merge)")
        return out
    with path.open() as f:
        rows = list(csv.reader(f))
    for r in rows[2:]:
        if len(r) < 48 or not r[0].strip():
            continue
        step = r[0].strip()
        try:
            out[step] = {
                "IG1_strict":  int(r[1]) if r[1].strip() else None,
                "IG2_strict":  int(r[2]) if r[2].strip() else None,
                "IG1_lenient": int(r[3]) if r[3].strip() else None,
                "IG2_lenient": int(r[4]) if r[4].strip() else None,
                "old_top1_picked_mwc": (float(r[47]) if r[47].strip() else None),
            }
        except ValueError:
            pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=0.7)
    ap.add_argument("--in", dest="inp", default=str(IN_JSON))
    ap.add_argument("--out", default=str(OUT_CSV))
    args = ap.parse_args()

    data = json.loads(Path(args.inp).read_text())
    human = load_human(HUMAN_CSV)
    print(f"loaded {len(data)} alignment records, {len(human)} human-annotated steps")

    # Lean schema: only the columns used in chemist-vs-metric analysis.
    # For the full 20-IG sorted-by-mwc breakdown, see out/bgcp_alignment_eval.json.
    # _top2_* is the BEST among verifier's top-1 and top-2 picks.
    # case = bucket the row falls into (see comment below in main loop).
    fields = [
        "step",
        "gt_score",
        "t1_mwc", "t1_Sr",
        "t2_mwc", "t2_Sr",
        "oracle_mwc",
        "h_IG1_strict",  "h_IG2_strict",
        "h_IG1_lenient", "h_IG2_lenient",
        "case",
    ]

    n_written = 0
    n_loose = 0              # human pass + mwc below threshold (loose-mwc cases)
    n_human_pass1_strict = 0
    n_top1_mwc_ge = 0
    n_top20_mwc_ge = 0
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rec in data:
            if rec.get("error"): continue
            step = rec["step"]
            igs = [g for g in rec.get("igs", []) if "beta" in g]
            if not igs: continue
            n_ig = len(igs)
            sorted_by_score = sorted(igs, key=lambda r: -full_score(r))
            sorted_by_mwc   = sorted(igs, key=lambda r: -r["mwc_to_GT"])
            t1 = sorted_by_score[0]
            t2 = sorted_by_score[1] if len(sorted_by_score) > 1 else None
            oracle = sorted_by_mwc[0]

            t1_mwc = t1["mwc_to_GT"]
            t2_mwc = t2["mwc_to_GT"] if t2 else 0.0
            best_t2_mwc = max(t1_mwc, t2_mwc) if t2 else t1_mwc
            oracle_mwc = oracle["mwc_to_GT"]
            top1_ge = 1 if t1_mwc      >= args.threshold else 0
            top2_ge = 1 if best_t2_mwc >= args.threshold else 0
            top20_ge = 1 if oracle_mwc >= args.threshold else 0

            h = human.get(step, {})
            ig1s = h.get("IG1_strict"); ig2s = h.get("IG2_strict")
            ig1l = h.get("IG1_lenient"); ig2l = h.get("IG2_lenient")
            human_p1s = ig1s if ig1s is not None else ""
            human_p2s = (1 if ((ig1s == 1) or (ig2s == 1)) else 0) if (ig1s is not None or ig2s is not None) else ""
            human_p1l = ig1l if ig1l is not None else ""
            human_p2l = (1 if ((ig1l == 1) or (ig2l == 1)) else 0) if (ig1l is not None or ig2l is not None) else ""

            # Discrepancy: human strict said pass@1 but mwc<threshold. NOT
            # human error. Tracks loose-mwc cases for threshold tuning.
            loose = ""
            if human_p1s != "":
                loose = 1 if (human_p1s == 1 and top1_ge == 0) else 0
                if human_p1s == 1: n_human_pass1_strict += 1
                if loose == 1: n_loose += 1
            if top1_ge == 1: n_top1_mwc_ge += 1
            if top20_ge == 1: n_top20_mwc_ge += 1

            S_gt = rec.get("gt_score", None)
            def s_ratio(ig):
                if S_gt is None or S_gt < 1e-9: return ""
                return round(full_score(ig) / S_gt, 4)

            t1_sr = s_ratio(t1)
            t2_sr = s_ratio(t2) if t2 else ""

            # case bucket. Priority: D > C > E > A > B > OK > unlabeled.
            # D = chemist=FAIL, verifier S_ratio>=2.0  (ranker wildly inflated)
            # C = chemist=FAIL, t1_mwc>=0.7 and S_ratio>=0.9  (verifier confident, chemist disagrees)
            # E = chemist split: IG1=fail, IG2=pass (rank inversion)
            # A = chemist=PASS, t1_mwc<0.5 and S_ratio<0.8 (real ranker miss)
            # B = chemist=PASS, t1_mwc<0.7 and S_ratio in [0.8, 1.5) (different mode, equiv chemistry)
            # OK = chemist=PASS, both metrics healthy
            # ""  = no human label
            sr_v = float(t1_sr) if t1_sr != "" else None
            case = ""
            if ig1s is not None:
                if ig1s == 0 and sr_v is not None and sr_v >= 2.0:
                    case = "D_inflation"
                elif ig1s == 0 and t1_mwc >= 0.7 and sr_v is not None and sr_v >= 0.9:
                    case = "C_chemist_fail_metric_pass"
                elif ig1s == 0 and ig2s == 1:
                    case = "E_rank_inversion"
                elif ig1s == 1 and t1_mwc < 0.5 and sr_v is not None and sr_v < 0.8:
                    case = "A_ranker_miss"
                elif ig1s == 1 and t1_mwc < 0.7 and sr_v is not None and 0.8 <= sr_v < 1.5:
                    case = "B_diff_mode_eq_chem"
                elif ig1s == 1:
                    case = "OK"

            row = {
                "step": step,
                "gt_score": round(S_gt, 6) if S_gt is not None else "",
                "t1_mwc":     round(t1_mwc, 4),
                "t1_Sr":      t1_sr,
                "t2_mwc":     round(best_t2_mwc, 4),
                "t2_Sr":      t2_sr,
                "oracle_mwc": round(oracle_mwc, 4),
                "h_IG1_strict":  ig1s if ig1s is not None else "",
                "h_IG2_strict":  ig2s if ig2s is not None else "",
                "h_IG1_lenient": ig1l if ig1l is not None else "",
                "h_IG2_lenient": ig2l if ig2l is not None else "",
                "case": case,
            }
            w.writerow(row)
            n_written += 1

    print(f"\nwrote {args.out}  ({n_written} steps, mwc threshold={args.threshold})")
    print(f"\nGround truth = human strict label.")
    print(f"  human pass@1 strict = {n_human_pass1_strict}/{n_written}  "
          f"({100*n_human_pass1_strict/n_written:.1f}%)")
    print(f"\nmwc-only flag (informational; mwc is loose, not a pass criterion):")
    print(f"  top-1 mwc >= {args.threshold}  : {n_top1_mwc_ge}/{n_written}  "
          f"({100*n_top1_mwc_ge/n_written:.1f}%)")
    print(f"  top-20 (oracle) mwc >= {args.threshold}  : {n_top20_mwc_ge}/{n_written}  "
          f"({100*n_top20_mwc_ge/n_written:.1f}%)")
    print(f"\nLoose-mwc cases (human said pass, mwc below {args.threshold}):")
    print(f"  count = {n_loose}/{n_human_pass1_strict}  "
          f"-- shows mwc threshold is too strict for these chemically-correct guesses")


if __name__ == "__main__":
    main()

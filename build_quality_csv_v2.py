"""Build merged quality CSV from v2 (multi-mechanism) eval JSON.

For each step, picks the "winning" mechanism (highest GT S) and reports
its top-1/top-2 IGs. Compares against human strict + lenient labels.

Schema (lean):
  step, n_mechs, winning_mech_cut,
  gt_S, gt_beta, gt_kappa,
  t1_label, t1_S, t1_beta, t1_kappa,
  t2_label, t2_S,
  union_top2_labels   (concat across mechs, dedup)
  h_IG1_strict, h_IG2_strict, h_IG1_lenient, h_IG2_lenient,
  case   (A_ranker_miss, B_diff_mech_eq_chem, etc. — same scheme as v1)
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

PROJECT = Path(__file__).resolve().parent
IN_JSON = PROJECT / "out" / "bgcp_alignment_eval_v2.json"
HUMAN_CSV = Path("/Users/yunhengz/Downloads/final_quality_measurement-strict-version.csv")
OUT_CSV = PROJECT / "out" / "bgcp_quality_merged_v2.csv"


def load_human(path):
    out = {}
    if not path.exists():
        print(f"  (no human CSV at {path})")
        return out
    with path.open() as f:
        rows = list(csv.reader(f))
    for r in rows[2:]:
        if len(r) < 48 or not r[0].strip(): continue
        step = r[0].strip()
        try:
            out[step] = {
                "IG1_strict":  int(r[1]) if r[1].strip() else None,
                "IG2_strict":  int(r[2]) if r[2].strip() else None,
                "IG1_lenient": int(r[3]) if r[3].strip() else None,
                "IG2_lenient": int(r[4]) if r[4].strip() else None,
            }
        except ValueError: pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=str(IN_JSON))
    ap.add_argument("--out", default=str(OUT_CSV))
    args = ap.parse_args()

    data = json.loads(Path(args.inp).read_text())
    human = load_human(HUMAN_CSV)
    print(f"loaded {len(data)} eval records, {len(human)} human-annotated steps")

    fields = ["step", "n_mechs", "winning_cut",
              "gt_S", "gt_beta", "gt_kappa",
              "t1_label", "t1_S", "t1_beta", "t1_kappa",
              "t2_label", "t2_S",
              "union_top_labels",
              "h_IG1_strict", "h_IG2_strict",
              "h_IG1_lenient", "h_IG2_lenient",
              "case"]
    n_written = 0
    counts_case = {}
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rec in data:
            if rec.get("error"): continue
            step = rec["step"]
            mechs = rec["mechanisms"]
            # Pick winning mech by GT S
            wmech = max(mechs, key=lambda m: m["gt"]["S"] if m.get("gt") else 0)
            gt = wmech.get("gt") or {}
            igs_sorted = sorted([ig for ig in wmech["igs"] if ig.get("S") is not None],
                                key=lambda x: -x["S"])
            t1 = igs_sorted[0] if igs_sorted else {}
            t2 = igs_sorted[1] if len(igs_sorted) > 1 else {}
            # union top-2 across mechs
            union = set()
            for m in mechs:
                ranked = sorted([ig for ig in m["igs"] if ig.get("S") is not None],
                                key=lambda x: -x["S"])
                for ig in ranked[:2]: union.add(ig["label"])
            h = human.get(step, {})
            ig1s = h.get("IG1_strict"); ig2s = h.get("IG2_strict")
            ig1l = h.get("IG1_lenient"); ig2l = h.get("IG2_lenient")

            # Case (compared against pass@1 strict expectation)
            t1_label = t1.get("label", "")
            t1_S = t1.get("S")
            case = ""
            if ig1s is not None:
                if ig1s == 1 and t1_S is not None and t1_S >= 0.5:
                    case = "OK_match"
                elif ig1s == 1 and t1_S is not None and t1_S < 0.5:
                    case = "A_low_S"
                elif ig1s == 0 and t1_S is not None and t1_S >= 1.0:
                    case = "C_high_S_chem_fail"
                elif ig1s == 0 and ig2s == 1:
                    case = "E_rank_inversion"
                elif ig1s == 1:
                    case = "B_diff"
                else:
                    case = "OK_fail"
            counts_case[case] = counts_case.get(case, 0) + 1

            w.writerow({
                "step": step, "n_mechs": len(mechs),
                "winning_cut": wmech.get("cut", ""),
                "gt_S": round(gt.get("S", 0), 4) if gt else "",
                "gt_beta": round(gt.get("beta", 0), 4) if gt else "",
                "gt_kappa": round(gt.get("kappa", 0), 4) if gt else "",
                "t1_label": t1_label, "t1_S": round(t1_S, 4) if t1_S is not None else "",
                "t1_beta": round(t1.get("beta", 0), 4) if t1 else "",
                "t1_kappa": round(t1.get("kappa", 0), 4) if t1 else "",
                "t2_label": t2.get("label", ""),
                "t2_S": round(t2["S"], 4) if t2.get("S") is not None else "",
                "union_top_labels": "|".join(sorted(union)),
                "h_IG1_strict": ig1s if ig1s is not None else "",
                "h_IG2_strict": ig2s if ig2s is not None else "",
                "h_IG1_lenient": ig1l if ig1l is not None else "",
                "h_IG2_lenient": ig2l if ig2l is not None else "",
                "case": case,
            })
            n_written += 1

    print(f"\nwrote {args.out}  ({n_written} steps)")
    print(f"\ncase distribution:")
    for k, v in sorted(counts_case.items(), key=lambda x: -x[1]):
        print(f"  {k or '(no human)':<30}  {v}")


if __name__ == "__main__":
    main()

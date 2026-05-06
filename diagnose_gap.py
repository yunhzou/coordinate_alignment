"""
For each step, decompose the gap between aggressive_v1 (k=2) and oracle:

  oracle = max over all (IG, mode) of gt_align
  oracle_imag = max over all (IG, imag-mode) of gt_align
  oracle_among_top2_IGs = max over (top-2 IG, all modes) of gt_align
  oracle_among_top2_IGs_imag = same but imag only
  ranker_top2 = max over (top-2 IG, picked-mode-by-bond) of gt_align

If oracle_among_top2_IGs == oracle, the gap is "wrong-mode-within-right-IG".
If oracle_among_top2_IGs < oracle, the gap is "wrong-IG-picked".
"""
from __future__ import annotations
import json, re, time
from pathlib import Path
import numpy as np

from improve_ranker import (load_step, imag_modes, mass_weighted_cos,
                            rk_aggressive_v1, evaluate_topk, cos_sim)

SRC = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_viewer')


def main():
    files = sorted(SRC.glob('*.html'))
    files = [f for f in files
             if f.name not in ('index.html', 'flat_view.html', 'guess_quality.html')]
    print(f"Loading {len(files)} steps...")
    t0 = time.time()
    scope = []
    for hp in files:
        sd = load_step(hp)
        if sd: scope.append(sd)
    print(f"  loaded {len(scope)} steps in {time.time()-t0:.0f}s\n")

    rows = []
    for sd in scope:
        gt_disp = sd['gt_disp']
        oracle_all = max((cos_sim(np.asarray(m['disp']), gt_disp)
                          for ts in sd['igs'] for m in ts['modes']), default=0)
        oracle_imag = max((cos_sim(np.asarray(m['disp']), gt_disp)
                           for ts in sd['igs']
                           for m in ts['modes'] if m['freq'] < 0), default=0)
        # aggressive_v1 top-2
        ranked = rk_aggressive_v1(sd, 0.7, 0.10, 1.0, 0.2)
        top_aligns = [cos_sim(np.asarray(p[1]['disp']), gt_disp) for p in ranked[:2]]
        ranker_top2 = max(top_aligns) if top_aligns else 0
        top_labels = set(p[2] for p in ranked[:2])
        # Find oracle-best IG label
        best_ig_label, best_ig_align = None, 0
        for ts in sd['igs']:
            for m in ts['modes']:
                a = cos_sim(np.asarray(m['disp']), gt_disp)
                if a > best_ig_align:
                    best_ig_align = a; best_ig_label = ts['label']
        # Oracle restricted to top-2 IGs
        oracle_top2_IGs = 0
        oracle_top2_IGs_imag = 0
        for ts in sd['igs']:
            if ts['label'] not in top_labels: continue
            for m in ts['modes']:
                a = cos_sim(np.asarray(m['disp']), gt_disp)
                oracle_top2_IGs = max(oracle_top2_IGs, a)
                if m['freq'] < 0:
                    oracle_top2_IGs_imag = max(oracle_top2_IGs_imag, a)
        rows.append(dict(
            step=sd['step'], oracle=oracle_all, oracle_imag=oracle_imag,
            oracle_top2=oracle_top2_IGs, oracle_top2_imag=oracle_top2_IGs_imag,
            ranker=ranker_top2, top_labels=','.join(sorted(top_labels)),
            best_ig=best_ig_label,
            wrong_ig=int(best_ig_label not in top_labels),
            wrong_mode=int(best_ig_label in top_labels and ranker_top2 < best_ig_align - 0.05),
        ))

    # Summary
    n = len(rows)
    arr_oracle = np.array([r['oracle'] for r in rows])
    arr_o_imag = np.array([r['oracle_imag'] for r in rows])
    arr_o_t2 = np.array([r['oracle_top2'] for r in rows])
    arr_o_t2_imag = np.array([r['oracle_top2_imag'] for r in rows])
    arr_rank = np.array([r['ranker'] for r in rows])

    print(f"{'metric':32s}  mean   ≥0.7    ≥0.5    ≥0.3")
    for name, a in (('oracle (all modes)', arr_oracle),
                    ('oracle_imag (imag modes only)', arr_o_imag),
                    ('oracle restricted to top-2 IGs', arr_o_t2),
                    ('oracle in top-2 IGs (imag only)', arr_o_t2_imag),
                    ('ranker (aggressive_v1 top-2)', arr_rank)):
        print(f"{name:32s}  {a.mean():.3f}  {(a>=0.7).mean()*100:5.1f}%  "
              f"{(a>=0.5).mean()*100:5.1f}%  {(a>=0.3).mean()*100:5.1f}%")

    # Decompose: where is the gap?
    wrong_ig = sum(r['wrong_ig'] for r in rows)
    wrong_mode = sum(r['wrong_mode'] for r in rows)
    correct = n - wrong_ig - wrong_mode
    print(f"\nGap decomposition (vs oracle):")
    print(f"  Both right (top-2 IGs include best, mode within 0.05 of best): {correct}/{n}")
    print(f"  Wrong IG (best IG not in top-2):                                 {wrong_ig}/{n}")
    print(f"  Right IG, wrong mode (best IG in top-2, but >0.05 mode gap):     {wrong_mode}/{n}")
    print()
    print(f"How far off is ranker when WRONG IG?")
    wig = [r for r in rows if r['wrong_ig']]
    if wig:
        ws = np.array([r['oracle']-r['ranker'] for r in wig])
        print(f"  N={len(wig)}, gap mean={ws.mean():.3f}, median={np.median(ws):.3f}")
    print(f"How far off when RIGHT IG, WRONG MODE?")
    wm = [r for r in rows if r['wrong_mode']]
    if wm:
        ws = np.array([r['oracle']-r['ranker'] for r in wm])
        print(f"  N={len(wm)}, gap mean={ws.mean():.3f}, median={np.median(ws):.3f}")

    # Save details for steps where best IG was missed
    out_csv = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_analysis/gap_diagnosis.csv')
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    import csv
    with out_csv.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\nCSV: {out_csv}")


if __name__ == '__main__':
    main()

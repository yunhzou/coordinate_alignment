"""
Per-IG oracle alignment CSV.

For every (step, initial_guess) pair compute:
  - best_align_any:  max gt_alignment over ALL modes of the IG (oracle, all)
  - best_align_imag: max gt_alignment restricted to imaginary modes
  - best_mode_idx, best_mode_freq, best_mode_is_imag
  - bond_overlap / rxn_overlap / core_fraction of the best-aligned mode

Also a per-step summary row (one per step) reporting:
  - oracle_step (max over all 20 IGs × all modes — the absolute ceiling)
  - oracle_step_imag (max over imaginary modes only)
  - n_ig
  - which IG label achieves oracle_step

Output:
  out/mode_analysis/oracle_alignment_per_ig.csv  (3200 rows)
  out/mode_analysis/oracle_alignment_per_step.csv (160 rows)
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import csv
import json
import re
import time
from pathlib import Path
import numpy as np

SRC = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_viewer')
OUT_PER_IG  = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_analysis/oracle_alignment_per_ig.csv')
OUT_PER_STEP = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_analysis/oracle_alignment_per_step.csv')
OUT_TOPK = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_analysis/oracle_topk_per_step.csv')


def cos_sim(a, b):
    a = np.asarray(a).reshape(-1); b = np.asarray(b).reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9: return 0.0
    return abs(float(a @ b)) / (na * nb)


def main():
    files = sorted(SRC.glob('*.html'))
    files = [f for f in files
             if f.name not in ('index.html', 'flat_view.html', 'oracle_view.html',
                               'guess_quality.html')]
    print(f"Reading {len(files)} steps...")
    t0 = time.time()
    per_ig_rows = []
    per_step_rows = []
    for hp in files:
        text = hp.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        gt = next((t for t in data['ts_list']
                   if t['label']=='groundtruth' and t.get('modes')), None)
        if gt is None: continue
        gt_disp = np.asarray(gt['modes'][gt['default_mode_idx']]['disp'])
        gt_freq = gt['modes'][gt['default_mode_idx']]['freq']

        igs = [t for t in data['ts_list']
               if t['label']!='groundtruth' and t.get('modes')]

        oracle_step_any = 0.0
        oracle_step_imag = 0.0
        oracle_ig_label = ''
        oracle_mode_idx = -1

        for ts in igs:
            best_any = (-1.0, None)
            best_imag = (-1.0, None)
            for m in ts['modes']:
                a = cos_sim(np.asarray(m['disp']), gt_disp)
                if a > best_any[0]:  best_any  = (a, m)
                if m['freq'] < 0 and a > best_imag[0]:
                    best_imag = (a, m)
            ba, bm = best_any
            bia, bim = best_imag
            per_ig_rows.append(dict(
                step=data['step'],
                ig_label=ts['label'],
                n_modes=ts['n_modes_total'],
                n_imag=ts['n_imag'],
                # Best ANY mode
                best_align_any=round(ba, 6),
                best_any_idx=bm['idx'] if bm else -1,
                best_any_freq=round(bm['freq'], 4) if bm else 0,
                best_any_is_imag=int(bm['freq']<0) if bm else 0,
                best_any_bond_ov=round(bm.get('bond_overlap',0), 4) if bm else 0,
                best_any_rxn_ov=round(bm.get('rxn_overlap',0), 4) if bm else 0,
                best_any_core_frac=round(bm.get('core_fraction',0), 4) if bm else 0,
                # Best IMAG mode
                best_align_imag=round(bia, 6) if bim else 0,
                best_imag_idx=bim['idx'] if bim else -1,
                best_imag_freq=round(bim['freq'], 4) if bim else 0,
                best_imag_bond_ov=round(bim.get('bond_overlap',0), 4) if bim else 0,
                best_imag_rxn_ov=round(bim.get('rxn_overlap',0), 4) if bim else 0,
                best_imag_core_frac=round(bim.get('core_fraction',0), 4) if bim else 0,
            ))
            if ba > oracle_step_any:
                oracle_step_any = ba
                oracle_ig_label = ts['label']
                oracle_mode_idx = bm['idx'] if bm else -1
            if bia > oracle_step_imag:
                oracle_step_imag = bia

        per_step_rows.append(dict(
            step=data['step'],
            n_ig=len(igs),
            gt_freq=round(gt_freq, 4),
            oracle_any=round(oracle_step_any, 6),
            oracle_imag=round(oracle_step_imag, 6),
            oracle_ig_label=oracle_ig_label,
            oracle_mode_idx=oracle_mode_idx,
        ))

    # Build top-k sheet: per step, rank IGs by best_align_any desc and
    # report the k-th best alignment. Both "any-mode" and "imag-only".
    K = 20
    by_step = {}
    for r in per_ig_rows:
        by_step.setdefault(r['step'], []).append(r)
    topk_rows = []
    for step in sorted(by_step.keys()):
        ig_rows = by_step[step]
        any_sorted  = sorted([r['best_align_any']  for r in ig_rows], reverse=True)
        imag_sorted = sorted([r['best_align_imag'] for r in ig_rows], reverse=True)
        row = {'step': step, 'n_ig': len(ig_rows)}
        for k in range(K):
            row[f'top{k+1}_any']  = round(any_sorted[k], 6)  if k < len(any_sorted)  else 0.0
            row[f'top{k+1}_imag'] = round(imag_sorted[k], 6) if k < len(imag_sorted) else 0.0
        topk_rows.append(row)

    OUT_PER_IG.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PER_IG.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(per_ig_rows[0].keys()))
        w.writeheader(); w.writerows(per_ig_rows)
    with OUT_PER_STEP.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(per_step_rows[0].keys()))
        w.writeheader(); w.writerows(per_step_rows)
    with OUT_TOPK.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(topk_rows[0].keys()))
        w.writeheader(); w.writerows(topk_rows)

    print(f"\nDone in {time.time()-t0:.1f}s")
    print(f"Per-IG rows : {len(per_ig_rows):>5}  ->  {OUT_PER_IG}")
    print(f"Per-step    : {len(per_step_rows):>5}  ->  {OUT_PER_STEP}")
    print(f"Top-k sheet : {len(topk_rows):>5}  ->  {OUT_TOPK}")
    a = np.array([r['oracle_any'] for r in per_step_rows])
    ai = np.array([r['oracle_imag'] for r in per_step_rows])
    print(f"\nPer-step oracle (any modes):  mean={a.mean():.3f}  ≥0.7={(a>=0.7).mean()*100:.1f}%  ≥0.5={(a>=0.5).mean()*100:.1f}%  ≥0.3={(a>=0.3).mean()*100:.1f}%")
    print(f"Per-step oracle (imag only):  mean={ai.mean():.3f}  ≥0.7={(ai>=0.7).mean()*100:.1f}%  ≥0.5={(ai>=0.5).mean()*100:.1f}%  ≥0.3={(ai>=0.3).mean()*100:.1f}%")
    # Print top-k headline
    print(f"\nTop-k oracle (any) — % steps reaching threshold by k-th best IG:")
    print(f"{'k':>3}  {'mean':>5}  {'≥0.7':>5}  {'≥0.5':>5}  {'≥0.3':>5}")
    for k in (1, 2, 3, 5, 10, 20):
        col = np.array([r[f'top{k}_any'] for r in topk_rows])
        print(f"{k:>3}  {col.mean():.3f}  {(col>=0.7).mean()*100:4.1f}%  "
              f"{(col>=0.5).mean()*100:4.1f}%  {(col>=0.3).mean()*100:4.1f}%")


if __name__ == '__main__':
    main()

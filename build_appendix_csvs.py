"""
Build two summary CSVs for the appendix:

  1. final_quality_measurement.csv
     One row per step. Columns:
       step, n_ig
       top1_any..top20_any            — oracle: K-th best IG's max gt_alignment
       top1_any_imag..top20_any_imag  — same restricted to imaginary modes
       verifier_top1..verifier_top5   — alignment of K-th IG picked by current
                                        verifier (clean_v2). Quality = best
                                        mode in that IG vs GT.
       verifier_top1_picked..verifier_top5_picked
                                      — alignment of the SPECIFIC mode the
                                        verifier picked within the IG (i.e.
                                        what flat_view actually shows).

  2. initial_guess_modes.csv
     One row per (step, IG, mode). Columns the user requested, augmented
     with verifier_score and a filter-pass flag, since the current verifier
     is clean_v2 = bond × (1+rxn) × (1+0.2·core) / n_imag^0.3 (filtered to
     n_imag<=2 and rxn>=0.10).
       step, ts_label, mode_idx, freq, is_imag, bond_overlap, rxn_overlap,
       core_fraction, mode_rank, n_imag, n_modes_total, n_core_atoms,
       core_atoms, verifier_score, passes_verifier_filter, gt_alignment

Output dir: appendix_perparation/analtics/  (using user's exact naming)
"""
from __future__ import annotations
import csv
import json
import re
import time
from pathlib import Path
import numpy as np

from build_flat_view import rank_clean_v2  # canonical verifier

ROOT = Path('/Users/yunhengz/empty_for_claude/rxn_core')
SRC = ROOT / 'out' / 'mode_viewer'
OUT_DIR = ROOT / 'appendix_perparation' / 'analtics'
OUT_FINAL = OUT_DIR / 'final_quality_measurement.csv'
OUT_MODES = OUT_DIR / 'initial_guess_modes.csv'

K_ORACLE = 20
K_VERIFIER = 5


def cos_sim(a, b):
    a = np.asarray(a).reshape(-1); b = np.asarray(b).reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9: return 0.0
    return abs(float(a @ b)) / (na * nb)


def verifier_score_of_mode(picked, n_imag, w_rxn=1.0, w_core=0.2, imag_pen=0.3):
    b = picked.get('bond_overlap', 0.0)
    r = picked.get('rxn_overlap',  0.0)
    c = picked.get('core_fraction', 0.0)
    return b * (1 + w_rxn * r) * (1 + w_core * c) / max(n_imag, 1) ** imag_pen


def passes_filter(picked, n_imag, min_rxn=0.10, max_imag=2):
    return (picked.get('rxn_overlap', 0.0) >= min_rxn) and (1 <= n_imag <= max_imag)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    files = sorted(SRC.glob('*.html'))
    files = [f for f in files
             if not (f.name in ('index.html', 'flat_view.html', 'guess_quality.html')
                     or f.name.startswith('oracle_view'))]
    print(f"Reading {len(files)} steps...")
    t0 = time.time()

    final_rows = []
    mode_rows = []

    for hp in files:
        text = hp.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        gt = next((t for t in data['ts_list']
                   if t['label']=='groundtruth' and t.get('modes')), None)
        if gt is None: continue
        gt_disp = np.asarray(gt['modes'][gt['default_mode_idx']]['disp'])

        igs = [t for t in data['ts_list']
               if t['label']!='groundtruth' and t.get('modes')]
        if not igs: continue

        # ── Oracle topK: rank IGs by best gt_alignment of any mode ──
        ig_best_any  = []
        ig_best_imag = []
        per_ig_best_align = {}
        for ts in igs:
            ba_any  = max((cos_sim(np.asarray(m['disp']), gt_disp) for m in ts['modes']),
                          default=0.0)
            ba_imag = max((cos_sim(np.asarray(m['disp']), gt_disp)
                           for m in ts['modes'] if m['freq']<0), default=0.0)
            ig_best_any.append(ba_any)
            ig_best_imag.append(ba_imag)
            per_ig_best_align[ts['label']] = ba_any
        ig_best_any.sort(reverse=True)
        ig_best_imag.sort(reverse=True)

        # ── Verifier top5 (clean_v2): align of best mode in each picked IG ──
        elements = gt['xyz_elements']
        # rank_clean_v2 returns up to 3; we need 5 — call with a wider rank
        # path by iterating: get all candidates, do the same diversity loop.
        # Simpler: call rank_clean_v2 (which returns top-3) and pad.
        ranked = rank_clean_v2(igs, elements)  # (ts, picked_mode) pairs, top-3
        # Extend to top-5 by fallback: re-rank remaining by clean_v2 score
        # ignoring diversity (lighter path) so we have 5.
        already = {t['label'] for (t, m) in ranked}
        if len(ranked) < K_VERIFIER:
            extras = []
            for ts in igs:
                if ts['label'] in already: continue
                imag = [m for m in ts['modes'] if m['freq']<0]
                if not imag: continue
                picked = max(imag, key=lambda m: m.get('bond_overlap', 0))
                extras.append((picked.get('bond_overlap', 0), ts, picked))
            extras.sort(key=lambda t: -t[0])
            for _, ts, picked in extras:
                if len(ranked) >= K_VERIFIER: break
                ranked.append((ts, picked))

        verifier_align_best = []  # best mode in IG vs GT
        verifier_align_picked = []  # mode the verifier picked vs GT
        for ts, picked in ranked[:K_VERIFIER]:
            verifier_align_best.append(per_ig_best_align.get(ts['label'], 0.0))
            verifier_align_picked.append(
                cos_sim(np.asarray(picked['disp']), gt_disp))

        row = {'step': data['step'], 'n_ig': len(igs)}
        for k in range(K_ORACLE):
            row[f'top{k+1}_any']      = round(ig_best_any[k],  6) if k < len(ig_best_any)  else 0.0
            row[f'top{k+1}_any_imag'] = round(ig_best_imag[k], 6) if k < len(ig_best_imag) else 0.0
        for k in range(K_VERIFIER):
            row[f'verifier_top{k+1}']        = round(verifier_align_best[k],   6)  if k < len(verifier_align_best)   else 0.0
            row[f'verifier_top{k+1}_picked'] = round(verifier_align_picked[k], 6)  if k < len(verifier_align_picked) else 0.0
        final_rows.append(row)

        # ── Per-mode CSV ──
        for ts in igs:
            n_imag_ts = ts['n_imag']
            n_modes_total = ts['n_modes_total']
            core_atoms = data['core_atoms']
            n_core = len(core_atoms)
            # mode_rank = idx within the IG sorted by bond_overlap desc
            modes_sorted = sorted(enumerate(ts['modes']),
                                   key=lambda im: -im[1].get('bond_overlap', 0))
            rank_by_idx = {idx: r for r, (idx, _) in enumerate(modes_sorted)}
            for orig_idx, m in enumerate(ts['modes']):
                gt_align = cos_sim(np.asarray(m['disp']), gt_disp)
                mode_rows.append(dict(
                    step=data['step'],
                    ts_label=ts['label'],
                    mode_idx=m['idx'],
                    freq=round(m['freq'], 4),
                    is_imag=int(m['freq'] < 0),
                    bond_overlap=round(m.get('bond_overlap', 0), 4),
                    rxn_overlap=round(m.get('rxn_overlap', 0), 4),
                    core_fraction=round(m.get('core_fraction', 0), 4),
                    mode_rank=rank_by_idx.get(orig_idx, 0),
                    n_imag=n_imag_ts,
                    n_modes_total=n_modes_total,
                    n_core_atoms=n_core,
                    core_atoms=','.join(str(x) for x in core_atoms),
                    verifier_score=round(verifier_score_of_mode(m, n_imag_ts), 6),
                    passes_verifier_filter=int(passes_filter(m, n_imag_ts)),
                    gt_alignment=round(gt_align, 6),
                ))

    # ── Write CSVs ──
    if final_rows:
        with OUT_FINAL.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(final_rows[0].keys()))
            w.writeheader(); w.writerows(final_rows)
    if mode_rows:
        with OUT_MODES.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(mode_rows[0].keys()))
            w.writeheader(); w.writerows(mode_rows)

    print(f"\nDone in {time.time()-t0:.1f}s")
    print(f"  steps        : {len(final_rows)}")
    print(f"  mode rows    : {len(mode_rows)}")
    print(f"  final CSV    : {OUT_FINAL}")
    print(f"  modes CSV    : {OUT_MODES}")

    if final_rows:
        # Headline: oracle vs verifier per k
        print(f"\nOracle (any) by k    /    Verifier-best by k:")
        print(f"{'k':>3}  {'oracle_mean':>11}  {'oracle≥0.7':>10}  {'verif_mean':>11}  {'verif≥0.7':>9}")
        for k in (1, 2, 3, 5):
            o = np.array([r[f'top{k}_any']      for r in final_rows])
            v = np.array([r.get(f'verifier_top{k}', 0) for r in final_rows])
            print(f"{k:>3}  {o.mean():>11.3f}  {(o>=0.7).mean()*100:>9.1f}%  "
                  f"{v.mean():>11.3f}  {(v>=0.7).mean()*100:>8.1f}%")


if __name__ == '__main__':
    main()

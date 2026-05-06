"""
Compare alternative per-IG ranking strategies against the GT-alignment
oracle. Pure data extraction from per-step HTMLs — no recompute.

Per IG TS, every ranker R picks one mode (its idea of "best") under
its criterion. Per step, R picks one IG (the IG whose R-picked mode
scores highest by R). We measure that IG's gt_alignment, where
gt_alignment = |cos similarity| of the R-picked mode against the GT's
bond_overlap-picked default mode (the gold reference for the step).

Oracle = upper bound: the best gt_alignment achievable across all
modes of all IGs in the step.

Output:
  out/mode_analysis/ranker_comparison.csv  — per (step, ranker) row
  prints summary table to stdout
"""
from __future__ import annotations
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


SRC_DIR = Path(__file__).parent / "out" / "mode_viewer"
OUT_CSV = Path(__file__).parent / "out" / "mode_analysis" / "ranker_comparison.csv"


def cos_sim(a, b):
    a = np.asarray(a).reshape(-1)
    b = np.asarray(b).reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9:
        return 0.0
    return abs(float(a @ b)) / (na * nb)


# Each ranker = (name, score_fn). score_fn(mode_dict) -> float; modes
# of the same TS are ranked by score_fn; the highest score wins. Imag
# modes are preferred (real modes get -inf score for these rankers
# unless explicitly chosen). For some rankers we override "imag-only".
def imag_only(score_fn):
    def f(m):
        return score_fn(m) if m['freq'] < 0 else -1e18
    return f


RANKERS = [
    ('bond_overlap',          imag_only(lambda m: m.get('bond_overlap', 0))),
    ('rxn_overlap',           imag_only(lambda m: m.get('rxn_overlap', 0))),
    ('core_fraction',         imag_only(lambda m: m.get('core_fraction', 0))),
    ('most_negative_freq',    imag_only(lambda m: -m['freq'])),  # higher score = more neg
    ('bond_x_core',           imag_only(lambda m: m.get('bond_overlap', 0) * m.get('core_fraction', 0))),
    ('bond_x_rxn',            imag_only(lambda m: m.get('bond_overlap', 0) * m.get('rxn_overlap', 0))),
    ('bond_plus_rxn',         imag_only(lambda m: m.get('bond_overlap', 0) + m.get('rxn_overlap', 0))),
    ('bond_x_absfreq',        imag_only(lambda m: m.get('bond_overlap', 0) * abs(m['freq']) / 1000.0)),
    # Rank-fusion: sum of ranks across bond_ov + rxn_ov + core_frac. Lower is
    # better, so we negate so higher score wins when picking max.
    # Implemented below per-TS.
]


def rank_fusion_score(modes_imag):
    """Sum-of-ranks across bond_ov + rxn_ov + core_frac for imag modes.
    Returns dict mode_idx_in_imag_list -> score (higher = better)."""
    n = len(modes_imag)
    if n == 0: return {}
    out = {i: 0 for i in range(n)}
    for key in ('bond_overlap', 'rxn_overlap', 'core_fraction'):
        order = sorted(range(n), key=lambda i: -modes_imag[i].get(key, 0))
        for rank, idx in enumerate(order):
            out[idx] += rank  # lower rank = higher score; we'll invert
    # Invert: highest sum-of-ranks → worst → lowest score
    return {i: -v for i, v in out.items()}


def evaluate_step(payload):
    """Returns dict: ranker_name → gt_alignment of ranker top-1 IG.
    Plus 'oracle' = best achievable gt_alignment for the step."""
    ts_by_label = {ts['label']: ts for ts in payload['ts_list']}
    gt = ts_by_label.get('groundtruth')
    if gt is None or not gt.get('modes'):
        return None
    gt_default = gt['modes'][gt.get('default_mode_idx', 0)]
    if gt_default['freq'] >= 0:
        return None
    gt_disp = gt_default['disp']

    # Pre-compute cos_sim of every mode in every IG vs gt_disp.
    # Store on the mode dict in-memory (cached).
    ig_tss = [ts for ts in payload['ts_list']
              if ts['label'] != 'groundtruth' and ts.get('modes')]
    if not ig_tss:
        return None
    for ts in ig_tss:
        for m in ts['modes']:
            m['_align'] = cos_sim(m['disp'], gt_disp)

    out = {'step': payload['step']}

    # Run each scoring ranker
    for name, score_fn in RANKERS:
        # Per-IG: pick mode with highest score_fn (only imag → score = -inf else)
        ig_picks = []
        for ts in ig_tss:
            scored = [(score_fn(m), m, ts['label']) for m in ts['modes']]
            scored.sort(key=lambda t: -t[0])
            top_score, top_mode, ts_label = scored[0]
            if top_score <= -1e17:
                # No imag mode in this TS — fall back to overall best by score_fn.replacement
                # Use a signal that doesn't penalise real modes
                fallback = max(ts['modes'], key=lambda m: score_fn.__wrapped__(m) if hasattr(score_fn, '__wrapped__') else -1e18)
                top_mode = fallback
            ig_picks.append((top_score, top_mode, ts_label))
        # Per-step: pick IG with highest top_score
        ig_picks.sort(key=lambda t: -t[0])
        chosen_score, chosen_mode, chosen_label = ig_picks[0]
        out[name] = chosen_mode['_align']
        out[name + '_label'] = chosen_label

    # Rank fusion (separate handling)
    ranker_picks = []
    for ts in ig_tss:
        imag_modes = [m for m in ts['modes'] if m['freq'] < 0]
        if not imag_modes:
            best_mode = max(ts['modes'], key=lambda m: m.get('bond_overlap', 0))
            best_score = -1e18
        else:
            scores = rank_fusion_score(imag_modes)
            best_i = max(scores.keys(), key=lambda i: scores[i])
            best_mode = imag_modes[best_i]
            best_score = scores[best_i]
        ranker_picks.append((best_score, best_mode, ts['label']))
    ranker_picks.sort(key=lambda t: -t[0])
    out['rank_fusion'] = ranker_picks[0][1]['_align']
    out['rank_fusion_label'] = ranker_picks[0][2]

    # Oracle: best _align across all modes of all IGs (any freq)
    oracle_align = 0.0
    oracle_label = '—'
    for ts in ig_tss:
        for m in ts['modes']:
            if m['_align'] > oracle_align:
                oracle_align = m['_align']
                oracle_label = ts['label']
    out['oracle'] = oracle_align
    out['oracle_label'] = oracle_label

    # Imag-restricted oracle (best _align among imag modes only)
    oracle_imag = 0.0
    oracle_imag_label = '—'
    for ts in ig_tss:
        for m in ts['modes']:
            if m['freq'] < 0 and m['_align'] > oracle_imag:
                oracle_imag = m['_align']
                oracle_imag_label = ts['label']
    out['oracle_imag'] = oracle_imag
    out['oracle_imag_label'] = oracle_imag_label

    return out


def main():
    files = sorted(SRC_DIR.glob("*.html"))
    files = [f for f in files
             if f.name not in ('index.html', 'flat_view.html', 'guess_quality.html')]

    rows = []
    for hp in files:
        text = hp.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        r = evaluate_step(data)
        if r is not None:
            rows.append(r)

    if not rows:
        print("no rows produced")
        return

    fieldnames = list(rows[0].keys())
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)

    # Summary
    rankers = [n for n, _ in RANKERS] + ['rank_fusion']
    print(f"Evaluated {len(rows)} steps.\n")
    print(f"{'ranker':25s}  {'mean':>7s}  {'median':>7s}  "
          f"{'≥0.7':>5s}  {'≥0.5':>5s}  {'gap_to_oracle':>14s}")
    oracle_vals = np.array([r['oracle'] for r in rows])
    oracle_imag_vals = np.array([r['oracle_imag'] for r in rows])
    print('-' * 72)
    for name in rankers:
        v = np.array([r[name] for r in rows])
        gap = oracle_vals.mean() - v.mean()
        print(f"{name:25s}  {v.mean():7.3f}  {np.median(v):7.3f}  "
              f"{(v >= 0.7).mean()*100:4.0f}%  {(v >= 0.5).mean()*100:4.0f}%  "
              f"{gap:14.3f}")
    print('-' * 72)
    print(f"{'oracle_imag (any imag)':25s}  {oracle_imag_vals.mean():7.3f}  "
          f"{np.median(oracle_imag_vals):7.3f}  "
          f"{(oracle_imag_vals >= 0.7).mean()*100:4.0f}%  "
          f"{(oracle_imag_vals >= 0.5).mean()*100:4.0f}%  "
          f"{oracle_vals.mean()-oracle_imag_vals.mean():14.3f}")
    print(f"{'oracle (any mode)':25s}  {oracle_vals.mean():7.3f}  "
          f"{np.median(oracle_vals):7.3f}  "
          f"{(oracle_vals >= 0.7).mean()*100:4.0f}%  "
          f"{(oracle_vals >= 0.5).mean()*100:4.0f}%  "
          f"{0.0:14.3f}")
    print(f"\nCSV: {OUT_CSV}")


if __name__ == "__main__":
    main()

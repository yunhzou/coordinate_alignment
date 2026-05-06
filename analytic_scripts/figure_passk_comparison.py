"""
Figure: pass@k comparison — verifier (clean_v2) vs uniform-random
baseline vs oracle ceiling, k = 1..5.

Reads appendix_perparation/analtics/final_quality_measurement.csv,
which contains for each step:

  oracle_pass<k>    — ceiling, constant in k (= best alignment over
                      the entire IG pool)
  verifier_pass<k>  — cumulative max alignment over verifier's first
                      k picks
  random_pass<k>    — closed-form E[max alignment over k IGs sampled
                      uniformly without replacement from the pool]

Output:
  appendix_perparation/figures/passk_comparison.{png,pdf}
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


CSV_PATH = PROJECT_ROOT / 'appendix_perparation' / 'analtics' / 'final_quality_measurement.csv'
OUT_DIR  = PROJECT_ROOT / 'appendix_perparation' / 'figures'

K_MAX = 5


def main():
    rows = list(csv.DictReader(CSV_PATH.open()))
    n = len(rows)
    print(f"Loaded {n} step rows from {CSV_PATH}")

    # Per-method per-k arrays
    methods = {
        'oracle':   ('oracle_pass{}',   '#7a7a7a',  '--', 'Oracle ceiling'),
        'verifier': ('verifier_pass{}', '#cc3333',  '-',  'clean_v2 verifier'),
        'random':   ('random_pass{}',   '#3355aa',  '-.', 'Uniform-random baseline'),
    }
    ks = list(range(1, K_MAX + 1))
    arrs = {}  # method → 2D array (n_steps × K)
    for name, (col_fmt, _, _, _) in methods.items():
        arr = np.zeros((n, K_MAX))
        for i, r in enumerate(rows):
            for k in ks:
                arr[i, k - 1] = float(r[col_fmt.format(k)])
        arrs[name] = arr

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 2x2 grid: mean alignment + ≥0.7 + ≥0.5 + ≥0.3 pass-rate
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    axes = axes.flatten()
    panels = [
        ('mean alignment',                  lambda a: a.mean(axis=0),                    None),
        ('% steps reaching alignment ≥ 0.7', lambda a: (a >= 0.7).mean(axis=0) * 100,    '%'),
        ('% steps reaching alignment ≥ 0.5', lambda a: (a >= 0.5).mean(axis=0) * 100,    '%'),
        ('% steps reaching alignment ≥ 0.3', lambda a: (a >= 0.3).mean(axis=0) * 100,    '%'),
    ]
    # Vertical-offset (in points) per method so labels don't collide.
    label_offsets = {'oracle': 12, 'verifier': 12, 'random': -18}

    def _label_ks(method, k_list):
        """Where to draw value annotations.
        Oracle is constant in k → only label endpoints (k=1, k=K)."""
        if method == 'oracle':
            return [k_list[0], k_list[-1]]
        return k_list

    for ax, (title, fn, units) in zip(axes, panels):
        # Pre-compute per-method curves to detect collisions for label placement.
        ys = {name: fn(arrs[name]) for name in methods}
        ymax_data = max(y.max() for y in ys.values())
        for name, (_, color, ls, label) in methods.items():
            y = ys[name]
            ax.plot(ks, y, marker='o', linestyle=ls, color=color, label=label, lw=2.2,
                    markersize=7)
            for k in _label_ks(name, ks):
                val = y[k - 1]
                # If this point is within ±2 % of another curve at same k,
                # nudge the label vertically to avoid overlap.
                offset_y = label_offsets[name]
                for other in methods:
                    if other == name: continue
                    if abs(ys[other][k - 1] - val) < (0.02 * (1 if units is None else 100)):
                        # collision: push verifier down if oracle is the same value;
                        # push random further down
                        if name == 'verifier' and ys.get('oracle', y)[k - 1] - val < 1e-6:
                            offset_y = -18
                        if name == 'random':
                            offset_y = -22
                fmt = f'{val:.2f}' if units is None else f'{val:.0f}'
                ax.annotate(fmt, (k, val), textcoords='offset points',
                            xytext=(0, offset_y), ha='center', fontsize=9, color=color)
        ax.set_xticks(ks)
        ax.set_xlabel('k  (top-k IGs picked)')
        ax.set_title(title, fontsize=12, pad=14)
        ax.grid(alpha=0.3)
        ax.set_xlim(0.55, K_MAX + 0.45)
        # Headroom so the top labels don't crowd the panel border
        if units == '%':
            ax.set_ylabel('% of steps')
            ax.set_ylim(0, max(ymax_data * 1.22, 5))
        else:
            ax.set_ylabel('mean cosine alignment')
            ax.set_ylim(0, max(ymax_data * 1.22, 0.05))
    axes[0].legend(loc='lower right', fontsize=10, framealpha=0.95)
    fig.suptitle(f'pass@k: clean_v2 verifier vs uniform-random baseline   '
                 f'(oracle ceiling shown for reference, N = 155 steps)',
                 fontsize=13, y=0.995)
    fig.subplots_adjust(top=0.90, hspace=0.42, wspace=0.25,
                        left=0.08, right=0.97, bottom=0.07)
    # NB: do NOT call tight_layout() here — it would override the
    # subplots_adjust spacing we set above.

    out_png = OUT_DIR / 'passk_comparison.png'
    out_pdf = OUT_DIR / 'passk_comparison.pdf'
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nSaved:\n  {out_png}\n  {out_pdf}")

    # Print headline numbers
    print(f"\n{'k':>3}  {'metric':28s}  {'oracle':>8s}  {'verif':>8s}  {'random':>8s}  {'Δ(verif−rd)':>11s}")
    for title, fn, _ in panels:
        for k in ks:
            o, v, rd = fn(arrs['oracle'])[k-1], fn(arrs['verifier'])[k-1], fn(arrs['random'])[k-1]
            print(f"{k:>3}  {title:28s}  {o:8.3f}  {v:8.3f}  {rd:8.3f}  {(v-rd):+11.3f}")


if __name__ == '__main__':
    main()

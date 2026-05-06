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

    # Four panels: mean alignment + ≥0.7 + ≥0.5 + ≥0.3 pass-rate
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.2))
    panels = [
        ('mean alignment',                  lambda a: a.mean(axis=0),                    None),
        ('% steps reaching alignment ≥ 0.7', lambda a: (a >= 0.7).mean(axis=0) * 100,    '%'),
        ('% steps reaching alignment ≥ 0.5', lambda a: (a >= 0.5).mean(axis=0) * 100,    '%'),
        ('% steps reaching alignment ≥ 0.3', lambda a: (a >= 0.3).mean(axis=0) * 100,    '%'),
    ]
    for ax, (title, fn, units) in zip(axes, panels):
        for name, (_, color, ls, label) in methods.items():
            y = fn(arrs[name])
            ax.plot(ks, y, marker='o', linestyle=ls, color=color, label=label, lw=2,
                    markersize=6)
            for k, val in zip(ks, y):
                fmt = f'{val:.2f}' if units is None else f'{val:.0f}'
                ax.annotate(fmt, (k, val), textcoords='offset points',
                            xytext=(0, 8 if name != 'random' else -14),
                            ha='center', fontsize=8, color=color)
        ax.set_xticks(ks)
        ax.set_xlabel('k (top-k IGs picked)')
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
        ax.set_xlim(0.7, K_MAX + 0.3)
        if units == '%':
            ax.set_ylim(bottom=0)
            ax.set_ylabel('% of steps')
        else:
            ax.set_ylabel('mean cosine alignment')
            ax.set_ylim(bottom=0)
    axes[0].legend(loc='lower right', fontsize=9)
    fig.suptitle(f'pass@k comparison on 155 elementary steps\n'
                 f'(verifier vs uniform-random baseline; oracle ceiling shown for reference)',
                 fontsize=11)
    fig.tight_layout()

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

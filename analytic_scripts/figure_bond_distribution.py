"""
Figure: distribution of (broken + formed) bond counts across steps.

Reads each per-step viewer HTML in
  appendix_perparation/viewer/mode_viewer/<step>.html
extracts `broken_bonds` and `formed_bonds_R` from the embedded JSON,
computes n_broken + n_formed per step, and saves a histogram with
integer bins (binwidth = 1).

Output:
  appendix_perparation/figures/bond_count_distribution.png
  appendix_perparation/figures/bond_count_distribution.pdf
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


SRC = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
OUT_DIR = PROJECT_ROOT / 'appendix_perparation' / 'figures'


def main():
    files = sorted(p for p in SRC.glob('*.html')
                   if p.name not in ('flat_view.html', 'guess_quality.html', 'index.html'))
    print(f"Reading {len(files)} per-step HTMLs from {SRC}")

    rows = []
    for p in files:
        text = p.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        n_b = len(data.get('broken_bonds', []))
        n_f = len(data.get('formed_bonds_R', []))
        rows.append((data['step'], n_b, n_f, n_b + n_f))

    n_steps = len(rows)
    n_broken = np.array([r[1] for r in rows])
    n_formed = np.array([r[2] for r in rows])
    n_total  = np.array([r[3] for r in rows])

    print(f"\n{n_steps} steps loaded")
    print(f"  broken: mean={n_broken.mean():.2f}  median={np.median(n_broken):.0f}  max={n_broken.max()}")
    print(f"  formed: mean={n_formed.mean():.2f}  median={np.median(n_formed):.0f}  max={n_formed.max()}")
    print(f"  total : mean={n_total.mean():.2f}   median={np.median(n_total):.0f}   max={n_total.max()}")

    # Histogram with bin=1
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    bin_max = int(max(n_broken.max(), n_formed.max(), n_total.max())) + 1
    bins = np.arange(-0.5, bin_max + 0.5, 1)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, arr, title, color in [
        (axes[0], n_broken, 'Broken bonds',          '#cc3333'),
        (axes[1], n_formed, 'Formed bonds',          '#2a8a2a'),
        (axes[2], n_total,  'Broken + formed total', '#3355aa'),
    ]:
        counts, _, _ = ax.hist(arr, bins=bins, color=color, edgecolor='white',
                               linewidth=0.7, alpha=0.85)
        ax.set_xlabel('count per step')
        ax.set_title(f'{title}\n(mean={arr.mean():.2f}, median={np.median(arr):.0f}, max={arr.max()})',
                     fontsize=10)
        ax.set_xticks(range(0, bin_max + 1))
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        # Annotate counts above each bar
        for b, c in zip(bins[:-1] + 0.5, counts):
            if c > 0:
                ax.text(b, c + 0.5, f'{int(c)}', ha='center', va='bottom', fontsize=8)
    axes[0].set_ylabel(f'number of steps (out of {n_steps})')
    fig.suptitle(f'Bond-count distribution per elementary step (N={n_steps}, binwidth=1)',
                 fontsize=12)
    fig.tight_layout()

    out_png = OUT_DIR / 'bond_count_distribution.png'
    out_pdf = OUT_DIR / 'bond_count_distribution.pdf'
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nSaved:\n  {out_png}\n  {out_pdf}")

    # Also print the cumulative breakdown for the report
    print(f"\n(broken, formed) → step count:")
    pair_counts = Counter((r[1], r[2]) for r in rows)
    for (b, f), n in sorted(pair_counts.items()):
        print(f"  ({b}, {f}): {n}")


if __name__ == '__main__':
    main()

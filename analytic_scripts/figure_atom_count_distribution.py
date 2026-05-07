"""
Figure: distribution of atom counts (n_atoms) across the 155-step
benchmark set.

Reads per-step viewer HTMLs in appendix_perparation/viewer/mode_viewer/
extracts `n_atoms`, plots a histogram with 10-atom bins.

Output:
  appendix_perparation/figures/atom_count_distribution.{png,pdf}
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

    n_atoms = []
    for p in files:
        text = p.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        n_atoms.append(data['n_atoms'])
    n_atoms = np.array(n_atoms)
    n = len(n_atoms)

    print(f"\n{n} steps")
    print(f"  min    = {n_atoms.min()}")
    print(f"  max    = {n_atoms.max()}")
    print(f"  mean   = {n_atoms.mean():.1f}")
    print(f"  median = {int(np.median(n_atoms))}")
    print(f"  IQR    = [{int(np.percentile(n_atoms,25))}, {int(np.percentile(n_atoms,75))}]")

    # 10-atom bins from 10 to 160
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    binwidth = 10
    bins = np.arange(10, n_atoms.max() + binwidth + 1, binwidth)
    fig, ax = plt.subplots(figsize=(10, 5))
    counts, _, patches = ax.hist(n_atoms, bins=bins, color='#3a6dbf',
                                  edgecolor='white', linewidth=0.8, alpha=0.9)
    # Annotate counts
    for b, c in zip(bins[:-1], counts):
        if c > 0:
            ax.text(b + binwidth/2, c + 0.3, f'{int(c)}',
                    ha='center', va='bottom', fontsize=10, color='#1f3d77')
    ax.set_xlabel('atom count per step', fontsize=11)
    ax.set_ylabel(f'number of steps (out of {n})', fontsize=11)
    ax.set_title(f'Atom-count distribution across the {n}-step benchmark   '
                 f'(binwidth = {binwidth} atoms)\n'
                 f'min = {n_atoms.min()},  median = {int(np.median(n_atoms))},  '
                 f'mean = {n_atoms.mean():.1f},  max = {n_atoms.max()}',
                 fontsize=11, pad=12)
    ax.set_xticks(bins)
    ax.set_xlim(bins[0] - 1, bins[-1] + 1)
    ax.grid(axis='y', linestyle=':', alpha=0.4)
    ax.set_axisbelow(True)
    # Mark mean and median
    ax.axvline(n_atoms.mean(),       color='#cc3333', linestyle='--', lw=1.4,
               label=f'mean = {n_atoms.mean():.1f}')
    ax.axvline(np.median(n_atoms),   color='#2a8a2a', linestyle='-.', lw=1.4,
               label=f'median = {int(np.median(n_atoms))}')
    ax.legend(loc='upper right', fontsize=10, framealpha=0.95)
    fig.tight_layout()

    out_png = OUT_DIR / 'atom_count_distribution.png'
    out_pdf = OUT_DIR / 'atom_count_distribution.pdf'
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nSaved:\n  {out_png}\n  {out_pdf}")


if __name__ == '__main__':
    main()

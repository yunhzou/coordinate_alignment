"""
Figure: human-judged verifier accuracy as a function of step size
(atom count). For each 10-atom bin we compute the pass rate
   accuracy = #(pass=1) / #(steps in bin)
where pass=1 iff at least one of the verifier's top-2 picks was
human-judged Good (column `pass` in the human-rated CSV).

Output:
  appendix_perparation/figures/accuracy_by_atom_count.{png,pdf}
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
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


VIEWER_DIR = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
HUMAN_CSV  = PROJECT_ROOT / 'appendix_perparation' / 'analtics' / 'final_quality_measurement-humanversion (1).csv'
OUT_DIR    = PROJECT_ROOT / 'appendix_perparation' / 'figures'

BINWIDTH = 10


def canonical(s):
    return re.sub(r'\s*\(.*?\)\s*$', '', s).strip()


def main():
    # Per-step n_atoms from per-step HTMLs
    n_atoms = {}
    for p in VIEWER_DIR.glob('*.html'):
        if p.name in ('flat_view.html', 'guess_quality.html', 'index.html'):
            continue
        text = p.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        n_atoms[data['step']] = data['n_atoms']
    print(f"Atom counts loaded for {len(n_atoms)} steps")

    # Per-step human pass labels
    lines = HUMAN_CSV.read_text().splitlines()
    header = lines[1].split(',')
    rows = [dict(zip(header, ln.split(','))) for ln in lines[2:] if ln.strip()]
    pass_label = {}
    for r in rows:
        step = canonical(r['step'])
        try: pass_label[step] = int(r['pass'])
        except (ValueError, KeyError): pass

    print(f"Pass labels loaded for {len(pass_label)} steps")

    # Join + bin
    joined = []
    for step, n in n_atoms.items():
        if step in pass_label:
            joined.append((step, n, pass_label[step]))
    print(f"Joined: {len(joined)} steps with both atom count and pass label")

    n_arr = np.array([j[1] for j in joined])
    p_arr = np.array([j[2] for j in joined])
    print(f"Overall pass rate: {p_arr.mean()*100:.1f}%   ({p_arr.sum()}/{len(p_arr)})")

    # Bins of width 10 from 10 up to 90, then a single open-ended "90+" bin.
    fixed_edges = np.arange(10, 90 + 1, BINWIDTH)
    bin_total = np.zeros(len(fixed_edges), dtype=int)  # last index = 90+
    bin_pass  = np.zeros(len(fixed_edges), dtype=int)
    bin_labels = [f"[{int(fixed_edges[i])}, {int(fixed_edges[i]+BINWIDTH)})"
                  for i in range(len(fixed_edges) - 1)] + ['90+']
    bin_centers = list(fixed_edges[:-1] + BINWIDTH/2) + [95]  # 90+ centered at 95
    for n, p in zip(n_arr, p_arr):
        if n >= 90:
            idx = len(fixed_edges) - 1
        else:
            idx = int((n - 10) // BINWIDTH)
            idx = max(0, min(idx, len(fixed_edges) - 2))
        bin_total[idx] += 1
        bin_pass[idx]  += int(p)
    bin_centers = np.array(bin_centers)

    # 95% Wilson confidence interval per bin (small-N appropriate)
    def wilson(k, n, z=1.96):
        if n == 0: return (0.0, 0.0, 0.0)
        p = k / n
        denom = 1 + z*z/n
        center = (p + z*z/(2*n)) / denom
        half = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
        return p, max(0.0, center - half), min(1.0, center + half)

    rates = []; lo = []; hi = []
    for k, n in zip(bin_pass, bin_total):
        p, l, h = wilson(int(k), int(n))
        rates.append(p); lo.append(l); hi.append(h)
    rates = np.array(rates); lo = np.array(lo); hi = np.array(hi)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, (ax_acc, ax_n) = plt.subplots(2, 1, figsize=(11, 7),
                                         sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    bar_x = bin_centers
    bar_h = rates * 100
    nonzero = bin_total > 0
    ax_acc.bar(bar_x[nonzero], bar_h[nonzero], width=BINWIDTH * 0.86,
                color='#3a6dbf', edgecolor='white', linewidth=0.8, alpha=0.9,
                label='accuracy = pass rate per bin')
    err_lo = (rates[nonzero] - lo[nonzero]) * 100
    err_hi = (hi[nonzero] - rates[nonzero]) * 100
    ax_acc.errorbar(bar_x[nonzero], bar_h[nonzero], yerr=[err_lo, err_hi],
                     fmt='none', ecolor='black', capsize=4, lw=1.2, alpha=0.7)
    for x, k, n_, r in zip(bar_x, bin_pass, bin_total, rates):
        if n_ > 0:
            ax_acc.text(x, r*100 + 2, f'{k}/{n_}', ha='center', va='bottom', fontsize=9)
    overall = p_arr.mean() * 100
    ax_acc.axhline(overall, color='#cc3333', linestyle='--', lw=1.4,
                    label=f'overall = {overall:.1f}%')
    ax_acc.set_ylabel('verifier accuracy (% steps with ≥ 1 good IG)', fontsize=11)
    ax_acc.set_ylim(0, 110)
    ax_acc.set_yticks(range(0, 101, 20))
    ax_acc.set_xticks(bar_x)
    ax_acc.set_xticklabels(bin_labels, fontsize=10)
    ax_acc.grid(axis='y', linestyle=':', alpha=0.4)
    ax_acc.set_axisbelow(True)
    ax_acc.legend(loc='lower left', fontsize=9, framealpha=0.95)
    ax_acc.set_title(f'Verifier accuracy vs step size   '
                     f'(N = {len(joined)}, binwidth = {BINWIDTH} atoms — last bin pools 90+ atoms; '
                     f'error bars = 95% Wilson CI)',
                     fontsize=12, pad=12)

    # Bottom: bin populations
    ax_n.bar(bar_x[nonzero], bin_total[nonzero], width=BINWIDTH * 0.86,
              color='#888', edgecolor='white', linewidth=0.6, alpha=0.85)
    for x, n_ in zip(bar_x, bin_total):
        if n_ > 0:
            ax_n.text(x, n_ + 0.4, f'{n_}', ha='center', va='bottom', fontsize=9, color='#444')
    ax_n.set_xlabel('atom count per step', fontsize=11)
    ax_n.set_ylabel('# steps in bin', fontsize=10)
    ax_n.set_xticks(bar_x)
    ax_n.set_xticklabels(bin_labels, fontsize=10)
    ax_n.set_xlim(bar_x[0] - BINWIDTH/2 - 1, bar_x[-1] + BINWIDTH/2 + 1)
    ax_n.grid(axis='y', linestyle=':', alpha=0.4)
    ax_n.set_axisbelow(True)
    fig.tight_layout()

    out_png = OUT_DIR / 'accuracy_by_atom_count.png'
    out_pdf = OUT_DIR / 'accuracy_by_atom_count.pdf'
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nSaved:\n  {out_png}\n  {out_pdf}")

    # Print the table
    print(f"\n{'bin':>10s}   {'pass':>4s} / {'total':>5s}    {'rate':>5s}    {'95% CI':>16s}")
    for i, lbl in enumerate(bin_labels):
        if bin_total[i] == 0: continue
        print(f"  {lbl:>10s}    {int(bin_pass[i]):>4} / {int(bin_total[i]):>5}    "
              f"{rates[i]*100:>5.1f}%    [{lo[i]*100:>5.1f}, {hi[i]*100:>5.1f}]")


if __name__ == '__main__':
    main()

"""
Figure: pass@1 and pass@2 accuracy by atom-count bin, with the
deterministic ReactOT baseline.

Three bars per bin (same style as accuracy_by_atom_count_with_reactot.py):
  - ours pass@1  : IG#1 column == 1
  - ours pass@2  : pass column == 1  (= IG#1 OR IG#2)
  - ReactOT      : single deterministic output per step; no pass@2

Bins: 10-atom width, last bin pools 90+. Top panel shows percentages
with k/N annotations and overall dashed lines. Bottom panel shows bin
populations and the union of elements per bin (organic first, then
alphabetical, wrapping at 6).

Output:
  appendix_perparation/figures/pass1_pass2_with_reactot.{png,pdf}
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


VIEWER_DIR = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
HUMAN_CSV  = PROJECT_ROOT / 'appendix_perparation' / 'analtics' / 'final_quality_measurement-humanversion (1).csv'
OUT_DIR    = PROJECT_ROOT / 'appendix_perparation' / 'figures'

BINWIDTH = 10

REACT_OT_PASS = {
    "Jackie_TS_19",
    "pr16.carbocation_ts2", "pr16.carbocation_ts3", "pr16.carbocation_ts4",
    "pr16.carbocation_ts5", "pr16.carbocation_ts7",
    "pr16.carbocation_ts11", "pr16.carbocation_ts12",
}


def canonical(s):
    return re.sub(r'\s*\(.*?\)\s*$', '', s).strip()


def main():
    # Per-step n_atoms + element set
    n_atoms = {}; elements_per_step = {}
    for p in VIEWER_DIR.glob('*.html'):
        if p.name in ('flat_view.html', 'guess_quality.html', 'index.html'): continue
        m = re.search(r"const DATA = (\{.*?\});\n", p.read_text(), re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        n_atoms[data['step']] = data['n_atoms']
        ts_list = data.get('ts_list', [])
        elements_per_step[data['step']] = (
            set(ts_list[0]['xyz_elements']) if ts_list and ts_list[0].get('xyz_elements') else set()
        )

    # Human IG#1 / IG#2 / pass labels
    lines = HUMAN_CSV.read_text().splitlines()
    header = lines[1].split(',')
    rows = [dict(zip(header, ln.split(','))) for ln in lines[2:] if ln.strip()]
    p1 = {}; p2 = {}
    for r in rows:
        s = canonical(r['step'])
        try: p1[s] = int(r['IG#1'])
        except (ValueError, KeyError): pass
        try: p2[s] = int(r['pass'])
        except (ValueError, KeyError): pass

    # Join
    joined = []
    for step, n in n_atoms.items():
        if step in p1 and step in p2:
            joined.append((step, n, p1[step], p2[step], 1 if step in REACT_OT_PASS else 0))

    n_arr  = np.array([j[1] for j in joined])
    p1_arr = np.array([j[2] for j in joined])
    p2_arr = np.array([j[3] for j in joined])
    ot_arr = np.array([j[4] for j in joined])
    print(f"{len(joined)} steps  |  pass@1={p1_arr.sum()}/{len(joined)}  "
          f"pass@2={p2_arr.sum()}/{len(joined)}  react_OT={ot_arr.sum()}/{len(joined)}")

    # Bins: [10,20)..[80,90), 90+
    fixed_edges = np.arange(10, 90 + 1, BINWIDTH)
    n_bins = len(fixed_edges)
    centers = list(fixed_edges[:-1] + BINWIDTH / 2) + [95]
    centers = np.array(centers)
    labels  = [f"[{int(e)}, {int(e+BINWIDTH)})" for e in fixed_edges[:-1]] + ['90+']

    bin_total = np.zeros(n_bins, dtype=int)
    bin_p1    = np.zeros(n_bins, dtype=int)
    bin_p2    = np.zeros(n_bins, dtype=int)
    bin_ot    = np.zeros(n_bins, dtype=int)
    bin_elems = [set() for _ in range(n_bins)]
    for step, n, p1v, p2v, otv in joined:
        idx = (n_bins - 1) if n >= 90 else int((n - 10) // BINWIDTH)
        idx = max(0, min(idx, n_bins - 1))
        bin_total[idx] += 1
        bin_p1[idx]    += int(p1v)
        bin_p2[idx]    += int(p2v)
        bin_ot[idx]    += int(otv)
        bin_elems[idx] |= elements_per_step.get(step, set())

    rates_p1 = np.where(bin_total > 0, bin_p1 / np.maximum(bin_total, 1), 0)
    rates_p2 = np.where(bin_total > 0, bin_p2 / np.maximum(bin_total, 1), 0)
    rates_ot = np.where(bin_total > 0, bin_ot / np.maximum(bin_total, 1), 0)

    # Element-list formatter (organic first then alphabetical, wrap at 6)
    PRIORITY = ['C', 'H', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']
    def elem_lines(eset, max_per_line=6):
        prio = [e for e in PRIORITY if e in eset]
        rest = sorted(e for e in eset if e not in PRIORITY)
        elems = prio + rest
        if not elems: return ''
        out = ','.join(elems[:max_per_line])
        if len(elems) > max_per_line:
            out += '\n' + ','.join(elems[max_per_line:])
        return out

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, (ax_acc, ax_n) = plt.subplots(2, 1, figsize=(15, 9),
                                          sharex=True,
                                          gridspec_kw={'height_ratios': [2.6, 1.4]})
    nz = bin_total > 0

    # Three side-by-side bars per bin
    w = BINWIDTH * 0.27
    x = centers
    # Color scheme: blues for ours (light=p@1, dark=p@2), pink for ReactOT
    C_OURS_P1 = '#7aa6dd'   # lighter blue
    C_OURS_P2 = '#3a6dbf'   # dark blue (matches the prior figure)
    C_OT      = '#cc3366'

    ax_acc.bar(x[nz] - w, rates_p1[nz]*100, w, color=C_OURS_P1,
                edgecolor='white', label='ours pass@1')
    ax_acc.bar(x[nz],      rates_p2[nz]*100, w, color=C_OURS_P2,
                edgecolor='white', label='ours pass@2')
    ax_acc.bar(x[nz] + w,  rates_ot[nz]*100, w, color=C_OT,
                edgecolor='white', label='ReactOT (deterministic)')
    # k/N annotations (just k as numerator for compactness; total is in
    # the bottom-panel bin label).
    for xi, n_, k1, k2, kot, r1, r2, rot in zip(
            x, bin_total, bin_p1, bin_p2, bin_ot, rates_p1, rates_p2, rates_ot):
        if n_ == 0: continue
        ax_acc.text(xi - w, r1*100 + 1.6, f'{k1}', ha='center', va='bottom',
                     fontsize=8.5, color='#1f3d77', fontweight='bold')
        ax_acc.text(xi,     r2*100 + 1.6, f'{k2}', ha='center', va='bottom',
                     fontsize=8.5, color='#0d2752', fontweight='bold')
        ax_acc.text(xi + w, rot*100 + 1.6, f'{kot}', ha='center', va='bottom',
                     fontsize=8.5, color='#7b1f3d', fontweight='bold')

    overall_p1 = p1_arr.mean() * 100
    overall_p2 = p2_arr.mean() * 100
    overall_ot = ot_arr.mean() * 100
    ax_acc.axhline(overall_p1, color=C_OURS_P1, linestyle='--', lw=1.2, alpha=0.85,
                    label=f'ours pass@1 overall = {overall_p1:.1f}%')
    ax_acc.axhline(overall_p2, color=C_OURS_P2, linestyle='--', lw=1.2, alpha=0.85,
                    label=f'ours pass@2 overall = {overall_p2:.1f}%')
    ax_acc.axhline(overall_ot, color=C_OT, linestyle='--', lw=1.2, alpha=0.85,
                    label=f'ReactOT overall = {overall_ot:.1f}%')
    ax_acc.set_ylabel('accuracy (% steps with $\\geq 1$ good IG)', fontsize=12)
    ax_acc.set_ylim(0, 110)
    ax_acc.set_yticks(range(0, 101, 20))
    ax_acc.set_xticks(x)
    ax_acc.set_xticklabels(labels, fontsize=11)
    ax_acc.grid(axis='y', linestyle=':', alpha=0.4)
    ax_acc.set_axisbelow(True)
    ax_acc.legend(loc='lower right', fontsize=9, framealpha=0.95, ncol=2)
    ax_acc.set_title(f'Pass@1 and Pass@2 accuracy by atom-count bin   '
                     f'(N = {len(joined)} steps; binwidth = {BINWIDTH} atoms; '
                     f'last bin pools 90+)',
                     fontsize=12, pad=12)

    # Bottom: bin populations + element annotation
    ax_n.bar(x[nz], bin_total[nz], width=BINWIDTH * 0.86, color='#7a7a7a',
              edgecolor='white', alpha=0.9)
    bin_max = bin_total.max() if bin_total.max() > 0 else 1
    count_offset = bin_max * 0.06
    elem_offset  = bin_max * 0.32
    for xi, n_, eset in zip(x, bin_total, bin_elems):
        if n_ > 0:
            ax_n.text(xi, n_ + count_offset, f'{n_}',
                      ha='center', va='bottom',
                      fontsize=11, color='#222', fontweight='bold')
            ax_n.text(xi, n_ + elem_offset, elem_lines(eset, max_per_line=6),
                      ha='center', va='bottom',
                      fontsize=10.5, color='#1c4e80', fontweight='bold',
                      linespacing=1.25,
                      bbox=dict(boxstyle='round,pad=0.25',
                                facecolor='white', edgecolor='#cfd9e8',
                                alpha=0.85))
    ax_n.set_xlabel('atom count per step', fontsize=12)
    ax_n.set_ylabel('# steps in bin', fontsize=11)
    ax_n.set_xticks(x)
    ax_n.set_xticklabels(labels, fontsize=11)
    ax_n.set_ylim(0, bin_max * 2.1)
    ax_n.grid(axis='y', linestyle=':', alpha=0.4)
    ax_n.set_axisbelow(True)
    fig.tight_layout()

    out_png = OUT_DIR / 'pass1_pass2_with_reactot.png'
    out_pdf = OUT_DIR / 'pass1_pass2_with_reactot.pdf'
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nSaved:\n  {out_png}\n  {out_pdf}")

    # Print breakdown
    print(f"\n{'bin':>10s}   {'pass@1':>14s}    {'pass@2':>14s}    {'ReactOT':>14s}")
    for i, lbl in enumerate(labels):
        if bin_total[i] == 0: continue
        print(f"  {lbl:>10s}    {bin_p1[i]:>3d}/{bin_total[i]:<3d} ({rates_p1[i]*100:5.1f}%)  "
              f"{bin_p2[i]:>3d}/{bin_total[i]:<3d} ({rates_p2[i]*100:5.1f}%)  "
              f"{bin_ot[i]:>3d}/{bin_total[i]:<3d} ({rates_ot[i]*100:5.1f}%)")


if __name__ == '__main__':
    main()

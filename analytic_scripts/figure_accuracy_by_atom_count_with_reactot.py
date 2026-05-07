"""
Figure: verifier accuracy + react_OT accuracy side-by-side per
atom-count bin.

react_OT pass set (user-supplied; all other 147 steps treated as fail):
  Jackie_TS_19, pr16.carbocation_ts{2,3,4,5,7,11,12}

Verifier pass = column `pass` in the human-judged final_quality CSV.

Same 10-atom bins as figure_accuracy_by_atom_count.py, last bin pools
90+ atoms.

Output:
  appendix_perparation/figures/accuracy_by_atom_count_with_reactot.{png,pdf}
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
    # n_atoms + element set per step
    n_atoms = {}
    elements_per_step = {}
    for p in VIEWER_DIR.glob('*.html'):
        if p.name in ('flat_view.html', 'guess_quality.html', 'index.html'): continue
        m = re.search(r"const DATA = (\{.*?\});\n", p.read_text(), re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        n_atoms[data['step']] = data['n_atoms']
        # Use the first TS's xyz_elements as the canonical element list
        ts_list = data.get('ts_list', [])
        elems = set()
        if ts_list and ts_list[0].get('xyz_elements'):
            elems = set(ts_list[0]['xyz_elements'])
        elements_per_step[data['step']] = elems

    # Verifier pass labels
    lines = HUMAN_CSV.read_text().splitlines()
    header = lines[1].split(',')
    rows = [dict(zip(header, ln.split(','))) for ln in lines[2:] if ln.strip()]
    verifier_pass = {canonical(r['step']): int(r['pass']) for r in rows
                     if r.get('pass') in ('0', '1')}

    # Join: only steps with all three (n_atoms, verifier label)
    joined = []
    for step, n in n_atoms.items():
        if step in verifier_pass:
            joined.append((step, n, verifier_pass[step], 1 if step in REACT_OT_PASS else 0))

    n_arr = np.array([j[1] for j in joined])
    v_arr = np.array([j[2] for j in joined])
    ot_arr = np.array([j[3] for j in joined])
    print(f"{len(joined)} steps  |  verifier pass: {v_arr.sum()}/{len(joined)}  |  "
          f"react_OT pass: {ot_arr.sum()}/{len(joined)}")

    # Bins: [10,20), [20,30), ..., [80,90), 90+
    fixed_edges = np.arange(10, 90 + 1, BINWIDTH)
    n_bins = len(fixed_edges)
    centers = list(fixed_edges[:-1] + BINWIDTH / 2) + [95]
    centers = np.array(centers)
    labels  = [f"[{int(e)}, {int(e+BINWIDTH)})" for e in fixed_edges[:-1]] + ['90+']

    bin_total = np.zeros(n_bins, dtype=int)
    bin_v     = np.zeros(n_bins, dtype=int)
    bin_ot    = np.zeros(n_bins, dtype=int)
    bin_elems = [set() for _ in range(n_bins)]
    for j in joined:
        step, n, v, ot = j
        idx = (n_bins - 1) if n >= 90 else int((n - 10) // BINWIDTH)
        idx = max(0, min(idx, n_bins - 1))
        bin_total[idx] += 1
        bin_v[idx]     += int(v)
        bin_ot[idx]    += int(ot)
        bin_elems[idx] |= elements_per_step.get(step, set())

    # Element list per bin — sorted with common organic first then by symbol.
    # Returns a 2-line string when the list is long.
    PRIORITY = ['C', 'H', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']
    def elem_lines(eset, max_per_line=6):
        prio = [e for e in PRIORITY if e in eset]
        rest = sorted(e for e in eset if e not in PRIORITY)
        elems = prio + rest
        if not elems: return ''
        line1 = elems[:max_per_line]
        line2 = elems[max_per_line:]
        out = ','.join(line1)
        if line2:
            out += '\n' + ','.join(line2)
        return out

    def wilson(k, n, z=1.96):
        if n == 0: return (0, 0, 0)
        p = k / n
        denom = 1 + z*z/n
        c = (p + z*z/(2*n)) / denom
        h = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
        return p, max(0.0, c-h), min(1.0, c+h)

    rates_v  = []; lo_v  = []; hi_v  = []
    rates_ot = []; lo_ot = []; hi_ot = []
    for k_v, k_ot, n_ in zip(bin_v, bin_ot, bin_total):
        p, l, h = wilson(int(k_v),  int(n_));  rates_v.append(p);  lo_v.append(l);  hi_v.append(h)
        p, l, h = wilson(int(k_ot), int(n_));  rates_ot.append(p); lo_ot.append(l); hi_ot.append(h)
    rates_v  = np.array(rates_v);  lo_v  = np.array(lo_v);  hi_v  = np.array(hi_v)
    rates_ot = np.array(rates_ot); lo_ot = np.array(lo_ot); hi_ot = np.array(hi_ot)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, (ax_acc, ax_n) = plt.subplots(2, 1, figsize=(13, 9),
                                          sharex=True,
                                          gridspec_kw={'height_ratios': [2.6, 1.4]})
    nz = bin_total > 0
    w = BINWIDTH * 0.40
    x = centers
    ax_acc.bar(x[nz] - w/2, rates_v[nz]*100,  w, color='#3a6dbf', edgecolor='white',
                label='clean_v2 verifier')
    ax_acc.bar(x[nz] + w/2, rates_ot[nz]*100, w, color='#cc3366', edgecolor='white',
                label='react_OT')
    # (CI error bars removed for visual cleanliness — k/N annotations
    # below give the sample size and Wilson CI is documented in the
    # accompanying CSV.)
    # Annotate k/N above each pair of bars
    for xi, n_, kv, kot, rv, rot in zip(x, bin_total, bin_v, bin_ot, rates_v, rates_ot):
        if n_ == 0: continue
        ax_acc.text(xi - w/2, rv*100 + 1.8, f'{kv}/{n_}', ha='center', va='bottom',
                     fontsize=9, color='#1f3d77', fontweight='bold')
        ax_acc.text(xi + w/2, rot*100 + 1.8, f'{kot}/{n_}', ha='center', va='bottom',
                     fontsize=9, color='#7b1f3d', fontweight='bold')
    overall_v  = v_arr.mean()  * 100
    overall_ot = ot_arr.mean() * 100
    ax_acc.axhline(overall_v,  color='#3a6dbf', linestyle='--', lw=1.2, alpha=0.7,
                    label=f'verifier overall = {overall_v:.1f}%')
    ax_acc.axhline(overall_ot, color='#cc3366', linestyle='--', lw=1.2, alpha=0.7,
                    label=f'react_OT overall = {overall_ot:.1f}%')
    ax_acc.set_ylabel('accuracy (% steps with ≥ 1 good IG)', fontsize=11)
    ax_acc.set_ylim(0, 110)
    ax_acc.set_yticks(range(0, 101, 20))
    ax_acc.set_xticks(x)
    ax_acc.set_xticklabels(labels, fontsize=10)
    ax_acc.grid(axis='y', linestyle=':', alpha=0.4)
    ax_acc.set_axisbelow(True)
    ax_acc.legend(loc='center left', fontsize=9, framealpha=0.95)
    ax_acc.set_title(f'Accuracy by atom-count bin   '
                     f'(N = {len(joined)} steps; binwidth = {BINWIDTH} atoms; '
                     f'last bin pools 90+)',
                     fontsize=12, pad=12)

    # Bin populations + element-set annotation above each bar (2 lines).
    # Layout: count number sits just above the bar; element list sits a
    # few units higher so it has breathing room. Larger font for the
    # element symbols and a soft white background for legibility.
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
    ax_n.set_ylim(0, bin_max * 2.1)   # extra headroom for taller element labels
    ax_n.grid(axis='y', linestyle=':', alpha=0.4)
    ax_n.set_axisbelow(True)
    fig.tight_layout()

    out_png = OUT_DIR / 'accuracy_by_atom_count_with_reactot.png'
    out_pdf = OUT_DIR / 'accuracy_by_atom_count_with_reactot.pdf'
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nSaved:\n  {out_png}\n  {out_pdf}")

    # Print
    print(f"\n{'bin':>10s}   {'verifier':>16s}   {'react_OT':>16s}")
    for i, lbl in enumerate(labels):
        if bin_total[i] == 0: continue
        v_pct  = rates_v[i]  * 100
        ot_pct = rates_ot[i] * 100
        print(f"  {lbl:>10s}    {bin_v[i]:>4d}/{bin_total[i]:<4d} ({v_pct:5.1f}%)   "
              f"{bin_ot[i]:>4d}/{bin_total[i]:<4d} ({ot_pct:5.1f}%)")


if __name__ == '__main__':
    main()

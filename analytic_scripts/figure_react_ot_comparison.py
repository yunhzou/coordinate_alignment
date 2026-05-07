"""
Figure: side-by-side per-step pass/fail comparison between the
clean_v2 verifier (human-judged) and the react_OT baseline (user
visual inspection on a 10-step subset).

Subset (the steps you visually inspected in alignment_view):
  Jackie_TS_19, pr16.carbocation_ts{2,3,4,5,7,8,9,11,12}

React OT labels are user-supplied:
  PASS: Jackie_TS_19, pr16.ts{2,3,4,5,7,11,12}
  FAIL: pr16.ts{8,9}

Verifier labels are read from the human-judged CSV
  appendix_perparation/analtics/final_quality_measurement-humanversion (1).csv
(`pass` column = at-least-one-of-IG1-IG2 judged Good).

Output:
  appendix_perparation/figures/react_ot_vs_verifier.{png,pdf}
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


HUMAN_CSV = PROJECT_ROOT / 'appendix_perparation' / 'analtics' / 'final_quality_measurement-humanversion (1).csv'
OUT_DIR   = PROJECT_ROOT / 'appendix_perparation' / 'figures'

# 10-step subset visually inspected in alignment_view, with react_OT labels
SUBSET = [
    ('Jackie_TS_19',           1),
    ('pr16.carbocation_ts2',   1),
    ('pr16.carbocation_ts3',   1),
    ('pr16.carbocation_ts4',   1),
    ('pr16.carbocation_ts5',   1),
    ('pr16.carbocation_ts7',   1),
    ('pr16.carbocation_ts8',   0),
    ('pr16.carbocation_ts9',   0),
    ('pr16.carbocation_ts11',  1),
    ('pr16.carbocation_ts12',  1),
]


def main():
    lines = HUMAN_CSV.read_text().splitlines()
    header = lines[1].split(',')
    rows = [dict(zip(header, ln.split(','))) for ln in lines[2:] if ln.strip()]
    canon = lambda s: re.sub(r'\s*\(.*?\)\s*$', '', s).strip()
    by_step = {canon(r['step']): r for r in rows}

    short_label = {
        'Jackie_TS_19':           'Jackie_TS_19',
        'pr16.carbocation_ts2':   'pr16.ts2',
        'pr16.carbocation_ts3':   'pr16.ts3',
        'pr16.carbocation_ts4':   'pr16.ts4',
        'pr16.carbocation_ts5':   'pr16.ts5',
        'pr16.carbocation_ts7':   'pr16.ts7',
        'pr16.carbocation_ts8':   'pr16.ts8',
        'pr16.carbocation_ts9':   'pr16.ts9',
        'pr16.carbocation_ts11':  'pr16.ts11',
        'pr16.carbocation_ts12':  'pr16.ts12',
    }

    steps, react_ot, verifier, ig1, ig2 = [], [], [], [], []
    for s, r_ot in SUBSET:
        steps.append(short_label[s])
        react_ot.append(r_ot)
        v = int(by_step[s]['pass']) if s in by_step else -1
        verifier.append(v)
        ig1.append(int(by_step[s]['IG#1']) if s in by_step else -1)
        ig2.append(int(by_step[s]['IG#2']) if s in by_step else -1)

    print(f"{'step':22s} {'verifier':>9s} {'react_OT':>9s}  IG#1 IG#2")
    for s, v, ot, a, b in zip(steps, verifier, react_ot, ig1, ig2):
        print(f"{s:22s} {v:>9} {ot:>9}    {a}    {b}")

    n = len(steps)
    v_pass = sum(verifier); ot_pass = sum(react_ot)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, (ax_step, ax_agg) = plt.subplots(1, 2, figsize=(14, 5.5),
                                           gridspec_kw={'width_ratios': [3, 1]})

    # Per-step grouped bars
    x = np.arange(n)
    w = 0.38
    b1 = ax_step.bar(x - w/2, verifier, w, color='#3a6dbf',
                      edgecolor='white', label='clean_v2 verifier (human-judged)')
    b2 = ax_step.bar(x + w/2, react_ot, w, color='#cc3366',
                      edgecolor='white', label='react_OT (visual inspection)')
    ax_step.set_xticks(x)
    ax_step.set_xticklabels(steps, rotation=30, ha='right', fontsize=10)
    ax_step.set_ylim(0, 1.18)
    ax_step.set_yticks([0, 1])
    ax_step.set_yticklabels(['fail', 'pass'])
    ax_step.set_title(f'Per-step pass/fail on the {n}-step inspected subset', pad=10)
    ax_step.legend(loc='lower right', fontsize=9, framealpha=0.95)
    ax_step.grid(axis='y', linestyle=':', alpha=0.4)
    ax_step.set_axisbelow(True)
    # Mark disagreement steps
    for i, (v, ot) in enumerate(zip(verifier, react_ot)):
        if v != ot:
            ax_step.axvspan(i - 0.5, i + 0.5, color='#fff2cc', alpha=0.5, zorder=0)
            ax_step.text(i, 1.10, '⚠ disagree', ha='center', fontsize=9, color='#aa6600')

    # Aggregate pass-rate bars
    rates = [v_pass / n * 100, ot_pass / n * 100]
    methods = ['verifier', 'react_OT']
    colors  = ['#3a6dbf', '#cc3366']
    bars = ax_agg.bar(methods, rates, color=colors, edgecolor='white', width=0.6)
    for b, r, k in zip(bars, rates, [v_pass, ot_pass]):
        ax_agg.text(b.get_x() + b.get_width()/2, r + 1, f'{k}/{n}\n({r:.0f}%)',
                     ha='center', va='bottom', fontsize=11)
    ax_agg.set_ylim(0, 110)
    ax_agg.set_yticks(range(0, 101, 20))
    ax_agg.set_ylabel('% of subset passed')
    ax_agg.set_title('Aggregate pass rate', pad=10)
    ax_agg.grid(axis='y', linestyle=':', alpha=0.4)
    ax_agg.set_axisbelow(True)

    fig.suptitle('clean_v2 verifier vs react_OT on the alignment-inspected subset',
                  fontsize=12.5, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_png = OUT_DIR / 'react_ot_vs_verifier.png'
    out_pdf = OUT_DIR / 'react_ot_vs_verifier.pdf'
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nSaved:\n  {out_png}\n  {out_pdf}")


if __name__ == '__main__':
    main()

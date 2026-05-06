"""
For each step, compute the GT's own picked-mode metrics:
  - bond_overlap (on the GT mode against its own broken/formed bonds)
  - rxn_overlap
  - core_fraction
  - n_imag (does GT have a clean single imaginary mode?)

If GT itself has LOW values, our broken/formed/core identification is
suspect for that step — the metric is essentially measuring against
a wrong target.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import json, re, time
from pathlib import Path
import numpy as np

from improve_ranker import load_step, imag_modes, cos_sim, rk_aggressive_v1

SRC = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_viewer')


def main():
    files = sorted(SRC.glob('*.html'))
    files = [f for f in files
             if f.name not in ('index.html', 'flat_view.html', 'guess_quality.html')]
    print(f"Loading {len(files)} steps...")
    t0 = time.time()
    rows = []
    for hp in files:
        text = hp.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        gt = next((t for t in data['ts_list'] if t['label']=='groundtruth'
                   and t.get('modes')), None)
        if gt is None: continue
        gt_mode = gt['modes'][gt['default_mode_idx']]
        gt_imag = [m for m in gt['modes'] if m['freq'] < 0]
        rows.append(dict(
            step=data['step'],
            gt_b=gt_mode.get('bond_overlap', 0),
            gt_r=gt_mode.get('rxn_overlap', 0),
            gt_c=gt_mode.get('core_fraction', 0),
            gt_freq=gt_mode['freq'],
            gt_n_imag=len(gt_imag),
            n_broken=len(data.get('broken_bonds', [])),
            n_formed=len(data.get('formed_bonds_R', [])),
            n_core=len(data.get('core_atoms', [])),
            n_atoms=data.get('n_atoms', 0),
        ))

    if not rows:
        print("no rows"); return
    print(f"loaded {len(rows)} steps in {time.time()-t0:.0f}s\n")

    arr_b = np.array([r['gt_b'] for r in rows])
    arr_c = np.array([r['gt_c'] for r in rows])
    arr_r = np.array([r['gt_r'] for r in rows])

    print("GT picked-mode metrics distribution:")
    print(f"{'metric':18s}  mean   median  <0.1   <0.2   <0.3   <0.5   <0.7")
    for name, a in (('bond_overlap', arr_b), ('rxn_overlap', arr_r),
                    ('core_fraction', arr_c)):
        print(f"{name:18s}  {a.mean():.3f}  {np.median(a):.3f}   "
              f"{(a<0.1).sum():3}    {(a<0.2).sum():3}    "
              f"{(a<0.3).sum():3}    {(a<0.5).sum():3}    {(a<0.7).sum():3}")

    print(f"\nGT n_imag distribution (clean TS = 1 imag mode):")
    from collections import Counter
    n_imag_dist = Counter(r['gt_n_imag'] for r in rows)
    for k in sorted(n_imag_dist.keys()):
        print(f"  n_imag={k}: {n_imag_dist[k]} steps")

    print(f"\nSteps where GT is suspicious (gt_b<0.1 OR gt_c<0.1):")
    sus = [r for r in rows if r['gt_b'] < 0.1 or r['gt_c'] < 0.1]
    print(f"  {len(sus)}/{len(rows)} steps")
    for r in sus[:20]:
        print(f"    {r['step']:50s}  b={r['gt_b']:.2f}  r={r['gt_r']:.2f}  "
              f"c={r['gt_c']:.2f}  freq={r['gt_freq']:.0f}  "
              f"n_imag={r['gt_n_imag']}  broken={r['n_broken']} formed={r['n_formed']}")

    # Cross-check: for steps where GT_c < 0.1, what's the oracle / ranker?
    print(f"\nLow-GT-c steps: how do they show up in ranker performance?")
    # Need to load step data with gt_disp + run ranker
    from improve_ranker import load_step
    for r in sus[:10]:
        hp = SRC / f"{r['step']}.html"
        sd = load_step(hp)
        if sd is None: continue
        gt_disp = sd['gt_disp']
        oracle = max((cos_sim(np.asarray(m['disp']), gt_disp)
                      for ts in sd['igs'] for m in ts['modes']), default=0)
        ranked = rk_aggressive_v1(sd, 0.7, 0.10, 1.0, 0.2)
        if ranked:
            top2 = max(cos_sim(np.asarray(p[1]['disp']), gt_disp)
                       for p in ranked[:2])
        else: top2 = 0
        print(f"    {r['step']:50s}  GT_b={r['gt_b']:.2f} GT_c={r['gt_c']:.2f}  "
              f"oracle={oracle:.2f}  ranker={top2:.2f}")

    # Save full csv
    out_csv = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_analysis/gt_quality.csv')
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    import csv
    with out_csv.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\nCSV: {out_csv}")


if __name__ == '__main__':
    main()

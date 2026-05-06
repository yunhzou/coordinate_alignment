"""
Re-evaluate the ranker excluding steps where GT itself is suspect:
  1. broken=0 AND formed=0 (no bond changes detected)
  2. GT picked-mode bond_overlap < 0.1 (GT mode doesn't project on bonds)
  3. GT picked-mode core_fraction < 0.1 (core_atoms misidentified)
  4. GT n_imag = 0 (no imaginary mode at all)

For each filter, report ranker mean and ≥thr stats vs oracle on the
remaining steps. If the gap shrinks, we know the metric/identification
is the problem — not the ranker.
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
        sd = load_step(hp)
        if sd is None: continue
        gt_disp = sd['gt_disp']
        oracle = max((cos_sim(np.asarray(m['disp']), gt_disp)
                      for ts in sd['igs'] for m in ts['modes']), default=0)
        oracle_imag = max((cos_sim(np.asarray(m['disp']), gt_disp)
                           for ts in sd['igs']
                           for m in ts['modes'] if m['freq'] < 0), default=0)
        ranked = rk_aggressive_v1(sd, 0.7, 0.10, 1.0, 0.2)
        if ranked:
            top2 = max(cos_sim(np.asarray(p[1]['disp']), gt_disp)
                       for p in ranked[:2])
        else: top2 = 0
        rows.append(dict(
            step=data['step'],
            gt_b=gt_mode.get('bond_overlap', 0),
            gt_c=gt_mode.get('core_fraction', 0),
            gt_r=gt_mode.get('rxn_overlap', 0),
            n_broken=len(data.get('broken_bonds', [])),
            n_formed=len(data.get('formed_bonds_R', [])),
            gt_freq=gt_mode['freq'],
            oracle=oracle, oracle_imag=oracle_imag, ranker=top2,
        ))
    print(f"loaded {len(rows)} steps\n")

    def stats(label, subset):
        n = len(subset)
        if n == 0:
            print(f"{label:50s}  N=0"); return
        o  = np.array([r['oracle'] for r in subset])
        oi = np.array([r['oracle_imag'] for r in subset])
        rk = np.array([r['ranker'] for r in subset])
        print(f"{label:50s}  N={n}")
        print(f"  oracle      mean={o.mean():.3f}  ≥0.7={(o>=0.7).mean()*100:.1f}%  "
              f"≥0.5={(o>=0.5).mean()*100:.1f}%  ≥0.3={(o>=0.3).mean()*100:.1f}%")
        print(f"  oracle_imag mean={oi.mean():.3f}  ≥0.7={(oi>=0.7).mean()*100:.1f}%  "
              f"≥0.5={(oi>=0.5).mean()*100:.1f}%  ≥0.3={(oi>=0.3).mean()*100:.1f}%")
        print(f"  ranker      mean={rk.mean():.3f}  ≥0.7={(rk>=0.7).mean()*100:.1f}%  "
              f"≥0.5={(rk>=0.5).mean()*100:.1f}%  ≥0.3={(rk>=0.3).mean()*100:.1f}%")
        print(f"  gap         mean={(o-rk).mean():.3f}")

    stats("ALL steps", rows)
    print()
    # Filter 1: at least 1 broken/formed bond identified
    stats("broken+formed >= 1", [r for r in rows if r['n_broken']+r['n_formed']>=1])
    print()
    # Filter 2: GT bond_overlap ≥ 0.1
    stats("GT bond_overlap >= 0.1", [r for r in rows if r['gt_b']>=0.1])
    print()
    # Filter 3: GT core_fraction ≥ 0.1
    stats("GT core_fraction >= 0.1", [r for r in rows if r['gt_c']>=0.1])
    print()
    # Filter 4: GT freq < 0
    stats("GT has imaginary mode", [r for r in rows if r['gt_freq']<0])
    print()
    # Combined: all reasonable
    stats("ALL filters (b>=0.1 AND c>=0.1 AND nbonds>=1 AND imag)",
          [r for r in rows
           if r['gt_b']>=0.1 and r['gt_c']>=0.1 and r['n_broken']+r['n_formed']>=1
           and r['gt_freq']<0])
    print()
    # Sus-only
    sus = [r for r in rows if not (r['gt_b']>=0.1 and r['gt_c']>=0.1
                                    and r['n_broken']+r['n_formed']>=1
                                    and r['gt_freq']<0)]
    stats("EXCLUDED ONLY (suspicious GT)", sus)


if __name__ == '__main__':
    main()

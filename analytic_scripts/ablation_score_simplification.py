"""
Ablation: does S(g) = β (just bond_overlap) give the same top-2
selection as the current S(g) = β × (1+w_r·ρ) × (1+w_c·κ) / n_imag^p?

Both rankers share:
  - the n_imag ≤ 2 filter
  - the rxn_overlap ≥ 0.10 filter
  - the greedy mass-weighted-cosine diversity penalty (α=0.7)
The only thing that changes is what we plug into the score.

For each step we report:
  - top-2 IG labels under each ranker
  - whether they agree as sets
Then aggregate: how many steps differ, and how many human-judged
pass-labels are affected.
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

from ranker import (mass_weighted_cos, imag_modes, ATOMIC_MASS,
                    rk_clean_v2)


VIEWER_DIR = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
HUMAN_CSV  = PROJECT_ROOT / 'appendix_perparation' / 'analtics' / 'final_quality_measurement-humanversion (1).csv'


def _diversity_select(candidates, alpha, elements, k=2):
    """Same greedy mass-weighted-cosine diversity rule used in clean_v2."""
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < k:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                mw = mass_weighted_cos(cand[1]['disp'], sel[1]['disp'], elements)
                score *= max(0.0, 1.0 - alpha * mw)
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


def rk_beta_only(igs, elements, alpha=0.7, min_rxn=0.10, max_imag=2, k=2):
    """Stripped ranker: S(g) = β. Same hard priors and same diversity rule."""
    cands = []
    for ts in igs:
        imag = imag_modes(ts)
        if not imag: continue
        n_imag = len(imag)
        if n_imag > max_imag: continue
        picked = max(imag, key=lambda m: m.get('bond_overlap', 0.0))
        if picked.get('rxn_overlap', 0.0) < min_rxn: continue
        score = picked.get('bond_overlap', 0.0)
        cands.append((score, picked, ts['label']))
    if not cands:
        # Fallback: bond_overlap top-k on unfiltered pool
        scored = []
        for ts in igs:
            imag = imag_modes(ts)
            modes = imag if imag else ts['modes']
            if not modes: continue
            picked = max(modes, key=lambda m: m.get('bond_overlap', 0.0))
            scored.append((picked.get('bond_overlap', 0.0), picked, ts['label']))
        scored.sort(key=lambda t: -t[0])
        return [(t[2], t[1]) for t in scored[:k]]
    selected = _diversity_select(cands, alpha, elements, k)
    return [(t[2], t[1]) for t in selected]


def get_clean_v2_labels(igs, elements):
    ranked = rk_clean_v2(igs, elements, k=2)
    return [t['label'] for t, _ in ranked[:2]]


def get_beta_only_labels(igs, elements):
    ranked = rk_beta_only(igs, elements, k=2)
    return [lbl for lbl, _ in ranked[:2]]


def main():
    files = sorted(VIEWER_DIR.glob('*.html'))
    files = [f for f in files
             if f.name not in ('flat_view.html', 'guess_quality.html', 'index.html')
             and not f.name.startswith('oracle_view')]
    print(f"Reading {len(files)} per-step HTMLs...")

    # Human pass labels (per-IG)
    human = {}
    if HUMAN_CSV.exists():
        lines = HUMAN_CSV.read_text().splitlines()
        header = lines[1].split(',')
        for ln in lines[2:]:
            if not ln.strip(): continue
            r = dict(zip(header, ln.split(',')))
            step = re.sub(r'\s*\(.*?\)\s*$', '', r['step']).strip()
            try:
                human[step] = (int(r['IG#1']), int(r['IG#2']), int(r['pass']))
            except Exception: pass

    same_set = 0; diff_set = 0
    same_order = 0; diff_order = 0
    rows = []
    for hp in files:
        text = hp.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        gt = next((t for t in data['ts_list']
                   if t['label']=='groundtruth' and t.get('modes')), None)
        if gt is None: continue
        elements = gt['xyz_elements']
        igs = [t for t in data['ts_list']
               if t['label']!='groundtruth' and t.get('modes')]
        if not igs: continue

        cv2 = get_clean_v2_labels(igs, elements)
        b   = get_beta_only_labels(igs, elements)
        cv2_set = set(cv2); b_set = set(b)
        if cv2_set == b_set:
            same_set += 1
            if tuple(cv2) == tuple(b): same_order += 1
            else: diff_order += 1
        else:
            diff_set += 1
        rows.append((data['step'], cv2, b))

    n = same_set + diff_set
    print(f"\nN = {n} steps with both rankers run successfully")
    print(f"Same {{IG#1, IG#2}} set:        {same_set}/{n}  ({100*same_set/n:.1f}%)")
    print(f"  -> same order (top-1, top-2): {same_order}/{n}  ({100*same_order/n:.1f}%)")
    print(f"  -> different order:           {diff_order}/{n}  ({100*diff_order/n:.1f}%)")
    print(f"Different set (different IGs): {diff_set}/{n}  ({100*diff_set/n:.1f}%)")

    # For the differing-set cases, show details
    if diff_set:
        print(f"\nSteps where the SET of top-2 IGs differs:")
        print(f"  {'step':46s}  {'clean_v2 top-2':22s}  {'β-only top-2':22s}")
        for step, cv2, b in rows:
            if set(cv2) == set(b): continue
            print(f"  {step:46s}  {','.join(cv2):22s}  {','.join(b):22s}")


if __name__ == '__main__':
    main()

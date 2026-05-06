"""
For each step, find the IG with the highest gt_alignment (= 'best IG').
Compute its rank under several candidate signals; the lower the rank,
the better the signal at identifying the best IG.

Aim: identify a signal where best-IG ranks ≤ 2 in many steps.
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np

from improve_ranker import load_step, imag_modes, mass_weighted_cos, cos_sim

SRC = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_viewer')


def kabsch_rmsd(A, B):
    """Heavy-atom-aware Kabsch RMSD between two coord arrays of same shape."""
    A = np.asarray(A); B = np.asarray(B)
    if A.shape != B.shape: return 1e9
    Ac = A - A.mean(axis=0); Bc = B - B.mean(axis=0)
    H = Ac.T @ Bc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Aaln = Ac @ R.T
    return float(np.sqrt(((Aaln - Bc) ** 2).sum() / len(A)))


def best_imag_align(ts, gt_disp):
    return max((cos_sim(np.asarray(m['disp']), gt_disp)
                for m in ts['modes'] if m['freq'] < 0), default=0)


def features_per_ig(sd):
    """Returns list of (label, feat-dict) for each IG."""
    feats = []
    igs = sd['igs']
    # TS structure consensus: mean Kabsch RMSD to other TSs
    ts_xyzs = [np.asarray(t['xyz_coords']) for t in igs]
    n = len(igs)
    rmsd_to_others = np.zeros(n)
    for i in range(n):
        rs = [kabsch_rmsd(ts_xyzs[i], ts_xyzs[j]) for j in range(n) if j != i]
        rmsd_to_others[i] = np.mean(rs) if rs else 0
    # Mode peer similarity (mass-weighted)
    elements = sd['elR']
    mode_peer = np.zeros(n)
    for i, ts in enumerate(igs):
        imag = imag_modes(ts)
        if not imag: continue
        m_i = max(imag, key=lambda m: m.get('bond_overlap', 0))
        s = 0; cnt = 0
        for j, ts_j in enumerate(igs):
            if i == j: continue
            imag_j = imag_modes(ts_j)
            if not imag_j: continue
            m_j = max(imag_j, key=lambda m: m.get('bond_overlap', 0))
            s += mass_weighted_cos(m_i['disp'], m_j['disp'], elements)
            cnt += 1
        mode_peer[i] = s / max(cnt, 1)

    for i, ts in enumerate(igs):
        imag = imag_modes(ts)
        if not imag:
            feats.append((ts['label'], {}))
            continue
        m_best = max(imag, key=lambda m: m.get('bond_overlap', 0))
        b = m_best.get('bond_overlap', 0)
        r = m_best.get('rxn_overlap', 0)
        c = m_best.get('core_fraction', 0)
        # sum of bond_overlap across imag modes
        bond_sum = sum(m.get('bond_overlap', 0) for m in imag)
        # max rxn across imag modes
        max_rxn = max((m.get('rxn_overlap', 0) for m in imag), default=0)
        feats.append((ts['label'], {
            'bond_overlap': b,
            'rxn_overlap': r,
            'core_fraction': c,
            'b_x_1plusr': b * (1 + r),
            'b_x_1plusr_x_1plus02c': b * (1 + r) * (1 + 0.2 * c),
            'mode_peer': mode_peer[i],
            'rmsd_consensus': -rmsd_to_others[i],  # negate so higher=better
            'bond_sum_imag': bond_sum,
            'max_rxn_imag': max_rxn,
            'b_x_peer': b * mode_peer[i],
            'b_x_1plusr_x_peer': b * (1 + r) * mode_peer[i],
            'imag_count': -len(imag),  # fewer imag = better (negate)
            'lowest_freq': -abs(min(m['freq'] for m in imag)),  # less imag = better, negate
            'b_x_1plusr_neg_rmsd': b * (1 + r) * (1.0 / (rmsd_to_others[i] + 0.5)),
        }))
    return feats


def rank_of_best(feats, best_label, key):
    """Rank (1=highest) of best-IG under feature `key`. Returns None if missing."""
    scored = [(lbl, f.get(key)) for lbl, f in feats if key in f]
    if not scored: return None
    scored.sort(key=lambda t: -(t[1] if t[1] is not None else -1e18))
    for r, (lbl, _) in enumerate(scored, 1):
        if lbl == best_label: return r
    return None


def main():
    files = sorted(SRC.glob('*.html'))
    files = [f for f in files
             if f.name not in ('index.html', 'flat_view.html', 'guess_quality.html')]
    print(f"Loading {len(files)} steps...")
    t0 = time.time()
    scope = []
    for hp in files:
        sd = load_step(hp)
        if sd: scope.append(sd)
    print(f"  loaded {len(scope)} steps in {time.time()-t0:.0f}s\n")

    all_keys = ['bond_overlap', 'rxn_overlap', 'core_fraction',
                'b_x_1plusr', 'b_x_1plusr_x_1plus02c',
                'mode_peer', 'rmsd_consensus', 'bond_sum_imag',
                'max_rxn_imag', 'b_x_peer', 'b_x_1plusr_x_peer',
                'imag_count', 'lowest_freq', 'b_x_1plusr_neg_rmsd']
    rank_dist = {k: [] for k in all_keys}
    for sd in scope:
        gt_disp = sd['gt_disp']
        # Find best IG (by best imag align)
        best = max(sd['igs'], key=lambda ts: best_imag_align(ts, gt_disp))
        best_label = best['label']
        feats = features_per_ig(sd)
        for k in all_keys:
            r = rank_of_best(feats, best_label, k)
            if r is not None: rank_dist[k].append(r)

    print(f"{'feature':32s}  rank=1   rank≤2   rank≤3   median")
    print('=' * 80)
    for k in all_keys:
        rs = rank_dist[k]
        if not rs: continue
        rs = np.array(rs)
        print(f"{k:32s}  {(rs==1).mean()*100:5.1f}%  {(rs<=2).mean()*100:5.1f}%  "
              f"{(rs<=3).mean()*100:5.1f}%  {np.median(rs):5.1f}")


if __name__ == '__main__':
    main()

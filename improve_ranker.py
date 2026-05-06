"""
Iterate on ranker improvements. All variants evaluated at k=2 (top-2 IGs)
on the same 160 BGCP steps.

Reads cached data from per-step HTMLs (no recompute of alignment, no xtb).
Loads R/P WBO matrices for dwbo-based variants. R↔P alignment runs once
per step (~0.5s) for variants needing mapping_RP.
"""
from __future__ import annotations
import json, re, time
from collections import Counter
from pathlib import Path
from typing import Callable

import numpy as np

from rxn_core_pq import align_from_arrays
from align_bgcp_coords import load_cached_xtb

SRC = Path('/Users/yunhengz/empty_for_claude/rxn_core/out/mode_viewer')
WORK_MODES = Path('/Users/yunhengz/empty_for_claude/rxn_core/work_modes')


def cos_sim(a, b):
    a = np.asarray(a).reshape(-1); b = np.asarray(b).reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9: return 0.0
    return abs(float(a @ b)) / (na * nb)


# ─── data loading ─────────────────────────────────────────────────────────

def load_step(step_html_path):
    """Returns dict: gt_disp, ig_tss (list of {label, modes_list, xyz}),
    broken_R, formed_R, core_R, n_atoms, mapping_RP, wboR, wboP."""
    text = step_html_path.read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m: return None
    data = json.loads(m.group(1))
    step = data['step']
    md = WORK_MODES / step
    if not (md / 'R' / 'wbo').exists() or not (md / 'P' / 'wbo').exists():
        return None
    elR, xyzR, wboR, _ = load_cached_xtb(md / 'R')
    elP, xyzP, wboP, _ = load_cached_xtb(md / 'P')
    rp = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP)
    mapping_RP = dict(rp['mapping'])
    gt = next((t for t in data['ts_list'] if t['label']=='groundtruth' and t.get('modes')), None)
    if gt is None: return None
    gt_disp = np.asarray(gt['modes'][gt['default_mode_idx']]['disp'])
    igs = [t for t in data['ts_list'] if t['label']!='groundtruth' and t.get('modes')]
    if not igs: return None
    return dict(
        step=step, gt_disp=gt_disp, igs=igs,
        broken_R=data['broken_bonds'], formed_R=data['formed_bonds_R'],
        core_R=data['core_atoms'], n_atoms=data['n_atoms'],
        mapping_RP=mapping_RP, wboR=wboR, wboP=wboP, xyzR=xyzR, elR=elR,
    )


# ─── shared helpers ───────────────────────────────────────────────────────

def imag_modes(ts):
    return [m for m in ts['modes'] if m['freq'] < 0]


def cos_modes(a_mode, b_mode):
    """Cosine similarity between two mode dicts' disps."""
    return cos_sim(np.asarray(a_mode['disp']), np.asarray(b_mode['disp']))


def evaluate_topk(ranker_fn, scope):
    """ranker_fn(step_data) -> list of (score, picked_mode, ig_label) sorted desc.
    Returns dict with top-1, top-2, top-3 distributions."""
    res = {'top1': [], 'top2': [], 'top3': [], 'oracle': []}
    for sd in scope:
        gt_disp = sd['gt_disp']
        # Oracle: best gt_align over all (IG, mode)
        oracle = max((cos_sim(np.asarray(m['disp']), gt_disp)
                      for ts in sd['igs'] for m in ts['modes']), default=0)
        res['oracle'].append(oracle)
        ranked = ranker_fn(sd)
        if not ranked:
            for k in ('top1','top2','top3'):
                res[k].append(0.0)
            continue
        aligns = [cos_sim(np.asarray(p[1]['disp']), gt_disp) for p in ranked[:3]]
        while len(aligns) < 3: aligns.append(0.0)
        res['top1'].append(aligns[0])
        res['top2'].append(max(aligns[:2]))
        res['top3'].append(max(aligns[:3]))
    for k in res: res[k] = np.array(res[k])
    return res


def summarise(name, results):
    """Cumulative coverage at thresholds for top-2 (the consumed output)."""
    t2 = results['top2']
    n = len(t2)
    cells = []
    for thr in (0.7, 0.5, 0.3):
        cells.append(f"{(t2>=thr).sum():>3}/{n}({100*(t2>=thr).mean():>4.1f}%)")
    return (f"{name:40s}  mean={t2.mean():.3f}  median={np.median(t2):.3f}  "
            f"≥0.7={cells[0]}  ≥0.5={cells[1]}  ≥0.3={cells[2]}")


# ─── BASELINE: bond_overlap ───────────────────────────────────────────────

def rk_bond_overlap(sd):
    picks = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        picks.append((picked.get('bond_overlap', 0), picked, ts['label']))
    picks.sort(key=lambda t: -t[0])
    return picks


# ─── CANDIDATE 1: bond_overlap with diversity-penalised top-2 ────────────
# Top-1 = highest bond_overlap as before. Top-2 onwards: penalise by
# similarity to already-selected modes.

def rk_bond_diverse(sd, alpha=1.0):
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates: return []
    # Greedy diverse selection
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                penalty = alpha * cos_modes(cand[1], sel[1])
                score *= max(0.0, 1.0 - penalty)
            if score > best_score:
                best_score = score
                best = cand
        selected.append(best)
        remaining.remove(best)
    return selected


# ─── CANDIDATE 2: frequency-band consensus ────────────────────────────────
# Pick top-1 by bond_overlap. For top-2, prefer IGs whose default mode
# frequency is close to the MEDIAN imag freq across the IG pool's picks.

def rk_freq_consensus(sd, freq_window=100.0):
    picks = []
    freqs = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        picks.append((picked.get('bond_overlap', 0), picked, ts['label']))
        if picked['freq'] < 0:
            freqs.append(picked['freq'])
    if not picks: return []
    if not freqs:
        picks.sort(key=lambda t: -t[0]); return picks
    median_freq = float(np.median(freqs))
    # Score = bond_overlap × Gaussian-similarity to median freq
    scored = []
    for s, m, lbl in picks:
        df = abs(m['freq'] - median_freq) if m['freq'] < 0 else freq_window * 2
        gauss = np.exp(-(df / freq_window) ** 2)
        scored.append((s * (0.5 + 0.5 * gauss), m, lbl))
    scored.sort(key=lambda t: -t[0])
    return scored


# ─── CANDIDATE 3: cross-IG mode centroid ──────────────────────────────────
# Compute centroid of all IG default-mode displacements (as-is; sign
# arbitrariness is handled by aligning each to the centroid before
# averaging via |cos|). Pick IG with highest |cos similarity| to centroid.

def rk_centroid(sd):
    picks = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        picks.append((picked, ts['label']))
    if not picks: return []
    # Build centroid in flattened 3N space; flip signs to align with first
    flats = [np.asarray(p[0]['disp']).reshape(-1) for p in picks]
    ref = flats[0]
    aligned = [ref]
    for f in flats[1:]:
        aligned.append(f * np.sign(np.dot(ref, f) + 1e-12))
    centroid = np.mean(aligned, axis=0)
    cn = np.linalg.norm(centroid)
    scored = []
    for (m, lbl), f in zip(picks, flats):
        score = abs(np.dot(f, centroid)) / max(np.linalg.norm(f) * cn, 1e-9)
        scored.append((score, m, lbl))
    scored.sort(key=lambda t: -t[0])
    return scored


# ─── CANDIDATE 4: hybrid — bond_overlap with cross-IG-centroid tiebreak ──

def rk_bond_with_centroid_tiebreak(sd, eps=0.05):
    # Compute candidates with bond_overlap
    picks = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        picks.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not picks: return []
    # Centroid
    flats = [np.asarray(p[1]['disp']).reshape(-1) for p in picks]
    ref = flats[0]
    aligned = [ref]
    for f in flats[1:]:
        aligned.append(f * np.sign(np.dot(ref, f) + 1e-12))
    centroid = np.mean(aligned, axis=0)
    cn = np.linalg.norm(centroid)
    # Tiebreak score: bond_ov + small * centroid_sim
    scored = []
    for (b, m, lbl), f in zip(picks, flats):
        c = abs(np.dot(f, centroid)) / max(np.linalg.norm(f) * cn, 1e-9)
        scored.append((b + eps * c, m, lbl))
    scored.sort(key=lambda t: -t[0])
    return scored


# ─── CANDIDATE 5: bond_overlap × |freq| × core_fraction ───────────────────

def rk_bond_freq_core(sd):
    picks = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        # Pick mode within IG by combined score
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0)
                     * m.get('core_fraction', 0))
        score = (picked.get('bond_overlap', 0)
                 * picked.get('core_fraction', 0)
                 * (abs(picked['freq']) / 1000.0 if picked['freq'] < 0 else 0))
        picks.append((score, picked, ts['label']))
    picks.sort(key=lambda t: -t[0])
    return picks


# ─── CANDIDATE 6: mode-mode consensus among IG picks ─────────────────────
# Score each IG by mean cosine to all OTHER IGs' picks. The IG most
# typical of the pool wins (peer-consensus).

def rk_peer_consensus(sd):
    picks = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        picks.append((picked, ts['label']))
    if not picks: return []
    flats = [np.asarray(p[0]['disp']).reshape(-1) for p in picks]
    norms = [np.linalg.norm(f) for f in flats]
    n = len(flats)
    scored = []
    for i, (m, lbl) in enumerate(picks):
        s = 0
        for j in range(n):
            if i == j: continue
            d = abs(float(flats[i] @ flats[j])) / max(norms[i] * norms[j], 1e-9)
            s += d
        s /= max(n - 1, 1)
        scored.append((s, m, lbl))
    scored.sort(key=lambda t: -t[0])
    return scored


# ─── CANDIDATE 7: peer + bond_ov filter ───────────────────────────────────
# Restrict to IGs with bond_ov ≥ threshold, then peer-consensus among them.

def rk_peer_filtered(sd, bond_min=0.3):
    picks = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('bond_overlap', 0) >= bond_min:
            picks.append((picked, ts['label']))
    if len(picks) < 2:
        # fall back to bond_overlap
        return rk_bond_overlap(sd)
    flats = [np.asarray(p[0]['disp']).reshape(-1) for p in picks]
    norms = [np.linalg.norm(f) for f in flats]
    n = len(flats)
    scored = []
    for i, (m, lbl) in enumerate(picks):
        s = 0
        for j in range(n):
            if i == j: continue
            d = abs(float(flats[i] @ flats[j])) / max(norms[i] * norms[j], 1e-9)
            s += d
        s /= max(n - 1, 1)
        scored.append((s, m, lbl))
    scored.sort(key=lambda t: -t[0])
    return scored


# ─── more diversity variants ──────────────────────────────────────────────

def rk_div_hard_cap(sd, cos_thr):
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    candidates.sort(key=lambda t: -t[0])
    if not candidates: return []
    selected = [candidates[0]]
    for cand in candidates[1:]:
        if all(cos_modes(cand[1], s[1]) < cos_thr for s in selected):
            selected.append(cand)
            if len(selected) >= 3: break
    while len(selected) < 3 and len(selected) < len(candidates):
        for c in candidates:
            if c not in selected:
                selected.append(c); break
    return selected


def rk_div_hard_cap_09(sd):  return rk_div_hard_cap(sd, 0.9)
def rk_div_hard_cap_085(sd): return rk_div_hard_cap(sd, 0.85)
def rk_div_hard_cap_08(sd):  return rk_div_hard_cap(sd, 0.8)


def rk_div_plus_freq(sd):
    """Diversity (alpha=0.5) but bond_overlap is weighted by freq-consensus."""
    candidates = []
    freqs = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
        if picked['freq'] < 0:
            freqs.append(picked['freq'])
    if not candidates: return []
    median_freq = float(np.median(freqs)) if freqs else -300.0
    boosted = []
    for s, m, lbl in candidates:
        df = abs(m['freq'] - median_freq) if m['freq'] < 0 else 1e6
        gauss = np.exp(-(df / 100.0) ** 2)
        boosted.append((s * (0.5 + 0.5 * gauss), m, lbl))
    selected = []
    remaining = list(boosted)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                penalty = 0.5 * cos_modes(cand[1], sel[1])
                score *= max(0.0, 1.0 - penalty)
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


def rk_div_squared(sd):
    """Diversity with bond_overlap squared (sharper preference for high score)."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates: return []
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0] ** 2
            for sel in selected:
                penalty = 0.5 * cos_modes(cand[1], sel[1])
                score *= max(0.0, 1.0 - penalty)
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


# ─── ATOMIC MASSES (for mass weighting) ───────────────────────────────────
ATOMIC_MASS = {
    'H':  1.008, 'He': 4.003,
    'Li': 6.94,  'Be': 9.012, 'B': 10.81, 'C': 12.011, 'N': 14.007,
    'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.99, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.948,
    'K': 39.098, 'Ca': 40.08, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942,
    'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693,
    'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.63, 'As': 74.922,
    'Se': 78.971, 'Br': 79.904, 'Kr': 83.798,
    'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906,
    'Mo': 95.95, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.906, 'Pd': 106.42,
    'Ag': 107.868, 'Cd': 112.414, 'In': 114.818, 'Sn': 118.710, 'Sb': 121.760,
    'Te': 127.60, 'I': 126.904, 'Xe': 131.293,
    'Cs': 132.905, 'Ba': 137.327, 'La': 138.906, 'Hf': 178.49, 'Ta': 180.948,
    'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084,
    'Au': 196.967, 'Hg': 200.592, 'Tl': 204.383, 'Pb': 207.2, 'Bi': 208.980,
}


def mass_weighted_cos(disp_a, disp_b, elements):
    """Cosine similarity with sqrt(mass) weighting per atom."""
    a = np.asarray(disp_a, dtype=float)
    b = np.asarray(disp_b, dtype=float)
    masses = np.array([ATOMIC_MASS.get(e, 12.0) for e in elements])
    sqm = np.sqrt(masses)[:, None]
    aw = (a * sqm).reshape(-1)
    bw = (b * sqm).reshape(-1)
    na = float(np.linalg.norm(aw)); nb = float(np.linalg.norm(bw))
    if na < 1e-9 or nb < 1e-9: return 0.0
    return abs(float(aw @ bw)) / (na * nb)


# ─── more candidates ──────────────────────────────────────────────────────

def rk_div_freq_filter(sd, freq_min=100, freq_max=3000, alpha=0.5):
    """Diversity but exclude IGs whose pick freq is outside [-freq_max, -freq_min]."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        f = picked['freq']
        if f >= 0 or abs(f) < freq_min or abs(f) > freq_max:
            # still include but with bond_overlap discount
            score = picked.get('bond_overlap', 0) * 0.5
        else:
            score = picked.get('bond_overlap', 0)
        candidates.append((score, picked, ts['label']))
    if not candidates: return []
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                score *= max(0.0, 1.0 - alpha * cos_modes(cand[1], sel[1]))
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


def rk_div_per_ig_top2(sd, alpha=0.5):
    """Each IG contributes its TOP 2 imag modes (not just top 1).
    Then global ranking with diversity penalty across all (IG, mode)."""
    candidates = []
    for ts in sd['igs']:
        imag = sorted(imag_modes(ts), key=lambda m: -m.get('bond_overlap', 0))
        if not imag:
            picked = max(ts['modes'], key=lambda m: m.get('bond_overlap', 0))
            candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
        else:
            for picked in imag[:2]:
                candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates: return []
    selected = []
    remaining = list(candidates)
    seen_labels = set()
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            # Skip if same IG already selected
            if cand[2] in seen_labels:
                continue
            score = cand[0]
            for sel in selected:
                score *= max(0.0, 1.0 - alpha * cos_modes(cand[1], sel[1]))
            if score > best_score:
                best_score = score; best = cand
        if best is None: break
        selected.append(best)
        seen_labels.add(best[2])
        remaining = [c for c in remaining if c[2] != best[2]]
    return selected


def rk_mass_weighted_div(sd, alpha=0.5):
    """Like diversity penalty but using mass-weighted cosine for similarity."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates: return []
    selected = []
    remaining = list(candidates)
    elements = sd['elR']
    while remaining and len(selected) < 3:
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


def rk_div_with_freq_match(sd, alpha=0.5, beta=0.3):
    """Diversity penalty + bonus for matching frequency to top-1 pick.
    Idea: top-1 sets a 'reference frequency', similar-freq IGs are
    more likely to be the same TS (= correct), but should still differ
    in mode (= diverse)."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates: return []
    candidates.sort(key=lambda t: -t[0])
    selected = [candidates[0]]
    ref_freq = candidates[0][1]['freq']
    remaining = candidates[1:]
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                score *= max(0.0, 1.0 - alpha * cos_modes(cand[1], sel[1]))
            # Freq-match bonus
            if cand[1]['freq'] < 0 and ref_freq < 0:
                df = abs(cand[1]['freq'] - ref_freq)
                score *= (1.0 + beta * np.exp(-(df / 100.0) ** 2))
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


def rk_div_normalized_pool(sd, alpha=0.5):
    """Diversity but normalize bond_overlap to [0, 1] within the pool of
    20 IGs (so the relative ranking matters more than absolute values)."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates: return []
    bond_max = max(c[0] for c in candidates) or 1.0
    candidates = [(s / bond_max, m, lbl) for s, m, lbl in candidates]
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                score *= max(0.0, 1.0 - alpha * cos_modes(cand[1], sel[1]))
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


# ─── CLUSTER / ENSEMBLE / RANK-FUSION CANDIDATES ──────────────────────────

def rk_cluster_consensus(sd, similarity_thr=0.7):
    """Cluster IGs by mode similarity (sign-blind). For each cluster,
    score = bond_ov_max × cluster_size. Pick representatives from the
    highest-scoring cluster, then second-highest, etc."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates: return []
    n = len(candidates)
    flats = [np.asarray(c[1]['disp']).reshape(-1) for c in candidates]
    norms = [np.linalg.norm(f) for f in flats]
    # Build similarity matrix
    sim = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if norms[i] > 1e-9 and norms[j] > 1e-9:
                sim[i, j] = abs(flats[i] @ flats[j]) / (norms[i] * norms[j])
    # Greedy clustering: each unvisited IG with highest bond_ov spawns a cluster
    order = sorted(range(n), key=lambda i: -candidates[i][0])
    cluster_id = [-1] * n
    cid = 0
    for i in order:
        if cluster_id[i] != -1: continue
        cluster_id[i] = cid
        for j in range(n):
            if cluster_id[j] == -1 and sim[i, j] >= similarity_thr:
                cluster_id[j] = cid
        cid += 1
    # Score each cluster
    clusters = {}
    for i, c in enumerate(cluster_id):
        clusters.setdefault(c, []).append(i)
    cluster_scores = []
    for c, members in clusters.items():
        max_bond = max(candidates[i][0] for i in members)
        size = len(members)
        cluster_scores.append((max_bond * np.sqrt(size), c, members))
    cluster_scores.sort(key=lambda t: -t[0])
    # Pick best from each cluster, in cluster-score order
    selected = []
    for _, c, members in cluster_scores:
        members.sort(key=lambda i: -candidates[i][0])
        selected.append(candidates[members[0]])
        if len(selected) >= 3: break
    return selected


def rk_rank_fusion(sd, alpha=0.5):
    """Sum-of-ranks across bond_ov, peer-similarity, centroid-alignment.
    Final IG ranking = sum of ranks (lower = better). Diversity penalty
    applied for top-2 selection."""
    picks = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        picks.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not picks: return []
    n = len(picks)
    flats = [np.asarray(p[1]['disp']).reshape(-1) for p in picks]
    norms = [np.linalg.norm(f) for f in flats]
    # Peer similarity (mean cos to others)
    peer = np.zeros(n)
    for i in range(n):
        s = 0; cnt = 0
        for j in range(n):
            if i == j or norms[i] < 1e-9 or norms[j] < 1e-9: continue
            s += abs(flats[i] @ flats[j]) / (norms[i] * norms[j])
            cnt += 1
        peer[i] = s / max(cnt, 1)
    # Centroid alignment
    ref = flats[0]
    aligned = [ref]
    for f in flats[1:]:
        aligned.append(f * np.sign(np.dot(ref, f) + 1e-12))
    centroid = np.mean(aligned, axis=0)
    cn = np.linalg.norm(centroid)
    cent = np.array([abs(np.dot(f, centroid)) / max(np.linalg.norm(f) * cn, 1e-9)
                     for f in flats])
    bond = np.array([p[0] for p in picks])
    # Ranks (lower = better)
    rank_bond = (-bond).argsort().argsort()
    rank_peer = (-peer).argsort().argsort()
    rank_cent = (-cent).argsort().argsort()
    fused = rank_bond + rank_peer + rank_cent
    # Convert to negative score for max-based diversity selection
    score = -fused.astype(float)
    candidates = [(score[i], picks[i][1], picks[i][2]) for i in range(n)]
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            sc = cand[0]
            for sel in selected:
                sc -= 5 * cos_modes(cand[1], sel[1])  # additive penalty
            if sc > best_score:
                best_score = sc; best = cand
        selected.append(best); remaining.remove(best)
    return selected


def rk_div_iterative(sd, alpha=0.5, n_iter=3):
    """Iterative refinement: start with bond_overlap top-2, recompute
    'reference centroid', re-rank by alignment with centroid, repeat."""
    picks = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        picks.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not picks: return []
    flats = [np.asarray(p[1]['disp']).reshape(-1) for p in picks]
    norms = np.array([np.linalg.norm(f) for f in flats])
    bond = np.array([p[0] for p in picks])
    # Initial: top-2 by bond_ov
    sel_idx = list(np.argsort(-bond)[:2])
    for _ in range(n_iter):
        # Centroid of current selection (sign-flipped to align)
        ref = flats[sel_idx[0]]
        accum = ref.copy()
        for k in sel_idx[1:]:
            accum = accum + flats[k] * np.sign(np.dot(ref, flats[k]) + 1e-12)
        centroid = accum / len(sel_idx)
        cn = np.linalg.norm(centroid) or 1.0
        cent = np.array([abs(np.dot(f, centroid)) / max(np.linalg.norm(f) * cn, 1e-9)
                         for f in flats])
        # Combined score: bond × centroid_align
        combo = bond * (0.5 + 0.5 * cent)
        # Diversity-aware top-2 selection
        order = np.argsort(-combo)
        new_sel = [int(order[0])]
        for j in order[1:]:
            j = int(j)
            if all(abs(np.dot(flats[j], flats[k])) /
                   max(norms[j] * norms[k], 1e-9) < 0.85 for k in new_sel):
                new_sel.append(j)
                if len(new_sel) >= 2: break
        if not new_sel:
            new_sel = [int(order[0]), int(order[1] if len(order) > 1 else order[0])]
        elif len(new_sel) < 2:
            for j in order:
                if int(j) not in new_sel:
                    new_sel.append(int(j)); break
        sel_idx = new_sel
    # Top 3 with same logic
    final = [(combo[i], picks[i][1], picks[i][2]) for i in sel_idx]
    # Add a 3rd by bond_ov diversity
    used = set(sel_idx)
    for j in np.argsort(-bond):
        j = int(j)
        if j in used: continue
        final.append((bond[j], picks[j][1], picks[j][2]))
        if len(final) >= 3: break
    return final


def rk_rxn_overlap_filter(sd, alpha=0.5, min_rxn=0.0):
    """Filter IGs by min_rxn (rxn_overlap on the picked mode), then
    diversity-rank by bond_overlap."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn:
            continue
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                score *= max(0.0, 1.0 - alpha * cos_modes(cand[1], sel[1]))
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


# ─── ITERATION 5: combine the best ideas ──────────────────────────────────

def rk_combined_v1(sd, alpha=0.5, min_rxn=0.10, min_freq=50):
    """Filter: rxn_ov ≥ min_rxn AND |freq| ≥ min_freq.
    Score: bond_ov × (1 + 0.3 × rxn_ov).
    Top-k: diversity penalty."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn:
            continue
        if abs(picked['freq']) < min_freq:
            continue
        score = picked.get('bond_overlap', 0) * (1 + 0.3 * picked.get('rxn_overlap', 0))
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                score *= max(0.0, 1.0 - alpha * cos_modes(cand[1], sel[1]))
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


def rk_combined_v2(sd, alpha=0.5, min_rxn=0.10, w_rxn=0.5, w_core=0.3):
    """Score: bond_ov × (1 + w_rxn × rxn_ov + w_core × core_frac).
    Filter: rxn_ov ≥ min_rxn. Diversity penalty."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn:
            continue
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r + w_core * c)
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            score = cand[0]
            for sel in selected:
                score *= max(0.0, 1.0 - alpha * cos_modes(cand[1], sel[1]))
            if score > best_score:
                best_score = score; best = cand
        selected.append(best); remaining.remove(best)
    return selected


def rk_combined_iterative(sd, alpha=0.5, min_rxn=0.10, n_iter=3):
    """rxn_ov filter + iterative centroid refinement."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn:
            continue
        candidates.append((picked.get('bond_overlap', 0), picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    flats = [np.asarray(p[1]['disp']).reshape(-1) for p in candidates]
    norms = np.array([np.linalg.norm(f) for f in flats])
    bond = np.array([p[0] for p in candidates])
    if len(candidates) < 2:
        return candidates
    sel_idx = list(np.argsort(-bond)[:2])
    for _ in range(n_iter):
        ref = flats[sel_idx[0]]
        accum = ref.copy()
        for k in sel_idx[1:]:
            accum = accum + flats[k] * np.sign(np.dot(ref, flats[k]) + 1e-12)
        centroid = accum / len(sel_idx)
        cn = np.linalg.norm(centroid) or 1.0
        cent = np.array([abs(np.dot(f, centroid)) / max(np.linalg.norm(f) * cn, 1e-9)
                         for f in flats])
        combo = bond * (0.5 + 0.5 * cent)
        order = np.argsort(-combo)
        new_sel = [int(order[0])]
        for j in order[1:]:
            j = int(j)
            if all(abs(np.dot(flats[j], flats[k])) /
                   max(norms[j] * norms[k], 1e-9) < 0.85 for k in new_sel):
                new_sel.append(j)
                if len(new_sel) >= 2: break
        if len(new_sel) < 2:
            for j in order:
                if int(j) not in new_sel:
                    new_sel.append(int(j)); break
        sel_idx = new_sel
    final = [(combo[i], candidates[i][1], candidates[i][2]) for i in sel_idx]
    used = set(sel_idx)
    for j in np.argsort(-bond):
        j = int(j)
        if j in used: continue
        final.append((bond[j], candidates[j][1], candidates[j][2]))
        if len(final) >= 3: break
    return final


# ─── ITERATION 6: aggressive combinations ─────────────────────────────────

def rk_aggressive_v1(sd, alpha=0.5, min_rxn=0.10, w_rxn=1.0, w_core=0.2):
    """Strongest combination: filter, weighted score, diversity, mass-aware."""
    candidates = []
    elements = sd['elR']
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn:
            continue
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r) * (1 + w_core * c)
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
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


def rk_two_stage(sd, alpha=0.5, min_rxn=0.10):
    """Stage 1: pick top-1 strict (highest score among rxn_ov≥0.10).
       Stage 2: pick top-2 by maxim score with hard exclusion of cos>0.85
                relative to top-1."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn:
            continue
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        score = b * (1 + r)
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    candidates.sort(key=lambda t: -t[0])
    selected = [candidates[0]]
    for cand in candidates[1:]:
        if cos_modes(cand[1], selected[0][1]) < 0.85:
            selected.append(cand)
            if len(selected) >= 3: break
    while len(selected) < 3 and len(selected) < len(candidates):
        for c in candidates:
            if c not in selected:
                selected.append(c); break
    return selected


def rk_anchor_then_complement(sd, min_rxn=0.10):
    """Top-1: best by bond × rxn (strong anchor).
       Top-2: among IGs with bond_ov ≥ 0.7×top1's bond_ov, pick the one
              MOST DIFFERENT from top-1 (1-cos)."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        if r < min_rxn: continue
        candidates.append((b * (1 + r), b, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, 0.5)
    candidates.sort(key=lambda t: -t[0])
    top1 = candidates[0]
    top1_bond = top1[1]
    bond_threshold = 0.7 * top1_bond
    diverse_pool = [c for c in candidates[1:] if c[1] >= bond_threshold]
    if not diverse_pool:
        diverse_pool = candidates[1:]
    diverse_pool.sort(key=lambda c: -1 * (1.0 - cos_modes(c[2], top1[2])))
    selected = [(top1[0], top1[2], top1[3])]
    for c in diverse_pool[:2]:
        selected.append((c[0], c[2], c[3]))
        if len(selected) >= 3: break
    while len(selected) < 3:
        for c in candidates:
            cand_t = (c[0], c[2], c[3])
            if cand_t not in selected:
                selected.append(cand_t); break
        else: break
    return selected


def rk_ensemble_vote(sd, min_rxn=0.10):
    """Each IG gets 3 votes from 3 rankers; sum votes; pick highest with
    diversity penalty."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn:
            continue
        candidates.append((picked, ts['label']))
    if not candidates: return rk_bond_diverse(sd, 0.5)
    n = len(candidates)
    bonds  = np.array([c[0].get('bond_overlap', 0) for c in candidates])
    rxns   = np.array([c[0].get('rxn_overlap', 0) for c in candidates])
    cores  = np.array([c[0].get('core_fraction', 0) for c in candidates])
    flats  = [np.asarray(c[0]['disp']).reshape(-1) for c in candidates]
    norms  = np.array([np.linalg.norm(f) for f in flats])
    # peer score
    peer = np.zeros(n)
    for i in range(n):
        s = 0
        for j in range(n):
            if i == j: continue
            s += abs(flats[i] @ flats[j]) / max(norms[i] * norms[j], 1e-9)
        peer[i] = s / max(n - 1, 1)
    # vote: top-3 in each metric
    votes = np.zeros(n)
    for arr in (bonds, rxns, cores, peer):
        for i in arr.argsort()[-5:]:  # top-5 in each
            votes[i] += 1
    score = votes + bonds  # tiebreak by bond_ov
    order = np.argsort(-score)
    selected_idx = []
    for i in order:
        i = int(i)
        ok = True
        for j in selected_idx:
            if cos_modes(candidates[i][0], candidates[j][0]) >= 0.85:
                ok = False; break
        if ok:
            selected_idx.append(i)
        if len(selected_idx) >= 3: break
    while len(selected_idx) < 3:
        for i in order:
            i = int(i)
            if i not in selected_idx:
                selected_idx.append(i); break
        else: break
    return [(score[i], candidates[i][0], candidates[i][1]) for i in selected_idx]


# ─── ITERATION 7: pool modes across IGs (within-IG mode selection) ───────
# Insight: aggressive_v1 collapses each IG to its bond_overlap-best mode.
# But the oracle ranges over ALL modes of ALL IGs. If the right basin is
# captured by an IG but the highest-bond_overlap mode in it isn't the right
# one, aggressive_v1 misses. Pool all imag modes from all IGs as candidates,
# then global diversity selection.

def rk_modepool(sd, alpha=0.7, min_rxn=0.10, w_rxn=1.0, w_core=0.2,
                topm_per_ig=3, one_per_ig=True):
    """Pool top-M imag modes per IG. Score same as aggressive_v1.
    one_per_ig: enforce at most one mode per IG in the top-3 (avoids
    duplicating obviously-similar modes)."""
    candidates = []
    elements = sd['elR']
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        # Sort imag modes by bond_overlap, take top-M
        imag_sorted = sorted(imag, key=lambda m: -m.get('bond_overlap', 0))
        for picked in imag_sorted[:topm_per_ig]:
            r = picked.get('rxn_overlap', 0)
            if r < min_rxn: continue
            b = picked.get('bond_overlap', 0)
            c = picked.get('core_fraction', 0)
            score = b * (1 + w_rxn * r) * (1 + w_core * c)
            candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    selected = []
    selected_labels = set()
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            if one_per_ig and cand[2] in selected_labels: continue
            score = cand[0]
            for sel in selected:
                mw = mass_weighted_cos(cand[1]['disp'], sel[1]['disp'], elements)
                score *= max(0.0, 1.0 - alpha * mw)
            if score > best_score:
                best_score = score; best = cand
        if best is None: break
        selected.append(best)
        selected_labels.add(best[2])
        remaining.remove(best)
    return selected


def rk_modepool_within_score(sd, alpha=0.7, min_rxn=0.10, w_rxn=1.0,
                              w_core=0.2, topm_per_ig=3, one_per_ig=True):
    """Same as modepool, but the within-IG mode is selected by the
    weighted score (bond × (1+rxn) × (1+core)), not by bond_overlap.
    Then ranking + diversity is on the same score."""
    candidates = []
    elements = sd['elR']
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        scored_modes = []
        for m in imag:
            r = m.get('rxn_overlap', 0)
            if r < min_rxn: continue
            b = m.get('bond_overlap', 0)
            c = m.get('core_fraction', 0)
            sc = b * (1 + w_rxn * r) * (1 + w_core * c)
            scored_modes.append((sc, m))
        scored_modes.sort(key=lambda t: -t[0])
        for sc, m in scored_modes[:topm_per_ig]:
            candidates.append((sc, m, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    selected = []
    selected_labels = set()
    remaining = list(candidates)
    while remaining and len(selected) < 3:
        best = None; best_score = -1e18
        for cand in remaining:
            if one_per_ig and cand[2] in selected_labels: continue
            score = cand[0]
            for sel in selected:
                mw = mass_weighted_cos(cand[1]['disp'], sel[1]['disp'], elements)
                score *= max(0.0, 1.0 - alpha * mw)
            if score > best_score:
                best_score = score; best = cand
        if best is None: break
        selected.append(best)
        selected_labels.add(best[2])
        remaining.remove(best)
    return selected


# ─── ITERATION 8: better within-IG mode selection ────────────────────────
# Maybe the issue is purely "which mode within each IG is best by bond_ov".
# Try: pick within IG by bond * (1 + rxn) or other combined score.

def rk_within_combo(sd, alpha=0.7, min_rxn=0.10, w_rxn=1.0, w_core=0.2,
                    within_w_rxn=1.0, within_w_core=0.0):
    """Within-IG mode pick by bond × (1 + within_w_rxn × rxn) ×
    (1 + within_w_core × core). IG ranking uses the global score."""
    candidates = []
    elements = sd['elR']
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        def within_score(m):
            return (m.get('bond_overlap', 0)
                    * (1 + within_w_rxn * m.get('rxn_overlap', 0))
                    * (1 + within_w_core * m.get('core_fraction', 0)))
        picked = max(imag, key=within_score)
        r = picked.get('rxn_overlap', 0)
        if r < min_rxn: continue
        b = picked.get('bond_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r) * (1 + w_core * c)
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
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


# ─── ITERATION 9: dwbo within-IG mode selection ──────────────────────────
# Compute dwbo_overlap for each mode and use it as part of within-IG pick.

def _dwbo_per_mode(sd, ts):
    """Compute dwbo_overlap for each mode in a TS. Cached lazily on the ts dict."""
    if '_dwbo_cache' in ts:
        return ts['_dwbo_cache']
    from eval_dwbo_overlap import dwbo_reaction_vector, dwbo_overlap_per_mode
    ts_xyz = np.asarray(ts['xyz_coords'])
    V = dwbo_reaction_vector(ts_xyz, sd['wboR'], sd['wboP'], sd['mapping_RP'])
    modes_R = np.asarray([mm['disp'] for mm in ts['modes']])
    out = dwbo_overlap_per_mode(modes_R, V)
    ts['_dwbo_cache'] = out
    return out


def rk_within_bond_x_dwbo(sd, alpha=0.7, min_rxn=0.10, w_rxn=1.0, w_core=0.2):
    """Within-IG: pick mode by bond × dwbo. Global rank: aggressive_v1 score."""
    candidates = []
    elements = sd['elR']
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        dwbo = _dwbo_per_mode(sd, ts)
        # Find imag modes' indices
        imag_idx = [i for i, m in enumerate(ts['modes']) if m['freq'] < 0]
        scored = [(ts['modes'][i].get('bond_overlap', 0) * float(dwbo[i]),
                   ts['modes'][i]) for i in imag_idx]
        if not scored: continue
        picked = max(scored, key=lambda x: x[0])[1]
        r = picked.get('rxn_overlap', 0)
        if r < min_rxn: continue
        b = picked.get('bond_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r) * (1 + w_core * c)
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
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


# ─── ITERATION 10: cluster-then-rank ──────────────────────────────────────
# Pool all imag modes across IGs. Cluster by mass-weighted cosine. For each
# cluster: score = max(score) × log(1+size). Pick top representatives.

def rk_cluster_pool(sd, alpha=0.7, min_rxn=0.10, w_rxn=1.0, w_core=0.2,
                    cluster_thr=0.7):
    candidates = []
    elements = sd['elR']
    for ts in sd['igs']:
        imag = imag_modes(ts)
        for m in imag:
            r = m.get('rxn_overlap', 0)
            if r < min_rxn: continue
            b = m.get('bond_overlap', 0)
            c = m.get('core_fraction', 0)
            sc = b * (1 + w_rxn * r) * (1 + w_core * c)
            candidates.append((sc, m, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    n = len(candidates)
    # Mass-weighted similarity matrix (only need for clustering — sparse path)
    order = sorted(range(n), key=lambda i: -candidates[i][0])
    cluster_id = [-1] * n
    cid = 0
    for i in order:
        if cluster_id[i] != -1: continue
        cluster_id[i] = cid
        for j in range(n):
            if cluster_id[j] == -1:
                mw = mass_weighted_cos(candidates[i][1]['disp'],
                                       candidates[j][1]['disp'], elements)
                if mw >= cluster_thr:
                    cluster_id[j] = cid
        cid += 1
    # Score each cluster: max(score) × log(1+size)  AND store representative
    clusters = {}
    for i, c in enumerate(cluster_id):
        clusters.setdefault(c, []).append(i)
    cluster_score = []
    for c, members in clusters.items():
        # representative: highest-scoring member with unique IG label
        members.sort(key=lambda i: -candidates[i][0])
        rep = members[0]
        size = len(set(candidates[i][2] for i in members))  # unique IG count
        sc = candidates[rep][0] * (1 + 0.3 * np.log1p(size))
        cluster_score.append((sc, rep, members))
    cluster_score.sort(key=lambda t: -t[0])
    # Pick reps from top clusters, with diversity penalty across reps
    selected = []
    remaining = list(cluster_score)
    while remaining and len(selected) < 3:
        best = None; best_sc = -1e18
        for cs, rep, members in remaining:
            score = cs
            for sel in selected:
                mw = mass_weighted_cos(candidates[rep][1]['disp'],
                                       candidates[sel[1]]['disp'] if False
                                       else candidates[sel[3]][1]['disp'],
                                       elements)
                score *= max(0.0, 1.0 - alpha * mw)
            if score > best_sc:
                best_sc = score; best = (cs, rep, members, rep)
        if best is None: break
        cs, rep, members, _ = best
        selected.append((cs, candidates[rep][1], candidates[rep][2], rep))
        remaining = [t for t in remaining if t[1] != rep]
    return [(s[0], s[1], s[2]) for s in selected]


# ─── ITERATION 12: core_fraction-dominant scores ─────────────────────────
# Diagnostic showed core_fraction alone has 30.8% rank-of-best-IG≤2
# vs 26.4% for bond_overlap. So core should be weighted more.

def _diversity_select(candidates, alpha, elements):
    selected = []
    remaining = list(candidates)
    while remaining and len(selected) < 3:
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


def rk_pure_core(sd, alpha=0.7, min_rxn=0.10):
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        candidates.append((picked.get('core_fraction', 0), picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    return _diversity_select(candidates, alpha, sd['elR'])


def rk_core_x_bond(sd, alpha=0.7, min_rxn=0.10):
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        sc = picked.get('bond_overlap', 0) * picked.get('core_fraction', 0)
        candidates.append((sc, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    return _diversity_select(candidates, alpha, sd['elR'])


def rk_sqrt_core_bond(sd, alpha=0.7, min_rxn=0.10):
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        sc = float(np.sqrt(max(0, picked.get('bond_overlap', 0)
                               * picked.get('core_fraction', 0))))
        candidates.append((sc, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    return _diversity_select(candidates, alpha, sd['elR'])


def rk_core_with_rxn(sd, alpha=0.7, min_rxn=0.10):
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        sc = picked.get('core_fraction', 0) * (1 + picked.get('rxn_overlap', 0))
        candidates.append((sc, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    return _diversity_select(candidates, alpha, sd['elR'])


# ─── ITERATION 13: Borda rank fusion + structural prior ──────────────────
# Idea: rank each IG by multiple signals, use sum-of-ranks (Borda).
# Add a "single-imaginary" prior — a clean TS has exactly one imag mode.

def rk_borda(sd, alpha=0.7, min_rxn=0.10, weights=None):
    """weights: dict of {feature: weight}. None = uniform."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        n_imag = len(imag)
        candidates.append({
            'b': picked.get('bond_overlap', 0),
            'r': picked.get('rxn_overlap', 0),
            'c': picked.get('core_fraction', 0),
            'n_imag': n_imag,
            'mode': picked, 'label': ts['label'],
        })
    if len(candidates) < 2:
        return rk_bond_diverse(sd, alpha)
    n = len(candidates)
    if weights is None:
        weights = {'b': 1.0, 'c': 1.0, 'r': 0.5, 'single_imag': 0.5}
    # Compute ranks per feature (1-based, lower = better)
    def ranks_by(key, descending=True):
        vals = [(i, c[key]) for i, c in enumerate(candidates)]
        vals.sort(key=lambda t: -t[1] if descending else t[1])
        out = [0] * n
        for r, (i, _) in enumerate(vals, 1):
            out[i] = r
        return out
    r_b = ranks_by('b')
    r_c = ranks_by('c')
    r_r = ranks_by('r')
    r_imag = ranks_by('n_imag', descending=False)  # fewer imag = better
    # Higher score = better; compute -weighted sum of ranks
    scores = []
    for i in range(n):
        s = -(weights['b'] * r_b[i] + weights['c'] * r_c[i]
              + weights['r'] * r_r[i] + weights['single_imag'] * r_imag[i])
        scores.append((s, candidates[i]['mode'], candidates[i]['label']))
    return _diversity_select(scores, alpha, sd['elR'])


def rk_single_imag_filter(sd, alpha=0.7, min_rxn=0.10, max_imag=2,
                           w_rxn=1.0, w_core=0.2):
    """Like aggressive_v1 but FILTER OUT IGs with > max_imag imaginary modes
    (those are not clean transition states)."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        if len(imag) > max_imag: continue  # not a clean TS
        picked = max(imag, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r) * (1 + w_core * c)
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    return _diversity_select(candidates, alpha, sd['elR'])


def rk_combined_clean(sd, alpha=0.7, min_rxn=0.10, max_imag=3):
    """aggressive_v1 + single-imag filter + Borda-style core/bond fusion."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        if len(imag) > max_imag: continue
        picked = max(imag, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        c = picked.get('core_fraction', 0)
        n_im = len(imag)
        score = b * (1 + r) * (1 + 0.5 * c) / max(n_im, 1) ** 0.3
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    return _diversity_select(candidates, alpha, sd['elR'])


def rk_clean_v2(sd, alpha=0.7, min_rxn=0.10, max_imag=2,
                 w_rxn=1.0, w_core=0.5, imag_pen=0.3):
    """combined_clean tuned + soft penalty for n_imag (rather than hard cutoff)."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        n_imag = len(imag)
        if n_imag > max_imag: continue
        picked = max(imag, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r) * (1 + w_core * c) / max(n_imag, 1) ** imag_pen
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    return _diversity_select(candidates, alpha, sd['elR'])


def rk_clean_max3(sd, alpha=0.7, min_rxn=0.10, w_rxn=1.0, w_core=0.5):
    """Soft max_imag (penalty), no hard cutoff."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        picked = max(imag, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        n_imag = len(imag)
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r) * (1 + w_core * c) / (1 + 0.3 * (n_imag - 1))
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, alpha)
    return _diversity_select(candidates, alpha, sd['elR'])


def rk_adaptive_alpha(sd, min_rxn=0.10, w_rxn=1.0, w_core=0.2):
    """Top-1 score is high (>=0.5 bond) → use alpha=0.8 for diversity.
       Top-1 score is low (<0.3 bond) → use alpha=0.3 (trust ranking more)."""
    candidates = []
    elements = sd['elR']
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r) * (1 + w_core * c)
        candidates.append((score, picked, ts['label']))
    if not candidates:
        return rk_bond_diverse(sd, 0.5)
    candidates.sort(key=lambda t: -t[0])
    top_b = candidates[0][1].get('bond_overlap', 0)
    if top_b >= 0.5:    alpha = 0.8
    elif top_b >= 0.3:  alpha = 0.6
    else:               alpha = 0.3
    return _diversity_select(candidates, alpha, elements)


def rk_normfuse(sd, alpha=0.7, min_rxn=0.10, w_b=1.0, w_c=0.5, w_r=0.3,
                 imag_pen=0.3, max_imag=3):
    """Normalize each signal to [0,1] within the step's IG pool, then sum."""
    rows = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        if not imag: continue
        if len(imag) > max_imag: continue
        picked = max(imag, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        rows.append({
            'b': picked.get('bond_overlap', 0),
            'r': picked.get('rxn_overlap', 0),
            'c': picked.get('core_fraction', 0),
            'n_imag': len(imag),
            'mode': picked, 'label': ts['label'],
        })
    if len(rows) < 2:
        return rk_bond_diverse(sd, alpha)
    bs = np.array([x['b'] for x in rows])
    cs = np.array([x['c'] for x in rows])
    rs = np.array([x['r'] for x in rows])
    bn = (bs - bs.min()) / max(bs.max() - bs.min(), 1e-9)
    cn = (cs - cs.min()) / max(cs.max() - cs.min(), 1e-9)
    rn = (rs - rs.min()) / max(rs.max() - rs.min(), 1e-9)
    candidates = []
    for i, x in enumerate(rows):
        s = (w_b * bn[i] + w_c * cn[i] + w_r * rn[i]) \
            / (1 + imag_pen * (x['n_imag'] - 1))
        candidates.append((float(s), x['mode'], x['label']))
    return _diversity_select(candidates, alpha, sd['elR'])


def rk_b_plus_c_rank(sd, alpha=0.7, min_rxn=0.10):
    """Score = bond_rank_pct + core_rank_pct (additive Borda on b+c)."""
    candidates = []
    for ts in sd['igs']:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        candidates.append({
            'b': picked.get('bond_overlap', 0),
            'c': picked.get('core_fraction', 0),
            'mode': picked, 'label': ts['label'],
        })
    if len(candidates) < 2:
        return rk_bond_diverse(sd, alpha)
    n = len(candidates)
    # Convert to rank percentiles
    bs = np.array([c['b'] for c in candidates])
    cs = np.array([c['c'] for c in candidates])
    rb = (-bs).argsort().argsort() + 1  # 1=best
    rc = (-cs).argsort().argsort() + 1
    scores = []
    for i in range(n):
        s = -(rb[i] + rc[i])
        scores.append((s, candidates[i]['mode'], candidates[i]['label']))
    return _diversity_select(scores, alpha, sd['elR'])


# ─── load all step data ───────────────────────────────────────────────────

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

    candidates = [
        ('BASELINE bond_overlap',           rk_bond_overlap),
        ('aggressive_v1 alpha=0.7 (prev best)',
                                            lambda s: rk_aggressive_v1(s, 0.7, 0.10, 1.0, 0.2)),
        # Iteration 12 — push core_fraction weight (since it's the best
        # rank=1-of-best-IG signal: 30.8% vs bond_overlap 26.4%)
        ('aggressive w_core=0.5',           lambda s: rk_aggressive_v1(s, 0.7, 0.10, 1.0, 0.5)),
        ('aggressive w_core=1.0',           lambda s: rk_aggressive_v1(s, 0.7, 0.10, 1.0, 1.0)),
        ('aggressive w_core=2.0',           lambda s: rk_aggressive_v1(s, 0.7, 0.10, 1.0, 2.0)),
        ('aggressive w_core=3.0',           lambda s: rk_aggressive_v1(s, 0.7, 0.10, 1.0, 3.0)),
        ('aggressive w_core=5.0 w_rxn=0.5', lambda s: rk_aggressive_v1(s, 0.7, 0.10, 0.5, 5.0)),
        ('aggressive w_core=10 w_rxn=0',    lambda s: rk_aggressive_v1(s, 0.7, 0.10, 0.0, 10.0)),
        # Pure core_fraction with diversity
        ('pure core_fraction + diversity',  lambda s: rk_pure_core(s, 0.7, 0.10)),
        ('core × bond filter',              lambda s: rk_core_x_bond(s, 0.7, 0.10)),
        ('sqrt(core × bond)',               lambda s: rk_sqrt_core_bond(s, 0.7, 0.10)),
        ('core_fraction × (1+r)',           lambda s: rk_core_with_rxn(s, 0.7, 0.10)),
        # Iter 13 — fusion + structural prior
        ('Borda b+c+r (uniform)',           lambda s: rk_borda(s, 0.7, 0.10, {'b':1,'c':1,'r':0.5,'single_imag':0})),
        ('Borda b+c+single_imag',           lambda s: rk_borda(s, 0.7, 0.10, {'b':1,'c':1,'r':0,'single_imag':0.7})),
        ('Borda all uniform',               lambda s: rk_borda(s, 0.7, 0.10, None)),
        ('aggressive + max_imag=2',         lambda s: rk_single_imag_filter(s, 0.7, 0.10, 2, 1.0, 0.2)),
        ('aggressive + max_imag=1',         lambda s: rk_single_imag_filter(s, 0.7, 0.10, 1, 1.0, 0.2)),
        ('combined_clean max_imag=3',       lambda s: rk_combined_clean(s, 0.7, 0.10, 3)),
        ('combined_clean max_imag=2',       lambda s: rk_combined_clean(s, 0.7, 0.10, 2)),
        ('rank-pct b+c',                    lambda s: rk_b_plus_c_rank(s, 0.7, 0.10)),
        # Iter 14 — clean_v2 sweep
        ('clean_v2 mi2 wc0.5 ip0.3',        lambda s: rk_clean_v2(s, 0.7, 0.10, 2, 1.0, 0.5, 0.3)),
        ('clean_v2 mi2 wc0.2 ip0.3',        lambda s: rk_clean_v2(s, 0.7, 0.10, 2, 1.0, 0.2, 0.3)),
        ('clean_v2 mi2 wc1.0 ip0.5',        lambda s: rk_clean_v2(s, 0.7, 0.10, 2, 1.0, 1.0, 0.5)),
        ('clean_v2 mi3 wc0.5 ip0.3',        lambda s: rk_clean_v2(s, 0.7, 0.10, 3, 1.0, 0.5, 0.3)),
        ('clean_max3 wc0.2',                lambda s: rk_clean_max3(s, 0.7, 0.10, 1.0, 0.2)),
        ('clean_max3 wc0.5',                lambda s: rk_clean_max3(s, 0.7, 0.10, 1.0, 0.5)),
        ('clean_max3 wc1.0',                lambda s: rk_clean_max3(s, 0.7, 0.10, 1.0, 1.0)),
        ('adaptive_alpha',                  lambda s: rk_adaptive_alpha(s, 0.10, 1.0, 0.2)),
        # Iter 15 — normalized fusion
        ('normfuse b1 c0.5 r0.3',           lambda s: rk_normfuse(s, 0.7, 0.10, 1.0, 0.5, 0.3, 0.3, 3)),
        ('normfuse b1 c1.0 r0.5',           lambda s: rk_normfuse(s, 0.7, 0.10, 1.0, 1.0, 0.5, 0.3, 3)),
        ('normfuse b1 c0.5 r0 mi3',         lambda s: rk_normfuse(s, 0.7, 0.10, 1.0, 0.5, 0.0, 0.3, 3)),
        ('normfuse b1 c1 r1 mi2',           lambda s: rk_normfuse(s, 0.7, 0.10, 1.0, 1.0, 1.0, 0.3, 2)),
        ('normfuse b2 c1 r0.5 mi3',         lambda s: rk_normfuse(s, 0.7, 0.10, 2.0, 1.0, 0.5, 0.3, 3)),
    ]
    print(f"{'method':40s}  k=2 results")
    print('=' * 100)
    for name, fn in candidates:
        res = evaluate_topk(fn, scope)
        print(summarise(name, res))


if __name__ == '__main__':
    main()

"""
Shared verifier (ranker) and small utilities used by viewers and
analytic scripts.

Two ranker variants are exported:

  * `rk_aggressive_v1` — first generation (rxn>=0.10 filter +
                          bond*(1+r)*(1+0.2c) score + diversity)
  * `rk_clean_v2`      — current verifier (max_imag<=2 filter +
                          bond*(1+r)*(1+0.2c)/n_imag^0.3 score +
                          mass-weighted-cosine diversity, alpha=0.7)

Plus small helpers: cos_sim, mass_weighted_cos, imag_modes,
load_step (per-step HTML payload reader), ATOMIC_MASS table.
"""
from __future__ import annotations
import json
import re
from pathlib import Path

import numpy as np


# ─── Atomic masses (used for mass-weighted cosine in diversity) ──────────
ATOMIC_MASS = {
    'H':1.008,'He':4.003,'Li':6.94,'Be':9.012,'B':10.81,'C':12.011,
    'N':14.007,'O':15.999,'F':18.998,'Ne':20.180,'Na':22.99,'Mg':24.305,
    'Al':26.982,'Si':28.085,'P':30.974,'S':32.06,'Cl':35.45,'Ar':39.948,
    'K':39.098,'Ca':40.08,'Sc':44.956,'Ti':47.867,'V':50.942,'Cr':51.996,
    'Mn':54.938,'Fe':55.845,'Co':58.933,'Ni':58.693,'Cu':63.546,'Zn':65.38,
    'Ga':69.723,'Ge':72.63,'As':74.922,'Se':78.971,'Br':79.904,'Kr':83.798,
    'Rb':85.468,'Sr':87.62,'Y':88.906,'Zr':91.224,'Nb':92.906,'Mo':95.95,
    'Tc':98.0,'Ru':101.07,'Rh':102.906,'Pd':106.42,'Ag':107.868,'Cd':112.414,
    'In':114.818,'Sn':118.710,'Sb':121.760,'Te':127.60,'I':126.904,'Xe':131.293,
    'Cs':132.905,'Ba':137.327,'La':138.906,'Hf':178.49,'Ta':180.948,'W':183.84,
    'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.084,'Au':196.967,'Hg':200.592,
    'Tl':204.383,'Pb':207.2,'Bi':208.980,
}


# ─── Tiny helpers ─────────────────────────────────────────────────────────

def cos_sim(a, b):
    """Sign-blind cosine similarity in flattened-3N space."""
    a = np.asarray(a).reshape(-1); b = np.asarray(b).reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9: return 0.0
    return abs(float(a @ b)) / (na * nb)


def mass_weighted_cos(disp_a, disp_b, elements):
    a = np.asarray(disp_a, dtype=float)
    b = np.asarray(disp_b, dtype=float)
    masses = np.array([ATOMIC_MASS.get(e, 12.0) for e in elements])
    sqm = np.sqrt(masses)[:, None]
    aw = (a * sqm).reshape(-1)
    bw = (b * sqm).reshape(-1)
    na = float(np.linalg.norm(aw)); nb = float(np.linalg.norm(bw))
    if na < 1e-9 or nb < 1e-9: return 0.0
    return abs(float(aw @ bw)) / (na * nb)


def imag_modes(ts):
    return [m for m in ts['modes'] if m['freq'] < 0]


def best_bond_overlap_imag(ts):
    """Default within-IG mode pick: highest bond_overlap among imag modes
    (or among all modes if no imag exists)."""
    imag = imag_modes(ts)
    key = lambda m: m.get('bond_overlap', m.get('rxn_overlap',
                                                m.get('core_fraction', 0)))
    if imag:
        return max(imag, key=key)
    return max(ts['modes'], key=key)


# ─── Per-step HTML payload reader (used by analytics) ────────────────────

def load_step_payload(html_path):
    """Parse the embedded JSON payload from a per-step viewer HTML."""
    text = Path(html_path).read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m:
        raise ValueError(f"could not find DATA= in {html_path}")
    return json.loads(m.group(1))


# ─── Diversity selection (greedy mass-weighted-cosine penalty) ────────────

def _diversity_select(candidates, alpha, elements, k=3):
    """`candidates`: list of (score, picked_mode_dict, ig_label).
    Greedy: each step pick the highest-score candidate after applying
    multiplicative penalty (1 - alpha * mass_cos) for similarity to
    each already-selected candidate."""
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


# ─── Rankers ──────────────────────────────────────────────────────────────

def rk_bond_overlap(igs, elements=None, k=3):
    """Baseline: rank IGs by bond_overlap of their best imag mode."""
    picks = []
    for ts in igs:
        imag = imag_modes(ts)
        modes = imag if imag else ts['modes']
        if not modes: continue
        picked = max(modes, key=lambda m: m.get('bond_overlap', 0))
        picks.append((picked.get('bond_overlap', 0), picked, ts['label']))
    picks.sort(key=lambda t: -t[0])
    return [(t, m) for _, m, t in picks[:k]]  # caller wants (ts, mode); legacy


def rk_aggressive_v1(igs, elements,
                     alpha=0.7, min_rxn=0.10, w_rxn=1.0, w_core=0.2, k=3):
    """First generation verifier."""
    cands = []
    for ts in igs:
        imag = imag_modes(ts)
        if not imag: continue
        picked = max(imag, key=lambda m: m.get('bond_overlap', 0))
        if picked.get('rxn_overlap', 0) < min_rxn: continue
        b = picked.get('bond_overlap', 0)
        r = picked.get('rxn_overlap', 0)
        c = picked.get('core_fraction', 0)
        score = b * (1 + w_rxn * r) * (1 + w_core * c)
        cands.append((score, picked, ts['label']))
    if not cands:
        return _bond_fallback(igs, k)
    selected = _diversity_select(cands, alpha, elements, k)
    return [(_find_ts(igs, lbl), m) for _, m, lbl in selected]


def rk_clean_v2(igs, elements,
                alpha=0.7, min_rxn=0.10, max_imag=2,
                w_rxn=1.0, w_core=0.2, imag_pen=0.3, k=3):
    """Current verifier (clean_v2 mi2 wc0.2 ip0.3, see improved_ranker.md).

    Filter: 1 <= n_imag <= max_imag AND best-imag-mode rxn_overlap >= min_rxn.
    Score: bond_overlap * (1 + w_rxn*rxn) * (1 + w_core*core)
            / n_imag^imag_pen.
    Selection: greedy with mass-weighted-cosine diversity penalty.

    Falls back to bond_overlap top-k (no filter, no diversity) if no IG
    passes the filter (rare)."""
    cands = []
    for ts in igs:
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
        cands.append((score, picked, ts['label']))
    if not cands:
        return _bond_fallback(igs, k)
    selected = _diversity_select(cands, alpha, elements, k)
    return [(_find_ts(igs, lbl), m) for _, m, lbl in selected]


def _bond_fallback(igs, k):
    scored = []
    for ts in igs:
        picked = best_bond_overlap_imag(ts)
        scored.append((picked.get('bond_overlap', 0), ts, picked))
    scored.sort(key=lambda t: -t[0])
    return [(t, m) for _, t, m in scored[:k]]


def _find_ts(igs, label):
    for ts in igs:
        if ts['label'] == label: return ts
    return None

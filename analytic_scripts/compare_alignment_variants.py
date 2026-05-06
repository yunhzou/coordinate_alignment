"""
Compare three alignment-score variants against human expert judgement.

Variants:
  1. cartesian      — current metric:  |d_a · d_b| / (||d_a|| ||d_b||)
  2. mass_weighted  — same but with sqrt(mass) per-atom weighting:
                       d̃_i = sqrt(m_i) · d_i ; cosine in d̃-space
  3. bond_stretch   — chemistry-aware: reduce each mode to a vector of
                       per-bond stretch coefficients
                         s_b(d) = (d_i - d_j) · û_ij     for each
                                                          broken/formed bond
                       Cosine in this (B+F)-dimensional reduced space.

Aggregation: per step, the alignment is the cumulative max over the
verifier's top-2 picks (matches the IG#1, IG#2 columns the human rated).

Reads:
  appendix_perparation/viewer/mode_viewer/<step>.html  (per-step payload)
  appendix_perparation/analtics/final_quality_measurement-humanversion (1).csv
                                                       (IG#1, IG#2, pass)

Output:
  out/mode_analysis/alignment_variants_vs_human.csv  (per-step values)
  prints headline agreement table
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import csv, json, re
from pathlib import Path
import numpy as np

from ranker import (rk_clean_v2, ATOMIC_MASS, cos_sim,
                    mass_weighted_cos, imag_modes)

VIEWER_DIR = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
HUMAN_CSV  = PROJECT_ROOT / 'appendix_perparation' / 'analtics' / 'final_quality_measurement-humanversion (1).csv'
OUT_CSV    = PROJECT_ROOT / 'out' / 'mode_analysis' / 'alignment_variants_vs_human.csv'


def load_step_payload(html_path):
    text = html_path.read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m: return None
    return json.loads(m.group(1))


def cartesian_align(d_a, d_b):
    return cos_sim(d_a, d_b)


def mass_weighted_align(d_a, d_b, elements):
    return mass_weighted_cos(d_a, d_b, elements)


def bond_stretch_vec(disp, ts_coords, bond_pairs):
    """Per-bond stretch coefficient s_b = (d_i - d_j) · u_ij for each
    (i, j) in bond_pairs. Returns a vector of length len(bond_pairs)."""
    disp = np.asarray(disp); coords = np.asarray(ts_coords)
    out = np.zeros(len(bond_pairs))
    for k, (i, j) in enumerate(bond_pairs):
        u = coords[j] - coords[i]
        n = float(np.linalg.norm(u))
        if n < 1e-9: continue
        u /= n
        out[k] = float((disp[i] - disp[j]) @ u)
    return out


def core_atoms_set(broken, formed):
    """Atoms touching any reactive bond — chemistry-aware focus set."""
    s = set()
    for (i, j) in broken: s.add(int(i)); s.add(int(j))
    for (i, j) in formed: s.add(int(i)); s.add(int(j))
    return sorted(s)


def core_cartesian_align(d_a, d_b, core):
    """Cartesian cosine restricted to core atoms only — discards
    spectator dilution but keeps full 3-vector richness."""
    if not core: return float('nan')
    a = np.asarray(d_a)[core].reshape(-1)
    b = np.asarray(d_b)[core].reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9: return float('nan')
    return abs(float(a @ b)) / (na * nb)


def core_mass_align(d_a, d_b, core, elements):
    """Mass-weighted cosine restricted to core atoms."""
    if not core: return float('nan')
    a = np.asarray(d_a)[core]; b = np.asarray(d_b)[core]
    masses = np.array([ATOMIC_MASS.get(elements[i], 12.0) for i in core])
    sqm = np.sqrt(masses)[:, None]
    aw = (a * sqm).reshape(-1); bw = (b * sqm).reshape(-1)
    na = float(np.linalg.norm(aw)); nb = float(np.linalg.norm(bw))
    if na < 1e-9 or nb < 1e-9: return float('nan')
    return abs(float(aw @ bw)) / (na * nb)


def bond_stretch_align(disp_a, disp_b, ts_coords_a, ts_coords_b,
                        broken, formed):
    """Cosine similarity in the reduced (B+F)-dimensional bond-stretch
    space. Each mode is reduced to one signed coefficient per reactive
    bond. Sign convention: positive means atoms moving apart.

    NB: the two TS coords differ between IG and GT, but the bond
    *indices* are R-frame — the same chemical bond in both. We use
    each side's own TS coords to compute its own unit vector."""
    bonds = list(broken) + list(formed)
    if not bonds: return float('nan')
    sa = bond_stretch_vec(disp_a, ts_coords_a, bonds)
    sb = bond_stretch_vec(disp_b, ts_coords_b, bonds)
    na = float(np.linalg.norm(sa)); nb = float(np.linalg.norm(sb))
    if na < 1e-9 or nb < 1e-9: return float('nan')
    return abs(float(sa @ sb)) / (na * nb)


def best_aligned_score(ts, gt_disp, gt_coords, broken, formed, score_kind, elements=None):
    """Best alignment across all modes of `ts` under the chosen metric.
    Returns the max alignment achievable from this IG."""
    best = -1.0
    core = core_atoms_set(broken, formed)
    for m in ts['modes']:
        d = np.asarray(m['disp'])
        if score_kind == 'cartesian':
            a = cartesian_align(d, gt_disp)
        elif score_kind == 'mass_weighted':
            a = mass_weighted_align(d, gt_disp, elements)
        elif score_kind == 'bond_stretch':
            a = bond_stretch_align(d, gt_disp,
                                   np.asarray(ts['xyz_coords']),
                                   gt_coords, broken, formed)
        elif score_kind == 'core_cartesian':
            a = core_cartesian_align(d, gt_disp, core)
        elif score_kind == 'core_mass':
            a = core_mass_align(d, gt_disp, core, elements)
        if a is None or (isinstance(a, float) and np.isnan(a)):
            a = -1.0
        if a > best: best = a
    return max(best, 0.0)


def main():
    # Read human ratings
    lines = HUMAN_CSV.read_text().splitlines()
    annotation = lines[0]
    header = lines[1].split(',')
    data_rows = [dict(zip(header, ln.split(','))) for ln in lines[2:] if ln.strip()]
    print(f"Loaded {len(data_rows)} human-rated rows")

    # Build canonical → step name lookup
    def canon(s): return re.sub(r'\s*\(.*?\)\s*$', '', s).strip()

    out_rows = []
    skipped = []
    for hr in data_rows:
        step_raw = hr['step']
        step = canon(step_raw)
        ig1 = hr.get('IG#1', ''); ig2 = hr.get('IG#2', '')
        pass_label = int(hr['pass']) if hr.get('pass') in ('0', '1') else None
        if pass_label is None:
            skipped.append((step, 'no_pass')); continue

        hp = VIEWER_DIR / f"{step}.html"
        if not hp.exists():
            skipped.append((step, 'no_html')); continue
        data = load_step_payload(hp)
        if data is None:
            skipped.append((step, 'parse')); continue
        gt = next((t for t in data['ts_list']
                   if t['label']=='groundtruth' and t.get('modes')), None)
        if gt is None:
            skipped.append((step, 'no_gt')); continue
        gt_disp = np.asarray(gt['modes'][gt['default_mode_idx']]['disp'])
        gt_coords = np.asarray(gt['xyz_coords'])
        elements = gt['xyz_elements']
        broken = data['broken_bonds']
        formed = data['formed_bonds_R']

        igs = [t for t in data['ts_list']
               if t['label']!='groundtruth' and t.get('modes')]
        if not igs:
            skipped.append((step, 'no_igs')); continue

        # Verifier picks (top-2 like the human rated)
        ranked = rk_clean_v2(igs, elements)[:2]
        if len(ranked) < 2:
            # pad from bond_overlap if needed
            already = {t['label'] for t, _ in ranked}
            for ts in igs:
                if ts['label'] in already: continue
                imag = imag_modes(ts)
                if not imag: continue
                pk = max(imag, key=lambda m: m.get('bond_overlap', 0))
                ranked.append((ts, pk))
                if len(ranked) >= 2: break
        if len(ranked) < 2:
            skipped.append((step, 'lt_2_picks')); continue

        # For each variant, compute pass@2 = max over the 2 IGs of
        # best-mode-in-IG alignment under that metric
        results = {}
        for kind in ('cartesian', 'mass_weighted', 'bond_stretch',
                     'core_cartesian', 'core_mass'):
            scores = []
            for ts, _ in ranked:
                s = best_aligned_score(ts, gt_disp, gt_coords, broken, formed,
                                       kind, elements=elements)
                scores.append(s)
            results[f'{kind}_pass2'] = round(max(scores), 6) if scores else 0.0
            results[f'{kind}_top1']  = round(scores[0], 6)
            results[f'{kind}_top2']  = round(scores[1], 6)

        out_rows.append(dict(
            step=step_raw, IG1=ig1, IG2=ig2, pass_label=pass_label,
            n_broken=len(broken), n_formed=len(formed),
            **results,
        ))

    print(f"Computed for {len(out_rows)} steps; skipped {len(skipped)}")
    if skipped: print(f"  skips: {skipped[:5]}{'...' if len(skipped)>5 else ''}")

    # Save full table
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader(); w.writerows(out_rows)
    print(f"\nWrote: {OUT_CSV}")

    # Headline: agreement with human pass label
    print(f"\n{'='*72}")
    print(f"AGREEMENT WITH HUMAN PASS LABEL  (N={len(out_rows)})")
    print(f"{'='*72}")
    for kind in ('cartesian', 'mass_weighted', 'bond_stretch',
                  'core_cartesian', 'core_mass'):
        col = f'{kind}_pass2'
        vals = np.array([r[col] for r in out_rows])
        labels = np.array([r['pass_label'] for r in out_rows])
        # Means within human classes
        m_pass = vals[labels==1].mean()
        m_fail = vals[labels==0].mean() if (labels==0).any() else float('nan')
        gap = m_pass - m_fail
        print(f"\n{kind}_pass2:")
        print(f"  mean | pass=1 (N={int((labels==1).sum())}): {m_pass:.3f}")
        print(f"  mean | pass=0 (N={int((labels==0).sum())}): {m_fail:.3f}")
        print(f"  Δ (pass − fail): {gap:+.3f}")
        # Threshold-based agreement
        for thr in (0.3, 0.5, 0.7):
            tp = ((labels==1) & (vals>=thr)).sum()
            fn = ((labels==1) & (vals< thr)).sum()
            fp = ((labels==0) & (vals>=thr)).sum()
            tn = ((labels==0) & (vals< thr)).sum()
            n = len(vals)
            print(f"  thr={thr}: agree={tp+tn}/{n}={(tp+tn)/n*100:.1f}%   "
                  f"TP={tp}  FP={fp}  TN={tn}  FN={fn}")


if __name__ == '__main__':
    main()

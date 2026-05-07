"""
Verify (or refute) the appendix claim:
  "the second IG is structurally distinct enough from the first to surface a
   different chemical hypothesis" -- i.e. the diversity penalty in the ranker
   is what actually produces the pass@2-over-pass@1 lift.

For every benchmark step we identify the verifier's top-1 and top-2 IGs, then
compute two structural-diversity metrics between them:
  - mass-weighted cosine on the picked imaginary-mode displacement vectors
  - Kabsch RMSD on the IG geometries (Angstrom)

We then split steps by their human-judged outcome:
  LIFT      : IG#1 wrong, IG#2 right  (= the steps that drive the 4.5pp lift)
  P1OK      : IG#1 already correct
  BOTH_BAD  : both wrong
  BOTH_GOOD : both correct
and compare distributions. If the diversity penalty is doing the work, LIFT
steps should have *higher* top1<->top2 diversity than BOTH_GOOD or BOTH_BAD.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent

import csv
import json
import re
import numpy as np

from ranker import rk_clean_v2, mass_weighted_cos


VIEWER_DIR = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
HUMAN_CSV  = PROJECT_ROOT / 'appendix_perparation' / 'analtics' / 'final_quality_measurement-humanversion (1).csv'


def kabsch_rmsd(P, Q):
    P = np.asarray(P, float); Q = np.asarray(Q, float)
    if P.shape != Q.shape: return float('nan')
    Pc = P.mean(0); Qc = Q.mean(0)
    H = (Q - Qc).T @ (P - Pc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Q_aln = (Q - Qc) @ R.T + Pc
    return float(np.sqrt(((P - Q_aln) ** 2).sum() / len(P)))


def canonical(s): return re.sub(r'\s*\(.*?\)\s*$', '', s).strip()


def main():
    # 1) Read human pass labels
    lines = HUMAN_CSV.read_text().splitlines()
    header = lines[1].split(',')
    human = {}
    for ln in lines[2:]:
        if not ln.strip(): continue
        r = dict(zip(header, ln.split(',')))
        try:
            human[canonical(r['step'])] = (int(r['IG#1']), int(r['IG#2']))
        except Exception: pass

    files = sorted(p for p in VIEWER_DIR.glob('*.html')
                   if p.name not in ('flat_view.html', 'guess_quality.html', 'index.html')
                   and not p.name.startswith('oracle_view'))

    rows = []  # (step, group, mwc, rmsd, ig1_label, ig2_label)
    for hp in files:
        text = hp.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if not m: continue
        data = json.loads(m.group(1))
        step = data['step']
        if step not in human: continue
        gt = next((t for t in data['ts_list']
                   if t['label']=='groundtruth' and t.get('modes')), None)
        if gt is None: continue
        elements = gt['xyz_elements']
        igs = [t for t in data['ts_list']
               if t['label']!='groundtruth' and t.get('modes')]
        ranked = rk_clean_v2(igs, elements, k=2)
        if len(ranked) < 2: continue
        (ts1, m1), (ts2, m2) = ranked[0], ranked[1]

        # Diversity metrics
        mwc = mass_weighted_cos(m1['disp'], m2['disp'], elements)
        coords1 = np.asarray(ts1['xyz_coords'], float)
        coords2 = np.asarray(ts2['xyz_coords'], float)
        rmsd = kabsch_rmsd(coords1, coords2)

        ig1_ok, ig2_ok = human[step]
        if   ig1_ok==0 and ig2_ok==1: group = 'LIFT'
        elif ig1_ok==1 and ig2_ok==1: group = 'BOTH_GOOD'
        elif ig1_ok==1 and ig2_ok==0: group = 'P1_ONLY'
        else:                          group = 'BOTH_BAD'
        rows.append((step, group, mwc, rmsd, ts1['label'], ts2['label']))

    print(f"Analyzed {len(rows)} steps with both rankings + human label\n")

    # Summary by group
    groups = {}
    for step, g, mwc, rmsd, l1, l2 in rows:
        groups.setdefault(g, []).append((mwc, rmsd, step, l1, l2))

    print(f"{'group':10s}  {'N':>3s}   {'MWC mean':>9s}  {'med':>5s}  {'p90':>5s}    "
          f"{'RMSD mean':>10s}  {'med':>5s}  {'p90':>5s}")
    print('-' * 80)
    for g in ['BOTH_GOOD', 'P1_ONLY', 'LIFT', 'BOTH_BAD']:
        arr = groups.get(g, [])
        if not arr:
            print(f"{g:10s}  {0:>3d}   (none)")
            continue
        mwcs  = np.array([x[0] for x in arr])
        rmsds = np.array([x[1] for x in arr])
        print(f"{g:10s}  {len(arr):>3d}   "
              f"{mwcs.mean():>9.3f}  {np.median(mwcs):>5.3f}  {np.percentile(mwcs,90):>5.3f}    "
              f"{rmsds.mean():>10.2f}  {np.median(rmsds):>5.2f}  {np.percentile(rmsds,90):>5.2f}")
    print()

    # Detail of LIFT steps
    print("LIFT steps (IG#1 wrong, IG#2 right) — top1<->top2 diversity:")
    print(f"  {'step':46s}  {'top1':10s}  {'top2':10s}  {'MWC':>5s}  {'RMSD':>5s}")
    for step, g, mwc, rmsd, l1, l2 in sorted(rows, key=lambda r: r[2] if r[1]=='LIFT' else 1e9):
        if g != 'LIFT': continue
        print(f"  {step[:46]:46s}  {l1:10s}  {l2:10s}  {mwc:>5.3f}  {rmsd:>5.2f}")

    # Counter-evidence: BOTH_GOOD steps with HIGHER diversity than the LIFT median
    print()
    print("Sanity: are LIFT steps actually MORE diverse than BOTH_GOOD?")
    if 'LIFT' in groups and 'BOTH_GOOD' in groups:
        lift_mwc = np.array([x[0] for x in groups['LIFT']])
        good_mwc = np.array([x[0] for x in groups['BOTH_GOOD']])
        print(f"  LIFT      MWC: median={np.median(lift_mwc):.3f}  (lower MWC = more diverse)")
        print(f"  BOTH_GOOD MWC: median={np.median(good_mwc):.3f}")
        from scipy.stats import mannwhitneyu
        try:
            U, p = mannwhitneyu(lift_mwc, good_mwc, alternative='less')
            print(f"  Mann-Whitney U (LIFT MWC < BOTH_GOOD MWC) p = {p:.3f}")
        except Exception as e:
            print(f"  (scipy unavailable: {e})")

        lift_r = np.array([x[1] for x in groups['LIFT']])
        good_r = np.array([x[1] for x in groups['BOTH_GOOD']])
        print(f"  LIFT      RMSD: median={np.median(lift_r):.2f} A")
        print(f"  BOTH_GOOD RMSD: median={np.median(good_r):.2f} A")


if __name__ == '__main__':
    main()

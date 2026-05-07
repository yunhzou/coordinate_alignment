"""
Figure: initial-guess diversity expressed as IG-vs-GT Kabsch RMSD.

Two side-by-side histograms, framed as a *spread* signal not a
quality signal:

  Left  : all (step, IG) pairs, $N \\approx 155 \\times 20 = 3100$.
          Frames how broadly the 20 IGs span around GT.
  Right : the ranker's top-2 picks per step, $N \\approx 310$.
          Frames that even after ranking the surfaced two are not
          duplicates (the diversity penalty pulls them apart).

Reads:
  appendix_perparation/Pure_Geometries_Elementary_Step/
    Benchmark_Guesses_Coordinate_Aligned_Version/<step>/
      groundtruth/*.xyz, initial_guess/<step>_..._iter<N>_<hash>.xyz
  appendix_perparation/viewer/mode_viewer/<step>.html
    (used to run rk_clean_v2 and identify the top-2 IG labels)

Output:
  appendix_perparation/figures/ig_diversity_rmsd.{png,pdf}
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
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ranker import rk_clean_v2, imag_modes


ALIGNED_DIR = (PROJECT_ROOT / 'appendix_perparation'
               / 'Pure_Geometries_Elementary_Step'
               / 'Benchmark_Guesses_Coordinate_Aligned_Version')
MODE_DIR = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
OUT_DIR  = PROJECT_ROOT / 'appendix_perparation' / 'figures'


def parse_xyz(path):
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines[0])
    coords = []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        coords.append([float(x) for x in parts[1:4]])
    return np.asarray(coords, dtype=float)


def kabsch_rmsd(P, Q):
    P = np.asarray(P, dtype=float); Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape: return np.nan
    Pc = P.mean(0); Qc = Q.mean(0)
    H = (Q - Qc).T @ (P - Pc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Q_aln = (Q - Qc) @ R.T + Pc
    return float(np.sqrt(((P - Q_aln) ** 2).sum() / len(P)))


def iter_label_from_filename(name):
    """`<step>_..._iter<N>_<hash>.xyz` -> 'iterN'."""
    m = re.search(r'_iter(\d+)_', name)
    return f'iter{m.group(1)}' if m else None


def get_top2_labels(step):
    """Run rk_clean_v2 on the per-step HTML payload and return the
    top-2 IG labels (e.g. ['iter5', 'iter12'])."""
    html = MODE_DIR / f"{step}.html"
    if not html.exists(): return []
    text = html.read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m: return []
    data = json.loads(m.group(1))
    gt = next((t for t in data['ts_list']
               if t['label'] == 'groundtruth' and t.get('modes')), None)
    if gt is None: return []
    elements = gt['xyz_elements']
    igs = [t for t in data['ts_list']
           if t['label'] != 'groundtruth' and t.get('modes')]
    ranked = rk_clean_v2(igs, elements)[:2]
    # Pad with bond_overlap fallback if the filter rejected too much
    already = {t['label'] for t, _ in ranked}
    if len(ranked) < 2:
        for ts in igs:
            if ts['label'] in already: continue
            imag = imag_modes(ts)
            if not imag: continue
            picked = max(imag, key=lambda m: m.get('bond_overlap', 0))
            ranked.append((ts, picked))
            if len(ranked) >= 2: break
    return [t['label'] for t, _ in ranked[:2]]


def main():
    steps = sorted(d.name for d in ALIGNED_DIR.iterdir() if d.is_dir())
    print(f"Processing {len(steps)} steps...")
    t0 = time.time()

    all_rmsd = []     # IG-vs-GT for every IG
    top2_rmsd = []    # IG-vs-GT for the verifier's top-2 only

    for i, step in enumerate(steps, 1):
        gt_dir = ALIGNED_DIR / step / 'groundtruth'
        ig_dir = ALIGNED_DIR / step / 'initial_guess'
        gt_files = sorted(gt_dir.glob('*.xyz')) if gt_dir.exists() else []
        ig_files = sorted(ig_dir.glob('*.xyz')) if ig_dir.exists() else []
        if not gt_files or not ig_files: continue
        gt = parse_xyz(gt_files[0])

        per_label_rmsd = {}
        for f in ig_files:
            xyz = parse_xyz(f)
            if xyz.shape != gt.shape: continue
            r = kabsch_rmsd(gt, xyz)
            if np.isnan(r): continue
            label = iter_label_from_filename(f.name)
            if label: per_label_rmsd[label] = r
            all_rmsd.append(r)

        top2 = get_top2_labels(step)
        for lbl in top2:
            if lbl in per_label_rmsd:
                top2_rmsd.append(per_label_rmsd[lbl])

        if i % 30 == 0:
            print(f"  [{i}/{len(steps)}]  ({time.time()-t0:.0f}s)")

    all_rmsd  = np.asarray(all_rmsd)
    top2_rmsd = np.asarray(top2_rmsd)
    print(f"\nAll IG-vs-GT  (N={len(all_rmsd)}):  "
          f"mean={all_rmsd.mean():.2f}Å, median={np.median(all_rmsd):.2f}Å, "
          f"5%-95%=[{np.percentile(all_rmsd,5):.2f}, {np.percentile(all_rmsd,95):.2f}]")
    print(f"Top-2 IG-vs-GT (N={len(top2_rmsd)}): "
          f"mean={top2_rmsd.mean():.2f}Å, median={np.median(top2_rmsd):.2f}Å, "
          f"5%-95%=[{np.percentile(top2_rmsd,5):.2f}, {np.percentile(top2_rmsd,95):.2f}]")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.2), sharey=False)

    upper = float(np.percentile(np.concatenate([all_rmsd, top2_rmsd]), 99.5))
    bins = np.linspace(0, np.ceil(upper * 2) / 2, 45)

    def panel(ax, data, color, title, n_label):
        ax.hist(data, bins=bins, color=color, edgecolor='white', alpha=0.9)
        ax.axvline(np.median(data), color='#cc3333', linestyle='--', lw=1.5,
                    label=f'median = {np.median(data):.2f} Å')
        ax.axvline(data.mean(), color='#dd8800', linestyle=':', lw=1.5,
                    label=f'mean = {data.mean():.2f} Å')
        ax.set_xlabel('Kabsch RMSD to ground-truth TS (Å)', fontsize=11)
        ax.set_ylabel(n_label, fontsize=11)
        ax.set_title(title, fontsize=12, pad=10)
        ax.legend(loc='upper right', fontsize=10, framealpha=0.95)
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        ax.set_axisbelow(True)
        ax.set_xlim(0, upper * 1.05)

    panel(ax_a, all_rmsd,
          '#3a6dbf',
          f'All 20 IGs vs. GT  (N = {len(all_rmsd)})',
          '# (step, IG) pairs')
    panel(ax_b, top2_rmsd,
          '#cc6699',
          f"Ranker's top-2 IGs vs. GT  (N = {len(top2_rmsd)})",
          '# (step, top-2) pairs')

    fig.suptitle('Initial-guess geometric diversity '
                 '(structural spread, not a quality metric)',
                 fontsize=13, y=1.0)
    fig.tight_layout()

    out_png = OUT_DIR / 'ig_diversity_rmsd.png'
    out_pdf = OUT_DIR / 'ig_diversity_rmsd.pdf'
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nSaved:\n  {out_png}\n  {out_pdf}")


if __name__ == '__main__':
    main()

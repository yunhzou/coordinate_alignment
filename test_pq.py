"""
Compare priority-queue alignment (rxn_core_pq) against the current
fragment-based alignment (rxn_core_frag) on a curated set of BGCP steps.

For each step:
  baseline: analyze() from rxn_core_frag    (graph_floor 0.5, no chirality)
  new:      analyze_pq() from rxn_core_pq   (graph_floor 0.2, branching, chirality)

Print a table: step | n_atoms | mapped_old/mapped_new | br/fm_old | br/fm_new | chir_new

Usage: python test_pq.py
"""
from __future__ import annotations

import re
import time
import traceback
from pathlib import Path

from rxn_core_frag import analyze
from rxn_core_pq import analyze_pq

ROOT = Path(__file__).parent
BGCP = ROOT / "Benchmark_Guesses_Collective_Package"
WORK_OLD = ROOT / "work_pq_compare_old"
WORK_NEW = ROOT / "work_pq_compare_new"
WORK_OLD.mkdir(exist_ok=True)
WORK_NEW.mkdir(exist_ok=True)


# (step_dir_name, charge, uhf)
STEPS = [
    ("pr1.tempo_ts2", 0, 0),
    ("pr11.cycloadditions_tsIa", 0, 0),
    ("pr12.Co_Silylation_JACS2015_TS_Dstar-Estar", 0, 0),
    ("pr14.Pd_hydroamination_JOC2025_TS14_step4_reductive_elimination", 0, 0),
    ("pr9.carbene.rearr_ts47a", 0, 0),
]


def concat_xyzs(dir_path):
    """Concat all xyz files in dir, translating fragments apart by 10 Å."""
    files = sorted(dir_path.glob("*.xyz"))
    if not files:
        return None
    chunks = []
    offset = 0.0
    for f in files:
        lines = f.read_text().splitlines()
        n = int(lines[0])
        coords = []
        elems = []
        for ln in lines[2:2 + n]:
            parts = ln.split()
            elems.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
        for i in range(n):
            coords[i][0] += offset
        for el, c in zip(elems, coords):
            chunks.append(f"{el} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}")
        offset += 10.0
    total = sum(int(f.read_text().splitlines()[0]) for f in files)
    return f"{total}\n\n" + "\n".join(chunks) + "\n"


def step_inputs(step_name):
    sd = BGCP / step_name
    rxyz = concat_xyzs(sd / "reactants")
    pxyz = concat_xyzs(sd / "products")
    if rxyz is None or pxyz is None:
        raise RuntimeError(f"missing R or P for {step_name}")
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step_name)
    rfile = WORK_OLD / f"{sanitized}_R.xyz"
    pfile = WORK_OLD / f"{sanitized}_P.xyz"
    rfile.write_text(rxyz)
    pfile.write_text(pxyz)
    return rfile, pfile


def fmt_bonds(bonds, elements):
    if not bonds:
        return "—"
    return ";".join(f"{i}({elements[i]})-{j}({elements[j]})"
                    for (i, j, _, _) in bonds[:6]) + (
        f" +{len(bonds)-6}" if len(bonds) > 6 else "")


def main():
    rows = []
    for step, chg, uhf in STEPS:
        sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step)
        try:
            rfile, pfile = step_inputs(step)
        except Exception as e:
            print(f"[skip] {step}: {e}")
            continue

        # baseline
        t0 = time.time()
        try:
            old = analyze(rfile, pfile, WORK_OLD / sanitized,
                          charge=chg, uhf=uhf)
            t_old = time.time() - t0
            old_n = len(old['elements_R'])
            old_m = len(old['mapping'])
            old_br = len(old['broken'])
            old_fm = len(old['formed'])
            old_br_str = fmt_bonds(old['broken'], old['elements_R'])
            old_fm_str = fmt_bonds(old['formed'], old['elements_P'])
        except Exception as e:
            print(f"[old fail] {step}: {e}")
            traceback.print_exc()
            continue

        # new
        t1 = time.time()
        try:
            new = analyze_pq(rfile, pfile, WORK_NEW / sanitized,
                             charge=chg, uhf=uhf)
            t_new = time.time() - t1
            new_n = len(new['elements_R'])
            new_m = new['n_mapped']
            new_br = new['n_broken']
            new_fm = new['n_formed']
            new_chir = new['chirality_violations']
            new_br_str = fmt_bonds(new['broken'], new['elements_R'])
            new_fm_str = fmt_bonds(new['formed'], new['elements_P'])
        except Exception as e:
            print(f"[new fail] {step}: {e}")
            traceback.print_exc()
            continue

        rows.append({
            'step': step, 'n': old_n,
            'm_old': old_m, 'm_new': new_m,
            'br_old': old_br, 'fm_old': old_fm,
            'br_new': new_br, 'fm_new': new_fm,
            'chir_new': new_chir,
            't_old': t_old, 't_new': t_new,
            'br_old_str': old_br_str, 'fm_old_str': old_fm_str,
            'br_new_str': new_br_str, 'fm_new_str': new_fm_str,
        })

    print()
    print("=" * 100)
    print(f"{'step':50s} {'N':>4} {'map old/new':>12} {'br/fm old':>10} "
          f"{'br/fm new':>10} {'chir':>5} {'Δt':>10}")
    print("-" * 100)
    for r in rows:
        print(f"{r['step'][:50]:50s} {r['n']:>4} "
              f"{r['m_old']:>5}/{r['m_new']:<6} "
              f"{r['br_old']}/{r['fm_old']:>7} "
              f"{r['br_new']}/{r['fm_new']:>7} "
              f"{r['chir_new']:>5} "
              f"{r['t_old']:>4.1f}s/{r['t_new']:>4.1f}s")
    print()
    print("Per-step bond detail:")
    for r in rows:
        print(f"  {r['step']}")
        print(f"    OLD br: {r['br_old_str']}")
        print(f"    NEW br: {r['br_new_str']}")
        print(f"    OLD fm: {r['fm_old_str']}")
        print(f"    NEW fm: {r['fm_new_str']}")


if __name__ == "__main__":
    main()

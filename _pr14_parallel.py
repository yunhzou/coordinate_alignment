"""Parallel cut_sweep on pr14 using all CPU cores.

Work unit: one (cut, seed_idx) pair = one find_islands_pq alignment.
Workers run independently, return chem signatures + mappings. Main process
deduplicates by chemistry signature.

Expected: 113 cuts × N seeds = many alignments, ~2s each serial, ~50s with 14 cores.
"""
from __future__ import annotations
import sys, time
sys.path.insert(0, 'src')
import multiprocessing as mp
import random
from pathlib import Path
import numpy as np

# Defer rxn_core imports to inside workers (faster process startup).
from rxn_core import parse_xyz, classify_bonds
from rxn_core.pq import build_graph

WORK = Path('appendix_perparation/xtb_frequency_calculations')
WBO_STRONG = 0.5
N_SEEDS_PER_CUT = 3


def load_step(d):
    xyz_path = next(p for p in d.glob('*.xyz') if 'xtbhess' not in p.name)
    el, xyz = parse_xyz(xyz_path)
    n = len(el); wbo = np.zeros((n, n))
    for ln in (d / 'wbo').read_text().splitlines():
        p = ln.split()
        if len(p) < 3: continue
        i, j = int(p[0])-1, int(p[1])-1
        wbo[i, j] = float(p[2]); wbo[j, i] = wbo[i, j]
    return el, np.asarray(xyz, float), wbo


# Globals set by worker init (one-time per process)
_W = {}

def _worker_init(elR, wboR, elT, wboT):
    """Cache (elR, wboR, elT, wboT) per worker; rebuild graphs once."""
    _W['elR'] = elR
    _W['wboR'] = wboR
    _W['elT'] = elT
    _W['wboT'] = wboT
    _W['g_P'] = build_graph(elT, wboT, bond_cut=0.2)
    _W['n'] = len(elR)


def _worker_run(args):
    """Run ONE alignment: (cut_tuple, seed_order). Return list of
    (chem_signature, mapping_items) for full bijections."""
    cut, order = args
    from rxn_core.pq import find_islands_pq, expand_mapping
    elR = _W['elR']; wboR = _W['wboR']; wboT = _W['wboT']
    g_R = build_graph(elR, wboR, bond_cut=0.2)
    for (i, j) in cut:
        if g_R.has_edge(i, j): g_R.remove_edge(i, j)
    g_P = _W['g_P']
    try:
        branches = find_islands_pq(g_R, g_P, list(order))
    except Exception as e:
        return [('__error__', str(e))]
    out = []
    for b in branches:
        mapping = expand_mapping(b.mapping, g_R, g_P)
        if len(mapping) < _W['n'] - 2: continue
        broken, formed, _, _ = classify_bonds(mapping, wboR, wboT)
        inv = {v: k for k, v in mapping.items()}
        br = tuple(sorted((min(a, b), max(a, b)) for (a, b, _, _) in broken))
        fm = tuple(sorted((min(inv.get(a, -1), inv.get(b, -1)),
                           max(inv.get(a, -1), inv.get(b, -1)))
                          for (a, b, _, _) in formed if a in inv and b in inv))
        out.append(((br, fm), tuple(sorted(mapping.items())), cut))
    return out


def main(step_name, n_workers=None):
    if n_workers is None:
        import os
        n_workers = os.cpu_count() or 1
    sd = WORK / step_name
    elR, xyzR, wboR = load_step(sd / 'R')
    elT, xyzT, wboT = load_step(sd / 'P')
    print(f'{step_name}: {len(elR)} atoms, {n_workers} workers', flush=True)

    # Generate work units: 1 baseline + one per strong R-bond cut
    n = len(elR)
    strong = [(i, j) for i in range(n) for j in range(i+1, n)
              if wboR[i, j] >= WBO_STRONG]
    print(f'  cuts: 1 baseline + {len(strong)} strong-edge cuts = {1+len(strong)} cut sets', flush=True)

    # Per cut: N_SEEDS_PER_CUT random seed orderings
    g_R = build_graph(elR, wboR, bond_cut=0.2)
    nodes = list(g_R.nodes())
    rng = random.Random(42)
    seed_orders = []
    for s in range(N_SEEDS_PER_CUT):
        order = list(nodes); rng.shuffle(order)
        seed_orders.append(tuple(order))

    work_units = []
    work_units.append(((), seed_orders[0]))  # baseline, single seed (cheap)
    for s in range(1, N_SEEDS_PER_CUT):
        work_units.append(((), seed_orders[s]))
    for (i, j) in strong:
        for s in range(N_SEEDS_PER_CUT):
            work_units.append((((i, j),), seed_orders[s]))
    print(f'  total work units: {len(work_units)}', flush=True)

    pool_chems = {}    # chem_sig -> (mapping_items, cut)
    n_errors = 0
    t0 = time.time()
    t_last = t0
    with mp.Pool(n_workers, initializer=_worker_init,
                  initargs=(elR, wboR, elT, wboT)) as pool:
        for i, results in enumerate(pool.imap_unordered(_worker_run, work_units, chunksize=4), 1):
            for entry in results:
                if entry[0] == '__error__':
                    n_errors += 1
                    continue
                chem_sig, mapping_items, cut = entry
                pool_chems.setdefault(chem_sig, (mapping_items, cut))
            if i % 50 == 0 or i == len(work_units):
                now = time.time()
                print(f'  [{i:4d}/{len(work_units)}]  pool={len(pool_chems)} chems  '
                      f'  elapsed={now-t0:6.1f}s  delta={now-t_last:5.1f}s  errors={n_errors}',
                      flush=True)
                t_last = now

    total = time.time() - t0
    print(f'\nDONE in {total:.1f}s.  {len(pool_chems)} chem classes, {n_errors} worker errors', flush=True)

    if not pool_chems:
        print('  NO CHEMS FOUND'); return
    mn = min(len(k[0]) + len(k[1]) for k in pool_chems)
    print(f'\nmin br+fm = {mn}, mechanisms at min:')
    for (br, fm), (mapping_items, cut) in sorted(pool_chems.items(),
                                                  key=lambda x: len(x[0][0])+len(x[0][1])):
        if len(br) + len(fm) > mn: continue
        br_str = ', '.join(f'{elR[a]}{a}-{elR[b]}{b}' for a, b in br)
        fm_str = ', '.join(f'{elR[a]}{a}-{elR[b]}{b}' for a, b in fm)
        cut_str = 'none' if not cut else ','.join(f'{elR[a]}{a}-{elR[b]}{b}' for a, b in cut)
        print(f'  br/fm={len(br)}/{len(fm)}  cut=[{cut_str}]')
        print(f'    broken: {br_str}')
        print(f'    formed: {fm_str}')


if __name__ == '__main__':
    step = sys.argv[1] if len(sys.argv) > 1 else 'pr14.Pd_hydroamination_JOC2025_TS3_step2_alkene_inserion'
    main(step)

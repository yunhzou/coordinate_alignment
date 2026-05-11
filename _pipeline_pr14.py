"""Full multi-mechanism ranking pipeline for pr14 with max core usage.

Per-step flow:
  1. parallel cut_sweep(R, P)   -> chem pool, pick min-bond mechanisms
  2. parallel cut_sweep(R, GT)  -> R↔GT bijection pool
  3. for each of 20 IGs: parallel cut_sweep(R, IG) -> R↔IG bijection pool
  4. for each mechanism m (with its own core_R):
       - dedup gt_pool by mapping restricted to core_R; score each rep;
         keep highest-S rep -> m.gt
       - for each IG: dedup ig_pool[ig] by core_R restriction; score reps;
         keep highest-S rep -> m.igs[ig]
  5. rank IGs per mech by S; flag top-2; write JSON + HTML

CORE-RESTRICTED DEDUP RATIONALE (for R↔TS / R↔IG, NOT R↔P):
  R↔P alignment dedups by (broken, formed) chemistry signature because
  spectator-permutation siblings produce identical chemistry — they're
  the SAME mechanism. R↔TS alignment instead dedups by the bijection's
  RESTRICTION TO core_R, because every term in the verifier score
  (beta, rho, kappa, V) only reads TS coords / mode displacements at
  R-frame core indices; spectator atom permutations don't affect any of
  those quantities. We use the *per-mechanism* core_R as the dedup key
  (option a) because different mechanisms can have different core_R, and
  per-mech dedup keeps the pool small without losing chemistry-relevant
  alternatives within each mechanism.
"""
from __future__ import annotations
import sys, json, time, re
sys.path.insert(0, 'src')
import multiprocessing as mp
import os, random
from pathlib import Path
import numpy as np

from rxn_core import (parse_xyz, classify_bonds, parse_g98_modes,
                      core_atoms_in_R_frame, fill_unmapped_greedy,
                      reaction_coord_delta, reindex_modes_to_R,
                      bond_overlap_per_mode, bond_reaction_vector,
                      rxn_overlap_per_mode)
from rxn_core.pq import build_graph, find_islands_pq, expand_mapping

WORK = Path('appendix_perparation/xtb_frequency_calculations')
OUT_ROOT = Path('out/bgcp_views')
STEP = 'pr14.Pd_hydroamination_JOC2025_TS3_step2_alkene_inserion'
WBO_STRONG = 0.5
N_SEEDS_PER_CUT = 3
W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3
N_WORKERS = os.cpu_count() or 14


def load(d):
    xyz_path = next(p for p in d.glob('*.xyz') if 'xtbhess' not in p.name)
    el, xyz = parse_xyz(xyz_path)
    n = len(el); wbo = np.zeros((n, n))
    for ln in (d / 'wbo').read_text().splitlines():
        p = ln.split()
        if len(p) < 3: continue
        i, j = int(p[0])-1, int(p[1])-1
        wbo[i, j] = float(p[2]); wbo[j, i] = wbo[i, j]
    return el, np.asarray(xyz, float), wbo


# ============================================================================
# Parallel cut_sweep
# ============================================================================
_W = {}
def _winit(elR, wboR, elT, wboT):
    _W['elR'] = elR; _W['wboR'] = wboR
    _W['elT'] = elT; _W['wboT'] = wboT
    _W['g_P'] = build_graph(elT, wboT, bond_cut=0.2)
    _W['n'] = len(elR)


def _wrun(args):
    cut, order = args
    g_R = build_graph(_W['elR'], _W['wboR'], bond_cut=0.2)
    for (i, j) in cut:
        if g_R.has_edge(i, j): g_R.remove_edge(i, j)
    try:
        branches = find_islands_pq(g_R, _W['g_P'], list(order))
    except Exception:
        return []
    out = []
    for b in branches:
        mapping = expand_mapping(b.mapping, g_R, _W['g_P'])
        if len(mapping) < _W['n'] - 2: continue
        broken, formed, _, _ = classify_bonds(mapping, _W['wboR'], _W['wboT'])
        inv = {v: k for k, v in mapping.items()}
        br = tuple(sorted((min(a, b), max(a, b)) for (a, b, _, _) in broken))
        fm_R = tuple(sorted((min(inv.get(a, -1), inv.get(b, -1)),
                              max(inv.get(a, -1), inv.get(b, -1)))
                             for (a, b, _, _) in formed if a in inv and b in inv))
        out.append(((br, fm_R), tuple(sorted(mapping.items())), cut))
    return out


def parallel_cut_sweep(elR, wboR, elT, wboT, n_workers=N_WORKERS, label=''):
    """Returns chem_sig -> {'mapping': dict, 'cut': tuple}."""
    t0 = time.time()
    n = len(elR)
    strong = [(i, j) for i in range(n) for j in range(i+1, n)
              if wboR[i, j] >= WBO_STRONG]
    g_R = build_graph(elR, wboR, bond_cut=0.2)
    nodes = list(g_R.nodes())
    rng = random.Random(42)
    seed_orders = []
    for _ in range(N_SEEDS_PER_CUT):
        order = list(nodes); rng.shuffle(order)
        seed_orders.append(tuple(order))
    work = [((), s) for s in seed_orders]
    for (i, j) in strong:
        for s in seed_orders: work.append((((i, j),), s))
    pool = {}
    with mp.Pool(n_workers, initializer=_winit,
                  initargs=(elR, wboR, elT, wboT)) as p:
        for results in p.imap_unordered(_wrun, work, chunksize=4):
            for chem_sig, mapping_items, cut in results:
                pool.setdefault(chem_sig, {
                    'mapping': dict(mapping_items),
                    'cut': cut,
                })
    print(f'  cut_sweep {label}: {len(pool)} chems in {time.time()-t0:.1f}s '
          f'({len(work)} units, {n_workers}w)', flush=True)
    return pool


def select_min(pool):
    if not pool: return {}
    m = min(len(k[0]) + len(k[1]) for k in pool)
    return {k: v for k, v in pool.items() if len(k[0]) + len(k[1]) == m}


# ============================================================================
# Scoring (single-mapping)
# ============================================================================
def score_one(elR, xyzR, elT, xyzT, mapping_RT, freqs, modes_TS,
              broken_R, formed_R, core_R, delta_RP):
    modes_R = reindex_modes_to_R(modes_TS, mapping_RT, len(elR))
    sq = (modes_R**2).sum(axis=2); total = sq.sum(axis=1)
    core_e = sq[:, core_R].sum(axis=1) if core_R else np.zeros(modes_R.shape[0])
    kappa = np.where(total > 1e-12, core_e / total, 0.0)
    rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R)
    ts_in_R = np.zeros_like(np.asarray(xyzR, float))
    for r, t in mapping_RT.items(): ts_in_R[r] = xyzT[t]
    V = bond_reaction_vector(ts_in_R, broken_R, formed_R)
    beta = bond_overlap_per_mode(modes_R, V)
    imag = list(np.where(freqs < 0)[0])
    if not imag: return None
    pk = max(imag, key=lambda k: beta[k])
    return {'S': float(beta[pk]*(1+W_RXN*rho[pk])*(1+W_CORE*kappa[pk])
                       / max(len(imag), 1)**IMAG_PEN),
            'beta': float(beta[pk]), 'rho': float(rho[pk]),
            'kappa': float(kappa[pk]),
            'freq': float(freqs[pk]), 'k': int(pk),
            'n_imag': len(imag)}


def best_with_core_dedup(rt_pool, elR, xyzR, elT, xyzT, freqs, modes_TS,
                          broken_R, formed_R, core_R, delta_RP):
    """Score every R↔T bijection in rt_pool that's UNIQUE on core_R.

    Two bijections that agree on { r → mapping[r] | r ∈ core_R } produce
    identical scores under this mechanism (the verifier only reads TS data
    at core indices), so we keep one rep per core-restricted key. Returns
    the highest-scoring rep's score dict.
    """
    core_R_set = frozenset(core_R)
    seen_core = set()
    best = None
    n_total = len(rt_pool); n_scored = 0
    for v in rt_pool.values():
        m_full = v['mapping']
        core_key = frozenset((r, m_full[r]) for r in core_R_set if r in m_full)
        if core_key in seen_core:
            continue
        seen_core.add(core_key)
        mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, m_full)
        s = score_one(elR, xyzR, elT, xyzT, mapping_RT, freqs, modes_TS,
                      broken_R, formed_R, core_R, delta_RP)
        n_scored += 1
        if s and (best is None or s['S'] > best['S']):
            best = s
    return best, n_total, n_scored


# ============================================================================
# Pipeline
# ============================================================================
def main():
    t_total = time.time()
    sd = WORK / STEP
    print(f'=== {STEP} ===', flush=True)
    elR, xyzR, wboR = load(sd / 'R')
    elP, xyzP, wboP = load(sd / 'P')
    print(f'  {len(elR)} atoms, {N_WORKERS} workers', flush=True)

    # 1. R↔P cut_sweep
    print('\n[1/5] R↔P cut_sweep', flush=True)
    rp_pool = parallel_cut_sweep(elR, wboR, elP, wboP, label='R↔P')
    rp_min = select_min(rp_pool)
    if not rp_min:
        print('  no min-bond mechanism found'); return
    print(f'  {len(rp_min)} mechanism(s) at min br+fm={min(len(k[0])+len(k[1]) for k in rp_pool)}',
          flush=True)

    # Build mechanism state
    mechanisms = []
    for mi, ((br_t, fm_t), info) in enumerate(rp_min.items(), 1):
        mapping_RP = info['mapping']
        inv_RP = {v: k for k, v in mapping_RP.items()}
        broken, formed, _, _ = classify_bonds(mapping_RP, wboR, wboP)
        broken_R = [(int(a), int(b)) for (a, b, _, _) in broken]
        formed_R = [(int(inv_RP[a]), int(inv_RP[b])) for (a, b, _, _) in formed
                    if a in inv_RP and b in inv_RP]
        core_R = list(core_atoms_in_R_frame(mapping_RP, broken, formed))
        full_RP = fill_unmapped_greedy(elR, xyzR, elP, xyzP, mapping_RP)
        delta_RP = reaction_coord_delta(np.asarray(xyzR), np.asarray(xyzP), full_RP)
        cut = info['cut']
        cut_name = ','.join(f'{elR[a]}{a}-{elR[b]}{b}' for a, b in cut) or 'none'
        br_label = ','.join(f'{elR[a]}{a}-{elR[b]}{b}' for a, b in br_t)
        mechanisms.append({
            'id': mi, 'cut': cut_name, 'label': f'#{mi}: {br_label} (cut: {cut_name})',
            'broken_R': broken_R, 'formed_R': formed_R, 'core_R': core_R,
            'delta_RP': delta_RP,
            'mapping_RP': mapping_RP,
        })
        print(f'    mech #{mi}: cut={cut_name}  core_R size={len(core_R)}', flush=True)

    # 2. R↔GT cut_sweep
    print('\n[2/5] R↔GT cut_sweep', flush=True)
    elT_gt, xyzT_gt, wboT_gt = load(sd / 'sp_groundtruth')
    freqs_gt, modes_gt = parse_g98_modes(sd / 'hess_groundtruth' / 'g98.out')
    gt_pool = parallel_cut_sweep(elR, wboR, elT_gt, wboT_gt, label='R↔GT')

    # 3. R↔IG cut_sweep for each of 20 IGs
    print('\n[3/5] R↔IG cut_sweep (×20)', flush=True)
    iter_dirs = sorted([d for d in sd.iterdir()
                        if d.is_dir() and re.match(r'hess_iter(\d+)$', d.name)],
                       key=lambda d: int(re.match(r'hess_iter(\d+)$', d.name).group(1)))
    ig_data = []   # list of (label, elI, xyzI, freqs, modes, ig_pool)
    for hess_dir in iter_dirs:
        label = hess_dir.name.replace('hess_', '')
        sp_dir = sd / f'sp_{label}'
        try:
            elI, xyzI, wboI = load(sp_dir)
            freqs_i, modes_i = parse_g98_modes(hess_dir / 'g98.out')
        except Exception as e:
            print(f'  {label}: missing data ({e})', flush=True)
            ig_data.append((label, None, None, None, None, None))
            continue
        ig_pool = parallel_cut_sweep(elR, wboR, elI, wboI, label=f'R↔{label}')
        ig_data.append((label, elI, xyzI, freqs_i, modes_i, ig_pool))

    # 4. Score every mech × (GT + 20 IGs) with core-restricted dedup
    print('\n[4/5] Score per mech with core-restricted dedup', flush=True)
    for mech in mechanisms:
        # GT under this mech
        s_gt, n_tot, n_scored = best_with_core_dedup(
            gt_pool, elR, xyzR, elT_gt, xyzT_gt, freqs_gt, modes_gt,
            mech['broken_R'], mech['formed_R'], mech['core_R'], mech['delta_RP'])
        mech['gt'] = s_gt
        print(f'  mech #{mech["id"]}: GT  pool={n_tot} -> {n_scored} unique on core  '
              f'S={s_gt["S"]:.3f}' if s_gt else
              f'  mech #{mech["id"]}: GT  no score', flush=True)
        # IGs under this mech
        mech['igs'] = []
        for label, elI, xyzI, freqs_i, modes_i, ig_pool in ig_data:
            if elI is None or ig_pool is None:
                mech['igs'].append({'label': label})
                continue
            s_ig, n_tot, n_scored = best_with_core_dedup(
                ig_pool, elR, xyzR, elI, xyzI, freqs_i, modes_i,
                mech['broken_R'], mech['formed_R'], mech['core_R'], mech['delta_RP'])
            entry = {'label': label}
            if s_ig: entry.update(s_ig)
            entry['core_unique'] = n_scored
            entry['pool_total'] = n_tot
            mech['igs'].append(entry)

    # 5. Rank per mech, dump
    print('\n[5/5] Rank + dump', flush=True)
    union_top = set()
    for mech in mechanisms:
        ranked = sorted([(i, ig) for i, ig in enumerate(mech['igs']) if ig.get('S') is not None],
                        key=lambda x: -x[1]['S'])
        top2 = {i for i, _ in ranked[:2]}
        for i, ig in enumerate(mech['igs']):
            ig['is_top2'] = (i in top2)
            if i in top2: union_top.add(ig['label'])
        print(f'\n  mech #{mech["id"]}: cut={mech["cut"]}  GT S='
              f'{mech["gt"]["S"]:.3f}' if mech['gt'] else
              f'  mech #{mech["id"]}: no GT score', flush=True)
        for rank, (i, ig) in enumerate(ranked[:5], 1):
            print(f'    #{rank}  {ig["label"]:<10}  S={ig["S"]:.3f}  '
                  f'beta={ig.get("beta", 0):.3f}  rho={ig.get("rho", 0):.3f}  '
                  f'kappa={ig.get("kappa", 0):.3f}  freq={ig.get("freq", 0):.0f}i  '
                  f'core_unique={ig["core_unique"]}/{ig["pool_total"]}',
                  flush=True)
    print(f'\nunion top-2 across mechs: {sorted(union_top)}', flush=True)

    # Dump JSON
    out_dir = OUT_ROOT / STEP
    out_dir.mkdir(parents=True, exist_ok=True)
    slim = {'step': STEP, 'n_atoms': len(elR), 'n_mechs': len(mechanisms),
            'wall_time_s': round(time.time() - t_total, 1),
            'mechanisms': []}
    for mech in mechanisms:
        slim['mechanisms'].append({
            'id': mech['id'], 'cut': mech['cut'],
            'broken_R': mech['broken_R'], 'formed_R': mech['formed_R'],
            'core_R': mech['core_R'],
            'gt': {k: mech['gt'][k] for k in ['S', 'beta', 'rho', 'kappa',
                                                'freq', 'n_imag']} if mech['gt'] else None,
            'igs': [{k: ig.get(k) for k in ['label', 'S', 'beta', 'rho', 'kappa',
                                              'freq', 'n_imag', 'is_top2',
                                              'core_unique', 'pool_total']}
                    for ig in mech['igs']],
        })
    (out_dir / 'pr14_pipeline.json').write_text(json.dumps(slim, indent=2))
    print(f'\nwrote {out_dir / "pr14_pipeline.json"}', flush=True)
    print(f'TOTAL wall-time: {time.time() - t_total:.1f}s', flush=True)


if __name__ == '__main__':
    main()

"""Build interactive multi-mechanism views for BGCP cached steps.

For each step:
  1. R<->P cut-sweep -> min-bondcount mechanisms.
  2. Under each mechanism: core-match GT + each IG from R and P endpoints,
     then pick best R-frame core witness.
  3. Rank IGs per mech, mark top-2; union across mechs.
  4. Write out/bgcp_views/<step>/view.html with mechanism switcher.
  5. Dump out/bgcp_alignment_eval.json for downstream CSV.

Parallelized via multiprocessing.
"""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import json
import multiprocessing as mp
import os
import re
import shutil
import time
import traceback
from dataclasses import dataclass
from pathlib import Path
import numpy as np

from rxn_core import (parse_xyz, classify_bonds, parse_g98_modes,
                      core_atoms_in_R_frame,
                      reindex_modes_to_R, bond_overlap_per_mode,
                      build_graph, cut_sweep, select_min_mechanisms,
                      ts_core_pool)
from rxn_core.chemistry_computations import (
    run_xtb, run_xtb_hess, write_xyz_str, xyz_with_disp,
)
from rxn_core.matcher import _nauty_orbits

PROJECT = Path(os.environ.get(
    "RXN_CORE_PROJECT",
    Path(__file__).resolve().parents[2],
))
WORK = Path(os.environ.get(
    "BGCP_WORK",
    PROJECT / "data" / "xtb_frequency_calculations",
))
OUT_ROOT = Path(os.environ.get(
    "BGCP_OUT_ROOT",
    PROJECT / "out" / "bgcp_views",
))
STAGE_ROOT = Path(os.environ.get(
    "BGCP_STAGE_ROOT",
    PROJECT / "out" / "bgcp_stages",
))
ALIGNMENT_OUT_ROOT = Path(os.environ.get(
    "BGCP_ALIGNMENT_OUT_ROOT",
    PROJECT / "out" / "bgcp_alignments",
))
EVAL_JSON = Path(os.environ.get(
    "BGCP_EVAL_JSON",
    PROJECT / "out" / "bgcp_alignment_eval.json",
))
CUT_FLOOR = float(os.environ.get("BGCP_CUT_FLOOR", "0.2"))
N_SEEDS_PER_RUN = 3  # cut + seed are orthogonal diversity sources; keep both modest
VIEW_MAX_BRANCHES = int(os.environ.get("BGCP_VIEW_MAX_BRANCHES", "100"))
CUTSWEEP_CHUNKSIZE = int(os.environ.get("BGCP_CUTSWEEP_CHUNKSIZE", "1"))
VIEW_ISO_TOL = float(os.environ.get("BGCP_ISO_TOL", "1.0"))
DWBO_THRESHOLD = float(os.environ.get("BGCP_DWBO_THRESHOLD", "0.5"))
METAL_DWBO_THRESHOLD = float(os.environ.get("BGCP_METAL_DWBO_THRESHOLD", "0.3"))
SYMMETRY_WBO_TOL = float(os.environ.get("BGCP_SYMMETRY_WBO_TOL", "0.2"))
BGCP_TIMING = os.environ.get("BGCP_TIMING", "0") == "1"
INCLUDE_GT = os.environ.get("BGCP_INCLUDE_GT", "0").lower() in {
    "1", "true", "yes", "on"
}
SYMMETRY_REPAIR = os.environ.get("BGCP_SYMMETRY_REPAIR", "1") != "0"
SYMMETRY_REPAIR_MIN_CHANGES = int(os.environ.get("BGCP_SYMMETRY_REPAIR_MIN_CHANGES", "1"))
SYMMETRY_REPAIR_MAX_EVALS = int(os.environ.get("BGCP_SYMMETRY_REPAIR_MAX_EVALS", "20000"))
TS_CORE_EDGE_FLOOR = float(os.environ.get("BGCP_TS_CORE_EDGE_FLOOR", "0.2"))
TS_CORE_MAX_CANDIDATES = int(os.environ.get("BGCP_TS_CORE_MAX_CANDIDATES", "20000"))
AUTO_INNER_WORKERS = int(os.environ.get("BGCP_AUTO_INNER_WORKERS", "8"))
XTB_CACHE_MODE = os.environ.get("BGCP_XTB_MODE", "auto").lower()
XTB_OMP_THREADS = os.environ.get("BGCP_XTB_OMP_THREADS", "auto")
XTB_MAX_THREADS = int(os.environ.get("BGCP_XTB_MAX_THREADS", "8"))
XTB_CHARGE = int(os.environ.get("BGCP_CHARGE", "0"))
XTB_MULTIPLICITY = int(os.environ.get("BGCP_MULTIPLICITY", "1"))
EVENT_WEIGHT_POWER = float(os.environ.get("BGCP_EVENT_WEIGHT_POWER", "1.0"))
WBO_PROGRESS_POWER = float(os.environ.get("BGCP_WBO_PROGRESS_POWER", "1.0"))


@dataclass
class StepInputs:
    """Loaded endpoint data for one R/P alignment problem.

    The name and directory fields are artifact labels, not a required input
    schema.  Build this object from arrays, direct XYZ files, or a benchmark
    step adapter; downstream stages only depend on the molecule data here.
    """
    step_name: str
    step_dir: Path
    elR: list
    xyzR: np.ndarray
    wboR: np.ndarray
    elP: list
    xyzP: np.ndarray
    wboP: np.ndarray


@dataclass
class TSTarget:
    """One GT/IG/TS candidate to verify under discovered mechanisms."""
    kind: str
    target_index: int
    label: str
    el: list
    xyz: np.ndarray
    wbo: np.ndarray
    freqs: np.ndarray
    modes: np.ndarray


@dataclass
class StagePaths:
    """On-disk artifacts for separately runnable stages."""
    root: Path
    rp_json: Path
    ts_json: Path
    view_dir: Path
    view_html: Path
    eval_slim_json: Path


def _default_worker_count():
    return max(1, int(os.cpu_count() or 2) - 1)


def pipeline_stage_paths(step_name, stage_root=None, out_root=None):
    """Return the standard artifact paths for one step."""
    stage_root = Path(stage_root or STAGE_ROOT)
    out_root = Path(out_root or OUT_ROOT)
    stage_dir = stage_root / step_name
    view_dir = out_root / step_name
    return StagePaths(
        root=stage_dir,
        rp_json=stage_dir / "rp_stage.json",
        ts_json=stage_dir / "ts_stage.json",
        view_dir=view_dir,
        view_html=view_dir / "view.html",
        eval_slim_json=view_dir / "_eval_slim.json",
    )


def alignment_output_dir(step_name, alignment_out_root=None):
    """Directory for clean Stage 1 aligned-coordinate exports."""
    return Path(alignment_out_root or ALIGNMENT_OUT_ROOT) / step_name


def rp_stage_config():
    """Current R-P mechanism-discovery hypotheses as a serializable dict."""
    return {
        'cut_floor': CUT_FLOOR,
        'graph_floor': 0.2,
        'iso_tol': VIEW_ISO_TOL,
        'dwbo_threshold': DWBO_THRESHOLD,
        'metal_dwbo_threshold': METAL_DWBO_THRESHOLD,
        'symmetry_wbo_tol': SYMMETRY_WBO_TOL,
        'n_seeds': N_SEEDS_PER_RUN,
        'max_branches': VIEW_MAX_BRANCHES,
        'chunksize': CUTSWEEP_CHUNKSIZE,
        'symmetry_repair': SYMMETRY_REPAIR,
        'symmetry_repair_min_changes': SYMMETRY_REPAIR_MIN_CHANGES,
        'symmetry_repair_max_evals': SYMMETRY_REPAIR_MAX_EVALS,
    }


def ts_stage_config():
    """Current TS/IG verification hypotheses as a serializable dict."""
    return {
        'iso_tol': VIEW_ISO_TOL,
        'edge_floor': TS_CORE_EDGE_FLOOR,
        'max_candidates': TS_CORE_MAX_CANDIDATES,
        'score': score_config(),
    }


def score_config():
    return {
        'formula': 'S = beta * wbo_progress^WBO_PROGRESS_POWER',
        'EVENT_WEIGHT_POWER': EVENT_WEIGHT_POWER,
        'WBO_PROGRESS_POWER': WBO_PROGRESS_POWER,
    }


def _resolve_xtb_threads(threads=None, max_threads=XTB_MAX_THREADS):
    raw = XTB_OMP_THREADS if threads is None else threads
    if isinstance(raw, str) and raw.lower() == "auto":
        requested = int(os.cpu_count() or 1)
    else:
        requested = int(raw)
    return max(1, min(max(1, requested), max(1, int(max_threads))))


def _normal_multiplicity(multiplicity):
    multiplicity = int(multiplicity)
    if multiplicity < 1:
        raise ValueError("multiplicity must be >= 1")
    return multiplicity


def _xtb_charge_uhf(charge=None, multiplicity=None):
    charge = XTB_CHARGE if charge is None else int(charge)
    multiplicity = (
        _normal_multiplicity(XTB_MULTIPLICITY) if multiplicity is None
        else _normal_multiplicity(multiplicity)
    )
    return charge, multiplicity - 1


def _xyz_path(d, include_xtbhess=False):
    d = Path(d)
    if not d.exists():
        return None
    xyzs = sorted(p for p in d.glob("*.xyz") if "xtbhess" not in p.name)
    if not xyzs and include_xtbhess:
        xyzs = sorted(d.glob("*.xyz"))
    return xyzs[0] if xyzs else None


def _xtb_available():
    return shutil.which("xtb") is not None


def _normal_xtb_mode(mode=None):
    mode = (mode or XTB_CACHE_MODE or "auto").lower()
    if mode in {"off", "false", "0", "cache", "cache_only"}:
        return "cache-only"
    if mode not in {"auto", "cache-only"}:
        raise ValueError("BGCP_XTB_MODE must be 'auto' or 'cache-only'")
    return mode


def _need_xtb(kind, path, xtb_mode=None):
    if _normal_xtb_mode(xtb_mode) == "cache-only":
        raise RuntimeError(
            f"missing {kind} cache at {path}; BGCP_XTB_MODE=cache-only")
    if not _xtb_available():
        raise RuntimeError(
            f"missing {kind} cache at {path}, and xtb is not on PATH")


def _ensure_sp_cache(d, label, xyz_fallback=None, xtb_mode=None,
                     charge=None, multiplicity=None):
    """Ensure one single-point cache directory has an XYZ and WBO file."""
    d = Path(d)
    local_xyz = _xyz_path(d)
    fallback_xyz = Path(xyz_fallback) if xyz_fallback else None
    if local_xyz is None and fallback_xyz is not None:
        d.mkdir(parents=True, exist_ok=True)
        name = fallback_xyz.name
        if "xtbhess" in name:
            safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label))
            name = f"{safe}.xyz"
        local_xyz = d / name
        if fallback_xyz.resolve() != local_xyz.resolve():
            shutil.copy(fallback_xyz, local_xyz)
    xyz = local_xyz or fallback_xyz
    wbo = d / "wbo"
    if xyz is None:
        raise RuntimeError(f"missing {label} xyz in {d}")
    if not wbo.exists():
        _need_xtb(f"{label} single-point WBO", wbo, xtb_mode)
        chrg, uhf = _xtb_charge_uhf(charge, multiplicity)
        run_xtb(xyz, d, charge=chrg, uhf=uhf,
                omp_threads=_resolve_xtb_threads())
    return d


def _ensure_hess_cache(hess_dir, label, xyz_fallback=None, xtb_mode=None,
                       charge=None, multiplicity=None):
    """Ensure one Hessian cache directory has g98.out.

    The hessian directory may contain its own XYZ. If not, the matching
    single-point XYZ is used as the source geometry and copied into the hessian
    cache before running xtb.
    """
    hess_dir = Path(hess_dir)
    g98 = hess_dir / "g98.out"
    if not g98.exists():
        hess_dir.mkdir(parents=True, exist_ok=True)
        xyz = _xyz_path(hess_dir, include_xtbhess=True) or (
            Path(xyz_fallback) if xyz_fallback else None)
        if xyz is None:
            raise RuntimeError(f"missing {label} hessian xyz in {hess_dir}")
        _need_xtb(f"{label} Hessian g98.out", g98, xtb_mode)
        chrg, uhf = _xtb_charge_uhf(charge, multiplicity)
        run_xtb_hess(
            xyz, hess_dir, charge=chrg, uhf=uhf,
            omp_threads=_resolve_xtb_threads())
    return parse_g98_modes(g98)


def load(d):
    xyz_path = _xyz_path(d)
    if xyz_path is None:
        raise RuntimeError(f"missing xyz in {d}")
    el, xyz = parse_xyz(xyz_path)
    n = len(el); wbo = np.zeros((n, n))
    for ln in (Path(d) / "wbo").read_text().splitlines():
        p = ln.split()
        if len(p) < 3: continue
        i, j = int(p[0])-1, int(p[1])-1
        wbo[i, j] = float(p[2]); wbo[j, i] = wbo[i, j]
    return el, np.asarray(xyz, float), wbo


def _load_sp(d, label, xyz_fallback=None, xtb_mode=None,
             charge=None, multiplicity=None):
    _ensure_sp_cache(d, label, xyz_fallback=xyz_fallback,
                     xtb_mode=xtb_mode,
                     charge=charge, multiplicity=multiplicity)
    return load(d)


def _load_hess(hess_dir, label, xyz_fallback=None, xtb_mode=None,
               charge=None, multiplicity=None):
    return _ensure_hess_cache(hess_dir, label, xyz_fallback=xyz_fallback,
                              xtb_mode=xtb_mode,
                              charge=charge, multiplicity=multiplicity)


def _iter_labels(sd):
    labels = set()
    if not Path(sd).exists():
        return []
    for d in Path(sd).iterdir():
        if not d.is_dir():
            continue
        m = re.match(r"(?:sp|hess)_iter(\d+)$", d.name)
        if m:
            labels.add(int(m.group(1)))
    return [f"iter{i}" for i in sorted(labels)]


def load_endpoint_from_xyz(xyz_path, workdir, label, *, charge=None,
                           multiplicity=None, xtb_mode=None):
    """Load one molecule from XYZ plus charge/multiplicity-backed WBO cache.

    `workdir` is caller-owned cache space.  If `workdir/wbo` is missing and
    `xtb_mode='auto'`, xtb is run from `xyz_path`; otherwise the existing cache
    is used.  This is the file-level primitive underneath Stage 1 and Stage 2.
    """
    cache = _ensure_sp_cache(
        workdir, label, xyz_fallback=xyz_path, xtb_mode=xtb_mode,
        charge=charge, multiplicity=multiplicity)
    return load(cache)


def step_inputs_from_arrays(step_name, elR, xyzR, wboR, elP, xyzP, wboP,
                            step_dir=None):
    """Build Stage 1 inputs directly from in-memory molecule arrays."""
    return StepInputs(
        step_name=str(step_name),
        step_dir=Path("." if step_dir is None else step_dir),
        elR=list(elR),
        xyzR=np.asarray(xyzR, float),
        wboR=np.asarray(wboR, float),
        elP=list(elP),
        xyzP=np.asarray(xyzP, float),
        wboP=np.asarray(wboP, float),
    )


def alignment_inputs_from_xyz(reactant_xyz, product_xyz, workdir=None, *,
                              name="alignment", charge=None,
                              multiplicity=None, xtb_mode=None,
                              reactant_workdir=None, product_workdir=None,
                              reactant_label="R", product_label="P"):
    """Build Stage 1 inputs from R/P XYZ files.

    This is the preferred file API for arbitrary molecules.  Pass explicit
    `reactant_workdir` and `product_workdir` to avoid any imposed directory
    layout, or pass `workdir` for the conventional `workdir/R` and `workdir/P`
    cache locations.
    """
    if reactant_workdir is None or product_workdir is None:
        if workdir is None:
            raise ValueError(
                "provide either workdir or both reactant_workdir/product_workdir")
        base = Path(workdir)
        reactant_workdir = base / reactant_label
        product_workdir = base / product_label
    elR, xyzR, wboR = load_endpoint_from_xyz(
        reactant_xyz, reactant_workdir, reactant_label,
        charge=charge, multiplicity=multiplicity, xtb_mode=xtb_mode)
    elP, xyzP, wboP = load_endpoint_from_xyz(
        product_xyz, product_workdir, product_label,
        charge=charge, multiplicity=multiplicity, xtb_mode=xtb_mode)
    step_dir = (
        Path(workdir) if workdir is not None
        else Path(reactant_workdir).parent
    )
    return step_inputs_from_arrays(
        name, elR, xyzR, wboR, elP, xyzP, wboP, step_dir=step_dir)


def load_step_inputs(step_name, *, charge=None, multiplicity=None,
                     xtb_mode=None):
    """Adapter: load/cache-fill R and P endpoints for one benchmark step."""
    sd = WORK / step_name
    if not sd.exists():
        raise RuntimeError(f"missing step directory: {sd}")
    return alignment_inputs_from_xyz(
        _xyz_path(sd / "R"), _xyz_path(sd / "P"),
        name=step_name, reactant_workdir=sd / "R", product_workdir=sd / "P",
        charge=charge, multiplicity=multiplicity, xtb_mode=xtb_mode)


def ts_target_from_arrays(kind, label, el, xyz, wbo, freqs, modes,
                          target_index=0):
    """Build a Stage 2 verification target directly from arrays."""
    return TSTarget(
        kind=str(kind),
        target_index=int(target_index),
        label=str(label),
        el=list(el),
        xyz=np.asarray(xyz, float),
        wbo=np.asarray(wbo, float),
        freqs=np.asarray(freqs, float),
        modes=np.asarray(modes, float),
    )


def _safe_cache_name(label):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label)).strip("_") or "target"


def ts_target_from_xyz(kind, label, xyz_path, workdir=None, *,
                       charge=None, multiplicity=None, xtb_mode=None,
                       target_index=0, sp_workdir=None, hess_workdir=None):
    """Build one Stage 2 target from TS/IG/GT XYZ plus charge/multiplicity.

    Pass explicit `sp_workdir` and `hess_workdir` when the caller owns the
    cache layout.  Passing `workdir` is just a convenience wrapper that creates
    `<workdir>/<label>_sp` and `<workdir>/<label>_hess`.
    """
    if sp_workdir is None or hess_workdir is None:
        if workdir is None:
            raise ValueError(
                "provide either workdir or both sp_workdir/hess_workdir")
        base = Path(workdir)
        safe = _safe_cache_name(label)
        sp_workdir = base / f"{safe}_sp"
        hess_workdir = base / f"{safe}_hess"
    el, xyz, wbo = load_endpoint_from_xyz(
        xyz_path, sp_workdir, label,
        charge=charge, multiplicity=multiplicity, xtb_mode=xtb_mode)
    freqs, modes = _load_hess(
        hess_workdir, label, xyz_fallback=_xyz_path(sp_workdir) or xyz_path,
        xtb_mode=xtb_mode, charge=charge, multiplicity=multiplicity)
    return ts_target_from_arrays(
        kind, label, el, xyz, wbo, freqs, modes, target_index=target_index)


def discover_mechanisms_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP,
                                    *, step_name="alignment",
                                    config=None, inner_workers=0):
    """Array-based R-P alignment / mechanism discovery entry point."""
    inputs = step_inputs_from_arrays(
        step_name, elR, xyzR, wboR, elP, xyzP, wboP)
    return run_rp_stage(inputs, config=config, inner_workers=inner_workers)


def discover_mechanisms_from_xyz(reactant_xyz, product_xyz, workdir=None, *,
                                 name="alignment", charge=None,
                                 multiplicity=None, xtb_mode=None,
                                 config=None, inner_workers=0,
                                 return_inputs=False,
                                 reactant_workdir=None,
                                 product_workdir=None):
    """File-based R-P alignment / mechanism discovery entry point."""
    inputs = alignment_inputs_from_xyz(
        reactant_xyz, product_xyz, workdir=workdir, name=name,
        charge=charge, multiplicity=multiplicity, xtb_mode=xtb_mode,
        reactant_workdir=reactant_workdir,
        product_workdir=product_workdir)
    result = run_rp_stage(inputs, config=config, inner_workers=inner_workers)
    return (inputs, result) if return_inputs else result


def load_ts_targets(inputs, include_gt=None):
    """Adapter: load/cache-fill GT and IG targets for one benchmark step.

    GT is optional and controlled separately from IG loading.  IGs are loaded
    from the conventional `sp_iter<N>` / `hess_iter<N>` cache pairs.
    """
    include_gt = INCLUDE_GT if include_gt is None else bool(include_gt)
    sd = inputs.step_dir
    targets = []
    if include_gt:
        gt_sp = sd / "sp_groundtruth"
        gt_hess = sd / "hess_groundtruth"
        xyz = _xyz_path(gt_sp) or _xyz_path(gt_hess, include_xtbhess=True)
        targets.append(ts_target_from_xyz(
            'gt', 'GT', xyz, target_index=-1,
            sp_workdir=gt_sp, hess_workdir=gt_hess))

    for label in _iter_labels(sd):
        hess_dir = sd / f"hess_{label}"
        sp_dir = sd / f"sp_{label}"
        try:
            xyz = _xyz_path(sp_dir) or _xyz_path(
                hess_dir, include_xtbhess=True)
            target = ts_target_from_xyz(
                'ig', label, xyz,
                target_index=len([t for t in targets if t.kind == 'ig']),
                sp_workdir=sp_dir, hess_workdir=hess_dir)
        except Exception:
            continue
        targets.append(target)
    return targets


def _estimate_cut_sweep_units(step_name):
    """Estimate R-P cut-sweep work units for scheduling.

    The actual work is `(no_cut + strong_R_edges) * N_SEEDS_PER_RUN`.  This is
    cheap to compute from the cached R WBO and lets auto mode start large
    steps early instead of discovering them as a slow tail.
    """
    try:
        _, _, wboR = load(WORK / step_name / "R")
    except Exception:
        return 1
    strong_edges = int(np.sum(np.triu(wboR >= CUT_FLOOR, 1)))
    return max(1, (strong_edges + 1) * N_SEEDS_PER_RUN)


def _bond_key(bonds, orbits=None):
    pairs = []
    for a, b in bonds:
        a = int(a); b = int(b)
        if orbits is not None:
            a = int(orbits[a]); b = int(orbits[b])
        if a > b:
            a, b = b, a
        pairs.append((a, b))
    return tuple(sorted(pairs))


def _mechanism_bond_key(mech, r_orbits):
    return (
        _bond_key(mech['broken_bonds_R'], r_orbits),
        _bond_key(mech['formed_bonds_R'], r_orbits),
    )


def _gt_score(mech):
    gt = mech.get('gt')
    return float(gt['S']) if gt and gt.get('S') is not None else float('-inf')


def dedupe_mechanisms_by_bond_changes(mechanisms, r_orbits):
    """Collapse final-view mechanisms with the same R-symmetry bond changes.

    The cut sweep may find multiple concrete alignments whose broken/formed
    R-index bonds differ only by swapping equivalent reactant atoms.  They are
    the same mechanism for the view, so keep the highest-GT-scoring
    representative and retain provenance for the slim JSON / button tooltip.
    """
    groups = {}
    for mech in mechanisms:
        key = _mechanism_bond_key(mech, r_orbits)
        groups.setdefault(key, []).append(mech)

    deduped = []
    for group in groups.values():
        rep = max(group, key=_gt_score)
        rep['dedup_count'] = sum(m.get('dedup_count', 1) for m in group)
        rep['dedup_source_ids'] = [int(m['id']) for m in group]
        rep['dedup_cuts'] = sorted({
            cut
            for m in group
            for cut in m.get('dedup_cuts', [m['cut']])
        })
        deduped.append(rep)

    for new_id, mech in enumerate(deduped, 1):
        suffix = re.sub(r"^#\d+:\s*", "", mech['label'])
        if mech['dedup_count'] > 1:
            suffix = f"{suffix} [dedup x{mech['dedup_count']}]"
        mech['id'] = new_id
        mech['label'] = f"#{new_id}: {suffix}"
    return deduped


def _core_pool_key(mapping, core_R):
    return (
        tuple((int(r), int(mapping[r])) for r in sorted(core_R)),
        (),
    )


def _add_core_pool_entry(pool, core_R, mapping, source, dedup_count=1):
    core_R = tuple(sorted(int(r) for r in core_R))
    mapping = {int(r): int(t) for r, t in dict(mapping).items()
               if int(r) in core_R}
    if len(mapping) != len(core_R):
        return
    sig = _core_pool_key(mapping, core_R)
    entry = pool.get(sig)
    if entry is None:
        pool[sig] = {
            'mapping': mapping,
            'cuts': frozenset(),
            'dedup_count': int(dedup_count),
            'sources': {source},
        }
    else:
        entry['dedup_count'] = entry.get('dedup_count', 1) + int(dedup_count)
        entry.setdefault('sources', set()).add(source)


def _product_core_pool_to_reactant(pool_PT, mapping_RP, core_R):
    """Convert a P->TS core pool into the R->TS frame used by scoring."""
    inv_RP = {int(p): int(r) for r, p in mapping_RP.items()}
    out = {}
    for entry in pool_PT.values():
        pulled = {}
        for p, t in entry.get('mapping', {}).items():
            r = inv_RP.get(int(p))
            if r is not None:
                pulled[r] = int(t)
        _add_core_pool_entry(
            out, core_R, pulled, 'P',
            dedup_count=entry.get('dedup_count', 1))
    return out


def _merge_endpoint_core_pools(core_R, r_pool, p_pool_as_r):
    """Union R->TS and pulled-back P->TS core alternatives.

    The merged object is still an R->TS pool because all scoring and normal-mode
    reindexing are R-frame operations.  `sources` records whether a core mapping
    was seen from the reactant endpoint, product endpoint, or both.
    """
    merged = {}
    for entry in r_pool.values():
        _add_core_pool_entry(
            merged, core_R, entry.get('mapping', {}), 'R',
            dedup_count=entry.get('dedup_count', 1))
    for entry in p_pool_as_r.values():
        _add_core_pool_entry(
            merged, core_R, entry.get('mapping', {}), 'P',
            dedup_count=entry.get('dedup_count', 1))
    return merged


def _pairs_to_product_frame(mapping_RP, pairs_R):
    out = []
    for a, b in pairs_R:
        if a in mapping_RP and b in mapping_RP:
            out.append((int(mapping_RP[a]), int(mapping_RP[b])))
    return out


def _ts_endpoint_pool_task(task):
    """Worker task for one endpoint-side TS core pool.

    A task is one `(target TS, mechanism, endpoint)` alignment.  The caller
    merges the returned R->TS and P->TS pools and then scores them in the main
    process so ranking stays deterministic.
    """
    t0 = time.time()
    pool = ts_core_pool(task['elS'], task['wboS'],
                        task['elT'], task['wboT'],
                        task['core_S'],
                        broken_R=task['broken_S'],
                        formed_R=task['formed_S'],
                        edge_floor=task['edge_floor'],
                        iso_tol=task['iso_tol'],
                        max_candidates=task['max_candidates'])
    return {
        'key': task['key'],
        'target_order': int(task['target_order']),
        'target_label': task['target_label'],
        'mech_id': int(task['mech_id']),
        'mech_pos': int(task['mech_pos']),
        'endpoint': task['endpoint'],
        'pool': pool,
        'n_pool': len(pool),
        'core_size': len(task['core_S']),
        'elapsed': time.time() - t0,
        'hit_cap': len(pool) >= task['max_candidates'],
    }


def _clip01(value):
    return max(0.0, min(1.0, float(value)))


def _event_progress(kind, wR, wP, wT):
    """How far the TS WBO has moved in the detected R->P event direction."""
    delta = abs(float(wP) - float(wR))
    if delta < 1e-12:
        return 1.0
    if kind == 'formed':
        return _clip01((float(wT) - float(wR)) / delta)
    return _clip01((float(wR) - float(wT)) / delta)


def _ranker_event_terms(wboR, wboP, wboT, mapping_RP, mapping_RT,
                        broken_R, formed_R, event_weight_power=1.0):
    """Per-event WBO terms used by the TS ranker.

    Event strength is the detected endpoint change |WBO_P - WBO_R|.  Progress
    is the mapped TS WBO's movement in the same direction: forming events
    require WBO_TS > WBO_R, broken events require WBO_TS < WBO_R.
    """
    terms = []
    for kind, pairs in (('broken', broken_R), ('formed', formed_R)):
        for a, b in pairs:
            if (a not in mapping_RP or b not in mapping_RP or
                    a not in mapping_RT or b not in mapping_RT):
                continue
            pa, pb = int(mapping_RP[a]), int(mapping_RP[b])
            ta, tb = int(mapping_RT[a]), int(mapping_RT[b])
            wR = float(wboR[a, b])
            wP = float(wboP[pa, pb])
            wT = float(wboT[ta, tb])
            delta = abs(wP - wR)
            weight = float(delta ** float(event_weight_power))
            terms.append({
                'kind': kind,
                'R_pair': [int(a), int(b)],
                'P_pair': [pa, pb],
                'T_pair': [ta, tb],
                'wbo_R': wR,
                'wbo_P': wP,
                'wbo_T': wT,
                'delta_wbo_RP': delta,
                'event_weight': weight,
                'ts_progress': _event_progress(kind, wR, wP, wT),
            })
    return terms


def _weighted_wbo_progress(event_terms):
    total = sum(float(t['event_weight']) for t in event_terms)
    if total < 1e-12:
        return 1.0
    return sum(float(t['event_weight']) * float(t['ts_progress'])
               for t in event_terms) / total


def _bond_reaction_vector_from_terms(xyz_TS_in_R, event_terms, *,
                                     weighted=True):
    xyz = np.asarray(xyz_TS_in_R, dtype=float)
    V = np.zeros_like(xyz)
    for term in event_terms:
        i, j = term['R_pair']
        v = xyz[j] - xyz[i]
        n = float(np.linalg.norm(v))
        if n < 1e-9:
            continue
        scale = float(term['event_weight']) if weighted else 1.0
        u = scale * (v / n)
        if term['kind'] == 'broken':
            V[i] -= u
            V[j] += u
        else:
            V[i] += u
            V[j] -= u
    return V


def score_one(elR, xyzR, elT, xyzT, mapping_RT, freqs, modes_TS,
              broken_R, formed_R, core_R, *,
              wboR, wboP, wboT, mapping_RP, score_weights=None):
    weights = dict(score_config() if score_weights is None else score_weights)
    event_weight_power = float(weights.get(
        'EVENT_WEIGHT_POWER', EVENT_WEIGHT_POWER))
    wbo_progress_power = float(weights.get(
        'WBO_PROGRESS_POWER', WBO_PROGRESS_POWER))
    mapping_RT = {int(r): int(t) for r, t in mapping_RT.items()}
    modes_R = reindex_modes_to_R(modes_TS, mapping_RT, len(elR))
    mode_norms = np.linalg.norm(modes_TS.reshape(modes_TS.shape[0], -1), axis=1)
    ts_in_R = np.asarray(xyzR, float).copy()
    for r, t in mapping_RT.items(): ts_in_R[r] = xyzT[t]
    event_terms = _ranker_event_terms(
        wboR, wboP, wboT, mapping_RP, mapping_RT, broken_R, formed_R,
        event_weight_power=event_weight_power)
    V = _bond_reaction_vector_from_terms(ts_in_R, event_terms, weighted=True)
    wbo_progress = _weighted_wbo_progress(event_terms)
    beta = bond_overlap_per_mode(modes_R, V, mode_norms=mode_norms)
    imag = list(np.where(freqs < 0)[0])
    if not imag: return None
    pk = max(imag, key=lambda k: beta[k])
    progress_factor = float(wbo_progress) ** wbo_progress_power

    def target_pairs(r_pairs):
        out = []
        for a, b in r_pairs:
            if a in mapping_RT and b in mapping_RT:
                out.append([int(mapping_RT[a]), int(mapping_RT[b])])
        return out

    if event_terms:
        v_flat = np.asarray(V).reshape(-1)
        v_norm = float(np.linalg.norm(v_flat))
        mode = modes_R[pk]
        mode_norm = float(mode_norms[pk])
        for term in event_terms:
            part = _bond_reaction_vector_from_terms(
                ts_in_R, [term], weighted=True)
            signed = 0.0
            if v_norm > 1e-12 and mode_norm > 1e-12:
                signed = float(mode.reshape(-1) @
                               (part.reshape(-1) / v_norm)) / mode_norm
            term['beta_contribution_signed'] = signed

    return {'S': float(beta[pk] * progress_factor),
            'beta': float(beta[pk]),
            'wbo_progress': float(wbo_progress),
            'wbo_progress_factor': progress_factor,
            'freq': float(freqs[pk]), 'k': int(pk),
            'score_formula': weights.get(
                'formula',
                'S = beta * wbo_progress^WBO_PROGRESS_POWER'),
            'event_terms': event_terms,
            'core_map': {str(int(r)): int(mapping_RT[r])
                         for r in core_R if r in mapping_RT},
            # Viewer fields stay in native target indexing.  Only annotations
            # are translated from R-core indices to target indices.
            'elements': list(elT),
            'xyz': np.asarray(xyzT, float).tolist(),
            'picked_disp': np.asarray(modes_TS[pk], float).tolist(),
            'core_atoms_T': [int(mapping_RT[r]) for r in core_R if r in mapping_RT],
            'broken_bonds_T': target_pairs(broken_R),
            'formed_bonds_T': target_pairs(formed_R),
            # Keep the R-indexed core-only materialization available for
            # debugging/scoring audits; the viewer no longer renders it.
            'xyz_in_R': ts_in_R.tolist(),
            'picked_disp_R': modes_R[pk].tolist()}


def best_under_mech_using_pool(elR, xyzR, wboR, wboP,
                                 elT, xyzT, wboT, freqs, modes_TS,
                                 rt_pool, mapping_RP,
                                 broken_R, formed_R, core_R,
                                 score_weights=None):
    """Score every R-frame core witness under one mech.

    `rt_pool` may contain native R->TS candidates and P->TS candidates that
    were pulled back through the R-P mechanism witness.  Two witnesses that
    agree on `core_R -> TS` are score-equivalent for this mechanism, so score
    one representative per exact core map and keep the highest S.
    """
    core_R_set = frozenset(core_R)
    seen_core = set()
    best = None
    for v in rt_pool.values():
        witness = v['mapping']
        # Per-mechanism core-restricted key
        core_key = frozenset((r, witness[r]) for r in core_R_set if r in witness)
        if core_key in seen_core:
            continue
        seen_core.add(core_key)
        mapping_RT = dict(witness)
        s = score_one(elR, xyzR, elT, xyzT, mapping_RT, freqs, modes_TS,
                      broken_R, formed_R, core_R,
                      wboR=wboR, wboP=wboP, wboT=wboT,
                      mapping_RP=mapping_RP,
                      score_weights=score_weights)
        if s:
            s['core_sources'] = sorted(v.get('sources', {'R'}))
            s['core_pool_dedup_count'] = int(v.get('dedup_count', 1))
        if s and (best is None or s['S'] > best['S']):
            best = s
    return best


def _json_ready(value):
    """Convert numpy/path/set-heavy stage records to plain JSON values."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return [_json_ready(v) for v in sorted(value)]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def write_stage_json(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_json_ready(record), indent=2))
    return path


def read_stage_json(path):
    return json.loads(Path(path).read_text())


def _int_mapping(mapping):
    return {int(k): int(v) for k, v in dict(mapping).items()}


def _int_pairs(pairs):
    return [(int(a), int(b)) for a, b in pairs]


def _mechanism_state(inputs, mech):
    mapping_RP = _int_mapping(mech['mapping_RP'])
    broken_R = _int_pairs(mech.get('broken_bonds_R', ()))
    formed_R = _int_pairs(mech.get('formed_bonds_R', ()))
    core_R = [int(r) for r in mech.get('core_atoms', ())]
    return mapping_RP, broken_R, formed_R, core_R


def _mechanism_for_view(mech):
    """Return a mechanism dict without private runtime-only fields."""
    out = dict(mech)
    out.pop('_state', None)
    return out


def _safe_file_stem(value):
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return stem or "item"


def write_rp_alignment_files(inputs, rp_result, out_dir=None):
    """Write clean Stage 1 alignment files for NEB/path-building use.

    Coordinates are only reindexed into the R atom frame.  No Kabsch or other
    spatial fitting is applied.  Each mechanism directory is self-contained and
    contains an R endpoint, the mechanism-specific aligned P endpoint, a
    two-frame XYZ, a mapping CSV, and a JSON metadata file.
    """
    out_dir = Path(out_dir) if out_dir is not None else alignment_output_dir(
        inputs.step_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    mech_root = out_dir / "mechanisms"
    mech_root.mkdir(exist_ok=True)

    (out_dir / "R.xyz").write_text(write_xyz_str(
        inputs.elR, inputs.xyzR, f"{inputs.step_name} reactant R frame"))
    (out_dir / "P_original.xyz").write_text(write_xyz_str(
        inputs.elP, inputs.xyzP, f"{inputs.step_name} product original order"))

    manifest_mechs = []
    for mech in rp_result.get('mechanisms', []):
        mech_id = int(mech['id'])
        name = f"mechanism_{mech_id:03d}"
        mdir = mech_root / name
        mdir.mkdir(exist_ok=True)
        mapping_RP = _int_mapping(mech['mapping_RP'])
        p_aligned = np.asarray(mech.get('product_xyz_in_R'), float)
        if p_aligned.shape != np.asarray(inputs.xyzR).shape:
            p_aligned = np.asarray(inputs.xyzR, float).copy()
            for i_R, i_P in mapping_RP.items():
                p_aligned[i_R] = inputs.xyzP[i_P]

        r_xyz = write_xyz_str(
            inputs.elR, inputs.xyzR,
            f"{inputs.step_name} {name} R endpoint")
        p_xyz = write_xyz_str(
            inputs.elR, p_aligned,
            f"{inputs.step_name} {name} P aligned to R atom order")
        (mdir / "R.xyz").write_text(r_xyz)
        (mdir / "P_aligned.xyz").write_text(p_xyz)
        (mdir / "neb_endpoints.xyz").write_text(r_xyz + p_xyz)

        rows = ["R_index,R_element,P_index,P_element"]
        for r in sorted(mapping_RP):
            p = mapping_RP[r]
            rows.append(f"{r},{inputs.elR[r]},{p},{inputs.elP[p]}")
        (mdir / "mapping_R_to_P.csv").write_text("\n".join(rows) + "\n")

        meta = {
            'step': inputs.step_name,
            'mechanism_id': mech_id,
            'label': mech.get('label'),
            'cut': mech.get('cut'),
            'mapping_RP': mapping_RP,
            'broken_bonds_R': mech.get('broken_bonds_R', []),
            'formed_bonds_R': mech.get('formed_bonds_R', []),
            'formed_bonds_P': mech.get('formed_bonds_P', []),
            'core_atoms_R': mech.get('core_atoms', []),
            'dedup_count': mech.get('dedup_count', 1),
            'dedup_source_ids': mech.get('dedup_source_ids', [mech_id]),
            'dedup_cuts': mech.get('dedup_cuts', [mech.get('cut')]),
            'files': {
                'reactant': 'R.xyz',
                'product_aligned': 'P_aligned.xyz',
                'neb_endpoints': 'neb_endpoints.xyz',
                'mapping_csv': 'mapping_R_to_P.csv',
            },
            'coordinate_policy': (
                'P_aligned.xyz is product geometry reindexed into R atom '
                'order using mapping_RP; no spatial/Kabsch alignment is used.'
            ),
        }
        write_stage_json(mdir / "mechanism.json", meta)
        manifest_mechs.append({
            'id': mech_id,
            'label': mech.get('label'),
            'directory': f"mechanisms/{name}",
            'reactant': f"mechanisms/{name}/R.xyz",
            'product_aligned': f"mechanisms/{name}/P_aligned.xyz",
            'neb_endpoints': f"mechanisms/{name}/neb_endpoints.xyz",
            'metadata': f"mechanisms/{name}/mechanism.json",
            'mapping_csv': f"mechanisms/{name}/mapping_R_to_P.csv",
            'broken_bonds_R': mech.get('broken_bonds_R', []),
            'formed_bonds_R': mech.get('formed_bonds_R', []),
            'core_atoms_R': mech.get('core_atoms', []),
        })

    manifest = {
        'stage': 'rp_alignment_files',
        'step': inputs.step_name,
        'n_atoms': len(inputs.elR),
        'source_stage': 'rp_stage.json',
        'root_files': {
            'reactant': 'R.xyz',
            'product_original_order': 'P_original.xyz',
        },
        'mechanisms': manifest_mechs,
        'coordinate_policy': (
            'All mechanism P endpoints are reindexed into the R atom frame. '
            'The files are intended for downstream path/NEB setup that needs '
            'matching atom order.'
        ),
    }
    write_stage_json(out_dir / "manifest.json", manifest)
    return {
        'stage': 'rp_alignment_files',
        'step': inputs.step_name,
        'out_dir': str(out_dir),
        'manifest': str(out_dir / "manifest.json"),
        'n_mechanisms': len(manifest_mechs),
    }


def _ts_score_records_for_export(mech):
    records = []
    if mech.get('gt') and mech['gt'].get('S') is not None:
        records.append(('gt', 'GT', mech['gt']))
    for item in mech.get('igs', []):
        if item.get('S') is not None:
            records.append(('ig', item.get('label', 'target'), item))
    return records


def write_ts_alignment_files(inputs, ts_result, out_dir=None):
    """Write selected best-S TS/IG/GT alignment files from Stage 2.

    Stage 2 scores mechanism-local core mappings.  The exported R-frame file is
    therefore a core-aligned materialization: mapped core atoms use the selected
    target coordinates, while unmapped spectator atoms remain at the reactant
    endpoint coordinates.  The native target coordinates are exported alongside
    it so downstream users can choose the representation they need.
    """
    out_dir = Path(out_dir) if out_dir is not None else (
        alignment_output_dir(inputs.step_name) / "ts_alignments")
    out_dir.mkdir(parents=True, exist_ok=True)
    mech_root = out_dir / "mechanisms"
    mech_root.mkdir(exist_ok=True)

    manifest_mechs = []
    for mech in ts_result.get('mechanisms', []):
        mech_id = int(mech['id'])
        mech_name = f"mechanism_{mech_id:03d}"
        mech_dir = mech_root / mech_name
        mech_dir.mkdir(exist_ok=True)
        target_entries = []

        for kind, label, score in _ts_score_records_for_export(mech):
            target_name = f"{kind}_{_safe_file_stem(label)}"
            tdir = mech_dir / target_name
            tdir.mkdir(exist_ok=True)
            native_elements = score.get('elements') or inputs.elR
            native_xyz = score.get('xyz')
            r_xyz = score.get('xyz_in_R')
            disp_R = score.get('picked_disp_R')

            if native_xyz is not None:
                (tdir / "TS_native.xyz").write_text(write_xyz_str(
                    native_elements, native_xyz,
                    f"{inputs.step_name} {mech_name} {label} native target"))
            if r_xyz is not None:
                (tdir / "TS_core_aligned_R_frame.xyz").write_text(
                    write_xyz_str(
                        inputs.elR, r_xyz,
                        f"{inputs.step_name} {mech_name} {label} "
                        "best-S core-aligned R frame"))
            if r_xyz is not None and disp_R is not None:
                (tdir / "picked_mode_R_frame.xyz").write_text(
                    xyz_with_disp(
                        inputs.elR, r_xyz, disp_R,
                        f"{inputs.step_name} {mech_name} {label} "
                        f"mode {score.get('k')} freq={score.get('freq')}"))

            score_doc = {
                'step': inputs.step_name,
                'mechanism_id': mech_id,
                'mechanism_label': mech.get('label'),
                'kind': kind,
                'label': label,
                'S': score.get('S'),
                'decomposition': {
                    'beta': score.get('beta'),
                    'wbo_progress': score.get('wbo_progress'),
                    'wbo_progress_factor': score.get('wbo_progress_factor'),
                    'freq': score.get('freq'),
                    'mode_index': score.get('k'),
                },
                'score_formula': score.get('score_formula'),
                'event_terms': score.get('event_terms', []),
                'core_map_R_to_target': score.get('core_map', {}),
                'core_sources': score.get('core_sources', []),
                'core_pool_dedup_count': score.get('core_pool_dedup_count'),
                'broken_bonds_R': mech.get('broken_bonds_R', []),
                'formed_bonds_R': mech.get('formed_bonds_R', []),
                'broken_bonds_T': score.get('broken_bonds_T', []),
                'formed_bonds_T': score.get('formed_bonds_T', []),
                'files': {
                    'target_native': 'TS_native.xyz',
                    'core_aligned_R_frame': 'TS_core_aligned_R_frame.xyz',
                    'picked_mode_R_frame': 'picked_mode_R_frame.xyz',
                },
                'coordinate_policy': (
                    'TS_core_aligned_R_frame.xyz is the selected best-S core '
                    'mapping materialized in R atom order. Only mapped core '
                    'atoms are replaced by target coordinates; spectators '
                    'remain at the reactant endpoint because Stage 2 does not '
                    'enumerate spectator bijections.'
                ),
            }
            write_stage_json(tdir / "score.json", score_doc)
            target_entries.append({
                'kind': kind,
                'label': label,
                'S': score.get('S'),
                'directory': f"mechanisms/{mech_name}/{target_name}",
                'target_native': (
                    f"mechanisms/{mech_name}/{target_name}/TS_native.xyz"),
                'core_aligned_R_frame': (
                    f"mechanisms/{mech_name}/{target_name}/"
                    "TS_core_aligned_R_frame.xyz"),
                'picked_mode_R_frame': (
                    f"mechanisms/{mech_name}/{target_name}/"
                    "picked_mode_R_frame.xyz"),
                'score_json': (
                    f"mechanisms/{mech_name}/{target_name}/score.json"),
            })

        manifest_mechs.append({
            'id': mech_id,
            'label': mech.get('label'),
            'directory': f"mechanisms/{mech_name}",
            'targets': target_entries,
        })

    manifest = {
        'stage': 'ts_alignment_files',
        'step': inputs.step_name,
        'n_atoms': len(inputs.elR),
        'source_stage': 'ts_stage.json',
        'mechanisms': manifest_mechs,
        'coordinate_policy': (
            'Stage 2 exports native target coordinates plus selected best-S '
            'core-aligned R-frame materializations. These are score/debug '
            'artifacts, not full spectator atom bijections.'
        ),
    }
    write_stage_json(out_dir / "manifest.json", manifest)
    return {
        'stage': 'ts_alignment_files',
        'step': inputs.step_name,
        'out_dir': str(out_dir),
        'manifest': str(out_dir / "manifest.json"),
        'n_mechanisms': len(manifest_mechs),
        'n_targets': sum(len(m['targets']) for m in manifest_mechs),
    }


def run_rp_stage(inputs, config=None, inner_workers=0):
    """Stage 1: discover mechanism-dependent R-P alignments.

    This is the reusable alignment/mechanism-discovery entry point.  It runs
    the no-cut plus one-edge cut sweep on R->P, selects minimum bond-change
    mechanisms, classifies broken/forming bonds, and stores the R-P witness
    needed for later TS verification.
    """
    cfg = dict(rp_stage_config() if config is None else config)
    t0 = time.time()
    pool = cut_sweep(
        inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
        n_workers=inner_workers,
        cut_floor=cfg.get('cut_floor', CUT_FLOOR),
        graph_floor=cfg.get('graph_floor', 0.2),
        iso_tol=cfg.get('iso_tol', VIEW_ISO_TOL),
        dwbo_threshold=cfg.get('dwbo_threshold', DWBO_THRESHOLD),
        metal_dwbo_threshold=cfg.get(
            'metal_dwbo_threshold', METAL_DWBO_THRESHOLD),
        symmetry_wbo_tol=cfg.get('symmetry_wbo_tol', SYMMETRY_WBO_TOL),
        n_seeds=cfg.get('n_seeds', N_SEEDS_PER_RUN),
        max_branches=cfg.get('max_branches', VIEW_MAX_BRANCHES),
        chunksize=cfg.get('chunksize', CUTSWEEP_CHUNKSIZE),
        symmetry_repair=cfg.get('symmetry_repair', SYMMETRY_REPAIR),
        symmetry_repair_min_changes=cfg.get(
            'symmetry_repair_min_changes', SYMMETRY_REPAIR_MIN_CHANGES),
        symmetry_repair_max_evals=cfg.get(
            'symmetry_repair_max_evals', SYMMETRY_REPAIR_MAX_EVALS),
    )
    if BGCP_TIMING:
        print(f"    {inputs.step_name} {'R-P':>12s} cut_sweep: "
              f"{len(pool):>4d} sigs in {time.time()-t0:.1f}s",
              flush=True)
    rp_min = select_min_mechanisms(pool)
    if not rp_min:
        raise RuntimeError("no min-bond mechanism")

    mechanisms = []
    for mi, (_sig, info) in enumerate(rp_min.items(), 1):
        mapping_RP = _int_mapping(info['mapping'])
        inv_RP = {v: k for k, v in mapping_RP.items()}
        broken, formed, _, _ = classify_bonds(
            mapping_RP, inputs.wboR, inputs.wboP,
            dwbo_threshold=cfg.get('dwbo_threshold', DWBO_THRESHOLD),
            elements_R=inputs.elR,
            elements_P=inputs.elP,
            metal_dwbo_threshold=cfg.get(
                'metal_dwbo_threshold', METAL_DWBO_THRESHOLD))
        broken_R = [(int(a), int(b)) for (a, b, _, _) in broken]
        formed_R = [(int(inv_RP[a]), int(inv_RP[b]))
                    for (a, b, _, _) in formed
                    if a in inv_RP and b in inv_RP]
        core_R = list(core_atoms_in_R_frame(mapping_RP, broken, formed))
        xyzP_in_R = np.asarray(inputs.xyzR, float).copy()
        for i_R, i_P in mapping_RP.items():
            xyzP_in_R[i_R] = inputs.xyzP[i_P]
        has_no_cut = bool(info.get('has_no_cut', False))
        cut = None if has_no_cut else next(iter(info['cuts']), None)
        cut_name = (
            f"{inputs.elR[cut[0]]}{cut[0]}-{inputs.elR[cut[1]]}{cut[1]}"
            if cut else "none"
        )
        br_label = ",".join(
            f"{inputs.elR[a]}{a}-{inputs.elR[b]}{b}" for a, b in broken_R)
        mechanisms.append({
            'id': mi,
            'cut': cut_name,
            'label': f"#{mi}: {br_label} (cut: {cut_name})",
            'dedup_count': int(info.get('dedup_count', 1)),
            'dedup_cuts': (["none"] if has_no_cut else []) + [
                f"{inputs.elR[a]}{a}-{inputs.elR[b]}{b}"
                for a, b in sorted(info['cuts'])
            ] or [cut_name],
            'mapping_RP': mapping_RP,
            'broken_bonds_R': broken_R,
            'formed_bonds_R': formed_R,
            'formed_bonds_P': [[int(a), int(b)] for (a, b, _, _) in formed],
            'core_atoms': [int(r) for r in core_R],
            'product_xyz_in_R': xyzP_in_R.tolist(),
        })

    r_orbits = _nauty_orbits(
        build_graph(inputs.elR, inputs.wboR, bond_cut=0.2),
        wbo_tol=cfg.get('symmetry_wbo_tol', SYMMETRY_WBO_TOL))
    mechanisms = dedupe_mechanisms_by_bond_changes(mechanisms, r_orbits)
    return {
        'stage': 'rp',
        'step': inputs.step_name,
        'n_atoms': len(inputs.elR),
        'config': cfg,
        'mechanisms': [_mechanism_for_view(m) for m in mechanisms],
        'timing': {'rp_seconds': time.time() - t0},
    }


def _add_ts_endpoint_tasks(tasks, inputs, key, target_order, target_label,
                           mech_pos, mech, target, config):
    mapping_RP, br_R, fm_R, core_R = _mechanism_state(inputs, mech)
    common = {
        'key': key,
        'target_order': target_order,
        'target_label': target_label,
        'mech_id': int(mech['id']),
        'mech_pos': int(mech_pos),
        'elT': target.el,
        'wboT': target.wbo,
        'edge_floor': config.get('edge_floor', TS_CORE_EDGE_FLOOR),
        'iso_tol': config.get('iso_tol', VIEW_ISO_TOL),
        'max_candidates': config.get(
            'max_candidates', TS_CORE_MAX_CANDIDATES),
    }
    tasks.append({
        **common,
        'endpoint': 'R',
        'elS': inputs.elR,
        'wboS': inputs.wboR,
        'core_S': list(core_R),
        'broken_S': br_R,
        'formed_S': fm_R,
    })
    core_P = [int(mapping_RP[r]) for r in core_R if r in mapping_RP]
    tasks.append({
        **common,
        'endpoint': 'P',
        'elS': inputs.elP,
        'wboS': inputs.wboP,
        'core_S': core_P,
        'broken_S': _pairs_to_product_frame(mapping_RP, br_R),
        'formed_S': _pairs_to_product_frame(mapping_RP, fm_R),
    })


def run_ts_endpoint_tasks(tasks, inner_workers=0):
    """Run the Stage 2 endpoint core-matching task list."""
    ts_workers = max(1, int(inner_workers or 1))
    if not tasks:
        return []
    if ts_workers <= 1 or len(tasks) == 1:
        return [_ts_endpoint_pool_task(t) for t in tasks]
    with cf.ProcessPoolExecutor(max_workers=min(ts_workers, len(tasks))) as ex:
        futs = [ex.submit(_ts_endpoint_pool_task, t) for t in tasks]
        return [f.result() for f in cf.as_completed(futs)]


def _selected_mechanisms(rp_result, mechanism_ids=None):
    mechanisms = [dict(m) for m in rp_result.get('mechanisms', [])]
    if mechanism_ids is None:
        return mechanisms
    wanted = {int(i) for i in mechanism_ids}
    return [m for m in mechanisms if int(m['id']) in wanted]


def run_ts_stage(inputs, rp_result, targets, config=None, inner_workers=0,
                 mechanism_ids=None):
    """Stage 2: verify GT/IG/TS targets under selected mechanisms."""
    cfg = dict(ts_stage_config() if config is None else config)
    mechanisms = _selected_mechanisms(rp_result, mechanism_ids)
    for mech in mechanisms:
        mech['gt'] = None
        mech['igs'] = []

    endpoint_tasks = []
    score_contexts = []
    target_order = 0

    for target in targets:
        if target.kind != 'gt':
            ig_index = len(mechanisms[0]['igs']) if mechanisms else 0
            for mech in mechanisms:
                mech['igs'].append({
                    'label': target.label,
                    'elements': list(target.el),
                    'xyz': np.asarray(target.xyz, float).tolist(),
                })
        else:
            ig_index = -1
        order = target_order
        target_order += 1
        for mech_pos, mech in enumerate(mechanisms):
            key = (target.kind, int(ig_index if target.kind != 'gt'
                                    else target.target_index), int(mech_pos))
            display_label = f"{target.label}:m{mech['id']}"
            _add_ts_endpoint_tasks(endpoint_tasks, inputs, key, order,
                                   display_label, mech_pos, mech, target, cfg)
            score_contexts.append({
                'key': key,
                'kind': target.kind,
                'target_index': int(ig_index if target.kind != 'gt'
                                    else target.target_index),
                'target_order': order,
                'display_label': display_label,
                'mech_pos': mech_pos,
                'target': target,
            })

    endpoint_results = run_ts_endpoint_tasks(
        endpoint_tasks, inner_workers=inner_workers)
    endpoint_by_key = {}
    for res in endpoint_results:
        endpoint_by_key.setdefault(res['key'], {})[res['endpoint']] = res

    if BGCP_TIMING:
        for res in sorted(endpoint_results,
                          key=lambda x: (x['target_order'],
                                         x['mech_id'],
                                         x['endpoint'])):
            if res['hit_cap']:
                print(f"    [warn] TS core pool hit cap={cfg.get('max_candidates')} "
                      f"{res['target_label']}:{res['endpoint']} "
                      f"core={res['core_size']}",
                      flush=True)
            print(f"    {inputs.step_name} {res['target_label'] + ':' + res['endpoint']:>12s} "
                  f"core_match: {res['n_pool']:>4d} sigs "
                  f"core={res['core_size']} in {res['elapsed']:.1f}s",
                  flush=True)

    for ctx in sorted(score_contexts,
                      key=lambda x: (x['target_order'], x['mech_pos'])):
        mech = mechanisms[ctx['mech_pos']]
        target = ctx['target']
        mapping_RP, br_R, fm_R, core_R = _mechanism_state(inputs, mech)
        parts = endpoint_by_key.get(ctx['key'], {})
        r_pool = parts.get('R', {}).get('pool', {})
        p_pool_native = parts.get('P', {}).get('pool', {})
        p_pool_as_r = _product_core_pool_to_reactant(
            p_pool_native, mapping_RP, core_R)
        merged = _merge_endpoint_core_pools(core_R, r_pool, p_pool_as_r)
        if BGCP_TIMING:
            print(f"    {inputs.step_name} {ctx['display_label']:>12s} core_union: "
                  f"R={len(r_pool)} P={len(p_pool_native)} "
                  f"merged={len(merged)}",
                  flush=True)
        s = best_under_mech_using_pool(
            inputs.elR, inputs.xyzR, inputs.wboR, inputs.wboP,
            target.el, target.xyz, target.wbo, target.freqs, target.modes,
            merged, mapping_RP, br_R, fm_R, core_R,
            score_weights=cfg.get('score', score_config()))
        if ctx['kind'] == 'gt':
            mech['gt'] = s
        elif s:
            mech['igs'][ctx['target_index']].update(s)

    union_top = set()
    for mech in mechanisms:
        ranked = sorted(
            [(i, ig) for i, ig in enumerate(mech['igs'])
             if ig.get('S') is not None],
            key=lambda x: -x[1]['S'])
        top2 = {i for i, _ in ranked[:2]}
        for i, ig in enumerate(mech['igs']):
            ig['is_top2'] = (i in top2)
            if i in top2:
                union_top.add(ig['label'])
    for mech in mechanisms:
        for ig in mech['igs']:
            ig['is_union_top'] = ig['label'] in union_top

    return {
        'stage': 'ts',
        'step': inputs.step_name,
        'config': cfg,
        'mechanisms': [_mechanism_for_view(m) for m in mechanisms],
        'endpoint_results': [
            {k: v for k, v in res.items() if k != 'pool'}
            for res in endpoint_results
        ],
        'union_top_labels': sorted(union_top),
    }


def build_view_data(inputs, rp_result, ts_result=None, include_gt=None):
    """Build the JSON object consumed by the HTML viewer."""
    include_gt = INCLUDE_GT if include_gt is None else bool(include_gt)
    view_score_config = score_config()
    if ts_result is not None:
        view_score_config = (
            ts_result.get('config', {}).get('score') or view_score_config)
    mechanisms = (
        [dict(m) for m in ts_result.get('mechanisms', [])]
        if ts_result is not None
        else [dict(m, gt=None, igs=[]) for m in rp_result.get('mechanisms', [])]
    )
    if include_gt and any(m.get('gt') for m in mechanisms):
        default_id = max(
            mechanisms, key=lambda m: m['gt']['S'] if m.get('gt') else 0
        )['id']
    else:
        default_id = mechanisms[0]['id'] if mechanisms else None
    return {
        'step': inputs.step_name,
        'n_atoms': len(inputs.elR),
        'reactant': {
            'elements': inputs.elR,
            'coords': np.asarray(inputs.xyzR).tolist(),
        },
        'product': {
            'elements': inputs.elP,
            'coords': np.asarray(inputs.xyzP).tolist(),
        },
        'mechanisms': mechanisms,
        'default_mech_id': default_id,
        'include_gt': bool(include_gt),
        'score_config': view_score_config,
    }


def build_eval_slim(view_data):
    slim = {
        'step': view_data['step'],
        'n_atoms': view_data['n_atoms'],
        'n_mechs': len(view_data.get('mechanisms', [])),
        'include_gt': bool(view_data.get('include_gt')),
        'score_config': view_data.get('score_config', score_config()),
        'mechanisms': [],
    }
    for mech in view_data.get('mechanisms', []):
        slim['mechanisms'].append({
            'id': mech['id'],
            'cut': mech['cut'],
            'dedup_count': mech.get('dedup_count', 1),
            'dedup_source_ids': mech.get('dedup_source_ids', [mech['id']]),
            'dedup_cuts': mech.get('dedup_cuts', [mech['cut']]),
            'broken_R': mech['broken_bonds_R'],
            'formed_R': mech['formed_bonds_R'],
            'core_R': mech['core_atoms'],
            'gt': ({k: mech['gt'].get(k) for k in [
                'S', 'beta', 'wbo_progress', 'wbo_progress_factor',
                'freq', 'core_map', 'core_sources',
            ]} if mech.get('gt') else None),
            'igs': [
                {k: ig.get(k) for k in [
                    'label', 'S', 'beta', 'wbo_progress',
                    'wbo_progress_factor', 'freq', 'core_map',
                    'core_sources', 'is_top2',
                ]}
                for ig in mech.get('igs', [])
            ],
        })
    return slim


def write_view_stage(inputs, rp_result, ts_result=None, out_root=None,
                     include_gt=None):
    """Stage 3: write viewer and slim eval artifacts."""
    data = build_view_data(inputs, rp_result, ts_result, include_gt=include_gt)
    slim = build_eval_slim(data)
    paths = pipeline_stage_paths(inputs.step_name, out_root=out_root)
    paths.view_dir.mkdir(parents=True, exist_ok=True)
    paths.view_html.write_text(HTML.format(
        title=f"BGCP &mdash; {inputs.step_name}  "
              f"({len(data['mechanisms'])} mechanisms)",
        data_json=json.dumps(_json_ready(data)),
    ))
    paths.eval_slim_json.write_text(json.dumps(_json_ready(slim)))
    return {
        'stage': 'view',
        'step': inputs.step_name,
        'view_html': str(paths.view_html),
        'eval_slim_json': str(paths.eval_slim_json),
        'slim': slim,
    }


HTML = r"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/jszip@3.10.1/dist/jszip.min.js"></script>
<style>
html,body{{margin:0;padding:0;font-family:-apple-system,sans-serif;background:#fafafa}}
body{{padding:14px;box-sizing:border-box}}
h2{{margin:0 0 4px;font-size:18px}}
.topbar{{display:flex;gap:10px;align-items:stretch;margin-bottom:14px}}
.mech-sel{{flex:1;background:white;padding:10px;border:1px solid #ccc;border-radius:6px}}
.mech-sel button{{padding:6px 12px;margin-right:6px;border:1px solid #aaa;background:#f0f0f0;border-radius:4px;cursor:pointer;font-family:ui-monospace,monospace;font-size:12px}}
.mech-sel button.active{{background:#ffd700;border-color:#a90;font-weight:600}}
.download-all{{padding:0 16px;border:1px solid #888;background:white;border-radius:6px;cursor:pointer;font-family:ui-monospace,monospace;font-size:12px;color:#024}}
.download-all:hover{{background:#eef}}
.index-toggle{{display:flex;align-items:center;gap:6px;padding:0 12px;border:1px solid #aaa;background:white;border-radius:6px;font-family:ui-monospace,monospace;font-size:12px;color:#024;white-space:nowrap}}
.index-toggle input{{margin:0}}
.ref-row{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-bottom:18px}}
.ref-row.no-gt{{grid-template-columns:repeat(2,1fr)}}
.ig-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}}
.panel{{background:white;border:1px solid #ddd;border-radius:6px;padding:6px 8px 8px}}
.panel.top2{{border:2px solid #d4af37}}
.panel.union{{box-shadow:0 0 0 2px #ff9}}
.ph{{display:flex;justify-content:space-between;align-items:baseline;font-size:12px;margin-bottom:4px}}
.ph .lbl{{font-weight:600;font-size:13px}}
.ph .rk{{font-family:ui-monospace,monospace;color:#024}}
.dl{{padding:3px 7px;margin-left:6px;border:1px solid #aaa;background:#fff;border-radius:4px;cursor:pointer;font-family:ui-monospace,monospace;font-size:11px;color:#024}}
.dl:hover{{background:#eef}}
.vw{{position:relative;width:100%;height:230px}}
.ref-row .vw{{height:300px}}
.vwbox{{position:absolute;inset:0}}
.meta{{font-family:ui-monospace,monospace;font-size:11px;color:#444;padding:3px 0 0;line-height:1.4}}
.meta b{{color:#024}}
</style></head><body>
<h2>{title}</h2>
<div class="topbar">
  <div class="mech-sel" id="mech-sel"></div>
  <label class="index-toggle"><input type="checkbox" id="showAtomIndices"> Aligned indices</label>
  <button class="download-all" id="downloadAllBtn">Download</button>
</div>
<div class="ref-row" id="ref-row">
  <div class="panel"><div class="ph"><span class="lbl">Reactant</span><span class="rk">static <button class="dl" onclick="downloadR()">XYZ</button></span></div>
    <div class="vw"><div id="vw_R" class="vwbox"></div></div></div>
  <div class="panel"><div class="ph"><span class="lbl">Product</span><span class="rk" id="prod_label">static <button class="dl" onclick="downloadP()">XYZ</button></span></div>
    <div class="vw"><div id="vw_P" class="vwbox"></div></div></div>
  <div class="panel" id="gt_panel"><div class="ph"><span class="lbl">Ground-truth TS</span><span class="rk" id="gt_S">S=? <button class="dl" onclick="downloadGT()">XYZ</button></span></div>
    <div class="vw"><div id="vw_GT" class="vwbox"></div></div>
    <div class="meta" id="gt_meta"></div></div>
</div>
<div class="ig-grid" id="grid"></div>
<script>
const DATA = {data_json};
let currentMechId = DATA.default_mech_id;
let showAtomIndices = false;
const elements = DATA.reactant.elements;
const xyzR_static = DATA.reactant.coords;
function findMech(id) {{ return DATA.mechanisms.find(m=>m.id===id); }}
function xyzText(els, xyz, comment) {{ let s = els.length+"\n"+(comment||"frame")+"\n"; for (let i=0;i<xyz.length;i++) s += els[i]+"  "+xyz[i][0].toFixed(6)+"  "+xyz[i][1].toFixed(6)+"  "+xyz[i][2].toFixed(6)+"\n"; return s; }}
function buildBody(els, xyz) {{ return xyzText(els, xyz, "frame"); }}
function buildBodyAt(els, xyz, disp, scale) {{ let out = []; for (let i=0;i<xyz.length;i++) {{ const x=xyz[i][0]+scale*disp[i][0], y=xyz[i][1]+scale*disp[i][1], z=xyz[i][2]+scale*disp[i][2]; out.push([x,y,z]); }} return xyzText(els, out, "frame"); }}
function xyzAt(xyz, disp, scale) {{ return xyz.map((p,i)=>[p[0]+scale*disp[i][0], p[1]+scale*disp[i][1], p[2]+scale*disp[i][2]]); }}
function safeName(s) {{ return String(s).replace(/[^A-Za-z0-9_.-]+/g, "_").replace(/^_+|_+$/g, "") || "item"; }}
function downloadBlob(name, blob) {{ const a = document.createElement("a"); const url = URL.createObjectURL(blob); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(url), 1000); }}
function downloadText(name, text) {{ downloadBlob(name, new Blob([text], {{type:"chemical/x-xyz"}})); }}
function downloadXYZ(name, els, xyz, comment) {{ if (!els || !xyz) return; downloadText(name, xyzText(els, xyz, comment)); }}
function downloadR() {{ downloadXYZ(safeName(DATA.step)+"_R.xyz", DATA.reactant.elements, DATA.reactant.coords, DATA.step+" R"); }}
function downloadP() {{ const mech = findMech(currentMechId); if (mech && mech.product_xyz_in_R) downloadXYZ(safeName(DATA.step)+"_P_aligned_mech"+mech.id+".xyz", DATA.reactant.elements, mech.product_xyz_in_R, DATA.step+" P aligned to R mech "+mech.id); else downloadXYZ(safeName(DATA.step)+"_P.xyz", DATA.product.elements, DATA.product.coords, DATA.step+" P"); }}
function downloadGT() {{ const mech = findMech(currentMechId); if (mech.gt) downloadXYZ(safeName(DATA.step)+"_GT_mech"+mech.id+".xyz", mech.gt.elements || elements, mech.gt.xyz || mech.gt.xyz_in_R, DATA.step+" GT mech "+mech.id); }}
function downloadIG(ig) {{ downloadXYZ(safeName(DATA.step)+"_"+safeName(ig.label)+".xyz", ig.elements || elements, ig.xyz || ig.xyz_in_R, DATA.step+" "+ig.label); }}
function scoreRecord(item) {{ if (!item || item.S === undefined || item.S === null) return null; return {{S:item.S, decomposition:{{beta:item.beta, wbo_progress:item.wbo_progress, wbo_progress_factor:item.wbo_progress_factor, freq:item.freq, mode_index:item.k}}, event_terms:item.event_terms || [], core_map:item.core_map, core_sources:item.core_sources, core_pool_dedup_count:item.core_pool_dedup_count}}; }}
function mechanismRecord(mech) {{ return {{id:mech.id, label:mech.label, cut:mech.cut, dedup_count:mech.dedup_count || 1, dedup_source_ids:mech.dedup_source_ids || [mech.id], dedup_cuts:mech.dedup_cuts || [mech.cut], broken_bonds_R:mech.broken_bonds_R, formed_bonds_R:mech.formed_bonds_R, formed_bonds_P:mech.formed_bonds_P || [], core_atoms_R:mech.core_atoms || [], gt:scoreRecord(mech.gt), igs:(mech.igs || []).map(ig => ({{label:ig.label, is_top2:!!ig.is_top2, is_union_top:!!ig.is_union_top, score:scoreRecord(ig)}}))}}; }}
function buildArchiveManifest() {{ return {{step:DATA.step, n_atoms:DATA.n_atoms, include_gt:!!DATA.include_gt, default_mech_id:DATA.default_mech_id, score_formula:"S = beta * wbo_progress^WBO_PROGRESS_POWER", score_config:DATA.score_config || null, mechanisms:(DATA.mechanisms || []).map(mechanismRecord), files:{{reactant:"R.xyz", product:"P.xyz", gt:"GT/GT.xyz if available", ig:"IG/<label>.xyz", per_mechanism:"mechanisms/mechanism_<id>.json", full_viewer_data:"viewer_data.json"}}}}; }}
function scoreMeta(item) {{ if (!item || item.beta === undefined) return "(no data)"; return "<b>&beta;</b>="+item.beta.toFixed(3)+" <b>wbo</b>="+item.wbo_progress.toFixed(3)+" <b>freq</b>="+item.freq.toFixed(0)+"i"; }}
async function downloadAll() {{ if (typeof JSZip === "undefined") {{ alert("Download library is not loaded"); return; }} const root = safeName(DATA.step); const zip = new JSZip(); zip.file(root+"/R.xyz", xyzText(DATA.reactant.elements, DATA.reactant.coords, DATA.step+" R")); zip.file(root+"/P.xyz", xyzText(DATA.product.elements, DATA.product.coords, DATA.step+" P")); const firstGT = (DATA.mechanisms || []).map(m => m.gt).find(gt => gt && (gt.xyz || gt.xyz_in_R)); if (firstGT) zip.file(root+"/GT/GT.xyz", xyzText(firstGT.elements || elements, firstGT.xyz || firstGT.xyz_in_R, DATA.step+" GT")); const seenIG = new Set(); for (const mech of DATA.mechanisms || []) {{ for (const ig of mech.igs || []) {{ if (seenIG.has(ig.label)) continue; const xyz = ig.xyz || ig.xyz_in_R; if (!xyz) continue; seenIG.add(ig.label); zip.file(root+"/IG/"+safeName(ig.label)+".xyz", xyzText(ig.elements || elements, xyz, DATA.step+" "+ig.label)); }} }} const manifest = buildArchiveManifest(); zip.file(root+"/mechanism.json", JSON.stringify(manifest, null, 2)); for (const mech of DATA.mechanisms || []) {{ zip.file(root+"/mechanisms/mechanism_"+String(mech.id).padStart(3,"0")+".json", JSON.stringify(mechanismRecord(mech), null, 2)); }} zip.file(root+"/viewer_data.json", JSON.stringify(DATA, null, 2)); const blob = await zip.generateAsync({{type:"blob"}}); downloadBlob(root+".zip", blob); }}
const animTimers = {{}};
function stopAnim(d) {{ if (animTimers[d]) {{ clearInterval(animTimers[d]); delete animTimers[d]; }} }}
function clearLabels(v) {{ if (v.removeAllLabels) v.removeAllLabels(); }}
function addAtomLabels(v, els, xyz) {{ if (!showAtomIndices || !els || !xyz) return; for (let i=0;i<xyz.length;i++) {{ v.addLabel(String(i), {{position:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, fontSize:10, fontColor:'black', backgroundColor:'white', backgroundOpacity:0.72, borderColor:'#666', borderThickness:0.5, inFront:true}}); }} }}
function drawBonds(v, xyz, pairs, color) {{ for (const [i,j] of pairs) {{ if (i>=xyz.length||j>=xyz.length) continue; v.addCylinder({{start:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, end:{{x:xyz[j][0],y:xyz[j][1],z:xyz[j][2]}}, color:color, radius:0.16, dashed:true}}); }} }}
function drawArrows(v, xyz, disp, core) {{ for (const i of core) {{ if (!disp||!disp[i]) continue; const d = disp[i]; const len = Math.hypot(d[0],d[1],d[2]); if (len<0.05) continue; v.addArrow({{start:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, end:{{x:xyz[i][0]+d[0]*1.5,y:xyz[i][1]+d[1]*1.5,z:xyz[i][2]+d[2]*1.5}}, color:'#0066cc', radius:0.07}}); }} }}
function makeStatic(divId, els, xyz, broken, formed) {{ stopAnim(divId); document.getElementById(divId).innerHTML=""; const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}}); v.addModel(buildBody(els, xyz), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); drawBonds(v, xyz, broken, 'red'); drawBonds(v, xyz, formed, 'green'); addAtomLabels(v, els, xyz); v.zoomTo(); v.render(); return v; }}
function makeAnimated(divId, els, xyz, disp, broken, formed, core) {{ stopAnim(divId); document.getElementById(divId).innerHTML=""; const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}}); v.addModel(buildBody(els, xyz), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); drawBonds(v, xyz, broken, 'red'); drawBonds(v, xyz, formed, 'green'); drawArrows(v, xyz, disp, core); addAtomLabels(v, els, xyz); v.zoomTo(); v.render(); let t=0; const period=30, amp=0.6; animTimers[divId] = setInterval(()=>{{ t=(t+1)%period; const scale = amp*Math.sin(2*Math.PI*t/period); const cur = xyzAt(xyz, disp, scale); v.removeAllModels(); v.removeAllShapes(); clearLabels(v); v.addModel(buildBodyAt(els, xyz, disp, scale), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); drawBonds(v, cur, broken, 'red'); drawBonds(v, cur, formed, 'green'); drawArrows(v, cur, disp, core); addAtomLabels(v, els, cur); v.render(); }}, 60); return v; }}
function render() {{ const mech = findMech(currentMechId); document.querySelectorAll('.mech-sel button[data-id]').forEach(b => {{ b.classList.toggle('active', parseInt(b.dataset.id)===currentMechId); }}); makeStatic('vw_R', DATA.reactant.elements, DATA.reactant.coords, mech.broken_bonds_R, []); const pAligned = !!mech.product_xyz_in_R; const pEls = pAligned ? DATA.reactant.elements : DATA.product.elements; const pXYZ = pAligned ? mech.product_xyz_in_R : DATA.product.coords; const pFormed = pAligned ? mech.formed_bonds_R : (mech.formed_bonds_P || []); makeStatic('vw_P', pEls, pXYZ, [], pFormed); document.getElementById('prod_label').innerHTML = (pAligned ? "aligned to R" : "static")+" (mech #"+mech.id+") <button class='dl' onclick='downloadP()'>XYZ</button>"; const showGT = !!(mech.gt && mech.gt.picked_disp); document.getElementById('ref-row').classList.toggle('no-gt', !showGT); document.getElementById('gt_panel').style.display = showGT ? "" : "none"; if (showGT) {{ makeAnimated('vw_GT', mech.gt.elements || elements, mech.gt.xyz || mech.gt.xyz_in_R, mech.gt.picked_disp, mech.gt.broken_bonds_T || mech.broken_bonds_R, mech.gt.formed_bonds_T || mech.formed_bonds_R, mech.gt.core_atoms_T || mech.core_atoms); document.getElementById('gt_S').innerHTML = "S = "+mech.gt.S.toFixed(3)+" <button class='dl' onclick='downloadGT()'>XYZ</button>"; document.getElementById('gt_meta').innerHTML = scoreMeta(mech.gt); }} else {{ stopAnim('vw_GT'); document.getElementById('vw_GT').innerHTML = ""; }} const grid = document.getElementById('grid'); grid.innerHTML = ""; const igs = [...mech.igs].sort((a,b) => (b.S||0) - (a.S||0)); igs.forEach((ig, idx) => {{ const div = document.createElement('div'); let cls = 'panel'; if (ig.is_top2) cls += ' top2'; if (ig.is_union_top && !ig.is_top2) cls += ' union'; div.className = cls; const sStr = ig.S !== undefined ? "S = "+ig.S.toFixed(3) : "no score"; const tag = ig.is_top2 ? '<span style="background:#d4af37;color:white;padding:1px 6px;border-radius:3px;font-size:11px;margin-left:4px">TOP2</span>' : (ig.is_union_top ? '<span style="background:#ff9;color:#660;padding:1px 6px;border-radius:3px;font-size:11px;margin-left:4px">union</span>' : ''); const dl = (ig.xyz || ig.xyz_in_R) ? '<button class="dl">XYZ</button> ' : ''; div.innerHTML = '<div class="ph"><span class="lbl">'+ig.label+tag+'</span><span class="rk">'+dl+sStr+'</span></div><div class="vw"><div id="vw_ig'+idx+'" class="vwbox"></div></div><div class="meta">'+scoreMeta(ig)+"</div>"; grid.appendChild(div); const btn = div.querySelector('button.dl'); if (btn) btn.onclick = () => downloadIG(ig); if (ig.picked_disp) makeAnimated("vw_ig"+idx, ig.elements || elements, ig.xyz || ig.xyz_in_R, ig.picked_disp, ig.broken_bonds_T || mech.broken_bonds_R, ig.formed_bonds_T || mech.formed_bonds_R, ig.core_atoms_T || mech.core_atoms); else if (ig.xyz || ig.xyz_in_R) makeStatic("vw_ig"+idx, ig.elements || elements, ig.xyz || ig.xyz_in_R, ig.broken_bonds_T || mech.broken_bonds_R, ig.formed_bonds_T || mech.formed_bonds_R); }}); }}
const ms = document.getElementById('mech-sel'); ms.innerHTML = "<span style='font-size:13px;margin-right:8px;color:#444'>Mechanism:</span>"; document.getElementById('downloadAllBtn').onclick = downloadAll; document.getElementById('showAtomIndices').onchange = (e) => {{ showAtomIndices = !!e.target.checked; render(); }}; DATA.mechanisms.forEach(m => {{ const b = document.createElement('button'); b.dataset.id = m.id; b.textContent = m.label + (m.gt ? "  GT S=" + m.gt.S.toFixed(3) : ""); if ((m.dedup_count||1) > 1) b.title = "Collapsed raw witnesses: "+m.dedup_count+"; source mechanisms: "+m.dedup_source_ids.join(", ")+"; cuts: "+m.dedup_cuts.join(", "); b.onclick = () => {{ currentMechId = m.id; render(); }}; ms.appendChild(b); }});
window.addEventListener('load', render);
</script>
</body></html>
"""


def run_full_pipeline_stage(step_name, inner_workers=0, mechanism_ids=None,
                            write_artifacts=True,
                            save_alignment_files=False, resume_rp=False):
    """Compose Stage 1, Stage 2, and Stage 3 for one step.

    If `resume_rp=True`, the R-P mechanism-discovery artifact is loaded from
    `rp_stage.json` and only the collective TS/IG validation plus view/export
    stages are run.  Missing Stage 1 artifacts are an error in resume mode.
    """
    inputs = load_step_inputs(step_name)
    paths = pipeline_stage_paths(step_name)
    if resume_rp:
        if not paths.rp_json.exists():
            raise RuntimeError(
                f"missing Stage 1 artifact: {paths.rp_json}; run --stage rp first")
        rp_result = read_stage_json(paths.rp_json)
    else:
        rp_result = run_rp_stage(inputs, inner_workers=inner_workers)
    targets = load_ts_targets(inputs, include_gt=INCLUDE_GT)
    ts_result = run_ts_stage(
        inputs, rp_result, targets, inner_workers=inner_workers,
        mechanism_ids=mechanism_ids)
    view_result = write_view_stage(
        inputs, rp_result, ts_result, include_gt=INCLUDE_GT)
    if write_artifacts:
        if not resume_rp:
            write_stage_json(paths.rp_json, rp_result)
        write_stage_json(paths.ts_json, ts_result)
    alignment_files = None
    ts_alignment_files = None
    if save_alignment_files:
        alignment_files = write_rp_alignment_files(inputs, rp_result)
        ts_alignment_files = write_ts_alignment_files(inputs, ts_result)
    return {
        'step': step_name,
        'rp': rp_result,
        'ts': ts_result,
        'view': view_result,
        'alignment_files': alignment_files,
        'ts_alignment_files': ts_alignment_files,
        'slim': view_result['slim'],
    }


def process_rp_step(step_name, inner_workers=0, write_view=True,
                    save_alignment_files=False):
    """Run only Stage 1 and write `rp_stage.json`.

    The optional view is an alignment/mechanism view: R/P only, no TS scores.
    """
    inputs = load_step_inputs(step_name)
    rp_result = run_rp_stage(inputs, inner_workers=inner_workers)
    paths = pipeline_stage_paths(step_name)
    write_stage_json(paths.rp_json, rp_result)
    view_result = None
    if write_view:
        view_result = write_view_stage(
            inputs, rp_result, ts_result=None, include_gt=False)
    alignment_files = None
    if save_alignment_files:
        alignment_files = write_rp_alignment_files(inputs, rp_result)
    return {
        'step': step_name,
        'rp': rp_result,
        'view': view_result,
        'alignment_files': alignment_files,
        'slim': view_result['slim'] if view_result else {
            'step': step_name,
            'n_mechs': len(rp_result.get('mechanisms', [])),
            'mechanisms': rp_result.get('mechanisms', []),
        },
    }


def process_ts_step(step_name, inner_workers=0, mechanism_ids=None,
                    write_view=True, save_alignment_files=False):
    """Run Stage 2 from a previously written `rp_stage.json`."""
    inputs = load_step_inputs(step_name)
    paths = pipeline_stage_paths(step_name)
    if not paths.rp_json.exists():
        raise RuntimeError(
            f"missing Stage 1 artifact: {paths.rp_json}; run --stage rp first")
    rp_result = read_stage_json(paths.rp_json)
    targets = load_ts_targets(inputs, include_gt=INCLUDE_GT)
    ts_result = run_ts_stage(
        inputs, rp_result, targets, inner_workers=inner_workers,
        mechanism_ids=mechanism_ids)
    write_stage_json(paths.ts_json, ts_result)
    view_result = None
    if write_view:
        view_result = write_view_stage(
            inputs, rp_result, ts_result, include_gt=INCLUDE_GT)
    ts_alignment_files = None
    if save_alignment_files:
        ts_alignment_files = write_ts_alignment_files(inputs, ts_result)
    return {
        'step': step_name,
        'rp': rp_result,
        'ts': ts_result,
        'view': view_result,
        'ts_alignment_files': ts_alignment_files,
        'slim': view_result['slim'] if view_result else build_eval_slim(
            build_view_data(inputs, rp_result, ts_result,
                            include_gt=INCLUDE_GT)),
    }


def process_view_step(step_name):
    """Run only Stage 3 from previously written stage artifacts."""
    inputs = load_step_inputs(step_name)
    paths = pipeline_stage_paths(step_name)
    if not paths.rp_json.exists():
        raise RuntimeError(
            f"missing Stage 1 artifact: {paths.rp_json}; run --stage rp first")
    rp_result = read_stage_json(paths.rp_json)
    ts_result = read_stage_json(paths.ts_json) if paths.ts_json.exists() else None
    view_result = write_view_stage(
        inputs, rp_result, ts_result, include_gt=INCLUDE_GT)
    return {
        'step': step_name,
        'rp': rp_result,
        'ts': ts_result,
        'view': view_result,
        'slim': view_result['slim'],
    }


def _target_specs_from_cli(target_xyzs, target_labels=None, target_kinds=None):
    specs = []
    target_labels = list(target_labels or [])
    target_kinds = list(target_kinds or [])
    for i, xyz in enumerate(target_xyzs or []):
        path = Path(xyz)
        specs.append({
            'xyz': path,
            'label': (
                target_labels[i] if i < len(target_labels)
                else path.stem
            ),
            'kind': (
                target_kinds[i] if i < len(target_kinds)
                else 'ig'
            ),
        })
    return specs


def _load_ts_target_from_spec_task(task):
    i, spec, workdir, charge, multiplicity, xtb_mode = task
    label = spec.get('label') or Path(spec['xyz']).stem
    kind = str(spec.get('kind', 'ig')).lower()
    target = ts_target_from_xyz(
        kind, label, spec['xyz'],
        workdir=Path(workdir),
        target_index=-(i + 1) if kind == 'gt' else i,
        charge=spec.get('charge', charge),
        multiplicity=spec.get('multiplicity', multiplicity),
        xtb_mode=spec.get('xtb_mode', xtb_mode),
    )
    return i, target


def load_ts_targets_from_specs(target_specs, workdir, *, charge=None,
                               multiplicity=None, xtb_mode=None,
                               inner_workers=0):
    """Load arbitrary Stage 2 targets from XYZ specs.

    Each spec is a dict with `xyz`, optional `label`, and optional `kind`
    (`ig` or `gt`).  Cache directories live under `workdir` unless the caller
    builds targets directly with `ts_target_from_xyz(..., sp_workdir=..., ...)`.
    """
    specs = list(target_specs or [])
    if not specs:
        return []

    tasks = [
        (i, spec, Path(workdir), charge, multiplicity, xtb_mode)
        for i, spec in enumerate(specs)
    ]
    workers = max(1, int(inner_workers or 1))
    if workers <= 1 or len(tasks) == 1:
        return [
            target
            for _i, target in (
                _load_ts_target_from_spec_task(task) for task in tasks
            )
        ]

    max_workers = min(workers, len(tasks))
    with cf.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(_load_ts_target_from_spec_task, task)
            for task in tasks
        ]
        results = [future.result() for future in cf.as_completed(futures)]
    return [target for _i, target in sorted(results, key=lambda item: item[0])]


def process_xyz_stage(name, reactant_xyz, product_xyz, *, workdir=None,
                      target_specs=None, stage='full', inner_workers=0,
                      mechanism_ids=None, save_alignment_files=False,
                      charge=None, multiplicity=None, xtb_mode=None,
                      rp_config=None, ts_config=None, resume_rp=False):
    """Run pipeline stages for arbitrary XYZ files, without a step schema."""
    workdir = Path(workdir or (PROJECT / "out" / "xyz_work" / name))
    inputs = alignment_inputs_from_xyz(
        reactant_xyz, product_xyz, workdir=workdir / "endpoints",
        name=name, charge=charge, multiplicity=multiplicity,
        xtb_mode=xtb_mode)
    paths = pipeline_stage_paths(name)
    include_gt = INCLUDE_GT or any(
        str(s.get('kind', '')).lower() == 'gt' for s in (target_specs or []))

    if stage == 'rp':
        rp_result = run_rp_stage(
            inputs, config=rp_config, inner_workers=inner_workers)
        write_stage_json(paths.rp_json, rp_result)
        view_result = write_view_stage(
            inputs, rp_result, ts_result=None, include_gt=False)
        alignment_files = (
            write_rp_alignment_files(inputs, rp_result)
            if save_alignment_files else None
        )
        return {
            'step': name,
            'rp': rp_result,
            'view': view_result,
            'alignment_files': alignment_files,
            'slim': view_result['slim'],
        }

    if stage == 'view':
        if not paths.rp_json.exists():
            raise RuntimeError(
                f"missing Stage 1 artifact: {paths.rp_json}; run --stage rp first")
        rp_result = read_stage_json(paths.rp_json)
        ts_result = read_stage_json(paths.ts_json) if paths.ts_json.exists() else None
        view_result = write_view_stage(
            inputs, rp_result, ts_result, include_gt=include_gt)
        return {
            'step': name,
            'rp': rp_result,
            'ts': ts_result,
            'view': view_result,
            'slim': view_result['slim'],
        }

    if stage in {'ts', 'post-rp'} or (stage == 'full' and resume_rp):
        if not paths.rp_json.exists():
            raise RuntimeError(
                f"missing Stage 1 artifact: {paths.rp_json}; run --stage rp first")
        rp_result = read_stage_json(paths.rp_json)
    elif stage == 'full':
        rp_result = run_rp_stage(
            inputs, config=rp_config, inner_workers=inner_workers)
        write_stage_json(paths.rp_json, rp_result)
    else:
        raise ValueError(f"unknown stage: {stage}")

    targets = load_ts_targets_from_specs(
        target_specs or [], workdir / "targets",
        charge=charge, multiplicity=multiplicity, xtb_mode=xtb_mode,
        inner_workers=inner_workers)
    ts_result = run_ts_stage(
        inputs, rp_result, targets, config=ts_config,
        inner_workers=inner_workers, mechanism_ids=mechanism_ids)
    write_stage_json(paths.ts_json, ts_result)
    view_result = write_view_stage(
        inputs, rp_result, ts_result, include_gt=include_gt)
    alignment_files = ts_alignment_files = None
    if save_alignment_files:
        if stage == 'full' and not resume_rp:
            alignment_files = write_rp_alignment_files(inputs, rp_result)
        ts_alignment_files = write_ts_alignment_files(inputs, ts_result)
    return {
        'step': name,
        'rp': rp_result,
        'ts': ts_result,
        'view': view_result,
        'alignment_files': alignment_files,
        'ts_alignment_files': ts_alignment_files,
        'slim': view_result['slim'],
    }


def process_step_stage(step_name, stage='full', inner_workers=0,
                       mechanism_ids=None, save_alignment_files=False,
                       resume_rp=False):
    """CLI-safe wrapper for one explicit stage."""
    try:
        if stage == 'rp':
            return process_rp_step(
                step_name, inner_workers=inner_workers,
                save_alignment_files=save_alignment_files)
        if stage in {'ts', 'post-rp'}:
            return process_ts_step(
                step_name, inner_workers=inner_workers,
                mechanism_ids=mechanism_ids,
                save_alignment_files=save_alignment_files)
        if stage == 'view':
            return process_view_step(step_name)
        if stage == 'full':
            return run_full_pipeline_stage(
                step_name, inner_workers=inner_workers,
                mechanism_ids=mechanism_ids,
                save_alignment_files=save_alignment_files,
                resume_rp=resume_rp)
        raise ValueError(f"unknown stage: {stage}")
    except Exception as e:
        return {"step": step_name, "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()}


def _process_step_stage_star(args):
    return process_step_stage(*args)


def process_step(step_name, inner_workers=0):
    """Compatibility wrapper for the historical full pipeline."""
    return process_step_stage(
        step_name, stage='full', inner_workers=inner_workers)


def main():
    global INCLUDE_GT, XTB_CACHE_MODE, XTB_OMP_THREADS, XTB_MAX_THREADS
    global XTB_CHARGE, XTB_MULTIPLICITY
    global VIEW_ISO_TOL, DWBO_THRESHOLD, METAL_DWBO_THRESHOLD, SYMMETRY_WBO_TOL
    global EVENT_WEIGHT_POWER, WBO_PROGRESS_POWER
    global STAGE_ROOT, ALIGNMENT_OUT_ROOT
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage",
                    choices=("full", "rp", "ts", "post-rp", "view"),
                    default=os.environ.get("BGCP_STAGE", "full"),
                    help="Pipeline stage to run. full composes rp+ts+view; "
                         "rp writes rp_stage.json; ts/post-rp resume from "
                         "rp_stage.json and run collective TS/IG validation; "
                         "view only rewrites HTML/eval.")
    ap.add_argument("--stage-root", default=str(STAGE_ROOT),
                    help="Directory for resumable rp_stage.json and "
                         "ts_stage.json artifacts. Default from "
                         "BGCP_STAGE_ROOT=out/bgcp_stages.")
    ap.add_argument("--alignment-out-root", default=str(ALIGNMENT_OUT_ROOT),
                    help="Directory for optional clean aligned "
                         "coordinate exports. Default from "
                         "BGCP_ALIGNMENT_OUT_ROOT=out/bgcp_alignments.")
    ap.add_argument("--save-alignment-files", action="store_true",
                    default=os.environ.get(
                        "BGCP_SAVE_ALIGNMENT_FILES", "0").lower() in {
                            "1", "true", "yes", "on"
                        },
                    help="Write clean mechanism-specific R/P-aligned XYZ "
                         "files during --stage rp or --stage full, and "
                         "best-S TS core-aligned files during --stage "
                         "ts/post-rp or --stage full.")
    ap.add_argument("--resume-rp", action="store_true",
                    default=os.environ.get(
                        "BGCP_RESUME_RP", "0").lower() in {
                            "1", "true", "yes", "on"
                        },
                    help="For --stage full, reuse existing rp_stage.json and "
                         "run only post-Stage-1 collective TS/IG validation "
                         "plus view/export. Missing rp_stage.json is an error.")
    ap.add_argument("--mechanism", type=int, action="append",
                    help="Restrict Stage 2 verification to a mechanism id. "
                         "Can be repeated. Default verifies all mechanisms.")
    ap.add_argument("--workers", type=int, default=_default_worker_count(),
                    help="Total CPU budget in auto mode, or outer step "
                         "parallelism in outer mode.")
    ap.add_argument("--inner-workers", type=int, default=0,
                    help="Explicit workers per step's inner R-P/TS work. In auto "
                         "mode, 0 means choose from --workers; >1 switches "
                         "to inner mode unless --parallel-mode is set.")
    ap.add_argument("--parallel-mode",
                    choices=("auto", "outer", "inner"),
                    default=os.environ.get("BGCP_PARALLEL_MODE", "auto"),
                    help="auto balances outer steps and inner workers; "
                         "outer is legacy many-steps/serial-inside-step mode; "
                         "inner runs steps serially with parallel inner work.")
    ap.add_argument("--auto-inner-workers", type=int,
                    default=AUTO_INNER_WORKERS,
                    help="Target inner workers per concurrent step in "
                         "auto mode. Default from BGCP_AUTO_INNER_WORKERS=8.")
    ap.add_argument("--iso-tol", type=float, default=VIEW_ISO_TOL,
                    help="WBO tolerance for active graph matching. "
                         "Default from BGCP_ISO_TOL=1.0.")
    ap.add_argument("--dwbo-threshold", type=float, default=DWBO_THRESHOLD,
                    help="Delta-WBO threshold for broken/formed bond "
                         "classification for non-metal pairs. Default from "
                         "BGCP_DWBO_THRESHOLD=0.5.")
    ap.add_argument("--metal-dwbo-threshold", type=float,
                    default=METAL_DWBO_THRESHOLD,
                    help="Delta-WBO threshold for broken/formed events where "
                         "either endpoint is a metal. Default from "
                         "BGCP_METAL_DWBO_THRESHOLD=0.3.")
    ap.add_argument("--symmetry-wbo-tol", type=float,
                    default=SYMMETRY_WBO_TOL,
                    help="WBO tolerance for symmetry-orbit bucketing. "
                         "Default from BGCP_SYMMETRY_WBO_TOL=0.2.")
    ap.add_argument("--event-weight-power", type=float,
                    default=EVENT_WEIGHT_POWER,
                    help="Exponent on each detected R-P event's |delta WBO| "
                         "when building the beta vector. Default from "
                         "BGCP_EVENT_WEIGHT_POWER=1.0.")
    ap.add_argument("--wbo-progress-power", type=float,
                    default=WBO_PROGRESS_POWER,
                    help="Exponent on the TS WBO progress factor in the "
                         "ranker. Default from BGCP_WBO_PROGRESS_POWER=1.0.")
    ap.add_argument("--xtb-mode",
                    choices=("auto", "cache-only"),
                    default=_normal_xtb_mode(XTB_CACHE_MODE),
                    help="auto fills missing WBO/g98 caches by running xtb "
                         "from available XYZ files; cache-only fails on "
                         "missing cache files.")
    ap.add_argument("--xtb-omp-threads",
                    default=XTB_OMP_THREADS,
                    help="Requested OMP_NUM_THREADS per xtb molecule. "
                         "Default 'auto' uses available CPUs capped by "
                         "--xtb-max-threads.")
    ap.add_argument("--xtb-max-threads", type=int,
                    default=XTB_MAX_THREADS,
                    help="Hard cap for OMP_NUM_THREADS per xtb molecule. "
                         "Default from BGCP_XTB_MAX_THREADS=8.")
    ap.add_argument("--charge", type=int, default=XTB_CHARGE,
                    help="Total molecular charge used only when auto-filling "
                         "missing xtb caches. Default from BGCP_CHARGE=0.")
    ap.add_argument("--multiplicity", type=int, default=XTB_MULTIPLICITY,
                    help="Spin multiplicity used only when auto-filling "
                         "missing xtb caches. Converted internally to "
                         "xtb --uhf=multiplicity-1. Default from "
                         "BGCP_MULTIPLICITY=1.")
    ap.add_argument("--include-gt", action="store_true",
                    default=INCLUDE_GT,
                    help="Load and score sp_groundtruth/hess_groundtruth. "
                         "Default is off for benchmark runs without GT.")
    ap.add_argument("--name", default=None,
                    help="Name for direct XYZ mode. Defaults to "
                         "<reactant-stem>_to_<product-stem>.")
    ap.add_argument("--reactant-xyz", default=None,
                    help="Direct Stage 1 reactant endpoint XYZ. Use with "
                         "--product-xyz instead of --steps.")
    ap.add_argument("--product-xyz", default=None,
                    help="Direct Stage 1 product endpoint XYZ. Use with "
                         "--reactant-xyz instead of --steps.")
    ap.add_argument("--workdir", default=None,
                    help="Direct XYZ mode cache work directory. Holds "
                         "endpoint, TS single-point, and TS Hessian caches.")
    ap.add_argument("--target-xyz", action="append", default=None,
                    help="Direct Stage 2 TS/IG/GT XYZ. Can be repeated.")
    ap.add_argument("--target-label", action="append", default=None,
                    help="Label for each --target-xyz. Defaults to file stem.")
    ap.add_argument("--target-kind", action="append",
                    choices=("ig", "gt"), default=None,
                    help="Kind for each --target-xyz. Defaults to ig.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--steps", nargs="+", default=None)
    args = ap.parse_args()

    XTB_CACHE_MODE = _normal_xtb_mode(args.xtb_mode)
    INCLUDE_GT = bool(args.include_gt)
    VIEW_ISO_TOL = float(args.iso_tol)
    DWBO_THRESHOLD = float(args.dwbo_threshold)
    METAL_DWBO_THRESHOLD = float(args.metal_dwbo_threshold)
    SYMMETRY_WBO_TOL = float(args.symmetry_wbo_tol)
    EVENT_WEIGHT_POWER = float(args.event_weight_power)
    WBO_PROGRESS_POWER = float(args.wbo_progress_power)
    STAGE_ROOT = Path(args.stage_root)
    ALIGNMENT_OUT_ROOT = Path(args.alignment_out_root)
    XTB_MAX_THREADS = max(1, int(args.xtb_max_threads))
    XTB_CHARGE = int(args.charge)
    XTB_MULTIPLICITY = _normal_multiplicity(args.multiplicity)
    XTB_OMP_THREADS = _resolve_xtb_threads(
        args.xtb_omp_threads, XTB_MAX_THREADS)
    os.environ["BGCP_XTB_MODE"] = XTB_CACHE_MODE
    os.environ["BGCP_INCLUDE_GT"] = "1" if INCLUDE_GT else "0"
    os.environ["BGCP_XTB_OMP_THREADS"] = str(XTB_OMP_THREADS)
    os.environ["BGCP_XTB_MAX_THREADS"] = str(XTB_MAX_THREADS)
    os.environ["BGCP_CHARGE"] = str(XTB_CHARGE)
    os.environ["BGCP_MULTIPLICITY"] = str(XTB_MULTIPLICITY)
    os.environ["BGCP_ISO_TOL"] = str(VIEW_ISO_TOL)
    os.environ["BGCP_DWBO_THRESHOLD"] = str(DWBO_THRESHOLD)
    os.environ["BGCP_METAL_DWBO_THRESHOLD"] = str(METAL_DWBO_THRESHOLD)
    os.environ["BGCP_SYMMETRY_WBO_TOL"] = str(SYMMETRY_WBO_TOL)
    os.environ["BGCP_EVENT_WEIGHT_POWER"] = str(EVENT_WEIGHT_POWER)
    os.environ["BGCP_WBO_PROGRESS_POWER"] = str(WBO_PROGRESS_POWER)
    os.environ["BGCP_STAGE_ROOT"] = str(STAGE_ROOT)
    os.environ["BGCP_ALIGNMENT_OUT_ROOT"] = str(ALIGNMENT_OUT_ROOT)
    os.environ["BGCP_SAVE_ALIGNMENT_FILES"] = (
        "1" if args.save_alignment_files else "0")
    os.environ["BGCP_RESUME_RP"] = "1" if args.resume_rp else "0"

    direct_mode = bool(args.reactant_xyz or args.product_xyz)
    if direct_mode:
        if args.steps:
            ap.error("use either --steps or --reactant-xyz/--product-xyz, not both")
        if not (args.reactant_xyz and args.product_xyz):
            ap.error("--reactant-xyz and --product-xyz must be provided together")
        n_targets = len(args.target_xyz or [])
        if args.target_label and len(args.target_label) > n_targets:
            ap.error("--target-label cannot be provided more times than --target-xyz")
        if args.target_kind and len(args.target_kind) > n_targets:
            ap.error("--target-kind cannot be provided more times than --target-xyz")
        name = args.name or (
            f"{Path(args.reactant_xyz).stem}_to_{Path(args.product_xyz).stem}")
        OUT_ROOT.mkdir(parents=True, exist_ok=True)
        STAGE_ROOT.mkdir(parents=True, exist_ok=True)
        if args.save_alignment_files:
            ALIGNMENT_OUT_ROOT.mkdir(parents=True, exist_ok=True)
        try:
            rec = process_xyz_stage(
                name,
                Path(args.reactant_xyz),
                Path(args.product_xyz),
                workdir=Path(args.workdir) if args.workdir else None,
                target_specs=_target_specs_from_cli(
                    args.target_xyz, args.target_label, args.target_kind),
                stage=args.stage,
                inner_workers=(
                    args.inner_workers if args.inner_workers > 0
                    else max(1, int(args.workers))),
                mechanism_ids=args.mechanism,
                save_alignment_files=args.save_alignment_files,
                charge=XTB_CHARGE,
                multiplicity=XTB_MULTIPLICITY,
                xtb_mode=XTB_CACHE_MODE,
                resume_rp=args.resume_rp,
            )
        except Exception as e:
            rec = {
                'step': name,
                'error': f"{type(e).__name__}: {e}",
                'trace': traceback.format_exc(),
            }
        if rec.get('error'):
            print(f"{name}: ERROR: {rec['error']}")
            EVAL_JSON.write_text(json.dumps([{
                'step': name,
                'error': rec['error'],
            }]))
        else:
            slim = rec['slim']
            print(f"{name}: mechs={slim.get('n_mechs', 0)} "
                  f"view={rec.get('view', {}).get('view_html')}")
            EVAL_JSON.write_text(json.dumps([slim]))
        print(f"wrote {EVAL_JSON}")
        return

    if not WORK.exists():
        print(f"No work directory: {WORK}")
        return
    all_steps = sorted(d.name for d in WORK.iterdir() if d.is_dir())
    if args.steps: steps = [s for s in all_steps if s in set(args.steps)]
    elif args.limit: steps = all_steps[:args.limit]
    else: steps = all_steps

    if not steps:
        print("No matching steps.")
        return

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    STAGE_ROOT.mkdir(parents=True, exist_ok=True)
    if args.save_alignment_files:
        ALIGNMENT_OUT_ROOT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    eval_records = []
    n_ok = n_err = 0

    def _record(i, rec):
        nonlocal n_ok, n_err
        if rec.get("error"):
            print(f"  [{i:>3d}/{len(steps)}] {rec['step']:60s}  ERROR: {rec['error'][:80]}", flush=True)
            n_err += 1
            eval_records.append({'step': rec['step'], 'error': rec['error']})
        else:
            slim = rec['slim']
            gt_scores = [
                m['gt']['S'] for m in slim.get('mechanisms', [])
                if m.get('gt')
            ]
            if gt_scores:
                gt_msg = f"best_GT_S={max(gt_scores):.3f}"
            else:
                gt_msg = "GT=skipped"
            print(f"  [{i:>3d}/{len(steps)}] {rec['step']:60s}  "
                  f"mechs={slim.get('n_mechs', len(slim.get('mechanisms', [])))}  "
                  f"{gt_msg}", flush=True)
            eval_records.append(slim)
            n_ok += 1

    mode = args.parallel_mode
    if mode == "auto" and args.inner_workers and args.inner_workers > 1:
        mode = "inner"
    worker_budget = max(1, int(args.workers))

    if mode == "inner":
        # Inner-parallel mode: steps run serially in main; each step's
        # cut_sweep and TS endpoint matching use inner_workers cores. Best for
        # a single step or a few large steps where inner work dominates cost.
        inner_workers = (args.inner_workers if args.inner_workers > 0
                         else worker_budget)
        print(f"Processing {len(steps)} steps serially; each step uses "
              f"{inner_workers} inner workers "
              f"(stage={args.stage}, "
              f"cut_sweep chunksize={CUTSWEEP_CHUNKSIZE}, "
              f"iso_tol={VIEW_ISO_TOL}, dwbo={DWBO_THRESHOLD}, "
              f"metal_dwbo={METAL_DWBO_THRESHOLD}, "
              f"sym_wbo_tol={SYMMETRY_WBO_TOL}, "
              f"event_weight_power={EVENT_WEIGHT_POWER}, "
              f"wbo_progress_power={WBO_PROGRESS_POWER}, "
              f"xtb_mode={XTB_CACHE_MODE}, "
              f"xtb_threads={XTB_OMP_THREADS}, charge={XTB_CHARGE}, "
              f"multiplicity={XTB_MULTIPLICITY}, include_gt={INCLUDE_GT}, "
              f"save_alignment_files={args.save_alignment_files}, "
              f"resume_rp={args.resume_rp})")
        for i, step in enumerate(steps, 1):
            rec = process_step_stage(
                step, stage=args.stage, inner_workers=inner_workers,
                mechanism_ids=args.mechanism,
                save_alignment_files=args.save_alignment_files,
                resume_rp=args.resume_rp)
            _record(i, rec)
    elif mode == "outer":
        # Outer-parallel mode: args.workers steps run concurrently; each
        # step's inner work is serial (no nested daemonic Pool). Best for
        # small/easy steps when nested process pools are undesirable.
        print(f"Processing {len(steps)} steps with {worker_budget} outer workers "
              f"(stage={args.stage}, legacy serial inner work inside each step, "
              f"iso_tol={VIEW_ISO_TOL}, dwbo={DWBO_THRESHOLD}, "
              f"metal_dwbo={METAL_DWBO_THRESHOLD}, "
              f"sym_wbo_tol={SYMMETRY_WBO_TOL}, "
              f"event_weight_power={EVENT_WEIGHT_POWER}, "
              f"wbo_progress_power={WBO_PROGRESS_POWER}, "
              f"xtb_mode={XTB_CACHE_MODE}, "
              f"xtb_threads={XTB_OMP_THREADS}, charge={XTB_CHARGE}, "
              f"multiplicity={XTB_MULTIPLICITY}, include_gt={INCLUDE_GT}, "
              f"save_alignment_files={args.save_alignment_files}, "
              f"resume_rp={args.resume_rp})")
        with mp.Pool(worker_budget) as pool:
            work = [
                (step, args.stage, 0, args.mechanism,
                 args.save_alignment_files, args.resume_rp)
                for step in steps
            ]
            for i, rec in enumerate(pool.imap_unordered(
                    _process_step_stage_star, work), 1):
                _record(i, rec)
    else:
        total_workers = worker_budget
        target_inner = (
            int(args.inner_workers) if args.inner_workers and args.inner_workers > 1
            else max(1, min(int(args.auto_inner_workers), total_workers))
        )
        outer_slots = max(1, min(len(steps), max(1, total_workers // target_inner)))
        inner_workers = max(1, total_workers // outer_slots)
        scheduled_steps = sorted(
            steps, key=_estimate_cut_sweep_units, reverse=True)
        print(f"Processing {len(steps)} steps in auto mode: "
              f"{outer_slots} concurrent steps x {inner_workers} "
              f"inner workers "
              f"(stage={args.stage}, total budget={total_workers}, "
              f"cut_sweep chunksize={CUTSWEEP_CHUNKSIZE}, "
              f"iso_tol={VIEW_ISO_TOL}, dwbo={DWBO_THRESHOLD}, "
              f"metal_dwbo={METAL_DWBO_THRESHOLD}, "
              f"sym_wbo_tol={SYMMETRY_WBO_TOL}, "
              f"event_weight_power={EVENT_WEIGHT_POWER}, "
              f"wbo_progress_power={WBO_PROGRESS_POWER}, "
              f"xtb_mode={XTB_CACHE_MODE}, "
              f"xtb_threads={XTB_OMP_THREADS}, charge={XTB_CHARGE}, "
              f"multiplicity={XTB_MULTIPLICITY}, include_gt={INCLUDE_GT}, "
              f"save_alignment_files={args.save_alignment_files}, "
              f"resume_rp={args.resume_rp})")
        with cf.ProcessPoolExecutor(max_workers=outer_slots) as executor:
            futures = {
                executor.submit(
                    process_step_stage, step, args.stage, inner_workers,
                    args.mechanism, args.save_alignment_files,
                    args.resume_rp): step
                for step in scheduled_steps
            }
            for i, fut in enumerate(cf.as_completed(futures), 1):
                step = futures[fut]
                try:
                    rec = fut.result()
                except Exception as e:
                    rec = {
                        "step": step,
                        "error": f"{type(e).__name__}: {e}",
                    }
                _record(i, rec)

    print(f"\n{n_ok} ok, {n_err} errors in {time.time()-t0:.0f}s")

    EVAL_JSON.write_text(json.dumps(eval_records))
    print(f"wrote {EVAL_JSON}  ({n_ok} step records)")


if __name__ == "__main__":
    main()

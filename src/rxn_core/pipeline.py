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
import copy
import concurrent.futures as cf
import json
import multiprocessing as mp
import os
import re
import shutil
import time
import traceback
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import numpy as np

from rxn_core import (parse_xyz, classify_bonds, parse_g98_modes,
                      core_atoms_in_R_frame,
                      reindex_modes_to_R, bond_overlap_per_mode,
                      build_graph, cut_sweep, cut_sweep_items,
                      merge_cut_sweep_pools, run_cut_sweep_chunk,
                      select_min_mechanisms, WeightedGraph,
                      match_weighted_subgraph)
from rxn_core.alignment.sweep import (
    attach_completed_candidate_groups,
    complete_chosen_automorphism_groups,
    run_no_cut_core_branch_records,
)
from rxn_core.alignment.index_chirality import (
    IndexChiralityConflict,
    analytical_family_static_context,
    analyze_group_chirality_branch,
    compile_analytical_mapping_family,
    fixed_mapping_aligned_rmsd,
    select_index_chirality_assignment,
)
from rxn_core.alignment.interpolation import (
    internal_coordinate_interpolation,
    proper_align_coordinates,
)
from rxn_core.alignment.post_aam import PostAAMMechanism
from rxn_core.chemistry_computations import (
    run_xtb, run_xtb_hess, write_xyz_str, xyz_with_disp,
)
from rxn_core.matcher import (
    _atom_tuple_orbit,
    _nauty_atom_generators,
    _nauty_orbits,
)
from rxn_core.smiles import smiles_to_formal_wbo

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
VIEWER_STATIC_DIR = Path(__file__).with_name("static")
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
# Compatibility name for serialized artifacts and older callers.  Pynauty
# must use the same tolerance as active-edge matching.
SYMMETRY_WBO_TOL = VIEW_ISO_TOL
INDEX_CHIRALITY = os.environ.get(
    "BGCP_INDEX_CHIRALITY", "off").strip().lower()
BGCP_TIMING = os.environ.get("BGCP_TIMING", "0") == "1"
INCLUDE_GT = os.environ.get("BGCP_INCLUDE_GT", "0").lower() in {
    "1", "true", "yes", "on"
}
SYMMETRY_REPAIR = os.environ.get("BGCP_SYMMETRY_REPAIR", "1") != "0"
SYMMETRY_REPAIR_MIN_CHANGES = int(os.environ.get("BGCP_SYMMETRY_REPAIR_MIN_CHANGES", "1"))
SYMMETRY_REPAIR_MAX_EVALS = int(os.environ.get("BGCP_SYMMETRY_REPAIR_MAX_EVALS", "20000"))
TS_ALIGN_GRAPH_FLOOR = float(os.environ.get("BGCP_TS_ALIGN_GRAPH_FLOOR", "0.2"))
TS_ALIGN_MAX_CORE_MAPS = int(os.environ.get("BGCP_TS_ALIGN_MAX_CORE_MAPS", "20000"))
PREFER_ENDPOINT_CONSENSUS = (
    os.environ.get("BGCP_PREFER_ENDPOINT_CONSENSUS", "1").lower()
    in {"1", "true", "yes", "on"}
)
AUTO_INNER_WORKERS = int(os.environ.get("BGCP_AUTO_INNER_WORKERS", "8"))
XTB_CACHE_MODE = os.environ.get("BGCP_XTB_MODE", "auto").lower()
XTB_OMP_THREADS = os.environ.get("BGCP_XTB_OMP_THREADS", "auto")
XTB_MAX_THREADS = int(os.environ.get("BGCP_XTB_MAX_THREADS", "8"))
XTB_WORKERS = os.environ.get("BGCP_XTB_WORKERS", "auto")
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
    energy_R: float | None = None
    energy_P: float | None = None


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
    energy: float | None = None


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
        'symmetry_wbo_tol': VIEW_ISO_TOL,
        'index_chirality': INDEX_CHIRALITY,
        'n_seeds': N_SEEDS_PER_RUN,
        'max_branches': VIEW_MAX_BRANCHES,
        'chunksize': CUTSWEEP_CHUNKSIZE,
        'symmetry_repair': SYMMETRY_REPAIR,
        'symmetry_repair_min_changes': SYMMETRY_REPAIR_MIN_CHANGES,
        'symmetry_repair_max_evals': SYMMETRY_REPAIR_MAX_EVALS,
        'anchor_map': {},
    }


def ts_stage_config():
    """Current TS/IG verification hypotheses as a serializable dict."""
    return {
        'iso_tol': VIEW_ISO_TOL,
        'graph_floor': TS_ALIGN_GRAPH_FLOOR,
        'dwbo_threshold': DWBO_THRESHOLD,
        'metal_dwbo_threshold': METAL_DWBO_THRESHOLD,
        'symmetry_wbo_tol': VIEW_ISO_TOL,
        'n_seeds': N_SEEDS_PER_RUN,
        'max_core_maps': TS_ALIGN_MAX_CORE_MAPS,
        'prefer_endpoint_consensus': PREFER_ENDPOINT_CONSENSUS,
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


def _available_cpus(default=1):
    for name in ("TSDISCO_SLURM_CPUS", "SLURM_CPUS_PER_TASK",
                 "SLURM_CPUS_ON_NODE"):
        value = os.environ.get(name)
        if value:
            try:
                return max(1, int(value))
            except ValueError:
                pass
    return max(1, int(default or os.cpu_count() or 1))


def _resolve_xtb_workers(inner_workers=0, workers=None):
    raw = XTB_WORKERS if workers is None else workers
    if isinstance(raw, str) and raw.lower() == "auto":
        threads = _resolve_xtb_threads()
        cpus = _available_cpus(default=inner_workers or os.cpu_count() or 1)
        cap_by_threads = max(1, cpus // max(1, threads))
        requested = int(inner_workers or cap_by_threads)
        return max(1, min(requested, cap_by_threads))
    return max(1, int(raw))


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


def _finite_float(value):
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


_FLOAT_RE = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
_ENERGY_PATTERNS = [
    re.compile(rf"\benergy\s*:\s*{_FLOAT_RE}", re.IGNORECASE),
    re.compile(rf"\btotal\s+energy\b\s*[:=]?\s*{_FLOAT_RE}", re.IGNORECASE),
    re.compile(rf"\bTOTAL\s+ENERGY\b\s*{_FLOAT_RE}"),
]


def _parse_energy_from_text(text):
    """Best-effort extraction of an xtb/XYZ energy in Hartree."""
    for pattern in _ENERGY_PATTERNS:
        match = pattern.search(text or "")
        if match:
            return _finite_float(match.group(1))
    return None


def _read_xyz_comment(path):
    try:
        with Path(path).open() as handle:
            handle.readline()
            return handle.readline().strip()
    except OSError:
        return ""


def _read_cache_energy(cache_dir):
    """Read cached xtb energy from XYZ comments or saved xtb output logs."""
    cache_dir = Path(cache_dir)
    if not cache_dir.exists():
        return None
    for xyz in sorted(cache_dir.glob("*.xyz")):
        energy = _parse_energy_from_text(_read_xyz_comment(xyz))
        if energy is not None:
            return energy
    preferred = [
        "xtb.stdout", "xtb.out", "xtb.log", "output.log",
        "xtbhess.stdout", "xtbhess.out", "g98.out",
    ]
    seen = set()
    candidates = []
    for name in preferred:
        path = cache_dir / name
        if path.exists():
            candidates.append(path)
            seen.add(path)
    for pattern in ("*.out", "*.log", "*.stdout"):
        for path in sorted(cache_dir.glob(pattern)):
            if path not in seen:
                candidates.append(path)
                seen.add(path)
    for path in candidates:
        try:
            energy = _parse_energy_from_text(path.read_text(errors="ignore"))
        except OSError:
            energy = None
        if energy is not None:
            return energy
    return None


def _frequency_list(freqs):
    if freqs is None:
        return None
    arr = np.asarray(freqs, float).reshape(-1)
    return [float(v) for v in arr if np.isfinite(v)]


def _frequency_summary(freqs):
    values = _frequency_list(freqs)
    if values is None:
        return None
    if not values:
        return {
            'n_modes': 0,
            'n_imaginary': 0,
            'imaginary_cm1': [],
            'lowest_cm1': None,
            'highest_cm1': None,
        }
    imaginary = [v for v in values if v < 0]
    return {
        'n_modes': len(values),
        'n_imaginary': len(imaginary),
        'imaginary_cm1': imaginary,
        'lowest_cm1': min(values),
        'highest_cm1': max(values),
    }


def _endpoint_metadata(label, energy):
    return {
        'label': label,
        'energy_hartree': _finite_float(energy),
        'energy_units': 'hartree',
    }


def _target_metadata(target):
    return {
        'energy_hartree': _finite_float(target.energy),
        'energy_units': 'hartree',
        'frequencies_cm1': _frequency_list(target.freqs),
        'frequency_summary': _frequency_summary(target.freqs),
        'frequency_units': 'cm^-1',
    }


def _attach_target_metadata(record, target):
    if record is None:
        return None
    record.update(_target_metadata(target))
    return record


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
                            step_dir=None, energy_R=None, energy_P=None):
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
        energy_R=_finite_float(energy_R),
        energy_P=_finite_float(energy_P),
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
    energy_R = _read_cache_energy(reactant_workdir)
    energy_P = _read_cache_energy(product_workdir)
    step_dir = (
        Path(workdir) if workdir is not None
        else Path(reactant_workdir).parent
    )
    return step_inputs_from_arrays(
        name, elR, xyzR, wboR, elP, xyzP, wboP, step_dir=step_dir,
        energy_R=energy_R, energy_P=energy_P)


def _weighted_graph_source_json(endpoint):
    return {
        "nodes": endpoint.nodes,
        "weights": np.asarray(endpoint.wbo, float).tolist(),
        "weight_name": "wbo",
        "coords": np.asarray(endpoint.coords, float).tolist(),
        "metadata": {
            "source": "smiles",
            "smiles": endpoint.smiles,
            "wbo_kind": "formal_bond_order",
            "hydrogen_policy": endpoint.hydrogen_policy,
            "atom_maps": {
                str(k): int(v) for k, v in sorted(endpoint.atom_maps.items())
            },
        },
    }


def write_smiles_source_files(workdir, reactant_endpoint, product_endpoint,
                              *, name="alignment"):
    """Write source/debug files for a formal-WBO SMILES R/P input pair."""
    source_dir = Path(workdir) / "source"
    source_dir.mkdir(parents=True, exist_ok=True)
    (source_dir / "R.xyz").write_text(write_xyz_str(
        reactant_endpoint.elements,
        reactant_endpoint.coords,
        f"{name} R 2D-from-CXSMILES"))
    (source_dir / "P.xyz").write_text(write_xyz_str(
        product_endpoint.elements,
        product_endpoint.coords,
        f"{name} P 2D-from-CXSMILES"))
    (source_dir / "R_wbo.json").write_text(json.dumps(
        np.asarray(reactant_endpoint.wbo, float).tolist(), indent=2))
    (source_dir / "P_wbo.json").write_text(json.dumps(
        np.asarray(product_endpoint.wbo, float).tolist(), indent=2))
    (source_dir / "R_graph.json").write_text(json.dumps(
        _weighted_graph_source_json(reactant_endpoint), indent=2))
    (source_dir / "P_graph.json").write_text(json.dumps(
        _weighted_graph_source_json(product_endpoint), indent=2))
    (source_dir / "atom_maps.json").write_text(json.dumps({
        "R_index_to_atom_map": {
            str(k): int(v)
            for k, v in sorted(reactant_endpoint.atom_maps.items())
        },
        "P_index_to_atom_map": {
            str(k): int(v)
            for k, v in sorted(product_endpoint.atom_maps.items())
        },
    }, indent=2))
    (source_dir / "smiles.json").write_text(json.dumps({
        "reactant_smiles": reactant_endpoint.smiles,
        "product_smiles": product_endpoint.smiles,
        "wbo_kind": "formal_bond_order",
        "hydrogen_policy": reactant_endpoint.hydrogen_policy,
        "coordinate_policy": "planar RDKit depiction for display only",
    }, indent=2))
    return source_dir


def smiles_inputs_from_strings(reactant_smiles, product_smiles, *,
                               name="alignment", workdir=None,
                               sanitize=True, component_spacing=3.0,
                               write_source=True,
                               expand_hydrogens=True):
    """Build Stage 1 inputs from SMILES/CXSMILES formal bond orders.

    Atom-map labels in CXSMILES are retained as metadata/source files only.
    They are not AAM anchors; pass ``anchor_map``/``--anchor`` separately when
    those labels should become hard mapping constraints.
    """
    r_endpoint = smiles_to_formal_wbo(
        reactant_smiles, sanitize=sanitize,
        component_spacing=component_spacing,
        expand_hydrogens=expand_hydrogens)
    p_endpoint = smiles_to_formal_wbo(
        product_smiles, sanitize=sanitize,
        component_spacing=component_spacing,
        expand_hydrogens=expand_hydrogens)
    step_dir = Path("." if workdir is None else workdir)
    if write_source and workdir is not None:
        write_smiles_source_files(
            step_dir, r_endpoint, p_endpoint, name=name)
    return step_inputs_from_arrays(
        name,
        r_endpoint.elements, r_endpoint.coords, r_endpoint.wbo,
        p_endpoint.elements, p_endpoint.coords, p_endpoint.wbo,
        step_dir=step_dir)


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
                          target_index=0, energy=None):
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
        energy=_finite_float(energy),
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
    energy = _read_cache_energy(sp_workdir)
    if energy is None:
        energy = _read_cache_energy(hess_workdir)
    return ts_target_from_arrays(
        kind, label, el, xyz, wbo, freqs, modes, target_index=target_index,
        energy=energy)


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


def discover_mechanisms_from_smiles(reactant_smiles, product_smiles, *,
                                    name="alignment", workdir=None,
                                    sanitize=True, component_spacing=3.0,
                                    config=None, inner_workers=0,
                                    return_inputs=False,
                                    expand_hydrogens=True):
    """SMILES/CXSMILES R-P mechanism discovery using formal bond orders."""
    inputs = smiles_inputs_from_strings(
        reactant_smiles, product_smiles,
        name=name, workdir=workdir, sanitize=sanitize,
        component_spacing=component_spacing,
        expand_hydrogens=expand_hydrogens)
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
    """Finalize ids after exact event-certificate grouping in AAM.

    Mechanism grouping now happens in the cut-sweep pool using a canonical
    colored-event certificate.  Re-grouping here by individual atom orbit IDs
    was weaker and discarded analytical branch families, so this outer/view
    stage is intentionally non-deduplicating.
    """
    deduped = list(mechanisms)
    for new_id, mech in enumerate(deduped, 1):
        mech.setdefault('dedup_source_ids', [int(mech.get('id', new_id))])
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


def _core_branch_record_signature(record, core_R):
    """Branch-state key for TS core records after conversion to R indexing."""
    core_R = tuple(sorted(int(r) for r in core_R))
    mapping = {int(r): int(t) for r, t in record.get('mapping', {}).items()}
    blocks = []
    block_r = set()
    for block in record.get('blocks', ()):
        r_atoms = tuple(sorted(
            int(r) for r in block.get('r_atoms', ())
            if int(r) in set(core_R)
        ))
        p_atoms = tuple(sorted(int(t) for t in block.get('p_atoms', ())))
        if not r_atoms:
            continue
        blocks.append((r_atoms, p_atoms))
        block_r.update(r_atoms)
    fixed = tuple(
        (int(r), int(mapping[r]))
        for r in core_R
        if r in mapping and r not in block_r
    )
    return fixed, tuple(sorted(blocks))


def _add_core_branch_record_entry(pool, core_R, record, source):
    core_R = tuple(sorted(int(r) for r in core_R))
    mapping = {int(r): int(t) for r, t in record.get('mapping', {}).items()
               if int(r) in set(core_R)}
    if len(mapping) != len(core_R):
        return
    blocks = []
    for block in record.get('blocks', ()):
        r_atoms = sorted(
            int(r) for r in block.get('r_atoms', ())
            if int(r) in set(core_R)
        )
        p_atoms = sorted(int(t) for t in block.get('p_atoms', ()))
        if r_atoms and len(p_atoms) >= len(r_atoms):
            blocks.append({'r_atoms': r_atoms, 'p_atoms': p_atoms})
    normalized = {
        'mapping': mapping,
        'blocks': blocks,
        'source': source,
        'dedup_count': int(record.get('dedup_count', 1)),
    }
    sig = _core_branch_record_signature(normalized, core_R)
    entry = pool.get(sig)
    if entry is None:
        pool[sig] = {
            'records': [normalized],
            'sources': {source},
            'dedup_count': int(normalized['dedup_count']),
        }
        return
    entry['sources'].add(source)
    entry['dedup_count'] = (
        int(entry.get('dedup_count', 1)) + int(normalized['dedup_count'])
    )
    record_key = (
        tuple(sorted(normalized['mapping'].items())),
        tuple(sorted(
            (tuple(block['r_atoms']), tuple(block['p_atoms']))
            for block in normalized['blocks']
        )),
        source,
    )
    existing_keys = {
        (
            tuple(sorted(rec.get('mapping', {}).items())),
            tuple(sorted(
                (tuple(block.get('r_atoms', ())),
                 tuple(block.get('p_atoms', ())))
                for block in rec.get('blocks', ())
            )),
            rec.get('source'),
        )
        for rec in entry['records']
    }
    if record_key not in existing_keys:
        entry['records'].append(normalized)


def _product_core_branch_records_to_reactant(records_PT, mapping_RP, core_R):
    """Convert P-indexed branch records to R-indexed core branch records."""
    inv_RP = {int(p): int(r) for r, p in mapping_RP.items()}
    core_set = {int(r) for r in core_R}
    converted = []
    for record in records_PT or ():
        pulled = {}
        for p, t in record.get('mapping', {}).items():
            r = inv_RP.get(int(p))
            if r in core_set:
                pulled[r] = int(t)
        if len(pulled) != len(core_set):
            continue
        blocks = []
        for block in record.get('blocks', ()):
            r_atoms = sorted({
                inv_RP[int(p)]
                for p in block.get('r_atoms', ())
                if int(p) in inv_RP and inv_RP[int(p)] in core_set
            })
            p_atoms = sorted(int(t) for t in block.get('p_atoms', ()))
            if r_atoms and len(p_atoms) >= len(r_atoms):
                blocks.append({'r_atoms': r_atoms, 'p_atoms': p_atoms})
        converted.append({
            'mapping': pulled,
            'blocks': blocks,
            'dedup_count': int(record.get('dedup_count', 1)),
        })
    return converted


def _merge_endpoint_core_branch_records(core_R, r_records, p_records_as_r):
    """Merge R-TS and pulled-back P-TS branch states before shuffling."""
    merged = {}
    for record in r_records or ():
        _add_core_branch_record_entry(merged, core_R, record, 'R')
    for record in p_records_as_r or ():
        _add_core_branch_record_entry(merged, core_R, record, 'P')
    return merged


def _ts_endpoint_pool_task(task):
    """Worker task for one endpoint-side TS no-cut core alignment.

    A task is one `(target TS, mechanism, endpoint)` alignment.  The caller
    merges the returned R->TS and P->TS pools and then scores them in the main
    process so ranking stays deterministic.
    """
    t0 = time.time()
    branch_records = run_no_cut_core_branch_records(
        task['elS'], task['wboS'],
        task['elT'], task['wboT'],
        task['core_S'],
        graph_floor=task['graph_floor'],
        iso_tol=task['iso_tol'],
        dwbo_threshold=task['dwbo_threshold'],
        metal_dwbo_threshold=task['metal_dwbo_threshold'],
        symmetry_wbo_tol=task['symmetry_wbo_tol'],
        n_seeds=task['n_seeds'],
        max_branches=task['max_core_maps'],
    )
    return {
        'key': task['key'],
        'target_order': int(task['target_order']),
        'target_label': task['target_label'],
        'mech_id': int(task['mech_id']),
        'mech_pos': int(task['mech_pos']),
        'endpoint': task['endpoint'],
        'branch_records': branch_records,
        'n_pool': len(branch_records),
        'core_size': len(task['core_S']),
        'elapsed': time.time() - t0,
        'hit_cap': len(branch_records) >= task['max_core_maps'],
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
            # Keep both native target indexing and the selected R-indexed
            # core-only materialization.  The viewer builds its own
            # geometry-preserving row permutation from the native target rows;
            # it must not render this core-only materialization as a molecule.
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


def _has_endpoint_consensus(score):
    return set(score.get('core_sources') or ()) >= {'R', 'P'}


def _better_score(candidate, incumbent):
    return incumbent is None or candidate['S'] > incumbent['S']


def _core_branch_record_automorphism_variants(
        record, core_R, g_T, symmetry_wbo_tol, generator_cache):
    """Materialize exact core maps by strict TS automorphism action.

    Fixed core images are singleton-colored.  Each explicit branch block's
    target pool gets a branch-local color, so nauty gives the subgroup that
    preserves this branch state.  The resulting tuple orbit can only move core
    atoms through valid automorphisms inside those branch pools.
    """
    core_R = tuple(sorted(int(r) for r in core_R))
    core_set = set(core_R)
    base = {int(r): int(t) for r, t in record.get('mapping', {}).items()
            if int(r) in core_set}
    if len(base) != len(core_R):
        return []

    blocks = []
    covered = set()
    for block in record.get('blocks', ()):
        r_atoms = tuple(sorted(
            int(r) for r in block.get('r_atoms', ())
            if int(r) in core_set and int(r) not in covered
        ))
        if not r_atoms:
            continue
        p_atoms = tuple(sorted(int(t) for t in block.get('p_atoms', ())))
        if len(p_atoms) < len(r_atoms):
            continue
        blocks.append((r_atoms, p_atoms))
        covered.update(r_atoms)

    fixed = {r: base[r] for r in core_R if r not in covered}
    if not blocks:
        return [base]

    tag_parts = {}

    def add_tag(atom, tag):
        tag_parts.setdefault(int(atom), []).append(tag)

    for r, t in fixed.items():
        add_tag(t, ('fixed', int(r), int(t)))
    block_sets = []
    for idx, (r_atoms, p_atoms) in enumerate(blocks):
        p_set = frozenset(int(t) for t in p_atoms)
        block_sets.append((tuple(r_atoms), p_set))
        for t in p_atoms:
            add_tag(t, ('block', int(idx)))

    atom_color_tags = {
        atom: tuple(parts)
        for atom, parts in tag_parts.items()
    }
    tag_key = tuple(sorted(atom_color_tags.items()))
    generators = generator_cache.get(tag_key)
    if generators is None:
        generators = _nauty_atom_generators(
            g_T, wbo_tol=float(symmetry_wbo_tol),
            atom_color_tags=atom_color_tags)
        generator_cache[tag_key] = generators

    seed = tuple(base[r] for r in core_R)
    out = {}
    for values in _atom_tuple_orbit(seed, generators):
        if len(set(values)) != len(values):
            continue
        mapping = {r: int(t) for r, t in zip(core_R, values)}
        if any(mapping[r] != t for r, t in fixed.items()):
            continue
        ok = True
        for r_atoms, p_set in block_sets:
            if any(mapping[r] not in p_set for r in r_atoms):
                ok = False
                break
        if not ok:
            continue
        out[_core_pool_key(mapping, core_R)] = mapping
    return list(out.values())


def best_under_mech_using_branch_pool(elR, xyzR, wboR, wboP,
                                      elT, xyzT, wboT, freqs, modes_TS,
                                      branch_pool, mapping_RP,
                                      broken_R, formed_R, core_R,
                                      score_weights=None,
                                      prefer_endpoint_consensus=True,
                                      symmetry_wbo_tol=0.2):
    """Score exact core maps after endpoint branch-state merge."""
    g_T = build_graph(elT, wboT, bond_cut=0.2)
    generator_cache = {}
    exact_pool = {}
    for entry in branch_pool.values():
        for record in entry.get('records', ()):
            source = record.get('source', 'R')
            variants = _core_branch_record_automorphism_variants(
                record, core_R, g_T, symmetry_wbo_tol, generator_cache)
            for mapping in variants:
                key = _core_pool_key(mapping, core_R)
                item = exact_pool.get(key)
                if item is None:
                    exact_pool[key] = {
                        'mapping': mapping,
                        'sources': {source},
                        'dedup_count': int(record.get('dedup_count', 1)),
                    }
                else:
                    item['sources'].add(source)
                    item['dedup_count'] = (
                        int(item.get('dedup_count', 1))
                        + int(record.get('dedup_count', 1))
                    )

    return best_under_mech_using_pool(
        elR, xyzR, wboR, wboP,
        elT, xyzT, wboT, freqs, modes_TS,
        exact_pool, mapping_RP, broken_R, formed_R, core_R,
        score_weights=score_weights,
        prefer_endpoint_consensus=prefer_endpoint_consensus)


def best_under_mech_using_pool(elR, xyzR, wboR, wboP,
                                 elT, xyzT, wboT, freqs, modes_TS,
                                 rt_pool, mapping_RP,
                                 broken_R, formed_R, core_R,
                                 score_weights=None,
                                 prefer_endpoint_consensus=True):
    """Score every R-frame core witness under one mech.

    `rt_pool` may contain native R->TS candidates and P->TS candidates that
    were pulled back through the R-P mechanism witness.  Two witnesses that
    agree on `core_R -> TS` are score-equivalent for this mechanism, so score
    one representative per exact core map.  When configured, prefer the best
    consensus map first; if no consensus map exists, fall back to the highest
    S over the endpoint-union pool.
    """
    core_R_set = frozenset(core_R)
    seen_core = set()
    best = None
    best_consensus = None
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
            s['endpoint_consensus'] = _has_endpoint_consensus(s)
        if s and _better_score(s, best):
            best = s
        if (s and prefer_endpoint_consensus and s['endpoint_consensus']
                and _better_score(s, best_consensus)):
            best_consensus = s
    return best_consensus or best


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

    Coordinates are only reindexed into the R atom frame.  No spatial fitting
    is applied.  Each mechanism directory is self-contained and
    contains an R endpoint, the mechanism-specific aligned P endpoint, a
    two-frame XYZ, a mapping CSV, and a JSON metadata file.
    """
    out_dir = Path(out_dir) if out_dir is not None else alignment_output_dir(
        inputs.step_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    mech_root = out_dir / "mechanisms"
    mech_root.mkdir(exist_ok=True)
    endpoint_metadata = {
        'reactant': _endpoint_metadata('R', inputs.energy_R),
        'product': _endpoint_metadata('P', inputs.energy_P),
    }

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
            'endpoint_metadata': endpoint_metadata,
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
                'order using mapping_RP; no spatial fitting is used.'
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
        'endpoint_metadata': endpoint_metadata,
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
    endpoint_metadata = {
        'reactant': _endpoint_metadata('R', inputs.energy_R),
        'product': _endpoint_metadata('P', inputs.energy_P),
    }

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
                'energy_hartree': score.get('energy_hartree'),
                'energy_units': score.get('energy_units', 'hartree'),
                'frequency_units': score.get('frequency_units', 'cm^-1'),
                'frequencies_cm1': score.get('frequencies_cm1'),
                'frequency_summary': score.get('frequency_summary'),
                'event_terms': score.get('event_terms', []),
                'core_map_R_to_target': score.get('core_map', {}),
                'core_sources': score.get('core_sources', []),
                'endpoint_consensus': score.get('endpoint_consensus'),
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
                'energy_hartree': score.get('energy_hartree'),
                'frequency_summary': score.get('frequency_summary'),
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
        'endpoint_metadata': endpoint_metadata,
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


def _rp_cfg(config=None):
    return dict(rp_stage_config() if config is None else config)


def _rp_cut_kwargs(cfg):
    return {
        'cut_floor': cfg.get('cut_floor', CUT_FLOOR),
        'graph_floor': cfg.get('graph_floor', 0.2),
        'iso_tol': cfg.get('iso_tol', VIEW_ISO_TOL),
        'dwbo_threshold': cfg.get('dwbo_threshold', DWBO_THRESHOLD),
        'metal_dwbo_threshold': cfg.get(
            'metal_dwbo_threshold', METAL_DWBO_THRESHOLD),
        'symmetry_wbo_tol': cfg.get('iso_tol', VIEW_ISO_TOL),
        'n_seeds': cfg.get('n_seeds', N_SEEDS_PER_RUN),
        'max_branches': cfg.get('max_branches', VIEW_MAX_BRANCHES),
        'chunksize': cfg.get('chunksize', CUTSWEEP_CHUNKSIZE),
        'symmetry_repair': cfg.get('symmetry_repair', SYMMETRY_REPAIR),
        'symmetry_repair_min_changes': cfg.get(
            'symmetry_repair_min_changes', SYMMETRY_REPAIR_MIN_CHANGES),
        'symmetry_repair_max_evals': cfg.get(
            'symmetry_repair_max_evals', SYMMETRY_REPAIR_MAX_EVALS),
        'anchor_map': cfg.get('anchor_map') or {},
    }


_ANALYTICAL_FAMILY_WORKER_CONTEXT = None
_ANALYTICAL_BRANCH_EVAL_CONTEXT = None


def _init_analytical_family_worker(context):
    global _ANALYTICAL_FAMILY_WORKER_CONTEXT
    _ANALYTICAL_FAMILY_WORKER_CONTEXT = context


def _init_analytical_branch_eval_worker(context):
    global _ANALYTICAL_BRANCH_EVAL_CONTEXT
    _ANALYTICAL_BRANCH_EVAL_CONTEXT = context


def _evaluate_analytical_branch_task(branch_index):
    """Evaluate one maximal family independently for chirality and RMSD."""
    context = _ANALYTICAL_BRANCH_EVAL_CONTEXT
    branch = context['branches'][branch_index]
    inputs = context['inputs']
    cfg = context['cfg']
    branch_mapping = _int_mapping(branch['mapping'])
    branch_hierarchy = branch.get('hierarchy') or {}
    branch_group_chirality = None
    branch_index_chirality = None
    try:
        if context['index_chirality_mode'] == 'preserve':
            branch_group_chirality = analyze_group_chirality_branch(
                branch_mapping,
                inputs.elR, inputs.xyzR, inputs.wboR,
                inputs.elP, inputs.xyzP, inputs.wboP,
                graph_floor=cfg.get('graph_floor', 0.2),
            )
            selection = select_index_chirality_assignment(
                branch_mapping, branch_hierarchy,
                inputs.elR, inputs.xyzR, inputs.wboR,
                inputs.elP, inputs.xyzP, inputs.wboP,
                graph_floor=cfg.get('graph_floor', 0.2),
                symmetry_wbo_tol=cfg.get('iso_tol', VIEW_ISO_TOL),
                dwbo_threshold=cfg.get('dwbo_threshold', DWBO_THRESHOLD),
                metal_dwbo_threshold=cfg.get(
                    'metal_dwbo_threshold', METAL_DWBO_THRESHOLD),
                anchor_map=cfg.get('anchor_map') or {},
                group_chirality_frames=(
                    branch_group_chirality.defined_frames),
                static_context=context['static_context'],
            )
            branch_mapping = selection.selected_mapping
            branch_index_chirality = selection.metadata
            branch_index_chirality['group_chirality_branch'] = (
                branch_group_chirality.metadata)
        rmsd = fixed_mapping_aligned_rmsd(
            branch_mapping, inputs.xyzR, inputs.xyzP)
        return ('ok', (
            float(rmsd), tuple(branch_mapping[r]
                               for r in sorted(branch_mapping)),
            branch_index, branch_mapping, branch_hierarchy,
            branch_group_chirality, branch_index_chirality,
        ))
    except IndexChiralityConflict as exc:
        return ('failure', {
            'branch_index': int(branch_index),
            'reason': str(exc),
            'source_mapping_RP': dict(branch_mapping),
            'diagnostics': getattr(exc, 'diagnostics', None),
        })


def _compile_analytical_family_task(payload):
    """Pickle-safe exact family compiler for post-AAM process batches."""
    mapping, hierarchy = payload
    context = _ANALYTICAL_FAMILY_WORKER_CONTEXT
    if context is None:
        raise RuntimeError("analytical family compiler context is unset")
    return compile_analytical_mapping_family(
        mapping, hierarchy,
        context['elements_R'], context['wbo_R'],
        context['elements_P'], context['wbo_P'],
        graph_floor=context['graph_floor'],
        symmetry_wbo_tol=context['symmetry_wbo_tol'],
        dwbo_threshold=context['dwbo_threshold'],
        metal_dwbo_threshold=context['metal_dwbo_threshold'],
        anchor_map=context['anchor_map'],
        static_context=context['static_context'],
    )


def _analytical_branch_payload_key(branch):
    """Exact relation input consumed by analytical-family compilation."""
    mapping = _int_mapping(branch['mapping'])
    hierarchy = branch.get('hierarchy') or {}
    fragments = tuple(sorted(
        tuple(sorted(int(atom) for atom in fragment.get('fragment', ())))
        for fragment in hierarchy.get('fragments') or ({
            'fragment': sorted(mapping),
        },)
    ))
    return tuple(sorted(mapping.items())), fragments


def _merge_stored_fragment_groups(kept_branch, incoming_branch):
    """Union exact generators while collapsing an identical branch input."""
    kept_fragments = {
        tuple(sorted(map(int, fragment.get('fragment') or ()))): fragment
        for fragment in (kept_branch.get('hierarchy') or {}).get(
            'fragments') or ()
    }
    for incoming in (incoming_branch.get('hierarchy') or {}).get(
            'fragments') or ():
        key = tuple(sorted(map(int, incoming.get('fragment') or ())))
        kept = kept_fragments.get(key)
        if kept is None:
            continue
        kept_symmetry = kept.setdefault('symmetry', {})
        incoming_symmetry = incoming.get('symmetry') or {}
        if ('automorph_generators' not in kept_symmetry
                and 'automorph_generators' not in incoming_symmetry):
            continue
        generators = {
            tuple(map(int, generator))
            for generator in kept_symmetry.get(
                'automorph_generators') or ()
        }
        generators.update(
            tuple(map(int, generator))
            for generator in incoming_symmetry.get(
                'automorph_generators') or ())
        kept_symmetry['automorph_generators'] = [
            list(generator) for generator in sorted(generators)]


def _compile_analytical_families(inputs, branches, cfg, static_context=None):
    """Compile independent exact cosets in parallel when CPUs are available."""
    context = {
        'elements_R': inputs.elR,
        'wbo_R': inputs.wboR,
        'elements_P': inputs.elP,
        'wbo_P': inputs.wboP,
        'graph_floor': cfg.get('graph_floor', 0.2),
        'symmetry_wbo_tol': cfg.get('iso_tol', VIEW_ISO_TOL),
        'dwbo_threshold': cfg.get('dwbo_threshold', DWBO_THRESHOLD),
        'metal_dwbo_threshold': cfg.get(
            'metal_dwbo_threshold', METAL_DWBO_THRESHOLD),
        'anchor_map': cfg.get('anchor_map') or {},
    }
    context['static_context'] = static_context or (
        analytical_family_static_context(
            inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
            graph_floor=context['graph_floor'],
            dwbo_threshold=context['dwbo_threshold'],
            metal_dwbo_threshold=context['metal_dwbo_threshold']))
    unique_payloads = []
    branch_payload_indices = []
    payload_index = {}
    for branch in branches:
        mapping = _int_mapping(branch['mapping'])
        hierarchy = branch.get('hierarchy') or {}
        # _masked_relation_data consumes exactly the mapping and fragment
        # partition.  Fragment order/indices are categorical names only, so
        # canonical sorting removes growth-history duplicates without
        # weakening any relational constraint.
        key = _analytical_branch_payload_key(branch)
        index = payload_index.get(key)
        if index is None:
            index = len(unique_payloads)
            payload_index[key] = index
            unique_payloads.append((mapping, hierarchy))
        branch_payload_indices.append(index)
    requested = int(cfg.get('post_aam_workers') or _available_cpus(default=1))
    # Large branch sets are relation-compilation bound and each payload is
    # independent.  Use the full 48-core node rather than leaving one third of
    # the allocation idle; smaller sets retain the lower process cap.
    worker_cap = 48 if len(unique_payloads) >= 128 else 8
    workers = min(len(unique_payloads), max(1, requested), worker_cap)
    # A process cannot create children when this pipeline itself is running as
    # a daemonic pool worker.  Sequential execution is the same exact
    # computation, only a scheduling difference.
    if (workers <= 1 or len(unique_payloads) < 16
            or mp.current_process().daemon):
        _init_analytical_family_worker(context)
        compiled = [
            _compile_analytical_family_task(item)
            for item in unique_payloads
        ]
    else:
        with mp.get_context('fork').Pool(
                processes=workers,
                initializer=_init_analytical_family_worker,
                initargs=(context,)) as pool:
            compiled = pool.map(
                _compile_analytical_family_task, unique_payloads)
    return [compiled[index] for index in branch_payload_indices]


def rp_cut_work_items(inputs, config=None):
    """Return independent R-P no-cut/one-edge cut work items."""
    cfg = _rp_cfg(config)
    return cut_sweep_items(inputs.wboR, cfg.get('cut_floor', CUT_FLOOR))


def run_rp_cut_chunk(inputs, cuts, config=None, inner_workers=0,
                     trace_path=None):
    """Run one R-P cut chunk and return a partial mechanism pool."""
    cfg = _rp_cfg(config)
    t0 = time.time()
    pool = run_cut_sweep_chunk(
        inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
        cuts, n_workers=max(1, int(inner_workers or 1)),
        trace_path=trace_path,
        **_rp_cut_kwargs(cfg))
    return {
        'stage': 'rp_cut_chunk',
        'step': inputs.step_name,
        'config': cfg,
        'cuts': cuts,
        'pool': pool,
        'trace_path': str(trace_path) if trace_path else None,
        'timing': {'rp_cut_chunk_seconds': time.time() - t0},
    }


def _dedupe_analytical_mapping_families(
        inputs, branches, cfg, static_context=None):
    """Keep maximal exact mapping cosets and retain subsumed provenance."""
    # Growth paths that have the same mapping and fragment partition compile
    # to the exact same colored relation and therefore the same coset.  The
    # former implementation compiled them once but expanded the shared object
    # back to every raw path before containment, recreating an O(raw*families)
    # loop.  Quotient those identical relation inputs first and carry all path
    # provenance on their single exact representative.
    payload_groups = {}
    for source_index, raw in enumerate(branches):
        key = _analytical_branch_payload_key(raw)
        group = payload_groups.get(key)
        branch = copy.deepcopy(raw)
        provenance = {
            'source_branch_index': int(source_index),
            'cuts': [list(map(int, cut))
                     for cut in branch.get('cuts') or ()],
            'encounter_count': int(branch.get('encounter_count', 1)),
            'fragment_count': len(
                (branch.get('hierarchy') or {}).get('fragments') or ()),
        }
        if group is None:
            payload_groups[key] = {
                'branch': branch,
                'provenance': [provenance],
            }
            continue
        kept = group['branch']
        _merge_stored_fragment_groups(kept, branch)
        kept['encounter_count'] = (
            int(kept.get('encounter_count', 1))
            + int(branch.get('encounter_count', 1)))
        cuts = {
            tuple(map(int, cut)) for cut in kept.get('cuts') or ()
        } | {
            tuple(map(int, cut)) for cut in branch.get('cuts') or ()
        }
        kept['cuts'] = [list(cut) for cut in sorted(cuts)]
        group['provenance'].append(provenance)

    grouped = list(payload_groups.values())
    grouped_branches = [group['branch'] for group in grouped]
    compiled = _compile_analytical_families(
        inputs, grouped_branches, cfg, static_context=static_context)
    entries = []
    for group, family in zip(grouped, compiled):
        branch = dict(group['branch'])
        branch['path_provenance'] = list(group['provenance'])
        branch['covered_path_count'] = len(branch['path_provenance'])
        branch['mapping_family'] = family.record()
        entries.append({'branch': branch, 'family': family})

    def merge_entry(kept, removed):
        kept_branch = kept['branch']
        removed_branch = removed['branch']
        kept_branch['encounter_count'] = (
            int(kept_branch.get('encounter_count', 1))
            + int(removed_branch.get('encounter_count', 1)))
        cuts = {
            tuple(map(int, cut))
            for cut in kept_branch.get('cuts') or ()
        } | {
            tuple(map(int, cut))
            for cut in removed_branch.get('cuts') or ()
        }
        kept_branch['cuts'] = [list(cut) for cut in sorted(cuts)]
        kept_branch.setdefault('path_provenance', []).extend(
            removed_branch.get('path_provenance') or ())
        kept_branch['covered_path_count'] = len(
            kept_branch['path_provenance'])

    # Phase 1: equal cosets normally share the same colored-relation
    # certificates and exact group invariant.  Resolve those dense duplicate
    # classes locally instead of comparing every new path with every maximal
    # family found so far.  A certificate mismatch is never treated as proof
    # of inequality; phase 2 still performs authoritative subset checks.
    equality_buckets = defaultdict(list)
    equality_unique = []
    for entry in entries:
        bucket = equality_buckets[entry['family'].equivalence_bucket]
        equivalent = next((
            candidate for candidate in bucket
            if entry['family'].equivalent(candidate['family'])
        ), None)
        if equivalent is not None:
            merge_entry(equivalent, entry)
            continue
        bucket.append(entry)
        equality_unique.append(entry)

    # Phase 2: visit larger groups first.  A later family cannot contain an
    # earlier strictly larger group, so maximal-coset reduction needs only one
    # directional subset test.  Equal groups missed by the fast certificate
    # bucket are also resolved here because coset inclusion at equal order is
    # equality.
    def descending_group_order(entry):
        mantissa, exponent = entry['family'].group_order
        return -(np.log10(mantissa) + exponent)

    maximal = []
    for entry in sorted(equality_unique, key=descending_group_order):
        covering = next((
            candidate for candidate in maximal
            if entry['family'].is_subset_of(candidate['family'])
        ), None)
        if covering is not None:
            merge_entry(covering, entry)
            continue
        maximal.append(entry)
    return [item['branch'] for item in maximal]


def run_rp_stage_from_pool(inputs, pool, config=None, elapsed=None):
    """Finalize Stage 1 from a full or merged cut-sweep pool."""
    post_start = time.time()
    phase_seconds = defaultdict(float)
    cfg = _rp_cfg(config)
    phase_start = time.time()
    rp_min = select_min_mechanisms(pool)
    phase_seconds['mechanism_selection_seconds'] += time.time() - phase_start
    if not rp_min:
        raise RuntimeError("no min-bond mechanism")

    mechanisms = []
    graph_floor = cfg.get('graph_floor', 0.2)
    g_R_full = build_graph(inputs.elR, inputs.wboR, bond_cut=graph_floor)
    g_P_full = build_graph(inputs.elP, inputs.wboP, bond_cut=graph_floor)
    index_chirality_mode = str(
        cfg.get('index_chirality', 'off')).lower()
    if index_chirality_mode not in {'off', 'preserve'}:
        raise ValueError(
            "index_chirality must be 'off' or 'preserve'")
    rejected_index_chirality = []
    phase_start = time.time()
    analytical_static_context = analytical_family_static_context(
        inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
        graph_floor=cfg.get('graph_floor', 0.2),
        dwbo_threshold=cfg.get('dwbo_threshold', DWBO_THRESHOLD),
        metal_dwbo_threshold=cfg.get(
            'metal_dwbo_threshold', METAL_DWBO_THRESHOLD))
    phase_seconds['analytical_context_seconds'] += time.time() - phase_start
    for mi, (_sig, info) in enumerate(rp_min.items(), 1):
        raw_analytical_branches = list(info.get('branches') or [{
            'mapping': info['mapping'],
            'hierarchy': info.get('branch_symmetry') or {},
            'encounter_count': int(info.get('dedup_count', 1)),
        }])
        raw_analytical_branches = attach_completed_candidate_groups(
            raw_analytical_branches, g_P_full,
            wbo_tol=cfg.get('iso_tol', VIEW_ISO_TOL))
        phase_start = time.time()
        analytical_branches = _dedupe_analytical_mapping_families(
            inputs, raw_analytical_branches, cfg,
            static_context=analytical_static_context)
        phase_seconds['analytical_family_dedupe_seconds'] += (
            time.time() - phase_start)
        analytical_info = dict(info)
        analytical_info['branches'] = analytical_branches
        analytical_info['mapping'] = dict(analytical_branches[0]['mapping'])
        analytical_info['branch_symmetry'] = (
            analytical_branches[0].get('hierarchy') or {})
        phase_start = time.time()
        post_aam = PostAAMMechanism.from_aam_graphs(
            _sig, analytical_info, g_R_full, g_P_full,
            symmetry_wbo_tolerance=cfg.get('iso_tol', VIEW_ISO_TOL))
        phase_seconds['post_aam_model_seconds'] += time.time() - phase_start
        mapping_RP = _int_mapping(analytical_info['mapping'])
        raw_branch_symmetry = analytical_info['branch_symmetry']
        group_chirality = None
        index_chirality = None
        evaluated_branches = []
        branch_failures = []
        eval_context = {
            'branches': analytical_branches,
            'inputs': inputs,
            'cfg': cfg,
            'index_chirality_mode': index_chirality_mode,
            'static_context': analytical_static_context,
        }
        eval_workers = min(
            len(analytical_branches), _available_cpus(default=1), 8)
        phase_start = time.time()
        if (eval_workers > 1 and len(analytical_branches) > 1
                and not mp.current_process().daemon):
            with mp.get_context('fork').Pool(
                    processes=eval_workers,
                    initializer=_init_analytical_branch_eval_worker,
                    initargs=(eval_context,)) as eval_pool:
                branch_results = eval_pool.map(
                    _evaluate_analytical_branch_task,
                    range(len(analytical_branches)))
        else:
            _init_analytical_branch_eval_worker(eval_context)
            branch_results = [
                _evaluate_analytical_branch_task(index)
                for index in range(len(analytical_branches))
            ]
        phase_seconds['chirality_rmsd_seconds'] += time.time() - phase_start
        phase_start = time.time()
        for status, record in branch_results:
            if status == 'ok':
                evaluated_branches.append(record)
            else:
                branch_failures.append(record)
        if not evaluated_branches:
            failure_reason = (
                branch_failures[0]['reason']
                if len(branch_failures) == 1
                else 'no analytical AAM branch satisfies chirality')
            rejected_index_chirality.append({
                'source_mechanism_id': int(mi),
                'reason': failure_reason,
                'source_mapping_RP': dict(mapping_RP),
                'analytical_branch_failures': branch_failures,
            })
            continue
        evaluated_branches.sort(key=lambda item: item[:3])
        (selected_rmsd, _mapping_key_for_rank, selected_branch_index,
         mapping_RP, raw_branch_symmetry, group_chirality,
         index_chirality) = evaluated_branches[0]
        selected_fragments = (
            (analytical_branches[selected_branch_index].get('hierarchy') or {})
            .get('fragments') or ())
        if any('automorph_generators' not in (fragment.get('symmetry') or {})
               for fragment in selected_fragments):
            raise RuntimeError(
                "selected completed AAM branch lacks fragment groups")
        exact_target_generators = [
            generator
            for fragment in selected_fragments
            for generator in ((fragment.get('symmetry') or {})
                              .get('automorph_generators') or ())
        ]
        branch_symmetry = complete_chosen_automorphism_groups(
            raw_branch_symmetry, mapping_RP,
            g_R_full, g_P_full, cfg.get('iso_tol', VIEW_ISO_TOL),
            exact_target_generators=exact_target_generators)
        if index_chirality is not None:
            index_chirality['selected_analytical_branch_index'] = int(
                selected_branch_index)
            index_chirality['selected_fixed_mapping_aligned_rmsd'] = float(
                selected_rmsd)
            index_chirality['analytical_branch_count'] = len(
                analytical_branches)
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
        post_aam_record = post_aam.symmetry_record()
        post_aam_record['selected_mapping'] = dict(mapping_RP)
        post_aam_record['selected_analytical_branch_index'] = int(
            selected_branch_index)
        post_aam_record['selected_fixed_mapping_aligned_rmsd'] = float(
            selected_rmsd)
        post_aam_record['branch_failures'] = branch_failures
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
            'post_aam': post_aam_record,
            'branch_symmetry': branch_symmetry,
            'index_chirality': index_chirality,
        })
        phase_seconds['mechanism_materialization_seconds'] += (
            time.time() - phase_start)

    if index_chirality_mode == 'preserve' and not mechanisms:
        reasons = "; ".join(
            f"#{item['source_mechanism_id']}: {item['reason']}"
            for item in rejected_index_chirality)
        raise IndexChiralityConflict(
            f"all {len(rp_min)} minimum-event mechanisms failed "
            f"index-direction consensus ({reasons})",
            diagnostics={
                'rejected_mechanisms': rejected_index_chirality,
            })
    phase_start = time.time()
    r_orbits = _nauty_orbits(
        build_graph(inputs.elR, inputs.wboR, bond_cut=0.2),
        wbo_tol=cfg.get('iso_tol', VIEW_ISO_TOL))
    mechanisms = dedupe_mechanisms_by_bond_changes(mechanisms, r_orbits)
    for mechanism in mechanisms:
        if mechanism.get('index_chirality'):
            mechanism['branch_symmetry']['index_chirality'] = dict(
                mechanism['index_chirality'])
    phase_seconds['final_symmetry_dedupe_seconds'] += time.time() - phase_start
    cut_sweep_seconds = float(elapsed or 0.0)
    post_aam_seconds = time.time() - post_start
    return {
        'stage': 'rp',
        'step': inputs.step_name,
        'n_atoms': len(inputs.elR),
        'config': cfg,
        'mechanisms': [_mechanism_for_view(m) for m in mechanisms],
        'rejected_index_chirality': rejected_index_chirality,
        'timing': {
            'rp_seconds': cut_sweep_seconds + post_aam_seconds,
            'cut_sweep_seconds': cut_sweep_seconds,
            'post_aam_seconds': post_aam_seconds,
            **{key: float(value) for key, value in phase_seconds.items()},
        },
    }


def merge_rp_cut_chunks(inputs, chunks, config=None):
    """Merge partial R-P cut chunks and finalize the Stage 1 result."""
    t0 = time.time()
    pool = merge_cut_sweep_pools([chunk.get('pool', {}) for chunk in chunks])
    elapsed = sum(
        float((chunk.get('timing') or {}).get('rp_cut_chunk_seconds', 0.0))
        for chunk in chunks
    ) or (time.time() - t0)
    return run_rp_stage_from_pool(
        inputs, pool, config=config, elapsed=elapsed)


def run_rp_stage(inputs, config=None, inner_workers=0):
    """Stage 1: discover mechanism-dependent R-P alignments.

    This is the reusable alignment/mechanism-discovery entry point.  It runs
    the no-cut plus one-edge cut sweep on R->P, selects minimum bond-change
    mechanisms, classifies broken/forming bonds, and stores the selected R-P
    bijection together with every retained analytical AAM branch
    needed for later TS verification.
    """
    cfg = _rp_cfg(config)
    t0 = time.time()
    pool = cut_sweep(
        inputs.elR, inputs.wboR, inputs.elP, inputs.wboP,
        n_workers=inner_workers, **_rp_cut_kwargs(cfg))
    if BGCP_TIMING:
        print(f"    {inputs.step_name} {'R-P':>12s} cut_sweep: "
              f"{len(pool):>4d} sigs in {time.time()-t0:.1f}s",
              flush=True)
    return run_rp_stage_from_pool(
        inputs, pool, config=cfg, elapsed=time.time() - t0)


def _add_ts_endpoint_tasks(tasks, inputs, key, target_order, target_label,
                           mech_pos, mech, target, config):
    mapping_RP, _br_R, _fm_R, core_R = _mechanism_state(inputs, mech)
    common = {
        'key': key,
        'target_order': target_order,
        'target_label': target_label,
        'mech_id': int(mech['id']),
        'mech_pos': int(mech_pos),
        'elT': target.el,
        'wboT': target.wbo,
        'graph_floor': config.get('graph_floor', TS_ALIGN_GRAPH_FLOOR),
        'iso_tol': config.get('iso_tol', VIEW_ISO_TOL),
        'dwbo_threshold': config.get('dwbo_threshold', DWBO_THRESHOLD),
        'metal_dwbo_threshold': config.get(
            'metal_dwbo_threshold', METAL_DWBO_THRESHOLD),
        'symmetry_wbo_tol': config.get('iso_tol', VIEW_ISO_TOL),
        'n_seeds': config.get('n_seeds', N_SEEDS_PER_RUN),
        'max_core_maps': config.get('max_core_maps', TS_ALIGN_MAX_CORE_MAPS),
    }
    tasks.append({
        **common,
        'endpoint': 'R',
        'elS': inputs.elR,
        'wboS': inputs.wboR,
        'core_S': list(core_R),
    })
    core_P = [int(mapping_RP[r]) for r in core_R if r in mapping_RP]
    tasks.append({
        **common,
        'endpoint': 'P',
        'elS': inputs.elP,
        'wboS': inputs.wboP,
        'core_S': core_P,
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
        target_view = {
            'label': target.label,
            'elements': list(target.el),
            'xyz': np.asarray(target.xyz, float).tolist(),
        }
        target_view.update(_target_metadata(target))
        if target.kind != 'gt':
            ig_index = len(mechanisms[0]['igs']) if mechanisms else 0
            for mech in mechanisms:
                mech['igs'].append(dict(target_view))
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
                print(f"    [warn] TS core alignment hit cap={cfg.get('max_core_maps')} "
                      f"{res['target_label']}:{res['endpoint']} "
                      f"core={res['core_size']}",
                      flush=True)
            print(f"    {inputs.step_name} {res['target_label'] + ':' + res['endpoint']:>12s} "
                  f"core_align: {res['n_pool']:>4d} sigs "
                  f"core={res['core_size']} in {res['elapsed']:.1f}s",
                  flush=True)

    for ctx in sorted(score_contexts,
                      key=lambda x: (x['target_order'], x['mech_pos'])):
        mech = mechanisms[ctx['mech_pos']]
        target = ctx['target']
        mapping_RP, br_R, fm_R, core_R = _mechanism_state(inputs, mech)
        parts = endpoint_by_key.get(ctx['key'], {})
        r_records = parts.get('R', {}).get('branch_records', [])
        p_records_native = parts.get('P', {}).get('branch_records', [])
        p_records_as_r = _product_core_branch_records_to_reactant(
            p_records_native, mapping_RP, core_R)
        merged = _merge_endpoint_core_branch_records(
            core_R, r_records, p_records_as_r)
        if BGCP_TIMING:
            print(f"    {inputs.step_name} {ctx['display_label']:>12s} core_union: "
                  f"R={len(r_records)} P={len(p_records_native)} "
                  f"merged={len(merged)}",
                  flush=True)
        s = best_under_mech_using_branch_pool(
            inputs.elR, inputs.xyzR, inputs.wboR, inputs.wboP,
            target.el, target.xyz, target.wbo, target.freqs, target.modes,
            merged, mapping_RP, br_R, fm_R, core_R,
            score_weights=cfg.get('score', score_config()),
            prefer_endpoint_consensus=cfg.get(
                'prefer_endpoint_consensus', True),
            symmetry_wbo_tol=cfg.get('iso_tol', VIEW_ISO_TOL))
        _attach_target_metadata(s, target)
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
            {k: v for k, v in res.items()
             if k not in {'pool', 'branch_records'}}
            for res in endpoint_results
        ],
        'union_top_labels': sorted(union_top),
    }


def _recompute_ts_top_flags(mechanisms):
    union_top = set()
    for mech in mechanisms:
        ranked = sorted(
            [(i, ig) for i, ig in enumerate(mech.get('igs', []))
             if ig.get('S') is not None],
            key=lambda x: -x[1]['S'])
        top2 = {i for i, _ in ranked[:2]}
        for i, ig in enumerate(mech.get('igs', [])):
            ig['is_top2'] = i in top2
            if i in top2:
                union_top.add(ig['label'])
    for mech in mechanisms:
        for ig in mech.get('igs', []):
            ig['is_union_top'] = ig.get('label') in union_top
    return sorted(union_top)


def merge_ts_stage_chunks(chunks):
    """Merge independent target/candidate TS chunks into one Stage 2 result."""
    chunks = [chunk for chunk in chunks if chunk]
    if not chunks:
        raise RuntimeError("no TS chunks to merge")
    merged_by_id = {}
    endpoint_results = []
    for chunk in chunks:
        endpoint_results.extend(chunk.get('endpoint_results', []))
        for mech in chunk.get('mechanisms', []):
            mech_id = int(mech['id'])
            current = merged_by_id.get(mech_id)
            if current is None:
                current = dict(mech)
                current['igs'] = []
                current['gt'] = None
                merged_by_id[mech_id] = current
            if mech.get('gt') and current.get('gt') is None:
                current['gt'] = mech.get('gt')
            current['igs'].extend(dict(ig) for ig in mech.get('igs', []))
    mechanisms = [
        merged_by_id[key] for key in sorted(merged_by_id)
    ]
    union_top = _recompute_ts_top_flags(mechanisms)
    return {
        'stage': 'ts',
        'step': chunks[0].get('step'),
        'config': chunks[0].get('config', ts_stage_config()),
        'mechanisms': mechanisms,
        'endpoint_results': endpoint_results,
        'union_top_labels': union_top,
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
    graph_floor = float(
        (rp_result.get('config') or {}).get('graph_floor', 0.2))
    for mechanism in mechanisms:
        mapping = _int_mapping(mechanism.get('mapping_RP') or {})
        product_in_R = np.asarray(
            mechanism.get('product_xyz_in_R'), dtype=float)
        if (len(mapping) != len(inputs.elR)
                or product_in_R.shape != np.asarray(inputs.xyzR).shape):
            mechanism['endpoint_interpolation'] = None
            continue
        mechanism['product_xyz_in_R_aligned'] = proper_align_coordinates(
            product_in_R, inputs.xyzR).tolist()
        inverse = {int(p): int(r) for r, p in mapping.items()}
        reactant_bonds = {
            (left, right)
            for left in range(len(inputs.elR))
            for right in range(left + 1, len(inputs.elR))
            if float(inputs.wboR[left, right]) >= graph_floor
        }
        product_bonds = set()
        for left_P in range(len(inputs.elP)):
            for right_P in range(left_P + 1, len(inputs.elP)):
                if float(inputs.wboP[left_P, right_P]) < graph_floor:
                    continue
                if left_P not in inverse or right_P not in inverse:
                    continue
                product_bonds.add(tuple(sorted(
                    (inverse[left_P], inverse[right_P]))))
        bonded_pairs = reactant_bonds | product_bonds
        persistent_bonds = reactant_bonds & product_bonds
        mechanism['endpoint_interpolation'] = internal_coordinate_interpolation(
            inputs.xyzR, product_in_R, inputs.elR,
            bonded_pairs=sorted(bonded_pairs),
            persistent_bonded_pairs=sorted(persistent_bonds),
            reactant_bonded_pairs=sorted(reactant_bonds),
            product_bonded_pairs=sorted(product_bonds), n_frames=101)
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
            'metadata': _endpoint_metadata('R', inputs.energy_R),
        },
        'product': {
            'elements': inputs.elP,
            'coords': np.asarray(inputs.xyzP).tolist(),
            'metadata': _endpoint_metadata('P', inputs.energy_P),
        },
        'mechanisms': mechanisms,
        'default_mech_id': default_id,
        'include_gt': bool(include_gt),
        'score_config': view_score_config,
        'metadata_units': {
            'energy': 'hartree',
            'frequency': 'cm^-1',
        },
    }


def build_eval_slim(view_data):
    slim = {
        'step': view_data['step'],
        'n_atoms': view_data['n_atoms'],
        'n_mechs': len(view_data.get('mechanisms', [])),
        'include_gt': bool(view_data.get('include_gt')),
        'score_config': view_data.get('score_config', score_config()),
        'reactant': {
            'metadata': view_data.get('reactant', {}).get('metadata', {}),
        },
        'product': {
            'metadata': view_data.get('product', {}).get('metadata', {}),
        },
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
                'freq', 'core_map', 'core_sources', 'endpoint_consensus',
                'energy_hartree', 'frequency_summary',
            ]} if mech.get('gt') else None),
            'igs': [
                {k: ig.get(k) for k in [
                    'label', 'S', 'beta', 'wbo_progress',
                    'wbo_progress_factor', 'freq', 'core_map',
                    'core_sources', 'endpoint_consensus', 'is_top2',
                    'energy_hartree', 'frequency_summary',
                ]}
                for ig in mech.get('igs', [])
            ],
        })
    return slim


def write_view_stage(inputs, rp_result, ts_result=None, out_root=None,
                     include_gt=None, return_data=False):
    """Stage 3: write viewer and slim eval artifacts."""
    data = build_view_data(inputs, rp_result, ts_result, include_gt=include_gt)
    slim = build_eval_slim(data)
    paths = pipeline_stage_paths(inputs.step_name, out_root=out_root)
    paths.view_dir.mkdir(parents=True, exist_ok=True)
    paths.view_html.write_text(HTML.format(
        title=f"BGCP &mdash; {inputs.step_name}  "
              f"({len(data['mechanisms'])} mechanisms)",
        data_json=json.dumps(_json_ready(data)),
        three_dmol_js=(VIEWER_STATIC_DIR / "3Dmol-min.js").read_text(),
        jszip_js=(VIEWER_STATIC_DIR / "jszip.min.js").read_text(),
    ))
    paths.eval_slim_json.write_text(json.dumps(_json_ready(slim)))
    result = {
        'stage': 'view',
        'step': inputs.step_name,
        'view_html': str(paths.view_html),
        'eval_slim_json': str(paths.eval_slim_json),
        'slim': slim,
    }
    if return_data:
        result['data'] = data
    return result


HTML = r"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<script>{three_dmol_js}</script>
<script>{jszip_js}</script>
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
.interp-panel{{margin:0 0 18px}}
.interp-panel .vw{{height:360px}}
.interp-controls{{display:flex;align-items:center;gap:10px;width:100%}}
.interp-controls input[type=range]{{flex:1}}
.clash-ok{{color:#16733b}}
.clash-bad{{color:#b00020;font-weight:700}}
</style></head><body>
<h2>{title}</h2>
<div class="topbar">
  <div class="mech-sel" id="mech-sel"></div>
  <label class="index-toggle"><input type="checkbox" id="showAtomIndices"> Atom labels</label>
  <label class="index-toggle"><input type="checkbox" id="showDegeneracy"> Degeneracy</label>
  <label class="index-toggle"><input type="checkbox" id="rOrdered" checked> R-order + spatial alignment</label>
  <button class="download-all" id="downloadAllBtn">Download</button>
</div>
<div class="ref-row" id="ref-row">
  <div class="panel"><div class="ph"><span class="lbl">Reactant</span><span class="rk">static <button class="dl" onclick="downloadR()">XYZ</button></span></div>
    <div class="vw"><div id="vw_R" class="vwbox"></div></div>
    <div class="meta" id="r_meta"></div></div>
  <div class="panel"><div class="ph"><span class="lbl">Product</span><span class="rk" id="prod_label">static <button class="dl" onclick="downloadP()">XYZ</button></span></div>
    <div class="vw"><div id="vw_P" class="vwbox"></div></div>
    <div class="meta" id="p_meta"></div></div>
  <div class="panel" id="gt_panel"><div class="ph"><span class="lbl">Ground-truth TS</span><span class="rk" id="gt_S">S=? <button class="dl" onclick="downloadGT()">XYZ</button></span></div>
    <div class="vw"><div id="vw_GT" class="vwbox"></div></div>
    <div class="meta" id="gt_meta"></div></div>
</div>
<div class="panel interp-panel">
  <div class="ph"><span class="lbl">R&rarr;P internal-coordinate interpolation</span><span class="rk" id="interp_method"></span></div>
  <div class="vw"><div id="vw_interp" class="vwbox"></div></div>
  <div class="interp-controls">
    <button class="dl" id="interpPlayBtn">Play</button>
    <input type="range" id="interpFrame" min="0" max="20" value="0" step="1">
    <span class="rk" id="interp_t">t=0.00</span>
    <label class="index-toggle"><input type="checkbox" id="showClashHighlights"> Highlight clashes</label>
  </div>
  <div class="meta" id="interp_meta"></div>
</div>
<div class="ig-grid" id="grid"></div>
<script>
const DATA = {data_json};
let currentMechId = DATA.default_mech_id;
let showAtomIndices = false;
let showDegeneracy = false;
let rOrdered = true;
let interpPlayTimer = null;
let interpViewer = null;
let interpViewerMechId = null;
let showClashHighlights = false;
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
function downloadP() {{ const mech = findMech(currentMechId); const aligned = mech && (mech.product_xyz_in_R_aligned || mech.product_xyz_in_R); if (rOrdered && aligned) downloadXYZ(safeName(DATA.step)+"_P_R_aligned_mech"+mech.id+".xyz", DATA.reactant.elements, aligned, DATA.step+" P R-ordered and spatially aligned mech "+mech.id); else downloadXYZ(safeName(DATA.step)+"_P_native.xyz", DATA.product.elements, DATA.product.coords, DATA.step+" P native"); }}
function targetFrame(item, mech) {{ if (!item) return {{aligned:false, elements:elements, xyz:null, disp:null, broken:[], formed:[], core:[], label:"missing"}}; const nativeEls = item.elements || elements; const nativeXYZ = item.xyz || null; const nativeDisp = item.picked_disp || null; if (!nativeXYZ) return {{aligned:false, elements:nativeEls, xyz:item.xyz_in_R || null, disp:item.picked_disp_R || null, broken:item.broken_bonds_T || mech.broken_bonds_R || [], formed:item.formed_bonds_T || mech.formed_bonds_R || [], core:item.core_atoms_T || mech.core_atoms || [], label:"native"}}; if (!rOrdered) return {{aligned:false, elements:nativeEls, xyz:nativeXYZ, disp:nativeDisp, broken:item.broken_bonds_T || mech.broken_bonds_R || [], formed:item.formed_bonds_T || mech.formed_bonds_R || [], core:item.core_atoms_T || mech.core_atoms || [], label:"native"}}; const n = nativeXYZ.length; const perm = Array(n).fill(null); const used = new Set(); const cmap = item.core_map || {{}}; for (const rk of Object.keys(cmap)) {{ const r = parseInt(rk), t = parseInt(cmap[rk]); if (Number.isInteger(r) && Number.isInteger(t) && r >= 0 && r < n && t >= 0 && t < n && !used.has(t)) {{ perm[r] = t; used.add(t); }} }} for (let r = 0; r < n; r++) {{ if (perm[r] !== null) continue; const want = elements[r]; let pick = -1; for (let t = 0; t < n; t++) {{ if (!used.has(t) && nativeEls[t] === want) {{ pick = t; break; }} }} if (pick < 0) for (let t = 0; t < n; t++) {{ if (!used.has(t)) {{ pick = t; break; }} }} if (pick >= 0) {{ perm[r] = pick; used.add(pick); }} }} const orderedXYZ = perm.map(t => nativeXYZ[t]); const orderedDisp = nativeDisp ? perm.map(t => nativeDisp[t]) : null; return {{aligned:true, elements:elements, xyz:orderedXYZ, disp:orderedDisp, broken:mech.broken_bonds_R || [], formed:mech.formed_bonds_R || [], core:mech.core_atoms || [], label:"R-ordered"}}; }}
function downloadTargetXYZ(name, item, mech, comment) {{ const frame = targetFrame(item, mech); downloadXYZ(name, frame.elements, frame.xyz, comment+" "+frame.label); }}
function targetSuffix() {{ return rOrdered ? "_R_ordered" : "_native"; }}
function downloadGT() {{ const mech = findMech(currentMechId); if (mech.gt) downloadTargetXYZ(safeName(DATA.step)+"_GT_mech"+mech.id+targetSuffix()+".xyz", mech.gt, mech, DATA.step+" GT mech "+mech.id); }}
function downloadIG(ig) {{ const mech = findMech(currentMechId); downloadTargetXYZ(safeName(DATA.step)+"_"+safeName(ig.label)+targetSuffix()+".xyz", ig, mech, DATA.step+" "+ig.label); }}
function isFiniteNumber(v) {{ return typeof v === "number" && isFinite(v); }}
function formatEnergy(v) {{ return isFiniteNumber(v) ? v.toFixed(8)+" Eh" : "n/a"; }}
function formatFreq(v) {{ return isFiniteNumber(v) ? v.toFixed(0)+" cm^-1" : "n/a"; }}
function structureMeta(item) {{ const meta = (item && item.metadata) ? item.metadata : item; if (!meta) return "<b>E</b>=n/a"; return "<b>E</b>="+formatEnergy(meta.energy_hartree); }}
function freqSummaryText(item) {{ if (!item || !item.frequency_summary) return ""; const fs = item.frequency_summary; const imag = fs.imaginary_cm1 || []; if (imag.length) return imag.map(formatFreq).join(", "); if (isFiniteNumber(fs.lowest_cm1)) return "min "+formatFreq(fs.lowest_cm1); return ""; }}
function targetAnalysisRecord(item) {{ if (!item) return null; return {{label:item.label || null, energy_hartree:item.energy_hartree ?? null, energy_units:item.energy_units || "hartree", frequency_units:item.frequency_units || "cm^-1", frequencies_cm1:item.frequencies_cm1 || null, frequency_summary:item.frequency_summary || null, picked_mode:{{index:item.k ?? null, frequency_cm1:item.freq ?? null}}}}; }}
function scoreRecord(item) {{ if (!item || item.S === undefined || item.S === null) return null; return {{S:item.S, energy_hartree:item.energy_hartree ?? null, energy_units:item.energy_units || "hartree", frequency_units:item.frequency_units || "cm^-1", frequencies_cm1:item.frequencies_cm1 || null, frequency_summary:item.frequency_summary || null, decomposition:{{beta:item.beta, wbo_progress:item.wbo_progress, wbo_progress_factor:item.wbo_progress_factor, freq:item.freq, mode_index:item.k}}, event_terms:item.event_terms || [], core_map:item.core_map, core_sources:item.core_sources, core_pool_dedup_count:item.core_pool_dedup_count, endpoint_consensus:item.endpoint_consensus}}; }}
function mechanismRecord(mech) {{ return {{id:mech.id, label:mech.label, cut:mech.cut, dedup_count:mech.dedup_count || 1, dedup_source_ids:mech.dedup_source_ids || [mech.id], dedup_cuts:mech.dedup_cuts || [mech.cut], broken_bonds_R:mech.broken_bonds_R, formed_bonds_R:mech.formed_bonds_R, formed_bonds_P:mech.formed_bonds_P || [], core_atoms_R:mech.core_atoms || [], branch_symmetry:mech.branch_symmetry || null, index_chirality:mech.index_chirality || null, gt:scoreRecord(mech.gt), igs:(mech.igs || []).map(ig => ({{label:ig.label, energy_hartree:ig.energy_hartree ?? null, frequency_summary:ig.frequency_summary || null, is_top2:!!ig.is_top2, is_union_top:!!ig.is_union_top, score:scoreRecord(ig)}}))}}; }}
function buildEnergyFrequencySummary() {{ return {{step:DATA.step, units:{{energy:"hartree", frequency:"cm^-1"}}, reactant:(DATA.reactant && DATA.reactant.metadata) || null, product:(DATA.product && DATA.product.metadata) || null, mechanisms:(DATA.mechanisms || []).map(mech => ({{id:mech.id, label:mech.label, gt:targetAnalysisRecord(mech.gt), igs:(mech.igs || []).map(targetAnalysisRecord)}}))}}; }}
function buildArchiveManifest() {{ return {{step:DATA.step, n_atoms:DATA.n_atoms, include_gt:!!DATA.include_gt, default_mech_id:DATA.default_mech_id, score_formula:"S = beta * wbo_progress^WBO_PROGRESS_POWER", score_config:DATA.score_config || null, endpoint_metadata:{{reactant:(DATA.reactant && DATA.reactant.metadata) || null, product:(DATA.product && DATA.product.metadata) || null}}, mechanisms:(DATA.mechanisms || []).map(mechanismRecord), files:{{reactant:"R.xyz", product:"P.xyz", gt:"GT/GT.xyz if available", ig:"IG/<label>.xyz", energy_frequency_summary:"energies_frequencies.json", per_mechanism:"mechanisms/mechanism_<id>.json", full_viewer_data:"viewer_data.json", view_html:"view.html"}}}}; }}
function scoreMeta(item) {{ if (!item) return "(no data)"; const parts = []; if (item.energy_hartree !== undefined && item.energy_hartree !== null) parts.push("<b>E</b>="+formatEnergy(item.energy_hartree)); if (item.beta !== undefined) parts.push("<b>&beta;</b>="+item.beta.toFixed(3)); if (item.wbo_progress !== undefined) parts.push("<b>wbo</b>="+item.wbo_progress.toFixed(3)); if (item.freq !== undefined && item.freq !== null) parts.push("<b>mode</b>="+formatFreq(item.freq)); const fs = freqSummaryText(item); if (fs) parts.push("<b>freqs</b>="+fs); return parts.length ? parts.join(" ") : "(no data)"; }}
function viewerHtmlForArchive() {{ return "<!doctype html>\n" + document.documentElement.outerHTML; }}
async function downloadAll() {{ if (typeof JSZip === "undefined") {{ alert("Download library is not loaded"); return; }} const root = safeName(DATA.step); const zip = new JSZip(); zip.file(root+"/R.xyz", xyzText(DATA.reactant.elements, DATA.reactant.coords, DATA.step+" R")); zip.file(root+"/P.xyz", xyzText(DATA.product.elements, DATA.product.coords, DATA.step+" P")); const firstGTMech = (DATA.mechanisms || []).find(m => m.gt && (m.gt.xyz || m.gt.xyz_in_R)); if (firstGTMech) {{ const frame = targetFrame(firstGTMech.gt, firstGTMech); zip.file(root+"/GT/GT"+targetSuffix()+".xyz", xyzText(frame.elements, frame.xyz, DATA.step+" GT "+frame.label)); }} const seenIG = new Set(); for (const mech of DATA.mechanisms || []) {{ for (const ig of mech.igs || []) {{ if (seenIG.has(ig.label)) continue; const frame = targetFrame(ig, mech); if (!frame.xyz) continue; seenIG.add(ig.label); zip.file(root+"/IG/"+safeName(ig.label)+targetSuffix()+".xyz", xyzText(frame.elements, frame.xyz, DATA.step+" "+ig.label+" "+frame.label)); }} }} const manifest = buildArchiveManifest(); zip.file(root+"/mechanism.json", JSON.stringify(manifest, null, 2)); zip.file(root+"/energies_frequencies.json", JSON.stringify(buildEnergyFrequencySummary(), null, 2)); for (const mech of DATA.mechanisms || []) {{ zip.file(root+"/mechanisms/mechanism_"+String(mech.id).padStart(3,"0")+".json", JSON.stringify(mechanismRecord(mech), null, 2)); }} zip.file(root+"/viewer_data.json", JSON.stringify(DATA, null, 2)); zip.file(root+"/view.html", viewerHtmlForArchive()); const blob = await zip.generateAsync({{type:"blob"}}); downloadBlob(root+".zip", blob); }}
const animTimers = {{}};
function stopAnim(d) {{ if (animTimers[d]) {{ clearInterval(animTimers[d]); delete animTimers[d]; }} }}
function clearLabels(v) {{ if (v.removeAllLabels) v.removeAllLabels(); }}
const DEG_PALETTE = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#17becf','#bcbd22','#3182bd','#e6550d','#31a354','#756bb1','#fdae6b','#74c476','#fd8d3c'];
function branchBlocks(mech) {{ return (mech && mech.branch_symmetry && mech.branch_symmetry.blocks) ? mech.branch_symmetry.blocks : []; }}
function branchColorGroups(mech) {{ return (mech && mech.branch_symmetry && mech.branch_symmetry.color_groups && mech.branch_symmetry.color_groups.length) ? mech.branch_symmetry.color_groups : branchBlocks(mech); }}
function branchRDegMap(mech) {{ const out = {{}}; branchColorGroups(mech).forEach((b, idx) => {{ for (const r of b.r_atoms || []) out[r] = idx; }}); return out; }}
function branchProductDegMap(mech, pAligned) {{ const out = {{}}; const mapping = mech.mapping_RP || {{}}; const pToR = {{}}; for (const r of Object.keys(mapping)) pToR[String(mapping[r])] = r; branchColorGroups(mech).forEach((b, idx) => {{ for (const p of b.p_atoms || []) {{ if (pAligned) {{ if (pToR[String(p)] !== undefined) out[pToR[String(p)]] = idx; }} else {{ out[p] = idx; }} }} }}); return out; }}
function rDegMap(mech) {{ return showDegeneracy ? branchRDegMap(mech) : null; }}
function mappedProductDegMap(mech, pAligned) {{ return showDegeneracy ? branchProductDegMap(mech, pAligned) : null; }}
function targetDegMap(frame, mech) {{ return frame && frame.aligned ? rDegMap(mech) : null; }}
function applyDegeneracyStyles(v, degMap) {{ if (!degMap) return; for (const idx of Object.keys(degMap)) {{ const group = parseInt(degMap[idx]); if (!Number.isFinite(group)) continue; const color = DEG_PALETTE[Math.abs(group) % DEG_PALETTE.length]; v.setStyle({{serial: parseInt(idx)}}, {{stick:{{radius:0.13,color:color}}, sphere:{{scale:0.34,color:color}}}}); }} }}
function reactantAtomLabels() {{ return DATA.reactant.elements.map((_, i) => "R"+i); }}
function productAtomLabels(mech, pAligned) {{ const els = pAligned ? DATA.reactant.elements : DATA.product.elements; return els.map((_, i) => String(i)); }}
function targetAtomLabels(frame) {{ if (frame && frame.aligned) return DATA.reactant.elements.map((_, i) => "R"+i); const els = (frame && frame.elements) || []; return els.map((_, i) => String(i)); }}
function addAtomLabels(v, els, xyz, labels=null) {{ if (!showAtomIndices || !els || !xyz) return; for (let i=0;i<xyz.length;i++) {{ const text = labels && labels[i] !== undefined ? String(labels[i]) : String(i); v.addLabel(text, {{position:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, fontSize:9, fontColor:'black', backgroundColor:'white', backgroundOpacity:0.72, borderColor:'#666', borderThickness:0.5, inFront:true}}); }} }}
function drawBonds(v, xyz, pairs, color) {{ for (const [i,j] of pairs) {{ if (i>=xyz.length||j>=xyz.length) continue; v.addCylinder({{start:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, end:{{x:xyz[j][0],y:xyz[j][1],z:xyz[j][2]}}, color:color, radius:0.16, dashed:true}}); }} }}
function drawArrows(v, xyz, disp, core) {{ for (const i of core) {{ if (!disp||!disp[i]) continue; const d = disp[i]; const len = Math.hypot(d[0],d[1],d[2]); if (len<0.05) continue; v.addArrow({{start:{{x:xyz[i][0],y:xyz[i][1],z:xyz[i][2]}}, end:{{x:xyz[i][0]+d[0]*1.5,y:xyz[i][1]+d[1]*1.5,z:xyz[i][2]+d[2]*1.5}}, color:'#0066cc', radius:0.07}}); }} }}
function makeStatic(divId, els, xyz, broken, formed, degMap=null, labels=null) {{ stopAnim(divId); document.getElementById(divId).innerHTML=""; const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}}); v.addModel(buildBody(els, xyz), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); applyDegeneracyStyles(v, degMap); drawBonds(v, xyz, broken, 'red'); drawBonds(v, xyz, formed, 'green'); addAtomLabels(v, els, xyz, labels); v.zoomTo(); v.render(); return v; }}
function makeAnimated(divId, els, xyz, disp, broken, formed, core, degMap=null, labels=null) {{ stopAnim(divId); document.getElementById(divId).innerHTML=""; const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}}); v.addModel(buildBody(els, xyz), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); applyDegeneracyStyles(v, degMap); drawBonds(v, xyz, broken, 'red'); drawBonds(v, xyz, formed, 'green'); drawArrows(v, xyz, disp, core); addAtomLabels(v, els, xyz, labels); v.zoomTo(); v.render(); let t=0; const period=30, amp=0.6; animTimers[divId] = setInterval(()=>{{ t=(t+1)%period; const scale = amp*Math.sin(2*Math.PI*t/period); const cur = xyzAt(xyz, disp, scale); v.removeAllModels(); v.removeAllShapes(); clearLabels(v); v.addModel(buildBodyAt(els, xyz, disp, scale), 'xyz'); v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}}); applyDegeneracyStyles(v, degMap); drawBonds(v, cur, broken, 'red'); drawBonds(v, cur, formed, 'green'); drawArrows(v, cur, disp, core); addAtomLabels(v, els, cur, labels); v.render(); }}, 60); return v; }}
function stopInterpolation() {{ if (interpPlayTimer) {{ clearInterval(interpPlayTimer); interpPlayTimer=null; }} document.getElementById('interpPlayBtn').textContent='Play'; }}
function renderInterpolationFrame() {{ const mech=findMech(currentMechId); const path=mech && mech.endpoint_interpolation; const div=document.getElementById('vw_interp'); if (!path || !path.frames || !path.frames.length) {{ div.innerHTML=''; interpViewer=null; interpViewerMechId=null; document.getElementById('interp_meta').textContent='Interpolation unavailable'; return; }} const slider=document.getElementById('interpFrame'); slider.max=String(path.frames.length-1); const index=Math.max(0,Math.min(path.frames.length-1,parseInt(slider.value)||0)); const frame=path.frames[index]; const sameViewer=!!(interpViewer && interpViewerMechId===mech.id); const savedView=(sameViewer && interpViewer.getView) ? interpViewer.getView() : null; if (!sameViewer) {{ div.innerHTML=''; interpViewer=$3Dmol.createViewer('vw_interp',{{backgroundColor:'white'}}); interpViewerMechId=mech.id; }} else {{ interpViewer.removeAllModels(); interpViewer.removeAllShapes(); clearLabels(interpViewer); }} const v=interpViewer; v.addModel(buildBody(DATA.reactant.elements,frame.coords),'xyz'); v.setStyle({{}},{{stick:{{radius:0.10}},sphere:{{scale:0.20}}}}); applyDegeneracyStyles(v,rDegMap(mech)); drawBonds(v,frame.coords,mech.broken_bonds_R,'red'); drawBonds(v,frame.coords,mech.formed_bonds_R,'green'); addAtomLabels(v,DATA.reactant.elements,frame.coords,reactantAtomLabels()); if (showClashHighlights) {{ const clashAtoms=new Set(); for (const pair of (frame.clashes.pairs||[])) for (const atom of pair.atoms) clashAtoms.add(atom); for (const atom of clashAtoms) v.setStyle({{serial:atom}},{{stick:{{radius:0.14,color:'#d000ff'}},sphere:{{scale:0.42,color:'#d000ff'}}}}); }} if (savedView && v.setView) v.setView(savedView); else v.zoomTo(); v.render(); document.getElementById('interp_t').textContent='t='+frame.t.toFixed(2); document.getElementById('interp_method').textContent=path.method; const ratio=frame.clashes.minimum_radius_ratio; const pairs=(frame.clashes.pairs||[]).map(p=>'R'+p.atoms[0]+'–R'+p.atoms[1]).join(', '); const residuals=frame.constraint_residuals||{{}}; const geometry=' <b>persistent bond error</b>='+((100*(residuals.max_persistent_bond_relative_error||0)).toFixed(2))+'% <b>angle error</b>='+((residuals.max_angle_error_degrees||0).toFixed(2))+'° <b>torsion error</b>='+((residuals.max_torsion_error_degrees||0).toFixed(2))+'°'; document.getElementById('interp_meta').innerHTML='<span class="'+(frame.clashes.count?'clash-bad':'clash-ok')+'"><b>severe overlaps</b>='+frame.clashes.count+'</span> <b>min clearance</b>='+(ratio===null?'n/a':ratio.toFixed(3)+' × covalent radii')+geometry+(pairs?' <b>pairs</b>='+pairs:''); }}
function toggleInterpolation() {{ if (interpPlayTimer) {{ stopInterpolation(); return; }} document.getElementById('interpPlayBtn').textContent='Pause'; interpPlayTimer=setInterval(()=>{{ const slider=document.getElementById('interpFrame'); let next=(parseInt(slider.value)||0)+1; if (next>parseInt(slider.max)) next=0; slider.value=String(next); renderInterpolationFrame(); }},180); }}
function render() {{ const mech = findMech(currentMechId); document.getElementById('r_meta').innerHTML = structureMeta(DATA.reactant); document.getElementById('p_meta').innerHTML = structureMeta(DATA.product); document.querySelectorAll('.mech-sel button[data-id]').forEach(b => {{ b.classList.toggle('active', parseInt(b.dataset.id)===currentMechId); }}); makeStatic('vw_R', DATA.reactant.elements, DATA.reactant.coords, mech.broken_bonds_R, [], rDegMap(mech), reactantAtomLabels()); const alignedP = mech.product_xyz_in_R_aligned || mech.product_xyz_in_R; const pAligned = !!(rOrdered && alignedP); const pEls = pAligned ? DATA.reactant.elements : DATA.product.elements; const pXYZ = pAligned ? alignedP : DATA.product.coords; const pFormed = pAligned ? mech.formed_bonds_R : (mech.formed_bonds_P || []); makeStatic('vw_P', pEls, pXYZ, [], pFormed, mappedProductDegMap(mech, pAligned), productAtomLabels(mech, pAligned)); document.getElementById('prod_label').innerHTML = (pAligned ? "R-ordered + spatially aligned" : "native P")+" (mech #"+mech.id+") <button class='dl' onclick='downloadP()'>XYZ</button>"; renderInterpolationFrame(); const gtFrame = mech.gt ? targetFrame(mech.gt, mech) : null; const showGT = !!(gtFrame && gtFrame.disp); document.getElementById('ref-row').classList.toggle('no-gt', !showGT); document.getElementById('gt_panel').style.display = showGT ? "" : "none"; if (showGT) {{ makeAnimated('vw_GT', gtFrame.elements, gtFrame.xyz, gtFrame.disp, gtFrame.broken, gtFrame.formed, gtFrame.core, targetDegMap(gtFrame, mech), targetAtomLabels(gtFrame)); document.getElementById('gt_S').innerHTML = "S = "+mech.gt.S.toFixed(3)+" <button class='dl' onclick='downloadGT()'>XYZ</button>"; document.getElementById('gt_meta').innerHTML = scoreMeta(mech.gt)+" <b>frame</b>="+gtFrame.label; }} else {{ stopAnim('vw_GT'); document.getElementById('vw_GT').innerHTML = ""; }} const grid = document.getElementById('grid'); grid.innerHTML = ""; const igs = [...mech.igs].sort((a,b) => (b.S||0) - (a.S||0)); igs.forEach((ig, idx) => {{ const frame = targetFrame(ig, mech); const div = document.createElement('div'); let cls = 'panel'; if (ig.is_top2) cls += ' top2'; if (ig.is_union_top && !ig.is_top2) cls += ' union'; div.className = cls; const sStr = ig.S !== undefined ? "S = "+ig.S.toFixed(3) : "no score"; const tag = ig.is_top2 ? '<span style="background:#d4af37;color:white;padding:1px 6px;border-radius:3px;font-size:11px;margin-left:4px">TOP2</span>' : (ig.is_union_top ? '<span style="background:#ff9;color:#660;padding:1px 6px;border-radius:3px;font-size:11px;margin-left:4px">union</span>' : ''); const dl = frame.xyz ? '<button class="dl">XYZ</button> ' : ''; div.innerHTML = '<div class="ph"><span class="lbl">'+ig.label+tag+'</span><span class="rk">'+dl+sStr+'</span></div><div class="vw"><div id="vw_ig'+idx+'" class="vwbox"></div></div><div class="meta">'+scoreMeta(ig)+" <b>frame</b>="+frame.label+"</div>"; grid.appendChild(div); const btn = div.querySelector('button.dl'); if (btn) btn.onclick = () => downloadIG(ig); if (frame.disp) makeAnimated("vw_ig"+idx, frame.elements, frame.xyz, frame.disp, frame.broken, frame.formed, frame.core, targetDegMap(frame, mech), targetAtomLabels(frame)); else if (frame.xyz) makeStatic("vw_ig"+idx, frame.elements, frame.xyz, frame.broken, frame.formed, targetDegMap(frame, mech), targetAtomLabels(frame)); }}); }}
const ms = document.getElementById('mech-sel'); ms.innerHTML = "<span style='font-size:13px;margin-right:8px;color:#444'>Mechanism:</span>"; document.getElementById('downloadAllBtn').onclick = downloadAll; document.getElementById('showAtomIndices').onchange = (e) => {{ showAtomIndices = !!e.target.checked; render(); }}; document.getElementById('showDegeneracy').onchange = (e) => {{ showDegeneracy = !!e.target.checked; render(); }}; document.getElementById('rOrdered').onchange = (e) => {{ rOrdered = !!e.target.checked; render(); }}; DATA.mechanisms.forEach(m => {{ const b = document.createElement('button'); b.dataset.id = m.id; b.textContent = m.label + (m.gt ? "  GT S=" + m.gt.S.toFixed(3) : ""); if ((m.dedup_count||1) > 1) b.title = "Analytical branch encounters: "+m.dedup_count+"; source mechanisms: "+m.dedup_source_ids.join(", ")+"; cuts: "+m.dedup_cuts.join(", "); b.onclick = () => {{ currentMechId = m.id; render(); }}; ms.appendChild(b); }});
document.getElementById('interpFrame').oninput = renderInterpolationFrame;
document.getElementById('interpPlayBtn').onclick = toggleInterpolation;
document.getElementById('showClashHighlights').onchange = (e) => {{ showClashHighlights=!!e.target.checked; renderInterpolationFrame(); }};
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


def _weighted_graph_from_json(path):
    path = Path(path)
    data = json.loads(path.read_text())
    weights = data.get("weights", data.get("wbo"))
    if weights is None:
        raise ValueError(f"{path} must contain 'weights' or 'wbo'")
    if "nodes" not in data:
        raise ValueError(f"{path} must contain 'nodes'")
    return WeightedGraph(
        nodes=data["nodes"],
        weights=np.asarray(weights, float),
        weight_name=data.get("weight_name", "wbo"),
        coords=(
            None if data.get("coords") is None
            else np.asarray(data.get("coords"), float)
        ),
        metadata=dict(data.get("metadata") or {}),
    )


def _parse_subgraph_anchor(raw):
    text = str(raw).strip()
    for sep in (":", "="):
        if sep in text:
            left, right = text.split(sep, 1)
            return int(left), int(right)
    raise ValueError(
        f"anchor must be formatted as query:target or query=target: {raw}")


def _subgraph_anchor_map_from_cli(values):
    out = {}
    for raw in values or ():
        q, t = _parse_subgraph_anchor(raw)
        if q in out and out[q] != t:
            raise ValueError(f"conflicting anchors for query node {q}")
        if t in out.values():
            raise ValueError(f"target node appears in multiple anchors: {t}")
        out[q] = t
    return out


def _subgraph_node_policy_from_cli(fields):
    fields = list(fields or [])
    if not fields:
        return None
    return fields[0] if len(fields) == 1 else tuple(fields)


def _subgraph_match_record(match):
    return {
        "mapping": {str(k): int(v) for k, v in sorted(match.mapping.items())},
        "query_nodes": [int(v) for v in match.query_nodes],
        "target_nodes": [int(v) for v in match.target_nodes],
        "deferred_edges": [list(map(int, e)) for e in match.deferred_edges],
        "symmetry_fragments": list(match.symmetry_fragments),
    }


def run_subgraph_cli(query_json, target_json, *, node_policy_fields=None,
                     anchor_values=None, graph_floor=0.2, iso_tol=1.0,
                     symmetry_wbo_tol=0.2, orbit_dedup=True,
                     seed_order=None):
    query = _weighted_graph_from_json(query_json)
    target = _weighted_graph_from_json(target_json)
    anchor_map = _subgraph_anchor_map_from_cli(anchor_values)
    node_policy = _subgraph_node_policy_from_cli(node_policy_fields)
    matches = match_weighted_subgraph(
        query,
        target,
        node_policy=node_policy,
        anchor_map=anchor_map,
        graph_floor=graph_floor,
        iso_tol=iso_tol,
        symmetry_wbo_tol=iso_tol,
        seed_order=seed_order,
        orbit_dedup=orbit_dedup,
    )
    return {
        "mode": "subgraph",
        "query_json": str(Path(query_json)),
        "target_json": str(Path(target_json)),
        "node_policy": list(node_policy) if isinstance(node_policy, tuple)
        else node_policy,
        "anchor_map": {str(k): int(v) for k, v in sorted(anchor_map.items())},
        "graph_floor": float(graph_floor),
        "iso_tol": float(iso_tol),
        "symmetry_wbo_tol": float(iso_tol),
        "orbit_dedup": bool(orbit_dedup),
        "n_matches": len(matches),
        "matches": [_subgraph_match_record(match) for match in matches],
    }


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
    workers = _resolve_xtb_workers(inner_workers)
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


def process_smiles_stage(name, reactant_smiles, product_smiles, *,
                         workdir=None, stage='full', inner_workers=0,
                         save_alignment_files=False, rp_config=None,
                         resume_rp=False, sanitize=True,
                         component_spacing=3.0,
                         expand_hydrogens=True):
    """Run R-P stages for SMILES/CXSMILES endpoints.

    This uses formal bond orders from the written SMILES graph.  It does not
    run xtb.  By default, RDKit hydrogen counts are materialized as H atoms.
    """
    workdir = Path(workdir or (PROJECT / "out" / "smiles_work" / name))
    inputs = smiles_inputs_from_strings(
        reactant_smiles, product_smiles,
        name=name, workdir=workdir, sanitize=sanitize,
        component_spacing=component_spacing,
        expand_hydrogens=expand_hydrogens)
    paths = pipeline_stage_paths(name)

    if stage in {'ts', 'post-rp'}:
        raise ValueError(
            "SMILES mode is an R-P formal-bond-order workflow; use XYZ/WBO "
            "inputs for TS/IG validation")

    if stage == 'view':
        if not paths.rp_json.exists():
            raise RuntimeError(
                f"missing Stage 1 artifact: {paths.rp_json}; run --stage rp first")
        rp_result = read_stage_json(paths.rp_json)
        view_result = write_view_stage(
            inputs, rp_result, ts_result=None, include_gt=False)
        return {
            'step': name,
            'rp': rp_result,
            'view': view_result,
            'slim': view_result['slim'],
        }

    if stage not in {'rp', 'full'}:
        raise ValueError(f"unknown stage: {stage}")

    if resume_rp:
        if not paths.rp_json.exists():
            raise RuntimeError(
                f"missing Stage 1 artifact: {paths.rp_json}; run --stage rp first")
        rp_result = read_stage_json(paths.rp_json)
    else:
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
        'ts': None,
        'view': view_result,
        'alignment_files': alignment_files,
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
    global XTB_WORKERS
    global XTB_CHARGE, XTB_MULTIPLICITY
    global VIEW_ISO_TOL, DWBO_THRESHOLD, METAL_DWBO_THRESHOLD, SYMMETRY_WBO_TOL
    global INDEX_CHIRALITY
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
    ap.add_argument("--symmetry-wbo-tol", type=float, default=None,
                    help=argparse.SUPPRESS)
    ap.add_argument("--index-chiral-mode", "--index-chirality",
                    dest="index_chirality",
                    choices=("off", "preserve"),
                    default=INDEX_CHIRALITY,
                    help="Post-process the selected final automorphism to "
                         "preserve endpoint index chirality. Default off.")
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
    ap.add_argument("--xtb-workers", default=XTB_WORKERS,
                    help="Concurrent xtb target-cache jobs. Default 'auto' "
                         "uses available CPUs divided by xtb OMP threads, "
                         "capped by inner workers.")
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
    ap.add_argument("--reactant-smiles", default=None,
                    help="Direct Stage 1 reactant endpoint SMILES/CXSMILES. "
                         "Uses formal bond orders instead of xtb WBO.")
    ap.add_argument("--product-smiles", default=None,
                    help="Direct Stage 1 product endpoint SMILES/CXSMILES. "
                         "Use with --reactant-smiles instead of --steps.")
    ap.add_argument("--smiles-component-spacing", type=float, default=3.0,
                    help="Display-only spacing between disconnected SMILES "
                         "components in generated 2D XYZ coordinates.")
    smiles_h = ap.add_mutually_exclusive_group()
    smiles_h.add_argument("--smiles-expand-hydrogens",
                          dest="smiles_expand_hydrogens",
                          action="store_true", default=True,
                          help="Materialize SMILES atom hydrogen counts as "
                               "explicit H atoms before formal-WBO graph "
                               "construction. This is the default.")
    smiles_h.add_argument("--smiles-preserve-explicit-only",
                          dest="smiles_expand_hydrogens",
                          action="store_false",
                          help="Keep only atoms explicitly present in the "
                               "parsed SMILES graph; atom hydrogen counts "
                               "remain implicit.")
    ap.add_argument("--workdir", default=None,
                    help="Direct XYZ mode cache work directory. Holds "
                         "endpoint, TS single-point, and TS Hessian caches. "
                         "In SMILES mode, holds source formal-WBO files.")
    ap.add_argument("--target-xyz", action="append", default=None,
                    help="Direct Stage 2 TS/IG/GT XYZ. Can be repeated.")
    ap.add_argument("--target-label", action="append", default=None,
                    help="Label for each --target-xyz. Defaults to file stem.")
    ap.add_argument("--target-kind", action="append",
                    choices=("ig", "gt"), default=None,
                    help="Kind for each --target-xyz. Defaults to ig.")
    ap.add_argument("--subgraph-query-json", default=None,
                    help="WeightedGraph JSON query for standalone subgraph "
                         "matching. Use with --subgraph-target-json.")
    ap.add_argument("--subgraph-target-json", default=None,
                    help="WeightedGraph JSON target for standalone subgraph "
                         "matching.")
    ap.add_argument("--subgraph-node-policy", action="append", default=None,
                    help="Node attribute/feature field for subgraph node "
                         "compatibility. Repeat for a multi-field key. "
                         "Default is same element.")
    ap.add_argument("--anchor", action="append", default=None,
                    help="Exact R:P anchor constraint for direct R-P AAM. "
                         "Can be repeated. Same format as --subgraph-anchor.")
    ap.add_argument("--subgraph-anchor", action="append", default=None,
                    help="Exact query:target anchor constraint for standalone "
                         "subgraph matching. In direct R-P mode this is also "
                         "accepted as an R:P AAM anchor. Can be repeated.")
    ap.add_argument("--subgraph-output", default=None,
                    help="Output JSON for standalone subgraph matching. "
                         "Defaults to BGCP_EVAL_JSON.")
    ap.add_argument("--subgraph-no-orbit-dedup", action="store_true",
                    help="Disable orbit dedupe for standalone subgraph "
                         "matching.")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--steps", nargs="+", default=None)
    args = ap.parse_args()

    XTB_CACHE_MODE = _normal_xtb_mode(args.xtb_mode)
    INCLUDE_GT = bool(args.include_gt)
    VIEW_ISO_TOL = float(args.iso_tol)
    DWBO_THRESHOLD = float(args.dwbo_threshold)
    METAL_DWBO_THRESHOLD = float(args.metal_dwbo_threshold)
    if (args.symmetry_wbo_tol is not None
            and not np.isclose(args.symmetry_wbo_tol, VIEW_ISO_TOL)):
        ap.error("--symmetry-wbo-tol must equal --iso-tol; use --iso-tol")
    SYMMETRY_WBO_TOL = VIEW_ISO_TOL
    INDEX_CHIRALITY = str(args.index_chirality)
    EVENT_WEIGHT_POWER = float(args.event_weight_power)
    WBO_PROGRESS_POWER = float(args.wbo_progress_power)
    STAGE_ROOT = Path(args.stage_root)
    ALIGNMENT_OUT_ROOT = Path(args.alignment_out_root)
    XTB_MAX_THREADS = max(1, int(args.xtb_max_threads))
    XTB_WORKERS = args.xtb_workers
    XTB_CHARGE = int(args.charge)
    XTB_MULTIPLICITY = _normal_multiplicity(args.multiplicity)
    XTB_OMP_THREADS = _resolve_xtb_threads(
        args.xtb_omp_threads, XTB_MAX_THREADS)
    os.environ["BGCP_XTB_MODE"] = XTB_CACHE_MODE
    os.environ["BGCP_INCLUDE_GT"] = "1" if INCLUDE_GT else "0"
    os.environ["BGCP_XTB_OMP_THREADS"] = str(XTB_OMP_THREADS)
    os.environ["BGCP_XTB_MAX_THREADS"] = str(XTB_MAX_THREADS)
    os.environ["BGCP_XTB_WORKERS"] = str(XTB_WORKERS)
    os.environ["BGCP_CHARGE"] = str(XTB_CHARGE)
    os.environ["BGCP_MULTIPLICITY"] = str(XTB_MULTIPLICITY)
    os.environ["BGCP_ISO_TOL"] = str(VIEW_ISO_TOL)
    os.environ["BGCP_DWBO_THRESHOLD"] = str(DWBO_THRESHOLD)
    os.environ["BGCP_METAL_DWBO_THRESHOLD"] = str(METAL_DWBO_THRESHOLD)
    os.environ["BGCP_SYMMETRY_WBO_TOL"] = str(SYMMETRY_WBO_TOL)
    os.environ["BGCP_INDEX_CHIRALITY"] = INDEX_CHIRALITY
    os.environ["BGCP_EVENT_WEIGHT_POWER"] = str(EVENT_WEIGHT_POWER)
    os.environ["BGCP_WBO_PROGRESS_POWER"] = str(WBO_PROGRESS_POWER)
    os.environ["BGCP_STAGE_ROOT"] = str(STAGE_ROOT)
    os.environ["BGCP_ALIGNMENT_OUT_ROOT"] = str(ALIGNMENT_OUT_ROOT)
    os.environ["BGCP_SAVE_ALIGNMENT_FILES"] = (
        "1" if args.save_alignment_files else "0")
    os.environ["BGCP_RESUME_RP"] = "1" if args.resume_rp else "0"

    subgraph_mode = bool(args.subgraph_query_json or args.subgraph_target_json)
    if subgraph_mode:
        if (args.steps or args.reactant_xyz or args.product_xyz
                or args.reactant_smiles or args.product_smiles):
            ap.error("subgraph mode cannot be combined with --steps or "
                     "direct endpoint inputs")
        if not (args.subgraph_query_json and args.subgraph_target_json):
            ap.error("--subgraph-query-json and --subgraph-target-json "
                     "must be provided together")
        try:
            result = run_subgraph_cli(
                args.subgraph_query_json,
                args.subgraph_target_json,
                node_policy_fields=args.subgraph_node_policy,
                anchor_values=list(args.subgraph_anchor or []) + list(
                    args.anchor or []),
                graph_floor=0.2,
                iso_tol=VIEW_ISO_TOL,
                symmetry_wbo_tol=SYMMETRY_WBO_TOL,
                orbit_dedup=not args.subgraph_no_orbit_dedup,
            )
        except Exception as e:
            result = {
                "mode": "subgraph",
                "query_json": args.subgraph_query_json,
                "target_json": args.subgraph_target_json,
                "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc(),
            }
        out_path = Path(args.subgraph_output or EVAL_JSON)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(result, indent=2))
        if result.get("error"):
            print(f"subgraph: ERROR: {result['error']}")
        else:
            print(f"subgraph: matches={result['n_matches']} "
                  f"out={out_path}")
        return

    smiles_mode = bool(args.reactant_smiles or args.product_smiles)
    if smiles_mode:
        if args.steps or args.reactant_xyz or args.product_xyz:
            ap.error("use only one input mode: --steps, XYZ, or SMILES")
        if args.target_xyz:
            ap.error("SMILES mode is R-P only; use XYZ/WBO inputs for TS targets")
        if not (args.reactant_smiles and args.product_smiles):
            ap.error("--reactant-smiles and --product-smiles must be provided together")
        name = args.name or "smiles_alignment"
        OUT_ROOT.mkdir(parents=True, exist_ok=True)
        STAGE_ROOT.mkdir(parents=True, exist_ok=True)
        if args.save_alignment_files:
            ALIGNMENT_OUT_ROOT.mkdir(parents=True, exist_ok=True)
        try:
            rp_config = None
            direct_anchor_values = list(args.anchor or []) + list(
                args.subgraph_anchor or [])
            if direct_anchor_values:
                rp_config = rp_stage_config()
                rp_config['anchor_map'] = _subgraph_anchor_map_from_cli(
                    direct_anchor_values)
            rec = process_smiles_stage(
                name,
                args.reactant_smiles,
                args.product_smiles,
                workdir=Path(args.workdir) if args.workdir else None,
                stage=args.stage,
                inner_workers=(
                    args.inner_workers if args.inner_workers > 0
                    else max(1, int(args.workers))),
                save_alignment_files=args.save_alignment_files,
                rp_config=rp_config,
                resume_rp=args.resume_rp,
                component_spacing=args.smiles_component_spacing,
                expand_hydrogens=args.smiles_expand_hydrogens,
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
            source_dir = Path(args.workdir or (
                PROJECT / "out" / "smiles_work" / name)) / "source"
            print(f"{name}: mechs={slim.get('n_mechs', 0)} "
                  f"view={rec.get('view', {}).get('view_html')} "
                  f"source={source_dir}")
            EVAL_JSON.write_text(json.dumps([slim]))
        print(f"wrote {EVAL_JSON}")
        return

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
            rp_config = None
            direct_anchor_values = list(args.anchor or []) + list(
                args.subgraph_anchor or [])
            if direct_anchor_values:
                rp_config = rp_stage_config()
                rp_config['anchor_map'] = _subgraph_anchor_map_from_cli(
                    direct_anchor_values)
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
                rp_config=rp_config,
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
              f"xtb_threads={XTB_OMP_THREADS}, "
              f"xtb_workers={_resolve_xtb_workers(inner_workers)}, "
              f"charge={XTB_CHARGE}, "
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
              f"xtb_threads={XTB_OMP_THREADS}, "
              f"xtb_workers={_resolve_xtb_workers(0)}, "
              f"charge={XTB_CHARGE}, "
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
              f"xtb_threads={XTB_OMP_THREADS}, "
              f"xtb_workers={_resolve_xtb_workers(inner_workers)}, "
              f"charge={XTB_CHARGE}, "
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

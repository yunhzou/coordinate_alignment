#!/usr/bin/env python3
"""Rerun BGCP ranking on the appendix_final benchmark layout.

The release layout stores reactants/products/groundtruth/initial_guess as raw
XYZ files.  This adapter prepares combined endpoint XYZs when a side has
multiple fragments, runs the existing direct-XYZ BGCP pipeline with the
manifest charge/multiplicity, and aggregates a compact GT + top-2 IG viewer.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import sys
import traceback
from collections import Counter
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
SRC = PROJECT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import numpy as np

from rxn_core.chemistry_computations.xyz import parse_xyz, write_xyz_str

DEFAULT_BENCHMARK_ROOT = Path(
    "/h/399/yunhengzou/appendix_final/benchmark"
)
DEFAULT_RUN_ROOT = Path(
    "/h/399/yunhengzou/appendix_final/bgcp_rerank"
)


def safe_name(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return name or "item"


def read_manifest(benchmark_root: Path) -> list[dict]:
    path = benchmark_root / "manifest.csv"
    rows = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append({
                "step_id": row["step_id"],
                "charge": int(row["charge"]),
                "multiplicity": int(row["multiplicity"]),
                "reactant_count": int(row["reactant_count"]),
                "product_count": int(row["product_count"]),
                "groundtruth_count": int(row["groundtruth_count"]),
            })
    return rows


def select_manifest_row(rows: list[dict], index: int | None,
                        step: str | None) -> dict:
    if step:
        for row in rows:
            if row["step_id"] == step:
                return row
        raise SystemExit(f"no manifest row for step: {step}")
    if index is None:
        raise SystemExit("provide --index or --step")
    if index < 0 or index >= len(rows):
        raise SystemExit(f"--index out of range: {index} for {len(rows)} rows")
    return rows[index]


def xyz_files(directory: Path) -> list[Path]:
    return sorted(directory.glob("*.xyz"))


def pairwise_min_distance(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0 or len(b) == 0:
        return math.inf
    diff = a[:, None, :] - b[None, :, :]
    return float(np.linalg.norm(diff, axis=2).min())


def combine_fragments(role_dir: Path, out_path: Path, *,
                      min_separation: float = 5.0) -> dict:
    files = xyz_files(role_dir)
    if not files:
        raise FileNotFoundError(f"no XYZ files in {role_dir}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    all_elements: list[str] = []
    placed_blocks: list[np.ndarray] = []
    fragment_meta = []
    min_inter_fragment = None

    for idx, path in enumerate(files):
        elements, coords = parse_xyz(path)
        coords = np.asarray(coords, float)
        shift = np.zeros(3, dtype=float)

        if placed_blocks:
            placed = np.vstack(placed_blocks)
            shift[0] = float(placed[:, 0].max() + min_separation -
                             coords[:, 0].min())
            shifted = coords + shift
            distance = pairwise_min_distance(placed, shifted)
            if distance < min_separation:
                shift[0] += float(min_separation - distance + 1.0e-6)
                shifted = coords + shift
                distance = pairwise_min_distance(placed, shifted)
            min_inter_fragment = (
                distance if min_inter_fragment is None
                else min(min_inter_fragment, distance)
            )
        else:
            shifted = coords

        placed_blocks.append(shifted)
        all_elements.extend(elements)
        fragment_meta.append({
            "file": str(path),
            "n_atoms": len(elements),
            "translation": [float(x) for x in shift],
        })

    combined = np.vstack(placed_blocks)
    comment = (
        f"combined from {len(files)} fragment(s); "
        f"min_inter_fragment_distance={min_inter_fragment}"
    )
    out_path.write_text(write_xyz_str(all_elements, combined, comment))
    return {
        "path": str(out_path),
        "files": [str(p) for p in files],
        "n_fragments": len(files),
        "n_atoms": len(all_elements),
        "element_counts": dict(Counter(all_elements)),
        "min_inter_fragment_distance": min_inter_fragment,
        "fragments": fragment_meta,
    }


def initial_guess_label(path: Path) -> str:
    match = re.search(r"_iter(\d+)_", path.name)
    if match:
        return f"iter{int(match.group(1))}"
    return safe_name(path.stem)


def iter_sort_key(path: Path) -> tuple[int, str]:
    match = re.search(r"_iter(\d+)_", path.name)
    if match:
        return int(match.group(1)), path.name
    return 10**9, path.name


def validate_compositions(step_dir: Path, reactant_xyz: Path,
                          product_xyz: Path, target_paths: list[Path]) -> None:
    reactant_elements, _ = parse_xyz(reactant_xyz)
    product_elements, _ = parse_xyz(product_xyz)
    reference = Counter(reactant_elements)
    errors = []
    if Counter(product_elements) != reference:
        errors.append("product composition differs from reactant")
    for path in target_paths:
        elements, _ = parse_xyz(path)
        if Counter(elements) != reference:
            errors.append(f"target composition differs: {path}")
    if errors:
        detail = "\n".join(errors)
        raise ValueError(f"{step_dir.name}: nonmatching compositions\n{detail}")


def configure_pipeline(args, charge: int, multiplicity: int):
    import rxn_core.pipeline as pipeline

    run_root = args.run_root
    pipeline.OUT_ROOT = run_root / "views"
    pipeline.STAGE_ROOT = run_root / "stages"
    pipeline.ALIGNMENT_OUT_ROOT = run_root / "alignments"
    pipeline.EVAL_JSON = run_root / "bgcp_alignment_eval.json"
    pipeline.INCLUDE_GT = True
    pipeline.XTB_CACHE_MODE = args.xtb_mode
    pipeline.XTB_MAX_THREADS = max(1, int(args.xtb_max_threads))
    pipeline.XTB_OMP_THREADS = pipeline._resolve_xtb_threads(
        args.xtb_omp_threads, pipeline.XTB_MAX_THREADS)
    pipeline.XTB_WORKERS = args.xtb_workers
    pipeline.XTB_CHARGE = int(charge)
    pipeline.XTB_MULTIPLICITY = pipeline._normal_multiplicity(multiplicity)

    os.environ.update({
        "BGCP_OUT_ROOT": str(pipeline.OUT_ROOT),
        "BGCP_STAGE_ROOT": str(pipeline.STAGE_ROOT),
        "BGCP_ALIGNMENT_OUT_ROOT": str(pipeline.ALIGNMENT_OUT_ROOT),
        "BGCP_EVAL_JSON": str(pipeline.EVAL_JSON),
        "BGCP_INCLUDE_GT": "1",
        "BGCP_XTB_MODE": pipeline.XTB_CACHE_MODE,
        "BGCP_XTB_OMP_THREADS": str(pipeline.XTB_OMP_THREADS),
        "BGCP_XTB_MAX_THREADS": str(pipeline.XTB_MAX_THREADS),
        "BGCP_XTB_WORKERS": str(pipeline.XTB_WORKERS),
        "BGCP_CHARGE": str(charge),
        "BGCP_MULTIPLICITY": str(multiplicity),
        "RXN_CORE_PROJECT": str(PROJECT),
    })
    return pipeline


def remove_step_outputs(run_root: Path, step: str) -> None:
    for rel in ("combined_endpoints", "work", "stages", "views",
                "alignments"):
        path = run_root / rel / step
        if path.exists():
            shutil.rmtree(path)
    summary = run_root / "summaries" / f"{safe_name(step)}.json"
    if summary.exists():
        summary.unlink()


def prepare_step(row: dict, args) -> tuple[Path, Path, list[dict], dict]:
    benchmark_root = args.benchmark_root
    run_root = args.run_root
    step = row["step_id"]
    step_dir = benchmark_root / step
    combined_dir = run_root / "combined_endpoints" / step

    reactant_meta = combine_fragments(
        step_dir / "reactants",
        combined_dir / "reactant_combined.xyz",
        min_separation=args.fragment_separation,
    )
    product_meta = combine_fragments(
        step_dir / "products",
        combined_dir / "product_combined.xyz",
        min_separation=args.fragment_separation,
    )

    gt_files = xyz_files(step_dir / "groundtruth")
    ig_files = sorted(xyz_files(step_dir / "initial_guess"), key=iter_sort_key)
    if not gt_files:
        raise FileNotFoundError(f"{step}: no groundtruth XYZ")
    if not ig_files:
        raise FileNotFoundError(f"{step}: no initial_guess XYZ files")

    gt_used = gt_files[0]
    targets = [{
        "xyz": gt_used,
        "label": "GT",
        "kind": "gt",
    }]
    for path in ig_files:
        targets.append({
            "xyz": path,
            "label": initial_guess_label(path),
            "kind": "ig",
        })

    validate_compositions(
        step_dir,
        Path(reactant_meta["path"]),
        Path(product_meta["path"]),
        [gt_used] + ig_files,
    )
    metadata = {
        "reactant": reactant_meta,
        "product": product_meta,
        "groundtruth_used": str(gt_used),
        "groundtruth_extra": [str(p) for p in gt_files[1:]],
        "initial_guess_count": len(ig_files),
        "initial_guess_files": [str(p) for p in ig_files],
    }
    return Path(reactant_meta["path"]), Path(product_meta["path"]), targets, metadata


def summarize_record(rec: dict, metadata: dict) -> dict:
    slim = rec.get("slim") or {}
    mechanisms = slim.get("mechanisms") or []
    best_gt = None
    best_ig = None
    for mech in mechanisms:
        gt = mech.get("gt")
        if gt and gt.get("S") is not None:
            if best_gt is None or gt["S"] > best_gt["S"]:
                best_gt = {"mechanism_id": mech.get("id"), **gt}
        for ig in mech.get("igs") or []:
            if ig.get("S") is None:
                continue
            if best_ig is None or ig["S"] > best_ig["S"]:
                best_ig = {"mechanism_id": mech.get("id"), **ig}
    return {
        "status": "ok",
        "step": rec.get("step"),
        "n_mechs": slim.get("n_mechs"),
        "view_html": rec.get("view", {}).get("view_html"),
        "eval_slim_json": rec.get("view", {}).get("eval_slim_json"),
        "best_gt": best_gt,
        "best_ig": best_ig,
        "metadata": metadata,
    }


def run_step(args) -> int:
    rows = read_manifest(args.benchmark_root)
    row = select_manifest_row(rows, args.index, args.step)
    step = row["step_id"]
    args.run_root.mkdir(parents=True, exist_ok=True)
    for rel in ("logs", "summaries"):
        (args.run_root / rel).mkdir(parents=True, exist_ok=True)

    if args.force:
        remove_step_outputs(args.run_root, step)

    summary_path = args.run_root / "summaries" / f"{safe_name(step)}.json"
    if summary_path.exists() and not args.force and not args.dry_run:
        old = json.loads(summary_path.read_text())
        if old.get("status") == "ok":
            print(f"{step}: already complete; summary={summary_path}")
            return 0

    try:
        reactant_xyz, product_xyz, targets, metadata = prepare_step(row, args)
        if args.dry_run:
            summary = {
                "status": "dry_run",
                "step": step,
                "charge": row["charge"],
                "multiplicity": row["multiplicity"],
                "metadata": metadata,
            }
            summary_path.write_text(json.dumps(summary, indent=2))
            print(f"{step}: dry-run prepared endpoints")
            return 0

        pipeline = configure_pipeline(args, row["charge"], row["multiplicity"])
        rec = pipeline.process_xyz_stage(
            step,
            reactant_xyz,
            product_xyz,
            workdir=args.run_root / "work" / step,
            target_specs=targets,
            stage="full",
            inner_workers=max(1, int(args.workers)),
            save_alignment_files=args.save_alignment_files,
            charge=row["charge"],
            multiplicity=row["multiplicity"],
            xtb_mode=args.xtb_mode,
        )
        if rec.get("error"):
            raise RuntimeError(rec["error"])
        summary = summarize_record(rec, metadata)
        summary.update({
            "charge": row["charge"],
            "multiplicity": row["multiplicity"],
            "workers": args.workers,
            "xtb_mode": args.xtb_mode,
            "xtb_omp_threads": pipeline.XTB_OMP_THREADS,
            "xtb_workers": pipeline._resolve_xtb_workers(args.workers),
        })
        summary_path.write_text(json.dumps(summary, indent=2))
        print(f"{step}: ok view={summary.get('view_html')}")
        return 0
    except Exception as exc:
        summary = {
            "status": "error",
            "step": step,
            "charge": row.get("charge"),
            "multiplicity": row.get("multiplicity"),
            "error": f"{type(exc).__name__}: {exc}",
            "trace": traceback.format_exc(),
        }
        summary_path.write_text(json.dumps(summary, indent=2))
        print(f"{step}: ERROR {summary['error']}", file=sys.stderr)
        return 1


def finite(value):
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def choose_mechanism(ts_result: dict) -> dict | None:
    mechanisms = ts_result.get("mechanisms") or []
    if not mechanisms:
        return None
    with_gt = [m for m in mechanisms if m.get("gt") and
               finite(m["gt"].get("S")) is not None]
    if with_gt:
        return max(with_gt, key=lambda m: finite(m["gt"].get("S")) or -1.0)

    def best_ig_score(mech):
        scores = [finite(ig.get("S")) for ig in mech.get("igs") or []]
        scores = [s for s in scores if s is not None]
        return max(scores) if scores else -1.0

    return max(mechanisms, key=best_ig_score)


def panel_from_score(role: str, item: dict | None) -> dict | None:
    if not item:
        return None
    freq_summary = item.get("frequency_summary") or {}
    return {
        "role": role,
        "label": item.get("label") or role,
        "S": finite(item.get("S")),
        "beta": finite(item.get("beta")),
        "wbo_progress": finite(item.get("wbo_progress")),
        "mode_idx": item.get("k"),
        "freq": finite(item.get("freq")),
        "n_imag": freq_summary.get("n_imaginary"),
        "n_modes_total": freq_summary.get("n_modes"),
        "xyz_elements": item.get("elements") or [],
        "xyz_coords": item.get("xyz") or [],
        "disp": item.get("picked_disp") or [],
        "broken_bonds": item.get("broken_bonds_T") or [],
        "formed_bonds": item.get("formed_bonds_T") or [],
        "core_atoms": item.get("core_atoms_T") or [],
    }


def flat_record_for_step(run_root: Path, step: str) -> dict | None:
    ts_path = run_root / "stages" / step / "ts_stage.json"
    if not ts_path.exists():
        return None
    ts_result = json.loads(ts_path.read_text())
    mech = choose_mechanism(ts_result)
    if not mech:
        return None
    igs = [
        ig for ig in mech.get("igs") or []
        if finite(ig.get("S")) is not None
    ]
    igs = sorted(igs, key=lambda ig: finite(ig.get("S")) or -1.0,
                 reverse=True)
    panels = [
        panel_from_score("GT", mech.get("gt")),
        panel_from_score("IG #1", igs[0] if len(igs) > 0 else None),
        panel_from_score("IG #2", igs[1] if len(igs) > 1 else None),
    ]
    if any(panel is None for panel in panels):
        return None
    return {
        "step": step,
        "mechanism_id": mech.get("id"),
        "cut": mech.get("cut"),
        "panels": panels,
    }


def html_template(data_json: str) -> str:
    return f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>BGCP appendix rerank - GT and top-2 IG modes</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
body {{ font-family: Arial, sans-serif; margin: 12px; background: #f7f7f7; color: #222; }}
.ctl {{ background: #fff; border: 1px solid #d0d0d0; border-radius: 6px; padding: 8px 10px; display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }}
select {{ min-width: 520px; max-width: 70vw; padding: 4px 6px; }}
input[type="text"] {{ padding: 4px 6px; width: 180px; }}
button {{ padding: 4px 8px; }}
.row {{ display: flex; gap: 12px; margin-top: 12px; }}
.pane {{ flex: 1; min-width: 0; background: #fff; border: 1px solid #d0d0d0; border-radius: 6px; padding: 8px; }}
.viewer {{ width: 100%; height: 540px; position: relative; }}
h2 {{ margin: 0 0 10px; font-size: 18px; }}
h3 {{ margin: 0 0 4px; font-size: 14px; }}
.stats {{ font-size: 12px; color: #444; line-height: 1.35; min-height: 34px; }}
.legend span {{ display: inline-block; border-radius: 4px; padding: 2px 7px; font-size: 12px; margin-left: 4px; }}
.broken {{ background: #ffd9d9; color: #7a0000; }}
.formed {{ background: #d9f0d9; color: #005a00; }}
.core {{ background: #dbe8ff; color: #003c86; }}
</style>
</head>
<body>
<h2>BGCP appendix rerank - ground truth and top-2 initial guesses</h2>
<div class="ctl">
  <label><b>Step</b> <select id="stepSel"></select></label>
  <input id="filter" type="text" placeholder="filter step">
  <button id="prevBtn">Prev</button>
  <button id="nextBtn">Next</button>
  <label>Amplitude <input id="amp" type="range" min="0.05" max="1.5" step="0.05" value="0.5"> <span id="ampVal">0.50</span></label>
  <label>Speed <input id="speed" type="range" min="50" max="600" step="50" value="200"> <span id="speedVal">200</span> ms</label>
  <button id="playBtn">Pause</button>
  <span class="legend"><span class="broken">broken</span><span class="formed">formed</span><span class="core">mode vector</span></span>
</div>
<div class="row" id="panes"></div>
<script>
const DATA = {data_json};
const stepNames = Object.keys(DATA).sort();
const sel = document.getElementById("stepSel");
const panes = document.getElementById("panes");
let viewers = [null, null, null];
let current = null;
let amp = 0.5;
let speed = 200;
let playing = true;
let timer = null;
let phase = 0;

function fmt(v, digits) {{
  return (typeof v === "number" && isFinite(v)) ? v.toFixed(digits) : "n/a";
}}
function rebuildOptions() {{
  const filter = document.getElementById("filter").value.toLowerCase();
  sel.innerHTML = "";
  for (const name of stepNames) {{
    if (!name.toLowerCase().includes(filter)) continue;
    const d = DATA[name];
    const gt = d.panels[0], a = d.panels[1], b = d.panels[2];
    const opt = document.createElement("option");
    opt.value = name;
    opt.textContent = `${{name}}  mech=${{d.mechanism_id}}  GT=${{fmt(gt.S, 3)}}  IG1=${{a.label}}:${{fmt(a.S, 3)}}  IG2=${{b.label}}:${{fmt(b.S, 3)}}`;
    sel.appendChild(opt);
  }}
  if (sel.options.length) render(sel.value);
}}
function xyzBody(elements, coords, disp, scale) {{
  let out = `${{elements.length}}\\nframe\\n`;
  for (let i = 0; i < elements.length; i++) {{
    const p = coords[i], d = disp[i] || [0, 0, 0];
    out += `${{elements[i]}}  ${{(p[0] + scale * d[0]).toFixed(6)}}  ${{(p[1] + scale * d[1]).toFixed(6)}}  ${{(p[2] + scale * d[2]).toFixed(6)}}\\n`;
  }}
  return out;
}}
function buildPanes(d) {{
  panes.innerHTML = "";
  for (let i = 0; i < 3; i++) {{
    const p = d.panels[i];
    const div = document.createElement("div");
    div.className = "pane";
    div.innerHTML = `<h3>${{p.role}}: ${{p.label}}</h3>
      <div class="stats">S=${{fmt(p.S, 4)}} beta=${{fmt(p.beta, 4)}} wbo=${{fmt(p.wbo_progress, 4)}} mode=${{p.mode_idx ?? "n/a"}} freq=${{fmt(p.freq, 2)}} cm^-1 imag=${{p.n_imag ?? "n/a"}}/${{p.n_modes_total ?? "n/a"}}</div>
      <div id="viewer_${{i}}" class="viewer"></div>`;
    panes.appendChild(div);
  }}
}}
function drawPane(i, scale) {{
  const p = DATA[current].panels[i];
  const v = viewers[i];
  if (!v) return;
  v.removeAllModels();
  v.removeAllShapes();
  v.addModel(xyzBody(p.xyz_elements, p.xyz_coords, p.disp, scale), "xyz");
  v.setStyle({{}}, {{stick: {{radius: 0.1}}, sphere: {{scale: 0.2}}}});
  const atoms = v.selectedAtoms({{}});
  for (const pair of p.broken_bonds || []) {{
    const a = atoms[pair[0]], b = atoms[pair[1]];
    if (a && b) v.addCylinder({{start: {{x:a.x,y:a.y,z:a.z}}, end: {{x:b.x,y:b.y,z:b.z}}, color: "red", radius: 0.10, dashed: true}});
  }}
  for (const pair of p.formed_bonds || []) {{
    const a = atoms[pair[0]], b = atoms[pair[1]];
    if (a && b) v.addCylinder({{start: {{x:a.x,y:a.y,z:a.z}}, end: {{x:b.x,y:b.y,z:b.z}}, color: "green", radius: 0.10, dashed: true}});
  }}
  for (const idx of p.core_atoms || []) {{
    const a = atoms[idx], d = p.disp[idx];
    if (!a || !d) continue;
    const len = Math.hypot(d[0], d[1], d[2]);
    if (len < 1e-4) continue;
    v.addArrow({{start: {{x:a.x,y:a.y,z:a.z}}, end: {{x:a.x + d[0] * 1.5, y:a.y + d[1] * 1.5, z:a.z + d[2] * 1.5}}, color: "#005fcc", radius: 0.08}});
  }}
  v.render();
}}
function drawAll(scale) {{
  for (let i = 0; i < 3; i++) drawPane(i, scale);
}}
function render(name) {{
  current = name;
  buildPanes(DATA[name]);
  for (let i = 0; i < 3; i++) {{
    viewers[i] = $3Dmol.createViewer(`viewer_${{i}}`, {{backgroundColor: "white"}});
    drawPane(i, 0);
    viewers[i].zoomTo();
  }}
}}
function restartTimer() {{
  if (timer) clearInterval(timer);
  if (!playing) return;
  timer = setInterval(() => {{
    phase += 0.18;
    drawAll(Math.sin(phase) * amp);
  }}, speed);
}}
document.getElementById("filter").addEventListener("input", rebuildOptions);
document.getElementById("prevBtn").addEventListener("click", () => {{
  if (sel.selectedIndex > 0) {{ sel.selectedIndex--; render(sel.value); }}
}});
document.getElementById("nextBtn").addEventListener("click", () => {{
  if (sel.selectedIndex < sel.options.length - 1) {{ sel.selectedIndex++; render(sel.value); }}
}});
sel.addEventListener("change", () => render(sel.value));
document.getElementById("amp").addEventListener("input", (e) => {{
  amp = Number(e.target.value);
  document.getElementById("ampVal").textContent = amp.toFixed(2);
}});
document.getElementById("speed").addEventListener("input", (e) => {{
  speed = Number(e.target.value);
  document.getElementById("speedVal").textContent = String(speed);
  restartTimer();
}});
document.getElementById("playBtn").addEventListener("click", () => {{
  playing = !playing;
  document.getElementById("playBtn").textContent = playing ? "Pause" : "Play";
  restartTimer();
}});
window.addEventListener("load", () => {{
  rebuildOptions();
  restartTimer();
}});
</script>
</body>
</html>
"""


def aggregate(args) -> int:
    rows = read_manifest(args.benchmark_root)
    run_root = args.run_root
    aggregate_dir = run_root / "aggregate"
    aggregate_dir.mkdir(parents=True, exist_ok=True)
    records = []
    flat_data = {}

    for row in rows:
        step = row["step_id"]
        summary_path = run_root / "summaries" / f"{safe_name(step)}.json"
        summary = (
            json.loads(summary_path.read_text())
            if summary_path.exists()
            else {"status": "missing", "step": step}
        )
        flat = flat_record_for_step(run_root, step)
        if flat:
            flat_data[step] = flat
        records.append((row, summary, flat))

    csv_path = aggregate_dir / "rankings.csv"
    with csv_path.open("w", newline="") as f:
        fieldnames = [
            "step", "status", "charge", "multiplicity", "n_mechs",
            "mechanism_id", "gt_S", "gt_freq", "ig1_label", "ig1_S",
            "ig1_freq", "ig2_label", "ig2_S", "ig2_freq", "view_html",
            "error",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row, summary, flat in records:
            out = {
                "step": row["step_id"],
                "status": summary.get("status"),
                "charge": row["charge"],
                "multiplicity": row["multiplicity"],
                "n_mechs": summary.get("n_mechs"),
                "view_html": summary.get("view_html"),
                "error": summary.get("error"),
            }
            if flat:
                panels = flat["panels"]
                out.update({
                    "mechanism_id": flat.get("mechanism_id"),
                    "gt_S": panels[0].get("S"),
                    "gt_freq": panels[0].get("freq"),
                    "ig1_label": panels[1].get("label"),
                    "ig1_S": panels[1].get("S"),
                    "ig1_freq": panels[1].get("freq"),
                    "ig2_label": panels[2].get("label"),
                    "ig2_S": panels[2].get("S"),
                    "ig2_freq": panels[2].get("freq"),
                })
            writer.writerow(out)

    data_path = aggregate_dir / "flat_view_data.json"
    data_path.write_text(json.dumps(flat_data, separators=(",", ":")))
    html_path = aggregate_dir / "flat_view.html"
    html_path.write_text(html_template(json.dumps(
        flat_data, separators=(",", ":"))))
    status_path = aggregate_dir / "aggregate_status.json"
    status_path.write_text(json.dumps({
        "n_manifest": len(rows),
        "n_summaries": sum(1 for _, s, _ in records
                           if s.get("status") in {"ok", "error", "dry_run"}),
        "n_ok": sum(1 for _, s, _ in records if s.get("status") == "ok"),
        "n_flat_view_steps": len(flat_data),
        "rankings_csv": str(csv_path),
        "flat_view_data": str(data_path),
        "flat_view_html": str(html_path),
    }, indent=2))
    print(f"aggregate: {len(flat_data)} flat-view steps")
    print(f"aggregate: wrote {html_path}")
    return 0


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--benchmark-root", type=Path,
                        default=DEFAULT_BENCHMARK_ROOT)
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run-step")
    add_common_args(run)
    run.add_argument("--index", type=int)
    run.add_argument("--step")
    run.add_argument("--workers", type=int,
                     default=int(os.environ.get("SLURM_CPUS_PER_TASK", "8")))
    run.add_argument("--xtb-mode", choices=("auto", "cache-only"),
                     default="auto")
    run.add_argument("--xtb-omp-threads", default="auto")
    run.add_argument("--xtb-max-threads", type=int, default=8)
    run.add_argument("--xtb-workers", default="auto")
    run.add_argument("--fragment-separation", type=float, default=5.0)
    run.add_argument("--save-alignment-files", action="store_true")
    run.add_argument("--force", action="store_true")
    run.add_argument("--dry-run", action="store_true")

    agg = sub.add_parser("aggregate")
    add_common_args(agg)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "run-step":
        return run_step(args)
    if args.command == "aggregate":
        return aggregate(args)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())

"""Local web app for interactive R-P AAM and weighted subgraph matching."""
from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
import time
import traceback
import uuid
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

import numpy as np

from .frag import WeightedGraph
from .pipeline import (
    build_view_data,
    rp_stage_config,
    run_rp_stage,
    smiles_inputs_from_strings,
    step_inputs_from_arrays,
    write_stage_json,
    write_view_stage,
)
from .smiles import smiles_to_formal_wbo, smiles_to_weighted_graph
from .subgraph import match_weighted_subgraph
from .session_store import MongoSessionStore


COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Li": 1.28, "Na": 1.66, "K": 2.03, "Mg": 1.41, "Ca": 1.76,
    "Al": 1.21, "Si": 1.11, "Fe": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ag": 1.45, "Au": 1.36, "Pt": 1.36, "Pd": 1.39, "V": 1.53,
    "Mo": 1.54, "W": 1.62,
}


def _safe_name(value, fallback="job"):
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("_")
    return text or fallback


def _normal_element(raw):
    text = str(raw or "").strip()
    if not text:
        return "X"
    return text[0].upper() + text[1:].lower()


def parse_xyz_text(text):
    """Parse one XYZ block from text."""
    lines = [line.rstrip() for line in str(text or "").splitlines()]
    while lines and not lines[0].strip():
        lines.pop(0)
    if len(lines) < 2:
        raise ValueError("XYZ input must contain an atom count and comment")
    try:
        n_atoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError("first XYZ line must be an integer atom count") from exc
    if n_atoms < 0:
        raise ValueError("XYZ atom count must be non-negative")
    if len(lines) < n_atoms + 2:
        raise ValueError(f"XYZ has {len(lines) - 2} atom rows, expected {n_atoms}")
    elements = []
    coords = np.zeros((n_atoms, 3), dtype=float)
    for idx in range(n_atoms):
        parts = lines[idx + 2].split()
        if len(parts) < 4:
            raise ValueError(f"XYZ atom row {idx} must have element x y z")
        elements.append(_normal_element(parts[0]))
        coords[idx] = [float(parts[1]), float(parts[2]), float(parts[3])]
    return elements, coords


def parse_wbo_text(text, n_atoms):
    """Parse a square matrix JSON/text block or an xTB-style ``i j w`` list."""
    raw = str(text or "").strip()
    if not raw:
        return None
    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        data = None
    if data is not None:
        matrix = np.asarray(data.get("weights", data.get("wbo"))
                            if isinstance(data, dict) else data, dtype=float)
        if matrix.shape != (n_atoms, n_atoms):
            raise ValueError(
                f"WBO matrix shape must be {(n_atoms, n_atoms)}, got {matrix.shape}")
        return matrix

    token_rows = [
        line.split()
        for line in raw.splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    if not token_rows:
        return None

    def _is_int_token(value):
        try:
            int(value)
            return True
        except ValueError:
            return False

    looks_like_triples = all(
        len(parts) >= 3
        and _is_int_token(parts[0])
        and _is_int_token(parts[1])
        and 1 <= int(parts[0]) <= n_atoms
        and 1 <= int(parts[1]) <= n_atoms
        for parts in token_rows
    )
    looks_like_matrix = (
        len(token_rows) == n_atoms
        and all(len(parts) == n_atoms for parts in token_rows)
    )
    if looks_like_matrix and not (
            all(len(parts) == 3 for parts in token_rows) and looks_like_triples):
        return np.asarray(
            [[float(v) for v in parts] for parts in token_rows], dtype=float)

    triples = []
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            triples.append((int(parts[0]) - 1, int(parts[1]) - 1, float(parts[2])))
            continue
        raise ValueError(f"cannot parse WBO line: {line!r}")
    matrix = np.zeros((n_atoms, n_atoms), dtype=float)
    for i, j, w in triples:
        if i < 0 or j < 0 or i >= n_atoms or j >= n_atoms:
            raise ValueError("WBO pair index is outside atom range")
        matrix[i, j] = matrix[j, i] = float(w)
    return matrix


def infer_single_bond_wbo(elements, coords, *, scale=1.25):
    """Distance-based single-bond guess for UI-only XYZ graph construction."""
    coords = np.asarray(coords, dtype=float)
    n_atoms = len(elements)
    matrix = np.zeros((n_atoms, n_atoms), dtype=float)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            ri = COVALENT_RADII.get(_normal_element(elements[i]), 0.77)
            rj = COVALENT_RADII.get(_normal_element(elements[j]), 0.77)
            dist = float(np.linalg.norm(coords[i] - coords[j]))
            if 0.25 < dist <= scale * (ri + rj):
                matrix[i, j] = matrix[j, i] = 1.0
    return matrix


def _bonds_from_wbo(wbo, floor=0.2):
    matrix = np.asarray(wbo, dtype=float)
    bonds = []
    for i in range(matrix.shape[0]):
        for j in range(i + 1, matrix.shape[1]):
            w = float(matrix[i, j])
            if w >= floor:
                bonds.append({"i": i, "j": j, "wbo": w})
    return bonds


def _endpoint_record(elements, coords, wbo=None, atom_maps=None, source=None):
    return {
        "elements": list(elements),
        "coords": np.asarray(coords, dtype=float).tolist(),
        "wbo": None if wbo is None else np.asarray(wbo, dtype=float).tolist(),
        "bonds": [] if wbo is None else _bonds_from_wbo(wbo),
        "atom_maps": {str(k): int(v) for k, v in sorted((atom_maps or {}).items())},
        "source": source,
    }


def _config_from_payload(payload):
    cfg = rp_stage_config()
    user_cfg = dict(payload.get("config") or {})
    for key in (
            "iso_tol", "dwbo_threshold", "metal_dwbo_threshold",
            "symmetry_wbo_tol"):
        if key in user_cfg and user_cfg[key] not in (None, ""):
            cfg[key] = float(user_cfg[key])
    for key in ("n_seeds", "max_branches"):
        if key in user_cfg and user_cfg[key] not in (None, ""):
            cfg[key] = int(user_cfg[key])
    anchors = payload.get("anchors") or {}
    cfg["anchor_map"] = {int(k): int(v) for k, v in dict(anchors).items()}
    return cfg


def _xyz_inputs_from_payload(payload, job_dir, name):
    r_el, r_xyz = parse_xyz_text(payload.get("reactant", {}).get("xyz", ""))
    p_el, p_xyz = parse_xyz_text(payload.get("product", {}).get("xyz", ""))
    mode = str(payload.get("xyz_wbo_mode") or "infer").lower()
    if mode == "provided":
        r_wbo = parse_wbo_text(payload.get("reactant", {}).get("wbo", ""), len(r_el))
        p_wbo = parse_wbo_text(payload.get("product", {}).get("wbo", ""), len(p_el))
        if r_wbo is None or p_wbo is None:
            raise ValueError("provided WBO mode requires both R and P WBO inputs")
        return step_inputs_from_arrays(
            name, r_el, r_xyz, r_wbo, p_el, p_xyz, p_wbo, step_dir=job_dir)
    if mode == "infer":
        return step_inputs_from_arrays(
            name,
            r_el, r_xyz, infer_single_bond_wbo(r_el, r_xyz),
            p_el, p_xyz, infer_single_bond_wbo(p_el, p_xyz),
            step_dir=job_dir)
    if mode == "xtb":
        from .pipeline import alignment_inputs_from_xyz
        r_path = job_dir / "R.xyz"
        p_path = job_dir / "P.xyz"
        r_path.write_text(payload.get("reactant", {}).get("xyz", ""))
        p_path.write_text(payload.get("product", {}).get("xyz", ""))
        xtb_mode = str(payload.get("xtb_mode") or "auto")
        charge = int(payload.get("charge") or 0)
        multiplicity = int(payload.get("multiplicity") or 1)
        return alignment_inputs_from_xyz(
            r_path, p_path, workdir=job_dir / "xtb", name=name,
            charge=charge, multiplicity=multiplicity, xtb_mode=xtb_mode)
    raise ValueError("xyz_wbo_mode must be 'infer', 'provided', or 'xtb'")


def _build_inputs_from_payload(payload, job_dir, name):
    mode = str(payload.get("mode") or "smiles").lower()
    if mode == "smiles":
        return smiles_inputs_from_strings(
            payload.get("reactant", {}).get("smiles", ""),
            payload.get("product", {}).get("smiles", ""),
            name=name, workdir=job_dir)
    if mode == "xyz":
        return _xyz_inputs_from_payload(payload, job_dir, name)
    raise ValueError("mode must be 'smiles' or 'xyz'")


def preview_payload(payload):
    mode = str(payload.get("mode") or "smiles").lower()
    if mode == "smiles":
        r = smiles_to_formal_wbo(payload.get("reactant", {}).get("smiles", ""))
        p = smiles_to_formal_wbo(payload.get("product", {}).get("smiles", ""))
        return {
            "mode": "smiles",
            "reactant": _endpoint_record(
                r.elements, r.coords, r.wbo, r.atom_maps, "smiles"),
            "product": _endpoint_record(
                p.elements, p.coords, p.wbo, p.atom_maps, "smiles"),
        }
    if mode == "xyz":
        r_el, r_xyz = parse_xyz_text(payload.get("reactant", {}).get("xyz", ""))
        p_el, p_xyz = parse_xyz_text(payload.get("product", {}).get("xyz", ""))
        wbo_mode = str(payload.get("xyz_wbo_mode") or "infer").lower()
        r_wbo = p_wbo = None
        if wbo_mode == "provided":
            r_wbo = parse_wbo_text(payload.get("reactant", {}).get("wbo", ""), len(r_el))
            p_wbo = parse_wbo_text(payload.get("product", {}).get("wbo", ""), len(p_el))
        elif wbo_mode == "infer":
            r_wbo = infer_single_bond_wbo(r_el, r_xyz)
            p_wbo = infer_single_bond_wbo(p_el, p_xyz)
        return {
            "mode": "xyz",
            "reactant": _endpoint_record(r_el, r_xyz, r_wbo, source="xyz"),
            "product": _endpoint_record(p_el, p_xyz, p_wbo, source="xyz"),
        }
    raise ValueError("mode must be 'smiles' or 'xyz'")


def _drawn_graph_to_weighted(spec):
    graph = spec.get("graph") or {}
    nodes_raw = list(graph.get("nodes") or [])
    nodes = []
    for node in nodes_raw:
        element = _normal_element(node.get("element", "C"))
        nodes.append({
            "element": element,
            "features": dict(node.get("features") or {}),
            "label": node.get("label"),
        })
    weights = np.zeros((len(nodes), len(nodes)), dtype=float)
    for edge in graph.get("edges") or []:
        i = int(edge.get("i"))
        j = int(edge.get("j"))
        order = float(edge.get("order", edge.get("wbo", 1.0)))
        if i == j:
            continue
        weights[i, j] = weights[j, i] = order
    coords = [
        [float(node.get("x", 0.0)), float(node.get("y", 0.0)), 0.0]
        for node in nodes_raw
    ]
    return WeightedGraph(nodes=nodes, weights=weights, coords=np.asarray(coords))


def _weighted_json_to_graph(spec):
    data = spec.get("graph_json")
    if isinstance(data, str):
        data = json.loads(data)
    if not isinstance(data, dict):
        raise ValueError("weighted JSON graph must be an object")
    weights = data.get("weights", data.get("wbo"))
    if weights is None:
        raise ValueError("weighted JSON graph must contain weights or wbo")
    return WeightedGraph(
        nodes=data.get("nodes") or [],
        weights=np.asarray(weights, dtype=float),
        weight_name=data.get("weight_name", "wbo"),
        coords=(None if data.get("coords") is None
                else np.asarray(data.get("coords"), dtype=float)),
        metadata=dict(data.get("metadata") or {}),
    )


def graph_from_subgraph_spec(spec):
    kind = str(spec.get("type") or "drawn").lower()
    if kind == "drawn":
        return _drawn_graph_to_weighted(spec)
    if kind == "smiles":
        return smiles_to_weighted_graph(spec.get("smiles", ""))
    if kind == "xyz":
        elements, coords = parse_xyz_text(spec.get("xyz", ""))
        wbo = parse_wbo_text(spec.get("wbo", ""), len(elements))
        if wbo is None:
            wbo = infer_single_bond_wbo(elements, coords)
        nodes = [{"element": element, "features": {}} for element in elements]
        return WeightedGraph(nodes=nodes, weights=wbo, coords=coords)
    if kind == "json":
        return _weighted_json_to_graph(spec)
    raise ValueError("subgraph source type must be drawn, smiles, xyz, or json")


def run_subgraph_payload(payload):
    query = graph_from_subgraph_spec(payload.get("query") or {})
    target = graph_from_subgraph_spec(payload.get("target") or {})
    node_policy = payload.get("node_policy") or None
    if isinstance(node_policy, list) and not node_policy:
        node_policy = None
    anchors = {int(k): int(v) for k, v in dict(payload.get("anchors") or {}).items()}
    iso_tol = float(payload.get("iso_tol") or 1.0)
    matches = match_weighted_subgraph(
        query,
        target,
        node_policy=node_policy,
        anchor_map=anchors,
        graph_floor=float(payload.get("graph_floor") or 0.2),
        iso_tol=iso_tol,
        symmetry_wbo_tol=iso_tol,
        orbit_dedup=not bool(payload.get("no_orbit_dedup")),
    )
    return {
        "n_matches": len(matches),
        "matches": [
            {
                "mapping": {str(k): int(v) for k, v in sorted(m.mapping.items())},
                "query_nodes": [int(v) for v in m.query_nodes],
                "target_nodes": [int(v) for v in m.target_nodes],
                "deferred_edges": [list(map(int, e)) for e in m.deferred_edges],
                "symmetry_fragments": list(m.symmetry_fragments),
            }
            for m in matches
        ],
    }


class AAMWebApp:
    """Stateful request handler target used by ``AAMRequestHandler``."""

    def __init__(self, work_root, *, mongo_uri=None, mongo_db="rxn_core",
                 session_store=None):
        self.work_root = Path(work_root)
        self.work_root.mkdir(parents=True, exist_ok=True)
        self.session_error = None
        self.mongo_uri = mongo_uri
        self.mongo_db = mongo_db
        self.session_store = session_store
        self._ensure_session_store()

    def _ensure_session_store(self):
        if self.session_store is not None or not self.mongo_uri:
            return self.session_store
        try:
            self.session_store = MongoSessionStore(
                self.mongo_uri, db_name=self.mongo_db)
            self.session_error = None
        except Exception as exc:
            self.session_error = f"{type(exc).__name__}: {exc}"
        return self.session_store

    def _require_session_store(self):
        store = self._ensure_session_store()
        if store is None:
            raise RuntimeError(
                self.session_error or "Mongo session store is not configured")
        return store

    def new_job_dir(self, name):
        stamp = time.strftime("%Y%m%d_%H%M%S")
        job_id = f"{stamp}_{_safe_name(name)}_{uuid.uuid4().hex[:8]}"
        job_dir = self.work_root / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        return job_id, job_dir

    def api_preview(self, payload):
        return preview_payload(payload)

    def api_rp(self, payload):
        name = _safe_name(payload.get("name") or "aam_alignment", "aam_alignment")
        job_id, job_dir = self.new_job_dir(name)
        inputs = _build_inputs_from_payload(payload, job_dir, name)
        cfg = _config_from_payload(payload)
        rp_result = run_rp_stage(
            inputs, config=cfg,
            inner_workers=max(1, int(payload.get("inner_workers") or 1)))
        write_stage_json(job_dir / "rp_stage.json", rp_result)
        view_result = write_view_stage(
            inputs, rp_result, ts_result=None, out_root=job_dir / "views",
            include_gt=False)
        view_path = Path(view_result["view_html"])
        rel_view = view_path.relative_to(job_dir)
        view_data = build_view_data(inputs, rp_result, ts_result=None,
                                    include_gt=False)
        reactant = _endpoint_record(
            inputs.elR, inputs.xyzR, inputs.wboR, source="reactant")
        reactant["metadata"] = view_data["reactant"].get("metadata", {})
        product = _endpoint_record(
            inputs.elP, inputs.xyzP, inputs.wboP, source="product")
        product["metadata"] = view_data["product"].get("metadata", {})
        return {
            "job_id": job_id,
            "name": name,
            "rp": rp_result,
            "view": view_result,
            "viewer_url": f"/artifacts/{job_id}/{rel_view.as_posix()}",
            "stage_url": f"/artifacts/{job_id}/rp_stage.json",
            "viewer_rel_path": rel_view.as_posix(),
            "stage_rel_path": "rp_stage.json",
            "reactant": reactant,
            "product": product,
            "work_dir": str(job_dir),
        }

    def api_subgraph(self, payload):
        return run_subgraph_payload(payload)

    def api_session_status(self):
        self._ensure_session_store()
        return {
            "enabled": self.session_store is not None,
            "error": self.session_error,
        }

    def api_list_sessions(self, limit=50):
        store = self._ensure_session_store()
        if store is None:
            return {
                "enabled": False,
                "error": self.session_error or "Mongo session store is not configured",
                "sessions": [],
            }
        return {
            "enabled": True,
            "sessions": store.list_sessions(limit=limit),
        }

    def api_save_session(self, payload):
        store = self._require_session_store()
        state = dict(payload.get("state") or {})
        name = payload.get("name") or state.get("name") or "session"
        session_id = payload.get("id") or state.get("session_id")
        job_id = payload.get("job_id") or (
            state.get("lastRun") or {}).get("job_id")
        job_dir = self.work_root / str(job_id) if job_id else None
        saved = store.save_session(
            name, state, job_dir=job_dir, session_id=session_id)
        return {"enabled": True, "session": saved}

    def api_load_session(self, session_id):
        store = self._require_session_store()
        return store.load_session(
            session_id, restore_root=self.work_root)

    def api_delete_session(self, session_id):
        store = self._require_session_store()
        return {
            "deleted": store.delete_session(session_id),
        }


class AAMRequestHandler(SimpleHTTPRequestHandler):
    """Small JSON API and artifact server."""

    server_version = "rxn-core-aam-webapp/0.1"

    def log_message(self, fmt, *args):  # noqa: D401 - stdlib hook
        if os.environ.get("RXN_CORE_WEBAPP_LOG", "0") == "1":
            super().log_message(fmt, *args)

    @property
    def app(self):
        return self.server.app

    def _send_json(self, data, status=200):
        raw = json.dumps(data).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def _send_html(self, html):
        raw = html.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def _read_payload(self):
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length else b"{}"
        return json.loads(raw.decode("utf-8") or "{}")

    def do_GET(self):  # noqa: N802 - stdlib hook
        parsed = urlparse(self.path)
        try:
            if parsed.path in {"/", "/index.html"}:
                self._send_html(APP_HTML)
                return
            if parsed.path == "/health":
                self._send_json({"ok": True})
                return
            if parsed.path == "/api/session-status":
                self._send_json(self.app.api_session_status())
                return
            if parsed.path == "/api/sessions":
                self._send_json(self.app.api_list_sessions())
                return
            if parsed.path.startswith("/api/sessions/"):
                session_id = unquote(parsed.path[len("/api/sessions/"):])
                self._send_json(self.app.api_load_session(session_id))
                return
            if parsed.path.startswith("/artifacts/"):
                self._serve_artifact(parsed.path)
                return
            self.send_error(404, "not found")
        except Exception as exc:
            self._send_json({
                "error": f"{type(exc).__name__}: {exc}",
                "trace": traceback.format_exc(),
            }, status=500)

    def _serve_artifact(self, path):
        rel = unquote(path[len("/artifacts/"):])
        parts = [p for p in rel.split("/") if p]
        if not parts:
            self.send_error(404, "missing artifact path")
            return
        job_id = parts[0]
        rest = Path(*parts[1:]) if len(parts) > 1 else Path("")
        root = (self.app.work_root / job_id).resolve()
        target = (root / rest).resolve()
        if root not in target.parents and target != root:
            self.send_error(403, "artifact path outside job")
            return
        if not target.exists() or not target.is_file():
            self.send_error(404, "artifact not found")
            return
        ctype = "text/html" if target.suffix == ".html" else "application/json"
        raw = target.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def do_POST(self):  # noqa: N802 - stdlib hook
        parsed = urlparse(self.path)
        try:
            payload = self._read_payload()
            if parsed.path == "/api/preview":
                self._send_json(self.app.api_preview(payload))
                return
            if parsed.path == "/api/rp":
                self._send_json(self.app.api_rp(payload))
                return
            if parsed.path == "/api/subgraph":
                self._send_json(self.app.api_subgraph(payload))
                return
            if parsed.path == "/api/sessions":
                self._send_json(self.app.api_save_session(payload))
                return
            self.send_error(404, "not found")
        except Exception as exc:
            self._send_json({
                "error": f"{type(exc).__name__}: {exc}",
                "trace": traceback.format_exc(),
            }, status=500)

    def do_DELETE(self):  # noqa: N802 - stdlib hook
        parsed = urlparse(self.path)
        try:
            if parsed.path.startswith("/api/sessions/"):
                session_id = unquote(parsed.path[len("/api/sessions/"):])
                self._send_json(self.app.api_delete_session(session_id))
                return
            self.send_error(404, "not found")
        except Exception as exc:
            self._send_json({
                "error": f"{type(exc).__name__}: {exc}",
                "trace": traceback.format_exc(),
            }, status=500)


class AAMHTTPServer(ThreadingHTTPServer):
    def __init__(self, server_address, app):
        super().__init__(server_address, AAMRequestHandler)
        self.app = app


APP_HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>rxn_core AAM Workbench</title>
  <script src="https://3dmol.org/build/3Dmol-min.js"></script>
  <style>
    :root{--line:#d8dde3;--ink:#17212b;--muted:#5d6875;--bg:#f6f7f9;--panel:#fff;--blue:#0f5f9f;--red:#c62828;--green:#16803c;--gold:#b88a00}
    *{box-sizing:border-box}
    html,body{margin:0;padding:0;background:var(--bg);color:var(--ink);font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}
    body{padding:12px}
    header{display:flex;align-items:center;gap:12px;margin-bottom:10px}
    h1{font-size:18px;margin:0;font-weight:650}
    .tabs{display:flex;gap:6px}
    .sessionbar{margin-left:auto}
    .sessionbar input{width:150px}
    .sessionbar select{width:190px}
    .tab{border:1px solid var(--line);background:#fff;border-radius:4px;padding:7px 11px;cursor:pointer;font-size:13px;color:#253342}
    .tab.active{border-color:#7aa7c7;background:#e9f3fb;color:#073c61;font-weight:600}
    .page{display:none}
    .page.active{display:block}
    .grid{display:grid;grid-template-columns:360px minmax(0,1fr);gap:12px;align-items:start}
    .panel{background:var(--panel);border:1px solid var(--line);border-radius:6px;padding:10px}
    .section{display:grid;gap:10px;margin-bottom:10px}
    .row{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
    .twocol{display:grid;grid-template-columns:1fr 1fr;gap:10px}
    label{font-size:12px;color:var(--muted);display:grid;gap:4px}
    label.toggle{display:flex;align-items:center;gap:6px}
    label.toggle input{margin:0}
    input,select,textarea,button{font:inherit}
    input,select,textarea{border:1px solid #c8d0d8;border-radius:4px;background:#fff;color:var(--ink);padding:6px}
    textarea{width:100%;min-height:86px;resize:vertical;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
    textarea.tall{min-height:138px}
    button{border:1px solid #9aa8b5;border-radius:4px;background:#fff;color:#073c61;padding:7px 10px;cursor:pointer}
    button.primary{background:#0f5f9f;color:white;border-color:#0f5f9f}
    button.warn{border-color:#bb8b00;color:#704f00;background:#fff8dd}
    button:disabled{opacity:.55;cursor:not-allowed}
    .metric{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;color:#314151}
    .viewergrid{display:grid;grid-template-columns:1fr 1fr;gap:10px}
    .viewergrid.compact .molbox{height:260px}
    .molbox{height:315px;border:1px solid var(--line);border-radius:6px;position:relative;background:#fff;overflow:hidden}
    .mol2d svg{width:100%;height:100%;display:block}
    .moltitle{position:absolute;left:8px;top:6px;z-index:2;background:rgba(255,255,255,.82);border:1px solid var(--line);border-radius:4px;padding:2px 6px;font-size:12px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
    .empty{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;color:#6d7885;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
    .anchors{max-height:162px;overflow:auto;border:1px solid var(--line);border-radius:4px}
    table{width:100%;border-collapse:collapse;font-size:12px}
    th,td{border-bottom:1px solid var(--line);padding:5px;text-align:left}
    th{background:#f9fafb;color:#52606d;font-weight:600}
    .pill{display:inline-flex;align-items:center;gap:4px;border:1px solid var(--line);border-radius:4px;padding:2px 5px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
    .swatch{display:inline-block;width:10px;height:10px;border-radius:50%;border:1px solid rgba(0,0,0,.24)}
    iframe{width:100%;height:640px;border:1px solid var(--line);border-radius:6px;background:#fff}
    pre{margin:0;white-space:pre-wrap;word-break:break-word;background:#f9fafb;border:1px solid var(--line);border-radius:4px;padding:8px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;max-height:270px;overflow:auto}
    .drawer{height:330px;border:1px solid var(--line);border-radius:6px;background:#fff}
    .small{font-size:12px;color:var(--muted)}
    .error{color:#a21919}
    .ok{color:#126b35}
    .subgrid{display:grid;grid-template-columns:minmax(280px,360px) minmax(280px,360px) minmax(0,1fr);gap:12px;align-items:start}
    .hidden{display:none!important}
    @media(max-width:1120px){.grid,.subgrid{grid-template-columns:1fr}.viewergrid{grid-template-columns:1fr}iframe{height:560px}}
  </style>
</head>
<body>
  <header>
    <h1>rxn_core AAM Workbench</h1>
    <div class="tabs">
      <button class="tab active" data-page="rpPage">R-P AAM</button>
      <button class="tab" data-page="subPage">Subgraph</button>
    </div>
    <div class="row sessionbar">
      <input id="sessionName" placeholder="session name">
      <button id="saveSessionBtn">Save</button>
      <button id="refreshSessionsBtn">Refresh</button>
      <select id="sessionSelect"><option value="">No saved sessions</option></select>
      <button id="loadSessionBtn">Load</button>
      <button id="deleteSessionBtn" class="warn">Delete</button>
    </div>
    <span id="status" class="metric"></span>
  </header>

  <main id="rpPage" class="page active">
    <div class="grid">
      <div class="panel">
        <div class="section">
          <div class="row">
            <label>Input
              <select id="rpMode">
                <option value="smiles">SMILES / CXSMILES</option>
                <option value="xyz">XYZ</option>
              </select>
            </label>
            <label>Name <input id="rpName" value="aam_alignment"></label>
            <label>Result
              <select id="rpRenderMode">
                <option value="2d">2D graph</option>
                <option value="3d">3D XYZ</option>
              </select>
            </label>
          </div>
          <div id="smilesInputs" class="twocol">
            <label>Reactant SMILES<textarea id="rSmiles" class="tall">[O:1][H:2]</textarea></label>
            <label>Product SMILES<textarea id="pSmiles" class="tall">[O-:1].[H+:2]</textarea></label>
          </div>
          <div id="xyzInputs" class="hidden">
            <div class="row">
              <input id="xyzWboMode" type="hidden" value="xtb">
              <label>Charge <input id="charge" type="number" value="0" style="width:74px"></label>
              <label>Multiplicity <input id="multiplicity" type="number" value="1" style="width:74px"></label>
            </div>
            <div class="twocol">
              <label>Reactant XYZ<input type="file" id="rXyzFile" accept=".xyz,text/plain"><textarea id="rXyz" class="tall"></textarea></label>
              <label>Product XYZ<input type="file" id="pXyzFile" accept=".xyz,text/plain"><textarea id="pXyz" class="tall"></textarea></label>
            </div>
            <div id="wboInputs" class="hidden">
              <textarea id="rWbo"></textarea>
              <textarea id="pWbo"></textarea>
            </div>
          </div>
          <div class="twocol">
            <label>iso_tol <input id="isoTol" type="number" value="1.0" step="0.1"></label>
            <label>dwbo_threshold <input id="dwboThreshold" type="number" value="0.5" step="0.1"></label>
            <label>n_seeds <input id="nSeeds" type="number" value="3" min="1"></label>
          </div>
          <div class="row">
            <button id="previewBtn">Preview</button>
            <button id="runRpBtn" class="primary">Run R-P AAM</button>
            <button id="clearAnchorsBtn" class="warn">Clear Anchors</button>
          </div>
        </div>
        <div class="section">
          <div class="row">
            <span class="metric">Selected R: <b id="selR">-</b></span>
            <span class="metric">Selected P: <b id="selP">-</b></span>
          </div>
          <div class="anchors">
            <table>
              <thead><tr><th>R</th><th>P</th><th></th></tr></thead>
              <tbody id="anchorRows"></tbody>
            </table>
          </div>
          <pre id="anchorJson">{}</pre>
        </div>
      </div>

      <div>
        <div class="viewergrid">
          <div class="molbox" id="rpRBox"><div class="moltitle">R</div><div class="empty">Preview R</div></div>
          <div class="molbox" id="rpPBox"><div class="moltitle">P</div><div class="empty">Preview P</div></div>
        </div>
        <div class="panel" style="margin-top:10px">
          <div class="row" style="justify-content:space-between">
            <span class="metric" id="rpSummary">No run yet</span>
            <a id="stageLink" class="small" target="_blank"></a>
          </div>
          <div id="rp2dResult" class="hidden" style="margin-top:8px">
            <div class="row" style="margin-bottom:8px">
              <div class="row" id="rp2dMechs"></div>
              <label class="toggle"><input type="checkbox" id="show2dDegeneracy"> Degeneracy</label>
            </div>
            <div class="viewergrid compact">
              <div class="molbox" id="rp2dR"><div class="moltitle">R 2D</div><div class="empty">Run R-P AAM</div></div>
              <div class="molbox" id="rp2dP"><div class="moltitle">P 2D</div><div class="empty">Run R-P AAM</div></div>
            </div>
          </div>
          <iframe id="rpViewer" class="hidden"></iframe>
          <pre id="rpJson" style="margin-top:8px"></pre>
        </div>
      </div>
    </div>
  </main>

  <main id="subPage" class="page">
    <div class="subgrid">
      <div class="panel">
        <div class="row" style="justify-content:space-between">
          <b>Query</b>
          <select id="qType">
            <option value="drawn">Drawn</option>
            <option value="smiles">SMILES</option>
            <option value="xyz">XYZ</option>
            <option value="json">Weighted JSON</option>
          </select>
        </div>
        <div id="qDrawControls" class="section">
          <div class="row">
            <label>Mode <select id="qDrawMode"><option value="atom">Atom</option><option value="bond">Bond</option><option value="delete">Delete</option></select></label>
            <label>Element <input id="qElement" value="C" style="width:70px"></label>
            <label>Bond <select id="qBond"><option value="1">1</option><option value="1.5">1.5</option><option value="2">2</option><option value="3">3</option></select></label>
            <button data-clear-draw="q">Clear</button>
          </div>
          <svg id="qDraw" class="drawer"></svg>
        </div>
        <label id="qSmilesWrap" class="hidden">Query SMILES<textarea id="qSmiles"></textarea></label>
        <div id="qXyzWrap" class="hidden">
          <label>Query XYZ<input type="file" id="qXyzFile"><textarea id="qXyz"></textarea></label>
          <label>Query WBO<textarea id="qWbo"></textarea></label>
        </div>
        <label id="qJsonWrap" class="hidden">Query WeightedGraph JSON<textarea id="qJson" class="tall"></textarea></label>
      </div>

      <div class="panel">
        <div class="row" style="justify-content:space-between">
          <b>Target</b>
          <select id="tType">
            <option value="drawn">Drawn</option>
            <option value="smiles">SMILES</option>
            <option value="xyz">XYZ</option>
            <option value="json">Weighted JSON</option>
          </select>
        </div>
        <div id="tDrawControls" class="section">
          <div class="row">
            <label>Mode <select id="tDrawMode"><option value="atom">Atom</option><option value="bond">Bond</option><option value="delete">Delete</option></select></label>
            <label>Element <input id="tElement" value="C" style="width:70px"></label>
            <label>Bond <select id="tBond"><option value="1">1</option><option value="1.5">1.5</option><option value="2">2</option><option value="3">3</option></select></label>
            <button data-clear-draw="t">Clear</button>
          </div>
          <svg id="tDraw" class="drawer"></svg>
        </div>
        <label id="tSmilesWrap" class="hidden">Target SMILES<textarea id="tSmiles"></textarea></label>
        <div id="tXyzWrap" class="hidden">
          <label>Target XYZ<input type="file" id="tXyzFile"><textarea id="tXyz"></textarea></label>
          <label>Target WBO<textarea id="tWbo"></textarea></label>
        </div>
        <label id="tJsonWrap" class="hidden">Target WeightedGraph JSON<textarea id="tJson" class="tall"></textarea></label>
      </div>

      <div class="panel">
        <div class="section">
          <div class="twocol">
            <label>Node policy <input id="nodePolicy" placeholder="element"></label>
            <label>iso_tol <input id="subIsoTol" type="number" value="1.0" step="0.1"></label>
            <label>graph_floor <input id="subGraphFloor" type="number" value="0.2" step="0.1"></label>
          </div>
          <button id="runSubBtn" class="primary">Run Subgraph Match</button>
          <span id="subSummary" class="metric">No match yet</span>
        </div>
        <pre id="subJson"></pre>
      </div>
    </div>
  </main>

  <script>
    const statusEl = document.getElementById("status");
    const rp = {preview:null, anchors:new Map(), selected:{R:null,P:null}, viewers:{R:null,P:null}};
    let currentSessionId = null;
    let lastRun = null;
    let selectedRpMechId = null;
    let lastSubResult = null;
    const palette = ["#e11d48","#2563eb","#16a34a","#d97706","#7c3aed","#0891b2","#be123c","#4d7c0f","#9333ea","#0f766e"];
    const drawState = {
      q:{nodes:[],edges:[],selected:null},
      t:{nodes:[],edges:[],selected:null}
    };

    function setStatus(text, cls=""){ statusEl.textContent=text; statusEl.className="metric "+cls; }
    function $(id){ return document.getElementById(id); }
    async function postJSON(url, payload){
      const res = await fetch(url,{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(payload)});
      const data = await res.json();
      if(!res.ok) throw new Error(data.error || "request failed");
      return data;
    }
    async function getJSON(url){
      const res = await fetch(url);
      const data = await res.json();
      if(!res.ok) throw new Error(data.error || "request failed");
      return data;
    }
    async function deleteJSON(url){
      const res = await fetch(url,{method:"DELETE"});
      const data = await res.json();
      if(!res.ok) throw new Error(data.error || "request failed");
      return data;
    }
    function readFileInto(inputId, textId){
      $(inputId).addEventListener("change", async e=>{
        const file=e.target.files && e.target.files[0]; if(!file) return;
        $(textId).value = await file.text();
      });
    }
    ["rXyzFile","pXyzFile","qXyzFile","tXyzFile"].forEach((id,i)=>{
      const targets=["rXyz","pXyz","qXyz","tXyz"]; readFileInto(id, targets[i]);
    });
    document.querySelectorAll(".tab").forEach(btn=>btn.onclick=()=>{
      document.querySelectorAll(".tab").forEach(b=>b.classList.remove("active"));
      document.querySelectorAll(".page").forEach(p=>p.classList.remove("active"));
      btn.classList.add("active"); $(btn.dataset.page).classList.add("active");
    });
    $("rpMode").onchange=()=>{
      const isXyz=$("rpMode").value==="xyz";
      $("smilesInputs").classList.toggle("hidden",isXyz);
      $("xyzInputs").classList.toggle("hidden",!isXyz);
    };
    $("xyzWboMode").onchange=()=> $("wboInputs").classList.toggle("hidden",$("xyzWboMode").value!=="provided");
    $("rpRenderMode").onchange=()=>{ renderRpOutput(); };
    $("show2dDegeneracy").onchange=()=>{ renderRp2DResult(); };

    function rpPayload(){
      const mode=$("rpMode").value;
      const payload={mode,name:$("rpName").value,anchors:Object.fromEntries(rp.anchors),config:{
        iso_tol:Number($("isoTol").value),dwbo_threshold:Number($("dwboThreshold").value),
        symmetry_wbo_tol:Number($("isoTol").value),n_seeds:Number($("nSeeds").value)
      }};
      if(mode==="smiles"){
        payload.reactant={smiles:$("rSmiles").value};
        payload.product={smiles:$("pSmiles").value};
      } else {
        payload.xyz_wbo_mode="xtb";
        payload.charge=Number($("charge").value); payload.multiplicity=Number($("multiplicity").value);
        payload.reactant={xyz:$("rXyz").value};
        payload.product={xyz:$("pXyz").value};
      }
      return payload;
    }
    $("previewBtn").onclick=async()=>{
      try{ setStatus("previewing");
        rp.preview=await postJSON("/api/preview",rpPayload());
        rp.selected={R:null,P:null}; renderRpPreview(); setStatus("preview ready","ok");
      }catch(e){ setStatus(e.message,"error"); }
    };
    $("runRpBtn").onclick=async()=>{
      try{ setStatus("running R-P AAM");
        const data=await postJSON("/api/rp",rpPayload());
        lastRun={job_id:data.job_id,viewer_url:data.viewer_url,stage_url:data.stage_url,viewer_rel_path:data.viewer_rel_path,stage_rel_path:data.stage_rel_path,name:data.name,reactant:data.reactant,product:data.product,mechanisms:data.rp.mechanisms || []};
        selectedRpMechId=(lastRun.mechanisms[0] && lastRun.mechanisms[0].id) || null;
        const mechs=data.rp.mechanisms || [];
        $("rpSummary").textContent=`mechanisms=${mechs.length} job=${data.job_id}`;
        $("rpJson").textContent=JSON.stringify(mechs.map(m=>({id:m.id,label:m.label,mapping_RP:m.mapping_RP,broken_bonds_R:m.broken_bonds_R,formed_bonds_R:m.formed_bonds_R,core_atoms:m.core_atoms,dedup_count:m.dedup_count})),null,2);
        renderRpOutput();
        setStatus("R-P complete","ok");
      }catch(e){ setStatus(e.message,"error"); }
    };
    $("clearAnchorsBtn").onclick=()=>{ rp.anchors.clear(); rp.selected={R:null,P:null}; renderAnchors(); renderRpPreview(); };

    function atomText(mol,i){ return (mol.elements[i]||"X")+i; }
    function anchorColor(r){ const keys=[...rp.anchors.keys()].sort((a,b)=>a-b); return palette[Math.max(0,keys.indexOf(r))%palette.length]; }
    function selectAtom(side,i){
      rp.selected[side]=i;
      if(rp.selected.R!==null && rp.selected.P!==null){
        for(const [r,p] of [...rp.anchors.entries()]) if(p===rp.selected.P && r!==rp.selected.R) rp.anchors.delete(r);
        rp.anchors.set(rp.selected.R,rp.selected.P);
      }
      renderAnchors(); renderRpPreview();
    }
    function renderAnchors(){
      $("selR").textContent=rp.selected.R===null?"-":String(rp.selected.R);
      $("selP").textContent=rp.selected.P===null?"-":String(rp.selected.P);
      const rows=[...rp.anchors.entries()].sort((a,b)=>a[0]-b[0]).map(([r,p])=>{
        return `<tr><td><span class="swatch" style="background:${anchorColor(r)}"></span>R${r}</td><td>P${p}</td><td><button data-del-anchor="${r}">x</button></td></tr>`;
      }).join("");
      $("anchorRows").innerHTML=rows || `<tr><td colspan="3" class="small">none</td></tr>`;
      document.querySelectorAll("[data-del-anchor]").forEach(b=>b.onclick=()=>{rp.anchors.delete(Number(b.dataset.delAnchor));renderAnchors();renderRpPreview();});
      $("anchorJson").textContent=JSON.stringify(Object.fromEntries(rp.anchors),null,2);
    }
    function clearBox(id,title){ $(id).classList.remove("mol2d"); $(id).innerHTML=`<div class="moltitle">${title}</div><div class="empty">Preview ${title}</div>`; }
    function renderRpPreview(){
      if(!rp.preview){ clearBox("rpRBox","R"); clearBox("rpPBox","P"); return; }
      const use2d=rp.preview.mode==="smiles";
      renderMolBox("rpRBox","R",rp.preview.reactant,use2d,(i)=>selectAtom("R",i));
      renderMolBox("rpPBox","P",rp.preview.product,use2d,(i)=>selectAtom("P",i));
    }
    function molColor(side,i){
      if(rp.selected[side]===i) return "#f1c40f";
      if(side==="R" && rp.anchors.has(i)) return anchorColor(i);
      if(side==="P"){ for(const [r,p] of rp.anchors.entries()) if(p===i) return anchorColor(r); }
      return null;
    }
    function renderMolBox(id,title,mol,use2d,onPick){
      const box=$(id); box.innerHTML=`<div class="moltitle">${title} atoms=${mol.elements.length}</div>`;
      box.classList.toggle("mol2d", !!use2d);
      if(use2d || !window.$3Dmol) render2DMol(box,title,mol,onPick);
      else render3DMol(id,title,mol,onPick);
    }
    function render2DMol(box,title,mol,onPick){
      box.classList.add("mol2d");
      const w=box.clientWidth||420,h=box.clientHeight||315,pad=34;
      const pts=mol.coords.map(p=>[Number(p[0]),Number(p[1])]);
      const xs=pts.map(p=>p[0]), ys=pts.map(p=>p[1]);
      const minX=Math.min(...xs,0), maxX=Math.max(...xs,1), minY=Math.min(...ys,0), maxY=Math.max(...ys,1);
      const sx=(w-2*pad)/Math.max(1e-6,maxX-minX), sy=(h-2*pad)/Math.max(1e-6,maxY-minY), s=Math.min(sx,sy);
      const map=p=>[pad+(p[0]-minX)*s, h-pad-(p[1]-minY)*s];
      const lines=(mol.bonds||[]).map(b=>{const a=map(pts[b.i]),c=map(pts[b.j]);return `<line x1="${a[0]}" y1="${a[1]}" x2="${c[0]}" y2="${c[1]}" stroke="#39434d" stroke-width="${Math.max(2,Number(b.wbo)||1)}"/>`;}).join("");
      const atoms=pts.map((p,i)=>{const q=map(p), color=molColor(title,i)||"#fff"; const stroke=color==="#fff"?"#2f3b46":color; return `<g data-atom="${i}" style="cursor:pointer"><circle cx="${q[0]}" cy="${q[1]}" r="15" fill="${color}" fill-opacity="${color==="#fff"?.92:.52}" stroke="${stroke}" stroke-width="2"/><text x="${q[0]}" y="${q[1]+4}" text-anchor="middle" font-size="12" font-weight="650">${mol.elements[i]}${i}</text></g>`;}).join("");
      box.innerHTML += `<svg viewBox="0 0 ${w} ${h}">${lines}${atoms}</svg>`;
      box.querySelectorAll("[data-atom]").forEach(g=>g.onclick=()=>onPick(Number(g.dataset.atom)));
    }
    function xyzText(mol){ let s=mol.elements.length+"\\nframe\\n"; mol.elements.forEach((e,i)=>{const p=mol.coords[i]; s+=`${e} ${p[0]} ${p[1]} ${p[2]}\\n`;}); return s; }
    function render3DMol(id,title,mol,onPick){
      const box=$(id); const div=document.createElement("div"); div.style.position="absolute"; div.style.inset="0"; box.appendChild(div);
      box.classList.remove("mol2d");
      const viewer=$3Dmol.createViewer(div,{backgroundColor:"white"}); rp.viewers[title]=viewer;
      viewer.addModel(xyzText(mol),"xyz"); viewer.setStyle({}, {stick:{radius:.1},sphere:{scale:.22}});
      viewer.setClickable({}, true, atom=>{ const idx=Number.isInteger(atom.index)?atom.index:atom.serial-1; onPick(idx); });
      for(let i=0;i<mol.elements.length;i++){
        const c=molColor(title,i); const p=mol.coords[i]; if(c){ viewer.addSphere({center:{x:p[0],y:p[1],z:p[2]},radius:.36,color:c,alpha:.48}); }
        viewer.addLabel(String(i),{position:{x:p[0],y:p[1],z:p[2]},fontSize:9,fontColor:"black",backgroundColor:"white",backgroundOpacity:.7,inFront:true});
      }
      viewer.zoomTo(); viewer.render();
    }
    function renderRpOutput(){
      const mode=$("rpRenderMode").value;
      const hasRun=!!(lastRun && lastRun.viewer_url);
      if(mode==="3d"){
        $("rp2dResult").classList.add("hidden");
        if(hasRun){ $("rpViewer").src=lastRun.viewer_url; $("rpViewer").classList.remove("hidden"); }
        else { $("rpViewer").classList.add("hidden"); $("rpViewer").src=""; }
      } else {
        $("rpViewer").classList.add("hidden");
        $("rp2dResult").classList.remove("hidden");
        renderRp2DResult();
      }
      $("stageLink").href=hasRun && lastRun.stage_url ? lastRun.stage_url : "";
      $("stageLink").textContent=hasRun && lastRun.stage_url ? "rp_stage.json" : "";
    }
    function pairKey(pair){ return pair.slice().sort((a,b)=>a-b).join("-"); }
    function mechanismById(){
      if(!lastRun || !lastRun.mechanisms) return null;
      return lastRun.mechanisms.find(m=>String(m.id)===String(selectedRpMechId)) || lastRun.mechanisms[0] || null;
    }
    function bondsWithEvents(mol, eventPairs, kind){
      const events=new Set((eventPairs || []).map(p=>pairKey(p.map(Number))));
      return (mol.bonds || []).map(b=>({...b,event:events.has(pairKey([b.i,b.j])) ? kind : null}));
    }
    function branchColorGroups(mech){
      const sym=mech && mech.branch_symmetry ? mech.branch_symmetry : null;
      if(!sym) return [];
      return (sym.color_groups && sym.color_groups.length) ? sym.color_groups : (sym.blocks || []);
    }
    function branchRDegMap(mech){
      const out={};
      branchColorGroups(mech).forEach((group,idx)=>{
        (group.r_atoms || []).forEach(r=>{ out[Number(r)] = idx; });
      });
      return out;
    }
    function branchProductDegMap(mech){
      const out={};
      branchColorGroups(mech).forEach((group,idx)=>{
        (group.p_atoms || []).forEach(p=>{ out[Number(p)] = idx; });
      });
      return out;
    }
    function degMapFor(side, mech){
      if(!$("show2dDegeneracy").checked) return null;
      return side==="R" ? branchRDegMap(mech) : branchProductDegMap(mech);
    }
    function renderRp2DResult(){
      if(!lastRun || !lastRun.reactant || !lastRun.product){
        clearBox("rp2dR","R 2D"); clearBox("rp2dP","P 2D"); $("rp2dMechs").innerHTML=""; return;
      }
      const mechs=lastRun.mechanisms || [];
      if(selectedRpMechId===null && mechs[0]) selectedRpMechId=mechs[0].id;
      $("rp2dMechs").innerHTML=mechs.map(m=>`<button data-rp-mech="${m.id}" class="${String(m.id)===String(selectedRpMechId)?"primary":""}">#${m.id}</button>`).join("");
      document.querySelectorAll("[data-rp-mech]").forEach(b=>b.onclick=()=>{ selectedRpMechId=b.dataset.rpMech; renderRp2DResult(); });
      const mech=mechanismById() || {};
      const rMol={...lastRun.reactant,bonds:bondsWithEvents(lastRun.reactant,mech.broken_bonds_R || [],"broken")};
      const pFormed=mech.formed_bonds_P || [];
      const pMol={...lastRun.product,bonds:bondsWithEvents(lastRun.product,pFormed,"formed")};
      render2DStaticBox("rp2dR","R 2D",rMol,mech.core_atoms || [],degMapFor("R",mech));
      render2DStaticBox("rp2dP","P 2D",pMol,[],degMapFor("P",mech));
    }
    function render2DStaticBox(id,title,mol,coreAtoms=[],degMap=null){
      const degCount=degMap ? new Set(Object.values(degMap)).size : 0;
      const suffix=degCount ? ` deg=${degCount}` : "";
      const box=$(id); box.innerHTML=`<div class="moltitle">${title} atoms=${mol.elements.length}${suffix}</div>`;
      box.classList.add("mol2d");
      const w=box.clientWidth||420,h=box.clientHeight||260,pad=34;
      const pts=mol.coords.map(p=>[Number(p[0]),Number(p[1])]);
      const xs=pts.map(p=>p[0]), ys=pts.map(p=>p[1]);
      const minX=Math.min(...xs,0), maxX=Math.max(...xs,1), minY=Math.min(...ys,0), maxY=Math.max(...ys,1);
      const sx=(w-2*pad)/Math.max(1e-6,maxX-minX), sy=(h-2*pad)/Math.max(1e-6,maxY-minY), s=Math.min(sx,sy);
      const map=p=>[pad+(p[0]-minX)*s, h-pad-(p[1]-minY)*s];
      const lines=(mol.bonds||[]).map(b=>{const a=map(pts[b.i]),c=map(pts[b.j]); const color=b.event==="broken"?"#c62828":(b.event==="formed"?"#16803c":"#39434d"); const dash=b.event?'stroke-dasharray="7 4"':""; return `<line x1="${a[0]}" y1="${a[1]}" x2="${c[0]}" y2="${c[1]}" stroke="${color}" stroke-width="${b.event?4:Math.max(2,Number(b.wbo)||1)}" ${dash}/>`;}).join("");
      const core=new Set((coreAtoms||[]).map(Number));
      const atoms=pts.map((p,i)=>{const q=map(p), rawGroup=degMap && Object.prototype.hasOwnProperty.call(degMap,i) ? Number(degMap[i]) : null; const hasDeg=Number.isFinite(rawGroup); const color=hasDeg ? palette[Math.abs(rawGroup)%palette.length] : null; const fill=hasDeg ? color : (core.has(i)?"#fff1b8":"#fff"); const fillOpacity=hasDeg ? ".22" : "1"; const stroke=hasDeg ? color : (core.has(i)?"#b88a00":"#2f3b46"); const ring=core.has(i)&&hasDeg ? `<circle cx="${q[0]}" cy="${q[1]}" r="18" fill="none" stroke="#b88a00" stroke-width="2"/>` : ""; return `<g>${ring}<circle cx="${q[0]}" cy="${q[1]}" r="15" fill="${fill}" fill-opacity="${fillOpacity}" stroke="${stroke}" stroke-width="${hasDeg?3:2}"/><text x="${q[0]}" y="${q[1]+4}" text-anchor="middle" font-size="12" font-weight="650">${mol.elements[i]}${i}</text></g>`;}).join("");
      box.innerHTML += `<svg viewBox="0 0 ${w} ${h}">${lines}${atoms}</svg>`;
    }

    function initDraw(which){
      const svg=$(which+"Draw");
      svg.addEventListener("click",e=>{
        const rect=svg.getBoundingClientRect(); const x=e.clientX-rect.left, y=e.clientY-rect.top;
        const hit=findDrawAtom(which,x,y); const mode=$(which+"DrawMode").value;
        if(mode==="atom" && hit===null){ drawState[which].nodes.push({element:$(which+"Element").value||"C",x,y}); }
        else if(mode==="delete" && hit!==null){ deleteDrawAtom(which,hit); }
        else if(mode==="bond" && hit!==null){ handleDrawBond(which,hit); }
        renderDraw(which);
      });
      renderDraw(which);
    }
    function findDrawAtom(which,x,y){ const ns=drawState[which].nodes; for(let i=ns.length-1;i>=0;i--){ const dx=ns[i].x-x,dy=ns[i].y-y; if(Math.hypot(dx,dy)<=18) return i; } return null; }
    function deleteDrawAtom(which,i){ const st=drawState[which]; st.nodes.splice(i,1); st.edges=st.edges.filter(e=>e.i!==i&&e.j!==i).map(e=>({i:e.i>i?e.i-1:e.i,j:e.j>i?e.j-1:e.j,order:e.order})); st.selected=null; }
    function handleDrawBond(which,i){ const st=drawState[which]; if(st.selected===null){ st.selected=i; return; } const a=Math.min(st.selected,i), b=Math.max(st.selected,i); if(a!==b){ const old=st.edges.find(e=>e.i===a&&e.j===b); if(old) old.order=Number($(which+"Bond").value); else st.edges.push({i:a,j:b,order:Number($(which+"Bond").value)}); } st.selected=null; }
    function renderDraw(which){
      const svg=$(which+"Draw"), st=drawState[which]; const w=svg.clientWidth||330,h=svg.clientHeight||330;
      const edges=st.edges.map(e=>{const a=st.nodes[e.i],b=st.nodes[e.j];return `<line x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}" stroke="#39434d" stroke-width="${Math.max(2,e.order)}"/><text x="${(a.x+b.x)/2+4}" y="${(a.y+b.y)/2-4}" font-size="11">${e.order}</text>`;}).join("");
      const nodes=st.nodes.map((n,i)=>`<g><circle cx="${n.x}" cy="${n.y}" r="16" fill="${st.selected===i?"#ffe08a":"#fff"}" stroke="#24313d" stroke-width="2"/><text x="${n.x}" y="${n.y+4}" text-anchor="middle" font-size="12" font-weight="650">${n.element}${i}</text></g>`).join("");
      svg.setAttribute("viewBox",`0 0 ${w} ${h}`); svg.innerHTML=edges+nodes;
    }
    ["q","t"].forEach(initDraw);
    document.querySelectorAll("[data-clear-draw]").forEach(b=>b.onclick=()=>{const w=b.dataset.clearDraw; drawState[w]={nodes:[],edges:[],selected:null}; renderDraw(w);});
    function switchSource(prefix){
      const type=$(prefix+"Type").value;
      [prefix+"DrawControls",prefix+"SmilesWrap",prefix+"XyzWrap",prefix+"JsonWrap"].forEach(id=>$(id).classList.add("hidden"));
      const map={drawn:"DrawControls",smiles:"SmilesWrap",xyz:"XyzWrap",json:"JsonWrap"};
      $(prefix+map[type]).classList.remove("hidden");
    }
    ["qType","tType"].forEach(id=>$(id).onchange=()=>switchSource(id[0]));
    switchSource("q"); switchSource("t");
    function subSpec(prefix){
      const type=$(prefix+"Type").value;
      if(type==="drawn") return {type,graph:drawState[prefix]};
      if(type==="smiles") return {type,smiles:$(prefix+"Smiles").value};
      if(type==="xyz") return {type,xyz:$(prefix+"Xyz").value,wbo:$(prefix+"Wbo").value};
      return {type,graph_json:$(prefix+"Json").value};
    }
    $("runSubBtn").onclick=async()=>{
      try{ setStatus("running subgraph");
        const policy=$("nodePolicy").value.trim();
        const data=await postJSON("/api/subgraph",{query:subSpec("q"),target:subSpec("t"),node_policy:policy||null,iso_tol:Number($("subIsoTol").value),graph_floor:Number($("subGraphFloor").value)});
        lastSubResult=data;
        $("subSummary").textContent=`matches=${data.n_matches}`;
        $("subJson").textContent=JSON.stringify(data,null,2);
        setStatus("subgraph complete","ok");
      }catch(e){ setStatus(e.message,"error"); }
    };

    function activePageId(){
      const active=document.querySelector(".page.active");
      return active ? active.id : "rpPage";
    }
    function setActivePage(pageId){
      document.querySelectorAll(".tab").forEach(b=>b.classList.toggle("active", b.dataset.page===pageId));
      document.querySelectorAll(".page").forEach(p=>p.classList.toggle("active", p.id===pageId));
    }
    function collectSessionState(){
      return {
        version:1,
        session_id:currentSessionId,
        name:$("sessionName").value || $("rpName").value || "aam_session",
        activePage:activePageId(),
        rp:{
          mode:$("rpMode").value,
          name:$("rpName").value,
          render_mode:$("rpRenderMode").value,
          show_2d_degeneracy:$("show2dDegeneracy").checked,
          selected_mech_id:selectedRpMechId,
          reactant_smiles:$("rSmiles").value,
          product_smiles:$("pSmiles").value,
          xyz_wbo_mode:"xtb",
          charge:$("charge").value,
          multiplicity:$("multiplicity").value,
          reactant_xyz:$("rXyz").value,
          product_xyz:$("pXyz").value,
          config:{iso_tol:$("isoTol").value,dwbo_threshold:$("dwboThreshold").value,n_seeds:$("nSeeds").value},
          anchors:Object.fromEntries(rp.anchors),
          selected:rp.selected,
          preview:rp.preview,
          lastRun,
          rpSummary:$("rpSummary").textContent,
          rpJson:$("rpJson").textContent
        },
        subgraph:{
          qType:$("qType").value,tType:$("tType").value,
          qDraw:drawState.q,tDraw:drawState.t,
          qDrawMode:$("qDrawMode").value,tDrawMode:$("tDrawMode").value,
          qElement:$("qElement").value,tElement:$("tElement").value,
          qBond:$("qBond").value,tBond:$("tBond").value,
          qSmiles:$("qSmiles").value,tSmiles:$("tSmiles").value,
          qXyz:$("qXyz").value,tXyz:$("tXyz").value,
          qWbo:$("qWbo").value,tWbo:$("tWbo").value,
          qJson:$("qJson").value,tJson:$("tJson").value,
          nodePolicy:$("nodePolicy").value,
          iso_tol:$("subIsoTol").value,
          graph_floor:$("subGraphFloor").value,
          result:lastSubResult,
          summary:$("subSummary").textContent,
          resultJson:$("subJson").textContent
        }
      };
    }
    function applySessionState(state, restored=null){
      state=state || {};
      currentSessionId=state.session_id || currentSessionId;
      $("sessionName").value=state.name || "";
      if(state.activePage) setActivePage(state.activePage);
      const rs=state.rp || {};
      if(rs.mode) $("rpMode").value=rs.mode;
      $("rpRenderMode").value=rs.render_mode || "2d";
      $("show2dDegeneracy").checked=!!rs.show_2d_degeneracy;
      selectedRpMechId=rs.selected_mech_id ?? null;
      $("rpName").value=rs.name || "aam_alignment";
      $("rSmiles").value=rs.reactant_smiles || "";
      $("pSmiles").value=rs.product_smiles || "";
      $("xyzWboMode").value="xtb";
      $("charge").value=rs.charge ?? "0";
      $("multiplicity").value=rs.multiplicity ?? "1";
      $("rXyz").value=rs.reactant_xyz || "";
      $("pXyz").value=rs.product_xyz || "";
      $("rWbo").value="";
      $("pWbo").value="";
      const cfg=rs.config || {};
      $("isoTol").value=cfg.iso_tol ?? "1.0";
      $("dwboThreshold").value=cfg.dwbo_threshold ?? "0.5";
      $("nSeeds").value=cfg.n_seeds ?? "3";
      $("rpMode").onchange(); $("xyzWboMode").onchange();
      rp.anchors=new Map(Object.entries(rs.anchors || {}).map(([k,v])=>[Number(k),Number(v)]));
      rp.selected=rs.selected || {R:null,P:null};
      rp.preview=rs.preview || null;
      lastRun=rs.lastRun || null;
      if(restored && restored.viewer_url){
        lastRun={...(lastRun || {}),job_id:restored.job_id,viewer_url:restored.viewer_url,stage_url:restored.stage_url};
      }
      if(lastRun && lastRun.viewer_url){
        $("stageLink").href=lastRun.stage_url || ""; $("stageLink").textContent=lastRun.stage_url ? "rp_stage.json" : "";
      } else {
        $("stageLink").href=""; $("stageLink").textContent="";
      }
      $("rpSummary").textContent=rs.rpSummary || "No run yet";
      $("rpJson").textContent=rs.rpJson || "";
      renderAnchors(); renderRpPreview(); renderRpOutput();

      const ss=state.subgraph || {};
      $("qType").value=ss.qType || "drawn"; $("tType").value=ss.tType || "drawn";
      drawState.q=ss.qDraw || {nodes:[],edges:[],selected:null};
      drawState.t=ss.tDraw || {nodes:[],edges:[],selected:null};
      $("qDrawMode").value=ss.qDrawMode || "atom"; $("tDrawMode").value=ss.tDrawMode || "atom";
      $("qElement").value=ss.qElement || "C"; $("tElement").value=ss.tElement || "C";
      $("qBond").value=ss.qBond || "1"; $("tBond").value=ss.tBond || "1";
      $("qSmiles").value=ss.qSmiles || ""; $("tSmiles").value=ss.tSmiles || "";
      $("qXyz").value=ss.qXyz || ""; $("tXyz").value=ss.tXyz || "";
      $("qWbo").value=ss.qWbo || ""; $("tWbo").value=ss.tWbo || "";
      $("qJson").value=ss.qJson || ""; $("tJson").value=ss.tJson || "";
      $("nodePolicy").value=ss.nodePolicy || "";
      $("subIsoTol").value=ss.iso_tol ?? "1.0";
      $("subGraphFloor").value=ss.graph_floor ?? "0.2";
      lastSubResult=ss.result || null;
      $("subSummary").textContent=ss.summary || "No match yet";
      $("subJson").textContent=ss.resultJson || (lastSubResult ? JSON.stringify(lastSubResult,null,2) : "");
      switchSource("q"); switchSource("t"); renderDraw("q"); renderDraw("t");
    }
    async function refreshSessions(){
      try{
        const data=await getJSON("/api/sessions");
        const sel=$("sessionSelect");
        if(!data.enabled){
          sel.innerHTML='<option value="">Mongo unavailable</option>';
          setStatus(data.error || "session store unavailable","error");
          return;
        }
        sel.innerHTML=(data.sessions || []).map(s=>`<option value="${s.id}">${s.name || s.id} (${s.artifact_count || 0})</option>`).join("") || '<option value="">No saved sessions</option>';
        if(currentSessionId) sel.value=currentSessionId;
      }catch(e){ setStatus(e.message,"error"); }
    }
    $("saveSessionBtn").onclick=async()=>{
      try{
        const state=collectSessionState();
        const data=await postJSON("/api/sessions",{id:currentSessionId,name:state.name,state,job_id:lastRun && lastRun.job_id});
        currentSessionId=data.session.id;
        $("sessionName").value=data.session.name || state.name;
        await refreshSessions();
        setStatus("session saved","ok");
      }catch(e){ setStatus(e.message,"error"); }
    };
    $("refreshSessionsBtn").onclick=refreshSessions;
    $("loadSessionBtn").onclick=async()=>{
      const id=$("sessionSelect").value; if(!id) return;
      try{
        const data=await getJSON("/api/sessions/"+encodeURIComponent(id));
        currentSessionId=data.id;
        const state=data.state || {};
        state.session_id=data.id;
        applySessionState(state,data.restored);
        setStatus("session loaded","ok");
      }catch(e){ setStatus(e.message,"error"); }
    };
    $("deleteSessionBtn").onclick=async()=>{
      const id=$("sessionSelect").value; if(!id) return;
      try{
        await deleteJSON("/api/sessions/"+encodeURIComponent(id));
        if(currentSessionId===id) currentSessionId=null;
        await refreshSessions();
        setStatus("session deleted","ok");
      }catch(e){ setStatus(e.message,"error"); }
    };
    getJSON("/api/session-status").then(data=>{
      if(data.enabled) refreshSessions();
      else if(data.error) setStatus("Mongo unavailable: "+data.error,"error");
    }).catch(()=>{});
    renderAnchors();
    renderRpOutput();
  </script>
</body>
</html>
"""


def run_server(host="127.0.0.1", port=8765, work_root=None,
               open_browser=False, mongo_uri=None, mongo_db="rxn_core"):
    root = Path(work_root or (Path(tempfile.gettempdir()) / "rxn_core_aam_webapp"))
    server = AAMHTTPServer(
        (host, int(port)),
        AAMWebApp(root, mongo_uri=mongo_uri, mongo_db=mongo_db))
    url = f"http://{host}:{server.server_port}/"
    if open_browser:
        webbrowser.open(url)
    print(f"rxn_core AAM webapp: {url}")
    print(f"artifact root: {root}")
    if server.app.session_store is not None:
        print(f"mongo sessions: enabled db={mongo_db}")
    elif server.app.session_error:
        print(f"mongo sessions: unavailable ({server.app.session_error})")
    else:
        print("mongo sessions: disabled")
    try:
        server.serve_forever()
    finally:
        server.server_close()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--work-root", default=None)
    parser.add_argument(
        "--mongo-uri",
        default=os.environ.get("RXN_CORE_MONGO_URI", "mongodb://localhost:27017"))
    parser.add_argument(
        "--mongo-db",
        default=os.environ.get("RXN_CORE_MONGO_DB", "rxn_core"))
    parser.add_argument("--no-mongo", action="store_true")
    parser.add_argument("--open", action="store_true")
    args = parser.parse_args(argv)
    run_server(
        args.host, args.port, args.work_root, args.open,
        mongo_uri=None if args.no_mongo else args.mongo_uri,
        mongo_db=args.mongo_db)


if __name__ == "__main__":
    main()

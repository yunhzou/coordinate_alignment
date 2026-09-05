"""Serialization and self-contained views for typed computational results."""
from __future__ import annotations

import html
import json
from pathlib import Path

import numpy as np

from .chemistry_computations import write_xyz_str
from .alignment.post_aam import AtomBijection
from .domain import (
    AAMProblem, AAMSearchConfig, AAMResult, AAMSearchMetrics, MolecularEndpoint, RPResult, ReactionContext,
    ResolvedMechanism, TSResult,
)
from .search_graph import AAMSearchGraph


def aam_record(result: AAMResult):
    """Persist the search graph before selecting mechanisms or representatives."""
    from dataclasses import asdict
    def endpoint(molecule):
        return {'elements': molecule.elements, 'coordinates': molecule.coordinates.tolist(),
                'wbo': molecule.wbo.tolist(), 'label': molecule.label,
                'energy': molecule.energy, 'metadata': dict(molecule.metadata)}
    return {'schema': 'rxn_core.aam/v1', 'name': result.problem.name,
            'reactant': endpoint(result.problem.reactant),
            'product': endpoint(result.problem.product),
            'config': asdict(result.config), 'metrics': asdict(result.metrics),
            'graph': result.graph.to_record()}


def aam_from_record(record):
    if record['schema'] != 'rxn_core.aam/v1':
        raise ValueError('unsupported AAM result schema')
    def endpoint(raw):
        return MolecularEndpoint(tuple(raw['elements']), raw['coordinates'], raw['wbo'], raw['label'],
                                  energy=raw.get('energy'), metadata=raw.get('metadata', {}))
    return AAMResult(AAMProblem(endpoint(record['reactant']), endpoint(record['product']), record['name']),
                     AAMSearchConfig(**record['config']), AAMSearchGraph.from_record(record['graph']),
                     AAMSearchMetrics(**record['metrics']))


def write_aam_bundle(result: AAMResult, output_directory):
    """Save raw AAM and an offline graph/path viewer; never rerun matching."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    record = aam_record(result)
    _json_dump(output / 'aam.json', record)
    assets = Path(__file__).parent / 'static'
    page = (assets / 'aam_search.html').read_text()
    page = page.replace('__LIBRARY__', (assets / '3Dmol-min.js').read_text())
    page = page.replace('__DATA__', json.dumps(record).replace('<', '\\u003c'))
    (output / 'search.html').write_text(page)
    return output


def _mapping_record(mapping):
    return {str(source): int(target)
            for source, target in mapping.as_dict().items()}


def rp_record(result: RPResult):
    """Return the stable JSON boundary representation of an ``RPResult``."""
    problem = result.analytical.aam.problem
    return {
        "schema": "rxn_core.rp/v2",
        "aam": aam_record(result.analytical.aam),
        "name": problem.name,
        "atom_count": problem.atom_count,
        "timing": {
            "aam_seconds": result.analytical.aam.metrics.elapsed_seconds,
            "family_seconds": result.analytical.elapsed_seconds,
            "selection_seconds": result.elapsed_seconds,
        },
        "search_metrics": vars(result.analytical.aam.metrics),
        "mechanisms": [{
            "id": index,
            "mapping": _mapping_record(mechanism.mapping),
            "broken_bonds": [list(bond) for bond in mechanism.broken_bonds],
            "formed_bonds": [list(bond) for bond in mechanism.formed_bonds],
            "core_atoms": list(mechanism.core_atoms),
            "fixed_mapping_rmsd": mechanism.fixed_mapping_rmsd,
            "selected_branch_index": mechanism.selected_branch_index,
            "chirality": dict(mechanism.chirality),
            "analytical_family_count": len(mechanism.analytical.branches),
        } for index, mechanism in enumerate(result.mechanisms, 1)],
    }


def ts_record(result: TSResult):
    """Return a compact JSON boundary representation of one typed TS run."""
    return {
        "schema": "rxn_core.ts/v2",
        "target": result.mechanisms[0].target.molecule.label
        if result.mechanisms else "",
        "elapsed_seconds": result.elapsed_seconds,
        "mechanisms": [{
            "id": index,
            "status": item.status,
            "reason": item.reason,
            "reactant_core_assignments": (
                0 if item.reactant_core_aam is None
                else len(item.reactant_core_aam.assignments)),
            "product_core_assignments": (
                0 if item.product_core_aam is None
                else len(item.product_core_aam.assignments)),
            "candidate_count": len(item.candidates),
            "selected": None if item.selected is None else {
                "assignment": {str(a): int(b)
                               for a, b in item.selected.assignment.pairs},
                "sources": sorted(item.selected.sources),
                "score": item.selected.score,
                "overlap": item.selected.overlap,
                "wbo_progress": item.selected.wbo_progress,
                "mode_index": item.selected.mode_index,
                "frequency": item.selected.frequency,
                "event_terms": [dict(term)
                                for term in item.selected.event_terms],
            },
        } for index, item in enumerate(result.mechanisms, 1)],
    }


def reaction_record(reaction: ReactionContext):
    """Serialize the exact R/P information required by the TS stage."""
    return {
        "schema": "rxn_core.resolved_reaction/v1",
        "name": reaction.problem.name,
        "atom_count": reaction.problem.atom_count,
        "mechanisms": [{
            "mapping": _mapping_record(item.mapping),
            "broken_bonds": [list(bond) for bond in item.broken_bonds],
            "formed_bonds": [list(bond) for bond in item.formed_bonds],
            "core_atoms": list(item.core_atoms),
        } for item in reaction.mechanisms],
    }


def reaction_from_record(record, problem: AAMProblem,
                         config: AAMSearchConfig | None = None):
    """Materialize a TS-stage context without executing R/P search."""
    if record.get("schema") != "rxn_core.resolved_reaction/v1":
        raise ValueError("unsupported resolved reaction schema")
    if int(record.get("atom_count", -1)) != problem.atom_count:
        raise ValueError("resolved reaction atom count differs from endpoints")
    mechanisms = []
    for raw in record.get("mechanisms") or ():
        mechanisms.append(ResolvedMechanism(
            mapping=AtomBijection.from_mapping(
                {int(a): int(b) for a, b in raw["mapping"].items()},
                degree=problem.atom_count),
            broken_bonds=tuple(tuple(bond)
                               for bond in raw.get("broken_bonds") or ()),
            formed_bonds=tuple(tuple(bond)
                               for bond in raw.get("formed_bonds") or ()),
            core_atoms=tuple(raw.get("core_atoms") or ()),
        ))
    if not mechanisms:
        raise ValueError("resolved reaction contains no mechanisms")
    return ReactionContext(
        problem, config or AAMSearchConfig(), tuple(mechanisms))


def _json_dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=float))


def _viewer_html(result: RPResult):
    problem = result.analytical.aam.problem
    reactant = write_xyz_str(
        problem.reactant.elements, problem.reactant.coordinates, "Reactant")
    mechanisms = []
    for index, mechanism in enumerate(result.mechanisms, 1):
        product = mechanism.mapping.product_in_reactant_order(
            problem.product.coordinates)
        mechanisms.append({
            "id": index,
            "product": write_xyz_str(
                problem.reactant.elements, product, f"Mechanism {index} product"),
            "broken": mechanism.broken_bonds,
            "formed": mechanism.formed_bonds,
            "rmsd": mechanism.fixed_mapping_rmsd,
        })
    library = (Path(__file__).parent / "static" / "3Dmol-min.js").read_text()
    payload = json.dumps({"reactant": reactant, "mechanisms": mechanisms})
    title = html.escape(problem.name or "R/P alignment")
    return f"""<!doctype html><html><head><meta charset=\"utf-8\">
<title>{title}</title><style>
body{{font:14px system-ui;margin:0;background:#111;color:#eee}}
header{{padding:12px 18px;background:#222}} select{{margin-left:8px}}
#grid{{display:grid;grid-template-columns:1fr 1fr;height:calc(100vh - 52px)}}
.panel{{position:relative;border-top:1px solid #444}} .label{{position:absolute;z-index:2;padding:8px}}
.view{{position:absolute;inset:0}}</style><script>{library}</script></head>
<body><header>{title}<select id=\"mechanism\"></select><span id=\"meta\"></span></header>
<div id=\"grid\"><div class=\"panel\"><b class=\"label\">R</b><div id=\"r\" class=\"view\"></div></div>
<div class=\"panel\"><b class=\"label\">P aligned</b><div id=\"p\" class=\"view\"></div></div></div>
<script>const data={payload};
const rv=$3Dmol.createViewer('r',{{backgroundColor:'#111'}}),pv=$3Dmol.createViewer('p',{{backgroundColor:'#111'}});
function draw(v,xyz){{v.removeAllModels();v.addModel(xyz,'xyz');v.setStyle({{}},{{stick:{{radius:.15}},sphere:{{scale:.28}}}});v.zoomTo();v.render();}}
draw(rv,data.reactant);const sel=document.getElementById('mechanism');
data.mechanisms.forEach(m=>{{const o=document.createElement('option');o.value=m.id;o.textContent=' mechanism '+m.id;sel.appendChild(o);}});
function update(){{const m=data.mechanisms.find(x=>x.id==sel.value)||data.mechanisms[0];if(!m)return;draw(pv,m.product);document.getElementById('meta').textContent=' RMSD '+m.rmsd.toFixed(4)+' Å';}}
sel.onchange=update;update();window.onresize=()=>{{rv.resize();pv.resize();}};</script></body></html>"""


def write_rp_bundle(result: RPResult, output_directory):
    """Write JSON, per-mechanism endpoint XYZ files, and an inline viewer."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    problem = result.analytical.aam.problem
    _json_dump(output / "rp.json", rp_record(result))
    from .ts import reaction_context_from_rp
    _json_dump(
        output / "reaction.json",
        reaction_record(reaction_context_from_rp(result)))
    (output / "R.xyz").write_text(write_xyz_str(
        problem.reactant.elements, problem.reactant.coordinates, "Reactant"))
    for index, mechanism in enumerate(result.mechanisms, 1):
        directory = output / f"mechanism_{index:03d}"
        directory.mkdir(exist_ok=True)
        (directory / "R.xyz").write_text(write_xyz_str(
            problem.reactant.elements, problem.reactant.coordinates, "Reactant"))
        aligned = mechanism.mapping.product_in_reactant_order(
            problem.product.coordinates)
        (directory / "P_aligned.xyz").write_text(write_xyz_str(
            problem.reactant.elements, aligned, "Product in reactant order"))
    (output / "view.html").write_text(_viewer_html(result))
    return output


def write_ts_record(result: TSResult, path):
    _json_dump(path, ts_record(result))
    return Path(path)

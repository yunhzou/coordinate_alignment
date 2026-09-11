"""Serialization and self-contained views for typed computational results."""
from __future__ import annotations

import html
import json
import gc
import gzip
import pickle
from pathlib import Path

import numpy as np

from .chemistry_computations import write_xyz_str
from .alignment.post_aam import AtomBijection
from .domain import (
    AAMProblem, AAMSearchConfig, AAMResult, AAMSearchMetrics, MolecularEndpoint, RPResult, ReactionContext,
    ResolvedMechanism, TSResult,
)
from .search_graph import AAMSearchGraph


def _aam_header(result):
    from dataclasses import asdict
    def endpoint(molecule):
        return {'elements': molecule.elements, 'coordinates': molecule.coordinates.tolist(),
                'wbo': molecule.wbo.tolist(), 'label': molecule.label,
                'energy': molecule.energy, 'metadata': dict(molecule.metadata)}
    return {'schema': 'rxn_core.aam/v1', 'name': result.problem.name,
            'reactant': endpoint(result.problem.reactant),
            'product': endpoint(result.problem.product),
            'config': asdict(result.config), 'metrics': asdict(result.metrics)}


def aam_record(result: AAMResult, *, copy_graph=True):
    """Persist the search graph before selecting mechanisms or representatives."""
    return {**_aam_header(result),'graph': result.graph.to_record(copy=copy_graph)}


def aam_from_record(record, *, copy_graph=True):
    if record['schema'] != 'rxn_core.aam/v1':
        raise ValueError('unsupported AAM result schema')
    return _aam_result(record,AAMSearchGraph.from_record(record['graph'],copy=copy_graph))


def _aam_result(record,graph):
    def endpoint(raw):
        return MolecularEndpoint(tuple(raw['elements']), raw['coordinates'], raw['wbo'], raw['label'],
                                  energy=raw.get('energy'), metadata=raw.get('metadata', {}))
    return AAMResult(AAMProblem(endpoint(record['reactant']), endpoint(record['product']), record['name']),
                     AAMSearchConfig(**record['config']), graph,
                     AAMSearchMetrics(**record['metrics']))


def write_aam_checkpoint(result,path):
    """Trusted internal snapshot preserving shared Python graph objects.

    Keep the producing engine with this file. JSON remains the interchange
    format; this pickle-based checkpoint must never be loaded from strangers.
    """
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    record={**_aam_header(result),'schema':'rxn_core.aam_checkpoint/v1','graph':result.graph}
    enabled=gc.isenabled()
    try:
        gc.disable()
        with gzip.open(temporary,'wb',compresslevel=1) as stream:
            pickle.dump(record,stream,protocol=5)
        temporary.replace(path)
    finally:
        if enabled:gc.enable()


def read_aam_checkpoint(path):
    """Read only trusted internally produced snapshots, never untrusted pickle."""
    enabled=gc.isenabled()
    try:
        gc.disable()
        with gzip.open(path,'rb') as stream:record=pickle.load(stream)
        if record['schema']!='rxn_core.aam_checkpoint/v1':
            raise ValueError('Unsupported AAM checkpoint schema')
        return _aam_result(record,record['graph'])
    finally:
        if enabled:gc.enable()


def write_graph_checkpoint(graph,path):
    """Save a trusted, finalized cut graph without expanding shared values."""
    _write_graph_checkpoint(graph,path,'rxn_core.finalized_cut/v1')


def _write_graph_checkpoint(graph,path,schema):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    with gzip.open(temporary,'wb',compresslevel=1) as stream:
        pickle.dump((schema,graph),stream,protocol=5)
    temporary.replace(path)


def read_graph_checkpoint(path):
    """Load only internally produced finalized-cut checkpoints."""
    with gzip.open(path,'rb') as stream:schema,graph=pickle.load(stream)
    if schema!='rxn_core.finalized_cut/v1':raise ValueError('Unsupported cut checkpoint schema')
    return graph


def read_aam(stream):
    """Load owned, acyclic JSON directly into the typed graph."""
    enabled=gc.isenabled()
    try:
        gc.disable()
        return aam_from_record(json.load(stream),copy_graph=False)
    finally:
        if enabled:gc.enable()


def raw_cut_paths(directory):
    """Completed raw cuts in either supported format, in canonical cut order."""
    directory=Path(directory);paths={}
    for path in (*directory.glob('cut_*.json'), *directory.glob('cut_*.raw.pkl.gz')):
        index=int(path.name.split('_')[1].split('.')[0])
        if index in paths:raise ValueError(f'Duplicate raw cut checkpoint {index}')
        paths[index]=path
    return [paths[index] for index in sorted(paths)]


def write_raw_cut(graph,path):
    """Persist the requested raw-cut format atomically; binary is trusted-only."""
    path=Path(path)
    if path.name.endswith('.raw.pkl.gz'):
        _write_graph_checkpoint(graph,path,'rxn_core.raw_cut/v1')
    elif path.suffix=='.json':
        temporary=path.with_suffix('.json.tmp')
        # Interchange JSON uses the fast C encoder. Large persistent searches
        # select compact checkpoints to avoid expanding shared graph payloads.
        temporary.write_text(json.dumps(graph.to_record(copy=False)))
        temporary.replace(path)
    else:raise ValueError(f'Unsupported raw cut format: {path}')


def read_raw_cut(path, *, tuple_pool=None):
    """Read only trusted internal binary cuts; JSON remains interchange-safe."""
    path=Path(path)
    if path.name.endswith('.raw.pkl.gz'):
        with gzip.open(path,'rb') as stream:schema,graph=pickle.load(stream)
        if schema!='rxn_core.raw_cut/v1':raise ValueError('Unsupported raw cut checkpoint schema')
        return graph
    if path.suffix=='.json':
        return AAMSearchGraph.from_record(json.loads(path.read_bytes()),copy=False,tuple_pool=tuple_pool)
    raise ValueError(f'Unsupported raw cut format: {path}')


def aam_json(result):
    """Encode the complete record while borrowing its immutable graph data."""
    enabled=gc.isenabled()
    try:
        gc.disable()
        return json.dumps(aam_record(result,copy_graph=False))+'\n'
    finally:
        if enabled:gc.enable()


def write_aam_bundle(result: AAMResult, output_directory):
    """Save raw AAM and an offline graph/path viewer; never rerun matching."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    serialized = aam_json(result)
    (output / 'aam.json').write_text(serialized)
    assets = Path(__file__).parent / 'static'
    page = (assets / 'aam_search.html').read_text()
    from .viewers import style_document
    page = style_document(page, layout='search_graph')
    page = page.replace('__LIBRARY__', (assets / '3Dmol-min.js').read_text())
    page = page.replace('__DATA__', serialized.replace('<', '\\u003c'))
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
    from .viewers import comparison_document, reaction_html
    problem = result.analytical.aam.problem
    endpoints = [dict(elements=e.elements, coordinates=e.coordinates.tolist(), wbo=e.wbo.tolist())
                 for e in (problem.reactant, problem.product)]
    records = []
    for index, mechanism in enumerate(result.mechanisms, 1):
        mapping = _mapping_record(mechanism.mapping)
        pairs = sorted((int(a), int(b)) for a, b in mapping.items())
        m = dict(pairs)
        events = []
        for kind, bonds in [('broken', mechanism.broken_bonds), ('formed', mechanism.formed_bonds)]:
            for a, b in bonds:
                events.append(dict(kind=kind, r=[a,b], p=[m[a],m[b]],
                                   wbo=[float(problem.reactant.wbo[a,b]),float(problem.product.wbo[m[a],m[b]])]))
        records.append(dict(name=f'Mechanism {index} · RMSD {mechanism.fixed_mapping_rmsd:.4f} Å',
                            mapping=pairs, events=events, provenance=dict(selected_branch=mechanism.selected_branch_index)))
    if not records:
        return '<!doctype html><p>No selected R/P mechanism.</p>'
    document = comparison_document(dict(index=0, name=problem.name or 'R/P alignment',
                                        endpoints=endpoints, records=records))
    return reaction_html(document)


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

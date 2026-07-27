"""Build a portable R/P_final deliverable from NEB orientation results."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import tarfile
import tempfile
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from rxn_core.chemistry_computations import parse_xyz, write_xyz_str
from .neb_interpolation import (
    audit_internal_coordinate_interpolation,
    internal_coordinate_images,
    write_interpolation_xyz,
)


PACKAGE_VERSION = "rxn_core.neb_orientation_deliverable/v4"
_DATA_PREFIX = "const DATA = "
_DATA_SUFFIX = "\nconst COLORS="
_REMOTE_RENDERER = (
    '<script src="https://3dmol.org/build/3Dmol-min.js"></script>'
)


class DeliverableError(ValueError):
    """The source data cannot produce a trustworthy deliverable."""


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: Mapping, *, compact: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compact:
        text = json.dumps(value, separators=(",", ":"), ensure_ascii=False)
    else:
        text = json.dumps(
            value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    path.write_text(text, encoding="utf-8")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _mapping(raw: Mapping) -> dict[int, int]:
    return {int(r): int(p) for r, p in dict(raw).items()}


def _mapping_list(raw: Mapping, atom_count: int) -> list[int]:
    mapping = _mapping(raw)
    if set(mapping) != set(range(atom_count)):
        raise DeliverableError("mapping does not cover every reactant atom")
    if set(mapping.values()) != set(range(atom_count)):
        raise DeliverableError("mapping is not bijective")
    return [mapping[r] for r in range(atom_count)]


def _same_coordinates(left: np.ndarray, right: np.ndarray) -> bool:
    return left.shape == right.shape and bool(
        np.allclose(left, right, atol=5.1e-7, rtol=0.0))


def _extract_original_viewer(
    path: Path,
) -> tuple[str, dict, str]:
    """Split the repository viewer around its embedded all-case JSON."""
    text = path.read_text(encoding="utf-8")
    start = text.find(_DATA_PREFIX)
    if start < 0:
        raise DeliverableError("original viewer has no embedded DATA object")
    value_start = start + len(_DATA_PREFIX)
    marker = text.find(_DATA_SUFFIX, value_start)
    if marker < 0:
        raise DeliverableError("original viewer DATA terminator is missing")
    raw = text[value_start:marker]
    if not raw.endswith(";"):
        raise DeliverableError("original viewer DATA is not semicolon-terminated")
    data = json.loads(raw[:-1])
    return text[:start], data, text[marker:]


def _replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise DeliverableError(
            f"original viewer template changed: expected one occurrence, "
            f"found {count}: {old[:80]!r}")
    return text.replace(old, new, 1)


def _viewer_shell(
    prefix: str,
    suffix: str,
    renderer_javascript: str,
    *,
    interpolation_image_count: int | None = None,
) -> tuple[str, str]:
    """Adapt the original all-case viewer without replacing its design."""
    inline_renderer = (
        "<script>\n"
        + renderer_javascript.replace("</script", r"<\/script")
        + "\n</script>"
    )
    prefix = _replace_once(prefix, _REMOTE_RENDERER, inline_renderer)

    replacements = {
        "R/P AAM + local Kabsch viewer":
            "R/P_final internal-coordinate path viewer",
        "R/P AAM + local Kabsch alignments":
            "R/P_final endpoint assignments",
        "AAM order + local Kabsch":
            "P_final (proper global fit)",
        "AAM order before Kabsch":
            "Final mapping before global fit",
        ">local fragment</label>":
            ">mutable shuffle atoms</label>",
        "mandatory local Kabsch":
            "AAM-constrained endpoint orientation",
        "AAM mappings preserved":
            "AAM-constrained final mappings",
        "local fragment=<b>":
            "mutable shuffle atoms=<b>",
        "Kabsch RMSD=<b>":
            "fixed-anchor RMSD=<b>",
        "fragment detector=":
            "orientation=",
        ">raw mechanism</a>":
            ">mechanism record</a>",
        ">local match</a>":
            ">orientation audit</a>",
    }
    prefix = prefix.replace(
        "</div>\n<div class=\"panels\">",
        "</div>\n"
        "<div class=\"legend\">"
        "<span><i class=\"dot\" style=\"background:#d62728\"></i>broken</span>"
        "<span><i class=\"dot\" style=\"background:#2ca02c\"></i>formed</span>"
        "<span><i class=\"dot\" style=\"background:#d4af37\"></i>core</span>"
        "<span><i class=\"dot\" style=\"background:#00bcd4\"></i>mutable</span>"
        "<span><i class=\"dot\" style=\"background:#e6550d\"></i>mapping change</span>"
        "<span>atom labels are zero-based R indices</span>"
        "</div>\n<div class=\"panels\">",
        1,
    )
    for old, new in replacements.items():
        prefix = prefix.replace(old, new)
        suffix = suffix.replace(old, new)

    suffix = _replace_once(
        suffix,
        "function styleViewer(viewer,elements,coords,degMap,"
        "fragmentIndices,bondPairs,bondColor,labels){",
        "function styleViewer(viewer,elements,coords,degMap,"
        "fragmentIndices,bondPairs,bondColor,labels,coreIndices,"
        "changedIndices){",
    )
    suffix = _replace_once(
        suffix,
        "  if($('fragment').checked)fragmentIndices.forEach(i=>{"
        "if(coords[i])viewer.addSphere({center:{x:coords[i][0],"
        "y:coords[i][1],z:coords[i][2]},radius:.28,color:'#00bcd4',"
        "alpha:.28})});",
        "  if($('fragment').checked)fragmentIndices.forEach(i=>{"
        "if(coords[i])viewer.addSphere({center:{x:coords[i][0],"
        "y:coords[i][1],z:coords[i][2]},radius:.28,color:'#00bcd4',"
        "alpha:.28})});\n"
        "  coreIndices.forEach(i=>{if(coords[i])viewer.addSphere({"
        "center:{x:coords[i][0],y:coords[i][1],z:coords[i][2]},"
        "radius:.40,color:'#d4af37',alpha:.38})});\n"
        "  changedIndices.forEach(i=>{if(coords[i])viewer.addSphere({"
        "center:{x:coords[i][0],y:coords[i][1],z:coords[i][2]},"
        "radius:.32,color:'#e6550d',alpha:.34})});",
    )
    suffix = _replace_once(
        suffix,
        "  const pBonds=native?m.formed_bonds_P:m.formed_bonds_R;",
        "  const pBonds=native?m.formed_bonds_P:m.formed_bonds_R;\n"
        "  const pCore=native?m.core_atoms_P:m.core_atoms;\n"
        "  const pChanged=native?m.mapping_change_atoms_P:"
        "m.mapping_change_atoms_R;",
    )
    suffix = _replace_once(
        suffix,
        "styleViewer(viewerR,c.reactant.elements,c.reactant.coords,"
        "makeDegMap(m,'R',false),m.local_fragment_R,m.broken_bonds_R,"
        "'red',rLabels);",
        "styleViewer(viewerR,c.reactant.elements,c.reactant.coords,"
        "makeDegMap(m,'R',false),m.local_fragment_R,m.broken_bonds_R,"
        "'red',rLabels,m.core_atoms,m.mapping_change_atoms_R);",
    )
    suffix = _replace_once(
        suffix,
        "styleViewer(viewerP,pElements,pCoords,makeDegMap(m,'P',native),"
        "pFragment,pBonds,'green',pLabels);",
        "styleViewer(viewerP,pElements,pCoords,makeDegMap(m,'P',native),"
        "pFragment,pBonds,'green',pLabels,pCore,pChanged);",
    )
    suffix = _replace_once(
        suffix,
        "$('productLabel').textContent=native?'native P order':"
        "(frame==='aam'?'AAM / R order':'AAM / R order + Kabsch');",
        "$('productLabel').textContent=(native?'native P order':"
        "(frame==='aam'?'final mapping / R order before global fit':"
        "'P_final / R order + proper global fit'))+"
        "` · mechanism #${m.id}`;",
    )
    suffix = _replace_once(
        suffix,
        "`AAM mapping: full R→P bijection (${m.mapping_RP.length} atoms). "
        "Symmetry payload: <span class=\"${m.symmetry_payload_present?"
        "'ok':'warn'}\">${m.symmetry_payload_present?'preserved':"
        "'none generated'}</span>; viewer groups=${m.symmetry_groups.length}. `+",
        "`Final AAM mapping: full R→P bijection "
        "(${m.mapping_RP.length} atoms); changes=<b>"
        "${m.mapping_change_atoms_R.length}</b>. "
        "Orientation: <span class=\"${m.final_orientation_violation_count"
        "===0?'ok':'warn'}\">${esc(m.orientation_status)}</span>; "
        "undefined frames=${m.undefined_frame_count}; "
        "viewer groups=${m.symmetry_groups.length}. `+",
    )
    suffix = _replace_once(
        suffix,
        "downloadText(`${safeName(c.step_id)}_R.xyz`,"
        "xyzBody(c.reactant.elements,c.reactant.coords));",
        "downloadText(`${safeName(c.step_id)}_mechanism_${m.id}_R.xyz`,"
        "xyzBody(c.reactant.elements,c.reactant.coords));",
    )
    suffix = _replace_once(
        suffix,
        "setTimeout(()=>downloadText(`${safeName(c.step_id)}_P_${frame}_"
        "mechanism_${m.id}.xyz`,xyzBody(pe,pc)),150)",
        "setTimeout(()=>downloadText(`${safeName(c.step_id)}_mechanism_"
        "${m.id}_${native?'P_native.xyz':"
        "(frame==='aam'?'P_final_unfitted.xyz':'P_final.xyz')}`,"
        "xyzBody(pe,pc)),150)",
    )
    if interpolation_image_count is not None:
        prefix, suffix = _add_interpolation_viewer(
            prefix,
            suffix,
            interpolation_image_count=interpolation_image_count,
        )
    return prefix, suffix


def _add_interpolation_viewer(
    prefix: str,
    suffix: str,
    *,
    interpolation_image_count: int,
) -> tuple[str, str]:
    """Add a scrub/play control to the preserved all-case viewer."""
    if interpolation_image_count < 2:
        raise DeliverableError(
            "interpolation image count must be at least two")
    prefix = _replace_once(
        prefix,
        "    </label>\n"
        "    <label><input type=\"checkbox\" id=\"atomIndices\">"
        "atom indices</label>",
        "    </label>\n"
        "    <button id=\"interpPlay\" type=\"button\">Play</button>\n"
        "    <label>internal-coordinate path\n"
        "      <input id=\"interpFrame\" type=\"range\" min=\"0\" "
        f"max=\"{interpolation_image_count - 1}\" step=\"1\" "
        f"value=\"{interpolation_image_count - 1}\">\n"
        "    </label>\n"
        f"    <output id=\"interpValue\">{interpolation_image_count - 1} "
        f"/ {interpolation_image_count - 1}</output>\n"
        "    <label><input type=\"checkbox\" id=\"atomIndices\">"
        "atom indices</label>",
    )
    prefix = _replace_once(
        prefix,
        ".controls input{vertical-align:middle;margin:0 4px 0 0}",
        ".controls input{vertical-align:middle;margin:0 4px 0 0}"
        ".controls input[type=range]{width:145px;padding:0}"
        ".controls output{min-width:74px;font:11px "
        "ui-monospace,SFMono-Regular,monospace}",
    )
    suffix = _replace_once(
        suffix,
        "let caseIndex=0,mechanismId=null,viewerR=null,viewerP=null;",
        "let caseIndex=0,mechanismId=null,viewerR=null,viewerP=null,"
        "playTimer=null;",
    )
    suffix = _replace_once(
        suffix,
        "function xyzBody(elements,coords){return elements.length+'\\nviewer"
        "\\n'+elements.map((e,i)=>`${e} ${coords[i][0]} ${coords[i][1]} "
        "${coords[i][2]}`).join('\\n')+'\\n'}",
        "function xyzBody(elements,coords){return elements.length+'\\nviewer"
        "\\n'+elements.map((e,i)=>`${e} ${coords[i][0]} ${coords[i][1]} "
        "${coords[i][2]}`).join('\\n')+'\\n'}\n"
        "function interpIndex(){return Number($('interpFrame').value)}\n"
        "function stopPlayback(){if(playTimer!==null){clearInterval("
        "playTimer);playTimer=null}$('interpPlay').textContent='Play'}\n"
        f"function resetInterpolation(){{$('interpFrame').value='"
        f"{interpolation_image_count - 1}';$('interpValue').textContent='"
        f"{interpolation_image_count - 1} / {interpolation_image_count - 1}'}}\n"
        "function updateInterpolationControls(){const enabled="
        "$('productFrame').value==='kabsch';$('interpFrame').disabled="
        "!enabled;$('interpPlay').disabled=!enabled;if(!enabled)"
        "stopPlayback()}\n"
        "function togglePlayback(){if(playTimer!==null){stopPlayback();"
        "return}if($('productFrame').value!=='kabsch')return;"
        f"if(Number($('interpFrame').value)>="
        f"{interpolation_image_count - 1})"
        "$('interpFrame').value='0';$('interpPlay').textContent='Pause';"
        "render(false);playTimer=setInterval(()=>{let i=Number("
        f"$('interpFrame').value)+1;if(i>={interpolation_image_count - 1})"
        f"{{i={interpolation_image_count - 1};"
        "$('interpFrame').value=String(i);render(false);stopPlayback();"
        "return}$('interpFrame').value=String(i);render(false)},180)}",
    )
    suffix = _replace_once(
        suffix,
        "function styleViewer(viewer,elements,coords,degMap,"
        "fragmentIndices,bondPairs,bondColor,labels,coreIndices,"
        "changedIndices){",
        "function styleViewer(viewer,elements,coords,degMap,"
        "fragmentIndices,bondPairs,bondColor,labels,coreIndices,"
        "changedIndices,fit){",
    )
    suffix = _replace_once(
        suffix,
        "  viewer.zoomTo();viewer.render()",
        "  if(fit)viewer.zoomTo();viewer.render()",
    )
    suffix = _replace_once(
        suffix,
        "function render(){\n"
        "  const c=currentCase(),m=currentMechanism(),"
        "frame=$('productFrame').value,native=frame==='native';",
        "function render(fit=true){\n"
        "  const c=currentCase(),m=currentMechanism(),"
        "frame=$('productFrame').value,native=frame==='native',"
        "pathIndex=frame==='kabsch'?interpIndex():"
        f"{interpolation_image_count - 1};",
    )
    suffix = _replace_once(
        suffix,
        "  const pCoords=native?c.product_native.coords:"
        "(frame==='aam'?m.product_aam_order:m.product_local_kabsch);",
        "  const pCoords=native?c.product_native.coords:"
        "(frame==='aam'?m.product_aam_order:"
        "m.interpolation.images[pathIndex]);",
    )
    suffix = _replace_once(
        suffix,
        "  $('viewerR').innerHTML='';$('viewerP').innerHTML='';\n"
        "  viewerR=$3Dmol.createViewer('viewerR',"
        "{backgroundColor:'white'});viewerP=$3Dmol.createViewer("
        "'viewerP',{backgroundColor:'white'});",
        "  if(!viewerR)viewerR=$3Dmol.createViewer('viewerR',"
        "{backgroundColor:'white'});if(!viewerP)viewerP=$3Dmol."
        "createViewer('viewerP',{backgroundColor:'white'});\n"
        "  [viewerR,viewerP].forEach(v=>{v.removeAllModels();"
        "v.removeAllShapes();v.removeAllLabels()});",
    )
    suffix = _replace_once(
        suffix,
        "styleViewer(viewerR,c.reactant.elements,c.reactant.coords,"
        "makeDegMap(m,'R',false),m.local_fragment_R,m.broken_bonds_R,"
        "'red',rLabels,m.core_atoms,m.mapping_change_atoms_R);",
        "styleViewer(viewerR,c.reactant.elements,c.reactant.coords,"
        "makeDegMap(m,'R',false),m.local_fragment_R,m.broken_bonds_R,"
        "'red',rLabels,m.core_atoms,m.mapping_change_atoms_R,fit);",
    )
    suffix = _replace_once(
        suffix,
        "styleViewer(viewerP,pElements,pCoords,makeDegMap(m,'P',native),"
        "pFragment,pBonds,'green',pLabels,pCore,pChanged);",
        "styleViewer(viewerP,pElements,pCoords,makeDegMap(m,'P',native),"
        "pFragment,pBonds,'green',pLabels,pCore,pChanged,fit);",
    )
    suffix = _replace_once(
        suffix,
        "$('productLabel').textContent=(native?'native P order':"
        "(frame==='aam'?'final mapping / R order before global fit':"
        "'P_final / R order + proper global fit'))+"
        "` · mechanism #${m.id}`;",
        "$('productLabel').textContent=(native?'native P order':"
        "(frame==='aam'?'final mapping / R order before global fit':"
        "`internal-coordinate preview image ${pathIndex}/"
        f"{interpolation_image_count - 1}`))"
        "+` · mechanism #${m.id}`;\n"
        "  $('interpValue').textContent=`${$('interpFrame').value} / "
        f"{interpolation_image_count - 1}`;"
        "updateInterpolationControls();",
    )
    suffix = _replace_once(
        suffix,
        "`Final AAM mapping: full R→P bijection "
        "(${m.mapping_RP.length} atoms); changes=<b>"
        "${m.mapping_change_atoms_R.length}</b>. "
        "Orientation: <span class=\"${m.final_orientation_violation_count"
        "===0?'ok':'warn'}\">${esc(m.orientation_status)}</span>; "
        "undefined frames=${m.undefined_frame_count}; "
        "viewer groups=${m.symmetry_groups.length}. `+",
        "`Final AAM mapping: full R→P bijection "
        "(${m.mapping_RP.length} atoms); changes=<b>"
        "${m.mapping_change_atoms_R.length}</b>. "
        "Internal-coordinate preview: <span class=\""
        "${m.interpolation.overall_status"
        "==='pass'?'ok':'warn'}\">${esc(m.interpolation.overall_status)}"
        "</span>; orientation-zero frames=<b>"
        "${m.interpolation.orientation_zero_frame_count}</b>; "
        "hard close approaches=<b>"
        "${m.interpolation.hard_collision_count}</b>. "
        "Native parity variables=<b>"
        "${m.native_index_chirality?"
        "m.native_index_chirality.candidate_search."
        "parity_variable_count:'n/a'}</b>; GF(2) equations=<b>"
        "${m.native_index_chirality?"
        "m.native_index_chirality.candidate_search."
        "gf2_equation_count:'n/a'}</b>; solved routes=<b>"
        "${m.native_index_chirality?"
        "m.native_index_chirality.candidate_search."
        "gf2_solved_route_count:'n/a'}</b>"
        " → chirality-safe candidates=<b>${m.native_index_chirality?"
        "m.native_index_chirality.allowed_candidate_count:'n/a'}</b>; "
        "immutable endpoint mismatches=<b>"
        "${m.native_index_chirality?"
        "m.native_index_chirality.immutable_source_mismatch_count:"
        "'n/a'}</b>. `+",
    )
    suffix = _replace_once(
        suffix,
        "`<a href=\"${esc(m.files.directory)}/${esc(m.files.mapping)}\">"
        "mapping CSV</a> · <a href=\"${esc(m.files.directory)}/"
        "${esc(m.files.mechanism)}\">mechanism record</a> · "
        "<a href=\"${esc(m.files.directory)}/"
        "${esc(m.files.local_match)}\">orientation audit</a>`;",
        "`<a href=\"${esc(m.files.directory)}/"
        "${esc(m.files.trajectory)}\">"
        f"{interpolation_image_count}-image XYZ</a> · "
        "<a href=\"${esc(m.files.directory)}/"
        "${esc(m.files.interpolation_report)}\">path audit</a> · "
        "<a href=\"${esc(m.files.directory)}/${esc(m.files.mechanism)}\">"
        "mapping/orientation metadata</a>`;",
    )
    suffix = _replace_once(
        suffix,
        "b.onclick=()=>{mechanismId=m.id;render()}",
        "b.onclick=()=>{stopPlayback();resetInterpolation();"
        "mechanismId=m.id;render(true)}",
    )
    suffix = _replace_once(
        suffix,
        "function setCase(index){caseIndex=(index+DATA.case_count)%"
        "DATA.case_count;$('caseSelect').value=String(caseIndex);"
        "buildMechanisms()}",
        "function setCase(index){stopPlayback();resetInterpolation();"
        "caseIndex=(index+DATA.case_count)%DATA.case_count;"
        "$('caseSelect').value=String(caseIndex);buildMechanisms()}",
    )
    suffix = _replace_once(
        suffix,
        "['productFrame','atomIndices','symmetry','fragment'].forEach("
        "id=>$(id).onchange=render);$('downloadPair').onclick=downloadPair;",
        "['atomIndices','symmetry','fragment'].forEach(id=>"
        "$(id).onchange=()=>render(false));"
        "$('productFrame').onchange=()=>{stopPlayback();"
        "resetInterpolation();updateInterpolationControls();render(true)};"
        "$('interpFrame').oninput=()=>{stopPlayback();render(false)};"
        "$('interpPlay').onclick=togglePlayback;"
        "$('downloadPair').onclick=downloadPair;"
        "document.addEventListener('visibilitychange',()=>{"
        "if(document.hidden)stopPlayback()});",
    )
    return prefix, suffix


def _compact_mechanism_record(
    *,
    step_id: str,
    source_relative_path: str,
    view_mechanism: Mapping,
    orientation: Mapping,
    native_index_chirality: Mapping | None = None,
) -> dict:
    record = {
        "schema_version": PACKAGE_VERSION,
        "step_id": step_id,
        "mechanism_id": int(view_mechanism["id"]),
        "label": view_mechanism.get("label"),
        "is_final_selected": bool(view_mechanism.get("selected")),
        "cut": view_mechanism.get("cut"),
        "dedup_count": int(view_mechanism.get("dedup_count") or 1),
        "dedup_cuts": view_mechanism.get("dedup_cuts") or [],
        "source_mechanism": source_relative_path,
        "source_mapping_RP": orientation["source_mapping_RP"],
        "selected_neb_mapping_RP": orientation[
            "selected_neb_mapping_RP"],
        "source_mapping_sha256": orientation["source_mapping_sha256"],
        "selected_mapping_sha256": orientation[
            "selected_mapping_sha256"],
        "core_atoms_R": list(view_mechanism.get("core_atoms") or ()),
        "source_broken_bonds_R": list(
            view_mechanism.get("broken_bonds_R") or ()),
        "source_formed_bonds_R": list(
            view_mechanism.get("formed_bonds_R") or ()),
        "source_formed_bonds_P": list(
            view_mechanism.get("formed_bonds_P") or ()),
        "broken_bonds_R": list(
            orientation["selected_broken_bonds_R"]),
        "formed_bonds_R": list(
            orientation["selected_formed_bonds_R"]),
        "formed_bonds_P": list(
            orientation["selected_formed_bonds_P"]),
        "orientation": {
            "status": orientation["status"],
            "selected_witness_index": int(
                orientation["selected_witness_index"]),
            "encoded_candidate_count": int(
                orientation["encoded_candidate_count"]),
            "selected_candidate_provenance": orientation[
                "selected_candidate_provenance"],
            "fixed_r_atoms": list(orientation["fixed_r_atoms"]),
            "mutable_r_atoms": list(orientation["mutable_r_atoms"]),
            "frame_count": len(orientation["orientation_frames"]),
            "defined_frames": orientation["orientation_frames"],
            "undefined_frame_count": orientation[
                "undefined_frame_count"],
            "source_violation_count": orientation[
                "source_orientation_violation_count"],
            "final_violation_count": orientation[
                "final_orientation_violation_count"],
            "mapping_changes": orientation["mapping_changes"],
            "allowed_shuffle_blocks": orientation[
                "allowed_shuffle_blocks"],
            "geometry_tiebreak": orientation["geometry_tiebreak"],
            "invariants": orientation["invariants"],
            "initial_path_certified": False,
        },
        "files": {
            "reactant": "R.xyz",
            "product_final": "P_final.xyz",
        },
    }
    if native_index_chirality:
        record["native_index_chirality"] = dict(native_index_chirality)
    return record


def _native_index_chirality_summary(mechanism: Mapping) -> dict | None:
    """Compact the native bounded-candidate audit for package metadata."""
    record = (
        (mechanism.get("branch_symmetry") or {})
        .get("index_chirality") or {}
    )
    if record.get("status") not in {"applied", "conflict"}:
        return None
    candidate_search = record.get("candidate_search") or {}
    candidate_search_summary = {
        key: candidate_search.get(key)
        for key in (
            "semantics",
            "seed_route_count",
            "rejected_seed_route_count",
            "explicit_witness_seed_count",
            "nested_alternate_seed_count",
            "fragment_parity_seed_count",
            "parity_variable_count",
            "gf2_equation_count",
            "gf2_solved_route_count",
            "unique_candidate_evaluation_count",
        )
        if candidate_search.get(key) is not None
    }
    candidate_search_summary.setdefault(
        "unique_candidate_evaluation_count",
        len(candidate_search.get("candidate_evaluations") or ()),
    )
    immutable_mismatches = [
        {
            key: frame.get(key)
            for key in (
                "id",
                "center_R",
                "neighbors_R_index_order",
                "reaction_event_incident",
                "reason",
                "source_mismatch_details",
            )
            if frame.get(key) is not None
        }
        for frame in record.get("immutable_frames") or ()
        if frame.get("source_index_chirality_mismatch")
    ]
    return {
        "schema_version": record.get("schema_version"),
        "policy": record.get("policy"),
        "status": record.get("status"),
        "source_mapping_sha256": record.get("source_mapping_sha256"),
        "selected_mapping_sha256": record.get("selected_mapping_sha256"),
        "source_violation_count": int(
            record.get("source_index_chirality_violation_count", 0)),
        "selected_violation_count": int(
            record.get("selected_index_chirality_violation_count", 0)),
        "defined_frame_count": int(record.get("defined_frame_count", 0)),
        "undefined_frame_count": int(record.get("undefined_frame_count", 0)),
        "switchable_r_atoms": [
            int(value)
            for value in record.get("switchable_r_atoms") or ()
        ],
        "candidate_search": candidate_search_summary,
        "allowed_candidate_count": int(
            record.get("allowed_candidate_count", 0)),
        "selected_candidate_id": record.get("selected_candidate_id"),
        "immutable_frame_count": int(
            record.get("immutable_frame_count", 0)),
        "immutable_source_mismatch_count": int(
            record.get("immutable_source_mismatch_count", 0)),
        "immutable_mismatch_frames": immutable_mismatches,
        "mapping_changes": list(record.get("mapping_changes") or ()),
        "selection_rule": record.get("selection_rule"),
        "invariants": dict(record.get("invariants") or {}),
    }


def _write_pair(
    directory: Path,
    elements: Sequence[str],
    coords_R: np.ndarray,
    coords_P_final: np.ndarray,
    step_id: str,
    mechanism_id: int,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "R.xyz").write_text(
        write_xyz_str(
            elements,
            coords_R,
            f"{step_id} mechanism {mechanism_id:03d}; R order",
        ),
        encoding="utf-8",
    )
    (directory / "P_final.xyz").write_text(
        write_xyz_str(
            elements,
            coords_P_final,
            f"{step_id} mechanism {mechanism_id:03d}; "
            "P_final in R order; proper global fit",
        ),
        encoding="utf-8",
    )


def _package_readme(
    case_count: int,
    mechanism_count: int,
    interpolation_image_count: int | None = None,
) -> str:
    interpolation = ""
    if interpolation_image_count is not None:
        interpolation = f"""
Each mechanism also contains:

```text
internal_coordinate_interpolation.xyz
interpolation_report.json
```

`internal_coordinate_interpolation.xyz` contains {interpolation_image_count}
endpoint-inclusive images in fixed R row order. The band is generated from
linearly interpolated all-pair distances with IDPP. Because distances cannot
distinguish mirror images, each native hard tetrahedron is then placed on its
recorded R-sign branch without changing its six internal pair distances.

The report checks every adjacent image segment for an exact orientation zero
and analytically checks every pair's closest approach on the piecewise path.
The viewer slider and Play button show these exact generated images; the HTML
does not recreate a Cartesian endpoint interpolation.

This remains an initial-path preview, not an optimized NEB path. Pair-distance
warnings and displacement are reported for review but do not alter the AAM
mapping.
"""
    return f"""# Endpoint-orientation and internal-coordinate-path deliverable

This package contains {case_count} cases and {mechanism_count} mechanisms.

For every mechanism:

```text
cases/<step>/mechanisms/mechanism_<NNN>/R.xyz
cases/<step>/mechanisms/mechanism_<NNN>/P_final.xyz
cases/<step>/mechanisms/mechanism_<NNN>/metadata.json
```

`R.xyz` and `P_final.xyz` have identical elements and identical zero-based
reactant atom order. `P_final.xyz` uses the selected AAM-constrained mapping
and one proper whole-product rigid fit. No reflection or independent local
rotation was applied.

Open `viewer.html` directly in a browser. It contains the renderer and all
case/mechanism coordinates, so viewing, navigation, and XYZ downloads do not
need a server or network connection.

The viewer defaults to `P_final`. It can also show the same selected mapping
before the global fit or native product order for diagnosis. Red dashed bonds
are broken bonds, green dashed bonds are formed bonds, gold marks reaction-core
atoms, cyan marks AAM-authorized mutable atoms, and orange marks mapping
changes.
{interpolation}

Endpoint orientation is audited for every defined local frame. Near-planar
frames without a stable sign are reported as undefined in `metadata.json`.
This package does not claim that an unconstrained later NEB optimization can
never undergo a physical inversion.
"""


def _write_checksums(root: Path) -> None:
    paths = sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.txt"
    )
    lines = [
        f"{_sha256_file(path)}  {path.relative_to(root).as_posix()}"
        for path in paths
    ]
    (root / "SHA256SUMS.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")


def _write_reproducible_archive(root: Path, archive_path: Path) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with archive_path.open("wb") as raw:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw,
            mtime=0,
        ) as compressed:
            with tarfile.open(
                fileobj=compressed,
                mode="w",
                format=tarfile.PAX_FORMAT,
            ) as archive:
                paths = [root, *sorted(root.rglob("*"))]
                for path in paths:
                    relative = path.relative_to(root)
                    arcname = (
                        Path(root.name)
                        if not relative.parts
                        else Path(root.name) / relative
                    )
                    info = archive.gettarinfo(
                        str(path), arcname=arcname.as_posix())
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mtime = 0
                    if path.is_file():
                        with path.open("rb") as handle:
                            archive.addfile(info, handle)
                    else:
                        archive.addfile(info)


def build_deliverable(
    source_root: str | Path,
    orientation_root: str | Path,
    output_root: str | Path,
    renderer_javascript_path: str | Path,
    *,
    renderer_license_path: str | Path | None = None,
    archive_path: str | Path | None = None,
    interpolation_image_count: int | None = None,
) -> dict:
    """Build one non-overwriting, self-contained all-case deliverable."""
    source_root = Path(source_root).resolve()
    orientation_root = Path(orientation_root).resolve()
    output_root = Path(output_root).resolve()
    renderer_path = Path(renderer_javascript_path).resolve()
    archive = None if archive_path is None else Path(archive_path).resolve()
    if output_root.exists():
        raise FileExistsError(
            f"refusing to overwrite existing output: {output_root}")
    if archive is not None and archive.exists():
        raise FileExistsError(
            f"refusing to overwrite existing archive: {archive}")

    source_index = _read_json(source_root / "index.json")
    orientation_index = _read_json(orientation_root / "index.json")
    prefix, viewer_data, suffix = _extract_original_viewer(
        source_root / "viewer.html")
    renderer_text = renderer_path.read_text(encoding="utf-8")
    prefix, suffix = _viewer_shell(
        prefix,
        suffix,
        renderer_text,
        interpolation_image_count=interpolation_image_count,
    )

    source_cases = {
        str(case["step_id"]): case
        for case in source_index["cases"]
    }
    view_cases = {
        str(case["step_id"]): case
        for case in viewer_data["cases"]
    }
    orientation_cases = orientation_index["cases"]
    if len(orientation_cases) != len(source_cases) or (
            set(source_cases) != set(view_cases)):
        raise DeliverableError("source, viewer, and orientation case sets differ")

    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(
        prefix=f"{output_root.name}.staging-",
        dir=output_root.parent,
    ))
    package_cases = []
    interpolation_summaries = []
    try:
        for case_position, orientation_case in enumerate(
                orientation_cases, start=1):
            step_id = str(orientation_case["step_id"])
            source_case = source_cases[step_id]
            view_case = view_cases[step_id]
            source_case_dir = source_root / source_case["directory"]
            elements_R, coords_R = parse_xyz(
                source_case_dir / source_case["files"]["reactant"])
            elements_P, coords_P = parse_xyz(
                source_case_dir
                / source_case["files"]["product_original_order"])
            if list(view_case["reactant"]["elements"]) != list(elements_R):
                raise DeliverableError(
                    f"{step_id}: viewer reactant elements differ")
            if list(view_case["product_native"]["elements"]) != list(
                    elements_P):
                raise DeliverableError(
                    f"{step_id}: viewer product elements differ")
            if not _same_coordinates(
                    np.asarray(view_case["reactant"]["coords"]), coords_R):
                raise DeliverableError(
                    f"{step_id}: viewer reactant coordinates differ")
            if not _same_coordinates(
                    np.asarray(view_case["product_native"]["coords"]), coords_P):
                raise DeliverableError(
                    f"{step_id}: viewer product coordinates differ")

            source_mechanisms = {
                int(item["mechanism_id"]): item
                for item in source_case["mechanisms"]
            }
            view_mechanisms = {
                int(item["id"]): item
                for item in view_case["mechanisms"]
            }
            case_mechanisms = []
            for summary in orientation_case["mechanisms"]:
                mechanism_id = int(summary["mechanism_id"])
                source_mechanism = source_mechanisms[mechanism_id]
                view_mechanism = view_mechanisms[mechanism_id]
                record_path = orientation_root / summary["record"]
                orientation = _read_json(record_path)
                mechanism_dir = record_path.parent
                final_elements, coords_P_final = parse_xyz(
                    mechanism_dir / "P_neb_ordered.xyz")
                result_elements, result_coords_R = parse_xyz(
                    mechanism_dir / "R.xyz")
                if list(result_elements) != list(elements_R) or (
                        list(final_elements) != list(elements_R)):
                    raise DeliverableError(
                        f"{step_id}/mechanism_{mechanism_id:03d}: "
                        "R/P_final element order differs")
                if not _same_coordinates(result_coords_R, coords_R):
                    raise DeliverableError(
                        f"{step_id}/mechanism_{mechanism_id:03d}: "
                        "packaged R differs from source")

                selected = _mapping_list(
                    orientation["selected_neb_mapping_RP"],
                    len(elements_R),
                )
                source_mapping = _mapping_list(
                    orientation["source_mapping_RP"],
                    len(elements_R),
                )
                for r, p in enumerate(selected):
                    if str(elements_R[r]) != str(elements_P[p]):
                        raise DeliverableError(
                            f"{step_id}/mechanism_{mechanism_id:03d}: "
                            f"element mismatch R{r}->P{p}")

                relative_dir = (
                    Path("cases") / step_id / "mechanisms"
                    / f"mechanism_{mechanism_id:03d}"
                )
                target_dir = staging / relative_dir
                _write_pair(
                    target_dir,
                    elements_R,
                    coords_R,
                    coords_P_final,
                    step_id,
                    mechanism_id,
                )
                interpolation_report = None
                interpolation_images = None
                if interpolation_image_count is not None:
                    interpolation_images, generation = (
                        internal_coordinate_images(
                            elements_R,
                            coords_R,
                            coords_P_final,
                            orientation["orientation_frames"],
                            interpolation_image_count,
                        )
                    )
                    write_interpolation_xyz(
                        target_dir / "internal_coordinate_interpolation.xyz",
                        elements_R,
                        interpolation_images,
                        step_id=step_id,
                        mechanism_id=mechanism_id,
                    )
                    interpolation_report = (
                        audit_internal_coordinate_interpolation(
                            elements_R,
                            interpolation_images,
                            orientation["orientation_frames"],
                            generation=generation,
                        )
                    )
                    interpolation_report.update({
                        "step_id": step_id,
                        "mechanism_id": mechanism_id,
                        "is_final_selected": bool(
                            summary["is_final_selected"]),
                        "endpoint_sha256": {
                            "R.xyz": _sha256_file(target_dir / "R.xyz"),
                            "P_final.xyz": _sha256_file(
                                target_dir / "P_final.xyz"),
                        },
                    })
                    _write_json(
                        target_dir / "interpolation_report.json",
                        interpolation_report,
                    )
                    interpolation_summaries.append({
                        "step_id": step_id,
                        "mechanism_id": mechanism_id,
                        "is_final_selected": bool(
                            summary["is_final_selected"]),
                        "overall_status": interpolation_report[
                            "overall_status"],
                        "orientation_zero_frame_count": (
                            interpolation_report["checks"]["orientation"][
                                "frames_with_interior_zero_count"]
                        ),
                        "hard_collision_count": (
                            interpolation_report["checks"]["pair_distance"][
                                "hard_collision_count"]
                        ),
                    })

                source_relative_path = (
                    Path(source_case["directory"])
                    / source_mechanism["directory"]
                    / "mechanism.json"
                ).as_posix()
                source_mechanism_record = _read_json(
                    source_case_dir
                    / source_mechanism["directory"]
                    / "mechanism.json"
                )
                native_index_chirality = (
                    _native_index_chirality_summary(
                        source_mechanism_record)
                )
                compact_record = _compact_mechanism_record(
                    step_id=step_id,
                    source_relative_path=source_relative_path,
                    view_mechanism=view_mechanism,
                    orientation=orientation,
                    native_index_chirality=native_index_chirality,
                )
                compact_record.update({
                    "atom_count": len(elements_R),
                    "coordinate_units": "angstrom",
                    "row_contract": {
                        "R_row_r": "native reactant atom r",
                        "P_final_row_r": (
                            "native product atom selected_neb_mapping_RP[r]"
                        ),
                        "same_neb_row_order": True,
                    },
                    "artifact_sha256": {
                        "R.xyz": _sha256_file(target_dir / "R.xyz"),
                        "P_final.xyz": _sha256_file(
                            target_dir / "P_final.xyz"),
                    },
                    "source_sha256": {
                        "mechanism.json": _sha256_file(
                            source_case_dir
                            / source_mechanism["directory"]
                            / "mechanism.json"
                        ),
                        "R.xyz": _sha256_file(
                            source_case_dir
                            / source_case["files"]["reactant"]
                        ),
                        "P_original.xyz": _sha256_file(
                            source_case_dir
                            / source_case["files"][
                                "product_original_order"]
                        ),
                        "neb_orientation.json": _sha256_file(record_path),
                    },
                })
                if interpolation_report is not None:
                    compact_record["interpolation"] = {
                        "path_model": interpolation_report["path_model"],
                        "image_count": interpolation_image_count,
                        "overall_status": interpolation_report[
                            "overall_status"],
                        "initial_internal_coordinate_path_certified": (
                            interpolation_report[
                                "initial_internal_coordinate_path_certified"]
                        ),
                        "optimized_neb_path_certified": False,
                        "orientation_zero_frame_count": (
                            interpolation_report["checks"]["orientation"][
                                "frames_with_interior_zero_count"]
                        ),
                        "hard_collision_count": (
                            interpolation_report["checks"]["pair_distance"][
                                "hard_collision_count"]
                        ),
                    }
                    compact_record["files"].update({
                        "internal_coordinate_interpolation": (
                            "internal_coordinate_interpolation.xyz"),
                        "interpolation_report": (
                            "interpolation_report.json"),
                    })
                    compact_record["artifact_sha256"].update({
                        "internal_coordinate_interpolation.xyz": (
                            _sha256_file(
                                target_dir
                                / "internal_coordinate_interpolation.xyz")
                        ),
                        "interpolation_report.json": _sha256_file(
                            target_dir / "interpolation_report.json"),
                    })
                _write_json(target_dir / "metadata.json", compact_record)

                blocks = list(orientation["allowed_shuffle_blocks"])
                mutable_R = sorted(
                    int(r) for r in orientation["mutable_r_atoms"])
                changed_R = sorted(
                    int(item["r_atom"])
                    for item in orientation["mapping_changes"]
                )
                if not changed_R and native_index_chirality:
                    changed_R = sorted(
                        int(item["r_atom"])
                        for item in native_index_chirality[
                            "mapping_changes"]
                    )
                core_R = [
                    int(value)
                    for value in view_mechanism.get("core_atoms") or ()
                ]
                view_mechanism["source_mapping_RP"] = source_mapping
                view_mechanism["mapping_RP"] = selected
                view_mechanism["broken_bonds_R"] = [
                    list(pair)
                    for pair in orientation["selected_broken_bonds_R"]
                ]
                view_mechanism["formed_bonds_R"] = [
                    list(pair)
                    for pair in orientation["selected_formed_bonds_R"]
                ]
                view_mechanism["formed_bonds_P"] = [
                    list(pair)
                    for pair in orientation["selected_formed_bonds_P"]
                ]
                view_mechanism["product_aam_order"] = [
                    list(map(float, coords_P[selected[r]]))
                    for r in range(len(selected))
                ]
                view_mechanism["product_local_kabsch"] = (
                    np.asarray(coords_P_final, dtype=float).tolist())
                view_mechanism["local_fragment_R"] = mutable_R
                view_mechanism["local_fragment_P"] = sorted(
                    selected[r] for r in mutable_R)
                view_mechanism["symmetry_groups"] = [
                    {
                        "r_atoms": [
                            int(value) for value in block["r_atoms"]],
                        "p_atoms": [
                            int(value) for value in block["p_atoms"]],
                        "assignments": f"{len(block['r_atoms'])}!",
                        "source": "exact_nested_aam_shuffle",
                    }
                    for block in blocks
                ]
                view_mechanism["symmetry_payload_present"] = bool(blocks)
                view_mechanism["fragment_detection"] = orientation["status"]
                geometry = orientation["geometry_tiebreak"]
                view_mechanism["kabsch"] = {
                    "rmsd_angstrom": float(
                        geometry["anchor_rmsd_angstrom"]),
                    "max_residual_angstrom": float(
                        geometry[
                            "maximum_mutable_displacement_angstrom"]),
                    "proper_rotation_determinant": float(
                        geometry["proper_rotation_determinant"]),
                    "fragment_atom_count": len(
                        geometry["anchor_r_atoms"]),
                }
                view_mechanism["core_atoms"] = core_R
                view_mechanism["core_atoms_P"] = [
                    selected[r] for r in core_R]
                view_mechanism["mapping_change_atoms_R"] = changed_R
                view_mechanism["mapping_change_atoms_P"] = [
                    selected[r] for r in changed_R]
                view_mechanism["orientation_status"] = orientation["status"]
                view_mechanism["native_index_chirality"] = (
                    native_index_chirality)
                view_mechanism["undefined_frame_count"] = int(
                    orientation["undefined_frame_count"])
                view_mechanism["final_orientation_violation_count"] = int(
                    orientation["final_orientation_violation_count"])
                view_mechanism["files"] = {
                    "directory": relative_dir.as_posix(),
                    "mapping": "metadata.json",
                    "mechanism": "metadata.json",
                    "local_match": "metadata.json",
                    "reactant": "R.xyz",
                    "product_final": "P_final.xyz",
                }
                if interpolation_report is not None:
                    view_mechanism["interpolation"] = {
                        "image_count": interpolation_image_count,
                        "images": np.asarray(
                            interpolation_images, dtype=float).tolist(),
                        "path_model": interpolation_report["path_model"],
                        "overall_status": interpolation_report[
                            "overall_status"],
                        "orientation_zero_frame_count": (
                            interpolation_report["checks"]["orientation"][
                                "frames_with_interior_zero_count"]
                        ),
                        "orientation_zero_count": (
                            interpolation_report["checks"]["orientation"][
                                "interior_zero_count"]
                        ),
                        "hard_collision_count": (
                            interpolation_report["checks"]["pair_distance"][
                                "hard_collision_count"]
                        ),
                        "minimum_pair_distance_angstrom": (
                            interpolation_report["checks"]["pair_distance"][
                                "closest_pair"][
                                    "minimum_distance_angstrom"]
                        ),
                        "maximum_displacement_angstrom": (
                            interpolation_report["checks"]["displacement"][
                                "all_atoms"]["maximum_angstrom"]
                        ),
                    }
                    view_mechanism["files"].update({
                        "trajectory": (
                            "internal_coordinate_interpolation.xyz"),
                        "interpolation_report": (
                            "interpolation_report.json"),
                    })

                mechanism_summary = {
                    "mechanism_id": mechanism_id,
                    "is_final_selected": bool(
                        summary["is_final_selected"]),
                    "status": orientation["status"],
                    "mapping_change_count": len(changed_R),
                    "undefined_frame_count": int(
                        orientation["undefined_frame_count"]),
                    "directory": relative_dir.as_posix(),
                    "files": {
                        "reactant": "R.xyz",
                        "product_final": "P_final.xyz",
                        "metadata": "metadata.json",
                    },
                }
                if interpolation_report is not None:
                    mechanism_summary["interpolation"] = {
                        "overall_status": interpolation_report[
                            "overall_status"],
                        "orientation_zero_frame_count": (
                            interpolation_report["checks"]["orientation"][
                                "frames_with_interior_zero_count"]
                        ),
                        "hard_collision_count": (
                            interpolation_report["checks"]["pair_distance"][
                                "hard_collision_count"]
                        ),
                        "files": {
                            "trajectory": (
                                "internal_coordinate_interpolation.xyz"),
                            "report": "interpolation_report.json",
                        },
                    }
                case_mechanisms.append(mechanism_summary)

            case_record = {
                "schema_version": PACKAGE_VERSION,
                "step_id": step_id,
                "selected_mechanism_id": int(
                    orientation_case["selected_mechanism_id"]),
                "mechanism_count": len(case_mechanisms),
                "mechanisms": case_mechanisms,
            }
            _write_json(
                staging / "cases" / step_id / "case.json",
                case_record,
            )
            package_cases.append(case_record)
            print(
                f"[{case_position:03d}/{len(orientation_cases):03d}] "
                f"{step_id} mechanisms={len(case_mechanisms)}",
                flush=True,
            )

        viewer_data["schema_version"] = (
            "rp-final-orientation-viewer/v3")
        viewer_data["title"] = (
            "R/P_final internal-coordinate path previews")
        viewer_data["case_count"] = len(package_cases)
        viewer_data["mechanism_count"] = sum(
            case["mechanism_count"] for case in package_cases)
        viewer_json = json.dumps(
            viewer_data, separators=(",", ":"), ensure_ascii=False
        ).replace("</script", r"<\/script")
        viewer_html = (
            prefix + _DATA_PREFIX + viewer_json + ";" + suffix)
        if _REMOTE_RENDERER in viewer_html or "<script src=" in prefix:
            raise DeliverableError("viewer still has an external script")
        (staging / "viewer.html").write_text(
            viewer_html, encoding="utf-8")

        index = {
            "schema_version": PACKAGE_VERSION,
            "case_count": len(package_cases),
            "mechanism_count": viewer_data["mechanism_count"],
            "source_archive": str(source_root),
            "orientation_result": str(orientation_root),
            "viewer": {
                "file": "viewer.html",
                "self_contained": True,
                "renderer": "3Dmol.js embedded inline",
                "renderer_sha256": _sha256_file(renderer_path),
                "all_case_data_embedded": True,
            },
            "coordinate_contract": {
                "reactant_file": "R.xyz",
                "product_file": "P_final.xyz",
                "both_in_reactant_atom_order": True,
                "product_proper_global_fit": True,
                "reflection_used": False,
            },
            "cases": package_cases,
        }
        if interpolation_image_count is not None:
            selected_interpolation = [
                item for item in interpolation_summaries
                if item["is_final_selected"]
            ]
            index["internal_coordinate_interpolation"] = {
                "path_model": (
                    "idpp_pair_distance_with_native_signed_frame_branch"),
                "image_count": interpolation_image_count,
                "mechanism_count": len(interpolation_summaries),
                "selected_mechanism_count": len(selected_interpolation),
                "failed_mechanism_count": sum(
                    item["overall_status"] == "fail"
                    for item in interpolation_summaries
                ),
                "failed_selected_mechanism_count": sum(
                    item["overall_status"] == "fail"
                    for item in selected_interpolation
                ),
                "selected_mechanisms_with_orientation_zero_count": sum(
                    item["orientation_zero_frame_count"] > 0
                    for item in selected_interpolation
                ),
                "selected_orientation_zero_frame_count": sum(
                    item["orientation_zero_frame_count"]
                    for item in selected_interpolation
                ),
                "selected_mechanisms_with_hard_collision_count": sum(
                    item["hard_collision_count"] > 0
                    for item in selected_interpolation
                ),
                "initial_internal_coordinate_path_certified_for_all_selected": all(
                    item["overall_status"] == "pass"
                    for item in selected_interpolation
                ),
                "optimized_neb_path_certified": False,
            }
        _write_json(staging / "deliverable.json", index)
        (staging / "README.md").write_text(
            _package_readme(
                len(package_cases),
                viewer_data["mechanism_count"],
                interpolation_image_count,
            ),
            encoding="utf-8",
        )
        licenses = staging / "LICENSES"
        licenses.mkdir(parents=True, exist_ok=True)
        if renderer_license_path is not None:
            shutil.copyfile(
                Path(renderer_license_path).resolve(),
                licenses / "3Dmol.js.txt",
            )
        else:
            (licenses / "3Dmol.js.txt").write_text(
                "3Dmol.js license was not supplied to the packager.\n",
                encoding="utf-8",
            )
        _write_checksums(staging)
        os.replace(staging, output_root)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    if archive is not None:
        _write_reproducible_archive(output_root, archive)
        archive.with_name(archive.name + ".sha256").write_text(
            f"{_sha256_file(archive)}  {archive.name}\n",
            encoding="utf-8",
        )
    return _read_json(output_root / "deliverable.json")

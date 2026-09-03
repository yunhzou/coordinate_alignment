#!/usr/bin/env python3
"""Build an interactive 3D viewer from merged blind catalog-search results."""
from __future__ import annotations

import argparse
from collections import Counter
from html import escape
import json
from pathlib import Path

from molecule_3d import mol_3d


COLORS = (
    "#2684ff", "#ff8b00", "#8b5cf6", "#00a896", "#e83e8c",
    "#795548", "#00acc1", "#7cb342", "#f4511e", "#5c6bc0",
    "#c0a000", "#6d4c91",
)
def _model(smiles, candidate, color, cache):
    if smiles not in cache:
        cache[smiles] = mol_3d(
            smiles, spread_ions=True, show_hydrogens=True)
    block, coords, elements = cache[smiles]
    return {
        "mol": block,
        "coords": coords,
        "elements": elements,
        "styles": [{"indices": candidate["retained_atoms"], "color": color}],
        "symmetryStyles": [{
            "indices": candidate["symmetry_retained_atoms"], "color": color,
        }],
        "labels": [
            {"atom": int(source), "text": f"P{int(target)}"}
            for source, target in candidate["mapping"]
        ],
        "broken": candidate["boundary_bonds"],
        "formed": [],
    }


def _group_precursors(precursors):
    """Collapse repeated stoichiometric copies into one visual reactant."""
    groups = []
    by_id = {}
    for item in precursors:
        precursor_id = item["precursor_id"]
        group = by_id.get(precursor_id)
        if group is None:
            group = {
                "precursor_id": precursor_id,
                "smiles": item["smiles"],
                "copies": [],
            }
            by_id[precursor_id] = group
            groups.append(group)
        group["copies"].append(item)

    for group in groups:
        copies = group["copies"]
        group["retained_atoms"] = sorted({
            atom for item in copies for atom in item["retained_atoms"]
        })
        group["symmetry_retained_atoms"] = sorted({
            atom for item in copies for atom in item["symmetry_retained_atoms"]
        })
        group["covered_target_atoms"] = sorted({
            atom for item in copies for atom in item["covered_target_atoms"]
        })
        group["boundary_bonds"] = sorted({
            tuple(sorted(bond))
            for item in copies for bond in item["boundary_bonds"]
        })
        group["leftover_atoms"] = sum(
            len(fragment)
            for item in copies for fragment in item["leftover_fragments"]
        )
        group["complete"] = all(item["complete"] for item in copies)
        group["multiplicity"] = len(copies)
        group["mapping"] = [
            pair for item in copies for pair in item["mapping"]
        ]
        group["symmetry_target_atoms"] = sorted({
            int(target)
            for item in copies
            for _source, targets in item["target_domains"]
            for target in targets
        })
    return groups


def _matches_expected(assembly, known_ids, expected_coverage):
    if {item["precursor_id"] for item in assembly["precursors"]} != known_ids:
        return False
    covered_by_id = {}
    for item in assembly["precursors"]:
        covered_by_id.setdefault(item["precursor_id"], set()).update(
            item["covered_target_atoms"])
    return all(
        len(covered_by_id.get(precursor_id, ())) == atom_count
        for precursor_id, atom_count in expected_coverage.items()
    )


def _payload(
        report, top_count, known_ids, expected_coverage, title, known_label,
        ground_truth_status="not-evaluated",
        ground_truth_note="No ground-truth metadata supplied."):
    target_smiles = report["target_smiles"]
    cache = {}
    target_block, target_coords, target_elements = mol_3d(
        target_smiles, show_hydrogens=True)
    ranked = report["assemblies"]
    known_rank = next(
         (index for index, assembly in enumerate(ranked, 1)
         if _matches_expected(assembly, known_ids, expected_coverage)),
        None,
    )
    selected_ranks = list(range(1, min(top_count, len(ranked)) + 1))
    if known_rank is not None and known_rank not in selected_ranks:
        selected_ranks.append(known_rank)

    selected = []
    if report.get("expected_assembly") is not None:
        selected.append(("ground truth", report["expected_assembly"], True))
    selected.extend(
        (rank, ranked[rank - 1], False) for rank in selected_ranks)

    expected_counts = Counter(report.get("expected_ids", ()))
    expected_found = report.get("expected_ids_found", {})

    assemblies = []
    for rank, assembly, ground_truth in selected:
        raw_precursors = assembly["precursors"]
        regions = [
            set(map(int, item["covered_target_atoms"]))
            for item in raw_precursors
        ]
        if set().union(*regions) != set(range(len(target_elements))):
            raise ValueError(
                f"assembly {rank} does not cover the full target")
        precursors = _group_precursors(raw_precursors)
        models = [
            _model(item["smiles"], item, COLORS[index], cache)
            for index, item in enumerate(precursors)
        ]
        displayed = set()
        product_styles = []
        for index, item in enumerate(precursors):
            owned = sorted(
                set(item["covered_target_atoms"]) - displayed)
            displayed.update(owned)
            product_styles.append({
                "indices": owned,
                "color": COLORS[index],
            })
        product_symmetry_styles = [
            {"indices": item["symmetry_target_atoms"], "color": COLORS[index]}
            for index, item in enumerate(precursors)
        ]
        models.append({
            "mol": target_block,
            "coords": target_coords,
            "elements": target_elements,
            "styles": product_styles,
            "symmetryStyles": product_symmetry_styles,
            "labels": [
                {"atom": atom, "text": f"P{atom}"}
                for atom in range(len(target_elements))
            ],
            "broken": [],
            "formed": assembly["formed_bonds"],
        })
        pattern = assembly.get("construction_pattern", "GT")
        assemblies.append({
            "rank": rank,
            "pattern": pattern,
            "known": ground_truth or _matches_expected(
                assembly, known_ids, expected_coverage),
            "ground_truth": ground_truth,
            "score": assembly["score"],
            "precursors": [{
                "id": item["precursor_id"],
                "smiles": item["smiles"],
                "retained": item["retained_atoms"],
                "unmatched": item["leftover_atoms"],
                "complete": item["complete"],
                "multiplicity": item["multiplicity"],
                "symmetry_positions": len(item["symmetry_target_atoms"]),
            } for item in precursors],
            "models": models,
        })
    return {
        "summary": {
            "title": title,
            "known_label": known_label,
            "catalog_rows": report["scan_counts"]["rows"],
            "searched": report["scan_counts"]["searched"],
            "matched_precursors": report["scan_counts"]["matched_precursors"],
            "fragment_candidates": report["scan_counts"]["fragment_candidates"],
            "capped": report["scan_counts"]["capped"],
            "assemblies": len(ranked),
            "known_rank": known_rank,
            "ground_truth_status": ground_truth_status,
            "ground_truth_note": ground_truth_note,
            "ground_truth_reactants": [
                {
                    "id": precursor_id,
                    "multiplicity": multiplicity,
                    "detected": bool(expected_found.get(precursor_id)),
                }
                for precursor_id, multiplicity in expected_counts.items()
            ],
            "explicit_hydrogens": True,
            "search_truncated": report["recommendation_search_truncated"],
        },
        "patterns": report["construction_patterns"],
        "assemblies": assemblies,
    }


def _html(payload):
    library = (Path(__file__).parents[1] / "src" / "rxn_core" / "static" /
               "3Dmol-min.js").read_text()
    data = json.dumps(payload, separators=(",", ":"))
    reactants = payload["summary"]["ground_truth_reactants"]
    ground_truth_reactants = " + ".join(
        escape(item["id"])
        + (f" ×{item['multiplicity']}" if item["multiplicity"] > 1 else "")
        + (" ✓ detected" if item["detected"] else " ✗ absent")
        for item in reactants
    ) or "No structured ground-truth reactants supplied"
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>Blind catalog retrosynthesis results</title><style>
:root{{--bg:#f3f5f8;--card:#fff;--ink:#172033;--muted:#64748b;--line:#dbe2ea;--blue:#2684ff;--orange:#ff8b00}}
*{{box-sizing:border-box}} body{{margin:0;height:100vh;overflow:hidden;font:13px system-ui;background:var(--bg);color:var(--ink)}}
header{{height:72px;padding:11px 18px;background:#101828;color:white;display:flex;align-items:center;justify-content:space-between}}
h1{{font-size:18px;margin:0 0 4px}} .muted{{color:#9fb0c8}} .metrics{{display:flex;gap:18px}} .metric b{{display:block;font-size:17px}} .metric small{{color:#9fb0c8}}
#layout{{display:grid;grid-template-columns:330px 1fr;height:calc(100vh - 72px)}}
aside{{background:white;border-right:1px solid var(--line);overflow:auto}} .intro{{padding:12px;border-bottom:1px solid var(--line);line-height:1.45}}
.result{{padding:10px 12px;border:0;border-bottom:1px solid var(--line);width:100%;text-align:left;background:white;cursor:pointer}}
.result:hover{{background:#f5f8fc}} .result.active{{background:#eaf3ff;box-shadow:inset 4px 0 var(--blue)}}
.patternhead{{padding:9px 12px;background:#e2e8f0;border-top:2px solid #94a3b8;border-bottom:1px solid var(--line);font-weight:800}}
.patternhead small{{display:block;color:#475569;font-weight:500;margin-top:2px}}
.rank{{font-weight:750;font-size:14px}} .badge{{background:#0f9d66;color:white;border-radius:10px;padding:2px 7px;margin-left:7px;font-size:10px}}
.truthbox{{margin-top:10px;padding:8px;border-radius:7px;background:#ecfdf3;border:1px solid #6ee7a8;color:#166534}} .truthbox.partial{{background:#fffbeb;border-color:#fbbf24;color:#92400e}} .truthbox.missing{{background:#fef2f2;border-color:#fca5a5;color:#991b1b}} .truthreactants{{display:block;margin-top:7px;padding-top:7px;border-top:1px solid currentColor;line-height:1.5}} .patternbadge{{background:#475569;color:white;border-radius:10px;padding:2px 7px;margin-left:7px;font-size:10px}}
.patternbadge{{background:#475569;color:white;border-radius:10px;padding:2px 7px;margin-left:7px;font-size:10px}}
.ids{{color:#475569;margin-top:4px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}} .score{{color:#64748b;font-size:11px;margin-top:4px}}
main{{display:grid;grid-template-rows:1fr 1fr;min-width:0}} #reactants{{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:8px;padding:8px 8px 4px}}
.panel{{position:relative;background:var(--card);border:1px solid var(--line);border-radius:9px;overflow:hidden;min-height:220px}}
#productWrap{{padding:4px 8px 8px}} #Ppanel{{height:100%}} .view{{position:absolute;inset:0}}
.label{{position:absolute;z-index:3;left:10px;top:9px;max-width:80%;background:#ffffffdf;border:1px solid var(--line);border-radius:7px;padding:6px 9px;pointer-events:none}}
.label b{{display:block}} .label small{{display:block;color:#64748b;margin-top:2px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}
.controls{{position:absolute;z-index:5;right:18px;top:82px;background:#ffffffdc;padding:5px 8px;border:1px solid var(--line);border-radius:7px}}
.dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin:0 5px 0 10px}} code{{font-size:11px}}
</style><script>{library}</script></head><body>
<header><div><h1>{payload['summary']['title']}</h1><div class="muted">Explicit-H fragment mappings · select an assembly to inspect it in 3D</div></div>
<div class="metrics"><span class="metric"><b>{payload['summary']['catalog_rows']:,}</b><small>catalog R</small></span><span class="metric"><b>{payload['summary']['matched_precursors']:,}</b><small>matched R</small></span><span class="metric"><b>{payload['summary']['fragment_candidates']:,}</b><small>fragments</small></span><span class="metric"><b>{payload['summary']['assemblies']}</b><small>ranked assemblies</small></span></div></header>
<div id="layout"><aside><div class="intro"><b>Each color is one unique precursor.</b><br>Repeated copies share one color and one R panel. Hydrogens are explicit. Unmatched atoms keep element colors. Symmetry mode colors every R position in the retained source orbits and every P position allowed by the compressed target domains. These are alternative matchable positions, not extra simultaneous assignments. <span style="color:#d33">Red = broken</span>; <span style="color:#159447">green = formed</span>.<div class="truthbox {payload['summary']['ground_truth_status']}"><b>Ground truth: {payload['summary']['ground_truth_status'].upper()}</b><br>{payload['summary']['ground_truth_note']}<span class="truthreactants"><b>Ground-truth raw ingredients</b><br>{ground_truth_reactants}</span></div><br>Cap-hit precursors: {payload['summary']['capped']:,}. <span style="color:#b45309;font-weight:700">Assembly search truncated: {'yes' if payload['summary']['search_truncated'] else 'no'}.</span></div><div id="list"></div></aside>
<main><div class="controls"><label><input id="fragments" type="checkbox" checked> color fragments</label><br><label><input id="symmetry" type="checkbox"> show symmetry domains</label><br><label><input id="labels" type="checkbox"> sampled P# identities</label></div><div id="reactants"></div>
<div id="productWrap"><section class="panel" id="Ppanel"><div class="label" id="LP"></div><div class="view" id="P"></div></section></div></main></div>
<script>const data={data}, colors={json.dumps(COLORS)}, viewers={{}};
function pt(m,i){{return {{x:m.coords[i][0],y:m.coords[i][1],z:m.coords[i][2]}}}}
function showModel(id,m){{let v=viewers[id];if(!v){{v=$3Dmol.createViewer(id,{{backgroundColor:'white'}});viewers[id]=v}}else{{v.removeAllModels();v.removeAllShapes();v.removeAllLabels()}}
 v.addModel(m.mol,'sdf');v.setStyle({{}},{{stick:{{radius:.12}},sphere:{{scale:.23}}}});if(document.getElementById('fragments').checked){{const symmetry=document.getElementById('symmetry').checked;const styles=symmetry&&(m.symmetryStyles||[]).length?m.symmetryStyles:m.styles;styles.forEach(s=>v.addStyle({{index:s.indices}},{{stick:{{color:s.color,radius:.19}},sphere:{{color:s.color,scale:.34}}}}))}}
 m.broken.forEach(b=>v.addCylinder({{start:pt(m,b[0]),end:pt(m,b[1]),radius:.10,color:'#e5484d'}}));m.formed.forEach(b=>v.addCylinder({{start:pt(m,b[0]),end:pt(m,b[1]),radius:.11,color:'#16a34a'}}));
 if(document.getElementById('labels').checked)(m.labels||[]).forEach(l=>v.addLabel(l.text,{{position:pt(m,l.atom),fontSize:10,fontColor:'#111',backgroundColor:'white',backgroundOpacity:.72,inFront:true}}));v.zoomTo();v.render()}}
function patternInfo(id){{return data.patterns.find(x=>x.pattern===id)}}
function patternText(id){{if(id==='GT')return 'independently evaluated known reactants';const p=patternInfo(id);return p?(p.fragment_sizes.length+' modules · atom sizes '+p.fragment_sizes.join(' + ')):''}}
function select(i){{const a=data.assemblies[i];document.querySelectorAll('.result').forEach(x=>x.classList.toggle('active',Number(x.dataset.index)===i));Object.keys(viewers).filter(k=>k.startsWith('R')).forEach(k=>delete viewers[k]);const wrap=document.getElementById('reactants');wrap.innerHTML='';a.precursors.forEach((r,j)=>{{const panel=document.createElement('section');panel.className='panel';panel.innerHTML='<div class="label" id="L'+j+'"></div><div class="view" id="R'+j+'"></div>';wrap.appendChild(panel);const mult=r.multiplicity>1?' ×'+r.multiplicity:'';document.getElementById('L'+j).innerHTML='<b><span style="color:'+colors[j%colors.length]+'">R'+(j+1)+'</span> · '+r.id+mult+'</b><small>'+r.smiles+'</small><small>retained '+r.retained.length+' atom positions; unmatched '+r.unmatched+' across copies</small>';showModel('R'+j,a.models[j])}});
 const retention=a.score.set_atom_retention===undefined?'':(' · direct retention '+(100*a.score.set_atom_retention).toFixed(1)+'%');const symmetryRetention=a.score.set_symmetry_atom_retention===undefined?'':(' · symmetry-adjusted '+(100*a.score.set_symmetry_atom_retention).toFixed(1)+'%');const chiral=a.score.chirality_violations===undefined?'':(' · chirality violations '+a.score.chirality_violations);const heading=a.ground_truth?'GROUND TRUTH':'P target · Pattern '+a.pattern;document.getElementById('LP').innerHTML='<b>'+heading+': '+patternText(a.pattern)+'</b><small>'+a.rank+' · '+a.score.broken_bonds+' broken, '+a.score.leftover_atoms+' unmatched, '+a.score.formed_bonds+' formed'+retention+symmetryRetention+chiral+'</small>';showModel('P',a.models[a.precursors.length])}}
const list=document.getElementById('list');let lastPattern=null;data.assemblies.forEach((a,i)=>{{if(a.pattern!==lastPattern){{const h=document.createElement('div');h.className='patternhead';h.innerHTML=(a.ground_truth?'GROUND TRUTH':'Pattern '+a.pattern)+'<small>'+patternText(a.pattern)+(a.ground_truth?'':' · colored regions on P define this construction')+'</small>';list.appendChild(h);lastPattern=a.pattern}}const b=document.createElement('button');b.className='result';b.dataset.index=i;b.innerHTML='<span class="rank">'+(a.ground_truth?'ground truth':'recommendation '+a.rank)+'</span>'+(a.known?'<span class="badge">GROUND TRUTH</span>':'')+'<div class="ids">'+a.precursors.map(x=>x.id+(x.multiplicity>1?' ×'+x.multiplicity:'')).join(' + ')+'</div><div class="score">direct retention '+(a.score.set_atom_retention===undefined?'n/a':(100*a.score.set_atom_retention).toFixed(1)+'%')+(a.score.set_symmetry_atom_retention===undefined?'':(' · symmetry-adjusted '+(100*a.score.set_symmetry_atom_retention).toFixed(1)+'%'))+' · broken '+a.score.broken_bonds+' · unmatched atoms '+a.score.leftover_atoms+' · formed '+a.score.formed_bonds+(a.score.chirality_violations===undefined?'':' · chirality '+a.score.chirality_violations)+'</div>';b.onclick=()=>select(i);list.appendChild(b)}});
function redraw(){{select([...document.querySelectorAll('.result')].findIndex(x=>x.classList.contains('active')))}}document.getElementById('labels').onchange=redraw;document.getElementById('fragments').onchange=redraw;document.getElementById('symmetry').onchange=redraw;if(data.assemblies.length){{const knownIndex=data.assemblies.findIndex(x=>x.known);select(knownIndex>=0?knownIndex:0)}}window.onresize=()=>Object.values(viewers).forEach(v=>v.resize());</script></body></html>"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top", type=int, default=24)
    parser.add_argument("--expected-id", action="append")
    parser.add_argument(
        "--expected-coverage", action="append", default=[], metavar="ID=ATOMS",
        help="require the selected known assembly to cover ATOMS target atoms with ID")
    parser.add_argument("--title", default="2-chlorobiphenyl")
    parser.add_argument("--known-label", default="Suzuki")
    parser.add_argument("--ground-truth-status", default="not evaluated")
    parser.add_argument("--ground-truth-note", default="No ground-truth metadata supplied.")
    args = parser.parse_args()
    report = json.loads(Path(args.results).read_text())
    known_ids = set(args.expected_id or ())
    expected_coverage = {}
    for value in args.expected_coverage:
        precursor_id, separator, atom_count = value.rpartition("=")
        if not separator or not precursor_id:
            parser.error(f"invalid --expected-coverage value: {value!r}")
        expected_coverage[precursor_id] = int(atom_count)
    payload = _payload(
        report, args.top, known_ids, expected_coverage,
        args.title, args.known_label, args.ground_truth_status,
        args.ground_truth_note)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_html(payload))
    print(output.resolve())


if __name__ == "__main__":
    main()

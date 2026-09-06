#!/usr/bin/env python3
"""Build an interactive 3D viewer from merged blind catalog-search results."""
from __future__ import annotations

import argparse
from collections import Counter
from html import escape
import json
import colorsys
from pathlib import Path

from molecule_3d import mol_3d


COLORS = (
    "#2684ff", "#ff8b00", "#8b5cf6", "#00a896", "#e83e8c",
    "#795548", "#00acc1", "#7cb342", "#f4511e", "#5c6bc0",
    "#c0a000", "#6d4c91",
)
def _color(index):
    if index < len(COLORS):
        return COLORS[index]
    rgb = colorsys.hsv_to_rgb((index * 0.61803398875) % 1, 0.7, 0.85)
    return "#" + "".join(f"{round(channel * 255):02x}" for channel in rgb)
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
        } | {int(atom) for item in copies for region in item.get("allowed_target_regions", ())
             for atom in region})
    return groups


def _matches_expected(assembly, expected_ids):
    return Counter(
        item["precursor_id"] for item in assembly["precursors"]
    ) == Counter(expected_ids)


def _fragment_colors(precursors, models):
    """Color saved matched partitions, never infer fragments from geometry.

    Identical source-atom fragments across copies share a color. Their separate
    mapped target occupations are retained. Conflicting colors on a merged R
    or overlapping P are explicit display alternatives, not new assignments.
    """
    fragments = []
    for source_index, group in enumerate(precursors):
        by_atoms = {}
        for copy_index, copy in enumerate(group['copies']):
            mapping = dict(copy['mapping'])
            parts = copy['retained_fragments']
            atoms = [a for part in parts for a in part]
            if len(atoms) != len(set(atoms)) or set(atoms) != set(mapping):
                raise ValueError('Saved matched fragments must partition the retained mapping')
            for part in parts:
                key = tuple(sorted(part))
                if key not in by_atoms:
                    index = len(fragments)
                    fragment = dict(index=index, source=source_index,
                        label=f'R{source_index + 1}.F{len(by_atoms) + 1}',
                        color=_color(index), source_atoms=list(key), occupations=[])
                    by_atoms[key] = fragment
                    fragments.append(fragment)
                by_atoms[key]['occupations'].append(dict(copy=copy_index + 1,
                    mapping=[[a, mapping[a]] for a in key]))
    for model_index, model in enumerate(models):
        target = model_index == len(precursors)
        claims = {}
        for fragment in fragments:
            if not target and fragment['source'] != model_index:
                continue
            atoms = (set(b for occupation in fragment['occupations']
                         for a, b in occupation['mapping']) if target else fragment['source_atoms'])
            for atom in atoms:
                claims.setdefault(atom, set()).add(fragment['index'])
        styles = []
        alternatives = {}
        for fragment in fragments:
            owned = sorted(a for a, owners in claims.items() if owners == {fragment['index']})
            if owned:
                styles.append(dict(indices=owned, color=fragment['color']))
        for atom, owners in sorted(claims.items()):
            if len(owners) > 1:
                alternatives.setdefault(tuple(sorted(owners)), []).append(atom)
        model['fragmentStyles'] = styles
        model['fragmentAlternatives'] = [dict(atoms=atoms, owners=list(owners), selected=owners[0])
                                        for owners, atoms in alternatives.items()]
    return fragments


def _payload(
        report, top_count, title,
        ground_truth_status="not-evaluated",
        ground_truth_note="No ground-truth metadata supplied."):
    target_smiles = report["target_smiles"]
    cache = {}
    target_block, target_coords, target_elements = mol_3d(
        target_smiles, show_hydrogens=True)
    ranked = report["assemblies"]
    expected_ids = report.get("expected_ids", ())
    known_index = next(
        (index for index, assembly in enumerate(ranked)
         if _matches_expected(assembly, expected_ids)), None)
    known_rank = report.get("expected_recommendation_rank")
    if known_rank is None and known_index is not None:
        known_rank = known_index + 1
    selected_ranks = list(range(1, min(top_count, len(ranked)) + 1))
    if known_index is not None and known_index + 1 not in selected_ranks:
        selected_ranks.append(known_index + 1)

    selected = [
        (rank, ranked[rank - 1], False) for rank in selected_ranks]
    selected.extend((f"validation {i}",assembly,True)
                    for i,assembly in enumerate(report.get('validation_assemblies',()),1))
    if known_rank is None and report.get("expected_mapping") is not None:
        selected.append(("ground truth", report["expected_mapping"], True))
    if report.get("diagnostic_assembly") is not None:
        selected.append(("diagnostic", report["diagnostic_assembly"], True))

    expected_counts = Counter(report.get("expected_ids", ()))
    expected_found = report.get("expected_ids_found", {})

    assemblies = []
    for rank, assembly, ground_truth in selected:
        raw_precursors = assembly["precursors"]
        regions = [
            set(map(int, item["covered_target_atoms"]))
            for item in raw_precursors
        ]
        complete_cover = set().union(*regions) == set(range(len(target_elements)))
        if not ground_truth and not complete_cover:
            raise ValueError(
                f"assembly {rank} does not cover the full target")
        precursors = _group_precursors(raw_precursors)
        models = [
            _model(item["smiles"], item, _color(index), cache)
            for index, item in enumerate(precursors)
        ]
        claims = {atom: [index for index, item in enumerate(precursors)
                         for copy in item["copies"] if atom in copy["covered_target_atoms"]]
                  for atom in range(len(target_elements))}
        owners = {atom: sorted(set(indices)) for atom, indices in claims.items()}
        alternatives = {}
        for atom, indices in owners.items():
            if len(indices) > 1:
                alternatives.setdefault(tuple(indices), []).append(atom)
        product_styles = []
        for index, item in enumerate(precursors):
            owned = sorted(atom for atom, indices in owners.items() if indices == [index])
            product_styles.append({
                "indices": owned,
                "color": _color(index),
            })
        product_symmetry_styles = [
            {"indices": item["symmetry_target_atoms"], "color": _color(index)}
            for index, item in enumerate(precursors)
        ]
        models.append({
            "mol": target_block,
            "coords": target_coords,
            "elements": target_elements,
            "styles": product_styles,
            "symmetryStyles": product_symmetry_styles,
            "supplierAlternatives": [dict(atoms=atoms, owners=list(indices), selected=indices[0])
                                     for indices, atoms in alternatives.items()],
            "labels": [
                {"atom": atom, "text": (f"P{atom} UNCOVERED" if not claims[atom] else
                    f"P{atom}" if len(claims[atom]) == 1 else
                    f"P{atom}: " + " or ".join(f"R{i + 1}" +
                        (f" ({claims[atom].count(i)} copies)" if claims[atom].count(i) > 1 else "")
                        for i in owners[atom])),
                 "always": len(claims[atom]) != 1}
                for atom in range(len(target_elements))
            ],
            "broken": [],
            "formed": assembly["formed_bonds"],
        })
        pattern = assembly.get("construction_pattern", "GT")
        fragments = _fragment_colors(precursors, models)
        assemblies.append({
            "fragments": fragments,
            "rank": rank,
            "pareto_layer": assembly.get('pareto_layer'),
            "pattern": pattern,
            "known": rank != "diagnostic" and (ground_truth or _matches_expected(
                assembly, expected_ids)),
            "ground_truth": ground_truth and rank != "diagnostic",
            "diagnostic": rank == "diagnostic",
            "complete_cover": complete_cover,
            "score": dict(assembly["score"],
                          covered_target_atoms=len(set().union(*regions)),
                          target_atom_count=len(target_elements)),
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
            "search_scope": report.get("search_scope", "Blind catalog search"),
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
        "unassembled_target": {
            "mol": target_block, "coords": target_coords, "elements": target_elements,
            "styles": [], "broken": [], "formed": [],
            "fragmentStyles": [], "fragmentAlternatives": [],
            "labels": [{"atom": atom, "text": f"P{atom}"}
                       for atom in range(len(target_elements))],
        } if not assemblies else None,
        "uncovered_target_atoms": report.get("uncovered_target_atoms", []),
    }


def _html(payload):
    score_plot = (Path(__file__).parents[1] / 'src/rxn_core/static/retro_score_plot.js').read_text()
    fragment_colors = (Path(__file__).parents[1] / 'src/rxn_core/static/retro_fragment_colors.js').read_text()
    library = (Path(__file__).parents[1] / "src" / "rxn_core" / "static" /
               "3Dmol-min.js").read_text()
    data = json.dumps(payload, separators=(",", ":"))
    reactants = payload["summary"]["ground_truth_reactants"]
    ground_truth_reactants = " + ".join(
        escape(item["id"])
        + (f" ×{item['multiplicity']}" if item["multiplicity"] > 1 else "")
        + (" ✓ detected" if item["detected"] else " ✗ no saved detection (bank absence not established)")
        for item in reactants
    ) or "No structured ground-truth reactants supplied"
    known_rank = payload["summary"]["known_rank"]
    recommendation_truth = (
        f"yes, candidate rank {known_rank} in the searched pool"
        if known_rank is not None else "no"
    )
    palette = [_color(i) for i in range(max(
        (len(a["precursors"]) for a in payload["assemblies"]), default=0))]
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>Geometric building-block results</title><style>
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
.truthbox.incomplete{{background:#fffbeb;border-color:#fbbf24;color:#92400e}} .truthbox.not-evaluated{{background:#f8fafc;border-color:#94a3b8;color:#334155}}
.ids{{color:#475569;margin-top:4px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}} .score{{color:#64748b;font-size:11px;margin-top:4px}}
.scoregrid{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:5px;margin-top:8px}} .scorecell{{background:#f1f5f9;border:1px solid #dbe2ea;border-radius:6px;padding:6px;min-width:0}} .scorecell b{{display:block;font-size:18px}} .scorecell small{{display:block;font-size:10px;color:#475569}} .scorecell.breaking b{{color:#c62828}} .scorecell.forming b{{color:#148343}}
#selectedScores .scoregrid{{grid-template-columns:repeat(6,minmax(0,1fr));margin:0 0 6px}} #selectedScores .scorecell b{{font-size:22px}} .sortbar{{padding:10px 12px;background:white;position:sticky;top:0;z-index:6;border-bottom:1px solid var(--line)}} .sortbar select{{width:100%;margin-top:4px;padding:6px}} .sortbar small{{display:block;color:#64748b;margin-top:4px}}
main{{display:grid;grid-template-rows:auto minmax(0,1fr);min-width:0;min-height:0;overflow:auto}}
#moleculeWorkspace{{display:grid;grid-template-columns:minmax(260px, .38fr) minmax(0, .62fr);min-height:540px;gap:8px;padding:0 8px 8px}}
#reactants{{display:flex;flex-direction:column;gap:8px;overflow:auto;min-height:0}} #reactants .panel{{flex:0 0 300px}}
#moleculeWorkspace.target-focus{{grid-template-columns:1fr}} #moleculeWorkspace.target-focus #reactants{{display:none}}
#scorePlotPanel{{margin:6px 8px;padding:6px 10px;background:white;border:1px solid var(--line);border-radius:8px}} #scorePlot{{display:block;width:100%;height:180px}} #scorePlotChoices{{display:flex;gap:5px;align-items:center;overflow:auto;white-space:nowrap;min-height:30px}} #scorePlotChoices button{{padding:4px 7px;background:white;border:1px solid #94a3b8;border-radius:5px;cursor:pointer}} #scorePlotChoices .chosen{{border:2px solid #2684ff;background:#eaf3ff}} .plotnote{{font-size:11px;color:#475569}}
.panel{{position:relative;background:var(--card);border:1px solid var(--line);border-radius:9px;overflow:hidden;min-height:220px}}
#productWrap{{display:grid;grid-template-rows:auto auto minmax(0,1fr);min-height:0;min-width:0}} #Ppanel{{height:100%;min-height:380px;display:flex;flex-direction:column}} .view{{position:absolute;inset:0}}
#Ppanel .label{{position:static;max-width:none;border:0;border-bottom:1px solid var(--line);border-radius:0;background:white;flex:none}}
#Ppanel .label small{{white-space:normal}} #P{{position:relative;flex:1;min-height:320px}}
#selectedScores .scoregrid{{grid-template-columns:repeat(4,minmax(0,1fr));gap:4px}} #selectedScores .scorecell{{padding:4px}} #selectedScores .scorecell b{{font-size:18px}}
.target-tools{{display:flex;gap:8px;align-items:center;padding:6px 0}} .target-tools button{{padding:6px 10px;border:1px solid #94a3b8;border-radius:6px;background:white;cursor:pointer}} .target-tools span{{color:var(--muted);font-size:11px}}
#fragmentLegend{{display:flex;flex-wrap:wrap;gap:5px;max-height:82px;overflow:auto;padding:4px 0}} #fragmentLegend>span:first-child{{width:100%;font-size:11px;color:var(--muted)}} .fragment-chip{{display:inline-flex;align-items:center;gap:5px;border:1px solid var(--line);border-radius:5px;padding:3px 6px;font-size:11px}} .fragment-chip i{{width:12px;height:12px;border-radius:3px;display:inline-block}}
#scorePlotPanel summary{{cursor:pointer;font-weight:700;padding:4px 0}} #scorePlotPanel .controls{{margin-bottom:2px}}
#scorePlotPanel,#moleculeWorkspace,#reactants{{min-width:0}} #scorePlotChoices{{max-width:100%;min-width:0}} header>div:first-child{{min-width:0}} header .muted{{max-height:32px;overflow:hidden}}
@media(max-width:1000px){{.metrics{{display:none}} #layout{{grid-template-columns:270px minmax(0,1fr)}} #moleculeWorkspace{{grid-template-columns:minmax(0,1fr);min-height:940px;grid-template-rows:minmax(620px,1fr) 300px}} #productWrap{{grid-row:1}} #reactants{{flex-direction:row;grid-row:2}} #reactants .panel{{flex:0 0 280px}} #moleculeWorkspace.target-focus{{grid-template-rows:1fr}}}}
@media(max-width:650px){{header{{height:auto;min-height:72px}} .metrics{{display:none}} #layout{{grid-template-columns:1fr;height:calc(100dvh - 72px);overflow:auto}} aside{{max-height:220px}} main{{overflow:visible}} #moleculeWorkspace{{min-height:900px}}}}
.label{{position:absolute;z-index:3;left:10px;top:9px;max-width:80%;background:#ffffffdf;border:1px solid var(--line);border-radius:7px;padding:6px 9px;pointer-events:none}}
.label b{{display:block}} .label small{{display:block;color:#64748b;margin-top:2px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}
.controls{{display:flex;gap:12px;flex-wrap:wrap;background:#ffffffdc;padding:4px 0}}
.dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin:0 5px 0 10px}} code{{font-size:11px}}
</style><script>{library}</script><script>{score_plot}</script><script>{fragment_colors}</script></head><body>
<header><div><h1>{payload['summary']['title']}</h1><div class="muted">{escape(payload['summary']['search_scope'])} · Explicit-H fragment mappings</div></div>
<div class="metrics"><span class="metric"><b>{payload['summary']['catalog_rows']:,}</b><small>catalog R</small></span><span class="metric"><b>{payload['summary']['matched_precursors']:,}</b><small>matched R</small></span><span class="metric"><b>{payload['summary']['fragment_candidates']:,}</b><small>fragments</small></span><span class="metric"><b>{payload['summary']['assemblies']}</b><small>ranked assemblies</small></span></div></header>
<div id="layout"><aside><div class="intro"><b>Each color is one matched fragment.</b><br>
R1.F1 and R1.F2 are different fragments of R1. The same color marks a fragment in R and its mapped atoms in P.
Repeated copies share one R panel; identical source-atom fragments share their color across copies. Hydrogens are explicit. Unmatched atoms keep element colors.
Colors are local to the selected assembly. For overlapping fragment colors,
use the buttons below to switch the displayed color; the legend lists the fragment assignments.
These are display alternatives, not independently validated new assemblies.
Symmetry mode uses precursor colors and shows alternative matchable positions, not extra simultaneous fragment assignments.<br>
<span style="color:#d33">Red = source cuts</span>; <span style="color:#159447">green = unsupported target connections</span>.
These are geometric connections, not validated reaction edits.
<div class="truthbox {payload['summary']['ground_truth_status']}"><b>Ground truth: {payload['summary']['ground_truth_status'].upper()}</b><br><b>Returned by recommender: {recommendation_truth}</b><br>{payload['summary']['ground_truth_note']}<span class="truthreactants"><b>Ground-truth raw ingredients</b><br>{ground_truth_reactants}</span></div><br>Cap-hit precursors: {payload['summary']['capped']:,}. <span style="color:#b45309;font-weight:700">Assembly search truncated: {'yes' if payload['summary']['search_truncated'] else 'no'}.</span></div><div id="suppliers" class="intro"></div><div id="list"></div></aside>
<main><section id="scorePlotPanel"><div class="controls"><label><input id="fragments" type="checkbox" checked> color matches</label><label>Color by <select id="colorBy"><option value="fragment">Matched fragment</option><option value="precursor">Precursor</option></select></label><label><input id="symmetry" type="checkbox"> symmetry domains (precursor colors)</label><label><input id="labels" type="checkbox"> sampled P# identities</label></div><div id="fragmentLegend"></div><details id="scoreDetails"><summary>Score plot · click to compare assemblies</summary><div class="plotnote">Blue: blind proposals · green: validation · ×N: coincident alternatives (no jitter). Upper left is better; trade-offs share a Pareto layer.</div><svg id="scorePlot" aria-label="Retention versus structural changes"></svg><div id="scorePlotChoices"></div></details></section><div id="moleculeWorkspace"><div id="reactants"></div>
<div id="productWrap"><div id="selectedScores"></div><div class="target-tools"><button id="fitTarget" type="button">Fit target</button><button id="expandTarget" type="button" aria-pressed="false">Expand target</button><span>Drag to rotate · scroll to zoom</span></div><section class="panel" id="Ppanel"><div class="label" id="LP"></div><div class="view" id="P"></div></section></div></div></main></div>
<script>const data={data}, colors={json.dumps(palette)}, viewers={{}};let activeAssembly=null;
function pt(m,i){{return {{x:m.coords[i][0],y:m.coords[i][1],z:m.coords[i][2]}}}}
function scoreCards(a){{const s=a.score;const cell=(label,value,detail,cls='')=>'<div class="scorecell '+cls+'"><small>'+label+'</small><b>'+value+'</b><small>'+detail+'</small></div>';return '<div class="scoregrid">'+cell('Pareto rank',a.ground_truth?'—':(a.pareto_layer||'—'),a.ground_truth?'validation only':'same rank = trade-offs')+cell('Bond breaking',s.broken_bonds===undefined?'—':s.broken_bonds,'source boundary cuts','breaking')+cell('Bond forming',s.formed_bonds===undefined?'—':s.formed_bonds,'target connections','forming')+cell('Total changes',s.broken_bonds===undefined||s.formed_bonds===undefined?'—':s.broken_bonds+s.formed_bonds,'breaking + forming')+cell('P coverage',(100*s.covered_target_atoms/s.target_atom_count).toFixed(1)+'%',s.covered_target_atoms+' / '+s.target_atom_count+' atoms incl. H')+cell('R retention',s.set_atom_retention===undefined?'—':(100*s.set_atom_retention).toFixed(1)+'%','unique P / all input atoms')+cell('Matched fragments',s.matched_fragment_count===undefined?'—':s.matched_fragment_count,'equal-score tie-break')+'</div>'}}
function supplierControls(m){{const wrap=document.getElementById('suppliers');wrap.replaceChildren();if(fragmentMode()&&activeAssembly){{fragmentSupplierControls(activeAssembly,wrap,showModel);return}}const groups=m.supplierAlternatives||[];if(!groups.length){{wrap.textContent='No overlap between different R species. Repeated copies use their R color.';return}}const title=document.createElement('b');title.textContent='Alternative R suppliers · display only';wrap.appendChild(title);groups.forEach(g=>{{const row=document.createElement('div');const atoms=document.createElement('div');atoms.textContent=g.atoms.map(a=>'P'+a).join(', ');row.appendChild(atoms);g.owners.forEach(r=>{{const button=document.createElement('button');button.textContent='R'+(r+1);button.title='Show this supplier’s color; all listed suppliers remain in the saved assembly';button.style.color=colors[r%colors.length];button.style.fontWeight=g.selected===r?'bold':'normal';button.style.border=g.selected===r?'2px solid currentColor':'1px solid #ccc';button.disabled=document.getElementById('symmetry').checked;button.onclick=()=>{{g.selected=r;showModel('P',m)}};row.appendChild(button)}});wrap.appendChild(row)}})}}
function showModel(id,m){{let v=viewers[id];if(!v){{v=$3Dmol.createViewer(id,{{backgroundColor:'white'}});viewers[id]=v}}else{{v.removeAllModels();v.removeAllShapes();v.removeAllLabels()}}
 v.addModel(m.mol,'sdf');v.setStyle({{}},{{stick:{{radius:.12}},sphere:{{scale:.23}}}});if(document.getElementById('fragments').checked){{const symmetry=document.getElementById('symmetry').checked;const styles=fragmentMode()?m.fragmentStyles:symmetry&&(m.symmetryStyles||[]).length?m.symmetryStyles:m.styles;styles.forEach(s=>v.addStyle({{index:s.indices}},{{stick:{{color:s.color,radius:.19}},sphere:{{color:s.color,scale:.34}}}}));if(fragmentMode())m.fragmentAlternatives.forEach(g=>{{const color=activeAssembly.fragments[g.selected].color;v.addStyle({{index:g.atoms}},{{stick:{{color,radius:.19}},sphere:{{color,scale:.34}}}})}})}}
 if(id==='P'||fragmentMode())supplierControls(activeAssembly?activeAssembly.models[activeAssembly.precursors.length]:m);
 if(id==='P'&&!fragmentMode()){{if(document.getElementById('fragments').checked&&!document.getElementById('symmetry').checked)(m.supplierAlternatives||[]).forEach(g=>v.addStyle({{index:g.atoms}},{{stick:{{color:colors[g.selected%colors.length],radius:.19}},sphere:{{color:colors[g.selected%colors.length],scale:.34}}}}))}}
 m.broken.forEach(b=>v.addCylinder({{start:pt(m,b[0]),end:pt(m,b[1]),radius:.10,color:'#e5484d'}}));m.formed.forEach(b=>v.addCylinder({{start:pt(m,b[0]),end:pt(m,b[1]),radius:.11,color:'#16a34a'}}));
 (m.labels||[]).filter(l=>l.always||document.getElementById('labels').checked).forEach(l=>v.addLabel(l.text,{{position:pt(m,l.atom),fontSize:10,fontColor:l.always?'#b42318':'#111',backgroundColor:'white',backgroundOpacity:.85,inFront:true}}));v.zoomTo();v.render()}}
function patternInfo(id){{return data.patterns.find(x=>x.pattern===id)}}
function patternText(id){{if(String(id).startsWith('GT-'))return 'Validation pattern '+id+' · complete supplier-set assembly; separate from blind ranks';if(id==='GT')return 'ground-truth AAM comparison; not returned by blind recommender';const p=patternInfo(id);return p?(p.fragment_sizes.length+' modules · atom sizes '+p.fragment_sizes.join(' + ')):''}}
function select(i){{const a=data.assemblies[i];activeAssembly=a;fragmentLegend(a);renderScorePlot(data.assemblies,i,j=>{{select(j);document.querySelector('.result.active')?.scrollIntoView({{block:'nearest'}})}});document.getElementById('selectedScores').innerHTML=scoreCards(a);document.querySelectorAll('.result').forEach(x=>x.classList.toggle('active',Number(x.dataset.index)===i));Object.keys(viewers).filter(k=>k.startsWith('R')).forEach(k=>delete viewers[k]);const wrap=document.getElementById('reactants');wrap.innerHTML='';a.precursors.forEach((r,j)=>{{const panel=document.createElement('section');panel.className='panel';panel.innerHTML='<div class="label" id="L'+j+'"></div><div class="view" id="R'+j+'"></div>';wrap.appendChild(panel);const mult=r.multiplicity>1?' ×'+r.multiplicity:'';document.getElementById('L'+j).innerHTML='<b><span style="color:'+(fragmentMode()?'#172033':colors[j%colors.length])+'">R'+(j+1)+'</span> · '+r.id+mult+'</b><small>'+r.smiles+'</small><small>retained '+r.retained.length+' atom positions; unmatched '+r.unmatched+' across copies</small>';showModel('R'+j,a.models[j])}});
 const retention=a.score.set_atom_retention===undefined?'':(' · direct retention '+(100*a.score.set_atom_retention).toFixed(1)+'%');const symmetryRetention=a.score.set_symmetry_atom_retention===undefined?'':(' · symmetry-adjusted '+(100*a.score.set_symmetry_atom_retention).toFixed(1)+'%');const chiral=a.score.chirality_violations===undefined?'':(' · stereochemistry not assessed');const coverage=a.complete_cover?'complete P cover':(a.score.covered_target_atoms+' / '+a.score.target_atom_count+' P atoms covered');const heading=a.diagnostic?(a.complete_cover?'KNOWN-SET COVER · saved AAM mappings':'FAILED COVER · saved-mapping diagnostic'):a.ground_truth?'GROUND TRUTH MATCHING':'P target · Pattern '+a.pattern;document.getElementById('LP').innerHTML='<b>'+heading+'</b><small>'+coverage+' · '+a.score.broken_bonds+' source cuts, '+a.score.leftover_atoms+' unmatched, '+a.score.formed_bonds+' target connections (not reaction edits)'+retention+symmetryRetention+chiral+'</small>';showModel('P',a.models[a.precursors.length])}}
const list=document.getElementById('list');const sortbar=document.createElement('div');sortbar.className='sortbar';sortbar.innerHTML='<label>Sort displayed results by<select id="scoreSort"><option value="pareto">Pareto rank (default)</option value="retention">R retention: highest first</option value="changes">Total changes: fewest first</option><option value="breaking">Bond breaking: fewest first</option><option value="forming">Bond forming: fewest first</option><option value="coverage">P coverage: highest first</option></select></label><small>Sorting changes display order, not Pareto rank. Validation entries stay separate. Bond counts are geometric proxies, not validated reaction events.</small>';list.before(sortbar);
function sortedAssemblies(assemblies,mode){{
const value=a=>{{const s=a.score;switch(mode){{case 'retention':return -(s.set_atom_retention??-Infinity);case 'changes':return s.broken_bonds+s.formed_bonds;case 'breaking':return s.broken_bonds;case 'forming':return s.formed_bonds;case 'coverage':return -s.covered_target_atoms/s.target_atom_count;default:return a.pareto_layer??Infinity}}}};
return assemblies.map((a,i)=>({{a,i}})).sort((x,y)=>Number(x.a.ground_truth)-Number(y.a.ground_truth)||(value(x.a)-value(y.a))||x.i-y.i)}}
function renderList(){{list.replaceChildren();let lastPattern=null;const ordered=sortedAssemblies(data.assemblies,document.getElementById('scoreSort').value);
ordered.forEach(({{a,i}})=>{{
 if(a.pattern!==lastPattern){{
  const h=document.createElement('div');h.className='patternhead';
  h.innerHTML=a.diagnostic?(a.complete_cover?'KNOWN-SET COVER':'FAILED COVER · partial diagnostic')+'<small>One copy of each known ingredient. Actual saved AAM occupations; separate from blind recommendations.</small>':(a.ground_truth?'GROUND TRUTH MATCHING':'Pattern '+a.pattern)+'<small>'+patternText(a.pattern)+' · overlapping suppliers listed explicitly</small>';
  list.appendChild(h);lastPattern=a.pattern;
 }}
 const b=document.createElement('button');b.className='result';b.dataset.index=i;
 const label=a.diagnostic?(a.complete_cover?'Complete known-ingredient combination':'Best incomplete combination'):a.ground_truth?'ground-truth comparison':a.pareto_layer?'Pareto layer '+a.pareto_layer+' · alternative '+a.rank:'recommendation '+a.rank;
 b.innerHTML='<span class="rank">'+label+'</span>'+(a.known?'<span class="badge">GROUND TRUTH</span>':'')+'<div class="ids">'+a.precursors.map(x=>x.id+(x.multiplicity>1?' ×'+x.multiplicity:'')).join(' + ')+'</div>'+scoreCards(a);
 b.onclick=()=>select(i);list.appendChild(b);
}});
}}
renderList();document.getElementById('scoreSort').onchange=()=>{{const active=document.querySelector('.result.active');const index=active?Number(active.dataset.index):0;renderList();if(data.assemblies.length)select(index)}};
document.getElementById('colorBy').onchange=()=>redraw();
function fitTarget(){{const v=viewers.P;if(v){{v.resize();v.zoomTo();v.render()}}}}
document.getElementById('fitTarget').onclick=fitTarget;
document.getElementById('expandTarget').onclick=function(){{const expanded=document.getElementById('moleculeWorkspace').classList.toggle('target-focus');this.textContent=expanded?'Show reactants':'Expand target';this.setAttribute('aria-pressed',String(expanded));requestAnimationFrame(fitTarget)}};
new ResizeObserver(()=>{{Object.values(viewers).forEach(v=>{{v.resize();v.render()}})}}).observe(document.getElementById('moleculeWorkspace'));
document.getElementById('scoreDetails').ontoggle=()=>{{if(data.assemblies.length){{const active=document.querySelector('.result.active');renderScorePlot(data.assemblies,active?Number(active.dataset.index):0,select)}}}};
function showUnassembled(){{document.getElementById('moleculeWorkspace').classList.add('target-focus');document.getElementById('LP').textContent='No complete assembly in saved detections. Target shown without an assigned mapping.';list.textContent='No recommendation. Unsupported target atoms: '+data.uncovered_target_atoms.map(a=>'P'+a).join(', ');showModel('P',data.unassembled_target)}}
function redraw(){{if(!data.assemblies.length){{showUnassembled();return}}select(Number(document.querySelector('.result.active').dataset.index))}}document.getElementById('labels').onchange=redraw;document.getElementById('fragments').onchange=redraw;document.getElementById('symmetry').onchange=redraw;if(data.assemblies.length){{const knownIndex=data.assemblies.findIndex(x=>x.known);select(knownIndex>=0?knownIndex:0)}}else{{showUnassembled()}}window.onresize=()=>Object.values(viewers).forEach(v=>v.resize());</script></body></html>"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top", type=int, default=24)
    parser.add_argument("--title", default="Geometric building-block assemblies")
    parser.add_argument("--ground-truth-status", default="not evaluated")
    parser.add_argument("--ground-truth-note", default="No ground-truth metadata supplied.")
    args = parser.parse_args()
    report = json.loads(Path(args.results).read_text())
    payload = _payload(
        report, args.top, args.title, args.ground_truth_status,
        args.ground_truth_note)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_html(payload))
    print(output.resolve())


if __name__ == "__main__":
    main()

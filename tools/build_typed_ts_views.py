#!/usr/bin/env python3
"""Build self-contained R/P plus ranked initial-guess viewers."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from rxn_core.chemistry_computations import parse_xyz, write_xyz_str


def _source_xyz(directory):
    files = sorted(path for path in Path(directory).glob("*.xyz")
                   if not path.name.startswith("xtb"))
    if len(files) != 1:
        raise RuntimeError(
            f"expected one source XYZ in {directory}, found {len(files)}")
    return parse_xyz(files[0])


def _mechanisms(document):
    rp = document["rp"]
    mechanisms = []
    for index, raw in enumerate(rp.get("mechanisms") or (), 1):
        if "mapping" in raw:
            mapping = raw["mapping"]
            broken = raw.get("broken_bonds") or ()
            formed = raw.get("formed_bonds") or ()
            core = raw.get("core_atoms") or ()
        else:
            raise ValueError(f"mechanism {index} has no resolved mapping")
        mechanisms.append({
            "id": index,
            "mapping": {int(a): int(b) for a, b in mapping.items()},
            "broken": [list(map(int, bond)) for bond in broken],
            "formed": [list(map(int, bond)) for bond in formed],
            "core": list(map(int, core)),
        })
    return mechanisms


def _html(case, elements_r, xyz_r, elements_p, xyz_p, document):
    mechanisms = _mechanisms(document)
    targets = {int(item["iteration"]): item
               for item in document.get("targets") or ()}
    for mechanism in mechanisms:
        product = [xyz_p[mechanism["mapping"][atom]].tolist()
                   for atom in range(len(elements_r))]
        mechanism["reactant_xyz"] = write_xyz_str(
            elements_r, xyz_r, f"{case} reactant")
        mechanism["product_xyz"] = write_xyz_str(
            elements_r, product, f"{case} product aligned")
        ranking = next((item["ranked_initial_guesses"]
                        for item in document.get("rankings") or ()
                        if int(item["mechanism_id"]) == mechanism["id"]), [])
        rows = []
        for rank, score in enumerate(ranking, 1):
            iteration = int(score["iteration"])
            target = targets[iteration]
            target_elements, target_xyz = _source_xyz(target["hessian_cache"])
            selected = target["mechanisms"][mechanism["id"] - 1]["selected"]
            rows.append({
                "rank": rank,
                **score,
                "xyz": write_xyz_str(
                    target_elements, target_xyz,
                    f"{case} initial guess {iteration}"),
                "core_target": ([] if selected is None else
                                list(map(int, selected["assignment"].values()))),
                "candidate_count": target["mechanisms"][
                    mechanism["id"] - 1]["candidate_count"],
            })
        mechanism["ranking"] = rows
    payload = json.dumps({"case": case, "mechanisms": mechanisms}).replace(
        "</", "<\\/")
    library = (Path(__file__).parents[1] / "src" / "rxn_core" / "static"
               / "3Dmol-min.js").read_text()
    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{case} — TS ranking</title><style>
*{{box-sizing:border-box}}html,body{{height:100%;margin:0;font:13px system-ui;color:#172033;background:#eef2f7}}
body{{display:grid;grid-template-rows:auto minmax(320px,52vh) 1fr;overflow:hidden}}
header{{background:#172033;color:#fff;padding:9px 14px;display:flex;gap:15px;align-items:center}}
h1{{font-size:16px;margin:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}select{{padding:5px}}
#events{{color:#dbe7ff;margin-left:auto}}#views{{display:grid;grid-template-columns:1fr 1fr 1fr;min-height:0}}
.panel{{position:relative;border-right:1px solid #cbd5e1;background:#111}}.panel:last-child{{border:0}}
.label{{position:absolute;z-index:2;color:#fff;background:#111b;padding:6px 9px;border-radius:0 0 5px 0}}
.viewer{{position:absolute;inset:0}}#lower{{min-height:0;overflow:auto;background:white;border-top:1px solid #cbd5e1}}
table{{border-collapse:collapse;width:100%}}th{{position:sticky;top:0;background:#e8eef7;z-index:2}}
th,td{{padding:6px 8px;border-bottom:1px solid #e2e8f0;text-align:right;white-space:nowrap}}
th:nth-child(2),td:nth-child(2){{text-align:left}}tbody tr{{cursor:pointer}}tbody tr:hover{{background:#edf5ff}}
tbody tr.active{{background:#cfe5ff}}.score{{font-weight:700}}.missing{{color:#9a3412}}
@media(max-width:850px){{body{{grid-template-rows:auto minmax(500px,65vh) 1fr}}#views{{grid-template-columns:1fr;grid-template-rows:repeat(3,1fr)}}}}
</style><script>{library}</script></head><body>
<header><h1>{case}</h1><label>Mechanism <select id="mechanism"></select></label><span id="events"></span></header>
<section id="views"><div class="panel"><b class="label">Reactant</b><div id="r" class="viewer"></div></div>
<div class="panel"><b class="label">Product aligned</b><div id="p" class="viewer"></div></div>
<div class="panel"><b class="label" id="iglabel">Initial guess</b><div id="ig" class="viewer"></div></div></section>
<section id="lower"><table><thead><tr><th>Rank</th><th>Guess</th><th>Score S</th><th>Mode overlap β</th><th>WBO progress</th><th>Imag. freq.</th><th>Endpoint support</th><th>Core candidates</th></tr></thead><tbody id="rows"></tbody></table></section>
<script>const DATA={payload};const viewers={{r:$3Dmol.createViewer('r',{{backgroundColor:'#111'}}),p:$3Dmol.createViewer('p',{{backgroundColor:'#111'}}),ig:$3Dmol.createViewer('ig',{{backgroundColor:'#111'}})}};
let mech,activeIteration=null;const select=document.getElementById('mechanism'),rows=document.getElementById('rows');
DATA.mechanisms.forEach(m=>{{const o=document.createElement('option');o.value=m.id;o.textContent=m.id;select.appendChild(o)}});
function draw(viewer,xyz,core){{viewer.removeAllModels();viewer.addModel(xyz,'xyz');viewer.setStyle({{}},{{stick:{{radius:.14}},sphere:{{scale:.25}}}});if(core&&core.length)viewer.addStyle({{index:core}},{{stick:{{color:'#ff8c00',radius:.22}},sphere:{{color:'#ff8c00',scale:.42}}}});viewer.zoomTo();viewer.render()}}
function showGuess(iteration){{const row=mech.ranking.find(x=>x.iteration===iteration);if(!row)return;activeIteration=iteration;draw(viewers.ig,row.xyz,row.core_target);document.getElementById('iglabel').textContent=`Initial guess ${{iteration}} · rank ${{row.rank}}`;[...rows.children].forEach(tr=>tr.classList.toggle('active',Number(tr.dataset.iteration)===iteration))}}
function update(){{mech=DATA.mechanisms.find(x=>x.id===Number(select.value))||DATA.mechanisms[0];draw(viewers.r,mech.reactant_xyz,mech.core);draw(viewers.p,mech.product_xyz,mech.core);document.getElementById('events').textContent=`broken ${{JSON.stringify(mech.broken)}} · formed ${{JSON.stringify(mech.formed)}}`;rows.innerHTML='';for(const row of mech.ranking){{const tr=document.createElement('tr');tr.dataset.iteration=row.iteration;tr.innerHTML=`<td>${{row.rank}}</td><td>iter${{row.iteration}}</td><td class="score">${{row.score.toFixed(6)}}</td><td>${{row.overlap.toFixed(6)}}</td><td>${{row.wbo_progress.toFixed(6)}}</td><td>${{row.frequency.toFixed(2)}}</td><td>${{row.sources.join(' + ')}}</td><td>${{row.candidate_count}}</td>`;tr.onclick=()=>showGuess(row.iteration);rows.appendChild(tr)}}if(mech.ranking.length)showGuess(mech.ranking[0].iteration);else{{viewers.ig.removeAllModels();viewers.ig.render();document.getElementById('iglabel').textContent='No scorable initial guess'}}}}
select.onchange=update;update();window.onresize=()=>Object.values(viewers).forEach(v=>v.resize());</script></body></html>"""


def build_case(case_directory, work_root):
    case_directory = Path(case_directory)
    document = json.loads((case_directory / "ts_scores.json").read_text())
    case = document["case"]
    elements_r, xyz_r = _source_xyz(Path(work_root) / case / "endpoints" / "R")
    elements_p, xyz_p = _source_xyz(Path(work_root) / case / "endpoints" / "P")
    output = case_directory / "view.html"
    output.write_text(_html(
        case, elements_r, xyz_r, elements_p, xyz_p, document))
    return output


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--case")
    args = parser.parse_args(argv)
    root = Path(args.results_root)
    directories = ([root / "cases" / args.case] if args.case else
                   sorted(path.parent for path in
                          (root / "cases").glob("*/ts_scores.json")))
    outputs = [str(build_case(directory, args.work_root))
               for directory in directories]
    print(json.dumps({"viewer_count": len(outputs), "outputs": outputs}))


if __name__ == "__main__":
    main()

import importlib.util
from pathlib import Path
import sys


TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))
SPEC = importlib.util.spec_from_file_location(
    "build_retro_db_viewer", TOOLS / "build_retro_db_viewer.py")
VIEWER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VIEWER)


def test_repeated_precursor_uses_union_of_compressed_symmetry_domains():
    common = {
        "precursor_id": "ligand",
        "smiles": "C",
        "retained_atoms": [0],
        "symmetry_retained_atoms": [0, 1, 2],
        "boundary_bonds": [],
        "leftover_fragments": [],
        "complete": True,
    }
    copies = [
        dict(common, covered_target_atoms=[1], mapping=[[0, 1]],
             target_domains=[[0, [1, 2, 3]]]),
        dict(common, covered_target_atoms=[2], mapping=[[0, 2]],
             target_domains=[[0, [1, 2, 3]]]),
    ]

    group, = VIEWER._group_precursors(copies)

    assert group["multiplicity"] == 2
    assert group["covered_target_atoms"] == [1, 2]
    assert group["symmetry_target_atoms"] == [1, 2, 3]
    assert group["symmetry_retained_atoms"] == [0, 1, 2]


def test_shared_target_support_is_not_assigned_to_first_precursor(monkeypatch):
    monkeypatch.setattr(VIEWER, "mol_3d", lambda *args, **kwargs:
        ("mock sdf", [[0., 0., 0.]] * 3, ["C"] * 3))
    def source(name, atoms):
        return {"precursor_id": name, "smiles": "CCC", "retained_atoms": [0, 1],
                "covered_target_atoms": atoms, "mapping": list(enumerate(atoms)),
                "symmetry_retained_atoms": [0, 1], "target_domains": list(enumerate([[a] for a in atoms])),
                "boundary_bonds": [], "leftover_fragments": [], "complete": True}
    report = {"target_smiles": "CCC", "assemblies": [{
        "precursors": [source("A", [0, 1]), source("B", [1, 2])],
        "formed_bonds": [], "score": {}}], "construction_patterns": [],
        "scan_counts": dict(rows=2, searched=2, matched_precursors=2, fragment_candidates=2, capped=0),
        "recommendation_search_truncated": False,
        "search_scope": "Known-ingredient check, not a blind bank scan"}
    payload = VIEWER._payload(report, 20, "Overlap test")
    target = payload["assemblies"][0]["models"][-1]
    assert target['supplierAlternatives'] == [dict(atoms=[1], owners=[0,1], selected=0)]
    assert all(s['color'] != '#9ca3af' for s in target['styles'])
    assert target["labels"][1]["text"] == "P1: R1 or R2"
    assert target['labels'][1]['always']
    html = VIEWER._html(payload)
    assert "not validated reaction edits" in html
    assert "Known-ingredient check, not a blind bank scan" in html
    assert "Returned by blind recommender" not in html

    report['expected_ids'] = ['A', 'B']
    matched = VIEWER._payload(report, 20, 'Known rank')
    assert matched['summary']['known_rank'] == 1
    assert 'yes, candidate rank 1 in the searched pool' in VIEWER._html(matched)
    del report['expected_ids']

    report['validation_assemblies'] = [dict(report['assemblies'][0], construction_pattern='GT-1')]
    checked = VIEWER._payload(report, 20, 'Validation separation')
    assert len(checked['assemblies']) == 2
    assert not checked['assemblies'][0]['ground_truth']
    assert checked['assemblies'][1]['ground_truth']
    assert checked['assemblies'][1]['pattern'] == 'GT-1'
    assert checked['assemblies'][1]['complete_cover']
    assert 'Validation pattern' in VIEWER._html(checked)
    report['assemblies'][0]['pareto_layer'] = 2
    checked = VIEWER._payload(report,20,'Pareto display')
    assert checked['assemblies'][0]['pareto_layer'] == 2
    assert 'Pareto layer' in VIEWER._html(checked)
    assert checked['assemblies'][0]['score']['covered_target_atoms'] == 3
    assert checked['assemblies'][0]['score']['target_atom_count'] == 3
    for label in ('Bond breaking', 'Bond forming', 'P coverage', 'R retention', 'Pareto rank'):
        assert label in VIEWER._html(checked)

    report['assemblies'][0]['precursors'][1]['precursor_id'] = 'A'
    target = VIEWER._payload(report, 20, 'Repeated source')['assemblies'][0]['models'][-1]
    assert target['supplierAlternatives'] == []
    assert target['styles'] == [dict(indices=[0,1,2],color=VIEWER._color(0))]
    assert target['labels'][1]['text'] == 'P1: R1 (2 copies)'


def test_no_cover_shows_unassigned_target_not_a_fabricated_assembly(monkeypatch):
    monkeypatch.setattr(VIEWER, "mol_3d", lambda *args, **kwargs:
        ("mock sdf", [[0., 0., 0.]], ["C"]))
    report = {"target_smiles": "C", "assemblies": [], "construction_patterns": [],
        "scan_counts": dict(rows=1, searched=1, matched_precursors=0,
                            fragment_candidates=0, capped=0),
        "recommendation_search_truncated": False, "uncovered_target_atoms": [0]}
    payload = VIEWER._payload(report, 20, "No cover")
    assert payload["assemblies"] == []
    assert payload["unassembled_target"]["styles"] == []
    assert payload["uncovered_target_atoms"] == [0]
    assert "No complete assembly in saved detections" in VIEWER._html(payload)


def test_score_sort_controls_keep_ranks_and_validation_separate(monkeypatch):
    import shutil
    import subprocess
    import pytest
    node = shutil.which('node')
    if not node:
        pytest.skip('JavaScript test requires Node')
    monkeypatch.setattr(VIEWER, 'mol_3d', lambda *args, **kwargs: ('sdf',[[0,0,0]],['C']))
    report = dict(target_smiles='C', assemblies=[], construction_patterns=[],
        scan_counts=dict(rows=0,searched=0,matched_precursors=0,fragment_candidates=0,capped=0),
        recommendation_search_truncated=False)
    html = VIEWER._html(VIEWER._payload(report,20,'Scores'))
    functions = html[html.index('function scoreCards'):html.index('function supplierControls')]
    functions += html[html.index('function sortedAssemblies'):html.index('function renderList')]
    code = functions + """
const assert=require('node:assert/strict');
const make=(rank,retention,broken,formed,covered,validation=false)=>({pareto_layer:rank,ground_truth:validation,score:{set_atom_retention:retention,broken_bonds:broken,formed_bonds:formed,covered_target_atoms:covered,target_atom_count:10}});
const rows=[make(2,.9,8,4,10),make(1,.6,3,6,9),make(null,1,0,0,10,true)];
for(const [mode,expected] of Object.entries({pareto:[1,0,2],retention:[0,1,2],changes:[1,0,2],breaking:[1,0,2],forming:[0,1,2],coverage:[0,1,2]})){
 assert.deepEqual(sortedAssemblies(rows,mode).map(x=>x.i),expected);
}
assert.deepEqual(rows.map(a=>a.pareto_layer),[2,1,null]);
const cards=scoreCards(rows[0]);
for(const text of ['90.0%','100.0%','10 / 10','>12<','>8<','>4<'])assert.ok(cards.includes(text),text);
assert.ok(scoreCards(rows[2]).includes('validation only'));
"""
    subprocess.run([node,'-e',code],check=True,capture_output=True,text=True)


def test_score_plot_coincident_choices_and_clicks():
    import shutil
    import subprocess
    import pytest
    node = shutil.which('node')
    if not node:
        pytest.skip('JavaScript test requires Node')
    source = (TOOLS.parent/'src/rxn_core/static/retro_score_plot.js').read_text()
    code = source + """
const assert=require('node:assert/strict');
const make=(fragments,retention=.9,validation=false)=>({rank:fragments,pattern:'p',ground_truth:validation,precursors:[{}],score:{matched_fragment_count:fragments,set_atom_retention:retention,broken_bonds:8,formed_bonds:4}});
const rows=[make(7),make(6),make(6,.9,true),make(6,.6)];
const groups=scorePlotGroups(rows);
assert.equal(groups.length,2);
assert.deepEqual(groups[0].members.map(m=>m.index),[1,2,0]);
function node(){return {children:[],attrs:{},style:{},appendChild(x){this.children.push(x)},replaceChildren(){this.children=[]},setAttribute(k,v){this.attrs[k]=v}}}
const svg=node(),choices=node();
global.document={getElementById:id=>id==='scorePlot'?svg:choices,createElement:node,createElementNS:node};
let selected=-1;
renderScorePlot(rows,0,i=>selected=i);
const circles=svg.children.filter(x=>x.attrs.role==='button');
assert.equal(circles.length,2);
circles[0].onclick();assert.equal(selected,1);
assert.equal(choices.children.length,4);
choices.children[2].onclick();assert.equal(selected,2);
assert.ok(choices.children[2].textContent.includes('Validation'));
circles[1].onkeydown({key:'Enter',preventDefault(){}});assert.equal(selected,3);
assert.deepEqual(rows.map(a=>a.score.matched_fragment_count),[7,6,6,6]);
"""
    subprocess.run([node,'-e',code],check=True,capture_output=True,text=True)

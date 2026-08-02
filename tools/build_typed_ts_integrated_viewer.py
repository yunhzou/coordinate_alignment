#!/usr/bin/env python3
"""Pack typed per-case TS viewers into one self-contained batch navigator."""
from __future__ import annotations

import argparse
import base64
import gzip
import json
from pathlib import Path


MARKER = b"<!-- RXN_CORE_SHARED_3DMOL -->"


def _gzip64(value):
    return base64.b64encode(
        gzip.compress(value, compresslevel=9, mtime=0)).decode("ascii")


def _split_library(document):
    start = document.find(b"<script>")
    end = document.find(b"</script>", start)
    if start < 0 or end < 0:
        raise ValueError("case viewer lacks its inline 3Dmol script")
    end += len(b"</script>")
    return document[start:end], document[:start] + MARKER + document[end:]


def build(results_root, output):
    root = Path(results_root).resolve()
    output = Path(output).resolve()
    records, common = [], None
    for score_path in sorted((root / "cases").glob("*/ts_scores.json")):
        case_dir = score_path.parent
        document = json.loads(score_path.read_text())
        viewer = case_dir / "view.html"
        source = viewer.read_bytes()
        library, body = _split_library(source)
        if common is None:
            common = library
        elif library != common:
            raise ValueError(f"3Dmol library differs for {document['case']}")
        rankings = document.get("rankings") or ()
        records.append({
            "id": document["case"],
            "atoms": int(document["atom_count"]),
            "mechanisms": len(rankings),
            "ranked": sum(len(item["ranked_initial_guesses"])
                          for item in rankings),
            "best": max((row["score"] for item in rankings
                         for row in item["ranked_initial_guesses"]),
                        default=None),
            "seconds": float(document.get("elapsed_seconds", 0.0)),
            "payload": _gzip64(body),
        })
    if not records or common is None:
        raise ValueError("no typed case viewers found")
    payload = json.dumps(records, separators=(",", ":")).replace("</", "<\\/")
    shared = _gzip64(common)
    marker = MARKER.decode()
    html = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Typed TS analysis — {len(records)} cases</title><style>
*{{box-sizing:border-box}}html,body{{height:100%;margin:0;font:14px system-ui;color:#172033;background:#eef2f7}}
body{{display:grid;grid-template-rows:auto 1fr;overflow:hidden}}header{{background:#172033;color:#fff;padding:9px 13px;display:flex;gap:15px;align-items:center}}
h1{{font-size:17px;margin:0}}header span{{color:#cbd8ee}}main{{display:grid;grid-template-columns:340px 1fr;min-height:0}}
aside{{background:white;border-right:1px solid #cbd5e1;display:grid;grid-template-rows:auto auto 1fr;min-height:0}}
#search{{margin:10px;padding:8px;border:1px solid #aebbd0;border-radius:6px}}#current{{padding:0 10px 9px;color:#526077;font-size:12px}}
#list{{overflow:auto;border-top:1px solid #e2e8f0}}.case{{width:100%;border:0;border-bottom:1px solid #edf0f5;background:#fff;padding:8px 10px;text-align:left;cursor:pointer}}
.case:hover{{background:#edf5ff}}.case.active{{background:#cfe5ff}}.name{{font-weight:700;display:block;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}
.meta{{font-size:11px;color:#66758d}}#stage{{position:relative;min-width:0;min-height:0}}iframe{{border:0;width:100%;height:100%;background:#fff}}
#loading{{display:none;position:absolute;inset:0;place-items:center;background:#fffE;z-index:2}}#loading.show{{display:grid}}
@media(max-width:850px){{main{{grid-template-columns:230px 1fr}}header span{{display:none}}}}
</style></head><body><header><h1>Typed TS analysis</h1><span>{len(records)} cases · R/P alignment · ranked initial guesses · reactive-mode scores</span></header>
<main><aside><input id="search" type="search" placeholder="Filter cases…"><div id="current">Select a case</div><div id="list"></div></aside>
<section id="stage"><div id="loading">Loading embedded case viewer…</div><iframe id="viewer"></iframe></section></main>
<script>const CASES={payload},COMMON="{shared}",MARKER="{marker}";const list=document.getElementById('list'),frame=document.getElementById('viewer'),loading=document.getElementById('loading'),current=document.getElementById('current');let filtered=CASES.slice(),selected=null,token=0;
function bytes64(s){{const b=atob(s),a=new Uint8Array(b.length);for(let i=0;i<b.length;i++)a[i]=b.charCodeAt(i);return a}}async function unzip(s){{const stream=new Blob([bytes64(s)]).stream().pipeThrough(new DecompressionStream('gzip'));return await new Response(stream).text()}}const common=unzip(COMMON);
function render(){{list.innerHTML='';for(const c of filtered){{const b=document.createElement('button');b.className='case'+(selected===c.id?' active':'');const best=c.best===null?'n/a':c.best.toFixed(5);b.innerHTML='<span class="name"></span><span class="meta"></span>';b.querySelector('.name').textContent=c.id;b.querySelector('.meta').textContent=`${{c.atoms}} atoms · ${{c.mechanisms}} mech · ${{c.ranked}} ranked · best ${{best}}`;b.onclick=()=>openCase(c);list.appendChild(b)}}}}
async function openCase(c){{const mine=++token;selected=c.id;render();loading.classList.add('show');current.textContent=`${{c.id}} — ${{c.atoms}} atoms, ${{c.mechanisms}} mechanisms, ${{c.ranked}} ranked rows`;location.hash=encodeURIComponent(c.id);try{{const [lib,body]=await Promise.all([common,unzip(c.payload)]);if(mine!==token)return;frame.srcdoc=body.replace(MARKER,lib)}}catch(e){{frame.srcdoc='<pre style="color:#b91c1c;padding:20px">'+String(e)+'</pre>'}}finally{{if(mine===token)loading.classList.remove('show')}}}}
document.getElementById('search').oninput=e=>{{const q=e.target.value.trim().toLowerCase();filtered=CASES.filter(c=>c.id.toLowerCase().includes(q)||String(c.atoms)===q);render()}};render();const wanted=decodeURIComponent(location.hash.slice(1));openCase(CASES.find(c=>c.id===wanted)||CASES[0]);</script></body></html>"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)
    return {"output": str(output), "cases": len(records),
            "bytes": output.stat().st_size}


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    print(json.dumps(build(args.results_root, args.output), indent=2))


if __name__ == "__main__":
    main()

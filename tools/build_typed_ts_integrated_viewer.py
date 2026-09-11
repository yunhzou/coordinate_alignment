#!/usr/bin/env python3
"""Pack typed per-case TS viewers into one self-contained batch navigator."""
from __future__ import annotations

from rxn_core.viewers import style_document

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
<title>Typed TS analysis — {len(records)} cases</title></head><body><header><h1>Typed TS analysis</h1><span>{len(records)} cases · R/P alignment · ranked initial guesses · reactive-mode scores</span></header>
<main><aside><input id="search" type="search" placeholder="Filter cases…"><div id="current">Select a case</div><div id="list"></div></aside>
<section id="stage"><div id="loading">Loading embedded case viewer…</div><iframe id="viewer"></iframe></section></main>
<script>const CASES={payload},COMMON="{shared}",MARKER="{marker}";const list=document.getElementById('list'),frame=document.getElementById('viewer'),loading=document.getElementById('loading'),current=document.getElementById('current');let filtered=CASES.slice(),selected=null,token=0;
function bytes64(s){{const b=atob(s),a=new Uint8Array(b.length);for(let i=0;i<b.length;i++)a[i]=b.charCodeAt(i);return a}}async function unzip(s){{const stream=new Blob([bytes64(s)]).stream().pipeThrough(new DecompressionStream('gzip'));return await new Response(stream).text()}}const common=unzip(COMMON);
function render(){{list.innerHTML='';for(const c of filtered){{const b=document.createElement('button');b.className='case'+(selected===c.id?' active':'');const best=c.best===null?'n/a':c.best.toFixed(5);b.innerHTML='<span class="name"></span><span class="meta"></span>';b.querySelector('.name').textContent=c.id;b.querySelector('.meta').textContent=`${{c.atoms}} atoms · ${{c.mechanisms}} mech · ${{c.ranked}} ranked · best ${{best}}`;b.onclick=()=>openCase(c);list.appendChild(b)}}}}
async function openCase(c){{const mine=++token;selected=c.id;render();loading.classList.add('show');current.textContent=`${{c.id}} — ${{c.atoms}} atoms, ${{c.mechanisms}} mechanisms, ${{c.ranked}} ranked rows`;location.hash=encodeURIComponent(c.id);try{{const [lib,body]=await Promise.all([common,unzip(c.payload)]);if(mine!==token)return;frame.srcdoc=body.replace(MARKER,lib)}}catch(e){{frame.srcdoc='<pre style="color:#b91c1c;padding:20px">'+String(e)+'</pre>'}}finally{{if(mine===token)loading.classList.remove('show')}}}}
document.getElementById('search').oninput=e=>{{const q=e.target.value.trim().toLowerCase();filtered=CASES.filter(c=>c.id.toLowerCase().includes(q)||String(c.atoms)===q);render()}};render();const wanted=decodeURIComponent(location.hash.slice(1));openCase(CASES.find(c=>c.id===wanted)||CASES[0]);</script></body></html>"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(style_document(html, 'catalog', 'typed_batch'))
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

#!/usr/bin/env python3
"""Build one offline HTML navigator containing every batch case viewer."""
from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path


COMMON_MARKER = b"<!-- RXN_CORE_SHARED_3DMOL -->"


def _read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _gzip64(data: bytes) -> str:
    return base64.b64encode(
        gzip.compress(data, compresslevel=9, mtime=0)).decode("ascii")


def _split_shared_script(document: bytes):
    start = document.find(b"<script>")
    if start < 0:
        raise ValueError("case viewer has no inline script")
    end = document.find(b"</script>", start)
    if end < 0:
        raise ValueError("case viewer has an unterminated inline script")
    end += len(b"</script>")
    return document[start:end], document[:start] + COMMON_MARKER + document[end:]


def _case_records(batch_root: Path):
    records = []
    for tier in ("small", "medium", "large"):
        manifest = _read_json(batch_root / "manifests" / f"{tier}.json")
        for case in manifest["cases"]:
            records.append((int(case["source_index"]), tier, case))
    return sorted(records)


def build(batch_root: Path, output: Path):
    batch_root = batch_root.resolve()
    output = output.resolve()
    inventory = _read_json(batch_root / "inventory.json")
    batch_summary = _read_json(batch_root / "batch_summary.json")
    cases = []
    common_script = None

    for _source_index, tier, case in _case_records(batch_root):
        step = str(case["step_id"])
        summary = _read_json(batch_root / "cases" / step / "summary.json")
        rp = _read_json(batch_root / "cases" / step / "rp_stage.json")
        view_path = batch_root / "views" / step / "view.html"
        document = view_path.read_bytes()
        shared, stripped = _split_shared_script(document)
        if common_script is None:
            common_script = shared
        elif shared != common_script:
            raise ValueError(f"shared 3Dmol script differs for {step}")
        reconstructed = stripped.replace(COMMON_MARKER, common_script, 1)
        if reconstructed != document:
            raise RuntimeError(f"viewer reconstruction failed for {step}")
        mechanisms = rp.get("mechanisms") or []
        cases.append({
            "id": step,
            "tier": tier,
            "atoms": int(case["atom_count"]),
            "mechanisms": len(mechanisms),
            "seconds": round(float(summary.get("elapsed_seconds", 0.0)), 3),
            "violations": sum(int(
                (mechanism.get("index_chirality") or {}).get(
                    "selected_index_chirality_violation_count", 0)
            ) for mechanism in mechanisms),
            "sha256": hashlib.sha256(document).hexdigest(),
            "payload": _gzip64(stripped),
        })

    if common_script is None:
        raise ValueError("batch contains no case viewers")
    data = json.dumps(cases, separators=(",", ":")).replace("</", "<\\/")
    common = _gzip64(common_script)
    title = f"AAM integrated viewer — {len(cases)} cases"
    total_mechanisms = sum(case["mechanisms"] for case in cases)
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<style>
*{{box-sizing:border-box}}html,body{{height:100%;margin:0;font:14px system-ui,-apple-system,sans-serif;color:#172033;background:#eef2f7}}
body{{display:grid;grid-template-rows:auto 1fr;overflow:hidden}}header{{background:#172033;color:white;padding:10px 14px;display:flex;align-items:center;gap:18px;box-shadow:0 2px 8px #0004}}
h1{{font-size:17px;margin:0;white-space:nowrap}}.stats{{display:flex;gap:12px;flex-wrap:wrap;color:#dbe7ff}}.stats b{{color:white}}.actions{{margin-left:auto;display:flex;gap:8px}}
button,input{{font:inherit}}button{{border:1px solid #b9c5d8;border-radius:6px;background:white;padding:6px 9px;cursor:pointer}}button:hover{{background:#edf4ff}}button:disabled{{opacity:.5;cursor:default}}
main{{min-height:0;display:grid;grid-template-columns:330px 1fr}}aside{{min-height:0;background:white;border-right:1px solid #cbd5e1;display:grid;grid-template-rows:auto auto 1fr}}
#search{{margin:10px;width:calc(100% - 20px);padding:8px;border:1px solid #aebbd0;border-radius:6px}}#current{{padding:0 10px 9px;color:#526077;font-size:12px}}
#list{{overflow:auto;border-top:1px solid #e2e8f0}}.case{{width:100%;border:0;border-bottom:1px solid #edf0f5;border-radius:0;text-align:left;padding:8px 10px;display:grid;grid-template-columns:1fr auto;gap:3px;background:white}}
.case.active{{background:#dcecff}}.name{{font-weight:650;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}.meta{{font-size:11px;color:#66758d;grid-column:1/-1}}.badge{{font-size:11px;border-radius:9px;padding:1px 6px;background:#edf2f7}}
#stage{{min-width:0;min-height:0;position:relative;background:#d9e0ea}}iframe{{border:0;width:100%;height:100%;background:white}}#loading{{position:absolute;inset:0;display:none;place-items:center;background:#fffE;font-size:16px;color:#334155;z-index:2}}#loading.show{{display:grid}}
@media(max-width:850px){{main{{grid-template-columns:240px 1fr}}header{{gap:8px}}.stats{{display:none}}}}
</style></head>
<body><header><h1>{title}</h1><div class="stats"><span><b>{inventory['case_count']}</b> cases</span><span><b>{total_mechanisms}</b> mechanisms</span><span><b>{batch_summary['error_count']}</b> errors</span><span><b>{batch_summary['total_elapsed_case_seconds']:.1f}</b> case-seconds</span></div><div class="actions"><button id="prev">Previous</button><button id="next">Next</button><button id="standalone" disabled>Open standalone</button></div></header>
<main><aside><input id="search" type="search" placeholder="Filter cases…"><div id="current">Select a case</div><div id="list"></div></aside><section id="stage"><div id="loading">Loading embedded viewer…</div><iframe id="viewer" title="Selected reaction viewer"></iframe></section></main>
<script>
const CASES={data};
const COMMON_B64="{common}";
const MARKER="{COMMON_MARKER.decode()}";
const list=document.getElementById('list'),frame=document.getElementById('viewer'),loading=document.getElementById('loading'),current=document.getElementById('current');
let filtered=CASES.slice(),selected=null,currentHTML=null,loadToken=0;
function decode64(value){{const raw=atob(value),out=new Uint8Array(raw.length);for(let i=0;i<raw.length;i++)out[i]=raw.charCodeAt(i);return out}}
async function gunzip64(value){{if(!('DecompressionStream' in window))throw new Error('This browser lacks DecompressionStream support. Use a current Chrome, Firefox, Edge, or Safari.');const stream=new Blob([decode64(value)]).stream().pipeThrough(new DecompressionStream('gzip'));return await new Response(stream).text()}}
const commonPromise=gunzip64(COMMON_B64);
function renderList(){{list.innerHTML='';for(const c of filtered){{const b=document.createElement('button');b.className='case'+(selected&&selected.id===c.id?' active':'');b.innerHTML='<span class="name"></span><span class="badge"></span><span class="meta"></span>';b.querySelector('.name').textContent=c.id;b.querySelector('.badge').textContent=c.tier;b.querySelector('.meta').textContent=`${{c.atoms}} atoms · ${{c.mechanisms}} mechanism${{c.mechanisms===1?'':'s'}} · ${{c.seconds.toFixed(1)}} s`;b.onclick=()=>openCase(c);list.appendChild(b)}}}}
async function openCase(c){{const token=++loadToken;selected=c;currentHTML=null;document.getElementById('standalone').disabled=true;renderList();loading.classList.add('show');current.textContent=`${{c.id}} — ${{c.atoms}} atoms, ${{c.mechanisms}} mechanisms, violations=${{c.violations}}`;location.hash=encodeURIComponent(c.id);try{{const [common,body]=await Promise.all([commonPromise,gunzip64(c.payload)]);if(token!==loadToken)return;const html=body.replace(MARKER,common);currentHTML=html;frame.srcdoc=html;document.getElementById('standalone').disabled=false}}catch(error){{frame.srcdoc='<pre style="color:#b91c1c;padding:20px">'+String(error)+'</pre>'}}finally{{if(token===loadToken)loading.classList.remove('show')}}}}
function move(delta){{if(!selected||!filtered.length)return;const i=filtered.findIndex(c=>c.id===selected.id);openCase(filtered[(i+delta+filtered.length)%filtered.length])}}
document.getElementById('prev').onclick=()=>move(-1);document.getElementById('next').onclick=()=>move(1);
document.getElementById('standalone').onclick=()=>{{if(!currentHTML)return;const url=URL.createObjectURL(new Blob([currentHTML],{{type:'text/html'}}));window.open(url,'_blank');setTimeout(()=>URL.revokeObjectURL(url),60000)}};
document.getElementById('search').oninput=e=>{{const q=e.target.value.trim().toLowerCase();filtered=CASES.filter(c=>c.id.toLowerCase().includes(q)||c.tier.includes(q)||String(c.atoms)===q);renderList()}};
renderList();const requested=decodeURIComponent(location.hash.slice(1));openCase(CASES.find(c=>c.id===requested)||CASES[0]);
</script></body></html>"""
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    return {
        "output": str(output),
        "cases": len(cases),
        "mechanisms": total_mechanisms,
        "bytes": output.stat().st_size,
        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.batch_root, args.output), indent=2))


if __name__ == "__main__":
    main()

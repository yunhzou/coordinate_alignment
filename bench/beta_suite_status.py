"""Build a download index and concise progress snapshot for a beta case suite."""
import argparse
from html import escape
import json
from pathlib import Path
import time


def update(root):
    rows=[]
    for item in json.loads((root/'submissions.json').read_text()):
        case_path=Path(item['manifest'])
        case=json.loads(case_path.read_text())
        status_path=case_path.with_name('status.json')
        status=json.loads(status_path.read_text()) if status_path.exists() else dict(stage='pending')
        run=Path(case['run'])
        queries=[]
        for query in sorted(run.glob('query_*')):
            if not (query/'manifest.json').exists():
                continue  # A newly submitted query may still be initializing.
            manifest=json.loads((query/'manifest.json').read_text())
            parts=[json.loads(p.read_text()) for p in (query/'parts').glob('*.json')]
            queries.append(dict(query=query.name,finished=len(parts),shards=manifest['shards'],
                rows=sum(p['rows'] for p in parts),capped=sum(p['capped'] for p in parts)))
        rows.append(dict(case=case['id'],job=item['job'],stage=status['stage'],queries=queries,
            seconds=status.get('elapsed',0) if status['stage'] in ('complete','failed') else
                time.time()-status['started'] if 'started' in status else 0,
            html=f"{case['id']}/viewer.html" if case_path.with_name('viewer.html').exists() else None))
    body=[]
    for row in rows:
        link=f'<a href="{escape(row["html"])}">Open / download HTML</a>' if row['html'] else 'Pending'
        progress='; '.join(f'{q["query"]}: {q["finished"]}/{q["shards"]} shards' for q in row['queries'])
        body.append(f'<tr><td>{escape(row["case"])}</td><td>{escape(row["stage"])}</td>'
            f'<td>{row["seconds"]/60:.1f} min</td><td>{escape(progress)}</td><td>{link}</td></tr>')
    html='''<!doctype html><meta charset="utf-8"><title>Eight beta case studies</title>
<style>body{font:15px system-ui;margin:32px;color:#172033}td,th{text-align:left;padding:12px;border-bottom:1px solid #ddd}a{color:#1767be}</style>
<h1>Eight case studies · current beta workflow</h1><p>155,305 bank structures. Explicit H · tolerance 1.0 · branch cap 100 · no sweep.
Pareto retention versus structural changes; equal scores prefer fewer fragments. Each viewer includes a clickable score plot and assembled R/P mappings.</p>
<p>Independent reference validation is labelled separately from blind recommendations. Case 5 uses the corrected nitrile isomer. Case 5 matching is reused; the other seven scans are fresh.</p>
<p style="color:#b42318">Case 1: blind processing exhausted its 12 GB allocation; its HTML shows only an independently validated reference assembly, not a successful blind recommendation.</p>
<table><tr><th>Target</th><th>Status</th><th>Elapsed (incl. queue after start)</th><th>Completed scan shards</th><th>Viewer</th></tr>'''+''.join(body)+'</table>'
    (root/'index.html').write_text(html)
    (root/'progress.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows),flush=True)
    return all(r['stage'] in ('complete','failed') for r in rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--watch',action='store_true')
    args=parser.parse_args()
    while True:
        done=update(args.root)
        if done or not args.watch:
            break
        time.sleep(30)


if __name__=='__main__':
    main()

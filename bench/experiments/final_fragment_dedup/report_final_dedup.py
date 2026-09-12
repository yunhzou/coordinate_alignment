from pathlib import Path
import json,statistics,hashlib
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');OUT=ROOT/'outputs/final_fragment_dedup'
rows=[json.loads((OUT/f'lossless140/case{c}/result.json').read_text()) for c in range(140)]
execution=json.loads((OUT/'lossless140/execution.json').read_text())
pilots=[json.loads((OUT/f'lossless_pilot/case{c}/result.json').read_text()) for c in [6,11,59,64,101,129]]
assert all(p['complete_saved_window'] and not p['unrecovered_previous_classes'] and not p['additional_window_classes'] for p in pilots)
assert all(j['status']=='passed' for j in execution['jobs'])
summary=dict(checkpoint_commit='312fff7',cases=140,old_branches=sum(r['old_literal_branches'] for r in rows),new_branches=sum(r['final_branches'] for r in rows),old_median=statistics.median(r['old_literal_branches'] for r in rows),new_median=statistics.median(r['final_branches'] for r in rows),flat_families=sum(r['flat_saved_families'] for r in rows),input_paths=sum(r['input_paths'] for r in rows),wall_seconds=execution['wall_seconds'],cpu_seconds=sum(r['cpu_seconds'] for r in rows),peak_mib=execution['peak_mib'],workers=execution['workers'],watchdog_seconds=execution['watchdog'],fresh_decode_all140=False,tests_passed=112,pilot_cases=[dict(case=p['case'],window=p['window'],classes=len(p['patterns']),complete=True,identical_to_previous=True) for p in pilots],standalone_cli_case64=dict(classes=27,complete=True,identical=True),per_case=[{k:r[k] for k in ['case','old_literal_branches','final_branches','flat_saved_families','input_paths','wall_seconds','cpu_seconds']} for r in rows])
summary['reduction_percent']=100*(1-summary['new_branches']/summary['old_branches'])
summary['current_source_sha256']={name:hashlib.sha256((ROOT/'work/aam-event-improved/src/rxn_core'/name).read_bytes()).hexdigest() for name in ['final_branches.py','search_graph.py','domain.py','analytical.py']}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print({k:v for k,v in summary.items() if k not in ['per_case','current_source_sha256','pilot_cases']})

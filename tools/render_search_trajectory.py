"""Render recorded fragment-growth events in the shared 3D reaction style."""
import argparse
import json
from pathlib import Path
from rxn_core.viewers import growth_trace_html


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('trace',type=Path)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    a.output.write_text(growth_trace_html(json.loads(a.trace.read_text())))
    print(a.output.resolve())

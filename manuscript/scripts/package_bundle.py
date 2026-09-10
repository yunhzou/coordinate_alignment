"""Package the reviewable manuscript without caches or build dependencies."""
import hashlib
import json
from pathlib import Path
import zipfile

MAN=Path(__file__).resolve().parents[1]
folders=['assets','includes','figs','evidence','animations','scripts','journals','examples']
names=['manuscript.pdf','preprint.tex','references.bib','README.md','EDITORIAL_NOTES.md',
       'requirements-build.txt','Makefile','.gitignore','.gitattributes','natbib.sty','fancyhdr.sty',
       'algorithmic.sty','editing_commands.tex']
files=[MAN/name for name in names]
for name in folders:
    files.extend(p for p in (MAN/name).rglob('*') if p.is_file() and '__pycache__' not in p.parts
                 and p.name!='phd072814s.pdf')
checks={str(p.relative_to(MAN)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
with zipfile.ZipFile(MAN/'manuscript_bundle.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,'manuscript/'+str(p.relative_to(MAN)))
    z.writestr('manuscript/SHA256SUMS.json',json.dumps(checks,indent=2)+'\n')
    for name in ['artifact-validation.json','viewer-validation.json']:
        z.write(MAN/'build'/name,'manuscript/validation/'+name)
with zipfile.ZipFile(MAN/'manuscript_bundle.zip') as z:
    assert z.testzip() is None
    assert 'manuscript/manuscript.pdf' in z.namelist()
print(f'Bundle: {len(files)} files; {(MAN/"manuscript_bundle.zip").stat().st_size/1e6:.2f} MB')

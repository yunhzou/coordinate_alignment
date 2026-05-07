"""
Build ts_discovery_paper/appendix_final/ — the supplementary release:
  - benchmark/         155 anonymized step folders (Jackie_TS_<n> -> TS_<n>)
  - benchmark.zip      zipped copy of benchmark/
  - flat_view.html     anonymized per-step three-panel mode-animation viewer
  - README.md          contents + intended use

Anonymization:
  - folder + filename: Jackie_TS_<n> -> TS_<n>
  - xyz comment line 2: strip /lustre/, /users/, yunhengz, "Jackie" tokens
  - flat_view.html string substitution: Jackie_TS_ -> TS_, "Jackie" -> ""

Run from repo root.
"""
from __future__ import annotations
import json
import re
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC_BENCH = (ROOT / 'appendix_perparation' / 'Pure_Geometries_Elementary_Step'
             / 'Benchmark_Guesses_Collective_Package')
SRC_VIEWER = ROOT / 'appendix_perparation' / 'viewer' / 'flat_view.html'
OUT = ROOT / 'ts_discovery_paper' / 'appendix_final'
OUT_BENCH = OUT / 'benchmark'

# Steps we keep: 19 Jackie_TS + 136 pr*.* (= 155). Skip dataset metadata files
# (manifest.csv / summary.json / missing_files.json / backfilled_initial_guesses.json)
# because they leak absolute paths and are not needed to use the benchmark.

PATH_PAT = re.compile(r'/lustre/[^\s|]+')


def anon_step(name: str) -> str:
    """Jackie_TS_07 -> TS_07. Pass-through for pr*.* and other names."""
    if name.startswith('Jackie_TS_'):
        return 'TS_' + name[len('Jackie_TS_'):]
    return name


def anon_filename(fname: str) -> str:
    """Drop a leading Jackie_ from a filename, e.g.
    Jackie_TS_07_benchmark_plain_iter5_xxxx.xyz -> TS_07_benchmark_plain_iter5_xxxx.xyz"""
    if fname.startswith('Jackie_'):
        return fname[len('Jackie_'):]
    return fname


def anon_comment(line: str, step_anon: str) -> str:
    """Rewrite an xyz comment line (line 2) for one of these patterns:
       1. '<step> <role> from /lustre/.../...'  (reactant/product/groundtruth)
       2. 'TS guess for Jackie TS_NN: ...'      (initial_guess descriptive)
       3. 'Jackie_TS_NN ...'                    (initial_guess plain)
       4. ' energy: ... gnorm: ... xtb: ...'    (xtb output, leave as-is)
       5. '<step> from xxx.pdf | charge=... multiplicity=...' (pr*.*)
    """
    s = line.rstrip('\n')
    # Strip any /lustre/... tail entirely.
    s = PATH_PAT.sub('', s)
    # 'Jackie_TS_<n>' tokens -> 'TS_<n>'
    s = re.sub(r'Jackie_TS_(\d+)', r'TS_\1', s)
    # 'TS guess for Jackie TS_07: ...' -> 'TS guess for TS_07: ...'
    s = re.sub(r'\bJackie TS_', 'TS_', s)
    # Stray 'Jackie' (e.g. 'Jackie/SM_01/...' would already be removed by path
    # strip; this catches anything else)
    s = re.sub(r'\bJackie\b', '', s)
    # Strip 'yunhengz' if it survived (shouldn't, after path strip)
    s = s.replace('yunhengz', '')
    # Tidy up dangling ' from ' or trailing whitespace from path removal
    s = re.sub(r'\s+from\s*$', '', s)
    s = re.sub(r'\s+\|\s*$', '', s)
    s = re.sub(r'  +', ' ', s).rstrip()
    if not s:
        s = step_anon
    return s + '\n'


def copy_xyz(src: Path, dst: Path, step_anon: str) -> None:
    text = src.read_text()
    lines = text.splitlines(keepends=True)
    if len(lines) >= 2:
        lines[1] = anon_comment(lines[1], step_anon)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(''.join(lines))


def build_benchmark() -> int:
    if OUT_BENCH.exists():
        shutil.rmtree(OUT_BENCH)
    OUT_BENCH.mkdir(parents=True)
    n_steps = 0
    n_files = 0
    step_dirs = sorted(p for p in SRC_BENCH.iterdir() if p.is_dir())
    manifest_rows = ['step_id,n_initial_guess,n_reactant,n_product,n_groundtruth\n']
    for sd in step_dirs:
        step_anon = anon_step(sd.name)
        out_step = OUT_BENCH / step_anon
        counts = {'initial_guess': 0, 'reactants': 0, 'products': 0, 'groundtruth': 0}
        for sub in ('reactants', 'products', 'groundtruth', 'initial_guess'):
            sd_sub = sd / sub
            if not sd_sub.exists():
                continue
            for f in sorted(sd_sub.iterdir()):
                if f.suffix.lower() != '.xyz':
                    continue
                out_name = anon_filename(f.name)
                copy_xyz(f, out_step / sub / out_name, step_anon)
                counts[sub] += 1
                n_files += 1
        n_steps += 1
        manifest_rows.append(
            f"{step_anon},{counts['initial_guess']},{counts['reactants']},"
            f"{counts['products']},{counts['groundtruth']}\n"
        )
    (OUT_BENCH / 'manifest.csv').write_text(''.join(manifest_rows))
    print(f"  benchmark/ : {n_steps} step folders, {n_files} xyz files")
    return n_steps


def build_viewer() -> None:
    if not SRC_VIEWER.exists():
        print(f"  WARN: flat viewer not found at {SRC_VIEWER}")
        return
    text = SRC_VIEWER.read_text()
    text = re.sub(r'Jackie_TS_(\d+)', r'TS_\1', text)
    text = re.sub(r'\bJackie TS_', 'TS_', text)
    text = re.sub(r'\bJackie\b', '', text)
    text = text.replace('yunhengz', '')
    text = PATH_PAT.sub('', text)
    (OUT / 'flat_view.html').write_text(text)
    print(f"  flat_view.html : {len(text):,} chars")


README = """# El Agente Disco — supplementary release

This archive accompanies the NeurIPS 2026 submission "El Agente Disco" and
contains the full elementary-step benchmark and the human-evaluation viewer
used to score pass@1 / pass@2.

## Contents

```
appendix_final/
  benchmark.zip      155 elementary steps; reactant + product + ground-truth TS
                     + 20 LLM-generated initial guesses per step (xyz)
  flat_view.html     standalone HTML viewer; per-step three-panel
                     mode-animation display used in pass@2 expert evaluation
  README.md          this file
```

## benchmark/ layout

```
benchmark/
  manifest.csv                     step_id, file counts per role
  TS_01/
    reactants/<name>.xyz           starting material (single conformer)
    products/<name>.xyz            product(s)
    groundtruth/<name>.xyz         DFT-optimized reference TS
    initial_guess/                 20 LLM-generated TS guesses
      TS_01_benchmark_plain_iter1_<hash>.xyz
      ...
  TS_02/
  ...
  pr1.tempo_ts1/
  pr1.tempo_ts2/
  ...
```

19 steps come from a curated in-house dataset (renumbered TS_01..TS_20, with
TS_14 absent in the original numbering — kept as-is). 136 steps come from
recent open-literature reactions, named `pr<paper>.<descriptor>`.

## flat_view.html

A standalone HTML page (no server, no build step). Open it in a browser and
each step exposes three panels: ground-truth TS mode animation, top-1
verifier-selected guess, top-2 verifier-selected guess. Used to score the
pass@2 metric reported in the paper.

## Anonymization

For double-blind review, all author identifiers and absolute file paths have
been stripped from xyz comment lines, filenames, and the viewer payload. The
underlying chemistry (atom orders, coordinates, bond lists, vibrational modes)
is unchanged.

## Intended use

- Reproducing the pass@1 / pass@2 numbers reported in the paper.
- Benchmarking new TS-prediction methods against a fixed expert-curated set.
- The 20 initial guesses per step are released as-is so that downstream
  rankers / verifiers can be evaluated without re-running an expensive
  LLM-agent sampler.

## License

To be set upon de-anonymization.
"""


def write_readme() -> None:
    (OUT / 'README.md').write_text(README)
    print(f"  README.md : {len(README):,} chars")


def make_zip() -> None:
    zip_path = OUT / 'benchmark.zip'
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for f in sorted(OUT_BENCH.rglob('*')):
            if f.is_file():
                z.write(f, f.relative_to(OUT))
    print(f"  benchmark.zip : {zip_path.stat().st_size/1e6:.1f} MB")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Building {OUT.relative_to(ROOT)} ...")
    build_benchmark()
    build_viewer()
    write_readme()
    make_zip()
    print("Done.")


if __name__ == '__main__':
    main()

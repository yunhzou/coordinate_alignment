# GRAFT manuscript

First draft, 10 September 2026, using the supplied Matter Lab preprint template.
Working title: **GRAFT: Symmetry-Aware Atom Mapping through Continuous Fragment Growth**.
GRAFT is the proposed method name; benchmark files retain their original AAM identifiers.
The paper centers on continuous growth, conditional branching, and context-preserving
deduplication. The SLAP single-edge sweep is a supporting cut-policy ablation.

The author order is Yunheng Zou, Olalla Nieto Faza, Shifa Hussain,
Varinia Bernales, and Alán Aspuru-Guzik, with Bernales and Aspuru-Guzik as principal
investigators. The supplied affiliation numbering is applied to Zou (1, 6),
Bernales (2, 7), and Aspuru-Guzik (1–9). Nieto Faza is affiliated with the
Department of Organic Chemistry, Universidade de Vigo, Spain (10).
Hussain's affiliation remains to be supplied. Corresponding-author designation,
the method name, funding, and contribution statements remain for confirmation.

Start with [the compiled paper](manuscript.pdf) and the
[offline interactive animation](animations/index.html). Download the HTML file
and open it in a browser if your file viewer does not render local HTML.
It contains its own data and needs no server or network connection.

The [review bundle](manuscript_bundle.zip) contains the PDF, LaTeX sources,
template assets, figure sources/data, figures, animations, and build scripts.

## Draft and evidence

| Material | Location |
| --- | --- |
| Main entry point, title and author fields | [preprint.tex](preprint.tex) |
| Abstract | [includes/include-abstract.tex](includes/include-abstract.tex) |
| Main text | [includes/paper.tex](includes/paper.tex) |
| Supporting information | [includes/supplement.tex](includes/supplement.tex) |
| Verified initial bibliography | [references.bib](references.bib) |
| Algorithm novelty assessment and closest prior work | [evidence/novelty-audit-20260910.md](evidence/novelty-audit-20260910.md) |
| Specific transition rules and implementation comparisons | [evidence/mechanism-comparison-20260910.md](evidence/mechanism-comparison-20260910.md) |
| Items to resolve before submission | [EDITORIAL_NOTES.md](EDITORIAL_NOTES.md) |
| Data provenance and figure-to-report mapping | [evidence/README.md](evidence/README.md) |
| Source paths and SHA-256 hashes | [evidence/sources.json](evidence/sources.json) |

The template's shared body/appendix includes now dispatch to this paper.
The original supplied ZIP remains outside this directory. The class retains
its supplied layout, font, and branding, with XeTeX/Tectonic compatibility fixes.
Journal wrappers in `journals/` remain unconfigured examples; only `preprint.tex`
is the compiled paper. Their submission metadata must be replaced before use.

## Figures and animations

Four main figures and two supporting figures have vector PDF/SVG and 240-dpi PNG
versions in `figs/`.
Figure numbering below follows the compiled paper.

| Figure | Content | Source basename |
| --- | --- | --- |
| 1 | Growth, candidate deduplication, branching and state sharing | `fig1_algorithm` |
| 2 | Recorded molecular growth and live candidate counts | `fig2_growth` |
| 3 | GRAFT Golden recovery and paired recorded search CPU | `fig3_golden` |
| 4 | 140-case scores, family membership and extra alternatives | `fig4_holdout` |
| S1 | Supporting SLAP single-edge constraint ablation | `figS1_slap_ablation` |
| S2 | Concrete recovery within event windows, earlier engine | `figS2_event_windows` |

- [Movie 1, MP4](animations/movie1_continuous_growth.mp4): 19.1 seconds,
  continuous growth and two retained placements on an 18-atom example.
- [Movie 2, MP4](animations/movie2_deferred_boundaries.mp4): 58.9 seconds,
  deferred boundaries and four successive fragment decisions on a 57-atom example.
- [Movie 1, GIF](animations/movie1_continuous_growth.gif) and
  [Movie 2, GIF](animations/movie2_deferred_boundaries.gif) provide presentation previews.

The HTML viewer adds event scrubbing, up to five sampled witnesses, terminal
selection, hydrogen display, rotation, and zoom. Atom numbers are zero-based.
Assigned product labels follow source identities. These are event-enabled Python
replays from frozen source `98b01b1`; they are illustrations, not benchmark reruns
or physical reaction trajectories. The schematic join in Figure 1 is conceptual;
the two movies do not claim a measured reconvergence or isolate deduplication's
contribution to performance.

## Rebuild

On the current workspace, dependencies are isolated in `.build-deps`:

```bash
bash scripts/build.sh
make animations
make verify
PYTHONPATH=.build-deps ../../.venv/bin/python scripts/validate_viewer.py
```

Tectonic 0.17.0 compiles the paper and resolves references. It must be on `PATH`;
the first build may download TeX packages. `build/` contains its cache, logs,
render checks and review screenshots. No benchmark computation is needed.

For another machine, install Python 3.12, Tectonic, and the pinned build packages
in a dedicated environment:

```bash
python3 -m venv .build-venv
.build-venv/bin/python -m pip install -r requirements-build.txt
MANUSCRIPT_PYTHON=.build-venv/bin/python bash scripts/build.sh
make animations PYTHON=.build-venv/bin/python
make verify PYTHON=.build-venv/bin/python
```

Optional browser validation uses Playwright Chromium. Install it with
`.build-venv/bin/python -m playwright install chromium`, then run
`scripts/validate_viewer.py` with that Python. `MANUSCRIPT_BROWSER` can specify an
existing Chromium executable. This cluster's missing browser libraries were
extracted locally under `build/browser-libs`, without changing system packages.

`scripts/prepare_evidence.py` refreshes snapshots and records an illustrative
trace (case 64 by default; case 1 as an argument). It requires the original
workspace's frozen engine and scientific environment. It is intentionally
excluded from the default build: rebuilding figures must not silently replace
the recorded scientific evidence.

`scripts/validate_outputs.py` checks report hashes, benchmark counts, element
compatibility/injectivity of displayed witnesses, figure formats, movie decoding,
and PDF references/layout warnings. `scripts/validate_viewer.py` exercises the UI.
These checks do not certify exhaustive search or chemical accuracy.

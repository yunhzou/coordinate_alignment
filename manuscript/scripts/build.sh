#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p build
manuscript_python="${MANUSCRIPT_PYTHON:-../../.venv/bin/python}"
if [[ ! -x "$manuscript_python" ]]; then manuscript_python=python3; fi
PYTHONPATH=.build-deps "$manuscript_python" scripts/build_figures.py
XDG_CACHE_HOME="$(pwd)/build/cache" tectonic --keep-logs --keep-intermediates --outdir build preprint.tex
cp build/preprint.pdf manuscript.pdf

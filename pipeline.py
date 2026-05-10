"""Thin shim so `python pipeline.py <step_dir>` still works in a clone
that hasn't been pip-installed. Adds src/ to sys.path then defers to
rxn_core.pipeline:main. After `pip install -e .`, the same entry point
is available as the `rxn-core-pipeline` console script."""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from rxn_core.pipeline import main

if __name__ == "__main__":
    main()

"""Thin shim for the packaged BGCP full-view pipeline.

Adds src/ to sys.path so `python pipeline.py ...` works in a clone that
hasn't been pip-installed, then defers to `rxn_core.pipeline:main`.
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from rxn_core.pipeline import main

if __name__ == "__main__":
    main()

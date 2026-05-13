"""Compatibility shim for the BGCP full-view pipeline."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from rxn_core.pipeline import main


if __name__ == "__main__":
    main()

"""Build the optional compiled matcher kernels in place.

    .venv/bin/python bench/build_fast.py [--annotate]

Compiles ``src/rxn_core/matcher/_fast.pyx`` to ``rxn_core/matcher/_fast.*.so``
next to its source.  Generated C, object files and the staging library go to
``build/fast/`` at the repository root.  Requires ``cython`` and
``setuptools`` in the environment; the package's declared dependencies are
unchanged and the extension is used only when ``RXN_CORE_FAST=1``.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup

ROOT = Path(__file__).resolve().parents[1]
PYX = Path("src") / "rxn_core" / "matcher" / "_fast.pyx"
BUILD = Path("build") / "fast"


def main(argv):
    os.chdir(ROOT)
    (BUILD / "cython").mkdir(parents=True, exist_ok=True)
    extension = Extension(
        "rxn_core.matcher._fast",
        [str(PYX)],
        extra_compile_args=["-O3"],
    )
    setup(
        name="rxn_core_fast_kernels",
        package_dir={"": "src"},
        ext_modules=cythonize(
            [extension],
            build_dir=str(BUILD / "cython"),
            compiler_directives={"language_level": "3", "binding": True},
            annotate="--annotate" in argv,
        ),
        script_args=[
            "build_ext", "--inplace",
            "--build-temp", str(BUILD / "temp"),
            "--build-lib", str(BUILD / "lib"),
        ],
    )


if __name__ == "__main__":
    main(sys.argv[1:])

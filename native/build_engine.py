"""Build the native AAM engine extension (rxn_core._engine).

    .venv/bin/python native/build_engine.py

Compiles the vendored nauty core (native/nauty, Apache-2.0) together with
native/src/engine.cpp, freeze.cpp, autgrp.cpp and repair.cpp into
src/rxn_core/_engine.<abi>.so.  Requires pybind11
and a C++17 compiler.  The extension is optional: rxn_core falls back to the
pure-Python engine when the extension is absent or RXN_CORE_NATIVE=0.
"""
import shutil
import sys
from pathlib import Path

import pybind11
from setuptools import Distribution, Extension
from setuptools.command.build_ext import build_ext

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "native"
NAUTY = NATIVE / "nauty"
BUILD = ROOT / "build" / "engine"

nauty_sources = [str(NAUTY / name) for name in (
    "nauty.c", "nautil.c", "naugraph.c", "schreier.c", "naurng.c")]

ext = Extension(
    "rxn_core._engine",
    sources=[str(NATIVE / "src" / "engine.cpp"),
             str(NATIVE / "src" / "freeze.cpp"),
             str(NATIVE / "src" / "autgrp.cpp"),
             str(NATIVE / "src" / "repair.cpp"),
             *nauty_sources],
    include_dirs=[pybind11.get_include(), str(NAUTY)],
    language="c++",
    extra_compile_args=["-O3", "-std=c++17", "-ffp-contract=off", "-fvisibility=hidden",
                        "-DMAXN=0", "-Wno-unused-result"],
)


class Build(build_ext):
    def build_extension(self, extension):
        # nauty is C; compile it with the C compiler flags, the .cpp files as C++
        c_flags = ["-O3", "-fPIC", "-DMAXN=0", "-Wno-unused-result"]
        objects = []
        compiler = self.compiler
        for source in extension.sources:
            is_c = source.endswith(".c")
            objects += compiler.compile(
                [source], output_dir=str(BUILD / "obj"),
                include_dirs=extension.include_dirs,
                extra_postargs=c_flags if is_c else extension.extra_compile_args)
        output = self.get_ext_fullpath(extension.name)
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        compiler.link_shared_object(
            objects, output, extra_postargs=["-std=c++17"],
            target_lang="c++")


def main():
    dist = Distribution({"name": "rxn_core_engine", "ext_modules": [ext]})
    cmd = Build(dist)
    cmd.build_lib = str(BUILD / "lib")
    cmd.build_temp = str(BUILD / "tmp")
    cmd.ensure_finalized()
    cmd.run()
    built = next((BUILD / "lib" / "rxn_core").glob("_engine*.so"))
    target = ROOT / "src" / "rxn_core" / built.name
    shutil.copy2(built, target)
    print("built", target)


if __name__ == "__main__":
    main()

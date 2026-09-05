from setuptools import Extension, setup
import pybind11


setup(
    ext_modules=[
        Extension(
            "rxn_core._group_ops",
            sources=["src/rxn_core/native/group_ops.cpp"],
            include_dirs=[pybind11.get_include()],
            language="c++",
            extra_compile_args=["-O3", "-std=c++17"],
        ),
        Extension(
            "rxn_core._native",
            sources=["src/rxn_core/native/paired_mapping.cpp"],
            language="c++",
            extra_compile_args=["-O3", "-std=c++17"],
        )
    ]
)

from setuptools import Extension, setup


setup(
    ext_modules=[
        Extension(
            "rxn_core._native",
            sources=["src/rxn_core/native/paired_mapping.cpp"],
            language="c++",
            extra_compile_args=["-O3", "-std=c++17"],
        )
    ]
)

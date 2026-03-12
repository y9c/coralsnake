from setuptools import setup, Extension
import platform


import os


def get_compile_args():
    args = [
        "-O3",
        "-funroll-loops",
        "-ffast-math",
        "-DNDEBUG",
        "-fomit-frame-pointer",
        "-flto",
        "-fPIC",
    ]
    machine = platform.machine().lower()

    # Do not use -march=native when building universal wheels for PyPI
    if os.environ.get("CIBUILDWHEEL", "0") != "1":
        if machine in ["x86_64", "amd64", "arm64", "aarch64"]:
            args.append("-march=native")
    return args


setup(
    name="seqops",
    ext_modules=[
        Extension(
            "coralsnake.seqops",
            ["coralsnake/seqops.c"],
            libraries=["z"],
            extra_compile_args=get_compile_args(),
            extra_link_args=["-flto"],
        )
    ],
)

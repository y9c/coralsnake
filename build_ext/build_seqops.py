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
    ]
    machine = platform.machine().lower()
    if machine in ["x86_64", "amd64"]:
        args.append("-march=native")
        args.append("-msse4.2")
    elif machine in ["arm64", "aarch64"]:
        args.append("-march=native")
    return args


# Path to BWA static library
bwa_lib = os.path.abspath("../bwamem/bwa/libbwa.a")

setup(
    name="seqops",
    ext_modules=[
        Extension(
            "coralsnake.seqops",
            ["coralsnake/seqops.c"],
            libraries=["z"],
            extra_compile_args=get_compile_args(),
            extra_link_args=["-flto", bwa_lib],  # Link directly to BWA
        )
    ],
)

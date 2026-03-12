from setuptools import setup, Extension
import platform


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

import os
import platform
from distutils.command.build_ext import build_ext

# from setuptools import Distribution, Extension
# from setuptools.command.build_ext import build_ext
from distutils.core import Distribution, Extension

from Cython.Build import cythonize

# Define the base path for minimap2 sources
minimap2_base = "coralsnake/minimap2"

# Define extra compile arguments based on the platform
extra_compile_args = [
    "-DHAVE_KALLOC",
    "-O3",
    "-Wno-sign-compare",
    "-Wno-unused-variable",
    "-Wno-unused-but-set-variable",
    "-Wno-unused-result",
]

include_dirs = [minimap2_base]

# Adjust compile args for ARM or x86 architectures
if platform.machine() in ["aarch64", "arm64"]:
    include_dirs.append(os.path.join(minimap2_base, "sse2neon"))
    extra_compile_args.extend(["-ftree-vectorize", "-DKSW_SSE2_ONLY", "-D__SSE2__"])
else:
    extra_compile_args.append("-msse4.1")  # Note: ancient x86_64 CPUs don't have SSE4

libraries = ["z", "m", "pthread"]


def build():
    sources = [
        "python/mappy.pyx",
        "align.c",
        "bseq.c",
        "lchain.c",
        "seed.c",
        "format.c",
        "hit.c",
        "index.c",
        "pe.c",
        "options.c",
        "ksw2_extd2_sse.c",
        "ksw2_exts2_sse.c",
        "ksw2_extz2_sse.c",
        "ksw2_ll_sse.c",
        "kalloc.c",
        "kthread.c",
        "map.c",
        "misc.c",
        "sdust.c",
        "sketch.c",
        "esterr.c",
        "splitidx.c",
    ]
    depends = [
        "minimap.h",
        "bseq.h",
        "kalloc.h",
        "kdq.h",
        "khash.h",
        "kseq.h",
        "ksort.h",
        "ksw2.h",
        "kthread.h",
        "kvec.h",
        "mmpriv.h",
        "sdust.h",
        "python/cmappy.h",
        "python/cmappy.pxd",
    ]
    sources = [os.path.join(minimap2_base, file) for file in sources]
    depends = [os.path.join(minimap2_base, file) for file in depends]

    extensions = [
        Extension(
            "mappy",
            sources=sources,
            depends=depends,
            extra_compile_args=extra_compile_args,
            include_dirs=include_dirs,
            libraries=libraries,
        )
    ]
    ext_modules = cythonize(
        extensions,
        include_path=include_dirs,
        compiler_directives={"binding": True, "language_level": 3},
    )

    distribution = Distribution(
        {
            "name": "coralsnake",
            "ext_modules": ext_modules,
            "cmdclass": {"build_ext": build_ext},
        }
    )
    distribution.package_dir = {"": "coralsnake"}

    # Define the output directory for the compiled extensions
    output_dir = "coralsnake"
    cmd = build_ext(distribution)
    cmd.build_lib = output_dir  # Direct output to coralsnake directory
    cmd.inplace = 1  # Build extensions in place
    cmd.ensure_finalized()
    cmd.run()

    # Ensure the outputs are set with the correct permissions
    for output in cmd.get_outputs():
        mode = os.stat(output).st_mode
        mode |= (
            mode & 0o444
        ) >> 2  # Make read-only by owner readable by group and others
        os.chmod(output, mode)


if __name__ == "__main__":
    build()

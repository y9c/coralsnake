"""Setup script for building C extensions."""
from setuptools import Extension, setup

# Define the C extension
seqops_module = Extension(
    "coralsnake.seqops",
    sources=["coralsnake/seqops.c"],
    extra_compile_args=["-O3", "-Wall"],
)

# Build the extension
setup(
    ext_modules=[seqops_module],
)


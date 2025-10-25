#!/usr/bin/env python3
import os
import sys
from distutils.core import setup, Extension

# Define the C extension
ext_modules = [
    Extension(
        "seqops",
        ["seqops.c"],
        include_dirs=[],
        libraries=[],
        library_dirs=[],
    ),
]

setup(
    name="seqops",
    ext_modules=ext_modules,
)

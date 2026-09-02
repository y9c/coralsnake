#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# coralsnake - exon-aware RNA analysis pipeline.

# Single source of truth is pyproject.toml; never hardcode here.
try:
    from importlib.metadata import version as _pkg_version

    __version__ = _pkg_version("coralsnake")
except Exception:  # running from a source checkout that is not installed
    __version__ = "0.0.0+source"


def __getattr__(name):
    """Lazy top-level exports to keep import light.

    - ``Mlogo`` comes from coralsnake.logo; its scoring engine is pure numpy,
      so importing it does NOT require the optional matplotlib extra.
    """
    if name == "Mlogo":
        from .logo import Mlogo

        return Mlogo
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

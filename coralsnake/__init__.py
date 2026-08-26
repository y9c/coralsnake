#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# coralsnake - transcriptome mapping in two colors.

__version__ = "0.0.212"


def __getattr__(name):
    """Lazy top-level exports to keep import light.

    - ``coralsnake.metagene`` is a subpackage (always available).
    - ``Mlogo`` comes from coralsnake.logo; its scoring engine is pure numpy,
      so importing it does NOT require the optional matplotlib extra.
    """
    if name == "Mlogo":
        from .logo import Mlogo

        return Mlogo
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

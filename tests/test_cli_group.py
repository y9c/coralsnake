#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regression tests for the CLI group error handling.

Verifies that `CoralsnakeGroup` converts genuine library errors into clean
``ClickException`` messages, but lets click's own control-flow exceptions
(``Exit``/``Abort``, e.g. from ``--help`` / ``--version``) propagate so they are
NOT rendered as bogus ``Error`` panels.
"""

import contextlib
import io

import pytest

from coralsnake.cli import cli

pytest.importorskip("rich_click")


def _run(args):
    """Run the CLI and capture stdout + stderr separately."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            cli.main(args, standalone_mode=False)
        except SystemExit:
            pass
    return out.getvalue(), err.getvalue()


class TestHelpNoErrorPanel:
    @pytest.mark.parametrize(
        "cmd",
        [
            ["--help"],
            ["annot", "--help"],
            ["effect", "--help"],
            ["annotate", "--help"],
            ["liftover", "--help"],
        ],
    )
    def test_subcommand_help_has_no_error_panel(self, cmd):
        """--help must not emit a bogus 'Error' panel (regression for Exit-as-RuntimeError)."""
        out, err = _run(cmd)
        combined = out + err
        assert "Error" not in combined, f"'--help' produced an Error panel for {cmd}"


class TestExitPropagation:
    def test_version_surfaces(self):
        """--version should print the version (Exit(0) propagates, no Error box)."""
        out, err = _run(["--version"])
        assert "version" in out
        assert "Error" not in (out + err)

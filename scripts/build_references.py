#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Build (or update) the built-in reference parquet files for coralsnake.

These are the files that ``coralsnake metagene --download <ref>`` fetches from
the ``data`` GitHub release of this repository (see ``coralsnake/config.py``).
Migrated from the standalone ``y9c/metagene`` repo (``scripts/
process_gtf_to_parquet.py``); the reference list, source GTF paths and output
names now come from ``coralsnake.config.BUILTIN_REFERENCES``, so this script
stays in sync with the package automatically.

Workflow
--------
1. Place the raw GTF files under ``--source-dir`` following the layout recorded
   in ``BUILTIN_REFERENCES[ref]["source_file"]``, e.g.::

       Homo_sapiens/raw/GRCh38.release110.gtf.gz

   (Ensembl: ``https://ftp.ensembl.org/pub/release-<N>/gtf/<species>/...``
    UCSC:    ``https://hgdownload.soe.ucsc.edu/goldenPath/<asm>/...``)

2. Build::

       python scripts/build_references.py --all            # everything
       python scripts/build_references.py --single <gtf> <name>
       python scripts/build_references.py --list

3. (Re)publish to the **fixed** ``data`` release of this repo::

       python scripts/build_references.py --all --publish

   Requires the ``gh`` CLI authenticated with write access to the repo.
   By design the release is fixed: if the tag already exists the script
   refuses to overwrite it (pass ``--force`` to clobber anyway). For a data
   update use a new tag (e.g. ``data-v2``) and bump ``GITHUB_DOWNLOAD_BASE``
   in ``coralsnake/config.py``.

The output schema is exactly ``coralsnake.gtf.prepare_exon_ref`` (the same
schema custom-GTF metagene runs produce locally), so hosted and local
references are interchangeable. A lightweight invariant check (exon offsets
vs. transcript length) runs on every build as a regression guard.
"""

import argparse
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path

import polars as pl

from coralsnake.config import BUILTIN_REFERENCES
from coralsnake.gtf import prepare_exon_ref

# Expected schema of the reference parquets (see prepare_exon_ref).
EXPECTED_COLUMNS = [
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "gene_id",
    "transcript_id",
    "exon_number",
    "transcript_length",
    "transcript_level",
    "Start_exon",
    "End_exon",
    "start_codon_pos",
    "stop_codon_pos",
]


def log(msg: str) -> None:
    print(msg, flush=True)


def check_invariants(df: pl.DataFrame, name: str) -> None:
    """Guard against corrupt rebuilds (e.g. wrong exon-order cumsum).

    For every row the exon offsets must be in-range, and the last exon of a
    transcript must end exactly at its transcript length.
    """
    bad = df.filter(
        (pl.col("End_exon") > pl.col("transcript_length"))
        | (pl.col("Start_exon") < 0)
        | (pl.col("End_exon") < 0)
    ).height
    if bad:
        raise ValueError(f"{name}: {bad} rows with out-of-range exon offsets")
    per_tx = (
        df.group_by("transcript_id", "transcript_length")
        .agg(max_end_exon=pl.col("End_exon").max())
        .filter(pl.col("max_end_exon") != pl.col("transcript_length"))
    ).height
    if per_tx:
        raise ValueError(
            f"{name}: {per_tx} transcripts whose last exon does not end at the "
            "transcript length (exon order / cumsum bug)"
        )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_one(gtf_path: str, out_path: str, compression_level: int) -> Path:
    """GTF -> prepare_exon_ref -> zstd parquet (atomic via tmp rename)."""
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")
    gtf_path_obj = Path(gtf_path)
    input_mb = gtf_path_obj.stat().st_size / (1024 * 1024)

    t0 = time.time()
    log(f"Processing {gtf_path} ({input_mb:.1f} MB) ...")
    df = prepare_exon_ref(gtf_path)
    if df.height == 0:
        raise ValueError(f"No features extracted from {gtf_path}")
    t1 = time.time()

    # Schema guard: a future prepare_exon_ref change must not silently
    # alter the hosted format. Column order is normalized to the hosted
    # layout (start_codon_pos before stop_codon_pos); readers select by
    # name, so order itself is not significant.
    if set(df.columns) != set(EXPECTED_COLUMNS):
        raise ValueError(
            f"prepare_exon_ref schema changed: {sorted(df.columns)} != "
            f"{sorted(EXPECTED_COLUMNS)}"
        )
    df = df.select(EXPECTED_COLUMNS)
    check_invariants(df, gtf_path_obj.name)

    tmp = Path(out_path).with_name(Path(out_path).name + ".tmp")
    df.write_parquet(tmp, compression="zstd", compression_level=compression_level)
    os.replace(tmp, out_path)
    out_mb = os.path.getsize(out_path) / (1024 * 1024)
    log(
        f"  OK {out_path}  {df.height} rows  {out_mb:.2f} MB  "
        f"(zstd level {compression_level}, {t1 - t0:.0f}s)"
    )
    return Path(out_path)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build built-in reference parquet files for coralsnake.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--single",
        nargs=2,
        metavar=("GTF_FILE", "NAME"),
        help="process a single GTF (NAME is a known reference or an output name)",
    )
    group.add_argument(
        "--all", action="store_true", help="process all known references"
    )
    group.add_argument("--list", action="store_true", help="list all known references")

    p.add_argument(
        "--source-dir",
        default=".",
        help="prefix for the source_file paths (default: current directory)",
    )
    p.add_argument(
        "--out",
        default="data/parquet",
        help="output directory for parquet files (default: data/parquet)",
    )
    p.add_argument(
        "--compression-level",
        type=int,
        default=22,
        help="zstd compression level 1-22 (default: 22, max compression)",
    )
    p.add_argument(
        "--publish",
        action="store_true",
        help="upload the built parquets to the GitHub data release (needs gh CLI)",
    )
    p.add_argument(
        "--repo",
        default="y9c/coralsnake",
        help="target repository for --publish (default: y9c/coralsnake)",
    )
    p.add_argument(
        "--tag",
        default="data",
        help="release tag for --publish (default: data)",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="allow overwriting assets on an existing (fixed) release",
    )
    return p.parse_args(argv)


def existing_release(repo: str, tag: str) -> bool:
    r = subprocess.run(
        ["gh", "api", f"repos/{repo}/releases/tags/{tag}"],
        capture_output=True,
        text=True,
    )
    return r.returncode == 0


def publish(parquets: list[Path], repo: str, tag: str, force: bool) -> None:
    if not existing_release(repo, tag):
        log(f"Creating release {tag} on {repo} (prerelease, fixed content) ...")
        notes = (
            "Built-in metagene reference parquets for `coralsnake metagene`.\n\n"
            "This release is **fixed**: asset content is immutable. Data updates "
            f"are published under a new tag (e.g. `{tag}-v2`) with a matching "
            "bump of `GITHUB_DOWNLOAD_BASE` in `coralsnake/config.py`.\n\n"
            "SHA-256:\n" + "\n".join(f"  {sha256(p)}  {p.name}" for p in parquets)
        )
        subprocess.run(
            [
                "gh",
                "release",
                "create",
                tag,
                "--repo",
                repo,
                "--target",
                "main",
                "--prerelease",
                "--title",
                "Built-in reference data (fixed)",
                "--notes",
                notes,
            ],
            check=True,
        )
    else:
        if not force:
            sys.exit(
                f"Error: release {tag} already exists on {repo} (fixed content). "
                "Use a new tag for updated data, or pass --force to clobber."
            )
        log(f"Release {tag} exists; clobbering assets (--force) ...")

    log(f"Uploading {len(parquets)} parquet file(s) ...")
    subprocess.run(
        [
            "gh",
            "release",
            "upload",
            tag,
            "--repo",
            repo,
            "--clobber",
            *[str(p) for p in parquets],
        ],
        check=True,
    )


def main(argv=None) -> None:
    args = parse_args(argv)
    source_dir = Path(args.source_dir)

    if args.list:
        log(f"{'#':>2}  reference / source file")
        for i, (name, info) in enumerate(BUILTIN_REFERENCES.items(), 1):
            src = source_dir / info["source_file"]
            mark = "present" if src.exists() else "missing"
            log(f"{i:2}. {name:<16} {src}  [{mark}]")
        return

    jobs: list[tuple[str, str, str]] = []  # (gtf path, out name, label)
    if args.single:
        gtf, name = args.single
        if name in BUILTIN_REFERENCES:
            out_name = BUILTIN_REFERENCES[name]["parquet_file"]
        else:
            out_name = f"{name}.parquet"
        jobs.append((gtf, out_name, name))
    else:  # --all
        for name, info in BUILTIN_REFERENCES.items():
            jobs.append(
                (str(source_dir / info["source_file"]), info["parquet_file"], name)
            )

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    built, failed = [], []
    for gtf, out_name, label in jobs:
        try:
            built.append(
                build_one(gtf, str(out_dir / out_name), args.compression_level)
            )
        except Exception as e:
            log(f"  FAILED {label}: {e}")
            failed.append((label, str(e)))

    if not args.publish:
        if failed:
            sys.exit(f"{len(failed)} reference(s) failed: {[f[0] for f in failed]}")
        log(f"\nBuilt {len(built)} parquet file(s) into {out_dir}/")
        return

    # publish only if the full set is present
    if failed:
        sys.exit(f"Refusing to publish: {len(failed)} reference(s) failed")
    publish(built, args.repo, args.tag, args.force)
    log("\nDone. Users pick up the files on next `coralsnake metagene --download`.")


if __name__ == "__main__":
    main()

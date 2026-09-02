import importlib.metadata

import rich_click as click

__VERSION__ = importlib.metadata.version("coralsnake")


click.rich_click.COMMAND_GROUPS = {
    "coralsnake": [
        {
            "name": "Reference Preparation",
            "commands": ["refine", "prepare"],
        },
        {
            "name": "Reference Data",
            "commands": ["reference"],
        },
        {
            "name": "Read Mapping",
            "commands": ["liftover"],
        },
        {
            "name": "Site & Variant Annotation",
            "commands": ["annotate", "motif", "coordinate"],
        },
        {
            "name": "Genomic Analysis",
            "commands": ["metagene", "group"],
        },
        {
            "name": "Visualization",
            "commands": ["logo"],
        },
        {
            "name": "Deprecated (use `annotate`)",
            "commands": ["annot", "effect"],
        },
    ]
}

click.rich_click.STYLE_OPTION = "bold green"
# click.rich_click.STYLE_COMMAND = "bold blue"


def _parse_csv(ctx, param, value, converter, error_msg):
    """Parse a comma-separated CLI value into a list via ``converter``."""
    if not value:
        return []
    try:
        return [converter(x.strip()) for x in value.split(",")]
    except (ValueError, TypeError):
        raise click.BadParameter(error_msg)


def _metagene_parse_ints(ctx, param, value):
    """Parse a comma-separated string into a list of ints (metagene CLI)."""
    return _parse_csv(ctx, param, value, int, "Must be comma-separated integers")


def _metagene_parse_strings(ctx, param, value):
    """Parse a comma-separated string into a list of strings (metagene CLI)."""
    return _parse_csv(ctx, param, value, str, "Must be comma-separated values")


def _metagene_parse_floats(ctx, param, value):
    """Parse a comma-separated string into a list of floats (metagene CLI)."""
    return _parse_csv(ctx, param, value, float, "Must be comma-separated numbers")


def _count_input_columns(path, with_header, separator):
    """Peek at the first data line of a delimited file (gz supported).

    Returns the number of columns, or None if it cannot be determined
    (empty file or unreadable).
    """
    import xopen

    try:
        with xopen.xopen(path, mode="r") as f:
            line = f.readline()
            if with_header:
                line = f.readline()
            if not line:
                return None
            return len(line.rstrip("\r\n").split(separator))
    except Exception:
        return None


class CoralsnakeGroup(click.RichGroup):
    """Click group that surfaces deliberate data/input errors cleanly.

    Subclasses ``rich_click.RichGroup`` so the standard rich-click styled /
    grouped help keeps working, while also converting library ``ValueError``
    / ``RuntimeError`` (bad input or data) into a concise ``ClickException``
    message instead of a raw traceback. Genuine bugs (other exception types)
    still raise normally.

    Note: ``click.exceptions.Exit`` and ``Abort`` are ``RuntimeError``
    subclasses used for normal CLI control flow (e.g. ``--help`` -> ``ctx.exit()``,
    ``--version`` -> ``Exit(0)``). Those must propagate untouched so the help /
    version paths don't get rendered as bogus ``Error`` panels.
    """

    def invoke(self, ctx):
        try:
            return super().invoke(ctx)
        except (click.exceptions.Exit, click.exceptions.Abort):
            raise
        except (ValueError, RuntimeError) as e:
            raise click.ClickException(str(e)) from e


@click.group(
    cls=CoralsnakeGroup,
    invoke_without_command=False,
    help="Coralsnake (transcriptome mapping utils)",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.version_option(__VERSION__, "--version", "-v")
@click.pass_context
def cli(ctx):
    pass


# ---------------------------------------------------------------------------
# Built-in reference management: `coralsnake reference <list|download|export>`
# ---------------------------------------------------------------------------


def _resolve_reference_artifact(value, kind: str):
    """Click callback: resolve a file argument that may be a reference name.

    kind: "table" (prepare-style TSV), "gtf", "fasta" (genome FASTA), or
    "faidx" (.fai index). A known reference name resolves to the cached
    artifact (downloading/deriving on first use); anything else passes
    through unchanged (existing file paths keep working as before).
    """
    from pathlib import Path

    from .config import BUILTIN_REFERENCES

    if value is None:
        return None

    def one(v: str) -> str:
        if v in BUILTIN_REFERENCES and not Path(v).exists():
            from .download import (
                download_genome,
                ensure_reference_gtf,
                ensure_reference_table,
            )

            if kind == "table":
                return str(ensure_reference_table(v))
            if kind == "gtf":
                return str(ensure_reference_gtf(v))
            fa = download_genome(v)
            return str(fa.with_name(fa.name + ".fai") if kind == "faidx" else fa)
        if not Path(v).exists():
            raise click.ClickException(
                f"{v!r} is not an existing file and not a built-in reference "
                f"name (see `coralsnake reference list`)."
            )
        return v

    if isinstance(value, (list, tuple)):
        return tuple(one(v) for v in value)
    return one(value)


def _ref_table_callback(ctx, param, value):
    return _resolve_reference_artifact(value, "table")


def _ref_gtf_callback(ctx, param, value):
    return _resolve_reference_artifact(value, "gtf")


def _ref_fasta_callback(ctx, param, value):
    return _resolve_reference_artifact(value, "fasta")


def _ref_faidx_callback(ctx, param, value):
    return _resolve_reference_artifact(value, "faidx")


@cli.group(
    "reference",
    help="Manage the built-in exon references (list, download, export, genome).",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
def reference():
    """Manage the built-in exon references.

    The references are exon-level parquets built from canonical GTFs (the
    same schema ``coralsnake prepare`` produces locally). One download
    serves every tool: ``metagene -r`` uses the parquet directly, and
    ``liftover`` / ``annotate`` / ``motif`` also accept a reference name
    (they reuse the cached object, deriving the table/GTF on demand; see
    ``reference export`` for the explicit text views). Genome FASTAs are
    too large to ship, so they are linked: ``reference genome <ref>``
    fetches the verified upstream genome on demand.
    """


@reference.command(
    "list",
    no_args_is_help=False,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
def reference_list():
    """List the built-in references with their download sizes."""
    from rich.console import Console

    from .download import list_references

    list_references(Console())


@reference.command(
    "download",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.argument("ref", metavar="REF|GROUP|all")
@click.option(
    "--with-genome",
    is_flag=True,
    help="Also fetch the linked genome FASTA(s) (large files, ~300-950 MB "
    "compressed for human; not part of the data release).",
)
def reference_download(ref, with_genome):
    """Download reference parquet(s) into the local cache.

    REF is a reference name (e.g. GRCh38), a group (human, mouse),
    or 'all'.
    """
    from .download import download_references

    try:
        download_references(ref, with_genome=with_genome)
    except (ValueError, RuntimeError) as e:
        raise click.ClickException(str(e))


@reference.command(
    "genome",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.argument("ref", metavar="REF")
def reference_genome(ref):
    """Download the linked genome FASTA for REF (+ .fa.fai index) into the cache.

    Genome FASTAs are too large to ship in the data release; they are
    fetched on demand from the verified upstream URL recorded for each
    reference (config.GENOME_URLS), decompressed, indexed, and
    cross-checked against the reference annotation's contig names.
    """
    from .download import download_genome

    try:
        download_genome(ref)
    except (ValueError, RuntimeError) as e:
        raise click.ClickException(str(e))


@reference.command(
    "export",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.argument("ref", metavar="REF")
@click.option(
    "--table",
    "table_path",
    type=click.Path(),
    default=None,
    help=(
        "Write the `prepare` annotation table (TSV); feeds `liftover -a`, "
        "`liftover --table`, `annotate --annotation`"
    ),
)
@click.option(
    "--gtf",
    "gtf_path",
    type=click.Path(),
    default=None,
    help="Write a GTF; feeds `annotate --reference-gtf`",
)
def reference_export(ref, table_path, gtf_path):
    """Export a reference as an annotation table and/or GTF.

    The reference is downloaded first if it is not in the cache.
    """
    if not table_path and not gtf_path:
        raise click.ClickException("nothing to export: give --table and/or --gtf")
    from .io import load_reference
    from .ref_export import export_gtf, export_table

    df = load_reference(ref)
    if table_path:
        export_table(df, table_path)
    if gtf_path:
        export_gtf(df, gtf_path)


@cli.command(
    help="Refine genome fasta and gtf files.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--fasta-file", "-f", "fasta_file", help="Fasta file.", default=None)
@click.option("--gtf-file", "-g", "gtf_file", help="GTF file.", default=None)
@click.option("--outdir", "-o", "outdir", help="Output directory.", default="./")
@click.option("--name", "-n", "name", help="Name of refined genome.", default=None)
@click.option(
    "--rename-mapper",
    "-m",
    "rename_mapper",
    help="Rename mapper file (TSV format, 1st column is old seqname, 2nd column is new seqname).",
    default=None,
)
@click.option(
    "--seqname-pattern",
    "-p",
    "seqname_pattern",
    help="Seqname pattern (regex).",
    default=None,
)
@click.option(
    "--canonical-transcripts",
    "-c",
    "canonical_transcripts",
    help="Canonical transcripts file (TSV; 1st column is the transcript ID).",
    default=None,
)
def refine(
    fasta_file,
    gtf_file,
    outdir,
    name,
    rename_mapper,
    seqname_pattern,
    canonical_transcripts,
):
    """Refine genome FASTA / GTF references for downstream coralsnake commands.

    Cleans seqnames (rename map / regex filter), normalizes gene/transcript
    names and types, creates missing gene/transcript/exon rows, merges
    overlapping exons, flags the canonical transcript (is_canonical), and keeps
    codon/UTR features. faidx and GTF indexing use pysam (no external tools).
    """
    if fasta_file is None and gtf_file is None:
        raise click.ClickException(
            "Nothing to refine: pass --fasta-file and/or --gtf-file."
        )
    from .refine import refine_genome_references

    refine_genome_references(
        input_fasta=fasta_file,
        input_gtf=gtf_file,
        outdir=outdir,
        name=name,
        rename_mapper=rename_mapper,
        seqname_pattern=seqname_pattern,
        canonical_transcripts=canonical_transcripts,
    )


@cli.command(
    help="Extract primary transcript from gtf/gff file.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--gtf-file", "-g", "gtf_file", help="GTF file.", required=True)
@click.option("--fasta-file", "-f", "fasta_file", help="Fasta file.")
@click.option("--output-file", "-o", "output_file", help="Output file.", required=True)
@click.option("--seq-file", "-s", "seq_file", help="Sequence file.")
@click.option(
    "--sanitize",
    "-z",
    "sanitize",
    help="Sanitize names to remove special characters.",
    is_flag=True,
)
@click.option(
    "--with-codon",
    "-c",
    "with_codon",
    help="Include codon in the output.",
    is_flag=True,
)
@click.option(
    "--with-genename",
    "-n",
    "with_genename",
    help="Include gene name in the output.",
    is_flag=True,
)
@click.option(
    "--with-biotype",
    "-t",
    "with_biotype",
    help="Include biotype in the output.",
    is_flag=True,
)
@click.option(
    "--with-txpos",
    "-x",
    "with_txpos",
    help="Include transcript position in the output.",
    is_flag=True,
)
@click.option(
    "--filter-biotype",
    "-b",
    "filter_biotype",
    help="Filter biotype.",
    default=None,
)
@click.option(
    "--seq-upper/--seq-lower",
    "-U/-u",
    "seq_upper",
    help="Convert sequence to uppercase.",
    is_flag=True,
    default=True,
)
@click.option(
    "--line-length",
    "-l",
    "line_length",
    help="Line length. (Default: 0, no wrap)",
    default=0,
)
def prepare(
    gtf_file,
    fasta_file,
    output_file,
    seq_file,
    sanitize,
    with_codon,
    with_genename,
    with_biotype,
    with_txpos,
    filter_biotype,
    seq_upper,
    line_length,
):
    from .gtf2tx import parse_file

    if seq_file is not None and fasta_file is None:
        raise click.ClickException(
            "Error: requires --fasta-file when --seq-file is used."
        )

    parse_file(
        gtf_file,
        fasta_file,
        output_file,
        seq_file,
        sanitize,
        with_codon,
        with_genename,
        with_biotype,
        with_txpos,
        filter_biotype,
        seq_upper,
        line_length,
    )


def _run_convert(
    direction, input_bam, output_bam, annotation_file, faidx_file, threads, sort
):
    """Dispatch a BAM coordinate conversion by direction.

    ``t2g`` = transcript -> genome; ``g2t`` = genome -> transcript.
    """
    if direction == "t2g" and not faidx_file:
        raise click.ClickException(
            "--faidx-file/-f is required for --direction t2g (genome FASTA index)."
        )
    if direction == "g2t":
        from .gbam2tbam import convert_bam as g2t

        g2t(input_bam, output_bam, annotation_file, threads, sort)
    elif direction == "t2g":
        from .tbam2gbam import convert_bam as t2g

        t2g(input_bam, output_bam, annotation_file, faidx_file, threads, sort)
    else:  # pragma: no cover
        raise ValueError(f"Unknown direction: {direction!r}")


@cli.command(
    help="Convert a BAM between genome and transcript coordinates (both directions).",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--direction",
    "-d",
    "direction",
    default="t2g",
    type=click.Choice(["t2g", "g2t"]),
    help="Conversion direction: 't2g' transcript->genome (default), or 'g2t' genome->transcript.",
)
@click.option("--input-bam", "-i", "input_bam", help="Input bam file.", required=True)
@click.option(
    "--output-bam", "-o", "output_bam", help="Output bam file.", required=True
)
@click.option(
    "--annotation-file",
    "-a",
    "annotation_file",
    callback=_ref_table_callback,
    help="Annotation file (prepare-style table), or a built-in reference "
    "name (e.g. GRCh38) for the cached table.",
    required=True,
)
@click.option(
    "--table",
    "table_mode",
    is_flag=True,
    help="Convert a sites TABLE between transcript(gene) and genome coordinates",
)
@click.option(
    "--gene-col",
    "gene_col",
    default="Chrom",
    help="gene column (t2g) / chrom column (g2t)",
)
@click.option("--pos-col", "pos_col", default="Pos", help="position column (1-based)")
@click.option(
    "--strand-col", "strand_col", default="Strand", help="strand column (g2t)"
)
@click.option(
    "--faidx-file",
    "-f",
    "faidx_file",
    callback=_ref_faidx_callback,
    help="Faidx file (required for 't2g'). A built-in reference name "
    "(e.g. GRCh38) resolves to the linked genome index.",
)
@click.option("--threads", "-t", "threads", help="Threads.", default=8)
@click.option("--sort", "-s", "sort", help="Sort.", is_flag=True)
def liftover(
    direction,
    input_bam,
    output_bam,
    annotation_file,
    faidx_file,
    threads,
    sort,
    table_mode,
    gene_col,
    pos_col,
    strand_col,
):
    """Convert reads (BAM) or a sites TABLE between genome and transcript
    coordinates (choose --direction).

    ``--direction t2g`` (default) remaps transcriptome-aligned reads back to
    genome coordinates; ``--direction g2t`` remaps genome-aligned reads back to
    transcript coordinates (this is the former standalone ``gbam2tbam`` command,
    now fully fused into ``liftover``).

    ``--table`` converts a tab-separated sites table instead:
      t2g: input has a gene column (--gene-col) and 1-based transcript
           position (--pos-col); appends GenomeChrom/GenomePos.
      g2t: input has chromosome (--gene-col), position (--pos-col) and strand
           (--strand-col); appends Gene/GenePos.
    """
    if table_mode:
        from .table_liftover import run_liftover_table

        run_liftover_table(
            input_bam,
            output_bam,
            annotation_file,
            direction,
            gene_col=gene_col,
            pos_col=pos_col,
            strand_col=strand_col,
        )
        return
    _run_convert(
        direction, input_bam, output_bam, annotation_file, faidx_file, threads, sort
    )


@cli.command(
    help="Deprecated: annotate sites from a table. Use 'annotate --annotation' instead.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--input-file", "-i", "input_file", help="Input file.", required=True)
@click.option("--output-file", "-o", "output_file", help="Output file.", required=True)
@click.option(
    "--annot-file", "-a", "annot_file", help="Annotation file.", required=True
)
@click.option(
    "--cols", "-c", "cols", help="Columns of Chrom,Pos,Strand", default="1,2,3"
)
@click.option("--keep-na", "-k", "keep_na", help="Keep NA.", is_flag=True)
@click.option(
    "--collapse-annot",
    "-l",
    "collapse_annot",
    help="Collapse annotation.",
    is_flag=True,
)
@click.option("--add-count", "-n", "add_count", help="Add count.", is_flag=True)
@click.option("--skip-header", "-H", "skip_header", help="Skip header.", is_flag=True)
def annot(
    input_file,
    output_file,
    annot_file,
    cols,
    keep_na,
    collapse_annot,
    add_count,
    skip_header,
):
    # --add-count is not compatible with --collapse-annot
    if add_count and collapse_annot:
        raise click.ClickException(
            "Error: --add-count is not compatible with --collapse-annot"
        )
    click.echo(
        "Warning: 'annot' is deprecated - use 'annotate --annotation' instead.",
        err=True,
    )
    from .annot import run_annot

    run_annot(
        input_file,
        output_file,
        annot_file,
        cols,
        keep_na,
        collapse_annot,
        add_count,
        skip_header,
    )


@cli.command(
    help="Group genes and build consensus sequences.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--fasta-file", "-f", "fasta_file", help="Fasta file.", required=True, multiple=True
)
@click.option(
    "--gtf-file",
    "-g",
    "gtf_file",
    help="GTF file.",
    required=True,
    multiple=True,
)
@click.option("--output-file", "-o", "output_file", help="Output file.", required=False)
@click.option(
    "--output-consensus",
    "-c",
    "output_consensus",
    help="Output artificial FASTA file containing consensus sequences.",
    required=False,
)
@click.option(
    "--gene-name-regex",
    "-r",
    "gene_name_regex",
    help="Gene name regex.",
    default=None,
    type=str,
)
@click.option(
    "--gene-biotype-list",
    "-b",
    "gene_biotype_list",
    help="Gene biotype list.",
    default=None,
)
@click.option(
    "--gene-length-limit",
    "-l",
    "gene_length_limit",
    help="Gene length limit.",
    default=300,
    type=int,
)
@click.option(
    "--cluster-threshold",
    "-s",
    "cluster_threshold",
    help="Clustering threshold (0-1).",
    default=0.1,
    type=float,
)
@click.option("--threads", "-t", "threads", help="Threads.", default=8)
def group(
    fasta_file,
    gtf_file,
    output_file,
    output_consensus,
    gene_name_regex,
    gene_biotype_list,
    gene_length_limit,
    cluster_threshold,
    threads,
):
    from .genegroup import group_genes

    group_genes(
        fasta_file,
        gtf_file,
        output_file,
        output_consensus,
        gene_name_regex,
        gene_biotype_list,
        gene_length_limit,
        cluster_threshold,
        threads,
    )


@cli.command(
    help="Run metagene profiling analysis on genomic sites.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.Path(exists=True),
    help="Input file path (BED, TSV or CSV, etc.)",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(),
    help="Output file path (TSV, CSV)",
)
@click.option(
    "--output-score",
    "-s",
    "output_score",
    type=click.Path(),
    help="Output file for binned score statistics",
)
@click.option(
    "--output-figure",
    "-p",
    "output_figure",
    type=click.Path(),
    help="Output file for metagene plot (requires 'coralsnake[plot]')",
)
@click.option(
    "--export-profile",
    "export_profile",
    type=click.Path(),
    help="Export the metagene profile matrix TSV (feature_type, feature_midpoint, "
    "count_*) - machine-readable input for downstream tools / report renderers.",
)
@click.option(
    "--reference",
    "-r",
    "reference",
    type=str,
    help="Built-in reference genome to use (e.g., GRCh38, GRCm39)",
)
@click.option(
    "--gtf",
    "-g",
    "gtf",
    type=click.Path(exists=True),
    help="GTF file path for custom reference",
)
@click.option(
    "--region",
    type=click.Choice(["all", "5utr", "cds", "3utr"]),
    default="all",
    help="Region to analyze (default: all)",
)
@click.option(
    "--bins",
    "-b",
    type=int,
    default=100,
    help="Number of bins for analysis (default: 100)",
)
@click.option(
    "--with-header",
    "-H",
    is_flag=True,
    help="Input file has header line",
)
@click.option(
    "--separator",
    "-S",
    type=str,
    default="\t",
    help="Separator for input file (default: tab)",
)
@click.option(
    "--meta-columns",
    "-m",
    "meta_columns",
    type=str,
    default="1,2,3,6",
    callback=_metagene_parse_ints,
    help="Input column indices (1-based) for genomic coordinates. The columns should contain Chromosome,Start,End,Strand or Chromosome,Site,Strand",
)
@click.option(
    "--weight-columns",
    "-w",
    "weight_columns",
    type=str,
    default="",
    callback=_metagene_parse_ints,
    help="Input column indices (1-based) for weight/score values",
)
@click.option(
    "--weight-names",
    "-n",
    "weight_names",
    type=str,
    default="",
    callback=_metagene_parse_strings,
    help="Names for weight columns",
)
@click.option(
    "--score-transform",
    "score_transform",
    type=click.Choice(["none", "log2", "log10"]),
    default="none",
    help="Transform to apply to scores (default: none)",
)
@click.option(
    "--normalize",
    is_flag=True,
    help="Normalize scores by transcript length",
)
@click.option(
    "--list",
    "list_references_flag",
    is_flag=True,
    help="List all available built-in references and exit",
)
@click.option(
    "--download",
    "download_ref",
    type=str,
    help=(
        "Download a reference (e.g., GRCh38), a group (human/mouse), "
        "or 'all' for every reference"
    ),
)
@click.option(
    "--export-table",
    "export_table",
    type=click.Path(),
    default=None,
    help=(
        "Export the reference as the `prepare` annotation table (TSV) and exit "
        "(use with -r/--reference); feeds `liftover -a`, `liftover --table` "
        "and `annotate --annotation`"
    ),
)
@click.option(
    "--export-gtf",
    "export_gtf",
    type=click.Path(),
    default=None,
    help=(
        "Export the reference as a GTF and exit (use with -r/--reference); "
        "feeds `annotate --reference-gtf`"
    ),
)
def metagene(
    input_file,
    output_file,
    output_score,
    output_figure,
    export_profile,
    reference,
    gtf,
    region,
    bins,
    with_header,
    separator,
    meta_columns,
    weight_columns,
    weight_names,
    score_transform,
    normalize,
    list_references_flag,
    download_ref,
    export_table,
    export_gtf,
):
    """Run metagene profiling analysis on genomic sites.

    Computes the distribution of genomic sites relative to gene regions
    (5'UTR, CDS, 3'UTR) and optionally produces publication-ready profiles.
    This is the full `metagene` package migrated into a coralsnake subcommand,
    built on the high-performance `polars` + `ruranges` stack.
    """
    from .annotation import map_to_transcripts, normalize_positions
    from .config import BUILTIN_REFERENCES
    from .download import download_references, list_references
    from .gtf import load_gtf
    from .io import load_reference, load_sites
    from .plotting import plot_profile

    # Handle list references option (deprecated: `coralsnake reference list`)
    if list_references_flag:
        click.echo(
            "Warning: 'metagene --list' is deprecated - "
            "use 'coralsnake reference list' instead.",
            err=True,
        )
        from rich.console import Console

        list_references(Console())
        return

    # Handle export options (deprecated: `coralsnake reference export`)
    if export_table or export_gtf:
        click.echo(
            "Warning: 'metagene --export-table/--export-gtf' is deprecated - "
            "use 'coralsnake reference export' instead.",
            err=True,
        )
        if not reference:
            raise click.ClickException(
                "--export-table/--export-gtf require -r/--reference"
            )
        if gtf:
            raise click.ClickException(
                "--export-table/--export-gtf export a built-in reference "
                "(-r/--reference), not a custom GTF"
            )
        from .ref_export import export_gtf as _export_gtf
        from .ref_export import export_table as _export_table

        ref_df = load_reference(reference)
        if export_table:
            _export_table(ref_df, export_table)
        if export_gtf:
            _export_gtf(ref_df, export_gtf)
        return

    # Handle download option (deprecated: `coralsnake reference download`)
    if download_ref:
        click.echo(
            "Warning: 'metagene --download' is deprecated - "
            "use 'coralsnake reference download' instead.",
            err=True,
        )
        try:
            click.echo(f"Downloading {download_ref}...")
            download_references(download_ref, silent=True)
            click.echo(f"✓ Downloaded {download_ref}")
        except Exception as e:
            click.echo(f"✗ Failed to download {download_ref}: {e}", err=True)
        return

    # Validate required options for analysis
    if not input_file:
        raise click.ClickException(
            "Input file is required for analysis (use -i/--input)"
        )
    if not output_file and not output_score and not export_profile:
        raise click.ClickException(
            "Output file is required for analysis (use -o/--output, "
            "-s/--output-score or --export-profile)"
        )

    if reference and gtf:
        raise click.ClickException("Cannot specify both --reference and --gtf options")
    if not reference and not gtf:
        raise click.ClickException("Must specify either --reference or --gtf option")

    # Pre-load reference data
    if reference:
        if reference not in BUILTIN_REFERENCES:
            raise click.ClickException(
                f"Unknown built-in reference: {reference}. "
                f"Available: {list(BUILTIN_REFERENCES.keys())}"
            )
        click.echo(f"Loading reference '{reference}'...")
        exon_ref = load_reference(reference)
        if exon_ref is None:
            raise click.ClickException(f"Failed to load reference '{reference}'")
        click.echo(f"✓ Reference '{reference}' ready")
    else:
        click.echo(f"Loading GTF file '{gtf}'...")
        exon_ref = load_gtf(gtf)
        click.echo("✓ GTF file loaded")

    # Convert 1-based to 0-based indices for meta_columns.
    # The default (1,2,3,6) assumes a >=6-column BED-like file; the common
    # 3-column site file (Chrom,Site,Strand) would crash on it, so auto-fall
    # back to 1,2,3 when the input has exactly 3 columns.
    meta_col_index = [col - 1 for col in meta_columns]
    if meta_columns == [1, 2, 3, 6]:
        ncols = _count_input_columns(input_file, with_header, separator)
        if ncols == 3:
            meta_col_index = [0, 1, 2]
            click.echo(
                "Note: input has 3 columns; using meta-columns 1,2,3 "
                "(Chrom,Site,Strand) instead of the default 1,2,3,6."
            )
    weight_col_index = [col - 1 for col in weight_columns]

    input_df = load_sites(
        input_file,
        with_header=with_header,
        meta_col_index=meta_col_index,
        separator=separator,
    )
    click.echo(f"Loaded {len(input_df)} input sites")

    if bins <= 0:
        raise click.ClickException("--bins must be a positive integer")
    if weight_col_index and max(weight_col_index) >= len(input_df.columns):
        raise click.ClickException(
            f"--weight-columns references column {max(weight_col_index) + 1} but the "
            f"input has only {len(input_df.columns)} column(s)"
        )

    annotated_df = map_to_transcripts(input_df, exon_ref)
    click.echo("✓ Annotated transcripts")

    import polars as pl

    # --- apply previously-ignored options: --normalize / --score-transform /
    #     --weight-names / --region ------------------------------------------
    weight_cols = [input_df.columns[i] for i in weight_col_index]
    if score_transform != "none" or normalize:
        for wc in weight_cols:
            expr = pl.col(wc).cast(pl.Float64, strict=False)
            if normalize:
                expr = expr / pl.col("transcript_length")
            if score_transform == "log2":
                expr = expr.log(2.0)
            elif score_transform == "log10":
                expr = expr.log(10.0)
            annotated_df = annotated_df.with_columns(expr.alias(wc))
    if weight_names:
        if len(weight_names) != len(weight_cols):
            raise click.ClickException(
                "--weight-names count must match --weight-columns"
            )
        for old, new in zip(weight_cols, weight_names):
            annotated_df = annotated_df.rename({old: new})
        # rename() preserves column positions, so re-locate the (renamed)
        # weight columns rather than assuming they are at 0..n-1.
        weight_col_index = [annotated_df.columns.index(n) for n in weight_names]
    if region != "all":
        target = {"5utr": "5UTR", "cds": "CDS", "3utr": "3UTR"}[region]
        annotated_df = (
            annotated_df.with_columns(
                transcript_pos=(pl.col("transcript_start") + pl.col("transcript_end"))
                // 2
            )
            .with_columns(
                feature_type=pl.when(
                    pl.col("transcript_pos") < pl.col("start_codon_pos")
                )
                .then(pl.lit("5UTR"))
                .when(pl.col("transcript_pos") > pl.col("stop_codon_pos"))
                .then(pl.lit("3UTR"))
                .otherwise(pl.lit("CDS"))
            )
            .filter(pl.col("feature_type") == target)
            .drop(["feature_type", "transcript_pos"])
        )

    if output_score or output_figure or export_profile:
        gene_bins, gene_stats, gene_splits = normalize_positions(
            annotated_df,
            split_strategy="median",
            bin_number=bins,
            weight_col_index=weight_col_index,
        )
        click.echo(
            f"Gene splits - 5'UTR: {gene_splits[0]:.3f}, "
            f"CDS: {gene_splits[1]:.3f}, 3'UTR: {gene_splits[2]:.3f}"
        )
        if gene_splits == (0.0, 0.0, 0.0):
            click.echo(
                "Warning: no mapped sites have CDS coordinates - the profile is "
                "empty. Check that the reference has start_codon/stop_codon features.",
                err=True,
            )

    # Save annotated data
    if output_file:
        annotated_df.write_csv(output_file, separator=separator)
        click.echo(f"✓ Saved annotated intervals to: {output_file}")

    # Save score statistics (if requested)
    if output_score:
        gene_bins.insert_column(
            0,
            pl.when(pl.col("feature_midpoint") < gene_splits[0])
            .then(pl.lit("5UTR"))
            .when((pl.col("feature_midpoint") > gene_splits[0] + gene_splits[1]))
            .then(pl.lit("3UTR"))
            .otherwise(pl.lit("CDS"))
            .alias("feature_type"),
        ).write_csv(output_score, separator=separator)
        click.echo(f"✓ Saved binned statistics to: {output_score}")

    # Save metagene profile matrix (machine-readable; same layout as output_score)
    if export_profile:
        gene_bins.insert_column(
            0,
            pl.when(pl.col("feature_midpoint") < gene_splits[0])
            .then(pl.lit("5UTR"))
            .when((pl.col("feature_midpoint") > gene_splits[0] + gene_splits[1]))
            .then(pl.lit("3UTR"))
            .otherwise(pl.lit("CDS"))
            .alias("feature_type"),
        ).write_csv(export_profile, separator=separator)
        click.echo(f"✓ Saved metagene profile to: {export_profile}")

    # Generate plot (optional matplotlib)
    if output_figure:
        plot_profile(gene_bins, gene_splits, output_figure)
        click.echo(f"✓ Saved plot to: {output_figure}")


@cli.command(
    help="Plot a DNA/RNA sequence-logo (requires 'coralsnake[plot]').",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--motifs",
    "-m",
    "motifs",
    type=str,
    multiple=True,
    help="Motif sequence(s). Repeatable, or comma-separated (e.g. -m ACGT -m ACGG)",
)
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.File("r"),
    help="Input file, one motif sequence per line (optionally 'seq<TAB>count' for weights). "
    "Use '-' to read from stdin (e.g. when piped from `coralsnake motif`).",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(),
    help="Output image file (e.g. logo.png, logo.svg). Optional if --matrix is given.",
)
@click.option(
    "--weights",
    "-w",
    "weights",
    type=str,
    default="",
    callback=_metagene_parse_floats,
    help="Comma-separated weights for the motifs (one per motif)",
)
@click.option("--t2u/--no-t2u", default=True, help="Convert T to U (default: T->U)")
@click.option(
    "--2bit/--no-2bit", "to2bit", default=True, help="Use 2-bit logo (default)"
)
@click.option("--normed", is_flag=True, default=False, help="Normalize letter heights")
@click.option(
    "--matrix",
    "matrix_file",
    type=click.Path(),
    help="Export the position x base score-matrix TSV (from Mlogo.scores) "
    "instead of / in addition to the figure (machine-readable output).",
)
def logo(motifs, input_file, output_file, weights, t2u, to2bit, normed, matrix_file):
    """Plot a sequence logo from a set of motif sequences.

    Uses the Mlogo engine migrated from the standalone `motiflogo` package.
    matplotlib is optional; it must be installed via ``pip install coralsnake[plot]``.
    """
    from .logo import Mlogo
    from .logo import _require_plotting

    # This raises a helpful error if matplotlib (optional 'plot' extra) is missing
    # - only needed when actually drawing the figure; the matrix export is pure numpy.
    if output_file:
        _require_plotting()

    motifs_list = []
    motif_weights = []
    if input_file is not None:
        # click.File already opened the path (or wired '-' to stdin)
        for raw in input_file:
            raw = raw.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split("\t")
            motifs_list.append(parts[0])
            # One weight per line (default 1.0) so motif_weights stays
            # aligned with motifs_list even when some lines carry no
            # explicit weight.
            if len(parts) > 1:
                try:
                    motif_weights.append(float(parts[1]))
                except ValueError:
                    motif_weights.append(1.0)
            else:
                motif_weights.append(1.0)
        if getattr(input_file, "name", "") != "-":
            input_file.close()
    for m in motifs:
        motifs_list.extend([x.strip() for x in m.split(",") if x.strip()])

    if len(motifs_list) == 0:
        raise click.ClickException("No motifs provided (use -m/--motifs or -i/--input)")

    if weights:
        if len(weights) != len(motifs_list):
            raise click.ClickException(
                f"--weights has {len(weights)} values but {len(motifs_list)} motifs given"
            )
        motif_weights = list(weights)
    elif not input_file:
        motif_weights = []

    mlogo = Mlogo(
        motifs=motifs_list,
        weights=motif_weights,
        t2u=t2u,
        to2bit=to2bit,
        normed=normed,
    )

    # machine-readable matrix export (preferred for downstream tools)
    if matrix_file:
        bases = ["A", "C", "G", "T", "U"]
        extra = []
        for col in mlogo.scores:
            for b, _ in col:
                if b not in bases and b not in extra:
                    extra.append(b)
        all_bases = bases + extra
        with open(matrix_file, "w") as fh:
            fh.write("position\t" + "\t".join(all_bases) + "\n")
            for i, col in enumerate(mlogo.scores, start=1):
                d = {b: v for b, v in col}
                fh.write(
                    f"{i}\t" + "\t".join(str(d.get(b, 0.0)) for b in all_bases) + "\n"
                )
        click.echo(f"✓ Saved logo matrix to: {matrix_file}")

    if output_file:
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(0.75 * len(mlogo.scores), 2.5))
        ax = fig.gca()
        mlogo.plot(ax=ax)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()
        click.echo(f"✓ Saved motif logo to: {output_file}")
    elif not matrix_file:
        raise click.ClickException("Provide --output and/or --matrix")


@cli.command(
    help="Fetch a genomic motif around the given sites.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--input", "-i", "input_file", default="-", help="Input position file.")
@click.option("--output", "-o", "output_file", default="-", help="Output file.")
@click.option(
    "--fasta",
    "-f",
    "fasta",
    required=True,
    callback=_ref_fasta_callback,
    help="Reference fasta file, or a built-in reference name (e.g. GRCh38): "
    "resolves to the linked genome FASTA (fetched on demand).",
)
@click.option(
    "--npad",
    "-n",
    "npad",
    default="10",
    help="Number of padding bases to call motif. Use comma for different left/right pads (eg. 2,3).",
)
@click.option("--with-header", "-H", is_flag=True, help="Input file has a header line.")
@click.option(
    "--columns",
    "-c",
    "columns",
    default="1,2,3",
    show_default=True,
    help="Sets columns for site info. (Chrom,Pos,Strand)",
)
@click.option("--to-upper", "-u", is_flag=True, help="Convert motif to upper case.")
@click.option("--wrap-site", "-w", is_flag=True, help="Wrap motif site.")
def motif(
    input_file, output_file, fasta, npad, with_header, columns, to_upper, wrap_site
):
    """Fetch a genomic motif around each variant site (strand-aware)."""
    from .motif import run_motif

    if "," in npad:
        parts = npad.split(",")
        if len(parts) != 2:
            raise click.ClickException(
                f"--npad should be 'N' or 'left,right', not {npad!r}"
            )
        lpad, rpad = parts
    else:
        lpad, rpad = npad, npad
    if not lpad.isdigit() or not rpad.isdigit():
        raise click.ClickException(f"--npad should be a positive integer, not {npad!r}")
    try:
        run_motif(
            input_file,
            output_file,
            fasta,
            int(lpad),
            int(rpad),
            with_header,
            columns,
            to_upper,
            wrap_site,
        )
    except ValueError as e:
        raise click.ClickException(str(e))


@cli.command(
    help="Map chromosome names between reference coordinate systems.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--input", "-i", "input_file", default="-", help="Input position file.")
@click.option("--output", "-o", "output_file", default="-", help="Output file.")
@click.option(
    "--reference-mapping",
    "-m",
    "reference_mapping",
    help="Mapping file (chrom in input \t chrom in reference db).",
)
@click.option(
    "--buildin-mapping",
    "-M",
    "buildin_mapping",
    help="Built-in mapping: U2E, E2U, U2E-hg38, E2U-hg38, U2E-mm39, E2U-mm39.",
)
@click.option(
    "--columns", "-c", "columns", default="1", help="Columns for site info (Chrom)."
)
@click.option("--with-header", "-H", is_flag=True, help="Input file has a header line.")
@click.option("--keep-original", "-k", is_flag=True, help="Keep original chrom name.")
def coordinate(
    input_file,
    output_file,
    reference_mapping,
    buildin_mapping,
    columns,
    with_header,
    keep_original,
):
    """Rename chromosome names between reference coordinate systems."""
    from .coordinate import run_coordinate

    try:
        run_coordinate(
            input_file,
            output_file,
            reference_mapping,
            buildin_mapping,
            columns,
            with_header,
            keep_original,
        )
    except ValueError as e:
        raise click.ClickException(str(e))


@cli.command(
    help=(
        "Annotate sites or variants (region + gene/transcript position + effect). "
        "Unifies 'annot' and 'effect'."
    ),
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--input", "-i", "input_file", default="-", help="Input position file.")
@click.option(
    "--output", "-o", "output_file", default="-", help="Output annotation file."
)
@click.option(
    "--reference-gtf",
    "-g",
    "reference_gtf",
    callback=_ref_gtf_callback,
    help="Reference GTF file, or a built-in reference name (e.g. GRCh38): "
    "uses the cached reference, deriving/caching the GTF on demand.",
)
@click.option(
    "--annotation",
    "annotation_table",
    callback=_ref_table_callback,
    help="Precomputed annotation table from `prepare` (fast table mode), "
    "or a built-in reference name for the cached table.",
)
@click.option(
    "--reference-transcript",
    "-f",
    "reference_transcript",
    multiple=True,
    callback=_ref_fasta_callback,
    help="Reference genome FASTA (needed for motif / codon / amino-acid). "
    "A built-in reference name resolves to the linked genome FASTA "
    "(fetched on demand, e.g. `annotate ... -f GRCh38`).",
)
@click.option("--strandness", "-s", is_flag=True, help="Use strand information.")
@click.option(
    "--npad", "-n", "npad", default=10, type=int, help="Padding bases for motif."
)
@click.option(
    "--all-effects",
    "-a",
    "all_effects",
    is_flag=True,
    help="Output all overlapping effects.",
)
@click.option("--with-header", "-H", is_flag=True, help="Input file has a header line.")
@click.option(
    "--columns",
    "-c",
    "columns",
    default="1,2,3,4,5",
    help="Columns for site info. (Chrom,Pos,Strand,Ref,Alt) Ref/Alt optional. Pos is 1-based.",
)
def annotate(
    input_file,
    output_file,
    reference_gtf,
    annotation_table,
    reference_transcript,
    npad,
    strandness,
    all_effects,
    with_header,
    columns,
):
    """Annotate genomic sites or variants (site + variant effect, unified)."""
    from .annotate import run_annotate

    try:
        run_annotate(
            input_file,
            output_file,
            reference_gtf,
            reference_transcript,
            npad,
            strandness,
            all_effects,
            with_header,
            columns,
            annotation_table,
        )
    except ValueError as e:
        raise click.ClickException(str(e))


@cli.command(
    help="Deprecated: annotate variant effects. Use 'annotate' (with a GTF + FASTA).",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--input", "-i", "input_file", default="-", help="Input position file.")
@click.option(
    "--output", "-o", "output_file", default="-", help="Output annotation file."
)
@click.option(
    "--reference-gtf",
    "-g",
    "reference_gtf",
    type=click.Path(exists=True),
    help="Reference GTF file.",
)
@click.option(
    "--reference-transcript",
    "reference_transcript",
    multiple=True,
    help="Reference transcript FASTA file(s).",
)
@click.option(
    "--reference-protein",
    "reference_protein",
    multiple=True,
    help="Reference protein FASTA file(s).",
)
@click.option("--strandness", "-s", is_flag=True, help="Use strand information.")
@click.option("--pU-mode", "-u", "pU_mode", is_flag=True, help="Prioritise RNA genes.")
@click.option(
    "--npad", "-n", "npad", default=10, type=int, help="Padding bases for motif."
)
@click.option("--all-effects", "-a", is_flag=True, help="Output all effects.")
@click.option("--with-header", "-H", is_flag=True, help="Input file has a header line.")
@click.option(
    "--columns",
    "-c",
    "columns",
    default="1,2,3,4,5",
    help="Columns for site info. (Chrom,Pos,Strand,Ref,Alt) Pos is 1-based.",
)
def effect(
    input_file,
    output_file,
    reference_gtf,
    reference_transcript,
    reference_protein,
    npad,
    strandness,
    all_effects,
    pU_mode,
    with_header,
    columns,
):
    """Annotate genomic variant effects."""
    click.echo(
        "Warning: 'effect' is deprecated - use 'annotate' (with a GTF + FASTA) instead.",
        err=True,
    )
    from .effect import run_effect

    try:
        run_effect(
            input_file,
            output_file,
            reference_gtf,
            reference_transcript,
            reference_protein,
            npad,
            strandness,
            all_effects,
            pU_mode,
            with_header,
            columns,
        )
    except ValueError as e:
        raise click.ClickException(str(e))


if __name__ == "__main__":
    cli()

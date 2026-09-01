import importlib.metadata

import rich_click as click

__VERSION__ = importlib.metadata.version("coralsnake")


click.rich_click.COMMAND_GROUPS = {
    "coralsnake": [
        {
            "name": "Read Mapping",
            "commands": ["prepare", "liftover"],
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


def _run_convert(direction, input_bam, output_bam, annotation_file, faidx_file, threads, sort):
    """Dispatch a BAM coordinate conversion by direction.

    ``t2g`` = transcript -> genome; ``g2t`` = genome -> transcript.
    """
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
    "--annotation-file", "-a", "annotation_file", help="Annotation file.", required=True
)
@click.option(
    "--faidx-file", "-f", "faidx_file", help="Faidx file (required for 't2g')."
)
@click.option("--threads", "-t", "threads", help="Threads.", default=8)
@click.option("--sort", "-s", "sort", help="Sort.", is_flag=True)
def liftover(direction, input_bam, output_bam, annotation_file, faidx_file, threads, sort):
    """Convert a BAM between genome and transcript coordinates (choose --direction).

    ``--direction t2g`` (default) remaps transcriptome-aligned reads back to
    genome coordinates; ``--direction g2t`` remaps genome-aligned reads back to
    transcript coordinates (this is the former standalone ``gbam2tbam`` command,
    now fully fused into ``liftover``).
    """
    _run_convert(direction, input_bam, output_bam, annotation_file, faidx_file, threads, sort)


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
    help="GTF/GFF file path for custom reference",
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
    help="Download a specific reference (e.g., GRCh38) or 'all' for all references",
)
def metagene(
    input_file,
    output_file,
    output_score,
    output_figure,
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

    # Handle list references option
    if list_references_flag:
        from rich.console import Console

        list_references(Console())
        return

    # Handle download option
    if download_ref:
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
    if not output_file:
        raise click.ClickException(
            "Output file is required for analysis (use -o/--output)"
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
            raise click.ClickException("--weight-names count must match --weight-columns")
        for old, new in zip(weight_cols, weight_names):
            annotated_df = annotated_df.rename({old: new})
        # rename() preserves column positions, so re-locate the (renamed)
        # weight columns rather than assuming they are at 0..n-1.
        weight_col_index = [annotated_df.columns.index(n) for n in weight_names]
    if region != "all":
        target = {"5utr": "5UTR", "cds": "CDS", "3utr": "3UTR"}[region]
        annotated_df = (
            annotated_df.with_columns(
                transcript_pos=(pl.col("transcript_start") + pl.col("transcript_end")) // 2
            )
            .with_columns(
                feature_type=pl.when(pl.col("transcript_pos") < pl.col("start_codon_pos"))
                .then(pl.lit("5UTR"))
                .when(pl.col("transcript_pos") > pl.col("stop_codon_pos"))
                .then(pl.lit("3UTR"))
                .otherwise(pl.lit("CDS"))
            )
            .filter(pl.col("feature_type") == target)
            .drop(["feature_type", "transcript_pos"])
        )

    if output_score or output_figure:
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
    type=click.Path(exists=True),
    help="Input file, one motif sequence per line (optionally 'seq<TAB>count' for weights)",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(),
    required=True,
    help="Output image file (e.g. logo.png, logo.svg)",
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
def logo(motifs, input_file, output_file, weights, t2u, to2bit, normed):
    """Plot a sequence logo from a set of motif sequences.

    Uses the Mlogo engine migrated from the standalone `motiflogo` package.
    matplotlib is optional; it must be installed via ``pip install coralsnake[plot]``.
    """
    from .logo import Mlogo
    from .logo import _require_plotting

    # This raises a helpful error if matplotlib (optional 'plot' extra) is missing
    _require_plotting()
    import matplotlib.pyplot as plt

    motifs_list = []
    motif_weights = []
    if input_file:
        with open(input_file) as f:
            for raw in f:
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

    fig = plt.figure(figsize=(0.75 * len(mlogo.scores), 2.5))
    ax = fig.gca()
    mlogo.plot(ax=ax)
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    click.echo(f"✓ Saved motif logo to: {output_file}")


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
    type=click.Path(exists=True),
    required=True,
    help="Reference fasta file.",
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
    type=click.Path(exists=True),
    help="Reference GTF file (GTF mode).",
)
@click.option(
    "--annotation",
    "annotation_table",
    type=click.Path(exists=True),
    help="Precomputed annotation table from `prepare` (fast table mode).",
)
@click.option(
    "--reference-transcript",
    "-f",
    "reference_transcript",
    multiple=True,
    help="Reference genome FASTA (needed for motif / codon / amino-acid).",
)
@click.option("--strandness", "-s", is_flag=True, help="Use strand information.")
@click.option(
    "--npad", "-n", "npad", default=10, type=int, help="Padding bases for motif."
)
@click.option("--all-effects", "-a", "all_effects", is_flag=True, help="Output all overlapping effects.")
@click.option("--with-header", "-H", is_flag=True, help="Input file has a header line.")
@click.option(
    "--columns",
    "-c",
    "columns",
    default="1,2,3,4,5",
    help="Columns for site info. (Chrom,Pos,Strand,Ref,Alt) Ref/Alt optional.",
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
    help="Columns for site info. (Chrom,Pos,Strand,Ref,Alt)",
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

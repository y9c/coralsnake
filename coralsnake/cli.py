import importlib.metadata

import rich_click as click

__VERSION__ = importlib.metadata.version("coralsnake")


class OptionEatAll(click.Option):
    """Custom Click option that consumes all arguments until the next flag.

    Based on: https://stackoverflow.com/questions/47631914/how-to-pass-several-list-of-arguments-to-click-option
    """

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop("save_other_options", True)
        nargs = kwargs.pop("nargs", -1)
        assert nargs == -1, f"nargs, if set, must be -1 not {nargs}"
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None
        # Ensure multiple is set to handle multiple values properly
        self.multiple = True

    def add_to_parser(self, parser, ctx):
        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            values = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                            break
                    if not done:
                        values.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                values += state.rargs
                state.rargs[:] = []

            # Process each value individually through Click's mechanism
            for v in values:
                self._previous_parser_process(v, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


click.rich_click.COMMAND_GROUPS = {
    "coralsnake": [
        {
            "name": "Read Mapping",
            "commands": ["prepare", "map", "liftover"],
        },
        {
            "name": "Site & Variant Annotation",
            "commands": ["annotate", "annot", "effect", "motif", "coordinate"],
        },
        {
            "name": "Genomic Analysis",
            "commands": ["metagene", "group"],
        },
        {
            "name": "Visualization",
            "commands": ["logo"],
        },
    ]
}

click.rich_click.OPTION_GROUPS = {
    "coralsnake map": [
        {
            "name": "Input/Output",
            "options": ["-1", "-2", "-r", "-o", "-u", "--report"],
        },
        {
            "name": "Strand-Specific Mapping",
            "options": [
                "--fwd-lib",
                "--rev-lib",
                "--fwd-ref",
                "--rev-ref",
                "--dbl-ref",
            ],
        },
        {
            "name": "Mapping Parameters",
            "options": [
                "-m",
                "--min-alignment-length",
                "--min-mapping-ratio",
                "--max-a2g-ratio",
                "--max-c2t-ratio",
                "--batch-size",
                "-t",
            ],
        },
        {
            "name": "Index Options",
            "options": ["--index-dir", "--index-only"],
        },
    ],
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
    help="Sanitize name to remove specical charaters.",
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


@cli.command(
    help="Remap transcriptome-aligned reads back to genome coordinates.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--input-bam", "-i", "input_bam", help="Input bam file.", required=True)
@click.option(
    "--output-bam", "-o", "output_bam", help="Output bam file.", required=True
)
@click.option(
    "--annotation-file", "-a", "annotation_file", help="Annotation file.", required=True
)
@click.option("--faidx-file", "-f", "faidx_file", help="Faidx file.", required=True)
@click.option("--threads", "-t", "threads", help="Threads.", default=8)
@click.option("--sort", "-s", "sort", help="Sort.", is_flag=True)
def liftover(input_bam, output_bam, annotation_file, faidx_file, threads, sort):
    from .tbam2gbam import convert_bam

    convert_bam(input_bam, output_bam, annotation_file, faidx_file, threads, sort)


@cli.command(
    help="""Map reads to reference genome using BWA-MEM.
    
    BWA-MEM parameters used by default:
    - Softclip supplementary alignments (softclip_supplementary=True)
    - Mark secondary alignments (mark_secondary=True)
    - Clipping penalties: -L 6,6 (clip_penalties=(6,6))
    - Unpaired read penalty: -U 24 (unpaired_penalty=24)
    - Minimum score threshold: -T 20 (min_score=20)
    - Insert size model: -I 80,60,450 (mean=80, std=60, max=450)
    
    These parameters optimize mapping for dual-base conversion chemistry (MK/KM).
    """,
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("-1", "--r1-file", help="r1 file", required=False)
@click.option("-2", "--r2-file", help="r2 file", required=False)
@click.option(
    "-r",
    "--ref-file",
    "ref_files",
    # multiple=True,
    cls=OptionEatAll,
    help="Reference file(s). Can be specified multiple times: -r ref1.fa -r ref2.fa. For multiple refs with single output, all mappings go to one BAM. For multiple refs with multiple outputs, maps ref1→out1, ref2→out2, etc.",
    required=False,
)
@click.option(
    "-o",
    "--output-file",
    "output_files",
    # multiple=True,
    cls=OptionEatAll,
    help="Output BAM file(s). Can be specified multiple times: -o out1.bam -o out2.bam. Must match the number of reference files unless multiple refs map to a single output.",
    required=False,
)
@click.option(
    "-u",
    "--unmap-file",
    "unmap_file",
    help="Output BAM file for unmapped reads. If specified, all unmapped reads will be written to this file instead of the regular output file(s).",
    required=False,
)
@click.option(
    "--report",
    "report_file",
    help="Output mapping summary report (HISAT2 style) to this file. Use '-' for stdout.",
    required=False,
)
@click.option(
    "--max-mismatches",
    "-m",
    type=int,
    default=10,
    help="Maximum allowed bad mismatches (wrong conversions + sequencing errors) for paired-end reads. Single-end uses half this value. (default: 10)",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    default=8,
    help="Number of worker processes for parallel mapping (default: 8). For multi-reference mapping with large references, consider reducing to 2-4 to avoid memory issues.",
)
@click.option(
    "--min-alignment-length",
    type=int,
    default=20,
    help="Minimum alignment length for filtering hits (default: 20)",
)
@click.option(
    "--min-mapping-ratio",
    type=float,
    default=0.5,
    help="Minimum mapping length ratio (mapped length / query length) for filtering hits (default: 0.5)",
)
@click.option(
    "--max-a2g-ratio",
    type=float,
    default=1.0,
    help="Maximum proportion of A to G mutations over total A (default: 1.0)",
)
@click.option(
    "--max-c2t-ratio",
    type=float,
    default=1.0,
    help="Maximum proportion of C to T mutations over total C (default: 1.0)",
)
@click.option(
    "--index-dir",
    type=click.Path(),
    multiple=True,
    help="Directory to store/load index files. Can be specified multiple times to match reference files. If not specified, uses temporary directory (default: None)",
)
@click.option(
    "--index-only",
    is_flag=True,
    default=False,
    help="Only build indices without mapping reads. Requires --index-dir to be specified. (default: False)",
)
@click.option(
    "--batch-size",
    type=int,
    default=5000,
    help="Number of reads per batch per worker (default: 5000)",
)
# Strand-specific mapping options (grouped)
@click.option(
    "--fwd-lib",
    "--forward-library",
    "library_type",
    flag_value="forward",
    default=True,
    help="[Library] Forward library orientation (default)",
)
@click.option(
    "--rev-lib",
    "--reverse-library",
    "library_type",
    flag_value="reverse",
    help="[Library] Reverse library orientation",
)
@click.option(
    "--fwd-ref",
    "--forward-reference",
    "reference_strand",
    flag_value="forward",
    help="[Reference] Only map to forward reference strand",
)
@click.option(
    "--rev-ref",
    "--reverse-reference",
    "reference_strand",
    flag_value="reverse",
    help="[Reference] Only map to reverse reference strand",
)
@click.option(
    "--dbl-ref",
    "--double-reference",
    "reference_strand",
    flag_value="double",
    default=True,
    help="[Reference] Map to both reference strands (default)",
)
def map(
    r1_file,
    r2_file,
    ref_files,
    output_files,
    unmap_file,
    report_file,
    max_mismatches,
    threads,
    min_alignment_length,
    min_mapping_ratio,
    max_a2g_ratio,
    max_c2t_ratio,
    index_dir,
    index_only,
    batch_size,
    library_type,
    reference_strand,
):
    import os

    from .mapping import map_file

    if ref_files is None:
        ref_files = []
    if output_files is None:
        output_files = []

    # Validate arguments
    if index_only:
        if not index_dir:
            click.echo(
                "❌ Error: --index-only requires --index-dir to be specified", err=True
            )
            raise click.Abort()
        if not r1_file and not r2_file and not output_files:
            # Index-only mode, these are not needed
            pass
    else:
        if not r1_file:
            click.echo("❌ Error: -1/--r1-file is required for mapping", err=True)
            raise click.Abort()
        if not output_files:
            click.echo("❌ Error: -o/--output-file is required for mapping", err=True)
            raise click.Abort()

        # Validate ref/output count matching
        # If multiple output files are specified, they must match the number of ref files
        # Single output file with multiple refs is allowed (priority-based mapping to one BAM)
        if len(output_files) > 1:
            if len(ref_files) != len(output_files):
                click.echo(
                    f"❌ Error: Number of reference files ({len(ref_files)}) must match "
                    f"number of output files ({len(output_files)}). Alternatively, use a "
                    f"single output file for priority-based mapping.",
                    err=True,
                )
                raise click.Abort()

        # Validate index-dir count matching
        if index_dir and len(index_dir) > 1 and len(ref_files) > 0:
            if len(ref_files) != len(index_dir):
                click.echo(
                    f"❌ Error: Number of reference files ({len(ref_files)}) must match "
                    f"number of index directories ({len(index_dir)}) when multiple index dirs are provided.",
                    err=True,
                )
                raise click.Abort()

        # If index-dir is provided and indices exist, ref-files can be omitted
        if not ref_files:
            if not index_dir:
                click.echo(
                    "❌ Error: --ref-file is required unless --index-dir is provided",
                    err=True,
                )
                raise click.Abort()

            # Check for BWA indices in all specified index directories
            for d in index_dir:
                # Check for standard 'ref' prefix
                if os.path.exists(os.path.join(d, "ref.orig.bwt")) and os.path.exists(
                    os.path.join(d, "ref.mk.bwt")
                ):
                    continue

                click.echo(
                    f"❌ Error: BWA indices ('ref.orig.bwt' and 'ref.mk.bwt') not found in --index-dir: {d}",
                    err=True,
                )
                raise click.Abort()

    # Determine forward library flag
    forward_library = library_type == "forward"

    # Determine orientation filter based on reference strand
    orientation_filter = None
    if reference_strand == "forward":
        orientation_filter = 1
    elif reference_strand == "reverse":
        orientation_filter = 2
    # else: reference_strand == "double", orientation_filter = None (map both)

    try:
        # Call map_file once with all refs and outputs
        # It will route reads to the appropriate output file based on which reference they mapped to
        map_file(
            r1_file,
            r2_file,
            ref_files,
            output_files,
            unmap_file,
            report_file,
            forward_library,
            max_mismatches,
            threads,
            min_alignment_length,
            min_mapping_ratio,
            max_a2g_ratio,
            max_c2t_ratio,
            index_dir,
            index_only,
            batch_size,
            orientation_filter,
        )
    except FileNotFoundError as e:
        click.echo(f"❌ {e}", err=True)
        raise click.Abort()

    if index_only:
        if len(index_dir) > 1:
            print(
                f"\n✅ Index building completed! Indices saved to: {', '.join(index_dir)}"
            )
        else:
            print(f"\n✅ Index building completed! Indices saved to: {index_dir[0]}")
    else:
        if len(output_files) > 1:
            print(
                f"\n✅ Mapping completed! Outputs saved to: {', '.join(output_files)}"
            )
        else:
            print(f"\n✅ Mapping completed! Output saved to: {output_files[0]}")


@cli.command(
    help="Annotate tsv file.",
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
    help="Group and find consenus of gene.",
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
    help="Output artifical fasta file containing consenus sequences.",
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

    # Convert 1-based to 0-based indices for meta_columns
    meta_col_index = [col - 1 for col in meta_columns]
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
        import polars as pl

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
                if len(parts) > 1:
                    try:
                        motif_weights.append(float(parts[1]))
                    except ValueError:
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
        lpad, rpad = npad.split(",")
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
    help="Reference GTF file.",
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
@click.option("--all-effects", "-a", is_flag=True, help="Output all overlapping effects.")
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
        )
    except ValueError as e:
        raise click.ClickException(str(e))


@cli.command(
    help="Annotate genomic variant effects.",
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

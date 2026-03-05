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
            "name": "Commands",
            "commands": ["prepare", "map", "liftover", "annot", "group"],
        },
    ]
}

click.rich_click.OPTION_GROUPS = {
    "coralsnake map": [
        {
            "name": "Input/Output",
            "options": ["-1", "-2", "-r", "-o", "-u"],
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


@click.group(
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
    help="Fetch genomic motif.",
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
    default=1000,
    help="Number of reads per batch per worker (default: 1000)",
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
    max_mismatches,
    threads,
    min_alignment_length,
    min_mapping_ratio,
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
        if index_dir and len(index_dir) > 1:
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
                idx_files_exist = False
                if os.path.exists(
                    os.path.join(d, "ref.orig.bwt")
                ) and os.path.exists(os.path.join(d, "ref.mk.bwt")):
                    idx_files_exist = True
                elif os.path.exists(
                    os.path.join(d, "ref1.orig.bwt")
                ) and os.path.exists(os.path.join(d, "ref1.mk.bwt")):
                    idx_files_exist = True

                if not idx_files_exist:
                    click.echo(
                        f"❌ Error: --ref-file is required because BWA indices not found in --index-dir: {d}",
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
            forward_library,
            max_mismatches,
            threads,
            min_alignment_length,
            min_mapping_ratio,
            index_dir,
            index_only,
            batch_size,
            orientation_filter,
        )
    except FileNotFoundError as e:
        click.echo(f"❌ {e}", err=True)
        raise click.Abort()

    if index_only:
        print(f"\n✅ Index building completed! Indices saved to: {index_dir}")
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


if __name__ == "__main__":
    cli()

import rich_click as click

click.rich_click.COMMAND_GROUPS = {
    "coralsnake": [
        {
            "name": "Commands",
            "commands": ["prepare", "map", "liftover"],
        },
    ]
}
click.rich_click.STYLE_OPTION = "bold green"
# click.rich_click.STYLE_COMMAND = "bold blue"


@click.group(
    invoke_without_command=False,
    help="Variant (genomic variant analysis in python)",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.version_option(None, "-v", "--version")
@click.pass_context
def cli(ctx):
    pass


@cli.command(
    help="Extract primary transcript from gtf/gff file.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("--gtf-file", "-g", "gtf_file", help="GTF file.", required=True)
@click.option("--fasta-file", "-f", "fasta_file", help="Fasta file.", required=True)
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
def prepare(
    gtf_file, fasta_file, output_file, seq_file, sanitize, with_codon, with_genename
):
    from .gtf2tx import parse_file

    parse_file(
        gtf_file, fasta_file, output_file, seq_file, sanitize, with_codon, with_genename
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
    help="Map reads to reference genome.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("-r", "--ref-file", help="reference file", required=True)
@click.option("-1", "--r1-file", help="r1 file", required=True)
@click.option("-2", "--r2-file", help="r2 file", required=True)
# if --fwd is passed (or leave empty), fwd_lib will be True, otherwise if --rev is passed, fwd_lib will be False
@click.option("--fwd/--rev", "fwd_lib", default=True, help="forward or reverse library")
def map(ref_file, r1_file, r2_file, fwd_lib):
    from .mapping import map_file

    map_file(ref_file, r1_file, r2_file, fwd_lib)


if __name__ == "__main__":
    cli()

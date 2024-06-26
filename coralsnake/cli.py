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
@click.option(
    "--gtf-file",
    "-g",
    "gtf_file",
    help="GTF file.",
    required=True,
)
@click.option(
    "--fasta-file",
    "-f",
    "fasta_file",
    help="Fasta file.",
    required=True,
)
@click.option(
    "--output-file",
    "-o",
    "output_file",
    help="Output file.",
    required=True,
)
@click.option(
    "--seq-file",
    "-s",
    "seq_file",
    help="Sequence file.",
    required=True,
)
def prepare(gtf_file, fasta_file, output_file, seq_file):
    from .gtf2tx import parse_file

    parse_file(gtf_file, fasta_file, output_file, seq_file)


@cli.command(
    help="Fetch genomic motif.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--input-bam",
    "-i",
    "input_bam",
    help="Input bam file.",
    required=True,
)
@click.option(
    "--output-bam",
    "-o",
    "output_bam",
    help="Output bam file.",
    required=True,
)
@click.option(
    "--annotation-file",
    "-a",
    "annotation_file",
    help="Annotation file.",
    required=True,
)
@click.option(
    "--faidx-file",
    "-f",
    "faidx_file",
    help="Faidx file.",
    required=True,
)
def liftover(input_bam, output_bam, annotation_file, faidx_file):
    from .tbam2gbam import convert_bam

    convert_bam(input_bam, output_bam, annotation_file, faidx_file)


@cli.command(
    help="Map reads to reference genome.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("-r", "--ref-file", help="reference file", required=True)
@click.option("-1", "--r1-file", help="r1 file", required=True)
@click.option("-2", "--r2-file", help="r2 file", required=True)
@click.option("-f", "--fwd-lib", is_flag=True, help="forward library")
def map(ref_file, r1_file, r2_file, fwd_lib):
    from .mapping import map_file

    map_file(ref_file, r1_file, r2_file, fwd_lib)


if __name__ == "__main__":
    cli()

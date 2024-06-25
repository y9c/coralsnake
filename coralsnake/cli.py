import rich_click as click


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
def extract(gtf_file, fasta_file, output_file, seq_file):
    from .gtf2tx import parse_file

    parse_file(gtf_file, fasta_file, output_file, seq_file)


@cli.command(
    help="Fetch genomic motif.",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--input",
    "-i",
    "input",
    default="-",
    help="Input position file.",
    required=False,
)
@click.option(
    "--output",
    "-o",
    "output",
    default="-",
    help="Output annotation file.",
    required=False,
)
@click.option(
    "--fasta",
    "-f",
    "fasta",
    help="reference fasta file.",
    required=True,
)
@click.option(
    "--npad",
    "-n",
    "npad",
    default="10",
    help="Number of padding base to call motif. "
    "If you want to set different left and right pads, "
    "use comma to separate them. (eg. 2,3)",
)
@click.option(
    "--with-header", "-H", help="With header line in input file.", is_flag=True
)
@click.option(
    "--columns",
    "-c",
    "columns",
    default="1,2,3",
    show_default=True,
    type=str,
    help="Sets columns for site info. (Chrom,Pos,Strand)",
)
@click.option("--to-upper", "-u", help="Convert motif to upper case.", is_flag=True)
@click.option("--wrap-site", "-w", help="Wrap motif site.", is_flag=True)
def motif(input, output, fasta, npad, with_header, columns, to_upper, wrap_site):
    from .motif import run_motif

    if "," in npad:
        lpad, rpad = npad.split(",")
    else:
        lpad, rpad = npad, npad
    # check if lpad and rpad are positive int
    # exit with error if not
    if not lpad.isdigit() or not rpad.isdigit():
        click.echo(f"Error: npad should be positive integer, not {npad}", err=True)
        exit(1)
    lpad = int(lpad)
    rpad = int(rpad)
    run_motif(
        input,
        output,
        fasta,
        lpad,
        rpad,
        with_header,
        columns,
        to_upper,
        wrap_site,
    )


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

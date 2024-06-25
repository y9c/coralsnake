rule combine_genes_fa:
    input:
        [
            os.path.expanduser(f)
            for f in config.get("spike", [])
            + [INTERNALDIR / "reference_file/transcript.fa"]
        ],
    output:
        fa=INTERNALDIR / "reference_file/genes.fa",
        fai=INTERNALDIR / "reference_file/genes.fa.fai",
    shell:
        """
        cat {input} > {output.fa}
        samtools faidx {output.fa} --fai-idx {output.fai}
        """


rule make_ref_conversion:
    input:
        INTERNALDIR / "reference_file/genes.fa",
    output:
        mk=INTERNALDIR / "reference_file/ref.MK.fa",
        km=INTERNALDIR / "reference_file/ref.KM.fa",
    shell:
        "python ../src/make_mk_km_conversion.py {input} {output.mk} {output.km}"


rule build_bwa_index:
    input:
        INTERNALDIR / "reference_file/ref.{type}.fa",
    output:
        INTERNALDIR / "reference_file/ref.{type}.amb",
    params:
        lambda wildcards: str(INTERNALDIR / f"reference_file/ref.{wildcards.type}"),
    shell:
        """
        ~/tools/bwa-mem2/bwa-mem2 index -p {params} {input}
        """


# cut adapters


rule cutadapt_SE:
    input:
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R1", "/"),
    output:
        c=temp(TEMPDIR / "trimmed_reads_SE/{sample}_{rn}_R1.fq.gz"),
        s=INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
        report="report_reads/trimming/{sample}_{rn}.json",
    params:
        library=lambda wildcards: SAMPLE2LIB[wildcards.sample],
    threads: 16
    shell:
        """
        cutseq -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule cutadapt_PE:
    input:
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R1", "/"),
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R2", "/"),
    output:
        c=[
            temp(TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_R1.fq.gz"),
            temp(TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_R2.fq.gz"),
        ],
        s=[
            INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R1.fq.gz",
            INTERNALDIR / "discarded_reads/{sample}_{rn}_tooshort_R2.fq.gz",
        ],
        report="report_reads/trimming/{sample}_{rn}.json",
    params:
        library=lambda wildcards: SAMPLE2LIB[wildcards.sample],
    threads: 16
    shell:
        """
        cutseq -t {threads} -A {params.library:q} -m 20 --auto-rc -o {output.c} -s {output.s} --json-file {output.report} {input} 
        """


rule report_cutadapt:
    input:
        [
            f"report_reads/trimming/{sample}_{rn}.json"
            for sample, d1 in SAMPLE2DATA.items()
            for rn, d2 in d1.items()
        ],
    output:
        "report_reads/report_cutadapt.html",
    shell:
        "multiqc -f -m cutadapt --no-data-dir -n {output} {input}"


# trimmed part qc


rule qc_trimmed:
    input:
        lambda wildcards: (
            TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_{rd}.fq.gz"
            if len(SAMPLE2DATA[wildcards.sample][wildcards.rn]) == 2
            else TEMPDIR / "trimmed_reads_SE/{sample}_{rn}_{rd}.fq.gz"
        ),
    output:
        html="report_reads/trimmed/{sample}_{rn}_{rd}/fastqc_report.html",
        text="report_reads/trimmed/{sample}_{rn}_{rd}/fastqc_data.txt",
        summary="report_reads/trimmed/{sample}_{rn}_{rd}/summary.txt",
    params:
        "report_reads/trimmed/{sample}_{rn}_{rd}",
    shell:
        "falco -o {params} {input}"


rule report_qc_trimmed:
    input:
        [
            f"report_reads/trimmed/{s}_{r}_{i}/fastqc_data.txt"
            for s, v in SAMPLE2DATA.items()
            for r, v2 in v.items()
            for i in v2.keys()
        ],
    output:
        "report_reads/report_trimmed.html",
    shell:
        "multiqc -f -m fastqc -n {output} {input}"


#########################
# mapping and converting
#########################


rule convert_reads:
    input:
        lambda wildcards: (
            TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_{rd}.fq.gz"
            if len(SAMPLE2DATA[wildcards.sample][wildcards.rn]) == 2
            else TEMPDIR / "trimmed_reads_SE/{sample}_{rn}_{rd}.fq.gz"
        ),
    output:
        mk=temp(TEMPDIR / "converted_reads/{sample}_{rn}_{rd}.MK.fq.gz"),
        km=temp(TEMPDIR / "converted_reads/{sample}_{rn}_{rd}.KM.fq.gz"),
    shell:
        "python ../src/make_mk_km_conversion.py {input} {output.mk} {output.km} YS ST"


rule map_to_fwd_PE:
    input:
        ref=INTERNALDIR / "reference_file/ref.MK.amb",
        f1=TEMPDIR / "converted_reads/{sample}_{rn}_R1.MK.fq.gz",
        f2=TEMPDIR / "converted_reads/{sample}_{rn}_R2.KM.fq.gz",
    output:
        temp(TEMPDIR / "stranded_map_PE/{sample}_{rn}.MK.bam"),
    threads: 20
    params:
        index=str(INTERNALDIR / "reference_file/ref.MK"),
    shell:
        """
        ~/tools/bwa-mem2/bwa-mem2 mem -t {threads} -C -Y -M -L 6,6 -U 24 {params.index} {input.f1} {input.f2} | \
            samtools view -e '!flag.read2==!flag.reverse && flag.proper_pair && !flag.unmap && !flag.munmap && cigar!~"H" && (![MC] || [MC]!~"H")' -Sb > {output}
        """


rule map_to_rev_PE:
    input:
        ref=INTERNALDIR / "reference_file/ref.KM.amb",
        f1=TEMPDIR / "converted_reads/{sample}_{rn}_R1.MK.fq.gz",
        f2=TEMPDIR / "converted_reads/{sample}_{rn}_R2.KM.fq.gz",
    output:
        temp(TEMPDIR / "stranded_map_PE/{sample}_{rn}.KM.bam"),
    params:
        index=str(INTERNALDIR / "reference_file/ref.KM"),
    threads: 20
    shell:
        """
        ~/tools/bwa-mem2/bwa-mem2 mem -t {threads} -C -Y -M -L 6,6 -U 24 {params.index} {input.f1} {input.f2} | \
            samtools view -e '!flag.read2!=!flag.reverse && flag.proper_pair && !flag.unmap && !flag.munmap && cigar!~"H" && (![MC] || [MC]!~"H")' -Sb > {output}
        """


rule map_to_fwd_SE:
    input:
        ref=INTERNALDIR / "reference_file/ref.MK.amb",
        fq=TEMPDIR / "converted_reads/{sample}_{rn}_R1.MK.fq.gz",
    output:
        temp(TEMPDIR / "stranded_map_SE/{sample}_{rn}.MK.bam"),
    threads: 20
    params:
        index=str(INTERNALDIR / "reference_file/ref.MK"),
    shell:
        """
        ~/tools/bwa-mem2/bwa-mem2 mem -t {threads} -C -Y -M -L 6,6 -U 24 {params.index} {input.fq} | samtools view -e '!flag.reverse && !flag.unmap && cigar!~"H" && (![MC] || [MC]!~"H")' -Sb > {output}
        """


rule map_to_rev_SE:
    input:
        ref=INTERNALDIR / "reference_file/ref.KM.amb",
        fq=TEMPDIR / "converted_reads/{sample}_{rn}_R1.MK.fq.gz",
    output:
        temp(TEMPDIR / "stranded_map_SE/{sample}_{rn}.KM.bam"),
    params:
        index=str(INTERNALDIR / "reference_file/ref.KM"),
    threads: 20
    shell:
        """
        ~/tools/bwa-mem2/bwa-mem2 mem -t {threads} -C -Y -M -L 6,6 -U 24 {params.index} {input.fq} | samtools view -e 'flag.reverse && !flag.unmap && cigar!~"H" && (![MC] || [MC]!~"H")' -Sb > {output}
        """


rule combine_map_join:
    input:
        lambda wildcards: [
            (
                TEMPDIR
                / (
                    "stranded_map_PE"
                    if len(SAMPLE2DATA[wildcards.sample][wildcards.rn]) == 2
                    else "stranded_map_SE"
                )
                / f"{wildcards.sample}_{wildcards.rn}.{t}.bam"
            )
            for t in ["MK"]
        ],
        # for t in (
        #     ["MK"]
        #     if SAMPLE2LIB[wildcards.sample]["strand"] == "+"
        #     else (
        #         ["KM"]
        #         if SAMPLE2LIB[wildcards.sample]["strand"] == "-"
        #         else ["MK", "KM"]
        #     )
        # )
    output:
        temp(TEMPDIR / "joined_bam/{sample}_{rn}.bam"),
    threads: 20
    shell:
        """
        samtools merge -@ {threads} -f -n -o {output} {input}
        """


rule replace_reads:
    input:
        fqs=lambda wildcards: (
            [
                TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_R1.fq.gz",
                TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_R2.fq.gz",
            ]
            if len(SAMPLE2DATA[wildcards.sample][wildcards.rn]) == 2
            else [TEMPDIR / "trimmed_reads_SE/{sample}_{rn}_R1.fq.gz"]
        ),
        bam=TEMPDIR / "joined_bam/{sample}_{rn}.bam",
    output:
        bam=temp(TEMPDIR / "replaced_bam/{sample}_{rn}.bam"),
    threads: 4
    shell:
        "python ../src/replace_reads_from_ys_tag.py -i {input.bam} -o {output.bam}"


rule combine_map_sortname:
    input:
        TEMPDIR / "replaced_bam/{sample}_{rn}.bam",
    output:
        temp(TEMPDIR / "name_sorted/{sample}_{rn}.bam"),
    threads: 20
    shell:
        "samtools sort -@ {threads} -n -o {output} {input}"


rule replace_reads_calmd:
    input:
        ref=INTERNALDIR / "reference_file/genes.fa",
        bam=TEMPDIR / "name_sorted/{sample}_{rn}.bam",
    output:
        temp(TEMPDIR / "calmd_bam/{sample}_{rn}.bam"),
    threads: 20
    shell:
        "samtools calmd -@ {threads} -b -Q {input.bam} {input.ref} > {output}"


rule add_tag:
    input:
        TEMPDIR / "calmd_bam/{sample}_{rn}.bam",
    output:
        TEMPDIR / "tagged_bam/{sample}_{rn}.bam",
    params:
        # is_reverse=lambda wildcards: "" if SAMPLE2LIB[wildcards.sample]["strand"] == "-" else "--forward-lib",
        is_reverse="--forward-lib",
    threads: 20
    shell:
        "python ../src/add_tag.py {input} {output} --threads {threads} {params.is_reverse}"


rule filter_and_sort_bam:
    input:
        TEMPDIR / "tagged_bam/{sample}_{rn}.bam",
    output:
        mapped=INTERNALDIR / "genes_bam/{sample}_{rn}.bam",
    threads: 18
    shell:
        """
        samtools sort --input-fmt-option='filter=[AP]<=5 && !flag.secondary' -@ {threads} -m 3G {input} -O BAM -o {output.mapped}
        """


# select unmapped reads


ruleorder: select_unmapped_reads_PE > select_unmapped_reads_SE


rule select_unmapped_reads_PE:
    input:
        fqs=lambda wildcards: [
            TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_R1.fq.gz",
            TEMPDIR / "trimmed_reads_PE/{sample}_{rn}_R2.fq.gz",
        ],
        bam=INTERNALDIR / "genes_bam/{sample}_{rn}.bam",
    output:
        INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
        INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R2.fq.gz",
    threads: 3
    shell:
        """
        python ../src/select_unmap.py -b {input.bam} -i {input.fqs} -o {output}
        """


rule select_unmapped_reads_SE:
    input:
        fqs=[TEMPDIR / "trimmed_reads_SE/{sample}_{rn}_R1.fq.gz"],
        bam=INTERNALDIR / "genes_bam/{sample}_{rn}.bam",
    output:
        INTERNALDIR / "discarded_reads/{sample}_{rn}_unmap_R1.fq.gz",
    threads: 14
    shell:
        """
        python ../src/select_unmap.py -b {input.bam} -i {input.fqs} -o {output}
        """


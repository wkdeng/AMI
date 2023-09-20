import os
import re

###**--- IN-FILE CONFIG ---**###
# SNAKEMAKE_FILE_DIR = config['SCRIPTS_DIR']
# GTF_PATH = config['gtf']
KALLISTO_IDX=config['kallisto_idx']
STAR_IDX=config['star_idx']
SAMPLES=config['samples']


rule summary:
    input:
        # exp=["data/kallisto_bs/{sample_name}/abundance.tsv".format(
        #     sample_name=x)
        #     for x in SAMPLES],
        # bam=["data/star/{sample_name}/Aligned.sortedByCoord.out.bam".format(
        #     sample_name=x)
        #     for x in SAMPLES],
        # exp_trimmed=["data/kallisto_trimmed/{sample_name}/abundance.tsv".format(
        #     sample_name=x)
        #     for x in SAMPLES],
        bam_trimmed=["data/star_trimmed/{sample_name}/Aligned.sortedByCoord.out.bam".format(
            sample_name=x)
            for x in SAMPLES]
        # fastqc=["data/fastqc/{sample_name}_R1_fastqc.zip".format(
        #     sample_name=x)
        #     for x in SAMPLES],
        # fastqc_trimmed=["data/fastqc_trimmed/{sample_name}_R1_fastqc.zip".format(
        #     sample_name=x)
        #     for x in SAMPLES]

rule adaptor_trim:
    input:
        r1="data/raw/{sample_name}_R1.fq.gz",
        r2="data/raw/{sample_name}_R2.fq.gz"
    output:
        r1="data/trimmed/{sample_name}_R1.fq.gz",
        r2="data/trimmed/{sample_name}_R2.fq.gz"
    log:
        "data/logs/trim/{sample_name}.log"
    params:
        out_dir = "data/trimmed",
        sample_name = "{sample_name}"
    shell:
        """
        cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o {params.out_dir}/{params.sample_name}_R1.fq.gz -p {params.out_dir}/{params.sample_name}_R2.fq.gz {input.r1} {input.r2} > {log} 2>&1"""

### alignment and preprocessing
rule star_map_trimmed:
    input:
        r1="data/trimmed/{sample_name}_R1.fq.gz",
        r2="data/trimmed/{sample_name}_R2.fq.gz"
    output:
        "data/star_trimmed/{sample_name}/Aligned.sortedByCoord.out.bam"
    log:
        "data/logs/star_trimmed/{sample_name}.log"
    params:
        reads = "data/raw/{sample_name}.fastq",
        prefix = "data/star_trimmed/{sample_name}/",
        index = STAR_IDX,
        max_hits = 100
    threads: 10
    shell:
        """
        STAR --genomeDir {params.index} \
        --readFilesIn {input.r1} {input.r2}  --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.prefix} \
        --outFilterMultimapNmax {params.max_hits} \
        --runThreadN {threads} \
        --twopassMode Basic \
        --readFilesCommand gunzip -c \
        --outStd Log >{log} 2>&1
        """

rule fastqc:
    input:
        r1="data/raw/{sample_name}_R1.fq.gz",
        r2="data/raw/{sample_name}_R2.fq.gz"
    output:
        "data/fastqc/{sample_name}_R1_fastqc.zip","data/fastqc/{sample_name}_R2_fastqc.zip"
    log:
        "data/fastqc/{sample_name}.log"
    params:
        out_dir = "data/fastqc"
    shell:
        """
        fastqc -o {params.out_dir} {input.r1} > {log} 2>&1
        fastqc -o {params.out_dir} {input.r2} > {log} 2>&1
        """

rule fastqc_trimmed:
    input:
        r1="data/trimmed/{sample_name}_R1.fq.gz",
        r2="data/trimmed/{sample_name}_R2.fq.gz"
    output:
        "data/fastqc_trimmed/{sample_name}_R1_fastqc.zip","data/fastqc_trimmed/{sample_name}_R2_fastqc.zip"
    log:
        "data/fastqc_trimmed/{sample_name}.log"
    params:
        out_dir = "data/fastqc_trimmed"
    shell:
        """
        fastqc -o {params.out_dir} {input.r1} > {log} 2>&1
        fastqc -o {params.out_dir} {input.r2} > {log} 2>&1
        """

### get gene expression
rule kallisto_quant:
    input:
        r1="data/raw/{sample_name}_R1.fq.gz",
        r2="data/raw/{sample_name}_R2.fq.gz"
    output:
        "data/kallisto_bs/{sample_name}/abundance.tsv"
    log:
        "data/logs/kallisto_bs/{sample_name}.log"
    params:
        outdir = "data/kallisto_bs/{sample_name}/",
        index = KALLISTO_IDX,
        out_log = "data/kallisto_bs/foo.txt",
        sample_name = "{sample_name}"
    threads: 5
    shell:
        """
        kallisto quant -b 100 -t {threads} -i {params.index} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1
        echo '{input} complete' >> {params.out_log}
        """

### get gene expression
rule kallisto_quant_trimmed:
    input:
        r1="data/trimmed/{sample_name}_R1.fq.gz",
        r2="data/trimmed/{sample_name}_R2.fq.gz"
    output:
        "data/kallisto_trimmed/{sample_name}/abundance.tsv"
    log:
        "data/logs/kallisto_trimmed/{sample_name}.log"
    params:
        outdir = "data/kallisto_trimmed/{sample_name}/",
        index = KALLISTO_IDX,
        out_log = "data/kallisto_trimmed/foo.txt",
        sample_name = "{sample_name}"
    threads: 5
    shell:
        """
        kallisto quant -b 100 -t {threads} -i {params.index} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1
        echo '{input} complete' >> {params.out_log}
        """

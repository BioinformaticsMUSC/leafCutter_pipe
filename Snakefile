# Minimal RNA-seq pipeline with QC and STAR alignment
# Author: Generated for RNA-seq transcript quantification

import pandas as pd
from pathlib import Path

# Configuration
configfile: "config/config.yaml"

# Load samples from TSV file
samples_df = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples_df["sample_id"].tolist()

# Determine STAR index path - use existing if provided, otherwise build new one
STAR_INDEX = config["star"]["existing_index"] if config["star"]["existing_index"] else "results/star_index"

# Helper function to get STAR index input
def get_star_index_input(wildcards):
    if config["star"]["existing_index"]:
        # If using existing index, no dependency needed (assume it exists)
        return []
    else:
        # If building new index, depend on the star_index rule
        return ["results/star_index"]

# Define the final target files
rule all:
    input:
        # FastQC reports
        expand("results/qc/fastqc/{sample}_{read}_fastqc.html", 
               sample=SAMPLES, read=["R1", "R2"]),
        expand("results/qc/fastqc/{sample}_{read}_fastqc.zip", 
               sample=SAMPLES, read=["R1", "R2"]),
        # MultiQC report
        "results/qc/multiqc/multiqc_report.html",
        # STAR alignment results
        expand("results/star/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("results/star/{sample}/Log.final.out", sample=SAMPLES),
        # BAM indices
        expand("results/star/{sample}/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES)

# FastQC quality control
rule fastqc:
    input:
        r1="data/fastq/{sample}_R1.fastq.gz",
        r2="data/fastq/{sample}_R2.fastq.gz"
    output:
        html_r1="results/qc/fastqc/{sample}_R1_fastqc.html",
        zip_r1="results/qc/fastqc/{sample}_R1_fastqc.zip",
        html_r2="results/qc/fastqc/{sample}_R2_fastqc.html",
        zip_r2="results/qc/fastqc/{sample}_R2_fastqc.zip"
    params:
        outdir="results/qc/fastqc"
    threads: 2
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastqc {input.r1} {input.r2} -o {params.outdir} -t {threads}
        """

# MultiQC to aggregate FastQC and STAR reports
rule multiqc:
    input:
        fastqc=expand("results/qc/fastqc/{sample}_{read}_fastqc.zip", 
                     sample=SAMPLES, read=["R1", "R2"]),
        star_logs=expand("results/star/{sample}/Log.final.out", sample=SAMPLES)
    output:
        "results/qc/multiqc/multiqc_report.html"
    params:
        outdir="results/qc/multiqc",
        search_dirs="results/qc/fastqc results/star"
    conda:
        "envs/qc.yaml"
    shell:
        """
        multiqc {params.search_dirs} -o {params.outdir} --force
        """

# STAR genome index (run once) - only if not using existing index
rule star_index:
    input:
        fasta=config["reference"]["genome_fasta"],
        gtf=config["reference"]["gtf"]
    output:
        directory("results/star_index")
    threads: 8
    params:
        overhang=config["star"]["read_length"] - 1
    conda:
        "envs/alignment.yaml"
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.overhang} \
             --runThreadN {threads}
        """

# STAR alignment
rule star_align:
    input:
        r1="data/fastq/{sample}_R1.fastq.gz",
        r2="data/fastq/{sample}_R2.fastq.gz",
        index=get_star_index_input
    output:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log="results/star/{sample}/Log.final.out",
        sj="results/star/{sample}/SJ.out.tab"
    params:
        prefix="results/star/{sample}/",
        index_path=STAR_INDEX,
        extra=config["star"]["extra_params"]
    threads: 8
    conda:
        "envs/alignment.yaml"
    shell:
        """
        STAR --runMode alignReads \
             --genomeDir {params.index_path} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard \
             --runThreadN {threads} \
             {params.extra}
        """

# Index BAM files for downstream analysis
rule index_bam:
    input:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "envs/alignment.yaml"
    shell:
        "samtools index {input}"
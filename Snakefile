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

# Helper function to detect FASTQ file extensions
def get_fastq_files(sample):
    """
    Automatically detect FASTQ file extensions for a given sample.
    Supports: .fastq.gz, .fq.gz, .fastq, .fq
    """
    fastq_dir = Path("data/fastq")
    
    # Possible extensions in order of preference (compressed first)
    extensions = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
    
    for ext in extensions:
        r1_file = fastq_dir / f"{sample}_R1{ext}"
        r2_file = fastq_dir / f"{sample}_R2{ext}"
        
        if r1_file.exists() and r2_file.exists():
            return {
                "r1": str(r1_file),
                "r2": str(r2_file),
                "extension": ext,
                "is_compressed": ext.endswith(".gz")
            }
    
    # If no files found, default to .fastq.gz (will cause error if files don't exist)
    return {
        "r1": f"data/fastq/{sample}_R1.fastq.gz",
        "r2": f"data/fastq/{sample}_R2.fastq.gz", 
        "extension": ".fastq.gz",
        "is_compressed": True
    }

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
        expand("results/star/{sample}/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        # Junction files
        expand("results/junctions/{sample}.junc", sample=SAMPLES),
        # Junction file list
        "results/junctions/juncfiles.txt",
        # LeafCutter clustering results
        "results/leafcutter/leafcutter_perind.counts.gz",
        "results/leafcutter/leafcutter_perind_numers.counts.gz"

# FastQC quality control
rule fastqc:
    input:
        lambda wildcards: [get_fastq_files(wildcards.sample)["r1"], 
                          get_fastq_files(wildcards.sample)["r2"]]
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
    envmodules:
        config["envmodules"]["biocontainers"],
        config["envmodules"]["fastqc"]
    shell:
        """
        fastqc {input} -o {params.outdir} -t {threads}
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
    envmodules:
        config["envmodules"]["biocontainers"],
        config["envmodules"]["multiqc"]
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
    envmodules:
        config["envmodules"]["biocontainers"],
        config["envmodules"]["star"]
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
        files=lambda wildcards: [get_fastq_files(wildcards.sample)["r1"], 
                                get_fastq_files(wildcards.sample)["r2"]],
        index=get_star_index_input
    output:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        log="results/star/{sample}/Log.final.out",
        sj="results/star/{sample}/SJ.out.tab"
    params:
        prefix="results/star/{sample}/",
        index_path=STAR_INDEX,
        extra=config["star"]["extra_params"],
        read_command=lambda wildcards: "zcat" if get_fastq_files(wildcards.sample)["is_compressed"] else "cat"
    threads: 8
    conda:
        "envs/alignment.yaml"
    envmodules:
        config["envmodules"]["biocontainers"],
        config["envmodules"]["star"]
    shell:
        """
        STAR --runMode alignReads \
             --genomeDir {params.index_path} \
             --readFilesIn {input.files} \
             --readFilesCommand {params.read_command} \
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
    envmodules:
        config["envmodules"]["biocontainers"],
        config["envmodules"]["samtools"]
    shell:
        "samtools index {input}"

# Extract splice junctions from BAM files using regtools
rule bam_to_junc:
    input:
        bam="results/star/{sample}/Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        junc="results/junctions/{sample}.junc"
    params:
        min_anchor_length=config["junction"]["min_anchor_length"],
        min_intron_size=config["junction"]["min_intron_size"],
        max_intron_size=config["junction"]["max_intron_size"]
    threads: 2
    conda:
        "envs/junction.yaml"
    envmodules:
        config["envmodules"]["biocontainers"],
        config["envmodules"]["regtools"]
    shell:
        """
        echo "Converting {input.bam} to {output.junc}"
        regtools junctions extract \
            -a {params.min_anchor_length} \
            -m {params.min_intron_size} \
            -M {params.max_intron_size} \
            {input.bam} \
            -o {output.junc}
        """

# Create a list of all junction files for downstream analysis
rule create_junc_list:
    input:
        expand("results/junctions/{sample}.junc", sample=SAMPLES)
    output:
        "results/junctions/juncfiles.txt"
    shell:
        """
        echo "Creating junction file list..."
        > {output}
        for junc_file in {input}; do
            echo "$junc_file" >> {output}
        done
        """

# LeafCutter clustering - groups junctions into intron clusters
rule leafcutter_cluster:
    input:
        juncfiles="results/junctions/juncfiles.txt"
    output:
        counts="results/leafcutter/leafcutter_perind.counts.gz",
        counts_numers="results/leafcutter/leafcutter_perind_numers.counts.gz",
        pooled="results/leafcutter/leafcutter_pooled",
        refined="results/leafcutter/leafcutter_refined",
        sortedlibs="results/leafcutter/leafcutter_sortedlibs"
    params:
        outprefix="leafcutter",
        rundir="results/leafcutter",
        max_intron_len=config["leafcutter"]["max_intron_len"],
        min_cluster_reads=config["leafcutter"]["min_cluster_reads"],
        min_cluster_ratio=config["leafcutter"]["min_cluster_ratio"],
        include_constitutive=config["leafcutter"]["include_constitutive"]
    threads: 4
    conda:
        "envs/leafcutter.yaml"
    shell:
        """
        mkdir -p {params.rundir}
        
        python scripts/leafcutter_cluster_regtools.py \
            --juncfiles {input.juncfiles} \
            --outprefix {params.outprefix} \
            --rundir {params.rundir} \
            --maxintronlen {params.max_intron_len} \
            --minclureads {params.min_cluster_reads} \
            --mincluratio {params.min_cluster_ratio} \
            --verbose \
            {params.include_constitutive}
        """
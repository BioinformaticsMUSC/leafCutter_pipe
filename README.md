# RNA-seq Pipeline: QC and STAR Alignment

A minimal Snakemake pipeline for RNA-seq transcript quantification with quality control and STAR alignment.

## Features

- **Quality Control**: FastQC for read quality assessment
- **Alignment**: STAR aligner for mapping reads to reference genome
- **Reporting**: MultiQC for aggregated quality reports
- **Reproducibility**: Conda environments for software management

## Directory Structure

```
leafCutter/
├── Snakefile                 # Main workflow file
├── config/
│   ├── config.yaml          # Pipeline configuration
│   └── samples.tsv          # Sample metadata
├── envs/
│   ├── qc.yaml             # Quality control environment
│   └── alignment.yaml       # Alignment environment
├── data/
│   └── fastq/              # Input FASTQ files
└── results/                # Pipeline outputs
    ├── qc/                 # Quality control results
    │   ├── fastqc/        # FastQC reports
    │   └── multiqc/       # MultiQC aggregated report
    └── star/              # STAR alignment results
```

## Setup

### 1. Install Dependencies

Ensure you have Snakemake and conda/mamba installed:

```bash
# Install mamba (faster than conda)
conda install -n base -c conda-forge mamba

# Install Snakemake
mamba install -c bioconda -c conda-forge snakemake
```

### 2. Configure the Pipeline

#### Update Reference Files

Edit `config/config.yaml` to specify paths to your reference genome and annotation:

```yaml
reference:
  genome_fasta: "/path/to/your/genome.fa"
  gtf: "/path/to/your/annotation.gtf"

star:
  existing_index: ""  # Leave empty to build new index
  # OR specify path to existing STAR index:
  # existing_index: "/path/to/existing/star_index"
```

**Using an Existing STAR Index**: If you already have a STAR index built for your reference genome, you can specify its path in the `existing_index` field. This will skip the time-consuming index building step and use your existing index directly.

#### Update Sample Information

Edit `config/samples.tsv` to list your samples:

```tsv
sample_id	condition
sample1	control
sample2	control
sample3	treatment
sample4	treatment
```

#### Adjust STAR Parameters

Update the read length in `config/config.yaml` to match your sequencing:

```yaml
star:
  read_length: 150  # Change this to your actual read length
```

### 3. Prepare Input Data

Place your paired-end FASTQ files in the `data/fastq/` directory following this naming convention:

```
{sample_id}_R1.{extension}
{sample_id}_R2.{extension}
```

**Supported file formats** (automatically detected):
- `.fastq.gz` (compressed FASTQ - preferred)
- `.fq.gz` (compressed FASTQ)
- `.fastq` (uncompressed FASTQ)
- `.fq` (uncompressed FASTQ)

**Examples:**
```
# Compressed files (recommended)
data/fastq/sample1_R1.fastq.gz
data/fastq/sample1_R2.fastq.gz

# Uncompressed files (also supported)  
data/fastq/sample2_R1.fastq
data/fastq/sample2_R2.fastq

# Short extensions (also supported)
data/fastq/sample3_R1.fq.gz
data/fastq/sample3_R2.fq.gz
```

The pipeline automatically detects the file format and adjusts the processing accordingly. Compressed files are preferred for storage efficiency.

## Running the Pipeline

### Dry Run (Recommended)

First, perform a dry run to check the workflow:

```bash
snakemake --dry-run --cores 1
```

### Full Pipeline

Run the complete pipeline:

```bash
# Local execution
snakemake --cores 8 --use-conda

# Or with more cores if available
snakemake --cores 16 --use-conda
```

### Specific Steps

Run only specific parts of the pipeline:

```bash
# Only QC
snakemake --cores 4 --use-conda results/qc/multiqc/multiqc_report.html

# Only STAR index
snakemake --cores 8 --use-conda results/star_index

# Only alignment for specific sample
snakemake --cores 8 --use-conda results/star/sample1/Aligned.sortedByCoord.out.bam
```

## Pipeline Rules

### 1. FastQC (`fastqc`)
- **Input**: Raw FASTQ files
- **Output**: FastQC HTML and ZIP reports
- **Purpose**: Quality assessment of raw sequencing reads

### 2. MultiQC (`multiqc`)
- **Input**: FastQC reports and STAR logs
- **Output**: Aggregated HTML report
- **Purpose**: Comprehensive quality overview

### 3. STAR Index (`star_index`)
- **Input**: Reference genome FASTA and GTF annotation
- **Output**: STAR genome index
- **Purpose**: One-time genome indexing for alignment

### 4. STAR Alignment (`star_align`)
- **Input**: FASTQ files and STAR index
- **Output**: Sorted BAM file and alignment logs
- **Purpose**: Map reads to reference genome

### 5. BAM Indexing (`index_bam`)
- **Input**: Sorted BAM file
- **Output**: BAM index file (.bai)
- **Purpose**: Enable rapid access to alignment data

## Outputs

### Quality Control
- `results/qc/fastqc/`: Individual FastQC reports for each sample
- `results/qc/multiqc/multiqc_report.html`: Comprehensive quality report

### Alignment
- `results/star/{sample}/Aligned.sortedByCoord.out.bam`: Aligned reads
- `results/star/{sample}/Log.final.out`: STAR alignment statistics
- `results/star/{sample}/SJ.out.tab`: Splice junction information

## Next Steps

This pipeline provides a solid foundation for RNA-seq analysis. You can extend it by adding:

- **Read trimming** (Trimmomatic/Cutadapt)
- **Transcript quantification** (featureCounts/RSEM/Salmon)
- **Differential expression analysis** (DESeq2/edgeR)
- **Functional enrichment analysis**
- **Alternative splicing analysis**

## Troubleshooting

### Common Issues

1. **Memory errors with STAR**: Reduce thread count or increase available memory
2. **Conda environment issues**: Clear conda cache with `conda clean --all`
3. **Path errors**: Ensure all file paths in config are absolute paths
4. **Using existing STAR index**: Make sure the existing index was built with the same STAR version and compatible parameters
5. **File format detection**: Ensure both R1 and R2 files for each sample use the same extension (e.g., both .fastq.gz or both .fastq)
6. **Missing files**: The pipeline will default to .fastq.gz if files are not found - check your file naming and extensions

### Performance Tips

- **Reuse STAR indices**: If you have multiple projects with the same reference genome, build the index once and reuse it by setting `star.existing_index` in your config
- **Index compatibility**: STAR indices are generally compatible across minor version updates, but may need rebuilding for major version changes

### Resource Requirements

- **STAR indexing**: ~30GB RAM, 30-60 minutes
- **STAR alignment**: ~8GB RAM per sample, 10-30 minutes per sample
- **FastQC**: ~1GB RAM per sample, 5-10 minutes per sample

## Version Control Best Practices

This repository includes a comprehensive `.gitignore` file that excludes:

- **Data files**: FASTQ, BAM, FASTA, GTF files (too large for git)
- **Results**: All pipeline outputs in `results/` directory
- **Temporary files**: Snakemake temporary files, logs, system files
- **Environment files**: Conda cache and Python bytecode

### What to Track in Git:
- Pipeline code (`Snakefile`)
- Configuration templates (`config/config.yaml`, `config/samples.tsv`)
- Environment definitions (`envs/*.yaml`)
- Documentation (`README.md`)

### What NOT to Track:
- Raw data files (FASTQ, reference genomes)
- Pipeline results and intermediate files
- Large binary files (STAR indices, BAM files)
- System-specific temporary files

**Note**: You may want to track your actual `config/samples.tsv` in a separate private repository if it contains sensitive sample information.

## Contact

For questions or issues, please check the Snakemake documentation or create an issue in the project repository.
# UMI C. elegans Analysis Pipeline

This repository contains scripts and tools for analyzing UMI (Unique Molecular Identifier) data from *C. elegans* samples.

## Example Results

This repository includes example results from running the full pipeline on complete FASTQ files. The following directories contain outputs that demonstrate what you can expect after completing the analysis:

- **`MultiQC/`** - Example aggregated quality control reports and metrics from the full dataset
- **`UMI_Metrics/`** - Example UMI family size histograms and summary statistics for all samples
- **`Analysis_Results_Combined/`** - Example downstream analysis results including:
  - Transcriptional error rate tables and statistics
  - Frequency distribution tables (all sites and variant sites only)
  - GO enrichment analysis results and visualizations
  - Global error rate comparisons between sample groups
- **`Variant_Analysis/`** - Example variant calling analysis results including statistical tests and visualizations

These folders are provided as references to help you understand the expected output structure and content when you run the pipeline on your own data.

## Getting Started

### Step 1: Clone the Repository

Clone this repository to your local machine:

```bash
git clone https://github.com/Prahlad-Lab/umi_celegans_analysis.git
```

### Step 2: Navigate to the Repository Directory

Change directory to the cloned repository:

```bash
cd umi_celegans_analysis
```

### Step 3: Download FASTQ Files

The sequencing data (FASTQ files) for this analysis are available from the European Nucleotide Archive (ENA) under project accession **PRJEB101605**. 

#### Option A: Use Sample Subset Files (Faster - Recommended for Testing)

For quick testing and validation of the pipeline, the repository includes pre-subsampled FASTQ files in the `Sample_Fastq/` directory. These subset files contain 100,000 reads per sample (~3-4MB per file) and allow you to run the entire pipeline in significantly less time without downloading large files from ENA.

**Benefits:**
- ✅ Already included in the repository (no download needed)
- ✅ Much faster pipeline execution (minutes instead of hours)
- ✅ Smaller disk space requirements (~50MB vs several GB)
- ✅ Ideal for testing, validation, and familiarization with the pipeline

**Files included in Sample_Fastq/:**
- N2.30min.HS.1_R1.subset.fastq.gz and N2.30min.HS.1_R2.subset.fastq.gz
- N2.30min.HS.2_R1.subset.fastq.gz and N2.30min.HS.2_R2.subset.fastq.gz
- N2.30min.HS.3_R1.subset.fastq.gz and N2.30min.HS.3_R2.subset.fastq.gz
- PRDE1.30min.HS.1_R1.subset.fastq.gz and PRDE1.30min.HS.1_R2.subset.fastq.gz
- PRDE1.30min.HS.2_R1.subset.fastq.gz and PRDE1.30min.HS.2_R2.subset.fastq.gz
- PRDE1.30min.HS.3_R1.subset.fastq.gz and PRDE1.30min.HS.3_R2.subset.fastq.gz



#### Option B: Download Full ENA Files (Required for Complete Analysis)

For complete, production-level analysis with full sequencing depth, download the full FASTQ files from ENA. The repository includes a download script that will retrieve all necessary FASTQ files:

```bash
# Make the download script executable
chmod +x ena-file-download-read_run-PRJEB101605-fastq_ftp-20251208-1709.sh

# Run the download script
./ena-file-download-read_run-PRJEB101605-fastq_ftp-20251208-1709.sh
```

This script will download 12 FASTQ files (paired-end reads for 6 samples) from the EBI FTP server. The files include:
- ERR15764772_1.fastq.gz and ERR15764772_2.fastq.gz
- ERR15764773_1.fastq.gz and ERR15764773_2.fastq.gz
- ERR15764774_1.fastq.gz and ERR15764774_2.fastq.gz
- ERR15764775_1.fastq.gz and ERR15764775_2.fastq.gz
- ERR15764776_1.fastq.gz and ERR15764776_2.fastq.gz
- ERR15764777_1.fastq.gz and ERR15764777_2.fastq.gz

**Note:** The download script uses `wget -nc` (no-clobber) which will skip files that have already been downloaded, making it safe to re-run if interrupted.

**Requirements:**
- `wget` must be installed on your system
- Sufficient disk space (FASTQ files can be several GB each)
- Stable internet connection for downloading from EBI FTP server

### Step 4: Extract Input Files

The repository includes compressed input files that need to be extracted before running the analysis:

```bash
cd input
tar -xzvf inputs.tar.gz
cd ..
```

This will extract the reference genome files and other required input data.

### Step 5: Set Up Conda Environments

Return to the main repository folder and set up the required conda environments. You have two options:

**Option 1: Combined Environment (Recommended for Simple Setup)**

```bash
./setup_conda_environments.sh all
```

This creates a single environment called `umi_celegans_analysis` with all tools. Activate it with:

```bash
conda activate umi_celegans_analysis
```

**Option 2: Separate Environments (Matches Pipeline Structure)**

```bash
./setup_conda_environments.sh separate
```

This creates individual environments for each tool (fastqc, bcftools, umi_tools, gatk4, fgbio, snpeff, samtools.v1.22, STAR).

### Step 6: Run the Analysis Pipeline

The main analysis pipeline is `UMI_analysis_pipeline_11.sh`. This script performs the complete UMI-based variant calling analysis using consensus reads.

**What the pipeline does:**
- Quality control (FastQC)
- Converts FastQ to unmapped BAM
- Extracts UMIs from reads
- Aligns reads using STAR
- Groups reads by UMI families
- **Generates consensus reads** from UMI families (using fgbio's CallMolecularConsensusReads)
- Re-aligns consensus reads to the reference genome
- Performs variant calling with GATK HaplotypeCaller on consensus BAMs
- Annotates variants with SnpEff
- Calculates allele proportions using bcftools mpileup
- Generates MultiQC summary report

**Key Features:**
The pipeline uses a consensus-based approach that:
- Groups reads by UMI sequences to identify molecular families
- Calls consensus sequences from each UMI family to reduce sequencing errors
- Produces high-quality consensus BAM files for more accurate variant calling
- Includes comprehensive quality metrics and MultiQC reporting

**Output folders created:**
The pipeline creates a base output directory with the following subdirectories. The default location is `/vscratch/grp-vprahlad/umi_celegans_consensus` for HPC environments, but this can be configured via the `BASE_OUTPUT_DIR` variable in the script.
- `fastqc/` - Quality control reports for raw sequencing data
- `FastqToUbam/` - Unmapped BAM files converted from FastQ
- `ExtractUmisFromBam/` - UMI-extracted BAM files
- `STARNoClip/` - Initial STAR alignment outputs (before consensus)
- `UMIAwareDuplicateMarkingGenomeNoClip/` - UMI-grouped BAM files
- `Consensus_Analysis/` - **Consensus read generation outputs including:**
  - Grouped BAMs by UMI families (`.fgbio_grouped.bam`)
  - Consensus unmapped BAMs (`.consensus.unmapped.bam`)
  - Re-aligned consensus BAMs (`.consensus.merged.bam`)
  - Alignment metrics for consensus reads
- `FilterBambyFamilySizeGenome/` - Family size statistics (kept for reference)
- `HaplotypeCaller/` - GATK variant calling outputs from consensus BAMs
- `Annotation/` - SnpEff annotated variants (VCF and BED formats)
- `Bcftools_Mpileup/` - bcftools pileup results
- `Allele_Proportions/` - Allele proportion calculations
- `MultiQC/` - Aggregated quality control reports

**Configuration:**

Before running the pipeline, you need to configure the following variables at the top of the `UMI_analysis_pipeline_11.sh` script:

- **`DIR_SEQS`** - Set this to the location where your FASTQ files are stored. 
  - For sample data: Use `./Sample_Fastq` directory in the repository
  - For full ENA data: Use the directory where you downloaded the FASTQ files in Step 3
  
- **`DIR_REF_INPUT`** - Set this to the location where the input files were unzipped (from Step 4).
  - Default: `./input` directory in the repository after extracting `inputs.tar.gz`
  
- **`PYTHON_SCRIPT_PATH`** - Set this to the path of the allele proportion calculation Python script.
  - A copy of this file is available in the `Python_Script/` folder as `calculate_allele_proportions_depth.py`
  - Example: `./Python_Script/calculate_allele_proportions_depth.py`

Example configuration:
```bash
DIR_SEQS="/path/to/your/fastq/files"  # or use "./Sample_Fastq" for sample data
DIR_REF_INPUT="/path/to/your/input"   # or use "./input" if extracted in repository
PYTHON_SCRIPT_PATH="./Python_Script/calculate_allele_proportions_depth.py"
```

**To run the pipeline:**

```bash
# Activate the conda environment
conda activate umi_celegans_analysis

# Run the pipeline
./UMI_analysis_pipeline_11.sh
```

**Note:** This pipeline is designed to run on an HPC cluster with SLURM. If running locally, you may need to:
- Remove or modify the SLURM directives at the top of the script (the `#SBATCH` lines)
- Adjust the input/output directory paths (such as `BASE_OUTPUT_DIR`, `DIR_SEQS`, and `DIR_REF_INPUT`) to match your local file system
- Ensure you have at least 100GB of memory and 32 CPU cores available for optimal performance

**Advanced Options:**
- The pipeline includes smart checkpoint/resume functionality - it will skip steps that have already completed successfully
- Set `CLEANUP_INTERMEDIATES="true"` (defined near the top of the script) to automatically remove intermediate files and save disk space
- Adjust parallelization parameters (e.g., `parallel -j 6`) based on your available resources

**Testing Information:** The pipeline was tested with UB CCR (University at Buffalo Center for Computational Research), using 32 cores and 100GB of memory with a 48-hour time limit.

### Step 7: Calculate UMI Metrics (Optional Quality Control)

After running the main analysis pipeline, you can calculate detailed UMI metrics to assess the quality and complexity of your UMI-based sequencing data. The `run_umi_metrics.sh` script generates family size histograms and summary statistics for each sample.

**What the script does:**
- Generates UMI family size histograms showing the distribution of reads per UMI family
- Calculates key quality metrics:
  - **Total Reads Processed**: Total number of sequenced reads
  - **Total UMI Families**: Number of unique molecular families identified
  - **Mean Family Size**: Average number of reads per UMI family
  - **Saturation**: Ratio of unique molecules to total reads (lower values indicate higher duplication/lower library complexity)
  - **Singleton Families (%)**: Percentage of UMI families with only one read (may indicate under-sequencing or technical artifacts)

**When to use it:**
- After completing the main pipeline (particularly after the UMI-aware duplicate marking step)
- To assess library complexity and sequencing depth
- To evaluate whether you have sufficient UMI diversity for your analysis
- To compare quality metrics across different samples or experimental conditions

**Prerequisites:**
- The main analysis pipeline (`UMI_analysis_pipeline_11.sh`) must have completed successfully
- The `UMIAwareDuplicateMarkingGenomeNoClip/` output directory must contain coordinate-sorted BAM files
- The `umi_celegans_analysis` conda environment must be activated

**Configuration:**

Before running the script, you may need to adjust the following variables at the top of `run_umi_metrics.sh`:

```bash
# Base output directory (should match your pipeline output location)
# Replace with your actual output directory path
BASE_OUTPUT_DIR="/path/to/your/output"  # e.g., "/vscratch/grp-vprahlad/umi_celegans_consensus"

# Sample names (update to match your samples)
# Example: sample_list=("sample1" "sample2" "sample3")
sample_list=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" 
             "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")
```

**To run the script:**

```bash
# Activate the conda environment
conda activate umi_celegans_analysis

# Run the metrics script
./run_umi_metrics.sh
```

**For HPC/SLURM environments:**
```bash
# Submit as a batch job
sbatch run_umi_metrics.sh
```

**Output files:**

The script creates a `UMI_Metrics/` directory in your output folder with the following files for each sample:

1. **`{sample}.family_size_histogram.txt`** - Tab-separated file showing:
   - Column 1: Family size (number of reads per UMI)
   - Column 2: Count (number of families with that size)

2. **`{sample}.summary_stats.txt`** - Text file with calculated metrics summary

**Example output interpretation:**

```
--- Metrics Summary for N2.30min.HS.1 ---
Total Reads Processed: 45823791
Total UMI Families:    12547823
Mean Family Size:      3.65
Saturation (Unique/Total): 0.274
Singleton Families (%): 35.2
```

**Interpreting the metrics:**
- **Higher Mean Family Size** (e.g., >3): Good UMI coverage, multiple reads per molecule
- **Lower Saturation** (e.g., <0.3): Higher library complexity, more unique molecules relative to reads
- **Lower Singleton %** (e.g., <40%): Better sequencing depth, most molecules have multiple reads
- **Higher Singleton %** (e.g., >60%): May indicate under-sequencing or need for deeper coverage

**Performance considerations:**
- The script processes samples in parallel (default: 6 samples at a time)
- Memory requirement: ~8GB per sample (controlled by `-Xmx8g` flag)
- Runtime: Typically 30 minutes to 2 hours depending on BAM file size
- Adjust parallelization in the script: Change `parallel -j 6` to match available cores
  - Example: For a system with 16 cores and 64GB RAM, you could use `parallel -j 8` (8 cores × 8GB = 64GB)
  - Rule of thumb: Set `-j` value such that `(number of parallel jobs × 8GB)` doesn't exceed your total available memory

**Troubleshooting:**

If the script fails:
1. **Check input files exist:**
   ```bash
   ls -lh ${BASE_OUTPUT_DIR}/UMIAwareDuplicateMarkingGenomeNoClip/*.MergeBamAlignment.coordinate_sorted.bam
   ```

2. **Verify conda environment:**
   ```bash
   conda activate umi_celegans_analysis
   fgbio --version  # Should show fgbio version 2.5.0
   ```

3. **Check available memory:**
   - If you get "OutOfMemoryError", increase the `-Xmx8g` value in the script (e.g., to `-Xmx16g`)

4. **For local runs (non-SLURM):**
   - Remove or comment out the `#SBATCH` directives at the top of the script
   - Ensure you have sufficient memory (at least 8GB per parallel job)

### Step 8: Set Up R Environment for Variant Analysis

To perform downstream variant analysis using R scripts, you need to set up an R environment with the required packages:

```bash
./setup_conda_r_variant_analysis.sh
```

This creates a conda environment called `r_variant_analysis` with all necessary R packages including:
- **Data manipulation and visualization**: tidyverse, scales, ggpubr, cowplot, extrafont, plotly, ggprism, patchwork, svglite
- **VCF and genomics file handling**: vcfR
- **Statistical analysis**: rstatix, broom
- **Document generation**: knitr, kableExtra, officer, flextable, readxl
- **Euler diagrams**: eulerr
- **Package installation tools**: remotes - For installing packages from GitHub
- **Bioconductor genomics packages**: 
  - rtracklayer, plyranges - For working with genomic ranges and tracks
  - GenomicFeatures, GenomicRanges - For genomic feature manipulation
  - DOSE, clusterProfiler - For functional enrichment analysis
  - org.Ce.eg.db - C. elegans genome annotation database

**Important:** After creating the environment, you need to install one additional package from GitHub:

```bash
# Activate the R environment
conda activate r_variant_analysis
```

### Step 9: Run R Analysis Scripts

The `R_Scripts/` directory contains scripts for downstream analysis of consensus-based variant calling results:

**Prerequisites:**
- The main UMI analysis pipeline must have completed successfully
- The `r_variant_analysis` conda environment must be activated
- You may need to update file paths in the R scripts to match your output directory structure

**1. Variant Genome Analysis** (`R_Scripts/Haplotypecaller.Variants.Genome.Analysis_consensus_12.02.05.2025.qmd` or `.R`)

This analysis is available in two formats:
- **Quarto/R Markdown notebook** (`.qmd`): For interactive analysis with rendered HTML output
- **Plain R script** (`.R`): For command-line execution without Quarto dependency

Both versions perform comprehensive analysis of variants from consensus reads called by GATK HaplotypeCaller. They compare variant profiles between N2 (wild-type) and PRDE1 (mutant) samples to assess transcriptional error rates.

**What it does:**
- Reads and processes annotated VCF files from both sample groups (N2 and PRDE1)
- Extracts variant information including allele depth (AD), read depth (DP), and SnpEff annotations (IMPACT)
- Filters variants based on depth and functional impact
- Performs statistical comparisons (Wilcoxon and t-tests) between N2 and PRDE1 samples
- Generates visualizations and tables for variant analysis

**Input files:** 
- Annotated VCF files from the `Annotation/` directory (generated by the main UMI analysis pipeline in Step 6):
  - `N2.30min.HS.1.consensus.ann.vcf`, `N2.30min.HS.2.consensus.ann.vcf`, `N2.30min.HS.3.consensus.ann.vcf`
  - `PRDE1.30min.HS.1.consensus.ann.vcf`, `PRDE1.30min.HS.2.consensus.ann.vcf`, `PRDE1.30min.HS.3.consensus.ann.vcf`
- **Note:** These files are automatically created by the variant annotation step of the main pipeline
- GTF annotation file (optional for gene quantile and intron/exon analysis): `Caenorhabditis_elegans.WBcel235.114.gtf` (from input files)

**Configuration:** Before running, you may want to update the GTF file path in the script. Open the script file (`.R` or `.qmd`) and modify the `gtf_path` variable in the "DATA LOADING: Define File Paths" section:
```r
# GTF File Path Configuration
# Set the path to your GTF annotation file
gtf_path <- "input/Caenorhabditis_elegans.WBcel235.114.gtf"
```

**How to run:**
```bash
# Activate the R environment
conda activate r_variant_analysis

# Option 1: Open in RStudio and run interactively (Quarto version)
# Open RStudio and then open the file:
# File -> Open File -> R_Scripts/Haplotypecaller.Variants.Genome.Analysis_consensus_12.02.05.2025.qmd

# Option 2: Render to HTML using Quarto (requires Quarto installation)
quarto render R_Scripts/Haplotypecaller.Variants.Genome.Analysis_consensus_12.02.05.2025.qmd

# Option 3: Run the plain R script using Rscript (no Quarto required)
Rscript R_Scripts/Haplotypecaller.Variants.Genome.Analysis_consensus_12.02.05.2025.R
```

**Output:** 
- HTML notebook with interactive visualizations
- Statistical test results comparing variant patterns between sample groups
- Tables showing variant impacts and functional classifications

---

**2. Transcriptional Error Analysis** (`R_Scripts/transcriptional_error_analysis_4.R`)

This comprehensive R script performs multi-level analysis of transcriptional errors using consensus read data. It integrates variant data with genomic features and 22G piRNA targets.

**What it does:**
- **Global transcriptional error analysis**: Calculates error rates across all samples with statistical tests (t-test and KS-test)
- **Gene-level analysis**: Identifies genes with elevated transcriptional errors
- **Frequency distribution analysis - All Sites**: Generates allele frequency tables for all genomic positions with 0.01 bin intervals
- **Frequency distribution analysis - Variant Sites Only**: Generates allele frequency tables specifically for variant positions with 0.05 bin intervals
- **22G piRNA overlap analysis**: Examines transcriptional errors in genes targeted by 22G piRNAs (optional, requires 22G data)
- **Functional enrichment**: Performs GO term and pathway enrichment for affected genes
- **Visualization**: Creates comprehensive plots including boxplots, Euler diagrams, and frequency distributions

**Input files:**
- Per-position allele proportion files from `Allele_Proportions/` directory (`.consensus.per_position_results.tsv`) - generated by the bcftools mpileup step in the main pipeline
- GTF annotation file: `Caenorhabditis_elegans.WBcel235.114.gtf` (from input files)

**Configuration:** Before running, update the following variables in the "CONFIGURATION & SETUP" section of the script:
```r
# Set your working directory (where the pipeline output is located)
setwd("/path/to/your/pipeline/output")

# Path to allele proportion files
consensus_path <- "/path/to/your/pipeline/output/Allele_Proportions"

# GTF annotation file (should be in your input directory)
gtf_file <- "Caenorhabditis_elegans.WBcel235.114.gtf"


```

**How to run:**
```bash
# Activate the R environment
conda activate r_variant_analysis

# Run the R script
Rscript R_Scripts/transcriptional_error_analysis_4.R

# Or run interactively in R
# First start R:
R
# Then inside R, run:
source("R_Scripts/transcriptional_error_analysis_4.R")
```

**Output:** 
- `Analysis_Results_Combined/` directory containing:
  - Transcriptional error rate tables (global and gene-level):
    - `Global_Error_Stats.csv` - Per-sample error rates
    - `Gene_Specific_Error_Rates.csv` - Error rates for individual genes
  - Frequency distribution tables (CSV format):
    - `Frequency_Binned_Allele_Proportions.csv` - All sites with 0.01 bin intervals
    - `Frequency_Binned_Allele_Proportions_Variant_Only.csv` - Variant sites only with 0.05 bin intervals
  - Statistical test results:
    - `Global_Error_Ttest_Result.csv` - T-test comparing N2 vs PRDE1
    - `Global_Error_KS_Test_Result.csv` - Kolmogorov-Smirnov test results
  - GO enrichment analysis results (`GO_Results_*.csv` and `GO_Enrichment_*.png`)
  - Visualizations:
    - `Global_Error_Rate.png` - Boxplot of error rates by group

**Key features:**
- Filters data by minimum depth (10 reads) and excludes fixed variants (Proportion = 1.0) and heterozygous positions (Proportion = 0.5)
- Performs both global and gene-specific transcriptional error analysis
- Generates two types of frequency distributions:
  - Fine-grained (0.01 bins) for all genomic positions including reference sites (Proportion = 0)
  - Coarser (0.05 bins) focusing only on variant sites (Proportion > 0)
- Integrates multiple data types (variants, genomic features, optional piRNA targets)
- Exports results in CSV format for easy analysis and publication

---

**Notes:**
- Both scripts require significant computational resources and may take several hours to complete depending on data size
- Ensure all input files exist and paths are correctly specified before running
- The scripts generate output directories automatically if they don't exist

## Prerequisites

Before running the analysis pipeline, you need to install the required bioinformatics tools using conda.

### Required Tools

The pipeline requires the following tools from the bioconda channel:

- **fastqc** = 0.12.1 - Quality control for sequencing data
- **bcftools** = 1.22 - Utilities for variant calling and manipulating VCFs/BCFs
- **umi_tools** = 1.1.6 - Tools for handling Unique Molecular Identifiers
- **gatk4** = 4.6.1.0 - Genome Analysis Toolkit for variant discovery
- **fgbio** = 2.5.0 - Tools for working with genomic and sequencing data
- **snpeff** = 5.2 - Genetic variant annotation and effect prediction
- **samtools** = 1.22.1 - Tools for manipulating SAM/BAM files
- **STAR** = 2.7.11b - RNA-seq aligner

## Installation

### Conda Installation

If you don't have conda installed, download and install Miniconda or Anaconda:

- **Miniconda**: https://docs.conda.io/en/latest/miniconda.html
- **Anaconda**: https://www.anaconda.com/products/distribution

For detailed setup instructions, see the [Getting Started](#getting-started) section above.

### Alternative Manual Installation

You can also manually create the environment using the provided YAML file:

```bash
conda env create -f environment.yml
conda activate umi_celegans_analysis
```

Or install tools individually:

```bash
conda create -n fastqc -c bioconda -c conda-forge fastqc=0.12.1 -y
conda create -n bcftools -c bioconda -c conda-forge bcftools=1.22 -y
conda create -n umi_tools -c bioconda -c conda-forge umi_tools=1.1.6 -y
conda create -n gatk4 -c bioconda -c conda-forge gatk4=4.6.1.0 -y
conda create -n fgbio -c bioconda -c conda-forge fgbio=2.5.0 -y
conda create -n snpeff -c bioconda -c conda-forge snpeff=5.2 -y
conda create -n samtools.v1.22 -c bioconda -c conda-forge samtools=1.22.1 -y
conda create -n STAR -c bioconda -c conda-forge star=2.7.11b -y
```

## Data Availability

The FASTQ files required to run this analysis pipeline are available from the European Nucleotide Archive (ENA) under project accession **PRJEB101605**. 

For detailed instructions on downloading the FASTQ files, see [Step 3: Download FASTQ Files](#step-3-download-fastq-files) in the Getting Started section. The repository includes a download script (`ena-file-download-read_run-PRJEB101605-fastq_ftp-20251208-1709.sh`) that automates the download of all required FASTQ files from the EBI FTP server.

## Usage

For complete setup and usage instructions, see the [Getting Started](#getting-started) section above.


## Pipeline Overview

The `UMI_analysis_pipeline_11.sh` performs comprehensive UMI-based variant calling using a consensus read approach. The pipeline workflow consists of the following steps:

### Pre-processing Steps (1-8):
1. **FastQC** - Quality control of raw sequencing data (R1 and R2 reads)
2. **FastqToUbam** - Convert paired-end FastQ files to unmapped BAM format with read group information
3. **ExtractUMIs** - Extract UMI sequences from reads using fgbio (8bp UMI on R1, 6bp soft-clipped spacer)
4. **Sort UMI BAM** - Sort extracted BAM files by query name for downstream processing
5. **STAR Alignment** - Initial alignment of UMI-extracted reads to the reference genome (two-pass mode)
6. **Sort Aligned BAM** - Sort aligned BAM files by query name
7. **MergeBamAlignment** - Merge unmapped and aligned BAMs to restore UMI tags and read metadata
8. **Sort Merged BAM** - Coordinate-sort merged BAMs with index creation

### Consensus Generation (Step 9):
9. **Consensus Analysis** - **Key innovation of this pipeline:**
   - **Group Reads by UMI** - Group aligned reads into molecular families based on UMI sequences (fgbio GroupReadsByUmi with adjacency strategy)
   - **Call Consensus Reads** - Generate consensus sequences from each UMI family to reduce sequencing errors (fgbio CallMolecularConsensusReads)
   - **Re-align Consensus Reads** - Align consensus reads back to the reference genome using STAR
   - **Merge Consensus BAMs** - Merge unmapped and aligned consensus BAMs to create final consensus BAM files
   - **Collect Metrics** - Generate alignment summary metrics for consensus reads

### Variant Calling & Analysis (Steps 10-16):
10. **SplitNCigarReads** - Split reads that span splice junctions for variant calling
11. **Base Quality Recalibration (BQSR)** - Recalibrate base quality scores using known variants
12. **HaplotypeCaller** - Call variants from consensus BAM files using GATK HaplotypeCaller
13. **Variant Annotation** - Annotate variants with functional effects using SnpEff (impact predictions, gene annotations)
14. **Bcftools Mpileup** - Generate pileup data for allele proportion analysis
15. **Allele Proportion Calculation** - Calculate allele frequencies and depth statistics using Python script
16. **MultiQC** - Aggregate quality control metrics and generate comprehensive HTML report

### Key Features:
- **Consensus-based approach**: Uses molecular consensus sequences from UMI families to minimize sequencing errors and improve variant calling accuracy
- **Parallel processing**: Leverages GNU parallel for efficient multi-sample processing
- **Smart checkpointing**: Automatically skips completed steps for efficient pipeline resumption
- **Comprehensive QC**: Includes FastQC reports, alignment metrics, and MultiQC aggregation

## Directory Structure

### Repository Structure:
```
.
├── UMI_analysis_pipeline_11.sh         # Main consensus-based UMI analysis pipeline
├── run_umi_metrics.sh                  # UMI family size metrics calculation script
├── setup_conda_environments.sh         # Main bioinformatics tools environment setup
├── setup_conda_r_variant_analysis.sh   # R environment setup for downstream analysis
├── environment.yml                     # Conda environment specification
├── README.md                           # This documentation file
├── QUICKSTART.md                       # Quick start guide
│
├── input/                              # Reference genome files and annotations
│   ├── Caenorhabditis_elegans.WBcel235.dna.toplevel.fa      # Reference genome
│   ├── Caenorhabditis_elegans.WBcel235.dna.toplevel.dict    # Sequence dictionary
│   ├── Caenorhabditis_elegans.WBcel235.114.gtf              # Gene annotations
│   └── random_subset.vcf.gz                                 # Known variants for BQSR
│
├── Sample_Fastq/                       # Pre-subsampled FASTQ files for testing (100K reads)
│   ├── N2.30min.HS.1_R1.subset.fastq.gz
│   ├── N2.30min.HS.1_R2.subset.fastq.gz
│   └── ... (additional sample pairs)
│
├── Python_Script/                      # Python utilities
│   └── calculate_allele_proportions_depth.py  # Allele proportion calculator
│
├── R_Scripts/                          # R analysis scripts for downstream variant analysis
│   ├── Haplotypecaller.Variants.Genome.Analysis_consensus_12.02.05.2025.qmd  # Quarto version
│   ├── Haplotypecaller.Variants.Genome.Analysis_consensus_12.02.05.2025.R    # Plain R script version
│   └── transcriptional_error_analysis_4.R
│
├── Analysis_Results_Combined/          # Example R analysis outputs
├── Variant_Analysis/                   # Example variant analysis results
├── MultiQC/                           # Example MultiQC reports
└── UMI_Metrics/                       # Example UMI family metrics
```

### Pipeline Output Structure:
The `UMI_analysis_pipeline_11.sh` creates the following output directories (default location: `/vscratch/grp-vprahlad/umi_celegans_consensus` or configured via `BASE_OUTPUT_DIR`):

```
BASE_OUTPUT_DIR/
├── fastqc/                            # FastQC quality control reports
│   ├── {sample}_R1_fastqc.html
│   └── {sample}_R2_fastqc.html
│
├── FastqToUbam/                       # Unmapped BAM files with metadata
│   └── {sample}.unmapped.bam
│
├── ExtractUmisFromBam/                # UMI-extracted BAM files
│   ├── {sample}.extractUMIs.out.bam
│   └── {sample}.querysort.extractUMIs.out.bam
│
├── STARNoClip/                        # Initial STAR alignment outputs
│   ├── {sample}.Aligned.out.bam
│   └── {sample}.querysort.Aligned.out.bam
│
├── UMIAwareDuplicateMarkingGenomeNoClip/  # Merged and sorted BAMs with UMI tags
│   ├── {sample}.MergeBamAlignment.bam
│   └── {sample}.MergeBamAlignment.coordinate_sorted.bam
│
├── Consensus_Analysis/                # ⭐ Consensus read generation outputs
│   ├── {sample}.fgbio_grouped.bam               # UMI-grouped reads
│   ├── {sample}.consensus.unmapped.bam          # Consensus unmapped reads
│   ├── {sample}.consensus.unmapped.querysort.bam
│   ├── {sample}.consensus.Aligned.out.bam       # Consensus aligned reads
│   ├── {sample}.consensus.querysort.Aligned.out.bam
│   ├── {sample}.consensus.merged.bam            # Final consensus BAM
│   └── {sample}.CONSENSUS_alignment_metrics.txt # Consensus QC metrics
│
├── FilterBambyFamilySizeGenome/       # Family size statistics (reference only)
│
├── HaplotypeCaller/                   # Variant calling from consensus BAMs
│   ├── {sample}.consensus.split.bam
│   ├── {sample}.consensus.recal.table
│   ├── {sample}.consensus.bqsr.bam             # Recalibrated consensus BAM
│   └── {sample}.consensus.vcf.gz               # Called variants
│
├── Annotation/                        # SnpEff variant annotations
│   ├── {sample}.consensus.ann.vcf              # Annotated VCF
│   ├── {sample}.consensus.ann.bed              # BED format annotations
│   ├── {sample}.consensus_summary.html         # Annotation statistics (HTML)
│   └── {sample}.consensus_summary.csv          # Annotation statistics (CSV)
│
├── Bcftools_Mpileup/                  # Pileup data for allele proportions
│   └── {sample}.consensus.pileup
│
├── Allele_Proportions/                # Allele frequency calculations
│   ├── {sample}.consensus.per_position_results.tsv  # Position-level data
│   └── {sample}.consensus.summary_totals.txt        # Summary statistics
│
├── UMI_Metrics/                       # UMI family size distributions (optional)
│   ├── {sample}.family_size_histogram.txt
│   └── {sample}.summary_stats.txt
│
└── MultiQC/                           # Aggregated QC reports
    └── multiqc_report.html
```

## Support

For questions or issues, please open an issue on the GitHub repository.

## License

Please refer to the repository license file for usage terms.

# UMI C. elegans Analysis Pipeline

This repository contains scripts and tools for analyzing UMI (Unique Molecular Identifier) data from *C. elegans* samples.

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

### Step 3: Extract Input Files

The repository includes compressed input files that need to be extracted before running the analysis:

```bash
cd input
tar -xzvf inputs.tar.gz
cd ..
```

This will extract the reference genome files and other required input data.

### Step 4: Set Up Conda Environments

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

### Step 5: Run the Analysis Pipeline

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

### Step 6: Set Up R Environment for Variant Analysis

To perform downstream variant analysis using R scripts, you need to set up an R environment with the required packages:

```bash
./setup_conda_r_variant_analysis.sh
```

This creates a conda environment called `r_variant_analysis` with all necessary R packages including:
- tidyverse, scales, ggpubr, cowplot for data manipulation and visualization
- vcfR for VCF file handling
- rstatix for statistical analysis
- Bioconductor packages (rtracklayer, plyranges)

Activate the R environment:

```bash
conda activate r_variant_analysis
```

### Step 7: Run R Analysis Scripts

The `R_Scripts/` directory contains scripts for downstream analysis:

**1. Allelic Proportion Analysis** (`R_Scripts/allelic_proportion_analysis_depth.r`)

This script analyzes allele proportions from variant calling results and creates visualization plots.

```bash
# Activate the R environment
conda activate r_variant_analysis

# Run the R script
Rscript R_Scripts/allelic_proportion_analysis_depth.r
```

Or run interactively in R/RStudio:
```r
source("R_Scripts/allelic_proportion_analysis_depth.r")
```

**Input:** Uses data from `Family_size1_allele_proportions_depth/` directory
**Output:** Histogram/density plots and violin plots showing allele proportion distributions

**2. Variant Genome Analysis** (`R_Scripts/Haplotypecaller.Variants.Genome.Analysis_singletons.R`)

This script performs comprehensive analysis of variants called by GATK HaplotypeCaller.

```bash
# Activate the R environment
conda activate r_variant_analysis

# Run the R script in R/RStudio
R
> source("R_Scripts/Haplotypecaller.Variants.Genome.Analysis_singletons.R")
```

**Input:** Uses annotated VCF files from `Family_size1_annotation/` directory
**Output:** Variant analysis results, statistical comparisons, and visualizations

**Note:** You may need to update file paths in the R scripts to match your output directory structure.

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

The fastq.gz files required to run this analysis pipeline will be made available through the Gene Expression Omnibus (GEO). The GEO accession number will be provided soon.

## Usage

For complete setup and usage instructions, see the [Getting Started](#getting-started) section above.

### Running Individual Pipeline Steps

Individual pipeline steps are available in the `Separate_Scripts/` directory if you need to run specific parts of the analysis:
- `task.FastqToUbam.sh` - Convert FastQ to unmapped BAM
- `task.ExtractUMIs.sh` - Extract UMIs from BAM files
- `task.STARNoClip.sh` - STAR alignment
- `task.GroupByUMIsGenomeneClip.sh` - Group reads by UMI
- And more...

## Pipeline Overview

The pipeline performs the following steps:

1. **FastQC** - Quality control of raw sequencing data
2. **FastqToUbam** - Convert FastQ files to unmapped BAM format
3. **ExtractUMIs** - Extract UMI sequences from reads
4. **STAR Alignment** - Align reads to the reference genome
5. **UMI Grouping** - Group reads by their UMI sequences
6. **BAM Filtering** - Filter BAM files by family size (singletons, families)
7. **GATK BQSR Pipeline** - Base quality score recalibration
8. **Variant Calling** - Call variants using GATK HaplotypeCaller
9. **Variant Annotation** - Annotate variants using SnpEff
10. **Allele Proportion Calculation** - Calculate allele frequencies

## Directory Structure

```
.
├── Complete_Analysis/           # Complete pipeline scripts
├── Separate_Scripts/            # Individual task scripts
├── R_Scripts/                   # R analysis scripts
├── input_broad/                 # Reference files
├── Family_size1_allele_proportions_depth10/  # Output data
├── Family_size1_annotation/     # Variant annotations
├── environment.yml              # Conda environment specification
├── setup_conda_environments.sh  # Setup script
└── README.md                    # This file
```

## Support

For questions or issues, please open an issue on the GitHub repository.

## License

Please refer to the repository license file for usage terms.

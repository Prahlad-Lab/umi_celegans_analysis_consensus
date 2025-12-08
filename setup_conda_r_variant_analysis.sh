#!/bin/bash

# ---
# Conda Environment Setup for R Variant Analysis
#
# This script creates a new conda environment named 'r_variant_analysis'
# with all the required R packages (from CRAN and Bioconductor)
# to run all R scripts in the R_Scripts directory.
# ---

# Define the name of the environment
ENV_NAME="r_variant_analysis"

echo "--- Creating Conda environment: $ENV_NAME ---"

# Create the environment with all packages
# We specify the channels in order of priority:
# 1. conda-forge: For R base and most CRAN packages
# 2. bioconda: For Bioconductor packages
# 3. defaults: For any remaining dependencies
conda create -n $ENV_NAME \
  -c conda-forge \
  -c bioconda \
  -c defaults \
  'r-base>=4.2' \
  'r-tidyverse' \
  'r-scales' \
  'r-ggpubr' \
  'r-cowplot' \
  'r-extrafont' \
  'r-plotly' \
  'r-ggprism' \
  'r-vcfr' \
  'r-rstatix' \
  'r-knitr' \
  'r-kableextra' \
  'r-officer' \
  'r-flextable' \
  'r-broom' \
  'r-eulerr' \
  'r-patchwork' \
  'r-readxl' \
  'r-svglite' \
  'bioconductor-rtracklayer' \
  'bioconductor-plyranges' \
  'bioconductor-genomicfeatures' \
  'bioconductor-genomicranges' \
  'bioconductor-dose' \
  'bioconductor-clusterprofiler' \
  'bioconductor-org.ce.eg.db' \
  -y

# Check if creation was successful
if [ $? -eq 0 ]; then
  echo ""
  echo "--- Environment '$ENV_NAME' created successfully. ---"
  echo ""
  echo "To activate this environment, run:"
  echo "  conda activate $ENV_NAME"
  echo ""
  echo "IMPORTANT: One additional package needs to be installed from Bioconductor:"
  echo "  - wbData (not available in conda)"
  echo ""
  echo "After activating the environment, install it with:"
  echo "  Rscript -e \"if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')\""
  echo "  Rscript -e \"BiocManager::install('wbData')\""
  echo ""
  echo "To deactivate the environment, run:"
  echo "  conda deactivate"
  echo ""
else
  echo ""
  echo "--- ERROR: Failed to create Conda environment '$ENV_NAME'. ---"
  echo "Please check your Conda installation and for any package conflicts."
  echo ""
fi

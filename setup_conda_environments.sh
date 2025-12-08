#!/bin/bash

# ========================================================================
# Setup script for creating conda environments for UMI C. elegans analysis
# ========================================================================
# This script creates conda environments with the required bioinformatics
# tools from the bioconda channel.
# 
# Usage:
#   ./setup_conda_environments.sh [option]
#
# Options:
#   all      - Create one combined environment with all tools (default)
#   separate - Create separate environments for each tool
#   help     - Show this help message
# ========================================================================

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show help
show_help() {
    echo "========================================================================="
    echo "Setup script for UMI C. elegans Analysis Conda Environments"
    echo "========================================================================="
    echo ""
    echo "Usage: ./setup_conda_environments.sh [option]"
    echo ""
    echo "Options:"
    echo "  all       - Create one combined environment with all tools (default)"
    echo "  separate  - Create separate environments for each tool"
    echo "  help      - Show this help message"
    echo ""
    echo "Required tools and versions:"
    echo "  - fastqc = 0.12.1"
    echo "  - bcftools = 1.22"
    echo "  - umi_tools = 1.1.6"
    echo "  - gatk4 = 4.6.1.0"
    echo "  - fgbio = 2.5.0"
    echo "  - snpeff = 5.2"
    echo "  - samtools = 1.22.1"
    echo "  - STAR = 2.7.11b"
    echo "  - multiqc = 1.32"
    echo ""
    echo "All tools will be installed from the bioconda channel."
    echo "========================================================================="
}

# Function to check if conda is installed
check_conda() {
    if ! command -v conda &> /dev/null; then
        print_error "conda is not installed or not in PATH"
        print_info "Please install Miniconda or Anaconda first:"
        print_info "https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    print_info "Found conda: $(conda --version)"
}

# Function to create combined environment
create_combined_environment() {
    print_info "Creating combined environment 'umi_celegans_analysis' with all tools..."
    
    if conda env list | grep -q "^umi_celegans_analysis "; then
        print_warning "Environment 'umi_celegans_analysis' already exists"
        read -p "Do you want to remove and recreate it? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_info "Removing existing environment..."
            conda env remove -n umi_celegans_analysis -y
        else
            print_info "Skipping environment creation"
            return
        fi
    fi
    
    print_info "Creating environment from environment.yml..."
    conda env create -f environment.yml
    
    print_info "✓ Combined environment created successfully!"
    print_info "Activate it with: conda activate umi_celegans_analysis"
}

# Function to create separate environments
create_separate_environments() {
    print_info "Creating separate environments for each tool..."
    
    # Array of tools with their versions and environment names
    declare -A tools=(
        ["fastqc"]="0.12.1"
        ["bcftools"]="1.22"
        ["umi_tools"]="1.1.6"
        ["gatk4"]="4.6.1.0"
        ["fgbio"]="2.5.0"
        ["snpeff"]="5.2"
        ["samtools.v1.22"]="1.22.1"
        ["STAR"]="2.7.11b"
        ["multiqc"]="1.32"
    )
    
    # Note: samtools.v1.22 is the environment name used in the pipeline scripts
    # STAR environment name matches the one used in the scripts
    
    for env_name in "${!tools[@]}"; do
        version="${tools[$env_name]}"
        # Get the package name (handle samtools.v1.22 -> samtools)
        if [[ "$env_name" == "samtools.v1.22" ]]; then
            package_name="samtools"
        else
            package_name="$env_name"
        fi
        
        print_info "Creating environment '$env_name' with $package_name=$version..."
        
        if conda env list | grep -q "^$env_name "; then
            print_warning "Environment '$env_name' already exists, skipping..."
            continue
        fi
        
        conda create -n "$env_name" -c bioconda -c conda-forge -c defaults \
            "$package_name=$version" -y
        
        print_info "✓ Environment '$env_name' created"
    done
    
    print_info "✓ All separate environments created successfully!"
    print_info "Environment names:"
    for env_name in "${!tools[@]}"; do
        echo "  - $env_name"
    done
}

# Main script logic
main() {
    print_info "UMI C. elegans Analysis - Conda Environment Setup"
    print_info "=================================================="
    
    check_conda
    
    # Parse command line argument
    option="${1:-all}"
    
    case "$option" in
        all)
            create_combined_environment
            ;;
        separate)
            create_separate_environments
            ;;
        help|--help|-h)
            show_help
            exit 0
            ;;
        *)
            print_error "Unknown option: $option"
            show_help
            exit 1
            ;;
    esac
    
    print_info "=================================================="
    print_info "Setup complete!"
}

# Run main function
main "$@"

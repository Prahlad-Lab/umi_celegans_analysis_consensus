#!/bin/bash

# Script to rename ENA downloaded files to descriptive sample names
# This script renames the ERR accession files to more meaningful sample names
# that indicate the strain, time point, condition, and replicate number.

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "ENA File Renaming Script"
echo "=========================================="
echo ""

# Define the mapping of ERR accessions to descriptive names
# Format: "ERR_ID:DESCRIPTIVE_NAME"
declare -a FILE_MAPPINGS=(
    "ERR15764772:N2.30min.HS.1"
    "ERR15764773:N2.30min.HS.2"
    "ERR15764774:N2.30min.HS.3"
    "ERR15764775:PRDE1.30min.HS.1"
    "ERR15764776:PRDE1.30min.HS.2"
    "ERR15764777:PRDE1.30min.HS.3"
)

# Function to rename a file
rename_file() {
    local source_file="$1"
    local target_file="$2"
    
    if [ ! -f "$source_file" ]; then
        echo -e "${YELLOW}  ⚠ Source file not found: $source_file${NC}"
        return 1
    fi
    
    if [ -f "$target_file" ]; then
        echo -e "${YELLOW}  ⚠ Target file already exists: $target_file (skipping)${NC}"
        return 1
    fi
    
    mv "$source_file" "$target_file"
    echo -e "${GREEN}  ✓ Renamed: $source_file → $target_file${NC}"
    return 0
}

# Counter for tracking
total_files=0
renamed_files=0
skipped_files=0
missing_files=0

# Process each mapping
for mapping in "${FILE_MAPPINGS[@]}"; do
    # Split the mapping into ERR ID and descriptive name
    err_id="${mapping%%:*}"
    descriptive_name="${mapping##*:}"
    
    echo "Processing: $err_id → $descriptive_name"
    
    # Process both read 1 and read 2
    for read_num in 1 2; do
        total_files=$((total_files + 1))
        
        source_file="${err_id}_${read_num}.fastq.gz"
        target_file="${descriptive_name}_R${read_num}.fastq.gz"
        
        if rename_file "$source_file" "$target_file"; then
            renamed_files=$((renamed_files + 1))
        else
            if [ ! -f "$source_file" ]; then
                missing_files=$((missing_files + 1))
            else
                skipped_files=$((skipped_files + 1))
            fi
        fi
    done
    echo ""
done

# Print summary
echo "=========================================="
echo "Summary:"
echo "=========================================="
echo "Total files processed: $total_files"
echo -e "${GREEN}Successfully renamed:  $renamed_files${NC}"
echo -e "${YELLOW}Skipped (existing):    $skipped_files${NC}"
echo -e "${RED}Missing source files:  $missing_files${NC}"
echo ""

if [ $renamed_files -eq 0 ] && [ $missing_files -gt 0 ]; then
    echo -e "${RED}❌ No files were renamed. Make sure you've downloaded the ENA files first.${NC}"
    echo "   Run: ./ena-file-download-read_run-PRJEB101605-fastq_ftp-20251208-1709.sh"
    exit 1
elif [ $renamed_files -gt 0 ]; then
    echo -e "${GREEN}✓ File renaming completed successfully!${NC}"
    exit 0
else
    echo -e "${YELLOW}⚠ All target files already exist. No renaming needed.${NC}"
    exit 0
fi

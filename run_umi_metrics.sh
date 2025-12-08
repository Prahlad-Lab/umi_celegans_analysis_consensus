#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="UMI_Metrics_Final"
#SBATCH --output=UMI_Metrics_Final.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --time=02:00:00

set -euo pipefail

# --- 1. Define Variables ---
BASE_OUTPUT_DIR="/vscratch/grp-vprahlad/umi_celegans_consensus"
DIR_UMI_MARKING="${BASE_OUTPUT_DIR}/UMIAwareDuplicateMarkingGenomeNoClip"
DIR_UMI_METRICS="${BASE_OUTPUT_DIR}/UMI_Metrics"

# Create output directory
mkdir -p ${DIR_UMI_METRICS}

# Define samples
sample_list=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")

# Activate Environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate umi_celegans_analysis

# --- 2. Define the Metrics Function ---
run_metrics_calc() {
    sample=$1
    # Input: Coordinate sorted BAM (Pre-consensus)
    input_bam="${DIR_UMI_MARKING}/${sample}.MergeBamAlignment.coordinate_sorted.bam"
    
    # Outputs
    hist_file="${DIR_UMI_METRICS}/${sample}.family_size_histogram.txt"
    summary_file="${DIR_UMI_METRICS}/${sample}.summary_stats.txt"

    # Check input
    if [ ! -f "$input_bam" ]; then
        echo "[ERROR] Input not found: $input_bam"
        return 1
    fi

    echo "[RUN] Generating Histogram for ${sample}..."
    
    # 1. Generate Family Size Histogram (Heavy lifting)
    # We send the BAM output to /dev/null because we only want the text report.
    if [ ! -s "${hist_file}" ]; then
        fgbio -Xmx8g GroupReadsByUmi \
            --input "${input_bam}" \
            --strategy Adjacency --edits 1 \
            --output /dev/null \
            --family-size-histogram "${hist_file}"
    fi

    # 2. Calculate Saturation & Means (Lightweight calculation)
    # We use awk to parse the histogram and calculate the metrics you need.
    # Logic: 
    #   Total Reads = Sum(FamilySize * Count)
    #   Total Molecules (Families) = Sum(Count)
    #   Saturation = Total Molecules / Total Reads (Higher is less saturated/more complex)
    #   Mean Family Size = Total Reads / Total Molecules
    
    echo "--- Metrics Summary for ${sample} ---" > "${summary_file}"
    awk '
    NR > 1 { 
        size = $1; 
        count = $2; 
        
        total_molecules += count; 
        total_reads += (size * count);
        
        if (size == 1) { singletons += count }
    }
    END {
        if (total_molecules > 0) {
            mean_fam = total_reads / total_molecules;
            saturation = total_molecules / total_reads;
            pct_singletons = (singletons / total_molecules) * 100;
            
            print "Total Reads Processed: " total_reads;
            print "Total UMI Families:    " total_molecules;
            print "Mean Family Size:      " mean_fam;
            print "Saturation (Unique/Total): " saturation;
            print "Singleton Families (%): " pct_singletons;
        } else {
            print "No data found.";
        }
    }' "${hist_file}" >> "${summary_file}"
    
    cat "${summary_file}"
}
export -f run_metrics_calc
export BASE_OUTPUT_DIR DIR_UMI_MARKING DIR_UMI_METRICS

# --- 3. Execute ---
echo "--- Starting Robust Metrics Generation ---"
parallel -j 6 run_metrics_calc {} ::: "${sample_list[@]}"
echo "--- Completed ---"
#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name="Consensus_UMI_Pipeline"
#SBATCH --output=Consensus_UMI_Pipeline.11.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --mem=100G
#SBATCH --time=48:00:00

#==============================================================================
# SAFETY & CONFIGURATION
#==============================================================================
set -euo pipefail
CLEANUP_INTERMEDIATES="false"

#==============================================================================
# DEFINE VARIABLES
#==============================================================================
START_TIME=$SECONDS
echo "--- Pipeline started at $(date) ---"

DIR_SEQS="/projects/rpci/vprahlad/umi_celegans/data/seqs"
DIR_REF_INPUT="/projects/rpci/vprahlad/umi_celegans/data/input"
REF_GENOME="${DIR_REF_INPUT}/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
REF_DICT="${DIR_REF_INPUT}/Caenorhabditis_elegans.WBcel235.dna.toplevel.dict"
REF_GTF_FILE="${DIR_REF_INPUT}/Caenorhabditis_elegans.WBcel235.114.gtf"
KNOWN_SITES_VCF="${DIR_REF_INPUT}/random_subset.vcf.gz"
STAR_INDEX="star_index"
PYTHON_SCRIPT_PATH="/vscratch/grp-vprahlad/umi_celegans_analysis/Separate_Scripts/calculate_allele_proportions_depth.py"

BASE_OUTPUT_DIR="/vscratch/grp-vprahlad/umi_celegans_consensus"
DIR_FASTQC="${BASE_OUTPUT_DIR}/fastqc"
DIR_FASTQ_TO_UBAM="${BASE_OUTPUT_DIR}/FastqToUbam"
DIR_EXTRACT_UMIS="${BASE_OUTPUT_DIR}/ExtractUmisFromBam"
DIR_STAR_ALIGN="${BASE_OUTPUT_DIR}/STARNoClip"
DIR_UMI_MARKING="${BASE_OUTPUT_DIR}/UMIAwareDuplicateMarkingGenomeNoClip"
DIR_CONSENSUS="${BASE_OUTPUT_DIR}/Consensus_Analysis"
DIR_FILTER_BAM="${BASE_OUTPUT_DIR}/FilterBambyFamilySizeGenome" # Kept for stats only
DIR_HAPLOTYPE_CALLER="${BASE_OUTPUT_DIR}/HaplotypeCaller"
DIR_ANNOTATION="${BASE_OUTPUT_DIR}/Annotation"
DIR_BCFTOOLS_PILEUP="${BASE_OUTPUT_DIR}/Bcftools_Mpileup"
DIR_ALLELE_PROPORTIONS="${BASE_OUTPUT_DIR}/Allele_Proportions"
DIR_MULTIQC="${BASE_OUTPUT_DIR}/MultiQC"

sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")
# 1. EXPORT VARIABLES so 'parallel' can see them
export REF_GENOME REF_DICT REF_GTF_FILE KNOWN_SITES_VCF STAR_INDEX PYTHON_SCRIPT_PATH
export BASE_OUTPUT_DIR DIR_FASTQC DIR_FASTQ_TO_UBAM DIR_EXTRACT_UMIS DIR_STAR_ALIGN \
       DIR_UMI_MARKING DIR_CONSENSUS DIR_FILTER_BAM DIR_HAPLOTYPE_CALLER \
       DIR_ANNOTATION DIR_BCFTOOLS_PILEUP DIR_ALLELE_PROPORTIONS DIR_MULTIQC
	   
mkdir -p ${DIR_FASTQC} ${DIR_FASTQ_TO_UBAM} ${DIR_EXTRACT_UMIS} ${DIR_STAR_ALIGN} \
         ${DIR_UMI_MARKING} ${DIR_CONSENSUS} ${DIR_FILTER_BAM} ${DIR_HAPLOTYPE_CALLER} \
         ${DIR_ANNOTATION} ${DIR_BCFTOOLS_PILEUP} ${DIR_ALLELE_PROPORTIONS} ${DIR_MULTIQC}
mkdir -p stdout

# Activate Environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate umi_celegans_analysis

#==============================================================================
# FUNCTIONS
#==============================================================================
run_fastqc() {
    sample=$1; DIR_IN=$2; DIR_OUT=$3
    
    # Check and run for Read 1
    if [ ! -s "${DIR_OUT}/${sample}_R1_fastqc.html" ]; then
        echo "[RUN] FastQC ${sample} R1"
        fastqc -q -o "${DIR_OUT}" "${DIR_IN}/${sample}_R1.fastq.gz"
    else
        echo "[SKIP] FastQC ${sample} R1"
    fi

    # Check and run for Read 2
    if [ ! -s "${DIR_OUT}/${sample}_R2_fastqc.html" ]; then
        echo "[RUN] FastQC ${sample} R2"
        fastqc -q -o "${DIR_OUT}" "${DIR_IN}/${sample}_R2.fastq.gz"
    else
        echo "[SKIP] FastQC ${sample} R2"
    fi
}
export -f run_fastqc
# --- Step 2: FastqToUbam ---
run_fastq_to_ubam() {
    sample=$1; DIR_IN=$2; DIR_OUT=$3
    output_file="${DIR_OUT}/${sample}.unmapped.bam"
    if [ -s "${output_file}" ]; then echo "[SKIP] FastqToUbam ${sample}"; return 0; fi
    
    echo "[RUN] FastqToUbam ${sample}"
    gatk FastqToSam SORT_ORDER=unsorted F1=${DIR_IN}/${sample}_R1.fastq.gz F2=${DIR_IN}/${sample}_R2.fastq.gz \
      SM=${sample} LB="L1" PL="ILLUMINA" PU="barcode1" RG="RG1" CN="BI" O=${output_file} ALLOW_EMPTY_FASTQ=True
}
export -f run_fastq_to_ubam

# --- Step 3: ExtractUMIs ---
run_extract_umis() {
    sample=$1; DIR_IN=$2; DIR_OUT=$3
    output_file="${DIR_OUT}/${sample}.extractUMIs.out.bam"
    if [ -s "${output_file}" ]; then echo "[SKIP] ExtractUMIs ${sample}"; return 0; fi

    echo "[RUN] ExtractUMIs ${sample}"
    fgbio ExtractUmisFromBam --input ${DIR_IN}/${sample}.unmapped.bam --read-structure +T --read-structure 8M6S+T --molecular-index-tags RX --output ${output_file}
}
export -f run_extract_umis

# --- Step 4: Sort UMI BAM ---
run_sort_umi_bam() {
    sample=$1; DIR_IN=$2
    output_file="${DIR_IN}/${sample}.querysort.extractUMIs.out.bam"
    if [ -s "${output_file}" ]; then echo "[SKIP] Sort UMI BAM ${sample}"; return 0; fi

    echo "[RUN] Sort UMI BAM ${sample}"
    gatk SortSam INPUT=${DIR_IN}/${sample}.extractUMIs.out.bam OUTPUT=${output_file} SORT_ORDER="queryname" MAX_RECORDS_IN_RAM=300000
}
export -f run_sort_umi_bam

# --- Step 5: STAR Alignment ---
run_star_align() {
    sample=$1; DIR_IN=$2; DIR_OUT=$3; IDX=$4
    output_file="${DIR_OUT}/${sample}.Aligned.out.bam"
    if [ -s "${output_file}" ]; then echo "[SKIP] STAR Alignment ${sample}"; return 0; fi

    echo "[RUN] STAR Alignment ${sample}"
    STAR --runThreadN 8 --genomeDir ${IDX} --readFilesIn ${DIR_IN}/${sample}.querysort.extractUMIs.out.bam \
         --outSAMtype BAM Unsorted --outFileNamePrefix ${DIR_OUT}/${sample}. --readFilesType SAM PE \
         --readFilesCommand samtools view -h --outSAMunmapped Within --twopassMode Basic
}
export -f run_star_align

# --- Step 6: Sort Aligned BAM ---
run_sort_aligned() {
    sample=$1; DIR_IN=$2; DIR_OUT=$3
    output_file="${DIR_OUT}/${sample}.querysort.Aligned.out.bam"
    if [ -s "${output_file}" ]; then echo "[SKIP] Sort Aligned BAM ${sample}"; return 0; fi

    echo "[RUN] Sort Aligned BAM ${sample}"
    gatk SortSam INPUT=${DIR_IN}/${sample}.Aligned.out.bam OUTPUT=${output_file} SORT_ORDER="queryname" MAX_RECORDS_IN_RAM=300000
}
export -f run_sort_aligned

# --- Step 7: MergeBamAlignment ---
run_merge_bam() {
    sample=$1; DIR_UNMAPPED=$2; DIR_ALIGNED=$3; DIR_OUT=$4; REF=$5
    output_file="${DIR_OUT}/${sample}.MergeBamAlignment.bam"
    if [ -s "${output_file}" ]; then echo "[SKIP] MergeBamAlignment ${sample}"; return 0; fi

    echo "[RUN] MergeBamAlignment ${sample}"
    gatk MergeBamAlignment --REFERENCE_SEQUENCE ${REF} --UNMAPPED_BAM ${DIR_UNMAPPED}/${sample}.querysort.extractUMIs.out.bam \
            --ALIGNED_BAM ${DIR_ALIGNED}/${sample}.querysort.Aligned.out.bam --OUTPUT ${output_file} \
            --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT
}
export -f run_merge_bam

# --- Step 8: Sort Merged BAM ---
run_sort_merged() {
    sample=$1; DIR_IN=$2
    output_file="${DIR_IN}/${sample}.MergeBamAlignment.coordinate_sorted.bam"
    if [ -s "${output_file}" ]; then echo "[SKIP] Sort Merged BAM ${sample}"; return 0; fi

    echo "[RUN] Sort Merged BAM ${sample}"
    gatk SortSam INPUT=${DIR_IN}/${sample}.MergeBamAlignment.bam OUTPUT=${output_file} SORT_ORDER="coordinate" CREATE_INDEX=true MAX_RECORDS_IN_RAM=300000
}
export -f run_sort_merged



# --- Step 9.5: Consensus Analysis (Corrected & Robust) ---
consensus_analysis() {
    sample=$1; DIR_UMI_MARKING=$2; DIR_CONSENSUS=$3; STAR_INDEX=$4; REF_GENOME=$5
    final_output="${DIR_CONSENSUS}/${sample}.consensus.merged.bam"
    
    # Global check: If the final file exists, skip everything
    if [ -s "${final_output}" ]; then echo "[SKIP] Consensus Gen ${sample}"; return 0; fi

    echo "[RUN] Generating Consensus for ${sample}"
    
    # 1. Group (fgbio)
    if [ ! -s "${DIR_CONSENSUS}/${sample}.fgbio_grouped.bam" ]; then
        fgbio -Xmx8g GroupReadsByUmi --input ${DIR_UMI_MARKING}/${sample}.MergeBamAlignment.coordinate_sorted.bam \
            --strategy Adjacency --edits 1 --output ${DIR_CONSENSUS}/${sample}.fgbio_grouped.bam
    else
        echo "[SKIP] Internal Step 1: GroupReadsByUmi (Found existing file)"
    fi
        
    # 2. Call Consensus
    if [ ! -s "${DIR_CONSENSUS}/${sample}.consensus.unmapped.bam" ]; then
        fgbio CallMolecularConsensusReads --input ${DIR_CONSENSUS}/${sample}.fgbio_grouped.bam \
            --output ${DIR_CONSENSUS}/${sample}.consensus.unmapped.bam --min-reads 1 --error-rate-post-umi 40 --error-rate-pre-umi 45
    else
        echo "[SKIP] Internal Step 2: CallMolecularConsensusReads (Found existing file)"
    fi

    # 2.5 Sort Unmapped Consensus BAM (NEW REQUIREMENT)
    # Creates: .consensus.unmapped.querysort.bam
    if [ ! -s "${DIR_CONSENSUS}/${sample}.consensus.unmapped.querysort.bam" ]; then
        echo "[RUN] Internal Step 2.5: Sorting Unmapped Consensus BAM"
        gatk SortSam \
            INPUT=${DIR_CONSENSUS}/${sample}.consensus.unmapped.bam \
            OUTPUT=${DIR_CONSENSUS}/${sample}.consensus.unmapped.querysort.bam \
            SORT_ORDER="queryname" \
            MAX_RECORDS_IN_RAM=300000
    fi
        
    # 3. Align Consensus
    # We use the ORIGINAL unmapped file here for alignment (order doesn't matter for STAR)
    if [ ! -s "${DIR_CONSENSUS}/${sample}.consensus.Aligned.out.bam" ]; then
        STAR --runThreadN 4 --genomeDir ${STAR_INDEX} --readFilesIn ${DIR_CONSENSUS}/${sample}.consensus.unmapped.bam \
             --outSAMtype BAM Unsorted --outFileNamePrefix ${DIR_CONSENSUS}/${sample}.consensus. \
             --readFilesType SAM PE --readFilesCommand samtools view -h --outSAMunmapped Within
    else
         echo "[SKIP] Internal Step 3: STAR Alignment (Found existing file)"
    fi
         
    # 3.5 Sort Aligned BAM by QueryName
    if [ ! -s "${DIR_CONSENSUS}/${sample}.consensus.querysort.Aligned.out.bam" ]; then
        echo "[RUN] Internal Step 3.5: Sorting Aligned Consensus BAM"
        gatk SortSam \
            INPUT=${DIR_CONSENSUS}/${sample}.consensus.Aligned.out.bam \
            OUTPUT=${DIR_CONSENSUS}/${sample}.consensus.querysort.Aligned.out.bam \
            SORT_ORDER="queryname" \
            MAX_RECORDS_IN_RAM=300000
    fi

    # 4. Merge Consensus
    # INPUT 1: .consensus.unmapped.querysort.bam (Sorted Unmapped)
    # INPUT 2: .consensus.querysort.Aligned.out.bam (Sorted Aligned)
    gatk MergeBamAlignment --REFERENCE_SEQUENCE ${REF_GENOME} \
        --UNMAPPED_BAM ${DIR_CONSENSUS}/${sample}.consensus.unmapped.querysort.bam \
        --ALIGNED_BAM ${DIR_CONSENSUS}/${sample}.consensus.querysort.Aligned.out.bam \
        --OUTPUT ${final_output} \
        --INCLUDE_SECONDARY_ALIGNMENTS false --VALIDATION_STRINGENCY SILENT
        
    # 5. Metrics
    gatk CollectAlignmentSummaryMetrics -R ${REF_GENOME} -I ${final_output} -O ${DIR_CONSENSUS}/${sample}.CONSENSUS_alignment_metrics.txt
}
export -f consensus_analysis

# --- Step 12/13: Variant Calling (Updated for Consensus Input) ---
run_variant_calling() {
    sample=$1
    type=$2       # Label (e.g., "consensus")
    input_bam=$3  # NOW TAKES CONSENSUS BAM
    D_HC=$4; D_ANN=$5; REF=$6; VCF_K=$7; REF_DICT=$8

    prefix="${sample}.${type}"
    final_output="${D_ANN}/${prefix}.ann.bed"
    if [ -s "${final_output}" ]; then echo "[SKIP] Variant Calling ${type} ${sample}"; return 0; fi

    echo "[RUN] Variant Calling (${type}) for ${sample}"
    
    # GATK Pipeline (SplitN -> BQSR -> HaplotypeCaller)
    gatk SplitNCigarReads -R ${REF} -I ${input_bam} -O ${D_HC}/${prefix}.split.bam
    gatk BaseRecalibrator -R ${REF} -I ${D_HC}/${prefix}.split.bam --use-original-qualities -O ${D_HC}/${prefix}.recal.table --known-sites ${VCF_K}
    gatk ApplyBQSR -R ${REF} -I ${D_HC}/${prefix}.split.bam -O ${D_HC}/${prefix}.bqsr.bam --bqsr-recal-file ${D_HC}/${prefix}.recal.table --sequence-dictionary ${REF_DICT}
    gatk HaplotypeCaller -R ${REF} -I ${D_HC}/${prefix}.bqsr.bam --dont-use-soft-clipped-bases true -O ${D_HC}/${prefix}.vcf.gz \
        --standard-min-confidence-threshold-for-calling 20 --dbsnp ${VCF_K}
        
    # Annotation
    snpEff WBcel235.105 -stats ${D_ANN}/${prefix}_summary.html ${D_HC}/${prefix}.vcf.gz > ${D_ANN}/${prefix}.ann.vcf
    snpEff ann WBcel235.105 -csvStats ${D_ANN}/${prefix}_summary.csv -o bedANN ${D_HC}/${prefix}.vcf.gz > ${final_output}
    
    rm -f ${D_HC}/${prefix}.split.bam ${D_HC}/${prefix}.split.bai
}
export -f run_variant_calling

# --- Step 14/15: Mpileup ---
run_mpileup_and_props() {
    sample=$1; type=$2; D_HC=$3; D_BCF=$4; D_PROP=$5; REF=$6; PY_SCRIPT=$7
    prefix="${sample}.${type}"
    final_output="${D_PROP}/${prefix}.summary_totals.txt"
    if [ -s "${final_output}" ]; then echo "[SKIP] Mpileup ${prefix}"; return 0; fi

    echo "[RUN] Mpileup for ${prefix}"
    # Note: Uses the BQSR BAM generated in the Variant Calling step
    bcftools mpileup -d 1000000 -f ${REF} ${D_HC}/${prefix}.bqsr.bam > ${D_BCF}/${prefix}.pileup
    python ${PY_SCRIPT} ${D_BCF}/${prefix}.pileup ${D_PROP}/${prefix}.per_position_results.tsv ${final_output} --min-depth 2
}
export -f run_mpileup_and_props

#==============================================================================
# EXECUTION
#==============================================================================

# --- Pre-Processing (Steps 1-9) ---
echo "--- Step 1: FastQC ---"
# Runs 6 samples in parallel. Each FastQC job is single-threaded (default), 
# which is efficient when running 6 at once.
parallel -j 6 run_fastqc {} "${DIR_SEQS}" "${DIR_FASTQC}" ::: "${sample[@]}"

echo "--- Step 2: FastqToUbam ---"
parallel -j 6 run_fastq_to_ubam {} "${DIR_SEQS}" "${DIR_FASTQ_TO_UBAM}" ::: "${sample[@]}"

echo "--- Step 3: ExtractUMIs ---"
parallel -j 6 run_extract_umis {} "${DIR_FASTQ_TO_UBAM}" "${DIR_EXTRACT_UMIS}" ::: "${sample[@]}"
if [ "$CLEANUP_INTERMEDIATES" = "true" ]; then rm -f ${DIR_FASTQ_TO_UBAM}/*.unmapped.bam; fi

echo "--- Step 4: Sort UMI BAM ---"
parallel -j 6 run_sort_umi_bam {} "${DIR_EXTRACT_UMIS}" ::: "${sample[@]}"
if [ "$CLEANUP_INTERMEDIATES" = "true" ]; then rm -f ${DIR_EXTRACT_UMIS}/*.extractUMIs.out.bam; fi

echo "--- Step 4.5: STAR Index ---"
if [ ! -d "${STAR_INDEX}" ] || [ ! "$(ls -A $STAR_INDEX)" ]; then
    echo "Generating STAR index..."
    mkdir -p ${STAR_INDEX}
    STAR --runThreadN ${SLURM_CPUS_PER_TASK:-32} --runMode genomeGenerate --genomeDir ${STAR_INDEX} \
         --genomeFastaFiles ${REF_GENOME} --sjdbGTFfile ${REF_GTF_FILE} --sjdbOverhang 100
fi

echo "--- Step 5: STAR Alignment ---"
parallel -j 4 run_star_align {} "${DIR_EXTRACT_UMIS}" "${DIR_STAR_ALIGN}" "${STAR_INDEX}" ::: "${sample[@]}"

echo "--- Step 6: Sort Aligned BAM ---"
parallel -j 6 run_sort_aligned {} "${DIR_STAR_ALIGN}" "${DIR_UMI_MARKING}" ::: "${sample[@]}"
if [ "$CLEANUP_INTERMEDIATES" = "true" ]; then rm -f ${DIR_STAR_ALIGN}/*.Aligned.out.bam; fi

echo "--- Step 7: MergeBamAlignment ---"
parallel -j 6 run_merge_bam {} "${DIR_EXTRACT_UMIS}" "${DIR_UMI_MARKING}" "${DIR_UMI_MARKING}" "${REF_GENOME}" ::: "${sample[@]}"
if [ "$CLEANUP_INTERMEDIATES" = "true" ]; then rm -f ${DIR_UMI_MARKING}/*.querysort.Aligned.out.bam ${DIR_EXTRACT_UMIS}/*.querysort.extractUMIs.out.bam; fi

echo "--- Step 8: Sort Merged BAM ---"
parallel -j 6 run_sort_merged {} "${DIR_UMI_MARKING}" ::: "${sample[@]}"
if [ "$CLEANUP_INTERMEDIATES" = "true" ]; then rm -f ${DIR_UMI_MARKING}/*.MergeBamAlignment.bam; fi



# --- Consensus Generation (Step 9.5) ---
echo "--- Step 9.5: Consensus Analysis ---"
parallel -j 4 consensus_analysis {} "${DIR_UMI_MARKING}" "${DIR_CONSENSUS}" "${STAR_INDEX}" "${REF_GENOME}" ::: "${sample[@]}"


# --- Downstream Analysis on CONSENSUS BAMs ---
echo "--- Steps 12 & 13: Variant Calling (Consensus) ---"
# NOTE: Changed input from filtered BAM to CONSENSUS BAM
parallel -j 6 run_variant_calling {1} "consensus" "${DIR_CONSENSUS}/{1}.consensus.merged.bam" \
    "${DIR_HAPLOTYPE_CALLER}" "${DIR_ANNOTATION}" "${REF_GENOME}" "${KNOWN_SITES_VCF}" "${REF_DICT}" ::: "${sample[@]}"

echo "--- Steps 14 & 15: Allele Proportions (Consensus) ---"
parallel -j 6 run_mpileup_and_props {1} "consensus" "${DIR_HAPLOTYPE_CALLER}" "${DIR_BCFTOOLS_PILEUP}" "${DIR_ALLELE_PROPORTIONS}" "${REF_GENOME}" "${PYTHON_SCRIPT_PATH}" ::: "${sample[@]}"


echo "--- Step 16: MultiQC ---"
if command -v multiqc &> /dev/null; then multiqc ${BASE_OUTPUT_DIR} -o ${DIR_MULTIQC}; fi

echo "--- Pipeline Completed ---"

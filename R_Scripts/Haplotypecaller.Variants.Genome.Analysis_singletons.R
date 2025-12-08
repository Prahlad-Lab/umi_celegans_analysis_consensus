# ---
# title: "Analysis Variants Genome"
# output:
#   html_notebook:
#     toc: TRUE
#     code_folding: hide
#     theme: flatly
# editor: visual
# ---
# 
# # Goal of the Project
# 
# The goal of the project was to assess whether a particular mutation in the nematode C. elegans caused more transcriptional errors. For this we used UMIs to tag mRNA that was extracted from 3 biological repeats of: 1. Control (wild type) C. elegans, called N2 and 2. Mutant C.elegans, called PRDE1
# 
# mRNA were tagged with UMIs using the SMART-Seq® Total RNA Pico Input with UMIs (ZapR™ Mammalian) at the GSR by Prashant and colleagues. RNA was tagged with the molecular barcodes during the first strand-cDNA synthesis. PCR amplification was conducted with 3 cycles of amplification for the initial library PCR and 12 cycles of amplification for the enrichment PCR. Sequencing was done on the Novaseq 6000; the number paired-end reads sequenced were in the range of 50 million-80million, with a read size of 100 bp.
# 
# # Analysis Steps
# 
# 1.  Fastqc was used to check the quality of the reads.
# 
# 2.  FGBIO Fastq to ubam was used to convert the fastq files to unmapped bam files and extract the umi-reads.
# 
# 3.  Align the reads with STAR. Two pass. Using Bulk RNA-seq Encode settings. BAM file aligned to the genome.
# 
# 4.  Intermediary steps Sort bam files first by coordinate then by query names.
# 
# 5.  Groups UMI with Umit-tools group
# 
# 6.  Markduplicates (Umi Aware) with GATK MarkDuplicates
# 
# 7.  Quantify BAM files and perform quality control RNAseqc
# 
# 8.  Call Variants with HaplotypeCaller
# 
# 9.  Annotate Variants with Snepeff
# 
# 10. Import variants to R to get differential variant analysis
# 
# ```{r,warning=F,message=FALSE}
rm(list = ls())
# ```
# 
# ```{r,warning=F,message=FALSE}
# 

#loadfonts(device = "win")
# plotting theme

# 
# ```
# 
# # VCF File Analysis and Comparison Script GLOBAL
# 
# ```{r}
# VCF File Analysis and Comparison Script
#
# This script reads VCF files for two samples (N2 and PRDE1), each with three replicates.
# It extracts variant information, including Allele Depth (AD), Read Depth (DP),
# and SNPeff annotations (IMPACT).
# It then filters variants based on DP and IMPACT, and performs several statistical
# comparisons between the two samples.
#
# Version 43: Changed SnpEff rate calculations to normalize by total read depth.

# --- 1. SETUP: Install and Load Required Packages ---
# --- 1. SETUP: Install and Load Required Packages ---

# Define required packages from CRAN
cran_packages <- c(
  "extrafont", 
  "cowplot", 
  "ggplot2", 
  "plotly", 
  "ggpubr", 
  "ggprism",
  "vcfR", 
  "tidyverse", 
  "broom", 
  "rstatix", 
  "knitr", 
  "kableExtra", 
  "officer", 
  "flextable"
)

# Define required packages from Bioconductor
bioc_packages <- c(
  "rtracklayer",  # For importing GTF files
  "plyranges",    # For genomic range manipulation
  "GenomicRanges" # Used explicitly with ::width
)

# --- Install CRAN Packages ---
# Check for and install any missing CRAN packages
missing_cran <- cran_packages[!cran_packages %in% installed.packages()[,"Package"]]
if (length(missing_cran) > 0) {
  cat("Installing missing CRAN packages:", paste(missing_cran, collapse=", "), "\n")
  # *** FIXED: Corrected the repos string (no markdown) ***
  install.packages(missing_cran, repos = "[https://cloud.r-project.org/](https://cloud.r-project.org/)")
}

# --- Install Bioconductor Packages ---
# First, ensure BiocManager itself is installed
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  # *** FIXED: Corrected the repos string (no markdown) ***
  install.packages("BiocManager", repos = "[https://cloud.r-project.org/](https://cloud.r-project.org/)")
}

# Check for and install any missing Bioconductor packages
missing_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[,"Package"]]
if (length(missing_bioc) > 0) {
  cat("Installing missing Bioconductor packages:", paste(missing_bioc, collapse=", "), "\n")
  BiocManager::install(missing_bioc)
}

# --- Load All Libraries ---
cat("Loading all required libraries...\n")
suppressPackageStartupMessages({
  library(extrafont)
  library(cowplot)
  library(ggplot2)
  library(plotly)
  library(ggpubr)
  library(ggprism)
  library(vcfR)
  library(tidyverse)
  library(broom)
  library(rstatix)
  library(knitr)
  library(kableExtra)
  library(officer)
  library(flextable)
  library(rtracklayer)
  library(plyranges)
  library(GenomicRanges)
})

cat("All packages installed and loaded successfully.\n")
old <- theme_set(theme_cowplot(font_size = 12,font_family = "Arial"))
# --- End of Setup ---


# --- 2. DATA LOADING: Define File Paths ---

# *** REMOVED problematic comment line that was causing a parsing error ***
n2_files <- c("../Family_size1_annotation/N2.30min.HS.1.singletons.ann.vcf", "../Family_size1_annotation/N2.30min.HS.2.singletons.ann.vcf", "../Family_size1_annotation/N2.30min.HS.3.singletons.ann.vcf")
prde1_files <-  c("../Family_size1_annotation/PRDE1.30min.HS.1.singletons.ann.vcf", "../Family_size1_annotation/PRDE1.30min.HS.2.singletons.ann.vcf", "../Family_size1_annotation/PRDE1.30min.HS.3.singletons.ann.vcf")

# Check if files exist to prevent errors
all_files <- c(n2_files, prde1_files)
if (!all(file.exists(all_files))) {
  stop("One or more VCF files were not found. Please check your file paths.")
}


# --- 3. PROCESSING FUNCTION: Extract Data from a VCF File ---

#' Process a single VCF file to extract relevant variant data.
process_vcf <- function(file_path, sample_name, replicate_num) {
  vcf <- read.vcfR(file_path, verbose = FALSE)
  tidy_vcf <- vcfR2tidy(vcf,
                        # Add LOF and NMD to the info fields to be extracted
                        info_fields = c("ANN", "LOF", "NMD"),
                        format_fields = c("AD", "DP"),
                        single_frame = TRUE)
  
  processed_data <- tidy_vcf$dat %>%
    as_tibble() %>%
    tidyr::separate_rows(ALT, sep = ",") %>%
    dplyr::group_by(CHROM, POS, REF, Indiv) %>%
    dplyr::mutate(allele_num = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      ad_list = str_split(gt_AD, ","),
      AD_Ref = as.numeric(map_chr(ad_list, 1)),
      AD_Alt = as.numeric(map2_chr(ad_list, allele_num, ~ .x[.y + 1])),
      ann_list = str_split(ANN, ",")
    ) %>%
    mutate(
      ANN_allele = map2_chr(ann_list, allele_num, ~ .x[.y]),
      Sample = sample_name,
      Replicate = replicate_num
    ) %>%
    tidyr::separate_rows(ANN_allele, sep = "&") %>%
    dplyr::mutate(
      # Safely extract SnpEff fields, returning NA if a field is missing
      Effect_Type = map_chr(str_split(ANN_allele, "\\|"), ~ if(length(.x) >= 2) .x[[2]] else NA_character_),
      IMPACT = map_chr(str_split(ANN_allele, "\\|"), ~ if(length(.x) >= 3) .x[[3]] else NA_character_),
      GENE = map_chr(str_split(ANN_allele, "\\|"), ~ if(length(.x) >= 4) .x[[4]] else NA_character_),
      Functional_Class = case_when(
        str_detect(Effect_Type, "missense_variant") ~ "MISSENSE",
        str_detect(Effect_Type, "synonymous_variant") ~ "SILENT",
        str_detect(Effect_Type, "stop_gained") ~ "NONSENSE",
        TRUE ~ "OTHER"
      ),
      # The LOF and NMD columns are now booleans based on the presence of the INFO tag
      LoF = ifelse(!is.na(LOF),TRUE,FALSE) # This assumes LOF is a boolean flag (TRUE/FALSE) from vcfR2tidy
      ,NMD = if_else(!is.na(NMD),TRUE,FALSE)  # This assumes NMD is a boolean flag (TRUE/FALSE) from vcfR2tidy
    ) %>%
    dplyr::mutate(
      AD_Ref = ifelse(is.na(AD_Ref), 0, AD_Ref),
      AD_Alt = ifelse(is.na(AD_Alt), 0, AD_Alt),
      gt_DP = as.numeric(gt_DP)
    ) %>%
    dplyr::select(Sample, Replicate, CHROM, POS, REF, ALT, GENE, IMPACT, Effect_Type, Functional_Class, LoF, NMD, AD_Ref, AD_Alt, gt_DP) %>%
    dplyr::rename(Read_Depth = gt_DP) %>%
    dplyr::mutate(Variant_Type = case_when(
      nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNP",
      nchar(REF) > nchar(ALT) ~ "Deletion",
      nchar(REF) < nchar(ALT) ~ "Insertion",
      TRUE ~ "Other"
    ))
  return(processed_data)
}


# --- 4. COMBINE DATA: Process All VCF Files ---

file_metadata <- tibble(
  path = c(n2_files, prde1_files),
  sample = rep(c("N2", "PRDE1"), each = 3),
  replicate = rep(1:3, times = 2)
)

all_variants_raw <- pmap_dfr(file_metadata, ~process_vcf(..1, ..2, ..3))
cat("Successfully processed", nrow(all_variants_raw), "total variants from all files.\n")

all_variants_raw =all_variants_raw |> dplyr::mutate(Variant_ratio = AD_Alt/Read_Depth)

# --- 5. FILTERING: Apply Filters for Impact and Read Depth ---

# Define filtering criteria
MIN_READ_DEPTH <- 10
DESIRED_IMPACTS <- c("HIGH", "MODERATE", "LOW")

# Apply filters sequentially
variants_filtered <- all_variants_raw %>% 
  dplyr::filter(Read_Depth >= MIN_READ_DEPTH) %>%
  dplyr::filter(IMPACT %in% DESIRED_IMPACTS)  |> 
  dplyr::filter(Variant_ratio < 1.0) 

cat("Filtered down to", nrow(variants_filtered), "total variant records matching criteria.\n")


# --- 6A. STATISTICAL ANALYSIS: Global Allele Proportions per Replicate ---
cat("\n--- Running Global Allele Proportion Analysis per Replicate ---\n")
# Calculate proportions for each replicate
replicate_proportions <- variants_filtered %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::summarise(
    Total_Alt = sum(AD_Alt),
    Total_Ref = sum(AD_Ref),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(
    Total_Reads = Total_Alt + Total_Ref,
    Proportion_Alt = Total_Alt / Total_Reads,
    Proportion_Ref = Total_Ref / Total_Reads
  )

# Perform Wilcoxon test on the alternate allele proportions
stat_test_alt_proportion <- replicate_proportions %>%
  rstatix::wilcox_test(Proportion_Alt ~ Sample, detailed = TRUE)
print(stat_test_alt_proportion)


# --- 6B. STATISTICAL ANALYSIS: Allele Read Proportions by Variant Type per Replicate ---
cat("\n--- Running Allele Read Proportion by Variant Type per Replicate ---\n")
# Calculate proportions of each variant type within the alternate alleles for each replicate
alt_allele_proportions_by_replicate <- variants_filtered %>%
  dplyr::filter(Variant_Type %in% c("SNP", "Insertion", "Deletion")) %>%
  dplyr::group_by(Sample, Replicate, Variant_Type) %>%
  dplyr::summarise(Total_Alt_Reads = sum(AD_Alt), .groups = 'drop') %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::mutate(Proportion = Total_Alt_Reads / sum(Total_Alt_Reads)) %>%
  dplyr::ungroup()

# Perform Wilcoxon test for each variant type's proportion
stat_test_alt_type_proportions <- alt_allele_proportions_by_replicate %>%
  group_by(Variant_Type) %>%
  rstatix::wilcox_test(Proportion ~ Sample, detailed = TRUE)
print(stat_test_alt_type_proportions)


# --- 6D. STATISTICAL ANALYSIS: Global Error Rate ---
cat("\n--- Running Wilcoxon Test on Global Error Rate ---\n")
error_rate_data <- variants_filtered %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::summarise(
    Total_Alt = sum(AD_Alt),
    Total_Ref = sum(AD_Ref),
    Error_Rate = Total_Alt / (Total_Alt + Total_Ref),
    .groups = 'drop'
  )
stat_test_error_rate <- error_rate_data %>%
  rstatix::wilcox_test(Error_Rate ~ Sample, detailed = TRUE)
print(stat_test_error_rate)


# --- Common Denominator for Rate Calculations ---
# Calculate the total read depth per replicate across all filtered variants
# This will be used as the denominator for all subsequent rate calculations
total_depth_per_replicate <- variants_filtered %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::summarise(Total_Replicate_Depth = sum(Read_Depth), .groups = 'drop')


# --- 6E. STATISTICAL ANALYSIS: Mismatch Error Rate ---
cat("\n--- Running Wilcoxon Test on Mismatch Error Rate ---\n")
# Calculate the rate for each substitution type
mismatch_error_data <- variants_filtered %>%
  dplyr::filter(Variant_Type == "SNP") %>%
  dplyr::mutate(Mismatch_Type = paste0(REF, ">", ALT)) %>%
  dplyr::group_by(Sample, Replicate, Mismatch_Type) %>%
  dplyr::summarise(Total_Alt_Mismatch = sum(AD_Alt), .groups = 'drop') %>%
  dplyr::left_join(total_depth_per_replicate, by = c("Sample", "Replicate")) %>%
  dplyr::mutate(Mismatch_Error_Rate = Total_Alt_Mismatch / Total_Replicate_Depth)

# Perform Wilcoxon test only on mismatch types with enough data
stat_test_mismatch_error_rate <- mismatch_error_data %>%
  group_by(Mismatch_Type) %>%
  filter(sum(Mismatch_Error_Rate[Sample == "N2"]) > 0 & sum(Mismatch_Error_Rate[Sample == "PRDE1"]) > 0) %>%
  rstatix::wilcox_test(Mismatch_Error_Rate ~ Sample, detailed = TRUE)
print(stat_test_mismatch_error_rate)


# --- 6F. STATISTICAL ANALYSIS: SnpEff Functional Class Rate ---
cat("\n--- Running Wilcoxon Test on SnpEff Functional Class Rate ---\n")
# Calculate rate for each functional class
functional_class_rate_data <- variants_filtered %>%
  dplyr::group_by(Sample, Replicate, Functional_Class) %>%
  dplyr::summarise(Total_Alt_Class = sum(AD_Alt), .groups = 'drop') %>% # Sum allele depths
  dplyr::left_join(total_depth_per_replicate, by = c("Sample", "Replicate")) %>%
  dplyr::mutate(Class_Rate = Total_Alt_Class / Total_Replicate_Depth)

# Perform Wilcoxon test only on classes with enough data
stat_test_functional_class_rate <- functional_class_rate_data %>%
  group_by(Functional_Class) %>%
  filter(sum(Class_Rate[Sample == "N2"], na.rm = TRUE) > 0 & sum(Class_Rate[Sample == "PRDE1"], na.rm = TRUE) > 0) %>%
  rstatix::wilcox_test(Class_Rate ~ Sample, detailed = TRUE)
print(stat_test_functional_class_rate)


# --- 6G. STATISTICAL ANALYSIS: SnpEff Effect Type Rate ---
cat("\n--- Running Wilcoxon Test on SnpEff Effect Type Rate ---\n")
# Calculate rate for each effect type
effect_type_rate_data <- variants_filtered %>%
  dplyr::group_by(Sample, Replicate, Effect_Type) %>%
  dplyr::summarise(Total_Alt_Effect = sum(AD_Alt), .groups = 'drop') %>% # Sum allele depths
  dplyr::left_join(total_depth_per_replicate, by = c("Sample", "Replicate")) %>%
  dplyr::mutate(Effect_Rate = Total_Alt_Effect / Total_Replicate_Depth)

# Perform Wilcoxon test only on effect types with enough data
stat_test_effect_type_rate <- effect_type_rate_data %>%
  group_by(Effect_Type) %>%
  filter(sum(Effect_Rate[Sample == "N2"], na.rm = TRUE) > 0 & sum(Effect_Rate[Sample == "PRDE1"], na.rm = TRUE) > 0) %>%
  rstatix::wilcox_test(Effect_Rate ~ Sample, detailed = TRUE)
print(stat_test_effect_type_rate)


# --- 6H. STATISTICAL ANALYSIS: SnpEff NMD Rate ---
cat("\n--- Running Wilcoxon Test on SnpEff NMD Rate ---\n")
# Sum alternate allele depths for NMD variants per replicate
nmd_alt_counts_per_replicate <- variants_filtered %>%
  dplyr::filter(NMD == TRUE) %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::summarise(Total_Alt_NMD = sum(AD_Alt), .groups = 'drop')

# Calculate rate by joining with total depth
nmd_rate_data <- total_depth_per_replicate %>%
  dplyr::full_join(nmd_alt_counts_per_replicate, by = c("Sample", "Replicate")) %>%
  dplyr::mutate(Total_Alt_NMD = coalesce(Total_Alt_NMD, 0L)) %>%
  dplyr::mutate(NMD_Rate = Total_Alt_NMD / Total_Replicate_Depth)

# Perform Wilcoxon test on NMD rates
stat_test_nmd_rate <- nmd_rate_data %>%
  rstatix::wilcox_test(NMD_Rate ~ Sample, detailed = TRUE)
print(stat_test_nmd_rate)


# --- 6I. STATISTICAL ANALYSIS: SnpEff LoF Rate ---
cat("\n--- Running Wilcoxon Test on SnpEff LoF Rate ---\n")
# Sum alternate allele depths for LoF variants per replicate
lof_alt_counts_per_replicate <- variants_filtered %>%
  dplyr::filter(LoF == TRUE) %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::summarise(Total_Alt_LoF = sum(AD_Alt), .groups = 'drop')

# Calculate rate by joining with total depth
lof_rate_data <- total_depth_per_replicate %>%
  dplyr::full_join(lof_alt_counts_per_replicate, by = c("Sample", "Replicate")) %>%
  dplyr::mutate(Total_Alt_LoF = coalesce(Total_Alt_LoF, 0L)) %>%
  dplyr::mutate(LoF_Rate = Total_Alt_LoF / Total_Replicate_Depth)

# Perform Wilcoxon test on LoF rates
stat_test_lof_rate <- lof_rate_data %>%
  rstatix::wilcox_test(LoF_Rate ~ Sample, detailed = TRUE)
print(stat_test_lof_rate)

# --- 8. ANALYSIS BY GENE LENGTH QUANTILE (Setup) ---
cat("\n--- Setting up Analysis by Gene Length Quantile ---\n")

# 8A. Load GTF and Define Gene Quantiles
gtf <- rtracklayer::import(con = "../input/Caenorhabditis_elegans.WBcel235.114.gtf")

protein_coding_genes <- gtf |> plyranges::filter(gene_biotype == "protein_coding") |> plyranges::filter(type == "gene")

# Define quantiles based on gene width (length)
quantiles <- quantile(GenomicRanges::width(protein_coding_genes), probs = c(0, 0.25, 0.5, 0.75, 1.0))

# Create new, descriptive labels with gene length ranges
q1_label <- paste0("Gene Length Q1: ", round(quantiles[1]), "-", round(quantiles[2]), "bp")
q2_label <- paste0("Gene Length Q2: ", round(quantiles[2] + 1), "-", round(quantiles[3]), "bp")
q3_label <- paste0("Gene Length Q3: ", round(quantiles[3] + 1), "-", round(quantiles[4]), "bp")
q4_label <- paste0("Gene Length Q4: ", round(quantiles[4] + 1), "-", round(quantiles[5]), "bp")
quantile_labels <- c(q1_label, q2_label, q3_label, q4_label)

genes_by_quantiles <- list(
  q1 = protein_coding_genes |> plyranges::filter(width <= quantiles[2]),
  q2 = protein_coding_genes |> plyranges::filter(width > quantiles[2] & width <= quantiles[3]),
  q3 = protein_coding_genes |> plyranges::filter(width > quantiles[3] & width <= quantiles[4]),
  q4 = protein_coding_genes |> plyranges::filter(width > quantiles[4])
)

# 8B. Filter Variants and Combine Quantile Dataframes
variants_filtered_quantile_1 <- variants_filtered |> dplyr::filter(GENE %in% genes_by_quantiles$q1$gene_name)
variants_filtered_quantile_2 <- variants_filtered |> dplyr::filter(GENE %in% genes_by_quantiles$q2$gene_name)
variants_filtered_quantile_3 <- variants_filtered |> dplyr::filter(GENE %in% genes_by_quantiles$q3$gene_name)
variants_filtered_quantile_4 <- variants_filtered |> dplyr::filter(GENE %in% genes_by_quantiles$q4$gene_name)

# Add a 'Quantile' identifier to each dataframe and combine them into one.
variants_by_quantile <- bind_rows(
  variants_filtered_quantile_1 %>% dplyr::mutate(Quantile = q1_label),
  variants_filtered_quantile_2 %>% dplyr::mutate(Quantile = q2_label),
  variants_filtered_quantile_3 %>% dplyr::mutate(Quantile = q3_label),
  variants_filtered_quantile_4 %>% dplyr::mutate(Quantile = q4_label)
) %>%
  # Make Quantile a factor to ensure correct plotting order
  dplyr::mutate(Quantile = factor(Quantile, levels = quantile_labels))


# --- 7. VISUALIZATIONS ---
cat("\n--- Generating Plots ---\n")

# Define the hardcoded n-value label
n_value_label <- tibble(label = "n = 2045-2054")
n_value_label_string <- "n = 2045-2054"


# 7A. Global Allele Proportions
# Summarize data for bar plot with SE bars
replicate_proportion_summary <- replicate_proportions %>%
  tidyr::pivot_longer(
    cols = c(Proportion_Alt, Proportion_Ref),
    names_to = "Allele_Type_Full",
    values_to = "Proportion"
  ) %>%
  dplyr::mutate(Allele_Type = ifelse(Allele_Type_Full == "Proportion_Ref", "Total_Ref", "Total_Alt")) %>%
  dplyr::group_by(Sample, Allele_Type) %>%
  dplyr::summarise(
    mean_prop = mean(Proportion),
    se_prop = sd(Proportion) / sqrt(dplyr::n()),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(Allele_Type = factor(Allele_Type, levels = c("Total_Ref", "Total_Alt")))

proportion_plot <- ggplot(replicate_proportion_summary, aes(x = Sample, y = mean_prop, fill = Allele_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
    width = 0.2,
    position = position_dodge(0.8)
  ) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(
    values = c("Total_Ref" = "#0000FF", "Total_Alt" = "#FF0000"),
    labels = c("Reference Alleles", "Alternate Alleles")
  ) +
  labs(
    subtitle = paste("Wilcoxon p-value (Alt Alleles):", format.pval(stat_test_alt_proportion$p, digits = 3)),
    x = "Sample",
    y = "Mean Proportion of Total Reads",
    fill = "Allele Type"
  ) +
  theme_cowplot(font_size = 12, font_family = "Arial")

print(ggplotly(proportion_plot))

# 7B. Allele Read Proportions by Variant Type
# Summarize data for bar plot with SE bars
alt_allele_proportion_summary <- alt_allele_proportions_by_replicate %>%
  dplyr::group_by(Sample, Variant_Type) %>%
  dplyr::summarise(
    mean_prop = mean(Proportion),
    se_prop = sd(Proportion) / sqrt(dplyr::n()),
    .groups = 'drop'
  )

allele_read_plot <- ggplot(alt_allele_proportion_summary, aes(x = Sample, y = mean_prop, fill = Variant_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
    width = 0.2,
    position = position_dodge(0.8)
  ) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(
    values = c("SNP" = "#00A087FF", "Deletion" = "#f39b7f", "Insertion" = "#4dbbd5"),
    name = "Alternate Allele Type"
  ) +
  labs(
    subtitle = "Comparing proportions of alternate allele types",
    x = "Sample",
    y = "Mean Proportion of Alternate Reads"
  ) +
  theme_cowplot(font_size = 12, font_family = "Arial")

print(ggplotly(allele_read_plot))

# 7D. Global Error Rate
error_rate_plot <- ggplot(error_rate_data, aes(x = Sample, y = Error_Rate, fill = Sample)) +
  geom_boxplot(staplewidth = 0.5, outlier.shape = NA, width = 0.3) +
  geom_jitter(width = 0.3, height = 0) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.15, 0.15))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(
    subtitle = paste("Wilcoxon Test p-value:", format.pval(stat_test_error_rate$p, digits = 3)),
    x = "Sample",
    y = "Percentage \nNon-Reference Reads/Total Reads"
  ) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(legend.position = "none")

print(ggplotly(error_rate_plot))

# 7E. Mismatch Error Rate (Consolidated Bar Plot)
mismatch_summary <- mismatch_error_data %>%
  dplyr::group_by(Sample, Mismatch_Type) %>%
  dplyr::summarise(
    mean_rate = mean(Mismatch_Error_Rate),
    se_rate = sd(Mismatch_Error_Rate) / sqrt(dplyr::n()),
    .groups = 'drop'
  )

mismatch_plots = mismatch_summary |> 
  dplyr::mutate(Mismatch_Type = fct_relevel(Mismatch_Type,c("A>G", "G>A", "G>C" ,"T>C", "A>C", "A>T", "C>A" ,"C>G" ,"C>T", "G>T", "T>A" ,"T>G"))) |> 
  ggplot( aes(x = Mismatch_Type, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(
    x = "Mismatch Type",
    y = "Mean Mismatch Error Rate"
  ) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))


# 7F. SnpEff Functional Class Rate Plot
functional_class_summary <- functional_class_rate_data %>%
  dplyr::group_by(Sample, Functional_Class) %>%
  dplyr::summarise(
    mean_rate = mean(Class_Rate, na.rm = TRUE),
    se_rate = sd(Class_Rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = 'drop'
  )

functional_class_plot <- ggplot(functional_class_summary, aes(x = Functional_Class, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(
    x = "SnpEff Functional Class",
    y = "Mean Error Rate \nper Total Reads"
  ) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))

# 7G. SnpEff Effect Type Rate Plot
# Choose top N effects to plot based on overall abundance
top_effects <- effect_type_rate_data %>%
  dplyr::group_by(Effect_Type) %>%
  dplyr::summarise(total_alt_ad = sum(Total_Alt_Effect, na.rm = TRUE)) %>%
  pull(Effect_Type)

effect_type_summary <- effect_type_rate_data %>%
  dplyr::filter(Effect_Type %in% top_effects) %>%
  dplyr::group_by(Sample, Effect_Type) %>%
  dplyr::summarise(
    mean_rate = mean(Effect_Rate, na.rm = TRUE),
    se_rate = sd(Effect_Rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = 'drop'
  )

effect_type_plot <- ggplot(effect_type_summary, aes(x = Effect_Type, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(
    x = "SnpEff Effect Type",
    y = "Mean Error Rate \nper Total Reads"
  ) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 1))


# --- [FIXED] 7H. SnpEff LoF Rate Plot ---
# This new section correctly creates the lof_rate_plot
lof_rate_summary <- lof_rate_data %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(
    mean_rate = mean(LoF_Rate, na.rm = TRUE),
    se_rate = sd(LoF_Rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = 'drop'
  )

lof_rate_plot <- ggplot(lof_rate_summary, aes(x = Sample, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  geom_jitter(
    data = lof_rate_data,
    aes(x = Sample, y = LoF_Rate),
    width = 0.1,
    height = 0
  ) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(
    subtitle = paste("Wilcoxon p-value:", format.pval(stat_test_lof_rate$p, digits = 3)),
    x = "Sample",
    y = "Mean Rate of LoF Alleles\nper Total Reads"
  ) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(legend.position = "none")


# --- [FIXED] 7I. SnpEff NMD Rate Plot ---
# This section was 7H and was mixing data. It now correctly plots NMD data.
nmd_rate_summary <- nmd_rate_data %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(
    mean_rate = mean(NMD_Rate, na.rm = TRUE),
    se_rate = sd(NMD_Rate, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = 'drop'
  )

nmd_rate_plot <- ggplot(nmd_rate_summary, aes(x = Sample, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate),
    width = 0.2,
    position = position_dodge(0.7)
  ) +
  geom_jitter(
    data = nmd_rate_data, # <-- FIXED
    aes(x = Sample, y = NMD_Rate), # <-- FIXED
    width = 0.1,
    height = 0
  ) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(
    subtitle = paste("Wilcoxon p-value:", format.pval(stat_test_nmd_rate$p, digits = 3)), # <-- FIXED
    x = "Sample",
    y = "Mean Rate of NMD Alleles\nper Total Reads" # <-- FIXED
  ) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(legend.position = "none")


# --- 8C. Statistical Analysis by Quantile ---
cat("\n--- Running Statistical Analysis by Gene Length Quantile ---\n")

# Allele Proportions by Quantile (per replicate)
quantile_replicate_proportions <- variants_by_quantile %>%
  dplyr::group_by(Quantile, Sample, Replicate) %>%
  dplyr::summarise(
    Total_Alt = sum(AD_Alt),
    Total_Ref = sum(AD_Ref),
    Total_Value = (Total_Alt + Total_Ref),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(Proportion_Alt = Total_Alt / (Total_Alt + Total_Ref))

quantile_stat_test_alt_proportion <- quantile_replicate_proportions %>%
  group_by(Quantile) %>%
  rstatix::wilcox_test(Proportion_Alt ~ Sample, detailed = TRUE)


# Chi-squared Test for Allele Read Proportions by Quantile
quantile_alt_allele_proportions_by_replicate <- variants_by_quantile %>%
  dplyr::filter(Variant_Type %in% c("SNP", "Insertion", "Deletion")) %>%
  dplyr::group_by(Quantile, Sample, Replicate, Variant_Type) %>%
  dplyr::summarise(Total_Alt_Reads = sum(AD_Alt), .groups = 'drop') %>%
  dplyr::group_by(Quantile, Sample, Replicate) %>%
  dplyr::mutate(Proportion = Total_Alt_Reads / sum(Total_Alt_Reads)) %>%
  dplyr::ungroup()

# Wilcoxon Test on Global Error Rate by Quantile
quantile_error_rate_data <- variants_by_quantile %>%
  dplyr::group_by(Quantile, Sample, Replicate) %>%
  dplyr::summarise(
    Error_Rate = sum(AD_Alt) / (sum(AD_Alt) + sum(AD_Ref)),
    .groups = 'drop'
  )
quantile_stat_test_error_rate <- quantile_error_rate_data %>%
  group_by(Quantile) %>%
  rstatix::wilcox_test(Error_Rate ~ Sample, detailed = TRUE)


# 8D. Visualizations by Quantile

# Faceted Global Allele Proportions
quantile_replicate_proportion_summary <- quantile_replicate_proportions %>%
  tidyr::pivot_longer(cols = c(Total_Alt, Total_Ref), names_to = "Allele_Type", values_to = "Count") %>%
  dplyr::group_by(Quantile, Sample) %>%
  dplyr::mutate(Proportion = Count / sum(Count)) %>%
  dplyr::group_by(Quantile, Sample, Allele_Type) %>%
  dplyr::summarise(
    mean_prop = mean(Proportion),
    se_prop = sd(Proportion) / sqrt(dplyr::n()),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(Allele_Type = factor(Allele_Type, levels = c("Total_Ref", "Total_Alt")))

p_values_fisher <- quantile_stat_test_alt_proportion %>% dplyr::mutate(p_label = paste("p =", format.pval(p, digits = 2)))

quantile_proportion_plot <- ggplot(quantile_replicate_proportion_summary, aes(x = Sample, y = mean_prop, fill = Allele_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
    width = 0.2,
    position = position_dodge(0.8)
  ) +
  facet_wrap(~ Quantile,axes = "all_y",labeller = label_wrap_gen(width=20)) +
  geom_text(data = p_values_fisher, aes(x = Inf, y = Inf, label = p_label), hjust = 1.1, vjust = 1.5, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Total_Ref" = "#0000FF", "Total_Alt" = "#FF0000"), labels = c("Reference", "Alternate")) +
  labs(x = "Sample", y = "Mean Proportion of Total Reads", fill = "Allele Type", caption = n_value_label_string) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(plot.caption = element_text(hjust = 1, size = 10))

print(ggplotly(quantile_proportion_plot))

# Faceted Allele Read Proportions by Variant Type
quantile_alt_allele_proportion_summary <- quantile_alt_allele_proportions_by_replicate %>%
  dplyr::group_by(Quantile, Sample, Variant_Type) %>%
  dplyr::summarise(
    mean_prop = mean(Proportion),
    se_prop = sd(Proportion) / sqrt(dplyr::n()),
    .groups = 'drop'
  )

quantile_allele_read_plot <- ggplot(quantile_alt_allele_proportion_summary, aes(x = Sample, y = mean_prop, fill = Variant_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop),
    width = 0.2,
    position = position_dodge(0.8)
  ) +
  facet_wrap(~ Quantile, axes = "all_y",labeller = label_wrap_gen(width=20)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("SNP" = "#00A087FF", "Deletion" = "#f39b7f", "Insertion" = "#4dbbd5"), name = "Alternate Allele Type") +
  labs(x = "Sample", y = "Mean Proportion of Alternate Reads", fill = "Allele Category", caption = n_value_label_string) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(plot.caption = element_text(hjust = 1, size = 10))

print(ggplotly(quantile_allele_read_plot))

# Faceted Global Error Rate
quantile_error_rate_plot <- ggplot(quantile_error_rate_data, aes(x = Sample, y = Error_Rate, fill = Sample)) +
  geom_boxplot(staplewidth = 0.5, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.1, height = 0) +
  facet_wrap(~ Quantile,axes = "all_y",labeller = label_wrap_gen(width=20)) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.15, 0.15))) +  scale_x_discrete(expand=c(0.3,0.3))+
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(x = "Sample", y = "Percentage \nNon-Reference Reads/Total Reads", caption = n_value_label_string) +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(legend.position = "none", plot.caption = element_text(hjust = 1, size = 10))

print(ggplotly(quantile_error_rate_plot))

# --- 8E. SUMMARY TABLE: Unique Variant Counts by Replicate ---
cat("\n--- Generating Summary Table of Variant Counts by Replicate ---\n")

# Total variants per replicate
total_counts_reps <- variants_filtered %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::summarise(Count = dplyr::n_distinct(CHROM, POS, REF, ALT), .groups = 'drop') %>%
  dplyr::mutate(Analysis = "Total Filtered Variants")

# Quantile variants per replicate
quantile_counts_reps <- variants_by_quantile %>%
  dplyr::group_by(Sample, Replicate, Quantile) %>%
  dplyr::summarise(Count = dplyr::n_distinct(CHROM, POS, REF, ALT), .groups = 'drop') %>%
  dplyr::rename(Analysis = Quantile)

# Combine and pivot
summary_table <- bind_rows(total_counts_reps, quantile_counts_reps) %>%
  dplyr::select(Analysis, Sample, Replicate, Count) %>%
  tidyr::pivot_wider(names_from = Replicate, values_from = Count, names_prefix = "Replicate_") %>%
  dplyr::arrange(Analysis, Sample)

# Display the formatted table
kable(summary_table, caption = "Summary of Unique Variant Counts by Replicate Used in Analyses") %>%
  kable_styling(bootstrap_options = "striped", full_width = FALSE) %>%
  scroll_box(width = "100%")

# --- 8F. INTERMEDIARY DATA TABLES ---
cat("\n--- Displaying and Saving Intermediary Data Tables ---\n")

# Define a function to display and save tables
display_and_save_outputs <- function(df, caption, filename_base, dir) {
  # Display kable in notebook
  cat("\n--- ", caption, " ---\n")
  print(kable(df, caption = caption) %>% kable_styling(bootstrap_options = "striped", full_width = FALSE))
  
  # Save as CSV
  write_csv(df, file.path(dir, paste0(filename_base, ".csv")))
  
  # Save as Word
  ft <- flextable(df)
  ft <- autofit(ft)
  doc <- read_docx()
  doc <- body_add_flextable(doc, value = ft)
  print(doc, target = file.path(dir, paste0(filename_base, ".docx")))
}

# --- Define Output Directories ---
output_dir <- "../output/Haplotype_Caller_variants"
tables_dir <- file.path(output_dir, "tables")
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

# Global analysis tables
display_and_save_outputs(replicate_proportions, "Global Allele Proportions", "global_allele_proportions", tables_dir)
display_and_save_outputs(alt_allele_proportions_by_replicate, "Allele Proportions by Variant Type", "alt_allele_proportions_by_type", tables_dir)
display_and_save_outputs(error_rate_data, "Global Error Rate", "global_error_rate", tables_dir)
display_and_save_outputs(mismatch_error_data, "Mismatch Error Rate", "mismatch_error_rate", tables_dir)
display_and_save_outputs(functional_class_rate_data, "Functional Class Rate", "functional_class_rate", tables_dir)
display_and_save_outputs(effect_type_rate_data, "Effect Type Rate", "effect_type_rate", tables_dir)
display_and_save_outputs(nmd_rate_data, "NMD Rate", "nmd_rate", tables_dir)
display_and_save_outputs(lof_rate_data, "LoF Rate", "lof_rate", tables_dir)

# Quantile analysis tables
display_and_save_outputs(quantile_replicate_proportions, "Quantile Allele Proportions", "quantile_allele_proportions", tables_dir)
display_and_save_outputs(quantile_alt_allele_proportions_by_replicate, "Quantile Allele Proportions by Type", "quantile_alt_allele_proportions_by_type", tables_dir)
display_and_save_outputs(quantile_error_rate_data, "Quantile Global Error Rate", "quantile_global_error_rate", tables_dir)


# --- 9. SAVE ALL RESULTS ---

# --- Define Output Directories ---
output_dir <- "../output/Haplotype_Caller_variants"
tables_dir <- file.path(output_dir, "tables")
stat_test_dir <- file.path(output_dir, "stat_test")
figures_dir <- file.path(output_dir, "figures")

# Create directories if they don't exist
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
if (!dir.exists(stat_test_dir)) dir.create(stat_test_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# --- Save Summary Table ---

write_csv(summary_table, file.path(tables_dir, "summary_variant_counts_for_analysis.csv"))
# Save summary table as Word document with landscape orientation
summary_ft <- flextable(summary_table)
summary_ft <- autofit(summary_ft)
# Define landscape page properties
landscape_prop <- prop_section(page_size = page_size(orient = "landscape"))
# Save the flextable to a docx file with landscape settings
save_as_docx(summary_ft, path = file.path(tables_dir, "summary_variant_counts_for_analysis.docx"), pr_section = landscape_prop)

cat("\nSuccessfully saved summary and intermediary tables to:", tables_dir, "\n")

# --- Save Plots ---
# Save Global Plots
ggsave(filename = file.path(figures_dir, "global_allele_proportions.svg"), plot = proportion_plot, width = 4.5, height = 4)
ggsave(filename = file.path(figures_dir, "alt_allele_proportions_by_type.svg"), plot = allele_read_plot, width = 5, height = 4)
ggsave(filename = file.path(figures_dir, "global_error_rate.svg"), plot = error_rate_plot, width = 2.625, height = 3)
ggsave(filename = file.path(figures_dir, "mismatch_error_rate_all.svg"), plot = mismatch_plots, width = 5, height = 3)
ggsave(filename = file.path(figures_dir, "functional_class_rate.svg"), plot = functional_class_plot, width = 5, height = 4)
ggsave(filename = file.path(figures_dir, "effect_type_rate.svg"), plot = effect_type_plot, width = 8, height = 6)
ggsave(filename = file.path(figures_dir, "nmd_rate.svg"), plot = nmd_rate_plot, width = 3.5, height = 4)
ggsave(filename = file.path(figures_dir, "lof_rate.svg"), plot = lof_rate_plot, width = 3.5, height = 4) # <-- This will now work

# Save Quantile Plots
ggsave(filename = file.path(figures_dir, "quantile_allele_proportions.svg"), plot = quantile_proportion_plot, width = 5, height = 5)
ggsave(filename = file.path(figures_dir, "quantile_alt_allele_proportions.svg"), plot = quantile_allele_read_plot, width = 5, height = 7)
ggsave(filename = file.path(figures_dir, "quantile_global_error_rate.svg"), plot = quantile_error_rate_plot, width = 4, height = 4)

cat("\nSuccessfully saved all plots to:", figures_dir, "\n")


# --- Save Statistical Results ---
write_csv(stat_test_alt_proportion, file.path(stat_test_dir, "summary_global_allele_proportions_wilcoxon.csv"))
write_csv(stat_test_alt_type_proportions, file.path(stat_test_dir, "summary_global_alt_allele_type_proportions_wilcoxon.csv"))
write_csv(stat_test_error_rate, file.path(stat_test_dir, "summary_global_error_rate_wilcoxon.csv"))
write_csv(stat_test_mismatch_error_rate, file.path(stat_test_dir, "summary_mismatch_error_rate_wilcoxon.csv"))
write_csv(stat_test_functional_class_rate, file.path(stat_test_dir, "summary_functional_class_rate_wilcoxon.csv"))
write_csv(stat_test_effect_type_rate, file.path(stat_test_dir, "summary_effect_type_rate_wilcoxon.csv"))
write_csv(stat_test_nmd_rate, file.path(stat_test_dir, "summary_nmd_rate_wilcoxon.csv"))
write_csv(stat_test_lof_rate, file.path(stat_test_dir, "summary_lof_rate_wilcoxon.csv"))
write_csv(quantile_stat_test_error_rate, file.path(stat_test_dir, "summary_quantile_error_rate_wilcoxon.csv"))

cat("Successfully saved all statistical summaries to:", stat_test_dir, "\n")


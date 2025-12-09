# ==============================================================================
# Title: Analysis Variants Genome (HaplotypeCaller)
# Converted from Quarto (.qmd) to R Script (.R)
#
# Goal of the Project:
# The goal of the project was to assess whether a particular mutation in the 
# nematode C. elegans caused more transcriptional errors.
#
# Experimental Context:
# - UMIs were used to tag mRNA extracted from 3 biological repeats of:
#   1. Control (wild type) C. elegans (N2)
#   2. Mutant C. elegans (PRDE1)
# - mRNA tagged with UMIs using SMART-Seq Total RNA Pico Input.
# - Sequencing: Novaseq 6000, paired-end, 50-80 million reads, 100 bp.
#
# Analysis Pipeline Summary:
# 1. QC (FastQC)
# 2. Fastq to uBAM (FGBIO)
# 3. Extract UMIs
# 4. Alignment (STAR)
# 5. Merge, Group, and De-duplicate (UMI-tools/FGBIO)
# 6. Call Variants (HaplotypeCaller)
# 7. Annotate (SnpEff)
# 8. Differential Variant Analysis (R - This script)
# ==============================================================================

# --- 1. SETUP: Initialization and Libraries ---

rm(list = ls()) # Clear environment

# Load plotting and analysis libraries
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
library(readxl) 
library(svglite) 

# Set fonts (Windows specific - adjust if on Linux/Mac)
loadfonts(device = "win")

# Set plotting theme
old <- theme_set(theme_cowplot(font_size = 12, font_family = "Arial"))


# ==============================================================================
# VCF File Analysis and Comparison Script GLOBAL
# ==============================================================================
# This script reads VCF files for two samples (N2 and PRDE1), extracts variant 
# information (AD, DP, IMPACT), filters variants, and performs statistical comparisons.

# --- 2. DATA LOADING: Define File Paths ---

# IMPORTANT: Replace these with the actual paths to your VCF files.
n2_files <- c("Annotation/N2.30min.HS.1.consensus.ann.vcf", 
              "Annotation/N2.30min.HS.2.consensus.ann.vcf", 
              "Annotation/N2.30min.HS.3.consensus.ann.vcf")
prde1_files <- c("Annotation/PRDE1.30min.HS.1.consensus.ann.vcf", 
                 "Annotation/PRDE1.30min.HS.2.consensus.ann.vcf", 
                 "Annotation/PRDE1.30min.HS.3.consensus.ann.vcf")

# Check if files exist to prevent errors
all_files <- c(n2_files, prde1_files)
if (!all(file.exists(all_files))) {
  warning("One or more VCF files were not found. Please check your file paths.")
}


# --- 3. PROCESSING FUNCTION: Extract Data from a VCF File ---

#' Process a single VCF file to extract relevant variant data.
process_vcf <- function(file_path, sample_name, replicate_num){
  
  if(!file.exists(file_path)) return(NULL) # Safety check
  
  vcf <- read.vcfR(file_path, verbose = FALSE)
  tidy_vcf <- vcfR2tidy(vcf,
                        info_fields = c("ANN", "LOF", "NMD","AC"),
                        format_fields = c("AD", "DP"),
                        single_frame = TRUE)
  
  processed_data <- tidy_vcf$dat %>%
    as_tibble() %>%
    separate_rows(ALT, sep = ",") %>%
    group_by(CHROM, POS, REF, Indiv) %>%
    mutate(allele_num = row_number()) %>%
    ungroup() %>%
    mutate(
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
    # --- STEP CAUSING DUPLICATES ---
    separate_rows(ANN_allele, sep = "&") %>% 
    # -------------------------------
  mutate(
    Effect_Type = map_chr(str_split(ANN_allele, "\\|"), ~ if(length(.x) >= 2) .x[[2]] else NA_character_),
    IMPACT = map_chr(str_split(ANN_allele, "\\|"), ~ if(length(.x) >= 3) .x[[3]] else NA_character_),
    GENE = map_chr(str_split(ANN_allele, "\\|"), ~ if(length(.x) >= 4) .x[[4]] else NA_character_),
    Functional_Class = case_when(
      str_detect(Effect_Type, "missense_variant") ~ "MISSENSE",
      str_detect(Effect_Type, "synonymous_variant") ~ "SILENT",
      str_detect(Effect_Type, "stop_gained") ~ "NONSENSE",
      TRUE ~ "OTHER"
    ),
    LoF = !is.na(LOF),
    NMD = !is.na(NMD)
  ) %>%
    mutate(
      AD_Ref = ifelse(is.na(AD_Ref), 0, AD_Ref),
      AD_Alt = ifelse(is.na(AD_Alt), 0, AD_Alt),
      gt_DP = as.numeric(gt_DP)
    ) %>%
    select(ANN, Sample, Replicate, CHROM, POS, REF, ALT, GENE, IMPACT, Effect_Type, Functional_Class, LoF, NMD, AD_Ref, AD_Alt, gt_DP, gt_AD, AC) %>%
    rename(Read_Depth = gt_DP) %>%
    mutate(Variant_Type = case_when(
      nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNP",
      nchar(REF) > nchar(ALT) ~ "Deletion",
      nchar(REF) < nchar(ALT) ~ "Insertion",
      TRUE ~ "Other"
    )) %>%
    
    # --- FIX: DEDUPLICATION LOGIC ---
    # Assign numeric rank to IMPACT to pick the "worst" one
    mutate(Impact_Rank = case_when(
      IMPACT == "HIGH" ~ 4,
      IMPACT == "MODERATE" ~ 3,
      IMPACT == "LOW" ~ 2,
      IMPACT == "MODIFIER" ~ 1,
      TRUE ~ 0
    )) %>%
    # Group by the unique genomic variant ID
    group_by(CHROM, POS, REF, ALT, Sample, Replicate) %>%
    # Sort by Impact (highest first) and take the top one
    arrange(desc(Impact_Rank)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-Impact_Rank) # Clean up helper column
  # --------------------------------
  
  return(processed_data)
}

# --- 4. COMBINE DATA: Process All VCF Files ---

file_metadata <- tibble(
  path = c(n2_files, prde1_files),
  sample = rep(c("N2", "PRDE1"), each = 3),
  replicate = rep(1:3, times = 2)
)

# Process files (using pmap to iterate over metadata)
# Note: Ensure files exist before running this chunk
if (all(file.exists(all_files))) {
  all_variants_raw <- pmap_dfr(file_metadata, ~process_vcf(..1, ..2, ..3))
  cat("Successfully processed", nrow(all_variants_raw), "total variants from all files.\n")
  
  all_variants_raw <- all_variants_raw |> mutate(Variant_ratio = AD_Alt/Read_Depth)
} else {
  stop("Files missing. Cannot proceed with variant processing.")
}


# --- 5. FILTERING: Apply Filters for Impact and Read Depth ---

# Define filtering criteria
MIN_READ_DEPTH <- 10
DESIRED_IMPACTS <- c("HIGH", "MODERATE", "LOW")

# Apply filters sequentially
variants_filtered <- all_variants_raw %>% 
  filter(Read_Depth >= MIN_READ_DEPTH) %>%
  filter(IMPACT %in% DESIRED_IMPACTS)  |> 
  filter(Variant_ratio < 1.0) |>  
  filter(Variant_ratio != 0.5) 

cat("Filtered down to", nrow(variants_filtered), "total variant records matching criteria.\n")

# Export filtered lists
list.variants.Filtered = variants_filtered |> as.data.frame() |>  group_by(Sample) |>  group_split()
# Note: Ensuring output directory exists before writing
if (!dir.exists("Variant_Analysis_3")) dir.create("Variant_Analysis_3")

final.list.varians.filtered = list("N2" = list.variants.Filtered[[1]], "PRDE1" = list.variants.Filtered[[2]])
writexl::write_xlsx(final.list.varians.filtered, path = "Variant_Analysis_3/Tables.Variants.Sites.Per.Samples.xlsx")


# --- 6. STATISTICAL ANALYSES ---

# 6A. Global Allele Proportions per Replicate
cat("\n--- Running Global Allele Proportion Analysis per Replicate ---\n")
replicate_proportions <- variants_filtered %>%
  group_by(Sample, Replicate) %>%
  summarise(
    Total_Alt = sum(AD_Alt),
    Total_Ref = sum(AD_Ref),
    .groups = 'drop'
  ) %>%
  mutate(
    Total_Reads = Total_Alt + Total_Ref,
    Proportion_Alt = Total_Alt / Total_Reads,
    Proportion_Ref = Total_Ref / Total_Reads
  )

# Tests
stat_test_alt_proportion <- replicate_proportions %>% rstatix::wilcox_test(Proportion_Alt ~ Sample, detailed = TRUE)
stat_t_test_alt_proportion <- replicate_proportions %>% rstatix::t_test(Proportion_Alt ~ Sample, detailed = TRUE)
print(stat_test_alt_proportion)


# 6B. Allele Read Proportions by Variant Type per Replicate
cat("\n--- Running Allele Read Proportion by Variant Type per Replicate ---\n")
alt_allele_proportions_by_replicate <- variants_filtered %>%
  filter(Variant_Type %in% c("SNP", "Insertion", "Deletion")) %>%
  group_by(Sample, Replicate, Variant_Type) %>%
  summarise(Total_Alt_Reads = sum(AD_Alt), .groups = 'drop') %>%
  group_by(Sample, Replicate) %>%
  mutate(Proportion = Total_Alt_Reads / sum(Total_Alt_Reads)) %>%
  ungroup()

stat_test_alt_type_proportions <- alt_allele_proportions_by_replicate %>%
  group_by(Variant_Type) %>%
  rstatix::wilcox_test(Proportion ~ Sample, detailed = TRUE)

stat_t_test_alt_type_proportions <- alt_allele_proportions_by_replicate %>%
  group_by(Variant_Type) %>%
  rstatix::t_test(Proportion ~ Sample, detailed = TRUE)


# 6D. Global Error Rate
cat("\n--- Running Wilcoxon Test on Global Error Rate ---\n")
error_rate_data <- variants_filtered %>%
  group_by(Sample, Replicate) %>%
  summarise(
    Total_Alt = sum(AD_Alt),
    Total_Ref = sum(AD_Ref),
    Error_Rate = Total_Alt / (Total_Alt + Total_Ref),
    .groups = 'drop'
  )
stat_test_error_rate <- error_rate_data %>% rstatix::wilcox_test(Error_Rate ~ Sample, detailed = TRUE)
stat_t_test_error_rate <- error_rate_data %>% rstatix::t_test(Error_Rate ~ Sample, detailed = TRUE)


# --- Common Denominator for Rate Calculations ---
total_depth_per_replicate <- variants_filtered %>%
  group_by(Sample, Replicate) %>%
  summarise(Total_Replicate_Depth = sum(Read_Depth), .groups = 'drop')


# 6E. Mismatch Error Rate
cat("\n--- Running Wilcoxon Test on Mismatch Error Rate ---\n")
mismatch_error_data <- variants_filtered %>%
  filter(Variant_Type == "SNP") %>%
  mutate(Mismatch_Type = paste0(REF, ">", ALT)) %>%
  group_by(Sample, Replicate, Mismatch_Type) %>%
  summarise(Total_Alt_Mismatch = sum(AD_Alt), .groups = 'drop') %>%
  left_join(total_depth_per_replicate, by = c("Sample", "Replicate")) %>%
  mutate(Mismatch_Error_Rate = Total_Alt_Mismatch / Total_Replicate_Depth)

stat_test_mismatch_error_rate <- mismatch_error_data %>%
  group_by(Mismatch_Type) %>%
  filter(sum(Mismatch_Error_Rate[Sample == "N2"]) > 0 & sum(Mismatch_Error_Rate[Sample == "PRDE1"]) > 0) %>%
  rstatix::wilcox_test(Mismatch_Error_Rate ~ Sample, detailed = TRUE)

stat_t_test_mismatch_error_rate <- mismatch_error_data %>%
  group_by(Mismatch_Type) %>%
  filter(sum(Mismatch_Error_Rate[Sample == "N2"]) > 0 & sum(Mismatch_Error_Rate[Sample == "PRDE1"]) > 0) %>%
  rstatix::t_test(Mismatch_Error_Rate ~ Sample, detailed = TRUE)


# 6F. SnpEff Functional Class Rate
cat("\n--- Running Wilcoxon Test on SnpEff Functional Class Rate ---\n")
functional_class_rate_data <- variants_filtered %>%
  group_by(Sample, Replicate, Functional_Class) %>%
  summarise(Total_Alt_Class = sum(AD_Alt), .groups = 'drop') %>%
  left_join(total_depth_per_replicate, by = c("Sample", "Replicate")) %>%
  mutate(Class_Rate = Total_Alt_Class / Total_Replicate_Depth)

stat_test_functional_class_rate <- functional_class_rate_data %>%
  group_by(Functional_Class) %>%
  filter(sum(Class_Rate[Sample == "N2"], na.rm = TRUE) > 0 & sum(Class_Rate[Sample == "PRDE1"], na.rm = TRUE) > 0) %>%
  rstatix::wilcox_test(Class_Rate ~ Sample, detailed = TRUE)

stat_t_test_functional_class_rate <- functional_class_rate_data %>%
  group_by(Functional_Class) %>%
  filter(sum(Class_Rate[Sample == "N2"], na.rm = TRUE) > 0 & sum(Class_Rate[Sample == "PRDE1"], na.rm = TRUE) > 0) %>%
  rstatix::t_test(Class_Rate ~ Sample, detailed = TRUE)


# 6G. SnpEff Effect Type Rate
cat("\n--- Running Wilcoxon Test on SnpEff Effect Type Rate ---\n")
effect_type_rate_data <- variants_filtered %>%
  group_by(Sample, Replicate, Effect_Type) %>%
  summarise(Total_Alt_Effect = sum(AD_Alt), .groups = 'drop') %>% 
  left_join(total_depth_per_replicate, by = c("Sample", "Replicate")) %>%
  mutate(Effect_Rate = Total_Alt_Effect / Total_Replicate_Depth)

stat_test_effect_type_rate <- effect_type_rate_data %>%
  group_by(Effect_Type) %>%
  filter(sum(Effect_Rate[Sample == "N2"], na.rm = TRUE) > 0 & sum(Effect_Rate[Sample == "PRDE1"], na.rm = TRUE) > 0) %>%
  rstatix::wilcox_test(Effect_Rate ~ Sample, detailed = TRUE)


# 6H. SnpEff NMD Rate
cat("\n--- Running Wilcoxon Test on SnpEff NMD Rate ---\n")
nmd_alt_counts_per_replicate <- variants_filtered %>%
  filter(NMD == TRUE) %>%
  group_by(Sample, Replicate) %>%
  summarise(Total_Alt_NMD = sum(AD_Alt), .groups = 'drop')

nmd_rate_data <- total_depth_per_replicate %>%
  full_join(nmd_alt_counts_per_replicate, by = c("Sample", "Replicate")) %>%
  mutate(Total_Alt_NMD = coalesce(Total_Alt_NMD, 0L)) %>%
  mutate(NMD_Rate = Total_Alt_NMD / Total_Replicate_Depth)

stat_test_nmd_rate <- nmd_rate_data %>% rstatix::wilcox_test(NMD_Rate ~ Sample, detailed = TRUE)
stat_t_test_nmd_rate <- nmd_rate_data %>% rstatix::t_test(NMD_Rate ~ Sample, detailed = TRUE)


# 6I. SnpEff LoF Rate
cat("\n--- Running Wilcoxon Test on SnpEff LoF Rate ---\n")
lof_alt_counts_per_replicate <- variants_filtered %>%
  filter(LoF == TRUE) %>%
  group_by(Sample, Replicate) %>%
  summarise(Total_Alt_LoF = sum(AD_Alt), .groups = 'drop')

lof_rate_data <- total_depth_per_replicate %>%
  full_join(lof_alt_counts_per_replicate, by = c("Sample", "Replicate")) %>%
  mutate(Total_Alt_LoF = coalesce(Total_Alt_LoF, 0L)) %>%
  mutate(LoF_Rate = Total_Alt_LoF / Total_Replicate_Depth)

stat_test_lof_rate <- lof_rate_data %>% rstatix::wilcox_test(LoF_Rate ~ Sample, detailed = TRUE)
stat_t_test_lof_rate <- lof_rate_data %>% rstatix::t_test(LoF_Rate ~ Sample, detailed = TRUE)


# 6J. Combined Deleterious Variants
cat("\n--- Running Wilcoxon Test on Combined Deleterious Variants Rate ---\n")
deleterious_effects <- c("missense_variant", "stop_gained", "frameshift_variant")

deleterious_alt_counts_per_replicate <- variants_filtered %>%
  filter(str_detect(Effect_Type, paste(deleterious_effects, collapse = "|"))) %>%
  group_by(Sample, Replicate) %>%
  summarise(Total_Alt_Deleterious = sum(AD_Alt), .groups = 'drop')

deleterious_rate_data <- total_depth_per_replicate %>%
  full_join(deleterious_alt_counts_per_replicate, by = c("Sample", "Replicate")) %>%
  mutate(Total_Alt_Deleterious = coalesce(Total_Alt_Deleterious, 0L)) %>%
  mutate(Deleterious_Rate = Total_Alt_Deleterious / Total_Replicate_Depth)

stat_test_deleterious_rate <- deleterious_rate_data %>% rstatix::wilcox_test(Deleterious_Rate ~ Sample, detailed = TRUE)
stat_t_test_deleterious_rate <- deleterious_rate_data %>% rstatix::t_test(Deleterious_Rate ~ Sample, detailed = TRUE)


# --- 7. ANALYSIS BY GENE LENGTH QUANTILE ---

cat("\n--- Setting up Analysis by Gene Length Quantile ---\n")
# 8A. Load GTF and Define Gene Quantiles
# Note: HARDCODED PATH. Update for your system.
gtf_path <- "Y:/Johnny/Roswell_Projects/Sehee_mRNA_piRNA/Github.Sehee.paper/Data/input/Caenorhabditis_elegans.WBcel235.111.gtf.gz"

if(file.exists(gtf_path)) {
  gtf <- rtracklayer::import(con = gtf_path)
  
  protein_coding_genes <- gtf |> plyranges::filter(gene_biotype == "protein_coding") |> plyranges::filter(type == "gene")
  
  # Define quantiles based on gene width (length)
  quantiles <- quantile(GenomicRanges::width(protein_coding_genes), probs = c(0, 0.25, 0.5, 0.75, 1.0))
  
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
  
  variants_by_quantile <- bind_rows(
    variants_filtered_quantile_1 %>% mutate(Quantile = q1_label),
    variants_filtered_quantile_2 %>% mutate(Quantile = q2_label),
    variants_filtered_quantile_3 %>% mutate(Quantile = q3_label),
    variants_filtered_quantile_4 %>% mutate(Quantile = q4_label)
  ) %>%
    mutate(Quantile = factor(Quantile, levels = quantile_labels))
  
  # 8C. Statistical Analysis by Quantile
  quantile_replicate_proportions <- variants_by_quantile %>%
    group_by(Quantile, Sample, Replicate) %>%
    summarise(
      Total_Alt = sum(AD_Alt),
      Total_Ref = sum(AD_Ref),
      Total_Value = (Total_Alt + Total_Ref),
      .groups = 'drop'
    ) %>%
    mutate(Proportion_Alt = Total_Alt / (Total_Alt + Total_Ref))
  
  quantile_stat_test_alt_proportion <- quantile_replicate_proportions %>%
    group_by(Quantile) %>%
    rstatix::wilcox_test(Proportion_Alt ~ Sample, detailed = TRUE)
  
  quantile_stat_t_test_alt_proportion <- quantile_replicate_proportions %>%
    group_by(Quantile) %>%
    rstatix::t_test(Proportion_Alt ~ Sample, detailed = TRUE)
  
  quantile_alt_allele_proportions_by_replicate <- variants_by_quantile %>%
    filter(Variant_Type %in% c("SNP", "Insertion", "Deletion")) %>%
    group_by(Quantile, Sample, Replicate, Variant_Type) %>%
    summarise(Total_Alt_Reads = sum(AD_Alt), .groups = 'drop') %>%
    group_by(Quantile, Sample, Replicate) %>%
    mutate(Proportion = Total_Alt_Reads / sum(Total_Alt_Reads)) %>%
    ungroup()
  
  quantile_error_rate_data <- variants_by_quantile %>%
    group_by(Quantile, Sample, Replicate) %>%
    summarise(
      Error_Rate = sum(AD_Alt) / (sum(AD_Alt) + sum(AD_Ref)),
      .groups = 'drop'
    )
  quantile_stat_test_error_rate <- quantile_error_rate_data %>%
    group_by(Quantile) %>%
    rstatix::wilcox_test(Error_Rate ~ Sample, detailed = TRUE)
  
  quantile_stat_t_test_error_rate <- quantile_error_rate_data %>%
    group_by(Quantile) %>%
    rstatix::t_test(Error_Rate ~ Sample, detailed = TRUE)
  
} else {
  warning("GTF file not found. Quantile analysis skipped.")
}

# --- 8. VISUALIZATIONS ---

cat("\n--- Generating Plots ---\n")
n_value_label <- tibble(label = "n = 2045-2054")
n_value_label_string <- "n = 2045-2054"

# 7A. Global Allele Proportions Plot
replicate_proportion_summary <- replicate_proportions %>%
  pivot_longer(cols = c(Proportion_Alt, Proportion_Ref), names_to = "Allele_Type_Full", values_to = "Proportion") %>%
  mutate(Allele_Type = ifelse(Allele_Type_Full == "Proportion_Ref", "Total_Ref", "Total_Alt")) %>%
  group_by(Sample, Allele_Type) %>%
  summarise(mean_prop = mean(Proportion), se_prop = sd(Proportion) / sqrt(n()), .groups = 'drop') %>%
  mutate(Allele_Type = factor(Allele_Type, levels = c("Total_Ref", "Total_Alt")))

proportion_plot <- ggplot(replicate_proportion_summary, aes(x = Sample, y = mean_prop, fill = Allele_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop), width = 0.2, position = position_dodge(0.8)) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("Total_Ref" = "#0000FF", "Total_Alt" = "#FF0000"), labels = c("Reference Alleles", "Alternate Alleles")) +
  labs(subtitle = paste("Wilcoxon p-value (Alt Alleles):", format.pval(stat_test_alt_proportion$p, digits = 3)), x = "Sample", y = "Mean Proportion of Total Reads", fill = "Allele Type") +
  theme_cowplot(font_size = 12, font_family = "Arial")

# 7B. Allele Read Proportions by Variant Type Plot
alt_allele_proportion_summary <- alt_allele_proportions_by_replicate %>%
  group_by(Sample, Variant_Type) %>%
  summarise(mean_prop = mean(Proportion), se_prop = sd(Proportion) / sqrt(n()), .groups = 'drop')

allele_read_plot <- ggplot(alt_allele_proportion_summary, aes(x = Sample, y = mean_prop, fill = Variant_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop), width = 0.2, position = position_dodge(0.8)) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("SNP" = "#00A087FF", "Deletion" = "#f39b7f", "Insertion" = "#4dbbd5"), name = "Alternate Allele Type") +
  labs(subtitle = "Comparing proportions of alternate allele types", x = "Sample", y = "Mean Proportion of Alternate Reads") +
  theme_cowplot(font_size = 12, font_family = "Arial")

# 7D. Global Error Rate Plot
error_rate_plot <- ggplot(error_rate_data, aes(x = Sample, y = Error_Rate, fill = Sample)) +
  geom_boxplot(staplewidth = 0.5, outlier.shape = NA, width = 0.3) +
  geom_jitter(width = 0.3, height = 0) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.15, 0.15))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(subtitle = paste("Wilcoxon Test p-value:", format.pval(stat_test_error_rate$p, digits = 3)), x = "Sample", y = "Percentage \nNon-Reference Reads/Total Reads") +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(legend.position = "none")

# 7E. Mismatch Error Rate Plot
mismatch_summary <- mismatch_error_data %>%
  group_by(Sample, Mismatch_Type) %>%
  summarise(mean_rate = mean(Mismatch_Error_Rate), se_rate = sd(Mismatch_Error_Rate) / sqrt(n()), .groups = 'drop')

mismatch_plots <- mismatch_summary |> 
  mutate(Mismatch_Type = fct_relevel(Mismatch_Type,c("A>G", "G>A", "G>C" ,"T>C", "A>C", "A>T", "C>A" ,"C>G" ,"C>T", "G>T", "T>A" ,"T>G"))) |> 
  ggplot(aes(x = Mismatch_Type, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate), width = 0.2, position = position_dodge(0.7)) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(x = "Mismatch Type", y = "Mean Mismatch Error Rate") +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))

# 7F. SnpEff Functional Class Rate Plot
functional_class_summary <- functional_class_rate_data %>%
  group_by(Sample, Functional_Class) %>%
  summarise(mean_rate = mean(Class_Rate, na.rm = TRUE), se_rate = sd(Class_Rate, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

functional_class_plot <- ggplot(functional_class_summary, aes(x = Functional_Class, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate), width = 0.2, position = position_dodge(0.7)) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(x = "SnpEff Functional Class", y = "Mean Error Rate \nper Total Reads") +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))

# 7G. SnpEff Effect Type Rate Plot
top_effects <- effect_type_rate_data %>%
  group_by(Effect_Type) %>%
  summarise(total_alt_ad = sum(Total_Alt_Effect, na.rm = TRUE)) %>%
  pull(Effect_Type)

effect_type_summary <- effect_type_rate_data %>%
  filter(Effect_Type %in% top_effects) %>%
  group_by(Sample, Effect_Type) %>%
  summarise(mean_rate = mean(Effect_Rate, na.rm = TRUE), se_rate = sd(Effect_Rate, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

effect_type_plot <- ggplot(effect_type_summary, aes(x = Effect_Type, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate), width = 0.2, position = position_dodge(0.7)) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(x = "SnpEff Effect Type", y = "Mean Error Rate \nper Total Reads") +
  theme_cowplot(font_size = 12, font_family = "Arial") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 1))

# 7H. NMD and LoF Plots
nmd_rate_summary <- nmd_rate_data %>%
  group_by(Sample) %>%
  summarise(mean_rate = mean(NMD_Rate, na.rm = TRUE), se_rate = sd(NMD_Rate, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

nmd_rate_plot <- ggplot(nmd_rate_summary, aes(x = Sample, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate), width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data = nmd_rate_data, aes(x = Sample, y = NMD_Rate), width = 0.1, height = 0) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(subtitle = paste("Wilcoxon p-value:", format.pval(stat_test_nmd_rate$p, digits = 3)), x = "Sample", y = "Mean Rate of NMD Alleles \nper Total Alleles") +
  theme_cowplot(font_size = 12, font_family = "Arial") + theme(legend.position = "none")

lof_rate_summary <- lof_rate_data %>%
  group_by(Sample) %>%
  summarise(mean_rate = mean(LoF_Rate, na.rm = TRUE), se_rate = sd(LoF_Rate, na.rm = TRUE) / sqrt(n()), .groups = 'drop')

lof_rate_plot <- ggplot(lof_rate_summary, aes(x = Sample, y = mean_rate, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate), width = 0.2, position = position_dodge(0.7)) +
  geom_jitter(data = lof_rate_data, aes(x = Sample, y = LoF_Rate), width = 0.1, height = 0) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::scientific, expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(subtitle = paste("Wilcoxon p-value:", format.pval(stat_test_lof_rate$p, digits = 3)), x = "Sample", y = "Mean Rate of LoF Alleles\nper Total Reads") +
  theme_cowplot(font_size = 12, font_family = "Arial") + theme(legend.position = "none")

deleterious_rate_plot <- ggplot(deleterious_rate_data, aes(x = Sample, y = Deleterious_Rate, fill = Sample)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0) +
  geom_text(data = n_value_label, aes(x = Inf, y = -Inf, label = label), hjust = 1.1, vjust = -1.1, inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
  labs(x = "Sample", y = "Percentage \nNon-Reference Reads/Total Reads") +
  theme_cowplot(font_size = 12, font_family = "Arial") + theme(legend.position = "none")


# 7I. Visualizations by Quantile
if(exists("variants_by_quantile")) {
  
  quantile_replicate_proportion_summary <- quantile_replicate_proportions %>%
    pivot_longer(cols = c(Total_Alt, Total_Ref), names_to = "Allele_Type", values_to = "Count") %>%
    group_by(Quantile, Sample) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    group_by(Quantile, Sample, Allele_Type) %>%
    summarise(mean_prop = mean(Proportion), se_prop = sd(Proportion) / sqrt(n()), .groups = 'drop') %>%
    mutate(Allele_Type = factor(Allele_Type, levels = c("Total_Ref", "Total_Alt")))
  
  p_values_fisher <- quantile_stat_test_alt_proportion %>% mutate(p_label = paste("p =", format.pval(p, digits = 2)))
  
  quantile_proportion_plot <- ggplot(quantile_replicate_proportion_summary, aes(x = Sample, y = mean_prop, fill = Allele_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop), width = 0.2, position = position_dodge(0.8)) +
    facet_wrap(~ Quantile,axes = "all_y",labeller = label_wrap_gen(width=20)) +
    geom_text(data = p_values_fisher, aes(x = Inf, y = Inf, label = p_label), hjust = 1.1, vjust = 1.5, inherit.aes = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Total_Ref" = "#0000FF", "Total_Alt" = "#FF0000"), labels = c("Reference", "Alternate")) +
    labs(x = "Sample", y = "Mean Proportion of Total Reads", fill = "Allele Type", caption = n_value_label_string) +
    theme_cowplot(font_size = 12, font_family = "Arial") + theme(plot.caption = element_text(hjust = 1, size = 10))
  
  quantile_alt_allele_proportion_summary <- quantile_alt_allele_proportions_by_replicate %>%
    group_by(Quantile, Sample, Variant_Type) %>%
    summarise(mean_prop = mean(Proportion), se_prop = sd(Proportion) / sqrt(n()), .groups = 'drop')
  
  quantile_allele_read_plot <- ggplot(quantile_alt_allele_proportion_summary, aes(x = Sample, y = mean_prop, fill = Variant_Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean_prop - se_prop, ymax = mean_prop + se_prop), width = 0.2, position = position_dodge(0.8)) +
    facet_wrap(~ Quantile, axes = "all_y",labeller = label_wrap_gen(width=20)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("SNP" = "#00A087FF", "Deletion" = "#f39b7f", "Insertion" = "#4dbbd5"), name = "Alternate Allele Type") +
    labs(x = "Sample", y = "Mean Proportion of Alternate Reads", fill = "Allele Category", caption = n_value_label_string) +
    theme_cowplot(font_size = 12, font_family = "Arial") + theme(plot.caption = element_text(hjust = 1, size = 10))
  
  quantile_error_rate_plot <- ggplot(quantile_error_rate_data, aes(x = Sample, y = Error_Rate, fill = Sample)) +
    geom_boxplot(staplewidth = 0.5, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.1, height = 0) +
    facet_wrap(~ Quantile,axes = "all_y",labeller = label_wrap_gen(width=20)) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.15, 0.15))) +  scale_x_discrete(expand=c(0.3,0.3))+
    scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
    labs(x = "Sample", y = "Percentage \nNon-Reference Reads/Total Reads", caption = n_value_label_string) +
    theme_cowplot(font_size = 12, font_family = "Arial") + theme(legend.position = "none", plot.caption = element_text(hjust = 1, size = 10))
}


# --- 9. ANALYSIS BY INTRON PRESENCE ---
cat("\n--- Running Analysis by Intron Presence ---\n")

if(exists("gtf")) {
  # 9A. Define Intron Status based on GTF
  gene_exon_counts <- gtf %>%
    as_tibble() %>%
    filter(type == "exon", gene_biotype == "protein_coding") %>%
    group_by(gene_name) %>%
    summarise(n_exons = n(), .groups = 'drop')
  
  genes_with_introns <- gene_exon_counts %>% filter(n_exons > 1) %>% pull(gene_name)
  genes_without_introns <- gene_exon_counts %>% filter(n_exons == 1) %>% pull(gene_name)
  
  # 9B. Filter Variants by Intron Status
  variants_intron_analysis <- variants_filtered %>%
    mutate(Intron_Status = case_when(
      GENE %in% genes_with_introns ~ "Genes with Introns",
      GENE %in% genes_without_introns ~ "Genes without Introns",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Intron_Status))
  
  # 9C. Calculate Error Rates
  intron_error_rates <- variants_intron_analysis %>%
    group_by(Intron_Status, Sample, Replicate) %>%
    summarise(
      Total_Alt = sum(AD_Alt),
      Total_Ref = sum(AD_Ref),
      Error_Rate = Total_Alt / (Total_Alt + Total_Ref),
      .groups = 'drop'
    )
  
  # 9D. Stats
  stat_test_intron_wilcox <- intron_error_rates %>% group_by(Intron_Status) %>% rstatix::wilcox_test(Error_Rate ~ Sample, detailed = TRUE)
  stat_test_intron_ttest <- intron_error_rates %>% group_by(Intron_Status) %>% rstatix::t_test(Error_Rate ~ Sample, detailed = TRUE)
  
  # 9E. Plot
  intron_plot <- ggplot(intron_error_rates, aes(x = Sample, y = Error_Rate, fill = Sample)) +
    geom_boxplot(staplewidth = 0.5, outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.1, height = 0) +
    facet_wrap(~ Intron_Status, axes = "all_y") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.15, 0.15))) +
    scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
    labs(
      subtitle = "Error Rates: Genes with Introns vs Genes without Introns",
      x = "Sample",
      y = "Percentage \nNon-Reference Reads/Total Reads",
      caption = n_value_label_string
    ) +
    theme_cowplot(font_size = 12, font_family = "Arial") +
    theme(legend.position = "none", plot.caption = element_text(hjust = 1, size = 10))
  
} else {
  warning("GTF object missing. Intron analysis skipped.")
}


# --- 10. SAVING RESULTS (TABLES & PLOTS) ---

cat("\n--- Saving Results ---\n")

# Define Output Directories
output_dir <- "Variant_Analysis_3/"
tables_dir <- file.path(output_dir, "tables")
stat_test_dir <- file.path(output_dir, "stat_test")
figures_dir <- file.path(output_dir, "figures")

# Create directories
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)
if (!dir.exists(stat_test_dir)) dir.create(stat_test_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Helper function to save tables
display_and_save_outputs <- function(df, caption, filename_base, dir) {
  write_csv(df, file.path(dir, paste0(filename_base, ".csv")))
  # Word doc generation
  ft <- flextable(df)
  ft <- autofit(ft)
  doc <- read_docx()
  doc <- body_add_flextable(doc, value = ft)
  print(doc, target = file.path(dir, paste0(filename_base, ".docx")))
}

# Save Global Tables
display_and_save_outputs(replicate_proportions, "Global Allele Proportions", "global_allele_proportions", tables_dir)
display_and_save_outputs(alt_allele_proportions_by_replicate, "Allele Proportions by Variant Type", "alt_allele_proportions_by_type", tables_dir)
display_and_save_outputs(error_rate_data, "Global Error Rate", "global_error_rate", tables_dir)
display_and_save_outputs(mismatch_error_data, "Mismatch Error Rate", "mismatch_error_rate", tables_dir)
display_and_save_outputs(functional_class_rate_data, "Functional Class Rate", "functional_class_rate", tables_dir)
display_and_save_outputs(effect_type_rate_data, "Effect Type Rate", "effect_type_rate", tables_dir)
display_and_save_outputs(nmd_rate_data, "NMD Rate", "nmd_rate", tables_dir)
display_and_save_outputs(lof_rate_data, "LoF Rate", "lof_rate", tables_dir)
display_and_save_outputs(deleterious_rate_data, "Combined Deleterious Rate", "deleterious_rate", tables_dir)

if(exists("variants_by_quantile")) {
  display_and_save_outputs(quantile_replicate_proportions, "Quantile Allele Proportions", "quantile_allele_proportions", tables_dir)
  display_and_save_outputs(quantile_alt_allele_proportions_by_replicate, "Quantile Allele Proportions by Type", "quantile_alt_allele_proportions_by_type", tables_dir)
  display_and_save_outputs(quantile_error_rate_data, "Quantile Global Error Rate", "quantile_global_error_rate", tables_dir)
}

if(exists("intron_error_rates")) {
  display_and_save_outputs(intron_error_rates, "Intron vs No-Intron Error Rates", "intron_error_rates", tables_dir)
}

# Save Plots
ggsave(filename = file.path(figures_dir, "global_allele_proportions.svg"), plot = proportion_plot, width = 4.5, height = 4)
ggsave(filename = file.path(figures_dir, "alt_allele_proportions_by_type.svg"), plot = allele_read_plot, width = 5, height = 4)
ggsave(filename = file.path(figures_dir, "global_error_rate.svg"), plot = error_rate_plot, width = 2.625, height = 3)
ggsave(filename = file.path(figures_dir, "mismatch_error_rate_all.svg"), plot = mismatch_plots, width = 5, height = 3)
ggsave(filename = file.path(figures_dir, "functional_class_rate.svg"), plot = functional_class_plot, width = 5, height = 4)
ggsave(filename = file.path(figures_dir, "effect_type_rate.svg"), plot = effect_type_plot, width = 8, height = 6)
ggsave(filename = file.path(figures_dir, "nmd_rate.svg"), plot = nmd_rate_plot, width = 3.5, height = 4)
ggsave(filename = file.path(figures_dir, "lof_rate.svg"), plot = lof_rate_plot, width = 3.5, height = 4)
ggsave(filename = file.path(figures_dir, "deleterious_rate.svg"), plot = deleterious_rate_plot, width = 3.5, height = 4)

if(exists("variants_by_quantile")) {
  ggsave(filename = file.path(figures_dir, "quantile_allele_proportions.svg"), plot = quantile_proportion_plot, width = 5, height = 5)
  ggsave(filename = file.path(figures_dir, "quantile_alt_allele_proportions.svg"), plot = quantile_allele_read_plot, width = 5, height = 7)
  ggsave(filename = file.path(figures_dir, "quantile_global_error_rate.svg"), plot = quantile_error_rate_plot, width = 4, height = 4)
}

if(exists("intron_plot")) {
  ggsave(filename = file.path(figures_dir, "intron_error_rate.svg"), plot = intron_plot, width = 5, height = 4)
}

# Save Stats (Wilcoxon and T-Test)
write_csv(stat_test_alt_proportion, file.path(stat_test_dir, "summary_global_allele_proportions_wilcoxon.csv"))
write_csv(stat_test_alt_type_proportions, file.path(stat_test_dir, "summary_global_alt_allele_type_proportions_wilcoxon.csv"))
write_csv(stat_test_error_rate, file.path(stat_test_dir, "summary_global_error_rate_wilcoxon.csv"))
write_csv(stat_test_mismatch_error_rate, file.path(stat_test_dir, "summary_mismatch_error_rate_wilcoxon.csv"))
write_csv(stat_test_functional_class_rate, file.path(stat_test_dir, "summary_functional_class_rate_wilcoxon.csv"))
write_csv(stat_test_effect_type_rate, file.path(stat_test_dir, "summary_effect_type_rate_wilcoxon.csv"))
write_csv(stat_test_nmd_rate, file.path(stat_test_dir, "summary_nmd_rate_wilcoxon.csv"))
write_csv(stat_test_lof_rate, file.path(stat_test_dir, "summary_lof_rate_wilcoxon.csv"))
write_csv(stat_test_deleterious_rate, file.path(stat_test_dir, "summary_deleterious_rate_wilcoxon.csv"))

write_csv(stat_t_test_alt_proportion, file.path(stat_test_dir, "summary_global_allele_proportions_ttest.csv"))
write_csv(stat_t_test_alt_type_proportions, file.path(stat_test_dir, "summary_global_alt_allele_type_proportions_ttest.csv"))
write_csv(stat_t_test_error_rate, file.path(stat_test_dir, "summary_global_error_rate_ttest.csv"))
write_csv(stat_t_test_mismatch_error_rate, file.path(stat_test_dir, "summary_mismatch_error_rate_ttest.csv"))
write_csv(stat_t_test_functional_class_rate, file.path(stat_test_dir, "summary_functional_class_rate_ttest.csv"))
write_csv(stat_t_test_nmd_rate, file.path(stat_test_dir, "summary_nmd_rate_ttest.csv"))
write_csv(stat_t_test_lof_rate, file.path(stat_test_dir, "summary_lof_rate_ttest.csv"))
write_csv(stat_t_test_deleterious_rate, file.path(stat_test_dir, "summary_deleterious_rate_ttest.csv"))

if(exists("variants_by_quantile")) {
  write_csv(quantile_stat_test_error_rate, file.path(stat_test_dir, "summary_quantile_error_rate_wilcoxon.csv"))
  write_csv(quantile_stat_test_alt_proportion, file.path(stat_test_dir, "summary_quantile_allele_proportions_wilcoxon.csv"))
  write_csv(quantile_stat_t_test_error_rate, file.path(stat_test_dir, "summary_quantile_error_rate_ttest.csv"))
  write_csv(quantile_stat_t_test_alt_proportion, file.path(stat_test_dir, "summary_quantile_allele_proportions_ttest.csv"))
}

if(exists("stat_test_intron_wilcox")) {
  write_csv(stat_test_intron_wilcox, file.path(stat_test_dir, "summary_intron_error_rate_wilcoxon.csv"))
  write_csv(stat_test_intron_ttest, file.path(stat_test_dir, "summary_intron_error_rate_ttest.csv"))
}


# --- 11. DEG SUBSET ANALYSIS ---
# Purpose: Analyze Global and Intron Error Rates specifically for genes identified as 
# differentially expressed (Up-regulated: padj < 0.05, log2FC > 0) in the DESeq2 analysis.

# HARDCODED PATH for DEG data
deg_file_path <- "Y:/Johnny/Roswell_Projects/Sehee_mRNA_piRNA/Github.Sehee.paper/Documents/Tables_Final/DESeq mRNA/relevant_deseq_tables/HeatShock_Recovery_vs_Control_DEG.xlsx"

analyze_deg_subset <- function(sheet_name, label, variants_df, intron_df) {
  
  if(!file.exists(deg_file_path)) {
    warning("DEG file not found: ", deg_file_path)
    return(NULL)
  }
  
  clean_name <- gsub("[^[:alnum:]]", "_", label)
  
  cat(paste0("\nProcessing DEG Subset: ", label, "\n"))
  
  # 1. Read and Filter DEG Data
  # Note: Requires correct sheet name
  tryCatch({
    deg_data <- read_excel(deg_file_path, sheet = sheet_name)
    target_genes <- deg_data %>%
      filter(padj < 0.05, log2FoldChange > 0) %>%
      pull(1) 
    
    cat(paste0("Identified ", length(target_genes), " up-regulated genes (padj < 0.05, log2FC > 0).\n"))
    
    if(length(target_genes) > 0) {
      
      # A. Global Error Rate for DEG Subset
      deg_variants <- variants_df %>% filter(GENE %in% target_genes)
      
      deg_error_rate_data <- deg_variants %>%
        group_by(Sample, Replicate) %>%
        summarise(
          Total_Alt = sum(AD_Alt),
          Total_Ref = sum(AD_Ref),
          Error_Rate = Total_Alt / (Total_Alt + Total_Ref),
          .groups = 'drop'
        )
      
      stat_test <- deg_error_rate_data %>% rstatix::t_test(Error_Rate ~ Sample, detailed = TRUE)
      
      p_global <- ggplot(deg_error_rate_data, aes(x = Sample, y = Error_Rate, fill = Sample)) +
        geom_boxplot(staplewidth = 0.5, outlier.shape = NA, width = 0.3) +
        geom_jitter(width = 0.3, height = 0) +
        scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.15, 0.15))) +
        scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
        labs(
          title = paste0("Global Error Rate\n(", label, ")"),
          subtitle = paste("T-test p-value:", format.pval(stat_test$p, digits = 3)),
          x = "Sample", y = "Error Rate (Non-Ref/Total)"
        ) +
        theme_cowplot(font_size = 12, font_family = "Arial") + theme(legend.position = "none")
      
      ggsave(filename = file.path(figures_dir, paste0("DEG_Global_Error_", clean_name, ".png")), plot = p_global, width = 4, height = 4)
      ggsave(filename = file.path(figures_dir, paste0("DEG_Global_Error_", clean_name, ".svg")), plot = p_global, width = 4, height = 4)
      write_csv(deg_error_rate_data, file.path(tables_dir, paste0("DEG_Global_Error_Data_", clean_name, ".csv")))
      
      # B. Intron Error Rate for DEG Subset
      deg_intron_variants <- intron_df %>% filter(GENE %in% target_genes)
      
      if(nrow(deg_intron_variants) > 0) {
        deg_intron_rates <- deg_intron_variants %>%
          group_by(Intron_Status, Sample, Replicate) %>%
          summarise(
            Total_Alt = sum(AD_Alt),
            Total_Ref = sum(AD_Ref),
            Error_Rate = Total_Alt / (Total_Alt + Total_Ref),
            .groups = 'drop'
          )
        
        write_csv(deg_intron_rates, file.path(tables_dir, paste0("DEG_Intron_Error_Data_", clean_name, ".csv")))
        
        # Plot: Genes WITH Introns
        intron_only_data <- deg_intron_rates %>% filter(Intron_Status == "Genes with Introns")
        if(nrow(intron_only_data) > 0) {
          p_intron_only <- ggplot(intron_only_data, aes(x = Sample, y = Error_Rate, fill = Sample)) +
            geom_boxplot(staplewidth = 0.5, outlier.shape = NA, width = 0.3) +
            geom_jitter(width = 0.3, height = 0) +
            scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.15, 0.15))) +
            scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
            labs(title = paste0("Error Rate: Genes WITH Introns\n(", label, ")"), x = "Sample", y = "Error Rate") +
            theme_cowplot(font_size = 12, font_family = "Arial") + theme(legend.position = "none")
          
          ggsave(filename = file.path(figures_dir, paste0("DEG_With_Intron_Error_", clean_name, ".png")), plot = p_intron_only, width = 4, height = 4)
          ggsave(filename = file.path(figures_dir, paste0("DEG_With_Intron_Error_", clean_name, ".svg")), plot = p_intron_only, width = 4, height = 4)
        }
        
        # Plot: Genes WITHOUT Introns
        no_intron_data <- deg_intron_rates %>% filter(Intron_Status == "Genes without Introns")
        if(nrow(no_intron_data) > 0) {
          p_no_intron <- ggplot(no_intron_data, aes(x = Sample, y = Error_Rate, fill = Sample)) +
            geom_boxplot(staplewidth = 0.5, outlier.shape = NA, width = 0.3) +
            geom_jitter(width = 0.3, height = 0) +
            scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.15, 0.15))) +
            scale_fill_manual(values = c("N2" = "#A0A0A4", "PRDE1" = "#FF0066")) +
            labs(title = paste0("Error Rate: Genes WITHOUT Introns\n(", label, ")"), x = "Sample", y = "Error Rate") +
            theme_cowplot(font_size = 12, font_family = "Arial") + theme(legend.position = "none")
          
          ggsave(filename = file.path(figures_dir, paste0("DEG_No_Intron_Error_", clean_name, ".png")), plot = p_no_intron, width = 4, height = 4)
          ggsave(filename = file.path(figures_dir, paste0("DEG_No_Intron_Error_", clean_name, ".svg")), plot = p_no_intron, width = 4, height = 4)
        }
      }
    } else {
      cat("No genes met the criteria (padj < 0.05, log2FC > 0).\n")
    }
  }, error = function(e) {
    warning("Error reading sheet '", sheet_name, "' from DEG file: ", e$message)
  })
}

# Run Analysis
if(exists("variants_filtered") && exists("variants_intron_analysis")) {
  analyze_deg_subset("N2_30mHS.vs.N2_Ctrl", "N2 HS Up-regulated Genes", variants_filtered, variants_intron_analysis)
  analyze_deg_subset("PRDE1_30mHS.vs.PRDE1_Ctrl", "PRDE1 HS Up-regulated Genes", variants_filtered, variants_intron_analysis)
  analyze_deg_subset("PRDE1_Ctrl.vs.N2_Ctrl", "PRDE1 Ctrl vs N2 Ctrl", variants_filtered, variants_intron_analysis)
  analyze_deg_subset("PRDE1_30mHS.vs.N2_30mHS", "PRDE1 HS vs N2 HS", variants_filtered, variants_intron_analysis)
}
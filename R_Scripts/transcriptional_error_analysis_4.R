# ==============================================================================
# COMBINED ANALYSIS SCRIPT
# - Transcriptional Errors (Global & Gene-level)
# - Frequency Tables (All Sites [0.01 bins] & Variant Only [0.05 bins])
# - 22G Overlap (Fixed Euler Plots)
# ==============================================================================

# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
if (!require("eulerr", quietly = TRUE)) install.packages("eulerr")
if (!require("flextable", quietly = TRUE)) install.packages("flextable") 
if (!require("officer", quietly = TRUE)) install.packages("officer")     

library(tidyverse)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(clusterProfiler) 
library(org.Ce.eg.db)    
library(DOSE)            
library(wbData)          
library(eulerr)
library(flextable)
library(officer)

# ==============================================================================
# 2. CONFIGURATION & SETUP
# ==============================================================================

# --- A. Directories & Files ---
setwd("Y:/Johnny/Roswell_Projects/Sehee_mRNA_piRNA/umi_celegans_consensus")

OUTPUT_DIR <- "Analysis_Results_Combined"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  message(paste("Created output directory:", OUTPUT_DIR))
}

consensus_path <- "Y:/Johnny/Roswell_Projects/Sehee_mRNA_piRNA/umi_celegans_consensus/Allele_Proportions"
gtf_file <- "Caenorhabditis_elegans.WBcel235.114.gtf" 
path_22G <- "Y:/Johnny/Roswell_Projects/Sehee_mRNA_piRNA/Github.Sehee.paper/Data/output/all_deseq_tables/norRNA_cond1_sx2499_30min_HS_cond2_N2_30min_HS_deseq_table.csv"

# --- B. Filtering Parameters ---
MIN_DEPTH <- 10         
MAX_VAF_SNP <- 1.0      
MIDDLE_VAF_SNP <- 0.5   

# --- C. Target Classes for 22G Analysis ---
target_classes_22G <- c("CSR Class 22G", "WAGO Class 22G", "unclassifed 22G")

# Load Wormbase Gene IDs
gids <- wbData::wb_load_gene_ids("WS269")

# ==============================================================================
# 3. DATA IMPORT & CLEANING
# ==============================================================================
message("--- Step 3: Loading Consensus Data ---")

process_sample <- function(filepath) {
  sample_name <- basename(filepath) %>% str_remove(".consensus.per_position_results.tsv")
  group <- ifelse(grepl("N2", sample_name), "N2", "PRDE1")
  
  df <- read_tsv(filepath, show_col_types = FALSE)
  
  df %>%
    filter(TotalReads >= MIN_DEPTH) %>%
    filter(Proportion != MAX_VAF_SNP) %>% 
    filter(Proportion != MIDDLE_VAF_SNP) %>%
    mutate(SampleID = sample_name, Group = group)
}

file_pattern <- glob2rx("*.consensus.per_position_results.tsv")
files <- list.files(pattern = file_pattern, full.names = TRUE, path = consensus_path)

if(length(files) == 0) stop("No input files found!")

all_data_list <- lapply(files, process_sample)
combined_data <- bind_rows(all_data_list)
message("Data loaded successfully.")

# ==============================================================================
# 3.5. FREQUENCY TABLE: ALL SITES (0.01 BINS)
# ==============================================================================
message("--- Step 3.5: Calculating Frequency All Sites (0.01 Bins) ---")

get_freq_data_all <- function(df, group_label) {
  d <- df %>% filter(Group == group_label)
  total <- nrow(d)
  
  # 1. Count Exact 0s
  n_zero <- sum(d$Proportion == 0)
  
  # 2. Process Bins > 0 with 0.01 step
  breaks_seq <- seq(0, 1, by = 0.01) 
  
  d_pos <- d %>% filter(Proportion > 0)
  d_pos$Bin <- cut(d_pos$Proportion, breaks = breaks_seq, right = TRUE, include.lowest = FALSE)
  
  counts_pos <- d_pos %>%
    group_by(Bin) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Template
  all_levels <- levels(cut(seq(0.001, 1, by=0.01), breaks = breaks_seq, right = TRUE))
  template <- data.frame(Bin = factor(all_levels, levels = all_levels))
  
  counts_pos_full <- template %>%
    left_join(counts_pos, by = "Bin") %>%
    mutate(Count = replace_na(Count, 0)) %>%
    mutate(Bin = as.character(Bin))
  
  # Combine [0,0] row
  row_zero <- data.frame(Bin = "[0,0]", Count = n_zero)
  final_df <- bind_rows(row_zero, counts_pos_full)
  
  # Percentages (Based on ALL sites)
  final_df$Percentage <- final_df$Count / total
  
  colnames(final_df) <- c("Bin", paste0("Sites_Counts_", group_label), paste0("Percentage_", group_label))
  return(final_df)
}

freq_N2 <- get_freq_data_all(combined_data, "N2")
freq_PRDE1 <- get_freq_data_all(combined_data, "PRDE1")
freq_table <- full_join(freq_N2, freq_PRDE1, by = "Bin")

# Format & Save Word Doc
freq_table_formatted <- freq_table %>%
  mutate(
    `Sites Counts N2` = format(Sites_Counts_N2, big.mark = ",", scientific = FALSE),
    `Percentage N2` = scales::percent(Percentage_N2, accuracy = 0.001),
    `Sites Counts PRDE1` = format(Sites_Counts_PRDE1, big.mark = ",", scientific = FALSE),
    `Percentage PRDE1` = scales::percent(Percentage_PRDE1, accuracy = 0.001)
  ) %>%
  select(Bin, `Sites Counts N2`, `Percentage N2`, `Sites Counts PRDE1`, `Percentage PRDE1`)

ft <- flextable(freq_table_formatted) %>% autofit() %>% theme_vanilla() %>% align(align = "center", part = "all")
save_as_docx(ft, path = file.path(OUTPUT_DIR, "Frequency_Binned_Allele_Proportions.docx"))
write_csv(freq_table_formatted, file.path(OUTPUT_DIR, "Frequency_Binned_Allele_Proportions.csv"))

# ==============================================================================
# 3.6. FREQUENCY TABLE: VARIANT SITES ONLY (0.05 BINS) [NEW]
# ==============================================================================
message("--- Step 3.6: Calculating Frequency Variant Sites Only (0.05 Bins) ---")

get_freq_data_variant <- function(df, group_label) {
  # Filter only variant sites (Proportion > 0)
  d_variant <- df %>% filter(Group == group_label & Proportion > 0)
  total_variant <- nrow(d_variant)
  
  # Define breaks: 0.05 step
  breaks_seq <- seq(0, 1, by = 0.05)
  
  # Binning (include.lowest=TRUE to catch exactly 0 if any slipped through, but we filtered >0)
  # For variants, usually the first bin is (0, 0.05]
  d_variant$Bin <- cut(d_variant$Proportion, breaks = breaks_seq, right = TRUE, include.lowest = FALSE)
  
  counts_pos <- d_variant %>%
    group_by(Bin) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Template for 0.05 bins
  all_levels <- levels(cut(seq(0.01, 1, by=0.05), breaks = breaks_seq, right = TRUE))
  template <- data.frame(Bin = factor(all_levels, levels = all_levels))
  
  counts_full <- template %>%
    left_join(counts_pos, by = "Bin") %>%
    mutate(Count = replace_na(Count, 0)) %>%
    mutate(Bin = as.character(Bin))
  
  # Note: No [0,0] row here because we are looking at Variants Only
  
  # Percentages (Based on VARIANT sites only)
  counts_full$Percentage <- counts_full$Count / total_variant
  
  colnames(counts_full) <- c("Bin", paste0("Sites_Counts_", group_label), paste0("Percentage_", group_label))
  return(counts_full)
}

var_N2 <- get_freq_data_variant(combined_data, "N2")
var_PRDE1 <- get_freq_data_variant(combined_data, "PRDE1")
var_table <- full_join(var_N2, var_PRDE1, by = "Bin")

# Format & Save Variant Word Doc
var_table_formatted <- var_table %>%
  mutate(
    `Sites Counts N2` = format(Sites_Counts_N2, big.mark = ",", scientific = FALSE),
    `Percentage N2` = scales::percent(Percentage_N2, accuracy = 0.01),
    `Sites Counts PRDE1` = format(Sites_Counts_PRDE1, big.mark = ",", scientific = FALSE),
    `Percentage PRDE1` = scales::percent(Percentage_PRDE1, accuracy = 0.01)
  ) %>%
  select(Bin, `Sites Counts N2`, `Percentage N2`, `Sites Counts PRDE1`, `Percentage PRDE1`)

ft_var <- flextable(var_table_formatted) %>% autofit() %>% theme_vanilla() %>% align(align = "center", part = "all") %>% 
  set_caption("Frequency of Binned Allele Proportions (Variant Sites Only)")

save_as_docx(ft_var, path = file.path(OUTPUT_DIR, "Frequency_Binned_Allele_Proportions_Variant_Only.docx"))
write_csv(var_table_formatted, file.path(OUTPUT_DIR, "Frequency_Binned_Allele_Proportions_Variant_Only.csv"))
message("Saved Variant-Only Frequency Table.")

# ==============================================================================
# 4. GLOBAL ERROR RATES
# ==============================================================================
message("--- Step 4: Global Error Rate Analysis ---")

global_stats <- combined_data %>%
  group_by(SampleID, Group) %>%
  summarise(
    Total_Errors = sum(AltReads),
    Total_Bases = sum(TotalReads),
    Error_Rate = Total_Errors / Total_Bases,
    .groups = "drop"
  )

write_csv(global_stats, file.path(OUTPUT_DIR, "Global_Error_Stats.csv"))

# Stats Tests
t_test_result <- t.test(Error_Rate ~ Group, data = global_stats)
t_test_df <- data.frame(
  statistic = t_test_result$statistic, p_value = t_test_result$p.value,
  mean_group1 = t_test_result$estimate[1], mean_group2 = t_test_result$estimate[2]
)
write_csv(t_test_df, file.path(OUTPUT_DIR, "Global_Error_Ttest_Result.csv"))

ks_result <- ks.test(combined_data$Proportion[combined_data$Group == "N2"], 
                     combined_data$Proportion[combined_data$Group == "PRDE1"])
write_csv(data.frame(stat=ks_result$statistic, p=ks_result$p.value), 
          file.path(OUTPUT_DIR, "Global_Error_KS_Test_Result.csv"))

# Plot
p1 <- ggplot(global_stats, aes(x = Group, y = Error_Rate, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, show.legend = F) +
  scale_y_continuous(labels = scales::scientific) +
  labs(y = "Error Rate (Errors/bp)", 
       subtitle = paste("KS p-val =", format.pval(ks_result$p.value, digits=3))) +
  theme_bw(base_size = 14) + theme(axis.title.x = element_blank())

ggsave(file.path(OUTPUT_DIR, "Global_Error_Rate.png"), p1, width = 5, height = 4)

# ==============================================================================
# 5. GENE-LEVEL ANALYSIS
# ==============================================================================
message("--- Step 5: Gene-Specific Error Rates ---")

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
genes_gr <- genes(txdb)

gr_errors <- makeGRangesFromDataFrame(combined_data, 
                                      seqnames.field = "CHROM", start.field = "POS", end.field = "POS", 
                                      keep.extra.columns = TRUE)

overlaps <- findOverlaps(gr_errors, genes_gr)
annotated_data <- combined_data[queryHits(overlaps), ]
annotated_data$GeneID <- names(genes_gr)[subjectHits(overlaps)]

gene_results <- annotated_data %>%
  group_by(GeneID, Group) %>%
  summarise(Gene_AltReads = sum(AltReads), Gene_TotalReads = sum(TotalReads), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = c(Gene_AltReads, Gene_TotalReads), values_fill = 0) %>%
  filter(Gene_TotalReads_N2 > 1000 & Gene_TotalReads_PRDE1 > 1000) %>%
  mutate(
    Rate_N2 = Gene_AltReads_N2 / Gene_TotalReads_N2,
    Rate_PRDE1 = Gene_AltReads_PRDE1 / Gene_TotalReads_PRDE1,
    Log2FC = log2(Rate_PRDE1 / Rate_N2)
  ) %>%
  mutate(GeneSymbol = wbData::i2s(GeneID, gids)) %>%
  relocate(GeneSymbol, .after = GeneID) %>% 
  arrange(desc(Log2FC))

write_csv(gene_results, file.path(OUTPUT_DIR, "Gene_Specific_Error_Rates.csv"))

genes_TE_UP <- gene_results %>% filter(Log2FC > 0) %>% pull(GeneID)
genes_TE_DOWN <- gene_results %>% filter(Log2FC < 0) %>% pull(GeneID)

# ==============================================================================
# 6. GO ENRICHMENT
# ==============================================================================
message("--- Step 6: GO Enrichment ---")

convert_go_ids <- function(ids) {
  sapply(ids, function(x) {
    split_ids <- unlist(strsplit(x, "/"))
    paste(wbData::i2s(split_ids, gids), collapse = "/")
  })
}

run_go <- function(gene_list, name, title) {
  if(length(gene_list) < 5) return(NULL)
  go <- enrichGO(gene = gene_list, OrgDb = org.Ce.eg.db, keyType = "WORMBASE", ont = "BP", pvalueCutoff = 0.05)
  if (!is.null(go) && nrow(go) > 0) {
    go@result <- mutate(go@result, FoldEnrichment = parse_ratio(GeneRatio)/parse_ratio(BgRatio))
    res_df <- filter(go@result, p.adjust < 0.05)
    if(nrow(res_df) > 0) {
      res_df$geneSymbol <- convert_go_ids(res_df$geneID)
      write_csv(res_df, file.path(OUTPUT_DIR, paste0("GO_Results_", name, ".csv")))
      p <- dotplot(go, x = "FoldEnrichment", showCategory = 15) + ggtitle(title)
      ggsave(file.path(OUTPUT_DIR, paste0("GO_Enrichment_", name, ".png")), p, width = 8, height = 6)
    }
  }
}

run_go(genes_TE_UP, "UP_Errors", "Higher Errors in PRDE1")
run_go(genes_TE_DOWN, "DOWN_Errors", "Lower Errors in PRDE1")

# ==============================================================================
# 7. 22G RNA OVERLAP (FIXED PLOTTING)
# ==============================================================================
message("--- Step 7: 22G RNA Overlap Analysis ---")

df_22G <- read_csv(path_22G, show_col_types = FALSE)
if(!"gene_id" %in% colnames(df_22G) && "baseMean" %in% colnames(df_22G)) {
  df_22G <- df_22G %>% rename(gene_id = 10) # Keeping your logic for column 10
}

genes_22G_UP <- df_22G %>% filter(Classifier %in% target_classes_22G, padj < 0.05, log2FoldChange > 0) %>% pull(gene_id)
genes_22G_DOWN <- df_22G %>% filter(Classifier %in% target_classes_22G, padj < 0.05, log2FoldChange < 0) %>% pull(gene_id)

plot_euler <- function(list_data, filename, plot_title) {
  fit <- euler(list_data)
  
  png(file.path(OUTPUT_DIR, paste0(filename, ".png")), width = 600, height = 600)
  p <- plot(fit, quantities = TRUE, main = plot_title)
  print(p) # ESSENTIAL for Euler plots inside functions
  dev.off()
  
  svg(file.path(OUTPUT_DIR, paste0(filename, ".svg")), width = 6, height = 6)
  p <- plot(fit, quantities = TRUE, main = plot_title)
  print(p) # ESSENTIAL
  dev.off()
}

plot_euler(list("Higher Errors (PRDE1)" = unique(genes_TE_UP), "Higher 22G (PRDE1)" = unique(genes_22G_UP)),
           "Euler_Overlap_UP_Errors_UP_22G", "Overlap: High Errors vs High 22G")

plot_euler(list("Lower Errors (PRDE1)" = unique(genes_TE_DOWN), "Lower 22G (PRDE1)" = unique(genes_22G_DOWN)),
           "Euler_Overlap_DOWN_Errors_DOWN_22G", "Overlap: Low Errors vs Low 22G")

plot_euler(list("Higher Errors (PRDE1)" = unique(genes_TE_UP), "Lower 22G (PRDE1)" = unique(genes_22G_DOWN)),
           "Euler_Overlap_UP_Errors_DOWN_22G", "Overlap: High Errors vs Low 22G")

plot_euler(list("Lower Errors (PRDE1)" = unique(genes_TE_DOWN), "Higher 22G (PRDE1)" = unique(genes_22G_UP)),
           "Euler_Overlap_DOWN_Errors_UP_22G", "Overlap: Low Errors vs High 22G")

message("--- Analysis Complete! Check Analysis_Results_Combined folder. ---")
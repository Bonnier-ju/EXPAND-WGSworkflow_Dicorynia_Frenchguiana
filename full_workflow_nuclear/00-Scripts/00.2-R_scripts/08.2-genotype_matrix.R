#!/usr/bin/env Rscript
# 08.2-genotype_matrix.R
# Post-process plink2 --export A output into clean genotype matrices and
# sample metadata for downstream GEA analyses.
#
# Inputs:
#   *.raw          — plink2 --export A output (individuals x SNPs, 0/1/2/NA)
#                    columns: #FID IID PAT MAT SEX PHENOTYPE <SNP_ID_1> ...
#   sample_sheet_complete_filtered.csv — individual-to-site mapping
#
# Outputs (out_dir/):
#   genotype_012.tsv          — clean matrix: rows=individuals, cols=SNPs (0/1/2)
#   genotype_012_imputed.tsv  — same with NA imputed to column mean (for RDA/LFMM)
#   sample_metadata.tsv       — individual_id, site, lat, long (for env join)
#   snp_ids.txt               — ordered SNP identifiers (CHR:POS:REF:ALT)
#   summary_stats.tsv         — per-SNP missingness, MAF, observed heterozygosity

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args         <- commandArgs(trailingOnly = TRUE)
raw_file     <- if (length(args) >= 1) args[1] else stop("raw_file required")
sample_sheet <- if (length(args) >= 2) args[2] else stop("sample_sheet required")
out_dir      <- if (length(args) >= 3) args[3] else stop("out_dir required")

plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# STEP 1: Load raw plink2 --export A output
# -------------------------------------------------------------------
cat("Loading plink2 raw matrix:", raw_file, "\n")
raw <- read.table(raw_file, header = TRUE, sep = "\t",
                  check.names = FALSE, stringsAsFactors = FALSE)

# plink2 --export A columns: #FID IID PAT MAT SEX PHENOTYPE <snp1> ...
meta_cols <- c("#FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
snp_cols  <- setdiff(names(raw), meta_cols)

cat(sprintf("INFO: %d individuals x %d SNPs loaded\n", nrow(raw), length(snp_cols)))

# Convert to numeric matrix
geno_mat <- as.matrix(raw[, snp_cols])
mode(geno_mat) <- "numeric"
rownames(geno_mat) <- raw[["IID"]]

# Strip plink allele suffix from column names: "CHR:POS:REF:ALT_ALT" -> "CHR:POS:REF:ALT"
colnames(geno_mat) <- sub("_[^_]+$", "", colnames(geno_mat))

# -------------------------------------------------------------------
# STEP 2: Per-SNP summary statistics
# -------------------------------------------------------------------
cat("Computing per-SNP statistics...\n")
n_ind       <- nrow(geno_mat)
n_missing   <- colSums(is.na(geno_mat))
miss_rate   <- n_missing / n_ind
alt_count   <- colSums(geno_mat, na.rm = TRUE)
n_called    <- n_ind - n_missing
maf_vec     <- alt_count / (2 * n_called)
maf_vec     <- pmin(maf_vec, 1 - maf_vec)
n_het       <- colSums(geno_mat == 1, na.rm = TRUE)
obs_het     <- n_het / n_called

summ_df <- data.frame(
  snp_id    = colnames(geno_mat),
  n_missing = n_missing,
  miss_rate = round(miss_rate, 4),
  maf       = round(maf_vec, 4),
  obs_het   = round(obs_het, 4),
  stringsAsFactors = FALSE
)
write.table(summ_df, file = file.path(out_dir, "snp_summary_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("INFO: mean missing rate: %.2f%%\n", mean(miss_rate) * 100))
cat(sprintf("INFO: mean MAF: %.3f\n", mean(maf_vec)))
cat(sprintf("INFO: SNPs with missing data: %d\n", sum(n_missing > 0)))

# -------------------------------------------------------------------
# STEP 3: Write raw 012 matrix (with NAs)
# -------------------------------------------------------------------
out_raw <- file.path(out_dir, "genotype_012.tsv")
write.table(cbind(individual_id = rownames(geno_mat), as.data.frame(geno_mat)),
            file = out_raw, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Raw 012 matrix written:", out_raw, "\n")

# -------------------------------------------------------------------
# STEP 4: Impute missing genotypes to column mean (integer round)
# Mean imputation is the standard approach for RDA and LFMM2.
# -------------------------------------------------------------------
cat("Imputing missing genotypes to column mean...\n")
geno_imp <- geno_mat
col_means <- colMeans(geno_mat, na.rm = TRUE)
for (j in seq_len(ncol(geno_imp))) {
  na_idx <- is.na(geno_imp[, j])
  if (any(na_idx)) geno_imp[na_idx, j] <- col_means[j]
}

out_imp <- file.path(out_dir, "genotype_012_imputed.tsv")
write.table(cbind(individual_id = rownames(geno_imp), as.data.frame(geno_imp)),
            file = out_imp, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Imputed 012 matrix written:", out_imp, "\n")

# -------------------------------------------------------------------
# STEP 5: Write ordered SNP ID list
# -------------------------------------------------------------------
writeLines(colnames(geno_mat), file.path(out_dir, "snp_ids.txt"))
cat(sprintf("SNP ID list written: %d SNPs\n", length(colnames(geno_mat))))

# -------------------------------------------------------------------
# STEP 6: Sample metadata (individual -> site -> lat/long)
# -------------------------------------------------------------------
cat("Building sample metadata from:", sample_sheet, "\n")
ss <- read.csv(sample_sheet, stringsAsFactors = FALSE)

# Keep only FG individuals present in the genotype matrix
ss_gea <- ss[ss$sample_id %in% rownames(geno_mat), ]

# Load geographic coordinates
GEOLOC_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/geoloc_site.csv"
if (file.exists(GEOLOC_FILE)) {
  geoloc <- read.csv(GEOLOC_FILE, stringsAsFactors = FALSE)
  names(geoloc)[names(geoloc) == "Sites"] <- "site"
  ss_gea <- merge(ss_gea[, c("sample_id", "site")], geoloc, by = "site", all.x = TRUE)
} else {
  warning("geoloc_site.csv not found — lat/long will be NA")
  ss_gea <- ss_gea[, c("sample_id", "site")]
  ss_gea$lat  <- NA_real_
  ss_gea$long <- NA_real_
}

# Reorder to match genotype matrix row order
ss_gea <- ss_gea[match(rownames(geno_mat), ss_gea$sample_id), ]

out_meta <- file.path(out_dir, "sample_metadata.tsv")
write.table(ss_gea, file = out_meta, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Sample metadata written:", out_meta, "\n")
cat(sprintf("INFO: %d individuals, %d sites\n",
            nrow(ss_gea), length(unique(ss_gea$site))))
print(table(ss_gea$site))

# -------------------------------------------------------------------
# STEP 7: QC plots
# -------------------------------------------------------------------
# MAF distribution
p_maf <- ggplot(summ_df, aes(x = maf)) +
  geom_histogram(bins = 50, fill = "#4575b4", color = "white", linewidth = 0.2) +
  labs(title = "MAF distribution - GEA SNP set",
       subtitle = sprintf("%d SNPs, 155 FG individuals", nrow(summ_df)),
       x = "Minor allele frequency", y = "Count") +
  geom_vline(xintercept = 0.05, color = "#d73027", linetype = "dashed", linewidth = 0.7)
ggsave(file.path(plot_dir, "snp_maf_distribution.png"), p_maf,
       width = 7, height = 4, dpi = 300)

# Missingness distribution
p_miss <- ggplot(summ_df[summ_df$miss_rate > 0, ], aes(x = miss_rate)) +
  geom_histogram(bins = 40, fill = "#d73027", color = "white", linewidth = 0.2) +
  labs(title = "Missing data rate per SNP",
       x = "Missing rate", y = "Count")
ggsave(file.path(plot_dir, "snp_missingness.png"), p_miss,
       width = 7, height = 4, dpi = 300)

# Individuals per site barplot
site_counts <- as.data.frame(table(site = ss_gea$site))
p_sites <- ggplot(site_counts, aes(x = reorder(site, -Freq), y = Freq)) +
  geom_col(fill = "#4575b4", width = 0.7) +
  geom_text(aes(label = Freq), vjust = -0.4, size = 3) +
  labs(title = "Individuals per site - GEA dataset (155 FG)",
       x = NULL, y = "N individuals") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
ggsave(file.path(plot_dir, "individuals_per_site.png"), p_sites,
       width = 10, height = 4, dpi = 300)

cat("QC plots written to:", plot_dir, "\n")

# -------------------------------------------------------------------
# Final summary
# -------------------------------------------------------------------
cat(sprintf("\nDONE 08.2 genotype matrix preparation completed\n"))
cat(sprintf("  Individuals           : %d\n", nrow(geno_mat)))
cat(sprintf("  SNPs                  : %d\n", ncol(geno_mat)))
cat(sprintf("  SNPs with any NA      : %d (%.1f%%)\n",
            sum(n_missing > 0), mean(n_missing > 0) * 100))
cat(sprintf("  Mean MAF              : %.3f\n", mean(maf_vec)))
cat(sprintf("  genotype_012.tsv      : %s\n", out_raw))
cat(sprintf("  genotype_012_imputed  : %s\n", out_imp))
cat(sprintf("  sample_metadata.tsv   : %s\n", out_meta))

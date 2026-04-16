#!/usr/bin/env Rscript
# 08.0.1-chelsa_extraction.R
# Post-process CHELSA values extracted by gdallocationinfo (SLURM step 3).
# No terra dependency — works with R 4.3.0.
#
# Inputs:
#   chelsa_raw_extracted.tsv  — raw integer values from gdallocationinfo
# Outputs:
#   chelsa_env_per_site.csv   — scaled to physical units
#   chelsa_env_scaled.csv     — z-score standardised (for LFMM/RDA)
#   chelsa_env_summary.tsv    — descriptive stats per variable
#   chelsa_high_correlations.tsv — pairs with |r| >= 0.8
#   plots/chelsa_bio_correlation.png/pdf

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
extracted_tsv <- if (length(args) >= 1) args[1] else stop("extracted_tsv required")
out_dir       <- if (length(args) >= 2) args[2] else stop("out_dir required")

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# STEP 1: Load extracted values
# -------------------------------------------------------------------
cat("Loading extracted CHELSA values:", extracted_tsv, "\n")
raw <- read.table(extracted_tsv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)

bio_cols <- paste0("BIO", 1:19)
for (col in bio_cols) {
  raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
}

cat(sprintf("INFO: %d sites x %d BIO variables loaded\n", nrow(raw), length(bio_cols)))

# -------------------------------------------------------------------
# STEP 2: Apply CHELSA v2.1 scale factors
# -------------------------------------------------------------------
# All UInt16 variables: actual = raw × 0.1 (+ offset where applicable)
# confirmed via gdalinfo on the downloaded GeoTIFFs.
#
# Absolute temperature (UInt16, Scale=0.1, Offset=-273.15): °C = raw×0.1 - 273.15
temp_abs_vars <- c("BIO1","BIO5","BIO6","BIO8","BIO9","BIO10","BIO11")
for (v in temp_abs_vars) raw[[v]] <- raw[[v]] / 10 - 273.15

# Temperature range/difference (UInt16, Scale=0.1, Offset=0): °C = raw×0.1
temp_diff_vars <- c("BIO2","BIO7")
for (v in temp_diff_vars) raw[[v]] <- raw[[v]] / 10

# BIO4 (temperature seasonality, SD×100): UInt16, Scale=0.1, Offset=0
raw[["BIO4"]] <- raw[["BIO4"]] / 10

# BIO3 (isothermality, %): Float32 — gdallocationinfo applies Scale=0.1 automatically,
# so the extracted value = actual_% × 0.1; multiply back by 10 to recover true %
raw[["BIO3"]] <- raw[["BIO3"]] * 10

# BIO12-BIO19 (precipitation, mm): UInt16, Scale=0.1, Offset=0
precip_vars <- paste0("BIO", 12:19)
for (v in precip_vars) raw[[v]] <- raw[[v]] / 10

env_df <- raw
colnames(env_df)[colnames(env_df) == "site"] <- "site"

out_raw <- file.path(out_dir, "chelsa_env_per_site.csv")
write.csv(env_df, file = out_raw, row.names = FALSE, quote = FALSE)
cat("Physical units values written:", out_raw, "\n")

# -------------------------------------------------------------------
# STEP 3: Z-score standardisation for LFMM/RDA
# -------------------------------------------------------------------
env_scaled <- env_df
env_scaled[bio_cols] <- as.data.frame(scale(env_df[bio_cols]))

out_scaled <- file.path(out_dir, "chelsa_env_scaled.csv")
write.csv(env_scaled, file = out_scaled, row.names = FALSE, quote = FALSE)
cat("Z-score scaled values written:", out_scaled, "\n")

# -------------------------------------------------------------------
# STEP 4: Summary statistics
# -------------------------------------------------------------------
summary_df <- data.frame(
  variable = bio_cols,
  mean     = round(colMeans(env_df[bio_cols], na.rm = TRUE), 3),
  sd       = round(apply(env_df[bio_cols], 2, sd, na.rm = TRUE), 3),
  min      = round(apply(env_df[bio_cols], 2, min, na.rm = TRUE), 3),
  max      = round(apply(env_df[bio_cols], 2, max, na.rm = TRUE), 3),
  n_sites  = colSums(!is.na(env_df[bio_cols])),
  stringsAsFactors = FALSE
)
out_summ <- file.path(out_dir, "chelsa_env_summary.tsv")
write.table(summary_df, file = out_summ, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary statistics written:", out_summ, "\n")
print(summary_df)

# -------------------------------------------------------------------
# STEP 5: Correlation matrix heatmap
# -------------------------------------------------------------------
cor_mat  <- cor(env_df[bio_cols], use = "pairwise.complete.obs")
cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Var1", "Var2", "Correlation")
cor_long$Var1 <- factor(cor_long$Var1, levels = bio_cols)
cor_long$Var2 <- factor(cor_long$Var2, levels = bio_cols)

p_cor <- ggplot(cor_long, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 1.9, color = "black") +
  scale_fill_gradientn(
    colours = c("#d73027", "#ffffff", "#4575b4"),
    limits  = c(-1, 1), name = "Pearson r"
  ) +
  labs(title = "CHELSA BIO1-BIO19 — Pearson correlation matrix",
       x = NULL, y = NULL) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    plot.title       = element_text(face = "bold"),
    legend.key.height = unit(0.8, "cm")
  )

ggsave(file.path(plot_dir, "chelsa_bio_correlation.png"), p_cor,
       width = 10, height = 9, dpi = 300)
ggsave(file.path(plot_dir, "chelsa_bio_correlation.pdf"), p_cor,
       width = 10, height = 9)
cat("Correlation heatmap written\n")

# -------------------------------------------------------------------
# STEP 6: Identify highly correlated pairs (|r| >= 0.8)
# -------------------------------------------------------------------
high_cor <- which(abs(cor_mat) >= 0.8 & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  hc_df <- data.frame(
    Var1 = rownames(cor_mat)[high_cor[, 1]],
    Var2 = colnames(cor_mat)[high_cor[, 2]],
    r    = round(cor_mat[high_cor], 3),
    stringsAsFactors = FALSE
  )
  hc_df <- hc_df[order(-abs(hc_df$r)), ]
  out_hc <- file.path(out_dir, "chelsa_high_correlations.tsv")
  write.table(hc_df, file = out_hc, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("INFO: %d highly correlated pairs (|r| >= 0.8): %s\n",
              nrow(hc_df), out_hc))
  print(hc_df)
} else {
  cat("INFO: no highly correlated pairs (|r| >= 0.8) detected\n")
}

# -------------------------------------------------------------------
# STEP 7: Barplots per BIO variable across sites
# -------------------------------------------------------------------
for (v in bio_cols) {
  df_v <- env_df[!is.na(env_df[[v]]), c("site", v)]
  df_v <- df_v[order(df_v[[v]]), ]
  df_v$site <- factor(df_v$site, levels = df_v$site)

  p <- ggplot(df_v, aes_string(x = "site", y = v)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_text(aes_string(label = sprintf('sprintf("%%.2f", %s)', v)),
              vjust = -0.4, size = 2.4) +
    labs(title = sprintf("CHELSA %s per site", v), x = NULL, y = v) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.title  = element_text(face = "bold"))

  fn <- file.path(plot_dir, sprintf("chelsa_%s_per_site.png", tolower(v)))
  ggsave(fn, p, width = max(8, nrow(df_v) * 0.45), height = 4, dpi = 300)
}
cat("Per-site barplots written\n")

cat("\nDONE 08.0.1 CHELSA post-processing completed\n")
cat("Outputs:\n")
cat("  Physical units :", out_raw,    "\n")
cat("  Scaled         :", out_scaled, "\n")
cat("  Summary        :", out_summ,   "\n")

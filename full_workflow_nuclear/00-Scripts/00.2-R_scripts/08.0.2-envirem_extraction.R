#!/usr/bin/env Rscript
# 08.0.2-envirem_extraction.R
# Post-process ENVIREM + elevation values extracted by gdallocationinfo (SLURM step 4).
#
# Inputs:
#   envirem_raw_extracted.tsv  — raw values from gdallocationinfo
#                                (site, lat, long, 18 ENVIREM vars, elevation)
# Outputs:
#   data/envirem/envirem_env_per_site.csv  — physical units
#   data/envirem/envirem_env_scaled.csv    — z-score standardised (for LFMM/RDA)
#   data/envirem/envirem_env_summary.tsv   — descriptive stats per variable
#   data/envirem/envirem_high_correlations.tsv — pairs with |r| >= 0.7
#   plots/envirem_correlation.png/pdf
#   plots/envirem_<var>_per_site.png
#
# Scale factors applied:
#   maxTempColdest, minTempWarmest : Int16 raw → ÷ 10 → °C
#   elevation                      : Float32, metres, no scaling
#   All other ENVIREM variables    : Float32, physical units, no scaling
#
# Variable sources:
#   ENVIREM (18 vars) — WorldClim v1.4 / ENVIREM v2 (Title & Bemmels 2019, Ecography)
#   elevation         — Copernicus DEM GLO-30, ~30 m resolution

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# Site color palette
SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else {
  warning("Site color file not found: ", SITE_COLORS_FILE)
  c()
}

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args          <- commandArgs(trailingOnly = TRUE)
extracted_tsv <- if (length(args) >= 1) args[1] else stop("extracted_tsv required")
out_dir       <- if (length(args) >= 2) args[2] else stop("out_dir required")

env_out_dir <- file.path(out_dir, "data/envirem")
plot_dir    <- file.path(out_dir, "plots")
dir.create(env_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,    recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# STEP 1: Load raw extracted values
# -------------------------------------------------------------------
cat("Loading raw ENVIREM values:", extracted_tsv, "\n")
raw <- read.table(extracted_tsv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)

env_cols <- c("annualPET", "aridityIndexThornthwaite", "climaticMoistureIndex",
              "continentality", "embergerQ", "growingDegDays0", "growingDegDays5",
              "maxTempColdest", "minTempWarmest", "monthCountByTemp10",
              "PETColdestQuarter", "PETDriestQuarter", "PETseasonality",
              "PETWarmestQuarter", "PETWettestQuarter", "thermicityIndex",
              "topoWet", "tri", "elevation")

for (col in env_cols) {
  raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
}
cat(sprintf("INFO: %d sites × %d variables loaded\n", nrow(raw), length(env_cols)))

# -------------------------------------------------------------------
# STEP 2: Apply scale factors
# -------------------------------------------------------------------
# maxTempColdest and minTempWarmest: Int16 raw value × 0.1 → °C
raw[["maxTempColdest"]] <- raw[["maxTempColdest"]] / 10
raw[["minTempWarmest"]] <- raw[["minTempWarmest"]] / 10

# All others: Float32, values already in physical units — no scaling needed.
# elevation: Copernicus DEM GLO-30, Float32, metres above EGM2008 geoid.

env_df <- raw

# Sanity checks
cat(sprintf("INFO: elevation range: %.0f – %.0f m\n",
            min(env_df$elevation, na.rm = TRUE),
            max(env_df$elevation, na.rm = TRUE)))
cat(sprintf("INFO: maxTempColdest range: %.1f – %.1f °C\n",
            min(env_df$maxTempColdest, na.rm = TRUE),
            max(env_df$maxTempColdest, na.rm = TRUE)))
cat(sprintf("INFO: annualPET range: %.0f – %.0f mm/yr\n",
            min(env_df$annualPET, na.rm = TRUE),
            max(env_df$annualPET, na.rm = TRUE)))

# -------------------------------------------------------------------
# STEP 3: Write physical-unit output
# -------------------------------------------------------------------
out_raw <- file.path(env_out_dir, "envirem_env_per_site.csv")
write.csv(env_df, file = out_raw, row.names = FALSE, quote = FALSE)
cat("Physical units written:", out_raw, "\n")

# -------------------------------------------------------------------
# STEP 4: Z-score standardisation for LFMM/RDA
# -------------------------------------------------------------------
env_scaled <- env_df
env_scaled[env_cols] <- as.data.frame(scale(env_df[env_cols]))

out_scaled <- file.path(env_out_dir, "envirem_env_scaled.csv")
write.csv(env_scaled, file = out_scaled, row.names = FALSE, quote = FALSE)
cat("Z-score scaled values written:", out_scaled, "\n")

# -------------------------------------------------------------------
# STEP 5: Summary statistics
# -------------------------------------------------------------------
summary_df <- data.frame(
  variable = env_cols,
  mean     = round(colMeans(env_df[env_cols], na.rm = TRUE), 3),
  sd       = round(apply(env_df[env_cols], 2, sd,  na.rm = TRUE), 3),
  min      = round(apply(env_df[env_cols], 2, min, na.rm = TRUE), 3),
  max      = round(apply(env_df[env_cols], 2, max, na.rm = TRUE), 3),
  n_sites  = colSums(!is.na(env_df[env_cols])),
  stringsAsFactors = FALSE
)
out_summ <- file.path(env_out_dir, "envirem_env_summary.tsv")
write.table(summary_df, file = out_summ, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary statistics written:", out_summ, "\n")
print(summary_df)

# -------------------------------------------------------------------
# STEP 6: Correlation matrix heatmap
# -------------------------------------------------------------------
cor_mat  <- cor(env_df[env_cols], use = "pairwise.complete.obs")
cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Var1", "Var2", "Correlation")
cor_long$Var1 <- factor(cor_long$Var1, levels = env_cols)
cor_long$Var2 <- factor(cor_long$Var2, levels = env_cols)

p_cor <- ggplot(cor_long, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 1.6, color = "black") +
  scale_fill_gradientn(
    colours = c("#d73027", "#ffffff", "#4575b4"),
    limits  = c(-1, 1), name = "Pearson r"
  ) +
  labs(title = "ENVIREM + elevation — Pearson correlation matrix (site level)",
       x = NULL, y = NULL) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 7.5),
    axis.text.y       = element_text(size = 7.5),
    plot.title        = element_text(face = "bold"),
    legend.key.height = unit(0.8, "cm")
  )

ggsave(file.path(plot_dir, "envirem_correlation.png"), p_cor,
       width = 12, height = 11, dpi = 300)
ggsave(file.path(plot_dir, "envirem_correlation.pdf"), p_cor,
       width = 12, height = 11)
cat("Correlation heatmap written\n")

# -------------------------------------------------------------------
# STEP 7: Identify highly correlated pairs (|r| >= 0.7)
# -------------------------------------------------------------------
# Threshold 0.7 following Legendre & Legendre (2012) and Zuur et al. (2010)
# as applied in the variable selection procedure (RDA + correlation filtering).
high_cor <- which(abs(cor_mat) >= 0.7 & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  hc_df <- data.frame(
    Var1 = rownames(cor_mat)[high_cor[, 1]],
    Var2 = colnames(cor_mat)[high_cor[, 2]],
    r    = round(cor_mat[high_cor], 3),
    stringsAsFactors = FALSE
  )
  hc_df <- hc_df[order(-abs(hc_df$r)), ]
  out_hc <- file.path(env_out_dir, "envirem_high_correlations.tsv")
  write.table(hc_df, file = out_hc, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("INFO: %d pairs with |r| >= 0.7: %s\n", nrow(hc_df), out_hc))
  print(hc_df)
} else {
  cat("INFO: no pairs with |r| >= 0.7 detected\n")
}

# -------------------------------------------------------------------
# STEP 8: Barplots per variable across sites
# -------------------------------------------------------------------
var_labels <- c(
  annualPET                = "Annual PET (mm/yr)",
  aridityIndexThornthwaite = "Thornthwaite aridity index",
  climaticMoistureIndex    = "Climatic moisture index",
  continentality           = "Conrad's continentality index",
  embergerQ                = "Emberger's pluviometric quotient",
  growingDegDays0          = "Growing degree days (base 0\u00b0C)",
  growingDegDays5          = "Growing degree days (base 5\u00b0C)",
  maxTempColdest           = "Max. temp. coldest month (\u00b0C)",
  minTempWarmest           = "Min. temp. warmest month (\u00b0C)",
  monthCountByTemp10       = "Months with mean temp > 10\u00b0C",
  PETColdestQuarter        = "PET coldest quarter (mm)",
  PETDriestQuarter         = "PET driest quarter (mm)",
  PETseasonality           = "PET seasonality (SD\u00d7100)",
  PETWarmestQuarter        = "PET warmest quarter (mm)",
  PETWettestQuarter        = "PET wettest quarter (mm)",
  thermicityIndex          = "Rivas-Mart\u00ednez thermicity index",
  topoWet                  = "Topographic wetness index",
  tri                      = "Terrain ruggedness index",
  elevation                = "Elevation (m, Copernicus DEM GLO-30)"
)

for (v in env_cols) {
  df_v <- env_df[!is.na(env_df[[v]]), c("site", v)]
  if (nrow(df_v) == 0) next
  df_v <- df_v[order(df_v[[v]]), ]
  df_v$site <- factor(df_v$site, levels = df_v$site)
  label_v <- if (v %in% names(var_labels)) var_labels[[v]] else v

  p <- ggplot(df_v, aes(x = site, y = .data[[v]], fill = site)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.3g", .data[[v]])),
              vjust = ifelse(df_v[[v]] >= 0, -0.4, 1.4), size = 2.4) +
    labs(title = label_v, x = NULL, y = label_v) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      plot.title      = element_text(face = "bold"),
      legend.position = "none"
    )

  if (length(site_colors) > 0)
    p <- p + scale_fill_manual(values = site_colors, na.value = "grey50")

  fn <- file.path(plot_dir, sprintf("envirem_%s_per_site.png", tolower(v)))
  ggsave(fn, p, width = max(8, nrow(df_v) * 0.45), height = 4, dpi = 300)
}
cat("Per-site barplots written\n")

cat("\nDONE 08.0.2 ENVIREM + elevation post-processing completed\n")
cat("Outputs:\n")
cat("  Physical units :", out_raw,    "\n")
cat("  Scaled         :", out_scaled, "\n")
cat("  Summary        :", out_summ,   "\n")

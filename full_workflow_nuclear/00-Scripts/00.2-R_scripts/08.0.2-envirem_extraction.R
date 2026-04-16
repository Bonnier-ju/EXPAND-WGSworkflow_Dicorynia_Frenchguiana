#!/usr/bin/env Rscript
# 08.0.2-envirem_extraction.R
# Compute ENVIREM bioclimatic variables from CHELSA v2.1 monthly data.
# Follows the ENVIREM variable definitions (Title et al. 2019, Ecography).
#
# Inputs:
#   envirem_raw_extracted.tsv  — raw CHELSA monthly values from gdallocationinfo
#                                columns: site, lat, long, tas_01..12, pr_01..12, pet_01..12
#   chelsa_env_per_site.csv    — CHELSA BIO1-19 in physical units (from 08.0.1)
#
# Outputs:
#   envirem_env_per_site.csv   — 14 ENVIREM variables in physical units
#   envirem_env_scaled.csv     — z-score standardised (for LFMM/RDA)
#   envirem_env_summary.tsv    — descriptive stats per variable
#   envirem_high_correlations.tsv — pairs with |r| >= 0.8
#   plots/envirem_correlation.png/pdf
#   plots/envirem_<var>_per_site.png
#
# ENVIREM variables computed:
#   From monthly tas + pr + pet:
#     annualPET            — annual sum of monthly PET (mm/year)
#     PETseasonality       — CV of monthly PET (SD/mean × 100, %)
#     PETColdestQuarter    — sum PET of coldest 3-month window (mm)
#     PETWarmestQuarter    — sum PET of warmest 3-month window (mm)
#     PETDriestQuarter     — sum PET of driest 3-month window (mm)
#     PETWettestQuarter    — sum PET of wettest 3-month window (mm)
#     aridityIndexThornthwaite — annualPET / annualPrecip (> 1 = arid)
#     climaticMoistureIndex — (P - PET) / (P + PET), range [-1, 1]
#     growingDegDays0      — GDD base 0°C (°C·days/year)
#     growingDegDays5      — GDD base 5°C (°C·days/year)
#     monthCountByTemp10   — number of months with mean temp > 10°C (0-12)
#   From CHELSA BIO variables (no monthly data needed):
#     continentality       — Conrad's continentality index (uses BIO7, latitude)
#     embergerQ            — Emberger's pluviometric quotient (uses BIO5, BIO6, BIO12)
#     thermicityIndex      — Rivas-Martínez thermicity index (uses BIO1, BIO5, BIO6)
#
# Reference: Title PQ & Bemmels JB (2019). ENVIREM: An expanded set of
#   bioclimatic and topographic variables increases flexibility and improves
#   performance of ecological niche modeling. Ecography 42:291-307.

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
args           <- commandArgs(trailingOnly = TRUE)
extracted_tsv  <- if (length(args) >= 1) args[1] else stop("extracted_tsv required")
chelsa_csv     <- if (length(args) >= 2) args[2] else stop("chelsa_env_per_site.csv required")
out_dir        <- if (length(args) >= 3) args[3] else stop("out_dir required")

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# STEP 1: Load raw monthly values and apply CHELSA v2.1 scale factors
# -------------------------------------------------------------------
cat("Loading raw monthly values:", extracted_tsv, "\n")
raw <- read.table(extracted_tsv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)

months   <- sprintf("%02d", 1:12)
tas_cols <- paste0("tas_", months)
pr_cols  <- paste0("pr_",  months)
pet_cols <- paste0("pet_", months)

for (col in c(tas_cols, pr_cols, pet_cols)) {
  raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
}

cat(sprintf("INFO: %d sites loaded\n", nrow(raw)))

# --- Scale factors ---
# tas: UInt16, Scale=0.1, Offset=-273.15 → °C = raw × 0.1 - 273.15
# pr:  UInt16, Scale=0.1, Offset=0       → mm/month = raw × 0.1
# pet: UInt16, Scale=0.1, Offset=0       → mm/month = raw × 0.1
# NOTE: if pet values appear 10× too large after running, the pet rasters
#       may be Float32 (scale already applied by gdallocationinfo) →
#       comment out the pet /10 line and set pet_scale to 1.
for (col in tas_cols) raw[[col]] <- raw[[col]] / 10 - 273.15
for (col in pr_cols)  raw[[col]] <- raw[[col]] / 10
for (col in pet_cols) raw[[col]] <- raw[[col]] / 10

# Sanity check: monthly mean temp in French Guiana should be ~24-28°C
tas_range <- range(unlist(raw[, tas_cols]), na.rm = TRUE)
pet_range <- range(unlist(raw[, pet_cols]), na.rm = TRUE)
cat(sprintf("INFO: tas range after scaling: %.1f to %.1f °C\n", tas_range[1], tas_range[2]))
cat(sprintf("INFO: pet range after scaling: %.1f to %.1f mm/month\n", pet_range[1], pet_range[2]))
if (any(abs(tas_range) > 60))
  warning("Monthly temperature out of expected range — check scale factors!")
if (any(pet_range > 500))
  warning("Monthly PET > 500 mm/month — pet scale factor may need adjustment (try dividing by 10 again)!")

# -------------------------------------------------------------------
# STEP 2: Load CHELSA BIO variables (needed for 3 index variables)
# -------------------------------------------------------------------
cat("Loading CHELSA BIO values:", chelsa_csv, "\n")
bio <- read.csv(chelsa_csv, stringsAsFactors = FALSE, check.names = FALSE)

# Merge on site column
df <- merge(raw, bio[, c("site", "BIO1", "BIO5", "BIO6", "BIO7", "BIO12")],
            by = "site", all.x = TRUE)
cat(sprintf("INFO: %d sites after merge with BIO variables\n", nrow(df)))

# -------------------------------------------------------------------
# Helper: find best 3-consecutive-month window (cyclic, months 1-12)
# -------------------------------------------------------------------
# criterion_vals: 12-length numeric vector (month values for assignment)
# target_vals:    12-length numeric vector (values to sum in winning window)
# fn: "max" or "min" to select window maximising or minimising criterion
best_quarter_sum <- function(criterion_vals, target_vals, fn = "max") {
  months12 <- 1:12
  windows  <- lapply(months12, function(m) c(m, (m %% 12) + 1, ((m + 1) %% 12) + 1))
  sums_c   <- sapply(windows, function(w) sum(criterion_vals[w], na.rm = TRUE))
  sums_t   <- sapply(windows, function(w) sum(target_vals[w],    na.rm = TRUE))
  best_w   <- if (fn == "max") which.max(sums_c) else which.min(sums_c)
  sums_t[best_w]
}

# Days per month (standard non-leap year, used for GDD)
days_per_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# -------------------------------------------------------------------
# STEP 3: Compute ENVIREM variables per site
# -------------------------------------------------------------------
cat("INFO: computing 14 ENVIREM variables per site...\n")

env_list <- lapply(seq_len(nrow(df)), function(i) {
  s   <- df[i, ]
  tas <- as.numeric(s[, tas_cols])
  pr  <- as.numeric(s[, pr_cols])
  pet <- as.numeric(s[, pet_cols])

  # Annual PET (mm/year)
  annualPET <- sum(pet, na.rm = TRUE)

  # PET seasonality: CV of monthly PET (%)
  PETseasonality <- if (!all(is.na(pet)) && mean(pet, na.rm = TRUE) != 0)
    (sd(pet, na.rm = TRUE) / mean(pet, na.rm = TRUE)) * 100
  else NA_real_

  # Quarter-based PET — 3 consecutive months, cyclic
  # Warmest quarter: window with max sum of monthly mean temp
  PETWarmestQuarter <- best_quarter_sum(tas, pet, fn = "max")
  # Coldest quarter: window with min sum of monthly mean temp
  PETColdestQuarter <- best_quarter_sum(tas, pet, fn = "min")
  # Wettest quarter: window with max sum of monthly precip
  PETWettestQuarter <- best_quarter_sum(pr, pet, fn = "max")
  # Driest quarter: window with min sum of monthly precip
  PETDriestQuarter  <- best_quarter_sum(pr, pet, fn = "min")

  # Thornthwaite aridity index: annualPET / annualPrecip
  # Values > 1 = water deficit (arid); < 1 = water surplus (humid tropical)
  annualPrecip <- sum(pr, na.rm = TRUE)
  aridityIndexThornthwaite <- if (annualPrecip > 0) annualPET / annualPrecip else NA_real_

  # Climatic moisture index: (P - PET) / (P + PET), range [-1 arid, +1 humid]
  climaticMoistureIndex <- if ((annualPrecip + annualPET) > 0)
    (annualPrecip - annualPET) / (annualPrecip + annualPET)
  else NA_real_

  # Growing degree days base 0°C (°C·days/year)
  growingDegDays0 <- sum(pmax(0, tas) * days_per_month, na.rm = TRUE)

  # Growing degree days base 5°C (°C·days/year)
  growingDegDays5 <- sum(pmax(0, tas - 5) * days_per_month, na.rm = TRUE)

  # Number of months with mean temp > 10°C (0-12)
  monthCountByTemp10 <- sum(tas > 10, na.rm = TRUE)

  # Conrad's continentality index: K = (1.7 × A / sin(|lat| + 10°)) - 14
  # A = amplitude of mean monthly temperatures
  # For tropical sites (all temps > 0), A is typically very small.
  A <- max(tas, na.rm = TRUE) - min(tas, na.rm = TRUE)
  lat_rad <- (abs(as.numeric(s$lat)) + 10) * pi / 180
  continentality <- (1.7 * A / sin(lat_rad)) - 14

  # Emberger's pluviometric quotient:
  # Q = (P × 100) / ((M + 273.15)² - (m + 273.15)²)
  # where P = annual precip (mm), M = mean max temp warmest month (BIO5, °C),
  #       m = mean min temp coldest month (BIO6, °C)
  BIO5  <- as.numeric(s$BIO5)
  BIO6  <- as.numeric(s$BIO6)
  BIO12 <- as.numeric(s$BIO12)
  M_K   <- BIO5 + 273.15
  m_K   <- BIO6 + 273.15
  embergerQ <- if (!is.na(BIO5) && !is.na(BIO6) && !is.na(BIO12) &&
                   (M_K^2 - m_K^2) != 0)
    (BIO12 * 100) / (M_K^2 - m_K^2)
  else NA_real_

  # Rivas-Martínez thermicity index: Ti = (T + M + m) × 10
  # T = annual mean temp (BIO1), M = max temp warmest month (BIO5),
  # m = min temp coldest month (BIO6)
  BIO1 <- as.numeric(s$BIO1)
  thermicityIndex <- if (!is.na(BIO1) && !is.na(BIO5) && !is.na(BIO6))
    (BIO1 + BIO5 + BIO6) * 10
  else NA_real_

  data.frame(
    site                     = s$site,
    lat                      = as.numeric(s$lat),
    long                     = as.numeric(s$long),
    annualPET                = round(annualPET,                3),
    PETseasonality           = round(PETseasonality,           3),
    PETWarmestQuarter        = round(PETWarmestQuarter,        3),
    PETColdestQuarter        = round(PETColdestQuarter,        3),
    PETWettestQuarter        = round(PETWettestQuarter,        3),
    PETDriestQuarter         = round(PETDriestQuarter,         3),
    aridityIndexThornthwaite = round(aridityIndexThornthwaite, 4),
    climaticMoistureIndex    = round(climaticMoistureIndex,    4),
    growingDegDays0          = round(growingDegDays0,          1),
    growingDegDays5          = round(growingDegDays5,          1),
    monthCountByTemp10       = monthCountByTemp10,
    continentality           = round(continentality,           3),
    embergerQ                = round(embergerQ,                4),
    thermicityIndex          = round(thermicityIndex,          2),
    stringsAsFactors = FALSE
  )
})

env_df <- do.call(rbind, env_list)
cat(sprintf("INFO: %d sites × 14 ENVIREM variables computed\n", nrow(env_df)))

# -------------------------------------------------------------------
# STEP 4: Write physical-unit output
# -------------------------------------------------------------------
env_cols <- c("annualPET", "PETseasonality", "PETWarmestQuarter", "PETColdestQuarter",
              "PETWettestQuarter", "PETDriestQuarter", "aridityIndexThornthwaite",
              "climaticMoistureIndex", "growingDegDays0", "growingDegDays5",
              "monthCountByTemp10", "continentality", "embergerQ", "thermicityIndex")

out_raw <- file.path(out_dir, "data/envirem/envirem_env_per_site.csv")
dir.create(dirname(out_raw), recursive = TRUE, showWarnings = FALSE)
write.csv(env_df, file = out_raw, row.names = FALSE, quote = FALSE)
cat("Physical units values written:", out_raw, "\n")

# -------------------------------------------------------------------
# STEP 5: Z-score standardisation for LFMM/RDA
# -------------------------------------------------------------------
env_scaled <- env_df
env_scaled[env_cols] <- as.data.frame(scale(env_df[env_cols]))

out_scaled <- file.path(out_dir, "data/envirem/envirem_env_scaled.csv")
write.csv(env_scaled, file = out_scaled, row.names = FALSE, quote = FALSE)
cat("Z-score scaled values written:", out_scaled, "\n")

# -------------------------------------------------------------------
# STEP 6: Summary statistics
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
out_summ <- file.path(out_dir, "data/envirem/envirem_env_summary.tsv")
write.table(summary_df, file = out_summ, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary statistics written:", out_summ, "\n")
print(summary_df)

# -------------------------------------------------------------------
# STEP 7: Correlation matrix heatmap
# -------------------------------------------------------------------
cor_mat  <- cor(env_df[env_cols], use = "pairwise.complete.obs")
cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Var1", "Var2", "Correlation")
cor_long$Var1 <- factor(cor_long$Var1, levels = env_cols)
cor_long$Var2 <- factor(cor_long$Var2, levels = env_cols)

p_cor <- ggplot(cor_long, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 1.7, color = "black") +
  scale_fill_gradientn(
    colours = c("#d73027", "#ffffff", "#4575b4"),
    limits  = c(-1, 1), name = "Pearson r"
  ) +
  labs(title = "ENVIREM variables — Pearson correlation matrix",
       x = NULL, y = NULL) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 7.5),
    axis.text.y       = element_text(size = 7.5),
    plot.title        = element_text(face = "bold"),
    legend.key.height = unit(0.8, "cm")
  )

ggsave(file.path(plot_dir, "envirem_correlation.png"), p_cor,
       width = 11, height = 10, dpi = 300)
ggsave(file.path(plot_dir, "envirem_correlation.pdf"), p_cor,
       width = 11, height = 10)
cat("Correlation heatmap written\n")

# -------------------------------------------------------------------
# STEP 8: Identify highly correlated pairs (|r| >= 0.8)
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
  out_hc <- file.path(out_dir, "data/envirem/envirem_high_correlations.tsv")
  write.table(hc_df, file = out_hc, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("INFO: %d highly correlated ENVIREM pairs (|r| >= 0.8): %s\n",
              nrow(hc_df), out_hc))
  print(hc_df)
} else {
  cat("INFO: no highly correlated ENVIREM pairs (|r| >= 0.8) detected\n")
}

# -------------------------------------------------------------------
# STEP 9: Barplots per ENVIREM variable across sites
# -------------------------------------------------------------------
var_labels <- c(
  annualPET                = "Annual PET (mm/year)",
  PETseasonality           = "PET seasonality (CV, %)",
  PETWarmestQuarter        = "PET warmest quarter (mm)",
  PETColdestQuarter        = "PET coldest quarter (mm)",
  PETWettestQuarter        = "PET wettest quarter (mm)",
  PETDriestQuarter         = "PET driest quarter (mm)",
  aridityIndexThornthwaite = "Thornthwaite aridity index (PET/P)",
  climaticMoistureIndex    = "Climatic moisture index",
  growingDegDays0          = "Growing degree days (base 0\u00b0C)",
  growingDegDays5          = "Growing degree days (base 5\u00b0C)",
  monthCountByTemp10       = "Months with mean temp > 10\u00b0C",
  continentality           = "Conrad's continentality index",
  embergerQ                = "Emberger's pluviometric quotient (Q)",
  thermicityIndex          = "Rivas-Mart\u00ednez thermicity index"
)

for (v in env_cols) {
  df_v <- env_df[!is.na(env_df[[v]]), c("site", v)]
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

cat("\nDONE 08.0.2 ENVIREM computation completed\n")
cat("Outputs:\n")
cat("  Physical units :", out_raw,    "\n")
cat("  Scaled         :", out_scaled, "\n")
cat("  Summary        :", out_summ,   "\n")

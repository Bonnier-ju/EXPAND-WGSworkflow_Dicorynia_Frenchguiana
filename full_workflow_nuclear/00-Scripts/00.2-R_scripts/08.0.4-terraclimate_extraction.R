#!/usr/bin/env Rscript
# 08.0.4-terraclimate_extraction.R
#
# Post-process raw monthly TerraClimate values extracted by
# 08.0.4-terraclimate_extraction.slurm.
#
# Input:
#   terraclimate_raw_monthly.tsv
#     Columns: site lat long variable year m01..m12
#     All values in physical units (GDAL applies scale_factor from NetCDF metadata).
#
# Annual aggregation rules:
#   SUM  : aet, def, ppt               (monthly fluxes → annual totals)
#   MEAN : PDSI, soil, srad, vpd       (stocks/indices → annual means)
#
# Derived output columns (11 — unique vs CHELSA/ENVIREM):
#   aet_mean, aet_sd     — actual evapotranspiration
#   def_mean, def_sd     — climatic water deficit
#   ppt_sd               — inter-annual precipitation variability (mean ≈ BIO12)
#   PDSI_mean, PDSI_sd   — Palmer Drought Severity Index
#   soil_mean, soil_sd   — soil moisture
#   srad_mean            — solar radiation (SD not prioritised)
#   vpd_mean             — vapour pressure deficit (SD not prioritised)
#
# Outputs:
#   data/terraclimate/terraclimate_env_per_site.csv
#     11 variables per site over 1991–2020
#   data/terraclimate/terraclimate_env_scaled.csv
#     z-score standardised (mean=0, sd=1) version of the above
#   data/terraclimate/terraclimate_env_summary.tsv
#     Descriptive statistics per variable
#   data/terraclimate/terraclimate_high_correlations.tsv
#     Variable pairs with |r| >= 0.70
#   plots/terraclimate_correlation.png/pdf
#   plots/terraclimate_{var}_per_site.png
#
# Args: raw_tsv  out_dir  year_start  year_end

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else {
  warning("sites_couleurs.csv not found — site labels will be black")
  c()
}

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
raw_tsv    <- if (length(args) >= 1) args[1] else stop("raw_tsv required")
out_dir    <- if (length(args) >= 2) args[2] else stop("out_dir required")
year_start <- as.integer(if (length(args) >= 3) args[3] else 1991)
year_end   <- as.integer(if (length(args) >= 4) args[4] else 2020)

data_dir <- file.path(out_dir, "data/terraclimate")
plot_dir <- file.path(out_dir, "plots")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Reading raw monthly data: %s\n", raw_tsv))
raw <- read.table(raw_tsv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, na.strings = c("NA", ""))

month_cols <- paste0("m", sprintf("%02d", 1:12))
cat(sprintf("  %d rows x %d columns\n", nrow(raw), ncol(raw)))

# Sanity check: all month columns present
missing_cols <- setdiff(month_cols, colnames(raw))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing month columns in raw TSV: %s", paste(missing_cols, collapse = ", ")))
}

# Convert month columns to numeric (NA-safe)
raw[, month_cols] <- lapply(raw[, month_cols], as.numeric)

# -------------------------------------------------------------------
# STEP 1: Annual aggregation
# Summation for fluxes; mean for stocks/indices.
# -------------------------------------------------------------------
cat("STEP 1: Annual aggregation ...\n")

SUM_VARS  <- c("aet", "def", "ppt")
MEAN_VARS <- c("PDSI", "soil", "srad", "vpd")

raw$annual_value <- vapply(seq_len(nrow(raw)), function(i) {
  vals <- as.numeric(raw[i, month_cols])
  if (all(is.na(vals))) return(NA_real_)
  var <- raw$variable[i]
  if (var %in% SUM_VARS)  return(sum(vals,  na.rm = TRUE))
  if (var %in% MEAN_VARS) return(mean(vals, na.rm = TRUE))
  mean(vals, na.rm = TRUE)
}, numeric(1))

# Keep only the 30-year analysis period
annual <- raw[raw$year >= year_start & raw$year <= year_end,
              c("site", "lat", "long", "variable", "year", "annual_value")]
annual$annual_value <- as.numeric(annual$annual_value)
cat(sprintf("  %d site × variable × year records after filtering to %d–%d\n",
            nrow(annual), year_start, year_end))

# -------------------------------------------------------------------
# STEP 2: 30-year climatological statistics (mean and SD)
# One row per site; after dropping redundant statistics → 11 columns.
# -------------------------------------------------------------------
cat("STEP 2: Computing 30-year statistics ...\n")

all_vars <- c(SUM_VARS, MEAN_VARS)

# Compute mean and SD per site × variable
clim_stats <- do.call(rbind, lapply(split(annual, list(annual$site, annual$variable)), function(d) {
  if (nrow(d) == 0) return(NULL)
  data.frame(
    site     = d$site[1],
    lat      = d$lat[1],
    long     = d$long[1],
    variable = d$variable[1],
    n_years  = sum(!is.na(d$annual_value)),
    mean_val = mean(d$annual_value, na.rm = TRUE),
    sd_val   = if (sum(!is.na(d$annual_value)) > 1)
                 sd(d$annual_value, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
clim_stats <- clim_stats[!is.null(clim_stats), ]
rownames(clim_stats) <- NULL

# Pivot wide: one row per site, columns = {var}_mean and {var}_sd
sites_meta <- unique(clim_stats[, c("site", "lat", "long")])
site_order <- sort(sites_meta$site)
sites_meta <- sites_meta[order(sites_meta$site), ]

wide <- sites_meta
for (v in all_vars) {
  sub  <- clim_stats[clim_stats$variable == v, c("site", "mean_val", "sd_val")]
  mean_col <- paste0(v, "_mean")
  sd_col   <- paste0(v, "_sd")
  names(sub)[2:3] <- c(mean_col, sd_col)
  wide <- merge(wide, sub, by = "site", all.x = TRUE, sort = FALSE)
}
wide <- wide[order(wide$site), ]

# Drop redundant statistics: ppt_mean ≈ CHELSA BIO12; srad_sd and vpd_sd
# not ecologically prioritised for GEA.
DROP_COLS <- c("ppt_mean", "srad_sd", "vpd_sd")
wide <- wide[, !colnames(wide) %in% DROP_COLS]

env_cols <- grep("_mean$|_sd$", colnames(wide), value = TRUE)
# env_cols (11): aet_mean, aet_sd, def_mean, def_sd, ppt_sd,
#                PDSI_mean, PDSI_sd, soil_mean, soil_sd, srad_mean, vpd_mean
cat(sprintf("  %d sites, %d derived variables\n", nrow(wide), length(env_cols)))

# -------------------------------------------------------------------
# STEP 3: Write per-site CSV and scaled version
# -------------------------------------------------------------------
cat("STEP 3: Writing output tables ...\n")

out_csv    <- file.path(data_dir, "terraclimate_env_per_site.csv")
out_scaled <- file.path(data_dir, "terraclimate_env_scaled.csv")
out_summ   <- file.path(data_dir, "terraclimate_env_summary.tsv")
out_hicorr <- file.path(data_dir, "terraclimate_high_correlations.tsv")

write.csv(wide, out_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_csv))

# z-score standardisation (scale each env variable)
scaled      <- wide
scaled[, env_cols] <- scale(wide[, env_cols])
write.csv(scaled, out_scaled, row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_scaled))

# Descriptive stats summary
summ_rows <- lapply(env_cols, function(col) {
  vals <- as.numeric(wide[[col]])
  data.frame(
    variable = col,
    n_valid  = sum(!is.na(vals)),
    mean     = round(mean(vals, na.rm = TRUE), 4),
    sd       = round(sd(vals, na.rm = TRUE), 4),
    min      = round(min(vals, na.rm = TRUE), 4),
    max      = round(max(vals, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
})
summ_df <- do.call(rbind, summ_rows)
write.table(summ_df, out_summ, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_summ))

# High correlations (|r| >= 0.70)
env_mat <- as.matrix(wide[, env_cols])
valid_cols <- env_cols[apply(!is.na(env_mat), 2, any)]
if (length(valid_cols) > 1) {
  cor_mat <- cor(wide[, valid_cols], use = "pairwise.complete.obs")
  hi_corr_rows <- lapply(seq_len(nrow(cor_mat)), function(i) {
    lapply(seq_len(ncol(cor_mat)), function(j) {
      if (j <= i) return(NULL)
      r <- cor_mat[i, j]
      if (!is.na(r) && abs(r) >= 0.70)
        data.frame(var1 = rownames(cor_mat)[i],
                   var2 = colnames(cor_mat)[j],
                   r    = round(r, 3), stringsAsFactors = FALSE)
      else NULL
    })
  })
  hi_corr_df <- do.call(rbind, unlist(hi_corr_rows, recursive = FALSE))
  if (!is.null(hi_corr_df) && nrow(hi_corr_df) > 0) {
    hi_corr_df <- hi_corr_df[order(abs(hi_corr_df$r), decreasing = TRUE), ]
    write.table(hi_corr_df, out_hicorr, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  %d high-correlation pairs (|r|>=0.70) written to %s\n",
                nrow(hi_corr_df), out_hicorr))
  }
}

# -------------------------------------------------------------------
# STEP 4: Correlation heatmap
# -------------------------------------------------------------------
cat("STEP 4: Correlation heatmap ...\n")

if (length(valid_cols) > 1) {
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("var1", "var2", "r")
  cor_df$var1 <- factor(cor_df$var1, levels = valid_cols)
  cor_df$var2 <- factor(cor_df$var2, levels = valid_cols)

  p_cor <- ggplot(cor_df, aes(x = var2, y = var1, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", r)),
              size = 2.2, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Pearson r") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(
      title = sprintf("TerraClimate variable correlations (n=%d sites, %d–%d)",
                      nrow(wide), year_start, year_end),
      x = NULL, y = NULL
    ) +
    theme(
      axis.text        = element_text(size = 7),
      legend.key.height = unit(0.8, "cm"),
      panel.grid       = element_blank()
    )

  ggsave(file.path(plot_dir, "terraclimate_correlation.png"), p_cor,
         width = 14, height = 12, dpi = 300)
  ggsave(file.path(plot_dir, "terraclimate_correlation.pdf"), p_cor,
         width = 14, height = 12)
  cat("  Correlation heatmap written\n")
}

# -------------------------------------------------------------------
# STEP 5: Per-site bar plots for retained variables
# Includes both _mean and _sd columns; subtitle adapts accordingly.
# -------------------------------------------------------------------
cat("STEP 5: Per-site variable plots ...\n")

plot_vars <- intersect(
  c("def_mean", "aet_mean", "soil_mean", "PDSI_mean", "vpd_mean", "srad_mean", "ppt_sd"),
  colnames(wide)
)

units_map <- c(
  def_mean  = "mm/yr",
  aet_mean  = "mm/yr",
  soil_mean = "mm",
  PDSI_mean = "index",
  vpd_mean  = "kPa",
  srad_mean = "W/m²/day",
  ppt_sd    = "mm/yr"
)

for (v in plot_vars) {
  plot_df <- wide[, c("site", v)]
  colnames(plot_df) <- c("site", "value")
  plot_df <- plot_df[!is.na(plot_df$value), ]
  plot_df <- plot_df[order(plot_df$value, decreasing = TRUE), ]
  plot_df$site <- factor(plot_df$site, levels = plot_df$site)

  unit_lbl   <- if (v %in% names(units_map)) units_map[[v]] else ""
  var_base   <- sub("_(mean|sd)$", "", v)
  var_suffix <- ifelse(grepl("_sd$", v), "inter-annual SD", "30-yr mean")

  p <- ggplot(plot_df, aes(x = site, y = value, fill = site)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = site_colors, na.value = "grey60", guide = "none") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(
      title    = sprintf("TerraClimate — %s per site", var_base),
      subtitle = sprintf("%s (%d–%d)", var_suffix, year_start, year_end),
      x = NULL,
      y = sprintf("%s (%s)", var_base, unit_lbl)
    ) +
    theme(
      axis.text.x   = element_text(size = 7, hjust = 1),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )

  png_path <- file.path(plot_dir, sprintf("terraclimate_%s_per_site.png", v))
  ggsave(png_path, p, width = 10, height = 5, dpi = 300)
}
cat(sprintf("  %d per-site plots written\n", length(plot_vars)))

# -------------------------------------------------------------------
# STEP 6: Interannual variability plot (SD vs mean) for variables
# that have both statistics in the output (aet, def, PDSI, soil).
# -------------------------------------------------------------------
cat("STEP 6: Variability overview ...\n")

var_pairs <- list(
  list(mean = "def_mean",  sd = "def_sd",  label = "CWD (def)"),
  list(mean = "aet_mean",  sd = "aet_sd",  label = "Actual ET (aet)"),
  list(mean = "PDSI_mean", sd = "PDSI_sd", label = "PDSI"),
  list(mean = "soil_mean", sd = "soil_sd", label = "Soil moisture")
)

for (vp in var_pairs) {
  if (!all(c(vp$mean, vp$sd) %in% colnames(wide))) next
  df_vp <- wide[, c("site", vp$mean, vp$sd)]
  colnames(df_vp) <- c("site", "mean_val", "sd_val")
  df_vp <- df_vp[!is.na(df_vp$mean_val), ]

  p_var <- ggplot(df_vp, aes(x = mean_val, y = sd_val,
                              color = site, label = site)) +
    geom_point(size = 3) +
    geom_text(nudge_y = 0.01 * diff(range(df_vp$sd_val, na.rm = TRUE)),
              size = 2.5, show.legend = FALSE) +
    scale_color_manual(values = site_colors, na.value = "grey60", guide = "none") +
    labs(
      title    = sprintf("TerraClimate %s — interannual variability", vp$label),
      subtitle = sprintf("%d–%d | each point = 1 site", year_start, year_end),
      x = sprintf("30-yr mean %s", vp$mean),
      y = sprintf("30-yr SD %s",   vp$sd)
    ) +
    theme(plot.subtitle = element_text(size = 8, color = "grey40"))

  vname <- sub("_mean$", "", vp$mean)
  ggsave(file.path(plot_dir, sprintf("terraclimate_%s_variability.png", vname)),
         p_var, width = 8, height = 6, dpi = 300)
}
cat("  Variability plots written\n")

cat(sprintf("\nDONE 08.0.4-terraclimate_extraction.R\n"))
cat(sprintf("  %d sites | %d variables | %d–%d\n",
            nrow(wide), length(env_cols), year_start, year_end))
cat(sprintf("  Per-site CSV   : %s\n", out_csv))
cat(sprintf("  Scaled CSV     : %s\n", out_scaled))
cat(sprintf("  Plots          : %s/terraclimate_*\n", plot_dir))

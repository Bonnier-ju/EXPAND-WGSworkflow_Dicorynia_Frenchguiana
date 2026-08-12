#!/usr/bin/env Rscript
# 08.1b.2-derived_variables.R
#
# Compute the 10 custom TerraClimate-derived environmental variables from the
# raw monthly data produced by 08.1b.1-terraclimate_monthly_extraction.slurm.
#
# IMPORTANT - scaling: values in the raw TSV are ALREADY in physical units.
# The raw monthly extraction (08.1b.1) uses gdallocationinfo on GDAL-opened
# NetCDF subdatasets, which applies the NetCDF scale_factor/add_offset
# automatically. Do NOT re-apply manual scaling factors (e.g. the x0.1 / x0.01
# used in Sylvain Schmitt's GEE-based pipeline, which pulls raw unscaled band
# values) - that would silently corrupt every value by 10x or 100x.
#
# Dry season definition (confirmed 2026-08-11, following Sylvain Schmitt's
# water-balance approach rather than the original plan's ppt<100mm threshold):
#   a month is DRY if aet > ppt (actual evapotranspiration exceeds rainfall).
#
# The 10 derived variables (per year, then mean across 1970-2000):
#   1. precip_annual         sum(ppt) over 12 months
#   2. dry_season_duration   count of dry months (aet > ppt)
#   3. dry_season_intensity  sum(aet - ppt) over dry months [deficit magnitude
#                            - confirmed 2026-08-11, replaces the original
#                            plan's "sum of ppt during dry months"]
#   4. wet_season_duration   count of wet months (aet <= ppt)
#   5. wet_season_intensity  sum(ppt - aet) over wet months [surplus magnitude
#                            - symmetric extension of the dry_season_intensity
#                            change, so that precip_contrast (below) subtracts
#                            two quantities in the same unit (water-balance mm)
#                            rather than mixing a deficit with a raw rainfall
#                            total. Not explicitly confirmed - flagged here,
#                            revisit if a pure-rainfall wet intensity is wanted.]
#   6. precip_contrast       wet_season_intensity - dry_season_intensity (same year)
#   7. cwd_annual            sum(def) over 12 months
#   8. tmax_dry_season       mean(tmax) over dry months; if no dry month that
#                            year, fall back to max(tmax) of that year
#   9. soil_moisture_min     min(soil) over 12 months -> mean of ANNUAL MINIMA
#  10. vpd_max_annual        max(vpd) over 12 months -> mean of ANNUAL MAXIMA
#                            [renamed from vpd_mean_annual, confirmed
#                            2026-08-11, per Sylvain's general advice to
#                            prefer mean-of-annual-maxima over mean-of-means]
#
# Input:
#   terraclimate_raw_monthly_19sites.tsv
#     Columns: site lat long variable year m01..m12 (physical units)
#
# Outputs:
#   data/terraclimate_derived_variables_19sites.csv   (site + 10 variables)
#   data/terraclimate_derived_variables_summary.tsv    (descriptive stats + CV)
#   plots/derived_{variable}_per_site.png
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
  warning("sites_couleurs.csv not found - site labels will be grey")
  c()
}

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args       <- commandArgs(trailingOnly = TRUE)
raw_tsv    <- if (length(args) >= 1) args[1] else stop("raw_tsv required")
out_dir    <- if (length(args) >= 2) args[2] else stop("out_dir required")
year_start <- as.integer(if (length(args) >= 3) args[3] else 1970)
year_end   <- as.integer(if (length(args) >= 4) args[4] else 2000)

data_dir <- file.path(out_dir, "data")
plot_dir <- file.path(out_dir, "plots")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Guard clause: verify input before running ----
if (!file.exists(raw_tsv) || file.info(raw_tsv)$size == 0) {
  stop(sprintf("ERROR: required input missing or empty: %s", raw_tsv))
}

cat(sprintf("Reading raw monthly data: %s\n", raw_tsv))
raw <- read.table(raw_tsv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, na.strings = c("NA", ""))
month_cols <- paste0("m", sprintf("%02d", 1:12))
missing_cols <- setdiff(month_cols, colnames(raw))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing month columns in raw TSV: %s", paste(missing_cols, collapse = ", ")))
}
raw[, month_cols] <- lapply(raw[, month_cols], as.numeric)
raw <- raw[raw$year >= year_start & raw$year <= year_end, ]
cat(sprintf("  %d rows after filtering to %d-%d\n", nrow(raw), year_start, year_end))

REQUIRED_VARS <- c("aet", "def", "pet", "ppt", "soil", "tmax", "tmin", "vpd")
missing_vars <- setdiff(REQUIRED_VARS, unique(raw$variable))
if (length(missing_vars) > 0) {
  stop(sprintf("Missing required variables in raw TSV: %s", paste(missing_vars, collapse = ", ")))
}

# -------------------------------------------------------------------
# STEP 1: Reshape to one row per site x year x month, one column per variable
# -------------------------------------------------------------------
cat("STEP 1: reshaping to site x year x month ...\n")

long <- do.call(rbind, lapply(month_cols, function(mc) {
  data.frame(
    site  = raw$site,
    year  = raw$year,
    variable = raw$variable,
    month = as.integer(sub("^m", "", mc)),
    value = raw[[mc]],
    stringsAsFactors = FALSE
  )
}))

wide <- reshape(long, idvar = c("site", "year", "month"),
                 timevar = "variable", direction = "wide")
colnames(wide) <- sub("^value\\.", "", colnames(wide))
wide <- wide[order(wide$site, wide$year, wide$month), ]
cat(sprintf("  %d site x year x month rows\n", nrow(wide)))

# -------------------------------------------------------------------
# STEP 2: Per-year metrics, per site
# -------------------------------------------------------------------
cat("STEP 2: computing per-year metrics ...\n")

compute_year <- function(d) {
  # d = 12 rows (months) for one site x year
  is_dry <- d$aet > d$ppt
  dry_n  <- sum(is_dry, na.rm = TRUE)
  wet_n  <- sum(!is_dry, na.rm = TRUE)

  dry_intensity <- if (dry_n > 0) sum((d$aet - d$ppt)[is_dry], na.rm = TRUE) else 0
  wet_intensity <- if (wet_n > 0) sum((d$ppt - d$aet)[!is_dry], na.rm = TRUE) else 0

  tmax_dry <- if (dry_n > 0) mean(d$tmax[is_dry], na.rm = TRUE) else max(d$tmax, na.rm = TRUE)

  data.frame(
    site                 = d$site[1],
    year                 = d$year[1],
    precip_annual        = sum(d$ppt, na.rm = TRUE),
    dry_season_duration  = dry_n,
    dry_season_intensity = dry_intensity,
    wet_season_duration  = wet_n,
    wet_season_intensity = wet_intensity,
    precip_contrast      = wet_intensity - dry_intensity,
    cwd_annual           = sum(d$def, na.rm = TRUE),
    tmax_dry_season      = tmax_dry,
    soil_moisture_min    = min(d$soil, na.rm = TRUE),
    vpd_max_annual       = max(d$vpd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

yearly <- do.call(rbind, lapply(split(wide, list(wide$site, wide$year), drop = TRUE), compute_year))
rownames(yearly) <- NULL
cat(sprintf("  %d site x year records\n", nrow(yearly)))

# -------------------------------------------------------------------
# STEP 3: Mean across years, per site (30-year normal, 1970-2000)
# -------------------------------------------------------------------
cat("STEP 3: averaging across years ...\n")

derived_vars <- c("precip_annual", "dry_season_duration", "dry_season_intensity",
                   "wet_season_duration", "wet_season_intensity", "precip_contrast",
                   "cwd_annual", "tmax_dry_season", "soil_moisture_min", "vpd_max_annual")

per_site <- do.call(rbind, lapply(split(yearly, yearly$site), function(d) {
  row <- data.frame(site = d$site[1], n_years = nrow(d), stringsAsFactors = FALSE)
  for (v in derived_vars) row[[v]] <- mean(d[[v]], na.rm = TRUE)
  row
}))
per_site <- per_site[order(per_site$site), ]
rownames(per_site) <- NULL

n_expected_years <- year_end - year_start + 1
if (any(per_site$n_years != n_expected_years)) {
  warning(sprintf("Some sites have != %d years of data - check for failed downloads:\n%s",
                  n_expected_years,
                  paste(sprintf("  %s: %d years", per_site$site[per_site$n_years != n_expected_years],
                                per_site$n_years[per_site$n_years != n_expected_years]),
                        collapse = "\n")))
}

# -------------------------------------------------------------------
# STEP 4: Write outputs
# -------------------------------------------------------------------
cat("STEP 4: writing outputs ...\n")

out_csv <- file.path(data_dir, "terraclimate_derived_variables_19sites.csv")
write.csv(per_site[, c("site", derived_vars)], out_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_csv))

summ_rows <- lapply(derived_vars, function(v) {
  vals <- per_site[[v]]
  m <- mean(vals, na.rm = TRUE)
  s <- sd(vals, na.rm = TRUE)
  data.frame(
    variable = v,
    n_sites  = sum(!is.na(vals)),
    mean     = round(m, 4),
    sd       = round(s, 4),
    min      = round(min(vals, na.rm = TRUE), 4),
    max      = round(max(vals, na.rm = TRUE), 4),
    cv_pct   = round(100 * s / m, 2),
    stringsAsFactors = FALSE
  )
})
summ_df <- do.call(rbind, summ_rows)
out_summ <- file.path(data_dir, "terraclimate_derived_variables_summary.tsv")
write.table(summ_df, out_summ, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_summ))

low_cv <- summ_df$variable[summ_df$cv_pct < 5]
if (length(low_cv) > 0) {
  cat(sprintf("  NOTE: low spatial variation (CV<5%%) - candidates for exclusion in 08.1b.5: %s\n",
              paste(low_cv, collapse = ", ")))
}

# -------------------------------------------------------------------
# STEP 5: Per-site bar plots
# -------------------------------------------------------------------
cat("STEP 5: per-site plots ...\n")

for (v in derived_vars) {
  plot_df <- per_site[, c("site", v)]
  colnames(plot_df) <- c("site", "value")
  plot_df <- plot_df[order(plot_df$value, decreasing = TRUE), ]
  plot_df$site <- factor(plot_df$site, levels = plot_df$site)

  p <- ggplot(plot_df, aes(x = site, y = value, fill = site)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = site_colors, na.value = "grey60", guide = "none") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(
      title    = sprintf("TerraClimate derived - %s per site", v),
      subtitle = sprintf("Mean across %d-%d", year_start, year_end),
      x = NULL, y = v
    ) +
    theme(
      axis.text.x   = element_text(size = 7, hjust = 1),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )

  ggsave(file.path(plot_dir, sprintf("derived_%s_per_site.png", v)), p,
         width = 10, height = 5, dpi = 300)
}
cat(sprintf("  %d per-site plots written\n", length(derived_vars)))

cat(sprintf("\nDONE 08.1b.2-derived_variables.R\n"))
cat(sprintf("  %d sites | %d variables | %d-%d\n",
            nrow(per_site), length(derived_vars), year_start, year_end))

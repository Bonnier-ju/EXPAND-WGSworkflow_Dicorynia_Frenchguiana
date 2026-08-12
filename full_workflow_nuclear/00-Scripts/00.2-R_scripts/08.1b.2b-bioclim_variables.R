#!/usr/bin/env Rscript
# 08.1b.2b-bioclim_variables.R
#
# Compute the 19 standard bioclimatic variables (BIO1-BIO19) from the
# TerraClimate monthly climatology (1970-2000), using dismo::biovars() -
# advised by Sylvain Schmitt (2026-08-11) so that BIO-equivalent variables
# (e.g. BIO6, BIO17, which mattered in the previous CHELSA-based run) can be
# recomputed from TerraClimate alone, instead of mixing CHELSA + TerraClimate
# as two different sources for the same kind of variable.
#
# biovars() expects one 30-year MONTHLY CLIMATOLOGY per site (12 values per
# variable: the mean of each calendar month across 1970-2000) - NOT the
# per-year-then-averaged approach used in 08.1b.2. This is the standard
# WorldClim/CHELSA convention for bioclimatic variables.
#
# IMPORTANT - scaling: as in 08.1b.2, the raw monthly TSV is already in
# physical units (GDAL applied scale_factor on extraction) - no manual
# rescaling needed here.
#
# Input:
#   terraclimate_raw_monthly_19sites.tsv (from 08.1b.1)
#     Columns: site lat long variable year m01..m12
#     Uses ppt, tmin, tmax only.
#
# Output:
#   data/bioclim_variables_19sites.csv       (site + BIO1..BIO19)
#   data/bioclim_variables_summary.tsv       (descriptive stats + CV)
#
# Args: raw_tsv  out_dir  year_start  year_end

suppressPackageStartupMessages({
  ok_dismo <- requireNamespace("dismo", quietly = TRUE)
})
if (!ok_dismo) {
  stop(paste(
    "Package 'dismo' not available in R_LIBS_USER.",
    "Check with: search_R_package dismo",
    "If absent, install into the project renv (see CLAUDE.md §5) before rerunning.",
    sep = "\n"
  ))
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
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Guard clause: verify input before running ----
if (!file.exists(raw_tsv) || file.info(raw_tsv)$size == 0) {
  stop(sprintf("ERROR: required input missing or empty: %s", raw_tsv))
}

cat(sprintf("Reading raw monthly data: %s\n", raw_tsv))
raw <- read.table(raw_tsv, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, na.strings = c("NA", ""))
month_cols <- paste0("m", sprintf("%02d", 1:12))
raw[, month_cols] <- lapply(raw[, month_cols], as.numeric)
raw <- raw[raw$year >= year_start & raw$year <= year_end, ]

NEEDED_VARS <- c("ppt", "tmin", "tmax")
missing_vars <- setdiff(NEEDED_VARS, unique(raw$variable))
if (length(missing_vars) > 0) {
  stop(sprintf("Missing required variables in raw TSV: %s", paste(missing_vars, collapse = ", ")))
}
raw <- raw[raw$variable %in% NEEDED_VARS, ]

# -------------------------------------------------------------------
# STEP 1: 30-year monthly climatology per site (mean of each calendar
# month across years) - one 12-value row per site x variable.
# -------------------------------------------------------------------
cat(sprintf("STEP 1: computing %d-%d monthly climatology ...\n", year_start, year_end))

sites <- sort(unique(raw$site))

clim_matrix <- function(varname) {
  sub <- raw[raw$variable == varname, ]
  m <- t(vapply(sites, function(s) {
    d <- sub[sub$site == s, month_cols]
    vapply(month_cols, function(mc) mean(d[[mc]], na.rm = TRUE), numeric(1))
  }, numeric(12)))
  rownames(m) <- sites
  colnames(m) <- month_cols
  m
}

prec_mat <- clim_matrix("ppt")
tmin_mat <- clim_matrix("tmin")
tmax_mat <- clim_matrix("tmax")

cat(sprintf("  %d sites x 12 months (ppt, tmin, tmax)\n", length(sites)))

# -------------------------------------------------------------------
# STEP 2: BIO1-BIO19 via dismo::biovars()
# -------------------------------------------------------------------
cat("STEP 2: computing BIO1-BIO19 (dismo::biovars) ...\n")

bio <- dismo::biovars(prec = prec_mat, tmin = tmin_mat, tmax = tmax_mat)
bio_df <- data.frame(site = sites, bio, stringsAsFactors = FALSE)
colnames(bio_df) <- sub("^bio", "BIO", colnames(bio_df))

bio_vars <- grep("^BIO", colnames(bio_df), value = TRUE)
cat(sprintf("  %d BIO variables computed for %d sites\n", length(bio_vars), nrow(bio_df)))

# -------------------------------------------------------------------
# STEP 3: Write outputs
# -------------------------------------------------------------------
cat("STEP 3: writing outputs ...\n")

out_csv <- file.path(data_dir, "bioclim_variables_19sites.csv")
write.csv(bio_df, out_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_csv))

summ_rows <- lapply(bio_vars, function(v) {
  vals <- bio_df[[v]]
  m <- mean(vals, na.rm = TRUE)
  s <- sd(vals, na.rm = TRUE)
  data.frame(
    variable = v,
    n_sites  = sum(!is.na(vals)),
    mean     = round(m, 4),
    sd       = round(s, 4),
    min      = round(min(vals, na.rm = TRUE), 4),
    max      = round(max(vals, na.rm = TRUE), 4),
    cv_pct   = round(100 * abs(s / m), 2),
    stringsAsFactors = FALSE
  )
})
summ_df <- do.call(rbind, summ_rows)
out_summ <- file.path(data_dir, "bioclim_variables_summary.tsv")
write.table(summ_df, out_summ, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_summ))

low_cv <- summ_df$variable[summ_df$cv_pct < 5]
if (length(low_cv) > 0) {
  cat(sprintf("  NOTE: low spatial variation (CV<5%%) - candidates for exclusion in 08.1b.5: %s\n",
              paste(low_cv, collapse = ", ")))
}

cat(sprintf("\nDONE 08.1b.2b-bioclim_variables.R\n"))
cat(sprintf("  %d sites | %d BIO variables | climatology %d-%d\n",
            nrow(bio_df), length(bio_vars), year_start, year_end))

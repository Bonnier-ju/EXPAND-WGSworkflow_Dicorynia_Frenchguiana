#!/usr/bin/env Rscript
# 08.1b.5-manual_selection_prep.R
#
# Semi-manual step: prepares the information needed for Julien to make an
# informed ecological selection of ~8-10 final variables from the 31
# candidates. Does NOT pick variables itself and does NOT invent an
# ecological "gradient" narrative - the PCA that would identify real
# gradients only runs in 08.1b.6, AFTER this manual selection, so any
# gradient-based recommendation here would be circular. Restructured versus
# the original plan on this point (flagged, not silent).
#
# What this script does provide, mechanically, from the correlation results
# of 08.1b.4:
#   - correlation groups: connected components of the |r|>=0.70 graph
#     (a "group" = a set of variables all inter-correlated above threshold;
#     singletons are variables with no correlated partner)
#   - CV% per variable (low CV = little spatial variation across the 19
#     sites = weak candidate for driving a GEA signal)
#   - a one-line description per variable (from the variable's own
#     definition - standard BIO1-19 definitions, or the formulas used in
#     08.1b.2/08.1b.3 for the custom/topo variables)
#
# Output:
#   correlation_analysis/variable_selection_template.txt
#     Editable template (# = excluded). Uncomment (remove #) the variables
#     you want to KEEP. Save in place. Feeds 08.1b.6 (PCA) and 08.1b.7 (RDA).
#   correlation_analysis/variable_selection_summary.tsv
#     Same information as a structured table.
#
# Args: out_dir

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else stop("out_dir required")

corr_dir <- file.path(out_dir, "correlation_analysis")

summary_tsv <- file.path(corr_dir, "variable_summary.tsv")
pairs_tsv   <- file.path(corr_dir, "correlated_pairs_table.tsv")

# ---- Guard clause: verify inputs before running ----
for (f in c(summary_tsv, pairs_tsv)) {
  if (!file.exists(f)) stop(sprintf("ERROR: required input missing: %s", f))
}

cat("STEP 1: loading correlation results ...\n")
summ <- read.table(summary_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pairs_df <- read.table(pairs_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                        colClasses = c("character", "character", "numeric"))

vars <- summ$variable

# -------------------------------------------------------------------
# STEP 2: Correlation groups = connected components of the |r|>=thr graph
# -------------------------------------------------------------------
cat("STEP 2: computing correlation groups (connected components) ...\n")

adj <- setNames(vector("list", length(vars)), vars)
for (v in vars) adj[[v]] <- character(0)
if (nrow(pairs_df) > 0) {
  for (i in seq_len(nrow(pairs_df))) {
    a <- pairs_df$variable_1[i]; b <- pairs_df$variable_2[i]
    if (a %in% vars) adj[[a]] <- c(adj[[a]], b)
    if (b %in% vars) adj[[b]] <- c(adj[[b]], a)
  }
}

visited  <- setNames(rep(FALSE, length(vars)), vars)
group_id <- setNames(rep(NA_integer_, length(vars)), vars)
gid <- 0
for (v in vars) {
  if (visited[[v]]) next
  gid <- gid + 1
  queue <- v
  while (length(queue) > 0) {
    cur <- queue[1]; queue <- queue[-1]
    if (visited[[cur]]) next
    visited[[cur]] <- TRUE
    group_id[[cur]] <- gid
    queue <- c(queue, adj[[cur]])
  }
}
group_size <- table(group_id)
cat(sprintf("  %d correlation groups (%d singletons, %d groups with >1 variable)\n",
            length(group_size), sum(group_size == 1), sum(group_size > 1)))

# -------------------------------------------------------------------
# STEP 3: Category + one-line description per variable
# -------------------------------------------------------------------
cat("STEP 3: assigning categories and descriptions ...\n")

descriptions <- c(
  precip_annual         = "Annual precipitation total (mm)",
  dry_season_duration   = "Number of dry months/year (aet>ppt)",
  dry_season_intensity  = "Cumulative water deficit during dry months (mm) - aet-ppt sum",
  wet_season_duration   = "Number of wet months/year (aet<=ppt)",
  wet_season_intensity  = "Cumulative water surplus during wet months (mm) - ppt-aet sum",
  precip_contrast       = "Wet minus dry season water-balance intensity (mm)",
  cwd_annual            = "Annual climatic water deficit total (mm)",
  tmax_dry_season       = "Mean max temperature during dry months (degC)",
  soil_moisture_min     = "Minimum monthly soil moisture (mm), mean of annual minima",
  vpd_max_annual        = "Maximum monthly VPD (kPa), mean of annual maxima",
  BIO1  = "Annual mean temperature",
  BIO2  = "Mean diurnal temperature range (mean of monthly max-min)",
  BIO3  = "Isothermality (BIO2/BIO7 x100)",
  BIO4  = "Temperature seasonality (SD x100)",
  BIO5  = "Max temperature of warmest month",
  BIO6  = "Min temperature of coldest month",
  BIO7  = "Temperature annual range (BIO5-BIO6)",
  BIO8  = "Mean temperature of wettest quarter",
  BIO9  = "Mean temperature of driest quarter",
  BIO10 = "Mean temperature of warmest quarter",
  BIO11 = "Mean temperature of coldest quarter",
  BIO12 = "Annual precipitation",
  BIO13 = "Precipitation of wettest month",
  BIO14 = "Precipitation of driest month",
  BIO15 = "Precipitation seasonality (coefficient of variation)",
  BIO16 = "Precipitation of wettest quarter",
  BIO17 = "Precipitation of driest quarter",
  BIO18 = "Precipitation of warmest quarter",
  BIO19 = "Precipitation of coldest quarter",
  elevation = "Elevation (m), Copernicus DEM GLO-30",
  tri       = "Terrain Ruggedness Index, Copernicus DEM GLO-30"
)

categorize <- function(v) {
  if (grepl("^BIO[0-9]+$", v)) return("BIO1-19 (biovars, TerraClimate climatology)")
  if (v %in% c("elevation", "tri")) return("Topographic (Copernicus DEM)")
  "Custom TerraClimate-derived"
}

sel <- data.frame(
  variable    = vars,
  category    = vapply(vars, categorize, character(1)),
  group_id    = as.integer(group_id[vars]),
  group_size  = as.integer(group_size[as.character(group_id[vars])]),
  cv_pct      = summ$cv_pct[match(vars, summ$variable)],
  n_correlated = summ$n_correlated[match(vars, summ$variable)],
  description = descriptions[vars],
  stringsAsFactors = FALSE
)
sel$low_variation <- sel$cv_pct < 5
sel <- sel[order(sel$category, sel$group_id, sel$variable), ]

out_summ <- file.path(corr_dir, "variable_selection_summary.tsv")
write.table(sel, out_summ, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_summ))

# -------------------------------------------------------------------
# STEP 4: Editable template (uncomment to KEEP), grouped by category
# -------------------------------------------------------------------
cat("STEP 4: writing editable selection template ...\n")

out_template <- file.path(corr_dir, "variable_selection_template.txt")
con <- file(out_template, "w")
writeLines(c(
  "# 08.1b.5 -- Variable selection template",
  "# ========================================",
  "# Instructions:",
  "#   1. Review the correlation heatmap: plots/env_correlation_heatmap.png",
  "#   2. Review correlated_pairs_table.tsv and variable_selection_summary.tsv",
  "#   3. For each correlation group (variables inter-correlated |r|>=0.70),",
  "#      keep ONE representative based on ecological relevance and/or lowest",
  "#      redundancy - the group_id/group_size columns identify these groups.",
  "#   4. Variables flagged low_variation=TRUE (CV<5%) are weak candidates -",
  "#      review before keeping.",
  "#   5. Edit this file: uncomment (remove '# ') the variables you want to KEEP.",
  "#   6. Save in place. This file feeds 08.1b.6 (PCA) and 08.1b.7 (RDA).",
  "#",
  "# No variable is pre-selected here - the PCA (08.1b.6) that would identify",
  "# real ecological gradients runs AFTER this manual step, so no gradient-",
  "# based recommendation is made at this stage.",
  "#",
  sprintf("# Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("# n = %d sites | %d candidate variables | threshold |r| >= 0.70", 19, nrow(sel)),
  "#"
), con)

for (cat_name in unique(sel$category)) {
  writeLines(sprintf("# --- %s ---", cat_name), con)
  block <- sel[sel$category == cat_name, ]
  for (i in seq_len(nrow(block))) {
    r <- block[i, ]
    flag <- if (r$low_variation) " [LOW VARIATION]" else ""
    line <- sprintf("# %s # group=%d (size=%d) | CV=%.1f%%%s | %s",
                    r$variable, r$group_id, r$group_size, r$cv_pct, flag, r$description)
    writeLines(line, con)
  }
  writeLines("#", con)
}
close(con)
cat(sprintf("  Wrote: %s\n", out_template))

cat(sprintf("\nDONE 08.1b.5-manual_selection_prep.R\n"))
cat(sprintf("  %d variables ready for manual review - edit %s\n",
            nrow(sel), out_template))
cat("  >>> STOP HERE. Wait for Julien's edited variable_selection_template.txt before running 08.1b.6. <<<\n")

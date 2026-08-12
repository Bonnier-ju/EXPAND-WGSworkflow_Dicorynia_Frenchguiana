#!/usr/bin/env Rscript
# 08.1b.8-final_export.R
#
# Final export of the refonte 08.1 pipeline: individual-level environmental
# matrix (input for all GEA methods, 08.2 onward) + a markdown summary
# report built from the REAL outputs of 08.1b.1-08.1b.7 (no fabricated
# numbers - every figure in the report is read from a prior step's output
# file; sections are left as "(not available)" if a file is missing rather
# than guessed).
#
# Input:
#   rda_vif/variables_final.txt
#   correlation_analysis/env_variables_merged_31vars.csv
#   correlation_analysis/variable_selection_summary.tsv
#   correlation_analysis/correlated_pairs_table.tsv
#   correlation_analysis/variable_selection_template.txt
#   pca/env_pca_eigenvalues.tsv
#   rda_vif/rda_forward_selection.tsv
#   rda_vif/rda_vif.tsv
#   08.2-genotype_matrix/sample_metadata.tsv   (excludes DBR_ADT, N=154)
#
# Output:
#   final/env_matrix_final_154ind.tsv
#   final/env_variable_selection_report.md
#
# Args: out_dir  meta_file

EXCLUDE_IDS <- c("DBR_ADT")

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
out_dir   <- if (length(args) >= 1) args[1] else stop("out_dir required")
meta_file <- if (length(args) >= 2) args[2] else stop("meta_file required")

corr_dir  <- file.path(out_dir, "correlation_analysis")
pca_dir   <- file.path(out_dir, "pca")
rda_dir   <- file.path(out_dir, "rda_vif")
final_dir <- file.path(out_dir, "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

final_vars_file <- file.path(rda_dir, "variables_final.txt")
env_merged_csv  <- file.path(corr_dir, "env_variables_merged_31vars.csv")

# ---- Guard clause: verify core inputs before running ----
for (f in c(final_vars_file, env_merged_csv, meta_file)) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    stop(sprintf("ERROR: required input missing or empty: %s", f))
  }
}

# -------------------------------------------------------------------
# STEP 1: Final individual-level environmental matrix
# -------------------------------------------------------------------
cat("STEP 1: building final individual-level environmental matrix ...\n")

final_vars <- readLines(final_vars_file)
final_vars <- final_vars[nchar(trimws(final_vars)) > 0]
cat(sprintf("  %d final variables: %s\n", length(final_vars), paste(final_vars, collapse = ", ")))

env_site <- read.csv(env_merged_csv, stringsAsFactors = FALSE, check.names = FALSE)
missing_vars <- setdiff(final_vars, colnames(env_site))
if (length(missing_vars) > 0) {
  stop(sprintf("Final variable(s) not found in %s: %s", env_merged_csv, paste(missing_vars, collapse = ", ")))
}
env_site <- env_site[, c("site", final_vars)]

meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[!(meta$sample_id %in% EXCLUDE_IDS), ]

env_ind <- merge(meta[, c("sample_id", "site")], env_site, by = "site", all.x = FALSE)
env_ind <- env_ind[!is.na(env_ind$site), c("sample_id", "site", final_vars)]
env_ind <- env_ind[order(env_ind$sample_id), ]

if (nrow(env_ind) != 154) {
  warning(sprintf("Expected 154 individuals in the final matrix, got %d - check %s and %s",
                  nrow(env_ind), meta_file, env_merged_csv))
}

out_matrix <- file.path(final_dir, sprintf("env_matrix_final_%dind.tsv", nrow(env_ind)))
write.table(env_ind, out_matrix, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s (%d individuals x %d variables)\n", out_matrix, nrow(env_ind), length(final_vars)))

# -------------------------------------------------------------------
# STEP 2: Summary report - built only from files that actually exist
# -------------------------------------------------------------------
cat("STEP 2: building summary report ...\n")

read_if_exists <- function(path, ...) if (file.exists(path)) tryCatch(read.table(path, ...), error = function(e) NULL) else NULL

sel_summary  <- read_if_exists(file.path(corr_dir, "variable_selection_summary.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pairs_tbl    <- read_if_exists(file.path(corr_dir, "correlated_pairs_table.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
eig_tbl      <- read_if_exists(file.path(pca_dir, "env_pca_eigenvalues.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
fwd_tbl      <- read_if_exists(file.path(rda_dir, "rda_forward_selection.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
vif_tbl      <- read_if_exists(file.path(rda_dir, "rda_vif.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
template_lines <- if (file.exists(file.path(corr_dir, "variable_selection_template.txt"))) readLines(file.path(corr_dir, "variable_selection_template.txt")) else character(0)
kept_manual <- sub("^([A-Za-z0-9_]+).*", "\\1", grep("^[A-Za-z0-9_]+\\s*#\\s*group=", template_lines, value = TRUE))

na_or <- function(x, fallback = "_(not available)_") if (is.null(x)) fallback else x

md <- c(
"# Environmental variable selection report -- 08.1b (TerraClimate refonte)",
"",
sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
"",
"## 1. Data sources and temporal interval",
"",
"- Climate: TerraClimate, monthly, 1970-2000 (31-year window), extracted via curl + gdallocationinfo (08.1b.1).",
"- Topography: Copernicus DEM GLO-30 (~30 m), tiles already available in the project (08.1b.3).",
"- No CHELSA/ENVIREM/manual habitat variables in this refonte - single-source climate + single-source topography (cf. CLAUDE.md §4).",
"",
"## 2. Candidate variables (31 total: 10 custom + 19 BIO + 2 topographic)",
"",
if (!is.null(sel_summary)) {
  c("| variable | category | CV% | n_correlated |",
    "|---|---|---|---|",
    sprintf("| %s | %s | %.1f | %d |", sel_summary$variable, sel_summary$category,
            sel_summary$cv_pct, sel_summary$n_correlated))
} else na_or(NULL),
"",
"## 3. Correlation analysis",
"",
sprintf("- %d variable pairs flagged at |Spearman r| >= 0.70 (see correlation_analysis/correlated_pairs_table.tsv and plots/env_correlation_heatmap.png).",
        if (!is.null(pairs_tbl)) nrow(pairs_tbl) else NA),
"",
"## 4. Manual selection (Julien, 08.1b.5/08.1b.6)",
"",
sprintf("- %d variables kept from the 31 candidates: %s",
        length(kept_manual), paste(kept_manual, collapse = ", ")),
"",
"## 5. Environmental PCA (08.1b.6)",
"",
if (!is.null(eig_tbl)) {
  c("| PC | variance % | cumulative % |",
    "|---|---|---|",
    sprintf("| %s | %.1f | %.1f |", eig_tbl$PC, eig_tbl$variance_pct, eig_tbl$cumulative_pct))
} else na_or(NULL),
"",
"See plots/env_pca_biplot_PC1PC2.png, plots/env_pca_biplot_PC1PC3.png, pca/env_pca_top_loadings.tsv.",
"",
"## 6. RDA forward selection (08.1b.7)",
"",
if (!is.null(fwd_tbl)) {
  c("| variable | R2adj (cumulative) |",
    "|---|---|",
    sprintf("| %s | %.4f |", fwd_tbl$variable, fwd_tbl$R2adj_cumul))
} else na_or(NULL),
"",
"## 7. VIF screening (08.1b.7, threshold < 10)",
"",
if (!is.null(vif_tbl)) {
  last_iter <- max(vif_tbl$iteration)
  vt <- vif_tbl[vif_tbl$iteration == last_iter, ]
  c(sprintf("Final iteration (%d):", last_iter),
    "| variable | VIF |", "|---|---|",
    sprintf("| %s | %.2f |", vt$variable, vt$vif))
} else na_or(NULL),
"",
"## 8. Final variable list",
"",
sprintf("- **%s** (N = %d individuals)", paste(final_vars, collapse = ", "), nrow(env_ind)),
sprintf("- Matrix: `%s`", basename(out_matrix)),
"",
"## 9. Notes",
"",
"- MC_87 and MC_88 are 1.6 km apart and fall in the same TerraClimate pixel (~4 km resolution) - identical TerraClimate values for both sites, documented and accepted (per the original 08.1 plan).",
sprintf("- DBR_ADT excluded from all individual-level tables (confirmed sampling error, cf. CLAUDE.md §4) - N = %d, not 155.", nrow(env_ind)),
"- No structure conditioning in the RDA forward selection step (variable selection, not candidate detection) - structure control is mandatory downstream, per method, in 08.2-08.9 (cf. CLAUDE.md §10-11).",
""
)

out_report <- file.path(final_dir, "env_variable_selection_report.md")
writeLines(md, out_report)
cat(sprintf("  Wrote: %s\n", out_report))

cat(sprintf("\nDONE 08.1b.8-final_export.R\n"))
cat(sprintf("  Final matrix: %s\n", out_matrix))
cat(sprintf("  Report      : %s\n", out_report))

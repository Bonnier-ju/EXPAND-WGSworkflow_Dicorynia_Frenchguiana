#!/usr/bin/env Rscript
# 08.2-lfmm2.R
# Latent Factor Mixed Model (LFMM2) genotype-environment association analysis.
#
# LFMM2 simultaneously estimates the effects of environmental variables (fixed
# effects) and K latent factors capturing unmeasured confounders (neutral
# population structure, cryptic relatedness). Because latent factors are
# estimated jointly from genotype and environment, effect size estimates have
# minimal confounding bias compared to methods that correct for structure
# post-hoc (Caye et al. 2019, Mol Ecol Res; Frichot & Francois 2015, MEE).
#
# Inputs:
#   genotype_012_imputed.tsv  — 155 ind x N SNPs, 0/1/2 dosage (imputed)
#   snp_ids.txt               — SNP identifiers scaffold:pos:ref:alt
#   sample_metadata.tsv       — sample_id, site, lat, long
#   variables_final.txt       — final env variables from 08.1d RDA+VIF
#   *_env_per_site.csv        — site-level environmental CSVs (08.0.x)
#
# Outputs (per environmental variable):
#   output/lfmm2_candidates_<VAR>.tsv  — candidate SNPs (scaffold, pos, p, q)
#   plots/lfmm2_pval_hist_<VAR>.png    — p-value distribution (raw + calibrated)
#   plots/lfmm2_manhattan_<VAR>.png    — Manhattan plot
#   output/lfmm2_candidates_all.tsv    — all candidates merged across variables
#   output/lfmm2_summary.tsv          — N candidates per variable

suppressPackageStartupMessages({
  library(LEA)
  library(qvalue)
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# STEP 1: Parse arguments
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 12) {
  stop("Usage: 08.2-lfmm2.R geno_matrix snp_ids sample_meta vars_final chelsa_csv envirem_csv manual_csv terra_csv out_dir plots_dir K FDR")
}
geno_file    <- args[1]
snp_ids_file <- args[2]
meta_file    <- args[3]
vars_file    <- args[4]
chelsa_csv   <- args[5]
envirem_csv  <- args[6]
manual_csv   <- args[7]
terra_csv    <- args[8]
out_dir      <- args[9]
plots_dir    <- args[10]
K            <- as.integer(args[11])
fdr_level    <- as.numeric(args[12])

cat(sprintf("K = %d | FDR = %.2f\n", K, fdr_level))

# -------------------------------------------------------------------
# STEP 2: Load final environmental variables
# -------------------------------------------------------------------
cat("Reading selected variables:", vars_file, "\n")
var_lines <- readLines(vars_file)
selected_vars <- var_lines[!grepl("^#", var_lines) & nchar(trimws(var_lines)) > 0]
selected_vars <- trimws(selected_vars)
cat(sprintf("INFO: %d final variables: %s\n", length(selected_vars),
            paste(selected_vars, collapse = ", ")))

# -------------------------------------------------------------------
# STEP 3: Build site-level environmental matrix
# -------------------------------------------------------------------
cat("Loading environmental data...\n")

chelsa  <- read.csv(chelsa_csv,  stringsAsFactors = FALSE, check.names = FALSE)
envirem <- read.csv(envirem_csv, stringsAsFactors = FALSE, check.names = FALSE)
manual  <- read.csv(manual_csv,  stringsAsFactors = FALSE, check.names = FALSE)
terra   <- read.csv(terra_csv,   stringsAsFactors = FALSE, check.names = FALSE)

env_site <- chelsa[, c("site", intersect(selected_vars, names(chelsa)))]

env_add <- function(env_site, src, selected_vars) {
  cols <- intersect(selected_vars, names(src))
  if (length(cols) == 0) return(env_site)
  merge(env_site, src[, c("site", cols)], by = "site", all.x = TRUE)
}
env_site <- env_add(env_site, envirem, selected_vars)
env_site <- env_add(env_site, manual,  selected_vars)
env_site <- env_add(env_site, terra,   selected_vars)

missing_vars <- setdiff(selected_vars, names(env_site))
if (length(missing_vars) > 0)
  warning("Variables not found in any env CSV: ", paste(missing_vars, collapse = ", "))

cat(sprintf("INFO: %d sites with environmental data\n", nrow(env_site)))

# -------------------------------------------------------------------
# STEP 4: Assign env values to individuals via site
# -------------------------------------------------------------------
cat("Loading sample metadata:", meta_file, "\n")
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)

ind_env <- merge(meta[, c("sample_id", "site")],
                 env_site[, c("site", selected_vars)],
                 by = "site", all.x = FALSE)
ind_env <- ind_env[complete.cases(ind_env), ]
cat(sprintf("INFO: %d individuals with complete environmental data\n", nrow(ind_env)))

# -------------------------------------------------------------------
# STEP 5: Load genotype matrix
# -------------------------------------------------------------------
cat("Loading genotype matrix:", geno_file, "\n")
cat("  (this may take several minutes for 433k SNPs)\n")
geno_raw <- read.table(geno_file, header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       row.names = 1)
cat(sprintf("INFO: %d individuals x %d SNPs loaded\n", nrow(geno_raw), ncol(geno_raw)))

# Keep only individuals present in both matrices, maintain same order
common_ind <- intersect(rownames(geno_raw), ind_env$sample_id)
cat(sprintf("INFO: %d individuals in common\n", length(common_ind)))

geno_mat <- as.matrix(geno_raw[common_ind, ])
ind_env  <- ind_env[match(common_ind, ind_env$sample_id), ]
env_mat  <- as.matrix(ind_env[, selected_vars])
rownames(env_mat) <- common_ind

cat(sprintf("INFO: genotype matrix: %d x %d\n", nrow(geno_mat), ncol(geno_mat)))
cat(sprintf("INFO: env matrix:      %d x %d\n", nrow(env_mat),  ncol(env_mat)))

# -------------------------------------------------------------------
# STEP 6: Load and parse SNP IDs
# -------------------------------------------------------------------
cat("Loading SNP IDs:", snp_ids_file, "\n")
snp_id_raw <- read.table(snp_ids_file, header = TRUE, stringsAsFactors = FALSE)[[1]]

# Parse scaffold:position:ref:alt
snp_parts <- strsplit(snp_id_raw, ":", fixed = TRUE)
snp_df <- data.frame(
  snp_id   = snp_id_raw,
  scaffold = sapply(snp_parts, `[`, 1),
  position = as.integer(sapply(snp_parts, `[`, 2)),
  ref      = sapply(snp_parts, `[`, 3),
  alt      = sapply(snp_parts, `[`, 4),
  stringsAsFactors = FALSE
)
cat(sprintf("INFO: %d SNPs parsed\n", nrow(snp_df)))

# Scaffold ordering for Manhattan plots (numeric sort)
scaffold_nums <- as.integer(gsub(".*_(\\d+)$", "\\1", snp_df$scaffold))
scaffold_order <- unique(snp_df$scaffold[order(scaffold_nums)])

# Compute cumulative position offset per scaffold
scaffold_sizes <- tapply(snp_df$position, snp_df$scaffold, max)
scaffold_offset <- c(0, cumsum(as.numeric(scaffold_sizes[scaffold_order[-length(scaffold_order)]])))
names(scaffold_offset) <- scaffold_order
snp_df$cum_pos <- snp_df$position + scaffold_offset[snp_df$scaffold]

# -------------------------------------------------------------------
# STEP 7: Fit LFMM2 model
# -------------------------------------------------------------------
cat(sprintf("\nSTEP 7: Fitting LFMM2 (K = %d)...\n", K))
cat("  This may take 30-90 minutes on 433k SNPs.\n")
mod_lfmm2 <- lfmm2(input = geno_mat, env = env_mat, K = K)
cat("  LFMM2 model fitted.\n")

# -------------------------------------------------------------------
# STEP 8: Association test (per variable, genomic control)
# -------------------------------------------------------------------
cat("STEP 8: LFMM2 association test (genomic.control = TRUE)...\n")
test_raw  <- lfmm2.test(object = mod_lfmm2, input = geno_mat, env = env_mat,
                        full = FALSE, genomic.control = FALSE)
test_cal  <- lfmm2.test(object = mod_lfmm2, input = geno_mat, env = env_mat,
                        full = FALSE, genomic.control = TRUE)

# pvalues: matrix [n_vars x n_snps]
pval_raw <- test_raw$pvalues
pval_cal <- test_cal$pvalues

# Ensure variable names as row names
if (is.null(rownames(pval_cal))) rownames(pval_cal) <- selected_vars
if (is.null(rownames(pval_raw))) rownames(pval_raw) <- selected_vars

cat(sprintf("  p-value matrix: %d variables x %d SNPs\n",
            nrow(pval_cal), ncol(pval_cal)))

# -------------------------------------------------------------------
# STEP 9: Per-variable FDR, candidates, plots
# -------------------------------------------------------------------
cat(sprintf("STEP 9: Per-variable FDR (threshold = %.2f)...\n", fdr_level))

all_candidates <- list()
summary_rows   <- list()

# Color palette for scaffolds in Manhattan plots
scaffold_colors <- rep(c("#2166AC", "#4DAC26"), length.out = length(scaffold_order))
names(scaffold_colors) <- scaffold_order

for (v in selected_vars) {
  cat(sprintf("  Processing: %s\n", v))

  pv_raw_v <- as.numeric(pval_raw[v, ])
  pv_cal_v <- as.numeric(pval_cal[v, ])

  # --- p-value histogram (raw vs calibrated) ---
  png(file.path(plots_dir, sprintf("lfmm2_pval_hist_%s.png", v)),
      width = 1400, height = 600, res = 150)
  par(mfrow = c(1, 2),
      bg = "white",
      mar = c(4, 4, 3, 1))
  hist(pv_raw_v, breaks = 40, col = "#4393C3", border = "white",
       main = paste("Raw p-values —", v),
       xlab = "p-value", ylab = "Frequency")
  hist(pv_cal_v, breaks = 40, col = "#D6604D", border = "white",
       main = paste("Calibrated p-values —", v),
       xlab = "p-value", ylab = "Frequency")
  dev.off()

  # --- FDR via q-values ---
  qv <- tryCatch(
    qvalue::qvalue(pv_cal_v, fdr.level = fdr_level),
    error = function(e) {
      warning(sprintf("  WARN: qvalue failed for %s (%s) — using BH correction\n", v, e$message))
      NULL
    }
  )

  if (!is.null(qv)) {
    candidate_idx <- which(qv$qvalues < fdr_level)
    q_vals        <- qv$qvalues
  } else {
    bh <- p.adjust(pv_cal_v, method = "BH")
    candidate_idx <- which(bh < fdr_level)
    q_vals        <- bh
  }

  n_cand <- length(candidate_idx)
  cat(sprintf("    %d candidate SNPs at FDR < %.2f\n", n_cand, fdr_level))

  # --- Candidate table ---
  if (n_cand > 0) {
    cand_df <- data.frame(
      snp_id   = snp_df$snp_id[candidate_idx],
      scaffold = snp_df$scaffold[candidate_idx],
      position = snp_df$position[candidate_idx],
      ref      = snp_df$ref[candidate_idx],
      alt      = snp_df$alt[candidate_idx],
      variable = v,
      pval_raw = pv_raw_v[candidate_idx],
      pval_cal = pv_cal_v[candidate_idx],
      qval     = q_vals[candidate_idx],
      stringsAsFactors = FALSE
    )
    cand_df <- cand_df[order(cand_df$qval), ]
    write.table(cand_df,
                file.path(out_dir, sprintf("lfmm2_candidates_%s.tsv", v)),
                sep = "\t", row.names = FALSE, quote = FALSE)
    all_candidates[[v]] <- cand_df
  } else {
    cat(sprintf("    No candidates for %s\n", v))
  }

  # --- Manhattan plot ---
  plot_df <- data.frame(
    cum_pos  = snp_df$cum_pos,
    scaffold = snp_df$scaffold,
    neg_log_p = -log10(pv_cal_v),
    is_cand  = seq_along(pv_cal_v) %in% candidate_idx,
    stringsAsFactors = FALSE
  )

  # Axis: midpoint per scaffold
  axis_df <- aggregate(cum_pos ~ scaffold, data = plot_df, FUN = mean)
  axis_df <- axis_df[match(scaffold_order, axis_df$scaffold), ]

  p_thresh <- -log10(max(pv_cal_v[candidate_idx], na.rm = TRUE))
  if (n_cand == 0) p_thresh <- NA

  g <- ggplot(plot_df, aes(x = cum_pos, y = neg_log_p, color = scaffold)) +
    geom_point(size = 0.3, alpha = 0.4, show.legend = FALSE) +
    scale_color_manual(values = scaffold_colors) +
    {if (!is.na(p_thresh))
       geom_hline(yintercept = p_thresh, linetype = "dashed", color = "#D6604D", linewidth = 0.5)
    } +
    {if (n_cand > 0)
       geom_point(data = subset(plot_df, is_cand),
                  aes(x = cum_pos, y = neg_log_p),
                  color = "#D6604D", size = 1.2, alpha = 0.8)
    } +
    scale_x_continuous(
      breaks = axis_df$cum_pos,
      labels = gsub("Super-Scaffold_", "SS", axis_df$scaffold)
    ) +
    labs(
      title = sprintf("LFMM2 Manhattan — %s  (%d candidates, FDR < %.2f)",
                      v, n_cand, fdr_level),
      x = "Scaffold",
      y = expression(-log[10](p))
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor  = element_blank()
    )

  ggsave(file.path(plots_dir, sprintf("lfmm2_manhattan_%s.png", v)),
         g, width = 14, height = 4, dpi = 150)
  ggsave(file.path(plots_dir, sprintf("lfmm2_manhattan_%s.pdf", v)),
         g, width = 14, height = 4)

  summary_rows[[v]] <- data.frame(variable = v, n_candidates = n_cand,
                                  stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------
# STEP 10: Cross-variable summary
# -------------------------------------------------------------------
cat("STEP 10: Writing summary outputs...\n")

# Merged candidate table
if (length(all_candidates) > 0) {
  all_cand_df <- do.call(rbind, all_candidates)
  rownames(all_cand_df) <- NULL
  write.table(all_cand_df,
              file.path(out_dir, "lfmm2_candidates_all.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  Total candidates (all variables): %d unique SNPs\n",
              length(unique(all_cand_df$snp_id))))
} else {
  cat("  No candidates across all variables.\n")
}

# Per-variable summary
summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL
write.table(summary_df,
            file.path(out_dir, "lfmm2_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== LFMM2 SUMMARY ===\n")
print(summary_df)
cat(sprintf("\nK = %d | FDR = %.2f\n", K, fdr_level))
cat(sprintf("Individuals: %d | SNPs: %d | Variables: %d\n",
            nrow(geno_mat), ncol(geno_mat), length(selected_vars)))

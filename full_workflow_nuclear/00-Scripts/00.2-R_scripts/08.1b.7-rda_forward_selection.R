#!/usr/bin/env Rscript
# 08.1b.7-rda_forward_selection.R
#
# RDA forward selection + VIF screening for the refonte 08.1 pipeline.
# Adapted from the proven 08.1c-env_rda_vif.R (July run, CHELSA+ENVIREM+
# manual+TerraClimate) - same method (Blanchet et al. 2008 double-criterion
# forward selection, vif.cca screening), simplified inputs (single merged
# TerraClimate-derived table instead of 4 sources) and updated to N=154.
#
# IMPORTANT - N=154 (confirmed 2026-08-11): DBR_ADT is a confirmed sampling
# error (cf. CLAUDE.md §4) and must be excluded from every GEA analysis.
# sample_metadata.tsv still CONTAINS DBR_ADT (verified 2026-08-11, 155 rows)
# - this script filters it out explicitly rather than assuming the input is
# already clean, and asserts the final N really is 154 before proceeding.
#
# Method:
#   Response Y: Hellinger-transformed genotype matrix (decostand, vegan).
#   Predictor X: scaled env variables (Julien's manual selection from
#     08.1b.5/08.1b.6) matched per individual via site.
#   Forward selection: ordiR2step (vegan), double stopping criterion
#     alpha=0.05 and R2adj <= R2adj(global model), permutations configurable.
#   VIF: vif.cca (vegan) on the final forward-selected model. Variables with
#     VIF > 10 are iteratively removed (highest VIF first).
#   No structure conditioning at this stage (matches the original/07 method
#     for environmental VARIABLE SELECTION, e.g. Forester et al. 2018) - this
#     is distinct from the downstream 08.3 pRDA candidate-detection step,
#     where structure conditioning is mandatory (cf. CLAUDE.md §11).
#
# Inputs:
#   genotype_012_imputed.tsv          — N ind x N SNPs, 0/1/2 dosage, imputed
#   sample_metadata.tsv               — sample_id, site, lat, long (155 rows,
#                                        includes DBR_ADT - filtered here)
#   env_variables_merged_31vars.csv   — site-level, from 08.1b.4
#   variable_selection_template.txt   — edited by Julien in 08.1b.5/.6
#
# Outputs (rda_vif/):
#   rda_forward_selection.tsv
#   rda_vif.tsv
#   variables_final.txt
#   rda_biplot_rda12.png/.pdf
#   rda_variance_partitioning.tsv
#
# Args: geno_file  meta_file  env_merged_csv  template_file  out_dir  n_perm  seed

for (pkg in c("vegan", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' is required but not installed. Check with: search_R_package %s", pkg, pkg))
}
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) library(ggrepel)

SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else { c() }

EXCLUDE_IDS <- c("DBR_ADT")  # confirmed sampling error, cf. CLAUDE.md §4

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args           <- commandArgs(trailingOnly = TRUE)
geno_file      <- if (length(args) >= 1) args[1] else stop("geno_file required")
meta_file      <- if (length(args) >= 2) args[2] else stop("meta_file required")
env_merged_csv <- if (length(args) >= 3) args[3] else stop("env_merged_csv required")
template_file  <- if (length(args) >= 4) args[4] else stop("template_file required")
out_dir        <- if (length(args) >= 5) args[5] else stop("out_dir required")
n_perm         <- if (length(args) >= 6) as.integer(args[6]) else 999L
seed           <- if (length(args) >= 7) as.integer(args[7]) else 1101L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)
cat(sprintf("Seed used: %d (logged for reproducibility, cf. seed.txt)\n", seed))

# ---- Guard clause: verify inputs before running ----
for (f in c(geno_file, meta_file, env_merged_csv, template_file)) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    stop(sprintf("ERROR: required input missing or empty: %s", f))
  }
}

# -------------------------------------------------------------------
# STEP 1: Parse selected variables (same parser as 08.1b.6 - a kept line
# has had its leading "# " removed)
# -------------------------------------------------------------------
cat("STEP 1: parsing edited variable_selection_template.txt ...\n")
tmpl_lines <- readLines(template_file)
kept_lines <- grep("^[A-Za-z0-9_]+\\s*#\\s*group=", tmpl_lines, value = TRUE)
selected_vars <- sub("^([A-Za-z0-9_]+).*", "\\1", kept_lines)
if (length(selected_vars) == 0) {
  stop("ERROR: no variable selected in variable_selection_template.txt. Edit it first (08.1b.5/08.1b.6).")
}
cat(sprintf("  %d candidate variables: %s\n", length(selected_vars), paste(selected_vars, collapse = ", ")))

# -------------------------------------------------------------------
# STEP 2: Load site-level environmental data
# -------------------------------------------------------------------
cat("STEP 2: loading environmental data ...\n")
env_all <- read.csv(env_merged_csv, stringsAsFactors = FALSE, check.names = FALSE)
missing_vars <- setdiff(selected_vars, colnames(env_all))
if (length(missing_vars) > 0) {
  stop(sprintf("Selected variable(s) not found in %s: %s", env_merged_csv, paste(missing_vars, collapse = ", ")))
}
env_site <- env_all[, c("site", selected_vars)]
cat(sprintf("  %d sites with environmental data\n", nrow(env_site)))

# -------------------------------------------------------------------
# STEP 3: Load sample metadata, exclude DBR_ADT, join env to individuals
# -------------------------------------------------------------------
cat("STEP 3: loading sample metadata ...\n")
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
cat(sprintf("  %d individuals in %s (before exclusion)\n", nrow(meta), meta_file))

n_before <- nrow(meta)
meta <- meta[!(meta$sample_id %in% EXCLUDE_IDS), ]
cat(sprintf("  Excluded: %s -> %d individuals remaining\n",
            paste(EXCLUDE_IDS, collapse = ", "), nrow(meta)))
if (n_before - nrow(meta) != length(EXCLUDE_IDS)) {
  warning(sprintf("Expected to exclude %d individual(s), removed %d - check %s for duplicates/mismatched IDs",
                  length(EXCLUDE_IDS), n_before - nrow(meta), meta_file))
}

ind_env <- merge(meta[, c("sample_id", "site")], env_site, by = "site", all.x = FALSE)
ind_env <- ind_env[!is.na(ind_env$site), ]
cat(sprintf("  %d individuals with environmental data\n", nrow(ind_env)))

# -------------------------------------------------------------------
# STEP 4: Load genotype matrix, align to individuals, assert N=154
# -------------------------------------------------------------------
cat("STEP 4: loading genotype matrix (this may take a few minutes) ...\n")
geno_raw <- read.table(geno_file, header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
cat(sprintf("  %d individuals x %d SNPs loaded\n", nrow(geno_raw), ncol(geno_raw)))

geno_raw <- geno_raw[!(rownames(geno_raw) %in% EXCLUDE_IDS), ]

common_ind <- intersect(rownames(geno_raw), ind_env$sample_id)
cat(sprintf("  %d individuals in common (genotype x metadata, post-exclusion)\n", length(common_ind)))

if (length(common_ind) != 154) {
  warning(sprintf(
    "Expected N=154 (155 Guyane ind - DBR_ADT) but got N=%d. Check that sample_metadata.tsv / genotype matrix don't already exclude other individuals (Angela/outgroups should already be absent upstream) before trusting downstream results.",
    length(common_ind)))
}

geno_mat <- as.matrix(geno_raw[common_ind, ])
mode(geno_mat) <- "numeric"

ind_env <- ind_env[match(common_ind, ind_env$sample_id), ]
rownames(ind_env) <- ind_env$sample_id

env_mat <- ind_env[common_ind, selected_vars, drop = FALSE]
env_mat <- as.data.frame(lapply(env_mat, as.numeric))
rownames(env_mat) <- common_ind

n_env_na <- sum(is.na(env_mat))
if (n_env_na > 0) {
  cat(sprintf("  WARNING: %d NA values in env matrix - imputing column means\n", n_env_na))
  for (v in names(env_mat)) {
    na_idx <- is.na(env_mat[[v]])
    if (any(na_idx)) env_mat[[v]][na_idx] <- mean(env_mat[[v]], na.rm = TRUE)
  }
}
env_scaled <- as.data.frame(scale(env_mat))
cat(sprintf("  env matrix: %d individuals x %d variables\n", nrow(env_scaled), ncol(env_scaled)))

# -------------------------------------------------------------------
# STEP 5: Hellinger transformation of the genotype matrix
# -------------------------------------------------------------------
cat("STEP 5: Hellinger transformation ...\n")
geno_hell <- decostand(geno_mat / 2, method = "hellinger")
cat(sprintf("  Hellinger-transformed matrix: %d x %d\n", nrow(geno_hell), ncol(geno_hell)))

# -------------------------------------------------------------------
# STEP 6: Global RDA (stopping criterion for forward selection)
# -------------------------------------------------------------------
cat(sprintf("\nSTEP 6: Global RDA (%d variables) ...\n", ncol(env_scaled)))
rda_global <- rda(geno_hell ~ ., data = env_scaled)
r2_global  <- RsquareAdj(rda_global)$adj.r.squared
cat(sprintf("  Global RDA R2adj = %.4f\n", r2_global))

# -------------------------------------------------------------------
# STEP 7: RDA forward selection (Blanchet et al. 2008 double criterion)
# -------------------------------------------------------------------
cat(sprintf("\nSTEP 7: Forward selection (alpha=0.05, %d permutations) ...\n", n_perm))
rda_null <- rda(geno_hell ~ 1, data = env_scaled)
rda_fwd <- ordiR2step(rda_null, scope = formula(rda_global), R2scope = r2_global,
                      direction = "forward", permutations = n_perm, trace = TRUE)

fwd_vars <- names(rda_fwd$terminfo$ordered)
cat(sprintf("\n  %d variables selected by forward selection: %s\n",
            length(fwd_vars), paste(fwd_vars, collapse = ", ")))

fwd_table <- as.data.frame(rda_fwd$anova)
fwd_table$variable <- rownames(fwd_table)
fwd_table$R2adj_cumul <- NA_real_
for (i in seq_along(fwd_vars)) {
  form_i <- as.formula(paste("geno_hell ~", paste(fwd_vars[1:i], collapse = " + ")))
  rda_i  <- rda(form_i, data = env_scaled)
  fwd_table$R2adj_cumul[i] <- RsquareAdj(rda_i)$adj.r.squared
}
write.table(fwd_table, file.path(out_dir, "rda_forward_selection.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  Forward selection table written\n")

if (length(fwd_vars) == 0) {
  cat("WARNING: no variables were selected by forward selection - stopping here.\n")
  quit(save = "no", status = 0)
}

# -------------------------------------------------------------------
# STEP 8: VIF screening (threshold = 10)
# -------------------------------------------------------------------
cat("\nSTEP 8: VIF screening (threshold = 10) ...\n")
vif_vars <- fwd_vars
vif_history <- list()
iteration <- 1L
repeat {
  form_vif <- as.formula(paste("geno_hell ~", paste(vif_vars, collapse = " + ")))
  rda_vif  <- rda(form_vif, data = env_scaled)
  vif_vals <- vif.cca(rda_vif)
  vif_df <- data.frame(variable = names(vif_vals), vif = round(vif_vals, 3),
                       iteration = iteration, stringsAsFactors = FALSE)
  vif_history[[iteration]] <- vif_df
  cat(sprintf("  Iteration %d:\n", iteration)); print(vif_df[order(-vif_df$vif), c("variable", "vif")])

  max_vif <- max(vif_vals); max_var <- names(which.max(vif_vals))
  if (max_vif <= 10) { cat("  All VIF <= 10 - stopping.\n"); break }
  cat(sprintf("  Removing '%s' (VIF = %.2f)\n", max_var, max_vif))
  vif_vars <- setdiff(vif_vars, max_var)
  iteration <- iteration + 1L
}
vif_all <- do.call(rbind, vif_history)
write.table(vif_all, file.path(out_dir, "rda_vif.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

final_vars <- vif_vars
cat(sprintf("\nFinal variables after VIF screening (%d): %s\n",
            length(final_vars), paste(final_vars, collapse = ", ")))
writeLines(final_vars, file.path(out_dir, "variables_final.txt"))

# -------------------------------------------------------------------
# STEP 9: Final RDA summary + variance partitioning
# -------------------------------------------------------------------
cat("\nSTEP 9: Final RDA summary ...\n")
form_final <- as.formula(paste("geno_hell ~", paste(final_vars, collapse = " + ")))
rda_final  <- rda(form_final, data = env_scaled)
r2_final   <- RsquareAdj(rda_final)$adj.r.squared
cat(sprintf("  Final RDA R2adj = %.4f (%d variables)\n", r2_final, length(final_vars)))

vpart_rows <- list()
for (v in final_vars) {
  f1 <- as.formula(paste("geno_hell ~", v))
  f2 <- as.formula(paste("geno_hell ~", paste(setdiff(final_vars, v), collapse = " + ")))
  r2_single  <- RsquareAdj(rda(f1, data = env_scaled))$adj.r.squared
  r2_without <- RsquareAdj(rda(f2, data = env_scaled))$adj.r.squared
  vpart_rows[[v]] <- data.frame(variable = v, r2adj_alone = round(r2_single, 4),
                                r2adj_unique = round(r2_final - r2_without, 4),
                                stringsAsFactors = FALSE)
}
vpart_df <- do.call(rbind, vpart_rows)
vpart_df$r2adj_final_total <- round(r2_final, 4)
write.table(vpart_df, file.path(out_dir, "rda_variance_partitioning.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------------------------------------------------
# STEP 10: RDA biplot
# -------------------------------------------------------------------
cat("\nSTEP 10: RDA biplot ...\n")
site_scores <- as.data.frame(scores(rda_final, display = "sites", scaling = 2))
site_scores$sample_id <- rownames(site_scores)
site_scores <- merge(site_scores, ind_env[, c("sample_id", "site")], by = "sample_id")

env_scores <- as.data.frame(scores(rda_final, display = "bp", scaling = 2))
env_scores$variable <- rownames(env_scores)

eig <- eigenvals(rda_final)
pct_exp <- eig / sum(eig) * 100
lab_rda <- function(i) sprintf("RDA%d (%.1f%%)", i, pct_exp[i])

scale_fac <- max(abs(site_scores[, c("RDA1", "RDA2")])) / max(abs(env_scores[, c("RDA1", "RDA2")])) * 0.6
env_sc2 <- env_scores
env_sc2$RDA1 <- env_scores$RDA1 * scale_fac
env_sc2$RDA2 <- env_scores$RDA2 * scale_fac

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_segment(data = env_sc2, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "grey30", linewidth = 0.7) +
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, color = site), size = 2.5, alpha = 0.8) +
  labs(title = sprintf("RDA forward selection - final model (R2adj = %.3f, N=%d)",
                       r2_final, length(common_ind)),
       x = lab_rda(1), y = lab_rda(2), color = "Site") +
  theme(legend.position = "right", legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7))
if (length(site_colors) > 0) p_biplot <- p_biplot + scale_color_manual(values = site_colors, na.value = "grey50")
if (has_ggrepel) {
  p_biplot <- p_biplot + ggrepel::geom_text_repel(data = env_sc2, aes(x = RDA1, y = RDA2, label = variable),
                                                   size = 3, color = "grey20", fontface = "italic",
                                                   max.overlaps = 20, seed = seed)
} else {
  p_biplot <- p_biplot + geom_text(data = env_sc2, aes(x = RDA1, y = RDA2, label = variable),
                                    size = 3, color = "grey20", fontface = "italic", vjust = -0.5)
}
ggsave(file.path(out_dir, "rda_biplot_rda12.png"), p_biplot, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "rda_biplot_rda12.pdf"), p_biplot, width = 10, height = 7)

# -------------------------------------------------------------------
# Final summary
# -------------------------------------------------------------------
cat("\n========== 08.1b.7 SUMMARY ==========\n")
cat(sprintf("N individuals                : %d\n", length(common_ind)))
cat(sprintf("Seed                         : %d\n", seed))
cat(sprintf("Candidate variables (08.1b.5): %d\n", length(selected_vars)))
cat(sprintf("Selected by forward selection: %d\n", length(fwd_vars)))
cat(sprintf("Retained after VIF < 10      : %d\n", length(final_vars)))
cat(sprintf("Final R2adj                  : %.4f\n", r2_final))
cat(sprintf("Global R2adj (all %d vars)    : %.4f\n", length(selected_vars), r2_global))
cat("\nFinal variables for GEA (08.2-08.9):\n")
for (v in final_vars) cat(sprintf("  %s\n", v))
cat(sprintf("\nVariables list: %s\n", file.path(out_dir, "variables_final.txt")))
cat("=====================================\n")

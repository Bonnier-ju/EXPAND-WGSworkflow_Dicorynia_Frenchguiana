#!/usr/bin/env Rscript
# 08.1c-env_rda_vif.R
# Final environmental variable selection for GEA:
#   (1) RDA forward selection — identifies variables significantly
#       associated with genome-wide SNP variation (Blanchet et al. 2008)
#   (2) VIF screening — removes residual multicollinearity among
#       selected variables (threshold VIF < 10)
#
# Inputs:
#   genotype_012_imputed.tsv  — 155 ind x N SNPs, 0/1/2 dosage, NAs imputed
#   sample_metadata.tsv       — individual_id, site, lat, long
#   chelsa_env_per_site.csv   — CHELSA BIO variables (temp in K/10, precip raw)
#   envirem_env_per_site.csv  — ENVIREM variables (physical units)
#   manual_variables_per_site.csv — elevation
#   variables_to_keep_template.txt — 11 pre-selected variables
#
# Outputs (rda_vif/):
#   rda_forward_selection.tsv — variables added, R2adj, F, p per step
#   rda_vif.tsv               — VIF values after forward selection
#   variables_final.txt       — final variable list for GEA (08.3-08.7)
#   rda_biplot_pc12.png/pdf   — RDA biplot (sites + env arrows)
#   rda_variance_partitioning.tsv — total R2adj explained by final set
#
# Method:
#   Response Y: Hellinger-transformed genotype matrix (decostand, vegan).
#   Predictor X: scaled env variables matched per individual via site.
#   Forward selection: ordiR2step (vegan), double stopping criterion
#     alpha=0.05 and R2adj <= R2adj(global model), 999 permutations.
#   VIF: vif.cca (vegan) on the final forward-selected model.
#     Variables with VIF > 10 are iteratively removed (highest VIF first).
#
# Reference: Forester et al. (2018) Mol Ecol Res; Blanchet et al. (2008) Ecology

# -------------------------------------------------------------------
# Check required packages
# -------------------------------------------------------------------
for (pkg in c("vegan", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package '%s' is required but not installed.\n", pkg),
         "Install with: install.packages('", pkg, "')\n", call. = FALSE)
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

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
geno_file   <- if (length(args) >= 1) args[1] else stop("geno_file required")
meta_file   <- if (length(args) >= 2) args[2] else stop("meta_file required")
chelsa_csv  <- if (length(args) >= 3) args[3] else stop("chelsa_csv required")
envirem_csv <- if (length(args) >= 4) args[4] else stop("envirem_csv required")
manual_csv  <- if (length(args) >= 5) args[5] else stop("manual_csv required")
vars_file   <- if (length(args) >= 6) args[6] else stop("vars_file required")
out_dir     <- if (length(args) >= 7) args[7] else stop("out_dir required")
n_perm      <- if (length(args) >= 8) as.integer(args[8]) else 999L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# -------------------------------------------------------------------
# STEP 1: Parse selected variables
# -------------------------------------------------------------------
cat("Parsing selected variables from:", vars_file, "\n")
lines <- readLines(vars_file)
selected_vars <- character(0)
for (ln in lines) {
  stripped <- trimws(ln)
  if (nchar(stripped) == 0 || grepl("^#", stripped)) next
  var_name <- trimws(sub("#.*", "", sub("\\s.*", "", stripped)))
  if (nchar(var_name) > 0) selected_vars <- c(selected_vars, var_name)
}
cat(sprintf("INFO: %d candidate variables: %s\n",
            length(selected_vars), paste(selected_vars, collapse = ", ")))

# -------------------------------------------------------------------
# STEP 2: Load and merge environmental data (site-level)
# -------------------------------------------------------------------
cat("Loading environmental data...\n")
chelsa  <- read.csv(chelsa_csv,  stringsAsFactors = FALSE, check.names = FALSE)
envirem <- read.csv(envirem_csv, stringsAsFactors = FALSE, check.names = FALSE)
manual  <- read.csv(manual_csv,  quote = "", stringsAsFactors = FALSE, check.names = FALSE)

# CHELSA scale corrections (see 08.0.1):
#   Absolute temperatures stored as Kelvin after /10 -> subtract 273.15
#   Precipitation stored raw (x10 of mm) -> divide by 10
temp_abs <- intersect(c("BIO1","BIO5","BIO6","BIO8","BIO9","BIO10","BIO11"), names(chelsa))
precip   <- intersect(paste0("BIO", 12:19), names(chelsa))
for (v in temp_abs) chelsa[[v]] <- chelsa[[v]] - 273.15
for (v in precip)   chelsa[[v]] <- chelsa[[v]] / 10

# Merge the three sources
env_site <- merge(
  chelsa[, c("site", intersect(selected_vars, names(chelsa)))],
  envirem[, c("site", intersect(selected_vars, names(envirem)))],
  by = "site", all = TRUE
)
if ("elevation" %in% selected_vars) {
  env_site <- merge(env_site,
                    manual[, c("site", "elevation")],
                    by = "site", all.x = TRUE)
}
env_site <- env_site[env_site$site != "Cameroun_Benin", ]

cat(sprintf("INFO: %d sites with environmental data\n", nrow(env_site)))

# -------------------------------------------------------------------
# STEP 3: Load sample metadata and match individuals to env variables
# -------------------------------------------------------------------
cat("Loading sample metadata:", meta_file, "\n")
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)

# Join env to each individual via site
ind_env <- merge(meta[, c("individual_id", "site")],
                 env_site, by = "site", all.x = FALSE)
ind_env <- ind_env[!is.na(ind_env$site), ]
cat(sprintf("INFO: %d individuals with environmental data\n", nrow(ind_env)))

# -------------------------------------------------------------------
# STEP 4: Load genotype matrix and align to individuals with env data
# -------------------------------------------------------------------
cat("Loading genotype matrix:", geno_file, "\n")
cat("  (this may take a few minutes for large SNP sets)\n")

geno_raw <- read.table(geno_file, header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       row.names = 1)
cat(sprintf("INFO: %d individuals x %d SNPs loaded\n", nrow(geno_raw), ncol(geno_raw)))

# Keep only individuals present in both matrices
common_ind <- intersect(rownames(geno_raw), ind_env$individual_id)
cat(sprintf("INFO: %d individuals in common\n", length(common_ind)))

geno_mat <- as.matrix(geno_raw[common_ind, ])
mode(geno_mat) <- "numeric"

ind_env <- ind_env[match(common_ind, ind_env$individual_id), ]
rownames(ind_env) <- ind_env$individual_id

env_vars_present <- intersect(selected_vars, names(ind_env))
missing_vars <- setdiff(selected_vars, names(ind_env))
if (length(missing_vars) > 0)
  warning("Variables not found in merged env data: ", paste(missing_vars, collapse = ", "))

env_mat <- ind_env[common_ind, env_vars_present, drop = FALSE]
env_mat <- as.data.frame(lapply(env_mat, as.numeric))
rownames(env_mat) <- common_ind

# Check for missing env values
n_env_na <- sum(is.na(env_mat))
if (n_env_na > 0) {
  cat(sprintf("WARNING: %d NA values in env matrix — imputing column means\n", n_env_na))
  for (v in names(env_mat)) {
    na_idx <- is.na(env_mat[[v]])
    if (any(na_idx)) env_mat[[v]][na_idx] <- mean(env_mat[[v]], na.rm = TRUE)
  }
}

# Scale environmental predictors (zero mean, unit variance)
env_scaled <- as.data.frame(scale(env_mat))

cat(sprintf("INFO: env matrix: %d individuals x %d variables\n",
            nrow(env_scaled), ncol(env_scaled)))

# -------------------------------------------------------------------
# STEP 5: Hellinger transformation of the genotype matrix
# Converts dosage counts to relative frequencies, then square-root-
# transforms — removes the arch effect and makes data suitable for RDA.
# Each individual's dosage vector is divided by its sum (total allele count),
# then square-rooted.
# -------------------------------------------------------------------
cat("Applying Hellinger transformation to genotype matrix...\n")
geno_hell <- decostand(geno_mat / 2, method = "hellinger")
cat(sprintf("INFO: Hellinger-transformed matrix: %d x %d\n",
            nrow(geno_hell), ncol(geno_hell)))

# -------------------------------------------------------------------
# STEP 6: Global RDA (full model with all candidate variables)
# Used to obtain the global R2adj as stopping criterion for forward selection
# -------------------------------------------------------------------
cat(sprintf("\nSTEP 6: Global RDA (%d variables)...\n", ncol(env_scaled)))
rda_global <- rda(geno_hell ~ ., data = env_scaled)
r2_global  <- RsquareAdj(rda_global)$adj.r.squared
cat(sprintf("INFO: Global RDA R2adj = %.4f\n", r2_global))

# -------------------------------------------------------------------
# STEP 7: RDA forward selection (Blanchet et al. 2008 double criterion)
# Stopping criteria: (1) p > alpha OR (2) R2adj exceeds global R2adj
# -------------------------------------------------------------------
cat(sprintf("\nSTEP 7: Forward selection (alpha=0.05, %d permutations)...\n", n_perm))

rda_null <- rda(geno_hell ~ 1, data = env_scaled)

rda_fwd <- ordiR2step(
  rda_null,
  scope   = formula(rda_global),
  R2scope = r2_global,
  direction = "forward",
  permutations = n_perm,
  trace   = TRUE
)

# Extract forward selection results
fwd_vars  <- names(rda_fwd$terminfo$ordered)
cat(sprintf("\nINFO: %d variables selected by forward selection: %s\n",
            length(fwd_vars), paste(fwd_vars, collapse = ", ")))

# Build step-by-step table
fwd_table <- as.data.frame(rda_fwd$anova)
fwd_table$variable <- rownames(fwd_table)
fwd_table$R2adj_cumul <- NA_real_

# Recompute cumulative R2adj at each step
for (i in seq_along(fwd_vars)) {
  form_i <- as.formula(paste("geno_hell ~", paste(fwd_vars[1:i], collapse = " + ")))
  rda_i  <- rda(form_i, data = env_scaled)
  fwd_table$R2adj_cumul[i] <- RsquareAdj(rda_i)$adj.r.squared
}

write.table(fwd_table, file = file.path(out_dir, "rda_forward_selection.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Forward selection table written\n")
print(fwd_table[, c("variable", "R2adj_cumul", "F", "Pr(>F)")])

if (length(fwd_vars) == 0) {
  cat("WARNING: no variables were selected by forward selection.\n")
  cat("         Check that the genotype matrix and env data are correctly formatted.\n")
  quit(save = "no", status = 0)
}

# -------------------------------------------------------------------
# STEP 8: VIF screening on forward-selected variables
# Iteratively removes the variable with the highest VIF until all VIF < 10.
# -------------------------------------------------------------------
cat("\nSTEP 8: VIF screening (threshold = 10)...\n")

vif_vars   <- fwd_vars
vif_history <- list()
iteration  <- 1L

repeat {
  form_vif <- as.formula(paste("geno_hell ~", paste(vif_vars, collapse = " + ")))
  rda_vif  <- rda(form_vif, data = env_scaled)
  vif_vals <- vif.cca(rda_vif)

  vif_df <- data.frame(
    variable  = names(vif_vals),
    vif       = round(vif_vals, 3),
    iteration = iteration,
    stringsAsFactors = FALSE
  )
  vif_history[[iteration]] <- vif_df
  cat(sprintf("  Iteration %d:\n", iteration))
  print(vif_df[order(-vif_df$vif), c("variable", "vif")])

  max_vif  <- max(vif_vals)
  max_var  <- names(which.max(vif_vals))

  if (max_vif <= 10) {
    cat(sprintf("  All VIF <= 10 — stopping.\n"))
    break
  }
  cat(sprintf("  Removing '%s' (VIF = %.2f)\n", max_var, max_vif))
  vif_vars <- setdiff(vif_vars, max_var)
  iteration <- iteration + 1L
}

# Write full VIF history
vif_all <- do.call(rbind, vif_history)
write.table(vif_all, file = file.path(out_dir, "rda_vif.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("VIF table written\n")

# Final variables retained
final_vars <- vif_vars
cat(sprintf("\nFinal variables after VIF screening (%d): %s\n",
            length(final_vars), paste(final_vars, collapse = ", ")))

# Write final variable list
writeLines(final_vars, file.path(out_dir, "variables_final.txt"))
cat("Final variable list written:", file.path(out_dir, "variables_final.txt"), "\n")

# -------------------------------------------------------------------
# STEP 9: Final RDA with retained variables — summary statistics
# -------------------------------------------------------------------
cat("\nSTEP 9: Final RDA summary...\n")
form_final <- as.formula(paste("geno_hell ~", paste(final_vars, collapse = " + ")))
rda_final  <- rda(form_final, data = env_scaled)
r2_final   <- RsquareAdj(rda_final)$adj.r.squared
cat(sprintf("INFO: Final RDA R2adj = %.4f (%d variables)\n", r2_final, length(final_vars)))

# Variance partitioning (each variable's marginal and individual contribution)
vpart_rows <- list()
for (v in final_vars) {
  f1 <- as.formula(paste("geno_hell ~", v))
  f2 <- as.formula(paste("geno_hell ~",
                         paste(setdiff(final_vars, v), collapse = " + ")))
  r2_single   <- RsquareAdj(rda(f1, data = env_scaled))$adj.r.squared
  r2_without  <- RsquareAdj(rda(f2, data = env_scaled))$adj.r.squared
  r2_unique   <- r2_final - r2_without
  vpart_rows[[v]] <- data.frame(
    variable      = v,
    r2adj_alone   = round(r2_single, 4),
    r2adj_unique  = round(r2_unique, 4),
    stringsAsFactors = FALSE
  )
}
vpart_df <- do.call(rbind, vpart_rows)
vpart_df$r2adj_final_total <- round(r2_final, 4)
write.table(vpart_df, file = file.path(out_dir, "rda_variance_partitioning.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Variance partitioning written\n")
print(vpart_df[order(-vpart_df$r2adj_unique), ])

# -------------------------------------------------------------------
# STEP 10: RDA biplot (sites + env arrows, RDA1 x RDA2)
# -------------------------------------------------------------------
cat("\nSTEP 10: RDA biplot...\n")

site_scores <- as.data.frame(scores(rda_final, display = "sites", scaling = 2))
site_scores$individual_id <- rownames(site_scores)
site_scores <- merge(site_scores, ind_env[, c("individual_id", "site")],
                     by = "individual_id")

env_scores  <- as.data.frame(scores(rda_final, display = "bp", scaling = 2))
env_scores$variable <- rownames(env_scores)

# Axis labels with % variance explained
eig     <- eigenvals(rda_final)
pct_exp <- eig / sum(eig) * 100
lab_rda <- function(i) sprintf("RDA%d (%.1f%%)", i, pct_exp[i])

scale_fac <- max(abs(site_scores[, c("RDA1","RDA2")])) /
             max(abs(env_scores[, c("RDA1","RDA2")])) * 0.6
env_sc2 <- env_scores
env_sc2$RDA1 <- env_scores$RDA1 * scale_fac
env_sc2$RDA2 <- env_scores$RDA2 * scale_fac

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_segment(data = env_sc2,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "grey30", linewidth = 0.7) +
  geom_point(data = site_scores,
             aes(x = RDA1, y = RDA2, color = site),
             size = 2.5, alpha = 0.8) +
  labs(title = sprintf("RDA forward selection - final model (R2adj = %.3f)", r2_final),
       x = lab_rda(1), y = lab_rda(2), color = "Site") +
  theme(legend.position    = "right",
        legend.key.size    = unit(0.35, "cm"),
        legend.text        = element_text(size = 7))

if (length(site_colors) > 0)
  p_biplot <- p_biplot +
    scale_color_manual(values = site_colors, na.value = "grey50")

if (has_ggrepel) {
  p_biplot <- p_biplot +
    ggrepel::geom_text_repel(data = env_sc2,
                             aes(x = RDA1, y = RDA2, label = variable),
                             size = 3, color = "grey20", fontface = "italic",
                             max.overlaps = 20, seed = 42)
} else {
  p_biplot <- p_biplot +
    geom_text(data = env_sc2,
              aes(x = RDA1, y = RDA2, label = variable),
              size = 3, color = "grey20", fontface = "italic", vjust = -0.5)
}

ggsave(file.path(out_dir, "rda_biplot_rda12.png"), p_biplot,
       width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "rda_biplot_rda12.pdf"), p_biplot,
       width = 10, height = 7)
cat("RDA biplot written\n")

# -------------------------------------------------------------------
# Final summary
# -------------------------------------------------------------------
cat("\n========== 08.1c SUMMARY ==========\n")
cat(sprintf("Candidate variables (08.1a/b): %d\n", length(selected_vars)))
cat(sprintf("Selected by forward selection: %d\n", length(fwd_vars)))
cat(sprintf("Retained after VIF < 10      : %d\n", length(final_vars)))
cat(sprintf("Final R2adj                  : %.4f\n", r2_final))
cat(sprintf("Global R2adj (all %d vars)    : %.4f\n",
            length(selected_vars), r2_global))
cat(sprintf("\nFinal variables for GEA (08.3-08.7):\n"))
for (v in final_vars) cat(sprintf("  %s\n", v))
cat(sprintf("\nVariables list: %s\n", file.path(out_dir, "variables_final.txt")))
cat("=====================================\n")

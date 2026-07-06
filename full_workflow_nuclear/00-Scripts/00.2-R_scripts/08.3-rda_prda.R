#!/usr/bin/env Rscript
# 08.3-rda_prda.R
# Redundancy Analysis (RDA) and partial RDA (pRDA) for genotype-environment
# association (GEA) analysis.
#
# Analysis is conducted at the individual level (137 individuals x 433k SNPs).
# Environmental and geographic predictors are site-level variables assigned to
# each individual based on their sampling site. Within-site genetic variation
# is captured in residuals, which is the expected behaviour: the analysis
# targets loci whose variation co-varies with between-site environmental
# gradients (Forester et al. 2018, Mol Ecol Res).
#
# Pipeline:
#   (1) Genomic PCA — first N_PCS axes used as neutral structure covariates
#   (2) Variance partitioning — quantifies the unique fractions of genetic
#       variance explained by environment, geography, and population structure
#   (3) pRDA genome scan (rdadapt) — identifies candidate SNPs with unusual
#       loadings on constrained RDA axes after conditioning on structure + geo
#
# Inputs:
#   genotype_012_imputed.tsv  — individuals x SNPs, 0/1/2 dosage (imputed)
#   snp_ids.txt               — SNP IDs scaffold:pos:ref:alt
#   sample_metadata.tsv       — sample_id, site, lat, long
#   variables_final.txt       — final env variables from 08.1d (RDA+VIF)
#   *_env_per_site.csv        — site-level environmental data (08.0.x)
#
# Reference: Forester et al. 2018, Mol Ecol Res; Capblancq & Forester 2021

suppressPackageStartupMessages({
  library(vegan)
  library(robust)
  library(qvalue)
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# Genome scan function (Capblancq & Forester 2021)
# Computes Mahalanobis distances on RDA locus scores, returns p- and
# q-values per SNP.
# -------------------------------------------------------------------
rdadapt <- function(rda, K) {
  zscores    <- rda$CCA$v[, 1:as.numeric(K), drop = FALSE]
  resscale   <- apply(zscores, 2, scale)
  resmaha    <- covRob(resscale, distance = TRUE, na.action = na.omit,
                       estim = "pairwiseGK")$dist
  lambda     <- median(resmaha) / qchisq(0.5, df = K)
  pvals      <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
  qv         <- qvalue(pvals)
  data.frame(p.values = pvals, q.values = qv$qvalues)
}

# -------------------------------------------------------------------
# STEP 1: Parse arguments
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 15) {
  stop("Usage: 08.3-rda_prda.R geno_matrix snp_ids sample_meta vars_final chelsa envirem manual terra out_rda out_prda out_vp out_plots N_PCS K_RDA FDR")
}
geno_file  <- args[1]
snp_file   <- args[2]
meta_file  <- args[3]
vars_file  <- args[4]
chelsa_csv <- args[5]
envirem_csv <- args[6]
manual_csv <- args[7]
terra_csv  <- args[8]
out_rda    <- args[9]
out_prda   <- args[10]
out_vp     <- args[11]
out_plots  <- args[12]
N_PCS      <- as.integer(args[13])
K_RDA      <- as.integer(args[14])
fdr_level  <- as.numeric(args[15])

cat(sprintf("N_PCS = %d | K_RDA = %d | FDR = %.2f\n", N_PCS, K_RDA, fdr_level))

# -------------------------------------------------------------------
# STEP 2: Load final environmental variables
# -------------------------------------------------------------------
cat("Reading selected variables:", vars_file, "\n")
var_lines     <- readLines(vars_file)
selected_vars <- trimws(var_lines[!grepl("^#", var_lines) & nchar(trimws(var_lines)) > 0])
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
for (src in list(envirem, manual, terra)) {
  cols <- intersect(selected_vars, names(src))
  if (length(cols) > 0)
    env_site <- merge(env_site, src[, c("site", cols)], by = "site", all.x = TRUE)
}
cat(sprintf("INFO: %d sites, %d env variables\n", nrow(env_site), length(selected_vars)))

# -------------------------------------------------------------------
# STEP 4: Assign env + geography to individuals
# -------------------------------------------------------------------
cat("Loading sample metadata:", meta_file, "\n")
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)

ind_data <- merge(meta[, c("sample_id", "site", "lat", "long")],
                  env_site[, c("site", selected_vars)],
                  by = "site", all.x = FALSE)
ind_data <- ind_data[complete.cases(ind_data), ]
cat(sprintf("INFO: %d individuals with complete env + geo data\n", nrow(ind_data)))

# -------------------------------------------------------------------
# STEP 5: Load genotype matrix
# -------------------------------------------------------------------
cat("Loading genotype matrix:", geno_file, "\n")
cat("  (this may take several minutes)\n")
geno_raw <- read.table(geno_file, header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       row.names = 1)
cat(sprintf("INFO: %d individuals x %d SNPs loaded\n", nrow(geno_raw), ncol(geno_raw)))

common_ind <- intersect(rownames(geno_raw), ind_data$sample_id)
cat(sprintf("INFO: %d individuals in common\n", length(common_ind)))

geno_mat <- as.matrix(geno_raw[common_ind, ])
ind_data  <- ind_data[match(common_ind, ind_data$sample_id), ]
rownames(ind_data) <- ind_data$sample_id

# Hellinger transformation
cat("Applying Hellinger transformation...\n")
geno_hell <- decostand(geno_mat, method = "hellinger")
cat(sprintf("INFO: Hellinger matrix: %d x %d\n", nrow(geno_hell), ncol(geno_hell)))

# -------------------------------------------------------------------
# STEP 6: Load and parse SNP IDs
# -------------------------------------------------------------------
cat("Loading SNP IDs:", snp_file, "\n")
snp_id_raw <- read.table(snp_file, header = TRUE, stringsAsFactors = FALSE)[[1]]
snp_parts  <- strsplit(snp_id_raw, ":", fixed = TRUE)
snp_df <- data.frame(
  snp_id   = snp_id_raw,
  scaffold = sapply(snp_parts, `[`, 1),
  position = as.integer(sapply(snp_parts, `[`, 2)),
  ref      = sapply(snp_parts, `[`, 3),
  alt      = sapply(snp_parts, `[`, 4),
  stringsAsFactors = FALSE
)

scaffold_nums  <- as.integer(gsub(".*_(\\d+)$", "\\1", snp_df$scaffold))
scaffold_order <- unique(snp_df$scaffold[order(scaffold_nums)])
scaffold_sizes <- tapply(snp_df$position, snp_df$scaffold, max)
scaffold_offset <- c(0, cumsum(as.numeric(scaffold_sizes[scaffold_order[-length(scaffold_order)]])))
names(scaffold_offset) <- scaffold_order
snp_df$cum_pos <- snp_df$position + scaffold_offset[snp_df$scaffold]

# -------------------------------------------------------------------
# STEP 7: Build predictor matrices
# -------------------------------------------------------------------
# Environmental predictors (scaled)
env_mat <- scale(as.matrix(ind_data[, selected_vars]))
env_df  <- as.data.frame(env_mat)

# Geographic predictors (scaled)
geo_mat <- scale(as.matrix(ind_data[, c("lat", "long")]))
geo_df  <- as.data.frame(geo_mat)

# -------------------------------------------------------------------
# STEP 8: Genomic PCA — neutral population structure axes
# -------------------------------------------------------------------
cat(sprintf("\nSTEP 8: Genomic PCA for population structure (N_PCS = %d)...\n", N_PCS))
pca_geno <- rda(geno_hell, scale = FALSE)

# Screeplot
png(file.path(out_plots, "genomic_pca_screeplot.png"),
    width = 900, height = 600, res = 150)
screeplot(pca_geno, npcs = 15, type = "barplot",
          main = "Genomic PCA — screeplot (neutral structure)")
dev.off()

# Extract first N_PCS axes
struct_scores <- scores(pca_geno, choices = 1:N_PCS,
                        display = "sites", scaling = 0)
struct_df <- as.data.frame(struct_scores)
colnames(struct_df) <- paste0("PC", 1:N_PCS)
rownames(struct_df) <- common_ind

# PCA biplot PC1 x PC2 coloured by site
pca_plot_df <- data.frame(struct_df,
                           site = ind_data$site,
                           stringsAsFactors = FALSE)

# Load site color palette
site_colors_file <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
if (file.exists(site_colors_file)) {
  sc <- read.csv(site_colors_file, stringsAsFactors = FALSE)
  site_colors <- setNames(sc[[2]], sc[[1]])
} else {
  site_colors <- NULL
}

g_pca <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = site)) +
  geom_point(size = 2.5, alpha = 0.85) +
  {if (!is.null(site_colors))
     scale_color_manual(values = site_colors, na.value = "grey50")
  } +
  labs(title = "Genomic PCA — PC1 x PC2 (neutral structure)",
       x = "PC1", y = "PC2", color = "Site") +
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"))

ggsave(file.path(out_plots, "genomic_pca_PC1_PC2.png"), g_pca,
       width = 10, height = 7, dpi = 150)
ggsave(file.path(out_plots, "genomic_pca_PC1_PC2.pdf"), g_pca,
       width = 10, height = 7)

cat(sprintf("  Genomic PCA done. PC1-PC%d used as structure covariates.\n", N_PCS))

# -------------------------------------------------------------------
# STEP 9: Variance partitioning
# -------------------------------------------------------------------
cat("\nSTEP 9: Variance partitioning (env | geo | structure)...\n")

# varpart requires data frames, not scaled matrices, for formula interface
vp <- varpart(geno_hell, env_df, geo_df, struct_df)

sink(file.path(out_vp, "variance_partitioning_results.txt"))
cat("=== Variance Partitioning: env | geography | structure ===\n\n")
print(vp)
sink()

# Venn diagram
png(file.path(out_plots, "variance_partitioning_venn.png"),
    width = 800, height = 800, res = 150)
plot(vp,
     Xnames = c("Environment", "Geography", "Structure"),
     bg     = c("#4393C3", "#74C476", "#F4A582"),
     alpha  = 80,
     main   = "Variance partitioning — genome-wide")
dev.off()

# Extract fraction table
vp_df <- data.frame(
  fraction    = c("[a] pure env", "[b] pure geo", "[c] pure struct",
                  "[ab] env+geo", "[ac] env+struct", "[bc] geo+struct",
                  "[abc] shared", "residuals"),
  adj_R2      = c(vp$part$indfract$Adj.R.squared,
                  vp$part$fract$Adj.R.squared[4]),
  stringsAsFactors = FALSE
)
write.table(vp_df, file.path(out_vp, "variance_partitioning_fractions.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("  Variance partitioning done.\n")
cat(sprintf("  [a] pure env   : %.4f\n", vp$part$indfract$Adj.R.squared[1]))
cat(sprintf("  [b] pure geo   : %.4f\n", vp$part$indfract$Adj.R.squared[2]))
cat(sprintf("  [c] pure struct: %.4f\n", vp$part$indfract$Adj.R.squared[3]))

# -------------------------------------------------------------------
# STEP 10: pRDA — full model + per-fraction models
# -------------------------------------------------------------------
cat("\nSTEP 10: pRDA models...\n")

all_covars <- cbind(env_df, geo_df, struct_df)

# Full model
pRDA_full <- rda(geno_hell ~ . , data = all_covars)
r2_full   <- RsquareAdj(pRDA_full)
cat(sprintf("  Full model R2adj = %.4f\n", r2_full$adj.r.squared))

# Pure environment | condition(geo + structure)
cond_geo_struct <- cbind(geo_df, struct_df)
pRDA_env <- rda(geno_hell ~ . + Condition(lat + long +
                  PC1 + PC2 + PC3),
                data = cbind(env_df, cond_geo_struct))

# Build formula dynamically from selected_vars and N_PCS
pc_names  <- paste0("PC", 1:N_PCS)
env_form  <- paste(selected_vars, collapse = " + ")
geo_form  <- "lat + long"
struct_form <- paste(pc_names, collapse = " + ")

pRDA_env <- rda(
  as.formula(sprintf("geno_hell ~ %s + Condition(%s + %s)",
                     env_form, geo_form, struct_form)),
  data = cbind(env_df, geo_df, struct_df)
)
r2_env <- RsquareAdj(pRDA_env)
cat(sprintf("  Pure env model R2adj = %.4f\n", r2_env$adj.r.squared))

pRDA_geo <- rda(
  as.formula(sprintf("geno_hell ~ %s + Condition(%s + %s)",
                     geo_form, env_form, struct_form)),
  data = cbind(env_df, geo_df, struct_df)
)
r2_geo <- RsquareAdj(pRDA_geo)
cat(sprintf("  Pure geo model R2adj = %.4f\n", r2_geo$adj.r.squared))

pRDA_struct <- rda(
  as.formula(sprintf("geno_hell ~ %s + Condition(%s + %s)",
                     struct_form, env_form, geo_form)),
  data = cbind(env_df, geo_df, struct_df)
)
r2_struct <- RsquareAdj(pRDA_struct)
cat(sprintf("  Pure struct model R2adj = %.4f\n", r2_struct$adj.r.squared))

# Significance tests (999 permutations)
cat("  Running permutation tests (999 perm)...\n")
anova_env    <- anova(pRDA_env,    permutations = 999)
anova_geo    <- anova(pRDA_geo,    permutations = 999)
anova_struct <- anova(pRDA_struct, permutations = 999)

# Save pRDA summary
sink(file.path(out_prda, "prda_model_summary.txt"))
cat("=== pRDA Model Summary ===\n\n")
cat("--- Full model ---\n"); print(summary(pRDA_full, display = NULL))
cat("\n--- Pure environment ---\n"); print(r2_env); print(anova_env)
cat("\n--- Pure geography ---\n");   print(r2_geo); print(anova_geo)
cat("\n--- Pure structure ---\n");   print(r2_struct); print(anova_struct)
sink()

prda_table <- data.frame(
  model   = c("full", "pure_env", "pure_geo", "pure_struct"),
  R2adj   = c(r2_full$adj.r.squared, r2_env$adj.r.squared,
              r2_geo$adj.r.squared,  r2_struct$adj.r.squared),
  p_value = c(NA,
              anova_env[1, "Pr(>F)"],
              anova_geo[1, "Pr(>F)"],
              anova_struct[1, "Pr(>F)"]),
  stringsAsFactors = FALSE
)
write.table(prda_table, file.path(out_prda, "prda_r2_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# -------------------------------------------------------------------
# STEP 11: RDA genome scan (rdadapt) on pure env pRDA
# -------------------------------------------------------------------
cat(sprintf("\nSTEP 11: RDA genome scan (rdadapt, K = %d)...\n", K_RDA))

# Screeplot of constrained axes
png(file.path(out_plots, "prda_env_screeplot.png"),
    width = 900, height = 600, res = 150)
screeplot(pRDA_env, main = "pRDA (pure env) — constrained axes screeplot")
dev.off()

cat("  Running rdadapt...\n")
rda_results <- rdadapt(pRDA_env, K_RDA)
cat(sprintf("  rdadapt done. %d SNPs tested.\n", nrow(rda_results)))

# --- Candidate identification ---
candidates_idx <- which(rda_results$q.values < fdr_level)
n_cand <- length(candidates_idx)
cat(sprintf("  %d candidate SNPs at FDR < %.2f\n", n_cand, fdr_level))

cand_df <- data.frame(
  snp_id   = snp_df$snp_id[candidates_idx],
  scaffold = snp_df$scaffold[candidates_idx],
  position = snp_df$position[candidates_idx],
  ref      = snp_df$ref[candidates_idx],
  alt      = snp_df$alt[candidates_idx],
  p_value  = rda_results$p.values[candidates_idx],
  q_value  = rda_results$q.values[candidates_idx],
  stringsAsFactors = FALSE
)
cand_df <- cand_df[order(cand_df$q_value), ]
write.table(cand_df, file.path(out_rda, "rdadapt_candidates.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Full results table
full_results <- data.frame(
  snp_id   = snp_df$snp_id,
  scaffold = snp_df$scaffold,
  position = snp_df$position,
  p_value  = rda_results$p.values,
  q_value  = rda_results$q.values,
  stringsAsFactors = FALSE
)
write.table(full_results, file.path(out_rda, "rdadapt_all_snps.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Manhattan plot ---
scaffold_colors <- rep(c("#2166AC", "#4DAC26"), length.out = length(scaffold_order))
names(scaffold_colors) <- scaffold_order

axis_df <- aggregate(cum_pos ~ scaffold, data = snp_df, FUN = mean)
axis_df <- axis_df[match(scaffold_order, axis_df$scaffold), ]

q_thresh <- ifelse(n_cand > 0,
                   -log10(max(rda_results$q.values[candidates_idx])),
                   NA)

manhattan_df <- data.frame(
  cum_pos  = snp_df$cum_pos,
  scaffold = snp_df$scaffold,
  neg_log_q = -log10(rda_results$q.values),
  is_cand  = seq_len(nrow(snp_df)) %in% candidates_idx,
  stringsAsFactors = FALSE
)

g_man <- ggplot(manhattan_df, aes(x = cum_pos, y = neg_log_q, color = scaffold)) +
  geom_point(size = 0.3, alpha = 0.4, show.legend = FALSE) +
  scale_color_manual(values = scaffold_colors) +
  {if (!is.na(q_thresh))
     geom_hline(yintercept = q_thresh, linetype = "dashed",
                color = "#D6604D", linewidth = 0.5)
  } +
  {if (n_cand > 0)
     geom_point(data = subset(manhattan_df, is_cand),
                aes(x = cum_pos, y = neg_log_q),
                color = "#D6604D", size = 1.5, alpha = 0.9)
  } +
  scale_x_continuous(
    breaks = axis_df$cum_pos,
    labels = gsub("Super-Scaffold_", "SS", axis_df$scaffold)
  ) +
  labs(
    title = sprintf("RDA genome scan (rdadapt) — %d candidates, FDR < %.2f",
                    n_cand, fdr_level),
    x = "Scaffold", y = expression(-log[10](q-value))
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(out_plots, "rdadapt_manhattan.png"), g_man,
       width = 14, height = 4, dpi = 150)
ggsave(file.path(out_plots, "rdadapt_manhattan.pdf"), g_man,
       width = 14, height = 4)

# --- RDA biplot (loci + env vectors) ---
locus_scores <- scores(pRDA_env, choices = 1:2,
                       display = "species", scaling = "none")
biplot_loci <- data.frame(
  snp_id = rownames(locus_scores),
  RDA1   = locus_scores[, 1],
  RDA2   = locus_scores[, 2],
  type   = ifelse(rownames(locus_scores) %in% cand_df$snp_id,
                  "Candidate", "Neutral"),
  stringsAsFactors = FALSE
)
biplot_loci$type <- factor(biplot_loci$type, levels = c("Neutral", "Candidate"))
biplot_loci <- biplot_loci[order(biplot_loci$type), ]

biplot_vars <- as.data.frame(scores(pRDA_env, choices = 1:2, display = "bp"))

g_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(data = biplot_loci,
             aes(x = RDA1 * 20, y = RDA2 * 20, color = type),
             size = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("Neutral" = "grey80", "Candidate" = "#D6604D")) +
  geom_segment(data = biplot_vars,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.02, "npc")),
               color = "black", linewidth = 0.4) +
  geom_text(data = biplot_vars,
            aes(x = 1.15 * RDA1, y = 1.15 * RDA2,
                label = rownames(biplot_vars)),
            size = 2.8) +
  labs(title = "pRDA biplot — loci and environmental vectors",
       x = "RDA1", y = "RDA2", color = "Locus") +
  theme(panel.grid = element_blank())

ggsave(file.path(out_plots, "rdadapt_biplot.png"), g_biplot,
       width = 10, height = 8, dpi = 150)
ggsave(file.path(out_plots, "rdadapt_biplot.pdf"), g_biplot,
       width = 10, height = 8)

# -------------------------------------------------------------------
# STEP 12: Summary
# -------------------------------------------------------------------
cat("\n=== RDA/pRDA SUMMARY ===\n")
cat(sprintf("Individuals : %d\n", nrow(geno_hell)))
cat(sprintf("SNPs        : %d\n", ncol(geno_hell)))
cat(sprintf("Env vars    : %s\n", paste(selected_vars, collapse = ", ")))
cat(sprintf("N_PCS struct: %d\n", N_PCS))
cat(sprintf("K_RDA       : %d\n", K_RDA))
cat(sprintf("FDR         : %.2f\n", fdr_level))
cat(sprintf("\nVariance partitioning:\n"))
cat(sprintf("  Pure env      : %.4f\n", r2_env$adj.r.squared))
cat(sprintf("  Pure geo      : %.4f\n", r2_geo$adj.r.squared))
cat(sprintf("  Pure structure: %.4f\n", r2_struct$adj.r.squared))
cat(sprintf("\nCandidates (rdadapt): %d SNPs at FDR < %.2f\n", n_cand, fdr_level))

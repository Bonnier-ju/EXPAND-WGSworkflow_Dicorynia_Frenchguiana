#!/usr/bin/env Rscript
# 08.1b.4-correlation_analysis.R
#
# Merge the 3 candidate variable tables produced by 08.1b.2 (10 custom
# TerraClimate-derived variables), 08.1b.2b (19 BIO variables) and 08.1b.3
# (2 topographic variables: elevation, tri) into a single 31-variable pool,
# compute pairwise Spearman correlations, and flag highly correlated pairs
# for manual review in 08.1b.5.
#
# Spearman confirmed 2026-08-11 (matches the threshold documented in
# CLAUDE.md §4; the original 08.1c run this replaces used Pearson, but that
# inconsistency is now resolved in favour of Spearman going forward).
#
# Input:
#   data/terraclimate_derived_variables_19sites.csv  (site + 10 vars)
#   data/bioclim_variables_19sites.csv                (site + 19 BIO vars)
#   data/topo_variables_19sites.csv                   (site + elevation, tri)
#
# Output:
#   correlation_analysis/env_variables_merged_31vars.csv
#   correlation_analysis/correlation_matrix_31vars.csv
#   correlation_analysis/correlated_pairs_table.tsv     (|r| >= 0.70)
#   correlation_analysis/variable_summary.tsv           (mean/sd/min/max/CV/n_corr)
#   plots/env_correlation_heatmap.png / .pdf
#
# Args: out_dir  [corr_threshold=0.70]

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
out_dir   <- if (length(args) >= 1) args[1] else stop("out_dir required")
corr_thr  <- if (length(args) >= 2) as.numeric(args[2]) else 0.70

data_dir <- file.path(out_dir, "data")
corr_dir <- file.path(out_dir, "correlation_analysis")
plot_dir <- file.path(out_dir, "plots")
dir.create(corr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

derived_csv <- file.path(data_dir, "terraclimate_derived_variables_19sites.csv")
bioclim_csv <- file.path(data_dir, "bioclim_variables_19sites.csv")
topo_csv    <- file.path(data_dir, "topo_variables_19sites.csv")

# ---- Guard clause: verify each input before running ----
for (f in c(derived_csv, bioclim_csv, topo_csv)) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    stop(sprintf("ERROR: required input missing or empty: %s", f))
  }
}

# -------------------------------------------------------------------
# STEP 1: Merge the 3 variable tables on site
# -------------------------------------------------------------------
cat("STEP 1: merging variable tables ...\n")

derived <- read.csv(derived_csv, stringsAsFactors = FALSE)
bioclim <- read.csv(bioclim_csv, stringsAsFactors = FALSE)
topo    <- read.csv(topo_csv, stringsAsFactors = FALSE)

merged <- merge(derived, bioclim, by = "site", all = TRUE)
merged <- merge(merged, topo, by = "site", all = TRUE)
merged <- merged[order(merged$site), ]

n_na_rows <- sum(!complete.cases(merged))
if (n_na_rows > 0) {
  warning(sprintf("%d site(s) have at least one missing variable after merge", n_na_rows))
}

env_vars <- setdiff(colnames(merged), "site")
cat(sprintf("  %d sites x %d candidate variables (10 custom + 19 BIO + 2 topo)\n",
            nrow(merged), length(env_vars)))

out_merged <- file.path(corr_dir, "env_variables_merged_31vars.csv")
write.csv(merged, out_merged, row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_merged))

# -------------------------------------------------------------------
# STEP 2: Spearman correlation matrix
# -------------------------------------------------------------------
cat("STEP 2: computing correlation matrix ...\n")

env_mat <- as.matrix(merged[, env_vars])
storage.mode(env_mat) <- "numeric"
valid_cols <- env_vars[apply(!is.na(env_mat), 2, sum) >= 3]
dropped <- setdiff(env_vars, valid_cols)
if (length(dropped) > 0) {
  warning(sprintf("Dropped from correlation (too few valid sites): %s",
                  paste(dropped, collapse = ", ")))
}

cor_mat <- cor(merged[, valid_cols], method = "spearman", use = "pairwise.complete.obs")

out_cor <- file.path(corr_dir, "correlation_matrix_31vars.csv")
write.csv(cor_mat, out_cor, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_cor))

# -------------------------------------------------------------------
# STEP 3: Flag correlated pairs (|r| >= threshold)
# -------------------------------------------------------------------
cat(sprintf("STEP 3: flagging pairs with |r| >= %.2f ...\n", corr_thr))

pairs_list <- list()
for (i in seq_len(nrow(cor_mat) - 1)) {
  for (j in (i + 1):ncol(cor_mat)) {
    r <- cor_mat[i, j]
    if (!is.na(r) && abs(r) >= corr_thr) {
      pairs_list[[length(pairs_list) + 1]] <- data.frame(
        variable_1  = rownames(cor_mat)[i],
        variable_2  = colnames(cor_mat)[j],
        spearman_r  = round(r, 3),
        stringsAsFactors = FALSE
      )
    }
  }
}
pairs_df <- if (length(pairs_list) > 0) {
  do.call(rbind, pairs_list)[order(-abs(do.call(rbind, pairs_list)$spearman_r)), ]
} else {
  data.frame(variable_1 = character(), variable_2 = character(), spearman_r = numeric())
}
out_pairs <- file.path(corr_dir, "correlated_pairs_table.tsv")
write.table(pairs_df, out_pairs, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  %d correlated pairs (|r|>=%.2f) written to %s\n",
            nrow(pairs_df), corr_thr, out_pairs))

# -------------------------------------------------------------------
# STEP 4: Variable summary (mean/sd/min/max/CV + n correlated partners)
# -------------------------------------------------------------------
cat("STEP 4: writing variable summary ...\n")

n_corr_partners <- sapply(valid_cols, function(v) {
  sum(pairs_df$variable_1 == v | pairs_df$variable_2 == v)
})

summ_rows <- lapply(valid_cols, function(v) {
  vals <- merged[[v]]
  m <- mean(vals, na.rm = TRUE)
  s <- sd(vals, na.rm = TRUE)
  data.frame(
    variable        = v,
    n_sites         = sum(!is.na(vals)),
    mean            = round(m, 4),
    sd              = round(s, 4),
    min             = round(min(vals, na.rm = TRUE), 4),
    max             = round(max(vals, na.rm = TRUE), 4),
    cv_pct          = round(100 * abs(s / m), 2),
    n_correlated    = n_corr_partners[[v]],
    stringsAsFactors = FALSE
  )
})
summ_df <- do.call(rbind, summ_rows)
summ_df <- summ_df[order(-summ_df$n_correlated), ]
out_summ <- file.path(corr_dir, "variable_summary.tsv")
write.table(summ_df, out_summ, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_summ))

low_cv <- summ_df$variable[summ_df$cv_pct < 5]
if (length(low_cv) > 0) {
  cat(sprintf("  NOTE: low spatial variation (CV<5%%): %s\n", paste(low_cv, collapse = ", ")))
}

# -------------------------------------------------------------------
# STEP 5: Correlation heatmap, variables ordered by hierarchical clustering
# -------------------------------------------------------------------
cat("STEP 5: correlation heatmap ...\n")

hc <- hclust(as.dist(1 - abs(cor_mat)), method = "average")
ord <- rownames(cor_mat)[hc$order]

cor_df <- as.data.frame(as.table(cor_mat))
names(cor_df) <- c("var1", "var2", "r")
cor_df$var1 <- factor(cor_df$var1, levels = ord)
cor_df$var2 <- factor(cor_df$var2, levels = ord)
cor_df$above_thr <- !is.na(cor_df$r) & abs(cor_df$r) >= corr_thr

p_cor <- ggplot(cor_df, aes(x = var2, y = var1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", r), fontface = ifelse(above_thr, "bold", "plain")),
            size = 1.6, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman r") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(
    title    = sprintf("Environmental variable correlations (Spearman, n=%d sites, %d variables)",
                       nrow(merged), length(valid_cols)),
    subtitle = sprintf("Bold labels: |r| >= %.2f | ordered by hierarchical clustering", corr_thr),
    x = NULL, y = NULL
  ) +
  theme(
    axis.text        = element_text(size = 6),
    legend.key.height = unit(0.8, "cm"),
    panel.grid       = element_blank()
  )

ggsave(file.path(plot_dir, "env_correlation_heatmap.png"), p_cor,
       width = 20, height = 18, dpi = 300, limitsize = FALSE)
ggsave(file.path(plot_dir, "env_correlation_heatmap.pdf"), p_cor,
       width = 20, height = 18, limitsize = FALSE)
cat("  Correlation heatmap written\n")

cat(sprintf("\nDONE 08.1b.4-correlation_analysis.R\n"))
cat(sprintf("  %d variables | %d correlated pairs (|r|>=%.2f)\n",
            length(valid_cols), nrow(pairs_df), corr_thr))

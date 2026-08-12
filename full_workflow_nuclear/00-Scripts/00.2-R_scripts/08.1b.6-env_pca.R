#!/usr/bin/env Rscript
# 08.1b.6-env_pca.R
#
# PCA on the environmental variables manually selected by Julien in
# correlation_analysis/variable_selection_template.txt (08.1b.5 output,
# edited). Exploratory step - does not select variables further, that
# happens statistically in 08.1b.7 (RDA forward selection).
#
# Site scores are coloured by the K=3 cluster assignment retained for the
# GEA pipeline (confirmed 2026-08-11, cf. CLAUDE.md §4 - distinct from the
# K=4 population-structure narrative in §3), using the cluster colours from
# the original 08.1 plan: Pop_1=#E07B3F, Pop_2=#3B7EC8, Pop_3=#58A45C.
#
# Loading interpretation table lists top contributors per PC mechanically
# (|loading| > 0.25) with an empty column left for Julien's ecological
# interpretation - no label is invented here (consistent with 08.1b.5).
#
# Input:
#   correlation_analysis/env_variables_merged_31vars.csv
#   correlation_analysis/variable_selection_template.txt  (edited by Julien)
#   metadata/sites_by_clusters.csv                          (column "K=3")
#
# Output:
#   pca/env_pca_eigenvalues.tsv
#   pca/env_pca_loadings.tsv
#   pca/env_pca_top_loadings.tsv
#   pca/env_pca_site_scores.tsv
#   plots/env_pca_screeplot.png/.pdf
#   plots/env_pca_biplot_PC1PC2.png/.pdf
#   plots/env_pca_biplot_PC1PC3.png/.pdf
#
# Args: out_dir  sites_by_clusters_csv

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

K3_COLORS <- c(Pop_1 = "#E07B3F", Pop_2 = "#3B7EC8", Pop_3 = "#58A45C")

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
out_dir     <- if (length(args) >= 1) args[1] else stop("out_dir required")
clusters_csv <- if (length(args) >= 2) args[2] else stop("sites_by_clusters_csv required")

corr_dir <- file.path(out_dir, "correlation_analysis")
pca_dir  <- file.path(out_dir, "pca")
plot_dir <- file.path(out_dir, "plots")
dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

merged_csv   <- file.path(corr_dir, "env_variables_merged_31vars.csv")
template_txt <- file.path(corr_dir, "variable_selection_template.txt")

# ---- Guard clause: verify inputs before running ----
for (f in c(merged_csv, template_txt, clusters_csv)) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    stop(sprintf("ERROR: required input missing or empty: %s", f))
  }
}

# -------------------------------------------------------------------
# STEP 1: Parse kept variables from the edited template
# A kept line has had its leading "# " removed, so it now starts directly
# with the variable name followed by " # group=..." - this pattern is
# unambiguous versus category headers / instructions (which never match it,
# commented or not).
# -------------------------------------------------------------------
cat("STEP 1: parsing edited variable_selection_template.txt ...\n")

lines <- readLines(template_txt)
kept_lines <- grep("^[A-Za-z0-9_]+\\s*#\\s*group=", lines, value = TRUE)
kept_vars <- sub("^([A-Za-z0-9_]+).*", "\\1", kept_lines)

if (length(kept_vars) == 0) {
  stop(paste(
    "ERROR: no variable selected in variable_selection_template.txt.",
    "Edit the file (uncomment - remove the leading '# ' - the variables to keep)",
    "before running 08.1b.6. See 08.1b.5 output for instructions.",
    sep = "\n"
  ))
}
if (length(kept_vars) < 2) {
  stop("ERROR: at least 2 variables are required for a PCA.")
}
cat(sprintf("  %d variables selected: %s\n", length(kept_vars), paste(kept_vars, collapse = ", ")))
if (length(kept_vars) < 5 || length(kept_vars) > 12) {
  cat(sprintf("  NOTE: plan target was ~8-10 variables - got %d. Not blocking, just flagging.\n",
              length(kept_vars)))
}

# -------------------------------------------------------------------
# STEP 2: Load data and cluster assignment
# -------------------------------------------------------------------
cat("STEP 2: loading data ...\n")

merged <- read.csv(merged_csv, stringsAsFactors = FALSE)
missing_vars <- setdiff(kept_vars, colnames(merged))
if (length(missing_vars) > 0) {
  stop(sprintf("Selected variable(s) not found in merged data: %s",
              paste(missing_vars, collapse = ", ")))
}

clusters <- read.csv(clusters_csv, stringsAsFactors = FALSE)
stopifnot("Sites" %in% colnames(clusters), "K.3" %in% colnames(clusters) || "K=3" %in% colnames(clusters))
k3_col <- if ("K=3" %in% colnames(clusters)) "K=3" else "K.3"
clusters <- clusters[, c("Sites", k3_col)]
colnames(clusters) <- c("site", "cluster")

pca_data <- merged[, c("site", kept_vars)]
pca_data <- merge(pca_data, clusters, by = "site", all.x = TRUE)
n_missing_cluster <- sum(is.na(pca_data$cluster))
if (n_missing_cluster > 0) {
  warning(sprintf("%d site(s) have no K=3 cluster assignment", n_missing_cluster))
}

n_na <- sum(!complete.cases(pca_data[, kept_vars]))
if (n_na > 0) {
  warning(sprintf("%d site(s) dropped from PCA due to missing values in selected variables", n_na))
  pca_data <- pca_data[complete.cases(pca_data[, kept_vars]), ]
}
cat(sprintf("  %d sites retained for PCA\n", nrow(pca_data)))

# -------------------------------------------------------------------
# STEP 3: PCA
# -------------------------------------------------------------------
cat("STEP 3: running PCA (centered, scaled) ...\n")

pca <- prcomp(pca_data[, kept_vars], center = TRUE, scale. = TRUE)
n_pc_report <- min(4, ncol(pca$x))

var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
eig_df <- data.frame(
  PC              = paste0("PC", seq_along(pca$sdev)),
  eigenvalue      = round(pca$sdev^2, 4),
  variance_pct    = round(100 * var_explained, 2),
  cumulative_pct  = round(100 * cumsum(var_explained), 2)
)
write.table(eig_df, file.path(pca_dir, "env_pca_eigenvalues.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  PC1-PC%d explain %.1f%% of variance\n", n_pc_report,
            eig_df$cumulative_pct[n_pc_report]))

loadings <- as.data.frame(pca$rotation[, seq_len(n_pc_report), drop = FALSE])
loadings$variable <- rownames(loadings)
loadings <- loadings[, c("variable", paste0("PC", seq_len(n_pc_report)))]
write.table(loadings, file.path(pca_dir, "env_pca_loadings.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

scores <- as.data.frame(pca$x[, seq_len(n_pc_report), drop = FALSE])
scores$site    <- pca_data$site
scores$cluster <- pca_data$cluster
scores <- scores[, c("site", "cluster", paste0("PC", seq_len(n_pc_report)))]
write.table(scores, file.path(pca_dir, "env_pca_site_scores.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Top contributing variables per PC (|loading| > 0.25), mechanical only.
top_rows <- list()
for (pc in paste0("PC", seq_len(n_pc_report))) {
  strong <- loadings[abs(loadings[[pc]]) > 0.25, c("variable", pc)]
  if (nrow(strong) == 0) next
  strong <- strong[order(-abs(strong[[pc]])), ]
  top_rows[[pc]] <- data.frame(
    PC = pc, variable = strong$variable, loading = round(strong[[pc]], 3),
    sign = ifelse(strong[[pc]] > 0, "+", "-"),
    ecological_interpretation = "",  # left blank for Julien - not invented
    stringsAsFactors = FALSE
  )
}
top_df <- do.call(rbind, top_rows)
write.table(top_df, file.path(pca_dir, "env_pca_top_loadings.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote eigenvalues, loadings, site scores, top loadings to %s\n", pca_dir))

# -------------------------------------------------------------------
# STEP 4: Scree plot
# -------------------------------------------------------------------
cat("STEP 4: plots ...\n")

p_scree <- ggplot(eig_df, aes(x = factor(PC, levels = PC), y = variance_pct, group = 1)) +
  geom_col(fill = "grey70") +
  geom_line(aes(y = cumulative_pct), color = "#B2182B", linewidth = 1) +
  geom_point(aes(y = cumulative_pct), color = "#B2182B", size = 2) +
  labs(title = "Environmental PCA - variance explained",
       subtitle = sprintf("%d variables, %d sites", length(kept_vars), nrow(pca_data)),
       x = NULL, y = "% variance (bars) / cumulative % (line)")
ggsave(file.path(plot_dir, "env_pca_screeplot.png"), p_scree, width = 8, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "env_pca_screeplot.pdf"), p_scree, width = 8, height = 5)

# -------------------------------------------------------------------
# STEP 5: Biplots PC1xPC2 and PC1xPC3
# -------------------------------------------------------------------
biplot_pair <- function(pc_x, pc_y) {
  sc <- scores[, c("site", "cluster", pc_x, pc_y)]
  colnames(sc) <- c("site", "cluster", "x", "y")

  ld <- loadings[, c("variable", pc_x, pc_y)]
  colnames(ld) <- c("variable", "x", "y")
  mult <- 0.8 * max(abs(sc[, c("x", "y")])) / max(abs(ld[, c("x", "y")]))
  ld$x <- ld$x * mult
  ld$y <- ld$y * mult

  xlab <- sprintf("%s (%.1f%%)", pc_x, eig_df$variance_pct[eig_df$PC == pc_x])
  ylab <- sprintf("%s (%.1f%%)", pc_y, eig_df$variance_pct[eig_df$PC == pc_y])

  ggplot() +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = x, yend = y),
                 arrow = arrow(length = unit(0.2, "cm")), color = "grey40") +
    geom_text(data = ld, aes(x = x, y = y, label = variable),
              size = 2.6, color = "grey20", fontface = "italic") +
    geom_point(data = sc, aes(x = x, y = y, color = cluster), size = 2.5) +
    geom_text(data = sc, aes(x = x, y = y, label = site, color = cluster),
              size = 2.4, vjust = -0.7, show.legend = FALSE) +
    scale_color_manual(values = K3_COLORS, na.value = "grey50", name = "Cluster (K=3)") +
    coord_equal() +
    labs(title = sprintf("Environmental PCA biplot - %s x %s", pc_x, pc_y), x = xlab, y = ylab)
}

p12 <- biplot_pair("PC1", "PC2")
ggsave(file.path(plot_dir, "env_pca_biplot_PC1PC2.png"), p12, width = 9, height = 7, dpi = 300)
ggsave(file.path(plot_dir, "env_pca_biplot_PC1PC2.pdf"), p12, width = 9, height = 7)

if (n_pc_report >= 3) {
  p13 <- biplot_pair("PC1", "PC3")
  ggsave(file.path(plot_dir, "env_pca_biplot_PC1PC3.png"), p13, width = 9, height = 7, dpi = 300)
  ggsave(file.path(plot_dir, "env_pca_biplot_PC1PC3.pdf"), p13, width = 9, height = 7)
}
cat("  Scree plot and biplot(s) written\n")

cat(sprintf("\nDONE 08.1b.6-env_pca.R\n"))
cat(sprintf("  %d variables | %d sites | PC1 explains %.1f%%\n",
            length(kept_vars), nrow(pca_data), eig_df$variance_pct[1]))

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

base_dir_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_test/07.1-pca_ibd_fst"
base_dir <- if (length(args) >= 1) args[1] else base_dir_default

eigenvec_file <- file.path(base_dir, "pca.eigenvec")
eigenval_file <- file.path(base_dir, "pca.eigenval")
pop_map_file <- file.path(base_dir, "sample_population_map.tsv")
out_dir <- file.path(base_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(eigenvec_file)) stop("Missing file: ", eigenvec_file)
if (!file.exists(eigenval_file)) stop("Missing file: ", eigenval_file)
if (!file.exists(pop_map_file)) stop("Missing file: ", pop_map_file)

pca <- read.table(eigenvec_file, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
eig <- scan(eigenval_file, what = numeric(), quiet = TRUE)
pop <- read.table(pop_map_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(pop) <- c("sample_id", "population")

pca$sample_id <- pca$IID
dat <- merge(pca, pop, by = "sample_id", all.x = TRUE, sort = FALSE)
dat$population[is.na(dat$population) | dat$population == ""] <- "UNKNOWN"

# Site/population centroids (one point per site)
dat_centroid <- aggregate(
  x = dat[, grep("^PC[0-9]+$", colnames(dat)), drop = FALSE],
  by = list(population = dat$population),
  FUN = mean
)
dat_centroid$n_individuals <- as.integer(table(dat$population)[dat_centroid$population])
dat_centroid$sample_id <- dat_centroid$population

pc_var <- eig / sum(eig) * 100
if (length(pc_var) < 4) {
  stop("Need at least 4 PCs in eigenvalues to plot requested combinations.")
}

make_pca_plot <- function(df, pcx, pcy, pc_var_vec, title_prefix = "PCA", label_col = "sample_id") {
  x_col <- paste0("PC", pcx)
  y_col <- paste0("PC", pcy)
  ggplot(df, aes_string(x = x_col, y = y_col, color = "population", label = label_col)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text(vjust = -0.6, size = 2.8, show.legend = FALSE) +
    labs(
      title = sprintf("%s - PC%d vs PC%d", title_prefix, pcx, pcy),
      x = sprintf("PC%d (%.2f%%)", pcx, pc_var_vec[pcx]),
      y = sprintf("PC%d (%.2f%%)", pcy, pc_var_vec[pcy]),
      color = "Population"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

pairs_to_plot <- list(
  c(1, 2),
  c(1, 3),
  c(2, 3),
  c(1, 4),
  c(2, 4),
  c(3, 4)
)

for (pp in pairs_to_plot) {
  pcx <- pp[1]
  pcy <- pp[2]

  # Version 1: all individuals
  plot_ind <- make_pca_plot(
    dat, pcx, pcy, pc_var,
    title_prefix = "PCA (all individuals)",
    label_col = "sample_id"
  )
  prefix_ind <- sprintf("pca_individuals_pc%d_pc%d", pcx, pcy)
  ggsave(file.path(out_dir, paste0(prefix_ind, ".png")), plot_ind, width = 8, height = 6, dpi = 300)
  ggsave(file.path(out_dir, paste0(prefix_ind, ".pdf")), plot_ind, width = 8, height = 6)

  # Version 2: one centroid per site/population
  plot_cent <- make_pca_plot(
    dat_centroid, pcx, pcy, pc_var,
    title_prefix = "PCA (site centroids)",
    label_col = "population"
  )
  prefix_cent <- sprintf("pca_site_centroids_pc%d_pc%d", pcx, pcy)
  ggsave(file.path(out_dir, paste0(prefix_cent, ".png")), plot_cent, width = 8, height = 6, dpi = 300)
  ggsave(file.path(out_dir, paste0(prefix_cent, ".pdf")), plot_cent, width = 8, height = 6)
}

summary_file <- file.path(out_dir, "pca_variance_explained.tsv")
write.table(
  data.frame(
    PC = paste0("PC", seq_along(pc_var)),
    variance_percent = round(pc_var, 6)
  ),
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("PCA plots written to:", out_dir, "\n")

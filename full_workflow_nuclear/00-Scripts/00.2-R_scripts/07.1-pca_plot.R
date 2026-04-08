#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# Site-specific color palette (consistent across all population graphics)
SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else {
  warning("Site color file not found: ", SITE_COLORS_FILE, " — using ggplot2 defaults")
  c()
}

args <- commandArgs(trailingOnly = TRUE)

base_dir_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_full/07.1-pca_ibd_fst"
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

# Outgroup samples excluded from the second group of plots.
# These are the 5 individuals from Cameroun/Benin used as outgroups.
OUTGROUP_SAMPLES <- c("DBR_AGN", "DBR_AGO", "DBR_AGP", "DBR_AGQ", "DBR_AGR")
dat_no_outgroup <- dat[!dat$sample_id %in% OUTGROUP_SAMPLES, ]

make_centroids <- function(d) {
  pc_cols <- grep("^PC[0-9]+$", colnames(d))
  ctr <- aggregate(
    x = d[, pc_cols, drop = FALSE],
    by = list(population = d$population),
    FUN = mean
  )
  ctr$n_individuals <- as.integer(table(d$population)[ctr$population])
  ctr$sample_id <- ctr$population
  ctr
}

dat_centroid            <- make_centroids(dat)
dat_centroid_no_outgroup <- make_centroids(dat_no_outgroup)

pc_var <- eig / sum(eig) * 100
if (length(pc_var) < 4) {
  stop("Need at least 4 PCs in eigenvalues to plot requested combinations.")
}

make_pca_plot <- function(df, pcx, pcy, pc_var_vec, title_prefix = "PCA", label_col = "sample_id",
                          color_values = site_colors) {
  x_col <- paste0("PC", pcx)
  y_col <- paste0("PC", pcy)
  p <- ggplot(df, aes_string(x = x_col, y = y_col, color = "population", label = label_col)) +
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
  if (length(color_values) > 0) {
    p <- p + scale_color_manual(values = color_values, na.value = "grey50")
  }
  p
}

pairs_to_plot <- list(c(1,2), c(1,3), c(2,3), c(1,4), c(2,4), c(3,4))

# Two groups of plots: all samples, then without outgroups.
# PC coordinates are those from the full PCA (not recomputed on the subset);
# the variance percentages therefore reflect the full 164-sample decomposition.
plot_groups <- list(
  list(
    subdir       = "all_samples",
    dat_ind      = dat,
    dat_cent     = dat_centroid,
    title_suffix = "all 164 ind."
  ),
  list(
    subdir       = "no_outgroup",
    dat_ind      = dat_no_outgroup,
    dat_cent     = dat_centroid_no_outgroup,
    title_suffix = sprintf("no outgroup, n=%d", nrow(dat_no_outgroup))
  )
)

for (grp in plot_groups) {
  sub_dir <- file.path(out_dir, grp$subdir)
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

  for (pp in pairs_to_plot) {
    pcx <- pp[1]
    pcy <- pp[2]

    # Individual-level plot
    plot_ind <- make_pca_plot(
      grp$dat_ind, pcx, pcy, pc_var,
      title_prefix = sprintf("PCA (%s)", grp$title_suffix),
      label_col    = "sample_id"
    )
    fn_ind <- sprintf("pca_individuals_pc%d_pc%d", pcx, pcy)
    ggsave(file.path(sub_dir, paste0(fn_ind, ".png")), plot_ind, width = 8, height = 6, dpi = 300)
    ggsave(file.path(sub_dir, paste0(fn_ind, ".pdf")), plot_ind, width = 8, height = 6)

    # Site centroid plot
    plot_cent <- make_pca_plot(
      grp$dat_cent, pcx, pcy, pc_var,
      title_prefix = sprintf("PCA centroids (%s)", grp$title_suffix),
      label_col    = "population"
    )
    fn_cent <- sprintf("pca_site_centroids_pc%d_pc%d", pcx, pcy)
    ggsave(file.path(sub_dir, paste0(fn_cent, ".png")), plot_cent, width = 8, height = 6, dpi = 300)
    ggsave(file.path(sub_dir, paste0(fn_cent, ".pdf")), plot_cent, width = 8, height = 6)
  }

  cat(sprintf("Plots written to: %s\n", sub_dir))
}

# -----------------------------------------------------------------------
# FST: pairwise heatmap and per-site mean FST barplot
# -----------------------------------------------------------------------
# Reads the pairwise FST summary produced by the SLURM step (plink --fst
# for each population pair). Produces:
#   - A heatmap of all pairwise weighted Weir-Cockerham FST values,
#     populations ordered by hierarchical clustering on the FST distance.
#   - A barplot of mean pairwise FST per site (average over all pairs
#     involving that site), colored with the site palette.
# Both are produced in two versions: all sites, and without outgroups.
# -----------------------------------------------------------------------
fst_pair_file <- file.path(base_dir, "fst_pairwise", "fst_pairwise_summary.tsv")
fst_dir <- file.path(out_dir, "fst")
dir.create(fst_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(fst_pair_file)) {
  fst_raw <- read.table(fst_pair_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, na.strings = c("NA", ""))
  fst_raw$fst_weighted <- suppressWarnings(as.numeric(fst_raw$fst_weighted))
  # Clamp small negative FST values (sampling artefact) to 0
  fst_raw$fst_weighted[!is.na(fst_raw$fst_weighted) & fst_raw$fst_weighted < 0] <- 0

  build_fst_matrix <- function(df) {
    pops <- sort(unique(c(df$pop1, df$pop2)))
    mat  <- matrix(0, nrow = length(pops), ncol = length(pops),
                   dimnames = list(pops, pops))
    for (i in seq_len(nrow(df))) {
      p1 <- df$pop1[i]; p2 <- df$pop2[i]; v <- df$fst_weighted[i]
      if (!is.na(v)) { mat[p1, p2] <- v; mat[p2, p1] <- v }
    }
    mat
  }

  make_fst_plots <- function(df, label, outgroup_samples = NULL) {
    if (!is.null(outgroup_samples)) {
      df <- df[!df$pop1 %in% outgroup_samples & !df$pop2 %in% outgroup_samples, ]
    }
    if (nrow(df) == 0) return(invisible(NULL))
    mat <- build_fst_matrix(df)

    # Order populations by hierarchical clustering on FST distance
    hc       <- hclust(as.dist(mat), method = "average")
    pop_ord  <- hc$labels[hc$order]

    # Long-format symmetric table (including diagonal = 0)
    all_pops <- rownames(mat)
    fst_long <- do.call(rbind, lapply(all_pops, function(p1) {
      data.frame(pop1 = p1, pop2 = all_pops, fst_weighted = mat[p1, all_pops],
                 stringsAsFactors = FALSE)
    }))
    fst_long$pop1 <- factor(fst_long$pop1, levels = pop_ord)
    fst_long$pop2 <- factor(fst_long$pop2, levels = pop_ord)
    fst_max <- max(fst_long$fst_weighted, na.rm = TRUE)

    # Heatmap
    p_heat <- ggplot(fst_long, aes(x = pop2, y = pop1, fill = fst_weighted)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(
        aes(label = ifelse(pop1 != pop2, sprintf("%.3f", fst_weighted), "")),
        size = 2.0, color = "black"
      ) +
      scale_fill_gradientn(
        colours  = c("#f7fbff", "#6baed6", "#08306b"),
        name     = "FST\n(weighted)",
        na.value = "grey90",
        limits   = c(0, fst_max)
      ) +
      labs(
        title = sprintf("Pairwise Weir-Cockerham FST — %s", label),
        x = NULL, y = NULL
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title  = element_text(face = "bold")
      )
    n_pops  <- length(all_pops)
    hw      <- max(6, n_pops * 0.5)
    ggsave(file.path(fst_dir, sprintf("fst_heatmap_%s.png", label)),
           p_heat, width = hw, height = hw * 0.9, dpi = 300)
    ggsave(file.path(fst_dir, sprintf("fst_heatmap_%s.pdf", label)),
           p_heat, width = hw, height = hw * 0.9)

    # Mean pairwise FST per site (barplot, sorted descending)
    off_diag    <- fst_long[fst_long$pop1 != fst_long$pop2, ]
    mean_fst    <- aggregate(fst_weighted ~ pop1, data = off_diag, FUN = mean, na.rm = TRUE)
    colnames(mean_fst) <- c("population", "mean_fst")
    mean_fst    <- mean_fst[order(-mean_fst$mean_fst), ]
    mean_fst$population <- factor(mean_fst$population, levels = mean_fst$population)

    p_bar <- ggplot(mean_fst, aes(x = population, y = mean_fst, fill = population)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = sprintf("%.3f", mean_fst)), vjust = -0.4, size = 2.8) +
      labs(
        title = sprintf("Mean pairwise FST per site — %s", label),
        x = NULL, y = "Mean weighted FST"
      ) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        plot.title   = element_text(face = "bold"),
        legend.position = "none"
      )
    if (length(site_colors) > 0) {
      p_bar <- p_bar + scale_fill_manual(values = site_colors, na.value = "grey50")
    }
    ggsave(file.path(fst_dir, sprintf("fst_mean_per_site_%s.png", label)),
           p_bar, width = max(8, n_pops * 0.5), height = 5, dpi = 300)
    ggsave(file.path(fst_dir, sprintf("fst_mean_per_site_%s.pdf", label)),
           p_bar, width = max(8, n_pops * 0.5), height = 5)

    # Save FST matrix as TSV
    write.table(as.data.frame(mat),
                file = file.path(fst_dir, sprintf("fst_matrix_%s.tsv", label)),
                sep = "\t", quote = FALSE, row.names = TRUE)

    cat(sprintf("FST plots written (%s): %s\n", label, fst_dir))
  }

  make_fst_plots(fst_raw, label = "all_sites")
  make_fst_plots(fst_raw, label = "no_outgroup", outgroup_samples = OUTGROUP_SAMPLES)

} else {
  cat("INFO: pairwise FST file not found, FST plots skipped:", fst_pair_file, "\n")
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

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    library(ggrepel)
    USE_GGREPEL <- TRUE
  } else {
    USE_GGREPEL <- FALSE
  }
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

# Arguments:
#   1: base_dir      — output directory of 07.1
#   2: pop_map_file  — sample_population_map.tsv
#   3: eigenvec_file — pca_<group>.eigenvec
#   4: eigenval_file — pca_<group>.eigenval
#   5: group_label   — e.g. "all_samples", "no_outgroup", "guiana_only"
base_dir_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_full/07.1-pca_ibd_fst"
base_dir      <- if (length(args) >= 1) args[1] else base_dir_default
pop_map_file  <- if (length(args) >= 2) args[2] else file.path(base_dir, "sample_population_map.tsv")
eigenvec_file <- if (length(args) >= 3) args[3] else file.path(base_dir, "pca.eigenvec")
eigenval_file <- if (length(args) >= 4) args[4] else file.path(base_dir, "pca.eigenval")
group_label   <- if (length(args) >= 5) args[5] else "all_samples"

out_dir <- file.path(base_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(eigenvec_file)) stop("Missing file: ", eigenvec_file)
if (!file.exists(eigenval_file)) stop("Missing file: ", eigenval_file)
if (!file.exists(pop_map_file)) stop("Missing file: ", pop_map_file)

pca <- read.table(eigenvec_file, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
eig <- scan(eigenval_file, what = numeric(), quiet = TRUE)
pop <- read.table(pop_map_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(pop) <- c("sample_id", "population")

# Each call to this script corresponds to one independent PCA group.
# dat contains only the individuals present in the eigenvec file
# (i.e. those included in this specific PCA computation).
pca$sample_id <- pca$IID
dat <- merge(pca, pop, by = "sample_id", all.x = TRUE, sort = FALSE)
dat$population[is.na(dat$population) | dat$population == ""] <- "UNKNOWN"

cat(sprintf("INFO: group=%s | %d individuals loaded from eigenvec\n",
            group_label, nrow(dat)))

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

dat_centroid <- make_centroids(dat)

pc_var <- eig / sum(eig) * 100
n_pc   <- length(pc_var)
if (n_pc < 2) stop("Need at least 2 PCs in eigenvalues.")
# Limit pairs to available PCs (max PC4 if available, else max available)
max_pc <- min(4, n_pc)

# Compute nudge vectors for centroid labels based on angle from centroid of all points.
# Each label is pushed outward from the center of the point cloud to reduce overlap.
compute_nudge <- function(x, y, dist = 0.015) {
  cx <- mean(x, na.rm = TRUE)
  cy <- mean(y, na.rm = TRUE)
  rx <- diff(range(x, na.rm = TRUE))
  ry <- diff(range(y, na.rm = TRUE))
  angle <- atan2((y - cy) / (ry + 1e-9), (x - cx) / (rx + 1e-9))
  list(
    nudge_x = cos(angle) * dist * rx,
    nudge_y = sin(angle) * dist * ry
  )
}

# Identify isolated individuals: those whose distance to their nearest neighbour
# in the (PCx, PCy) space exceeds `isolation_quantile` of the distribution of
# all nearest-neighbour distances. Only isolated individuals receive a label.
find_isolated <- function(df, x_col, y_col, isolation_quantile = 0.80) {
  x <- df[[x_col]]
  y <- df[[y_col]]
  n <- nrow(df)
  if (n <= 2) return(rep(TRUE, n))
  # Scale axes to unit variance before computing distances
  xs <- (x - mean(x, na.rm = TRUE)) / (sd(x, na.rm = TRUE) + 1e-9)
  ys <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE) + 1e-9)
  nn_dist <- vapply(seq_len(n), function(i) {
    min(sqrt((xs[-i] - xs[i])^2 + (ys[-i] - ys[i])^2))
  }, numeric(1))
  threshold <- quantile(nn_dist, probs = isolation_quantile, na.rm = TRUE)
  nn_dist >= threshold
}

make_pca_plot <- function(df, pcx, pcy, pc_var_vec, title_prefix = "PCA", label_col = "sample_id",
                          color_values = site_colors, is_centroid = FALSE) {
  x_col <- paste0("PC", pcx)
  y_col <- paste0("PC", pcy)

  p <- ggplot(df, aes_string(x = x_col, y = y_col, color = "population")) +
    geom_point(size = 3, alpha = 0.9) +
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

  if (is_centroid) {
    # Centroid plot: label all points but nudge outward to reduce overlap
    nudge <- compute_nudge(df[[x_col]], df[[y_col]], dist = 0.06)
    if (USE_GGREPEL) {
      p <- p + ggrepel::geom_text_repel(
        aes_string(label = label_col),
        size = 3, show.legend = FALSE,
        max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3
      )
    } else {
      p <- p + geom_text(
        aes_string(label = label_col),
        nudge_x = nudge$nudge_x, nudge_y = nudge$nudge_y,
        size = 3, show.legend = FALSE
      )
    }
  } else {
    # Individual plot: label only isolated individuals
    isolated <- find_isolated(df, x_col, y_col, isolation_quantile = 0.80)
    df_label <- df[isolated, ]
    if (nrow(df_label) > 0) {
      if (USE_GGREPEL) {
        p <- p + ggrepel::geom_text_repel(
          data = df_label,
          aes_string(label = label_col),
          size = 2.5, show.legend = FALSE,
          max.overlaps = Inf, box.padding = 0.3
        )
      } else {
        p <- p + geom_text(
          data = df_label,
          aes_string(label = label_col),
          vjust = -0.6, size = 2.5, show.legend = FALSE
        )
      }
    }
  }

  if (length(color_values) > 0) {
    p <- p + scale_color_manual(values = color_values, na.value = "grey50")
  }
  p
}

# Build list of PC pairs to plot (only pairs within available PCs)
pairs_to_plot <- Filter(function(p) p[1] <= max_pc && p[2] <= max_pc,
                        list(c(1,2), c(1,3), c(2,3), c(1,4), c(2,4), c(3,4)))

# One group per script call — write plots into plots/<group_label>/
sub_dir <- file.path(out_dir, group_label)
dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

title_suffix <- sprintf("%s (n=%d)", group_label, nrow(dat))

for (pp in pairs_to_plot) {
  {
    pcx <- pp[1]
    pcy <- pp[2]

    # Individual-level plot (labels only for isolated individuals)
    plot_ind <- make_pca_plot(
      dat, pcx, pcy, pc_var,
      title_prefix = sprintf("PCA (%s)", title_suffix),
      label_col    = "sample_id",
      is_centroid  = FALSE
    )
    fn_ind <- sprintf("pca_individuals_pc%d_pc%d", pcx, pcy)
    ggsave(file.path(sub_dir, paste0(fn_ind, ".png")), plot_ind, width = 8, height = 6, dpi = 300)
    ggsave(file.path(sub_dir, paste0(fn_ind, ".pdf")), plot_ind, width = 8, height = 6)

    # Site centroid plot (all labels, nudged to avoid overlap)
    plot_cent <- make_pca_plot(
      dat_centroid, pcx, pcy, pc_var,
      title_prefix = sprintf("PCA centroids (%s)", title_suffix),
      label_col    = "population",
      is_centroid  = TRUE
    )
    fn_cent <- sprintf("pca_site_centroids_pc%d_pc%d", pcx, pcy)
    ggsave(file.path(sub_dir, paste0(fn_cent, ".png")), plot_cent, width = 8, height = 6, dpi = 300)
    ggsave(file.path(sub_dir, paste0(fn_cent, ".pdf")), plot_cent, width = 8, height = 6)
  }
}

cat(sprintf("Plots written to: %s\n", sub_dir))


summary_file <- file.path(sub_dir, "pca_variance_explained.tsv")
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

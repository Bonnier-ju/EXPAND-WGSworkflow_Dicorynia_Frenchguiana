#!/usr/bin/env Rscript
# 07b.1-pca_plot.R
# PCA plots for purity-filtered individuals (07b.1).
# Produces two coloring schemes per PC pair:
#   - by cluster (Pop_1/2/3, from original K=3 ADMIXTURE Q)
#   - by site    (geographic origin, west→east)
#
# Inputs (6 positional args):
#   1  out_dir       output directory for this threshold (T70|T80|T90|T95)
#   2  cluster_map   TSV: sample_id <tab> cluster (Pop_1/2/3)
#   3  site_map      TSV: sample_id <tab> site
#   4  eigenvec      PLINK PCA eigenvec file (tab-separated, header)
#   5  eigenval      PLINK PCA eigenval file
#   6  threshold_tag label for plot titles (e.g. "T80")

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

POP_COLORS <- c(Pop_1 = "#EE7600", Pop_2 = "#458B00", Pop_3 = "#CD2626")

SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  site_colors <- setNames(pal$couleur_hex, pal$site)
} else {
  warning("Site color file not found — using ggplot2 defaults for site coloring")
  site_colors <- c()
}

# -------------------------------------------------------------------
# STEP 1: Parse arguments
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6)
  stop("Usage: 07b.1-pca_plot.R out_dir cluster_map site_map eigenvec eigenval threshold_tag")

out_dir       <- args[1]
cluster_file  <- args[2]
site_file     <- args[3]
eigenvec_file <- args[4]
eigenval_file <- args[5]
threshold_tag <- args[6]

plots_dir <- file.path(out_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

for (f in c(cluster_file, site_file, eigenvec_file, eigenval_file))
  if (!file.exists(f)) stop("Missing input file: ", f)

# -------------------------------------------------------------------
# STEP 2: Load data
# -------------------------------------------------------------------
pca  <- read.table(eigenvec_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
eig  <- scan(eigenval_file, what = numeric(), quiet = TRUE)
pc_var <- eig / sum(eig) * 100
n_pc   <- length(pc_var)
max_pc <- min(4, n_pc)

cluster_map <- read.table(cluster_file, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, col.names = c("sample_id", "cluster"))
site_map    <- read.table(site_file, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, col.names = c("sample_id", "site"))

pca$sample_id <- pca$IID
dat <- merge(pca, cluster_map, by = "sample_id", all.x = TRUE)
dat <- merge(dat,  site_map,   by = "sample_id", all.x = TRUE)

dat$cluster[is.na(dat$cluster)] <- "UNKNOWN"
dat$site[is.na(dat$site)]       <- "UNKNOWN"

n_ind <- nrow(dat)
cat(sprintf("INFO: threshold=%s | %d individuals | %d PCs\n",
            threshold_tag, n_ind, n_pc))

variance_df <- data.frame(PC = paste0("PC", seq_along(pc_var)),
                          variance_pct = round(pc_var, 4))
write.table(variance_df, file.path(plots_dir, "pca_variance_explained.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("PC1:", round(pc_var[1], 2), "% | PC2:", round(pc_var[2], 2),
    "% | PC3:", round(pc_var[3], 2), "%\n")

# -------------------------------------------------------------------
# STEP 3: Helper — identify isolated individuals for labeling
# -------------------------------------------------------------------
find_isolated <- function(df, x_col, y_col, q = 0.80) {
  x <- df[[x_col]]; y <- df[[y_col]]
  n <- nrow(df)
  if (n <= 2) return(rep(TRUE, n))
  xs <- (x - mean(x, na.rm = TRUE)) / (sd(x, na.rm = TRUE) + 1e-9)
  ys <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE) + 1e-9)
  nn <- vapply(seq_len(n), function(i)
    min(sqrt((xs[-i] - xs[i])^2 + (ys[-i] - ys[i])^2)), numeric(1))
  nn >= quantile(nn, probs = q, na.rm = TRUE)
}

# -------------------------------------------------------------------
# STEP 4: Helper — build one PCA scatter plot
# -------------------------------------------------------------------
make_plot <- function(df, pcx, pcy, color_col, color_vals, title, label_col = "sample_id") {
  x_col <- paste0("PC", pcx)
  y_col <- paste0("PC", pcy)

  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                       color = .data[[color_col]])) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(
      title = title,
      x = sprintf("PC%d (%.2f%%)", pcx, pc_var[pcx]),
      y = sprintf("PC%d (%.2f%%)", pcy, pc_var[pcy])
    ) +
    theme(
      panel.grid.minor  = element_blank(),
      legend.title      = element_blank(),
      legend.text       = element_text(size = 8)
    )

  # Label only isolated individuals
  isolated <- find_isolated(df, x_col, y_col)
  df_label <- df[isolated, ]

  if (nrow(df_label) > 0) {
    if (USE_GGREPEL) {
      p <- p + ggrepel::geom_text_repel(
        data = df_label,
        aes(label = .data[[label_col]]),
        size = 2.2, show.legend = FALSE,
        max.overlaps = Inf, box.padding = 0.3
      )
    } else {
      p <- p + geom_text(
        data = df_label,
        aes(label = .data[[label_col]]),
        vjust = -0.6, size = 2.2, show.legend = FALSE
      )
    }
  }

  if (length(color_vals) > 0)
    p <- p + scale_color_manual(values = color_vals, na.value = "grey50")

  p
}

# -------------------------------------------------------------------
# STEP 5: Generate plots for each PC pair
# -------------------------------------------------------------------
pc_pairs <- Filter(
  function(p) p[1] <= max_pc && p[2] <= max_pc,
  list(c(1,2), c(1,3), c(2,3), c(1,4), c(2,4), c(3,4))
)

n_total <- length(pc_pairs) * 2  # cluster + site per pair
n_done  <- 0

for (pp in pc_pairs) {
  pcx <- pp[1]; pcy <- pp[2]
  pair_tag <- sprintf("pc%d_pc%d", pcx, pcy)

  # --- Cluster-colored ---
  p_cl <- make_plot(
    df         = dat,
    pcx        = pcx, pcy = pcy,
    color_col  = "cluster",
    color_vals = POP_COLORS,
    title      = sprintf("PCA %s — %s — by cluster (n=%d)", threshold_tag, pair_tag, n_ind),
    label_col  = "sample_id"
  )
  ggsave(file.path(plots_dir, sprintf("pca_cluster_%s_%s.png", threshold_tag, pair_tag)),
         p_cl, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plots_dir, sprintf("pca_cluster_%s_%s.pdf", threshold_tag, pair_tag)),
         p_cl, width = 8, height = 6)

  # --- Site-colored ---
  p_si <- make_plot(
    df         = dat,
    pcx        = pcx, pcy = pcy,
    color_col  = "site",
    color_vals = site_colors,
    title      = sprintf("PCA %s — %s — by site (n=%d)", threshold_tag, pair_tag, n_ind),
    label_col  = "sample_id"
  )
  ggsave(file.path(plots_dir, sprintf("pca_site_%s_%s.png", threshold_tag, pair_tag)),
         p_si, width = 8, height = 6, dpi = 300)

  n_done <- n_done + 2
}

cat(sprintf("Plots written: %d files → %s\n", n_done, plots_dir))

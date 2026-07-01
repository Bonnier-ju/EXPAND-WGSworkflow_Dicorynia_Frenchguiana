#!/usr/bin/env Rscript
# 07.6-fastStructure_plot.R
# Ancestry bar plots from fastStructure best-replicate .named.Q files.
# Layout mirrors 07.2-admixture_plot.R: single ggplot per K, white
# population-separator lines, site labels at 45° coloured by site palette.
#
# Args: out_dir  map_file  k_min  k_max  bed_base  [log_file]
#   out_dir  : directory containing *.named.Q files
#   map_file : two-column TSV (sample_id  site), no header
#   k_min    : minimum K value
#   k_max    : maximum K value
#   bed_base : prefix of the *.named.Q files (e.g. "fastStructure_input")
#   log_file : optional path to pipeline log (chooseK extraction)

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# Cluster fill palette (up to 10 components)
CLUSTER_COLORS <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#66C2A5", "#FFD700", "#999999"
)

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1) args[1] else stop("out_dir required")
map_file <- if (length(args) >= 2) args[2] else stop("map_file required")
k_min    <- as.integer(if (length(args) >= 3) args[3] else 1)
k_max    <- as.integer(if (length(args) >= 4) args[4] else 5)
bed_base <- if (length(args) >= 5) args[5] else "fastStructure_input"
log_file <- if (length(args) >= 6) args[6] else {
  log_dir <- file.path(dirname(out_dir), "logs")
  logs    <- list.files(log_dir, pattern = "\\.log$", full.names = TRUE)
  if (length(logs) > 0) logs[1] else NA_character_
}

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Site colour palette
# -------------------------------------------------------------------
BASE    <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
col_csv <- file.path(BASE, "metadata/sites_couleurs.csv")
if (file.exists(col_csv)) {
  col_df      <- read.csv(col_csv, stringsAsFactors = FALSE)
  site_colors <- setNames(col_df$couleur_hex, col_df$site)
} else {
  site_colors <- c()
  cat("WARN: sites_couleurs.csv not found — site labels will be black\n")
}

# -------------------------------------------------------------------
# Sample -> site map
# -------------------------------------------------------------------
map_df <- read.table(map_file, header = FALSE, sep = "\t",
                     col.names = c("sample_id", "site"),
                     stringsAsFactors = FALSE)
site_order   <- sort(unique(map_df$site))
map_df$site  <- factor(map_df$site, levels = site_order)

# -------------------------------------------------------------------
# chooseK result from log file
# -------------------------------------------------------------------
choosek_label <- ""
if (!is.na(log_file) && file.exists(log_file) && file.info(log_file)$size > 0) {
  log_lines <- readLines(log_file, warn = FALSE)
  hit_marg  <- grep("maximizes marginal likelihood\\s*=\\s*([0-9]+)",
                    log_lines, value = TRUE)
  hit_comp  <- grep("explain structure\\s*.*=\\s*([0-9]+)",
                    log_lines, value = TRUE)
  km <- if (length(hit_marg) > 0)
          as.integer(gsub(".*=\\s*([0-9]+).*", "\\1", hit_marg[1])) else NA
  kc <- if (length(hit_comp) > 0)
          as.integer(gsub(".*=\\s*([0-9]+).*", "\\1", hit_comp[1])) else NA
  if (!is.na(km))
    choosek_label <- sprintf("chooseK: marginal=%d, structure=%d", km, kc)
}
cat(sprintf("INFO: chooseK result: %s\n",
            if (nzchar(choosek_label)) choosek_label else "not found"))

# -------------------------------------------------------------------
# One barplot per K  (mirrors build_admixture_plot in 07.2-admixture_plot.R)
# -------------------------------------------------------------------
for (K in seq(k_min, k_max)) {
  q_file <- file.path(out_dir, sprintf("%s.%d.named.Q", bed_base, K))
  if (!file.exists(q_file)) {
    cat(sprintf("WARN: %s not found — skipping K=%d\n", q_file, K))
    next
  }

  # Read named.Q (IID + K proportion columns, tab-separated, no header)
  cluster_cols <- paste0("cluster_", seq_len(K))
  q_raw <- read.table(q_file, header = FALSE, sep = "",
                      stringsAsFactors = FALSE)
  colnames(q_raw) <- c("sample_id", cluster_cols)

  # Merge with site map
  q_df <- merge(q_raw, map_df, by = "sample_id", all.x = TRUE, sort = FALSE)
  q_df$site <- as.character(q_df$site)
  q_df$site[is.na(q_df$site)] <- "Unknown"
  q_df$site <- factor(q_df$site, levels = c(site_order, "Unknown"))

  # Sort: site, then dominant component descending within site
  q_df$dominant_cluster <- apply(q_df[, cluster_cols, drop = FALSE], 1, which.max)
  q_df$max_prop         <- apply(q_df[, cluster_cols, drop = FALSE], 1, max)
  q_df <- q_df[order(q_df$site, q_df$dominant_cluster, -q_df$max_prop), ]
  q_df$ind_order <- seq_len(nrow(q_df))

  # Site boundary positions (separator lines + label positions)
  site_order_present <- unique(as.character(q_df$site))
  site_bounds <- do.call(rbind, lapply(site_order_present, function(s) {
    idx <- which(as.character(q_df$site) == s)
    data.frame(site = s,
               xmin = min(idx) - 0.5,
               xmid = mean(idx),
               stringsAsFactors = FALSE)
  }))

  # Long format
  q_long <- pivot_longer(q_df,
                         cols      = all_of(cluster_cols),
                         names_to  = "cluster",
                         values_to = "proportion")
  q_long$cluster <- factor(q_long$cluster, levels = cluster_cols)

  n_ind <- nrow(q_df)

  p <- ggplot(q_long, aes(x = ind_order, y = proportion, fill = cluster)) +
    geom_col(width = 1, position = "stack") +
    geom_vline(data    = site_bounds[-1, ],
               mapping = aes(xintercept = xmin),
               color   = "white", linewidth = 0.6) +
    annotate("text",
             x      = site_bounds$xmid,
             y      = -0.05,
             label  = site_bounds$site,
             colour = "black",
             angle  = 45, hjust = 1, vjust = 1, size = 2.4) +
    scale_fill_manual(
      values = setNames(CLUSTER_COLORS[seq_len(K)], cluster_cols),
      labels = paste0("K", seq_len(K)),
      name   = "Component"
    ) +
    scale_x_continuous(breaks = NULL, expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.35, 1.02)) +
    labs(
      title    = sprintf("fastStructure — K=%d intra-cluster  (n=%d)", K, n_ind),
      subtitle = if (nzchar(choosek_label)) choosek_label else NULL,
      x = NULL, y = "Ancestry proportion"
    ) +
    coord_cartesian(clip = "off") +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      panel.grid      = element_blank(),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      plot.subtitle   = element_text(size = 8, color = "grey40"),
      plot.margin     = margin(t = 5, r = 10, b = 55, l = 5)
    )

  png_path <- file.path(plot_dir, sprintf("fastStructure_K%d.png", K))
  ggsave(png_path, p,
         width  = max(12, n_ind * 0.07),
         height = 4,
         dpi    = 300)
  cat(sprintf("Plot written: %s\n", png_path))
}

cat("DONE 07.6-fastStructure_plot.R\n")

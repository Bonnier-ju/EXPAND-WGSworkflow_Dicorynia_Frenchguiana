#!/usr/bin/env Rscript
# 07.6-fastStructure_plot.R
# Produces ancestry barplots from fastStructure .named.Q files.
# Individuals are sorted by site; sites are annotated with a colour strip.
# chooseK model-selection result is shown as a title annotation.
# Mirrors the visual output of 07.2-admixture_plot.R.
#
# Args: out_dir  map_file  k_min  k_max  bed_base  [log_file]
#   out_dir  : directory containing *.named.Q and the log file
#   map_file : two-column TSV (sample_id  site), no header
#   k_min    : minimum K value
#   k_max    : maximum K value
#   bed_base : prefix of the *.named.Q files (e.g. "fastStructure_input")
#   log_file : optional path to pipeline log; used to extract chooseK result

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1) args[1] else stop("out_dir required")
map_file <- if (length(args) >= 2) args[2] else stop("map_file required")
k_min    <- as.integer(if (length(args) >= 3) args[3] else 1)
k_max    <- as.integer(if (length(args) >= 4) args[4] else 5)
bed_base <- if (length(args) >= 5) args[5] else "fastStructure_input"
log_file <- if (length(args) >= 6) args[6] else
              file.path(dirname(out_dir), "logs",
                        list.files(file.path(dirname(out_dir), "logs"),
                                   pattern = "\\.log$")[1])

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Site colour palette
# -------------------------------------------------------------------
BASE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
col_csv <- file.path(BASE, "metadata/sites_couleurs.csv")
if (file.exists(col_csv)) {
  col_df    <- read.csv(col_csv, stringsAsFactors = FALSE)
  site_colors <- setNames(col_df$couleur_hex, col_df$site)
} else {
  site_colors <- c()
  cat("WARN: sites_couleurs.csv not found — sites will use default colours\n")
}

# Cluster fill palette (up to 5 components)
CLUSTER_COLORS <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# -------------------------------------------------------------------
# Sample -> site map
# -------------------------------------------------------------------
map_df <- read.table(map_file, header = FALSE, sep = "\t",
                     col.names = c("sample_id", "site"),
                     stringsAsFactors = FALSE)

# Site ordering: alphabetical within sites, sites ordered by their
# first occurrence in the map file (preserves geographic grouping)
site_order <- unique(map_df$site)
map_df$site <- factor(map_df$site, levels = site_order)

# -------------------------------------------------------------------
# chooseK result from log file
# -------------------------------------------------------------------
choosek_marg <- NA_integer_
choosek_comp <- NA_integer_
if (!is.null(log_file) && !is.na(log_file) &&
    file.exists(log_file) && file.info(log_file)$size > 0) {
  log_lines <- readLines(log_file, warn = FALSE)
  hit_marg <- grep("maximizes marginal likelihood\\s*=\\s*([0-9]+)",
                   log_lines, value = TRUE)
  hit_comp <- grep("explain structure\\s*.*=\\s*([0-9]+)",
                   log_lines, value = TRUE)
  if (length(hit_marg) > 0)
    choosek_marg <- as.integer(gsub(".*=\\s*([0-9]+).*", "\\1", hit_marg[1]))
  if (length(hit_comp) > 0)
    choosek_comp <- as.integer(gsub(".*=\\s*([0-9]+).*", "\\1", hit_comp[1]))
}
choosek_label <- ""
if (!is.na(choosek_marg))
  choosek_label <- sprintf("chooseK: marginal=%d, structure=%d",
                            choosek_marg, choosek_comp)
cat(sprintf("INFO: chooseK result: %s\n",
            if (nzchar(choosek_label)) choosek_label else "not found"))

# -------------------------------------------------------------------
# One barplot per K
# -------------------------------------------------------------------
for (K in seq(k_min, k_max)) {
  q_file <- file.path(out_dir, sprintf("%s.%d.named.Q", bed_base, K))
  if (!file.exists(q_file)) {
    cat(sprintf("WARN: %s not found — skipping K=%d\n", q_file, K))
    next
  }

  # Read named.Q (IID + K proportion columns, no header, space-separated)
  col_names <- c("sample_id", paste0("C", seq_len(K)))
  q_df <- read.table(q_file, header = FALSE, sep = "",
                     col.names = col_names, stringsAsFactors = FALSE)

  # Join with site map
  q_df <- merge(q_df, map_df, by = "sample_id", all.x = TRUE)
  q_df$site[is.na(q_df$site)] <- "Unknown"
  q_df$site <- factor(q_df$site, levels = c(site_order, "Unknown"))

  # Sort: by site first, then by dominant component descending within site
  dom_col <- paste0("C", which.max(colMeans(q_df[, paste0("C", seq_len(K)),
                                                  drop = FALSE])))
  q_df <- q_df[order(q_df$site, -q_df[[dom_col]]), ]
  q_df$ind_order <- seq_len(nrow(q_df))

  # Reshape to long format
  long_df <- reshape(q_df[, c("sample_id", "ind_order", "site",
                               paste0("C", seq_len(K)))],
                     varying   = paste0("C", seq_len(K)),
                     v.names   = "proportion",
                     timevar   = "component",
                     times     = paste0("C", seq_len(K)),
                     direction = "long")
  long_df$component <- factor(long_df$component,
                               levels = paste0("C", seq_len(K)))

  # ---- barplot ----
  p_bar <- ggplot(long_df,
                  aes(x = ind_order, y = proportion, fill = component)) +
    geom_col(width = 1, linewidth = 0) +
    scale_fill_manual(values = setNames(CLUSTER_COLORS[seq_len(K)],
                                        paste0("C", seq_len(K))),
                      name = "Component") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title    = sprintf("fastStructure — K=%d intra-cluster", K),
         subtitle = if (nzchar(choosek_label)) choosek_label else NULL,
         x = NULL, y = "Ancestry proportion") +
    theme(axis.text.x    = element_blank(),
          axis.ticks.x   = element_blank(),
          panel.grid     = element_blank(),
          plot.subtitle  = element_text(size = 8, color = "grey40"),
          legend.position = "right")

  # ---- site colour strip (below barplot) ----
  # One row per individual, coloured by site
  strip_df <- unique(q_df[, c("ind_order", "site")])
  strip_colors <- if (length(site_colors) > 0) site_colors else NULL

  p_strip <- ggplot(strip_df, aes(x = ind_order, y = 1, fill = site)) +
    geom_col(width = 1, linewidth = 0) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    (if (!is.null(strip_colors))
       scale_fill_manual(values = strip_colors, na.value = "grey50")
     else
       scale_fill_discrete()) +
    labs(x = NULL, y = NULL, fill = "Site") +
    theme(axis.text    = element_blank(),
          axis.ticks   = element_blank(),
          panel.grid   = element_blank(),
          legend.text  = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.35, "cm"))

  # ---- site label positions (tick marks at site boundaries) ----
  site_bounds <- aggregate(ind_order ~ site, data = q_df, FUN = mean)
  p_labels <- ggplot(site_bounds, aes(x = ind_order, y = 0.5, label = site)) +
    geom_text(size = 2.5, angle = 45, hjust = 1, vjust = 0.5) +
    scale_x_continuous(limits = range(q_df$ind_order), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void()

  # ---- assemble with gridExtra ----
  n_ind  <- nrow(q_df)
  width  <- max(8, n_ind * 0.07)

  png_path <- file.path(plot_dir, sprintf("fastStructure_K%d.png", K))
  png(png_path, width = width, height = 5.5, units = "in", res = 300)
  grid.arrange(
    p_bar,
    p_strip  + theme(legend.position = "right"),
    p_labels,
    ncol    = 1,
    heights = c(6, 0.8, 1.8)
  )
  dev.off()
  cat(sprintf("Plot written: %s\n", png_path))
}

cat("DONE 07.6-fastStructure_plot.R\n")

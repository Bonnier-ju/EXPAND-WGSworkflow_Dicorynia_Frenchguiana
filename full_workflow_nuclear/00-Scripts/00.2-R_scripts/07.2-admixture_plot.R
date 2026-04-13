#!/usr/bin/env Rscript
# 07.2-admixture_plot.R
# Visualize ADMIXTURE results:
#   1. Cross-validation error plot — best K selection
#   2. STRUCTURE-style ancestry bar plots for each K

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# Cluster color palette — 10 visually distinct colors for K=1..10
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
k_min    <- if (length(args) >= 3) as.integer(args[3]) else 1L
k_max    <- if (length(args) >= 4) as.integer(args[4]) else 10L
bed_base <- if (length(args) >= 5) args[5] else "admixture_input"

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Load population map
# -------------------------------------------------------------------
map <- read.table(map_file, header = FALSE, sep = "\t",
                  col.names = c("sample_id", "population"),
                  stringsAsFactors = FALSE)

# -------------------------------------------------------------------
# STEP 1: Cross-validation error plot
# -------------------------------------------------------------------
cv_file <- file.path(out_dir, "cross_validation.txt")
if (!file.exists(cv_file)) stop("cross_validation.txt not found: ", cv_file)

cv_lines <- readLines(cv_file)
cv_df <- do.call(rbind, lapply(cv_lines, function(line) {
  k_match  <- regmatches(line, regexpr("K=([0-9]+)", line))
  cv_match <- regmatches(line, regexpr("[0-9]+\\.[0-9]+$", line))
  if (length(k_match) == 0 || length(cv_match) == 0) return(NULL)
  data.frame(K = as.integer(sub("K=", "", k_match)),
             cv_error = as.numeric(cv_match),
             stringsAsFactors = FALSE)
}))
cv_df <- cv_df[order(cv_df$K), ]

best_k <- cv_df$K[which.min(cv_df$cv_error)]
cat(sprintf("INFO: best K = %d (CV error = %.5f)\n",
            best_k, min(cv_df$cv_error)))

cv_df$is_best <- cv_df$K == best_k

p_cv <- ggplot(cv_df, aes(x = K, y = cv_error)) +
  geom_line(color = "grey40", linewidth = 0.8) +
  geom_point(aes(color = is_best), size = 3.5) +
  geom_text(aes(label = ifelse(is_best, sprintf("K=%d", K), "")),
            vjust = -1, hjust = 0.5, size = 3.5, color = "red") +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "red"),
                     guide  = "none") +
  scale_x_continuous(breaks = cv_df$K) +
  labs(
    title    = "ADMIXTURE - cross-validation error",
    subtitle = sprintf("Best K = %d  (CV = %.5f)", best_k, min(cv_df$cv_error)),
    x = "K", y = "CV error"
  ) +
  theme(plot.subtitle = element_text(color = "red"))

ggsave(file.path(plot_dir, "admixture_cv_error.png"), p_cv,
       width = 7, height = 4, dpi = 300)
ggsave(file.path(plot_dir, "admixture_cv_error.pdf"), p_cv,
       width = 7, height = 4)
cat("CV error plot written\n")

# -------------------------------------------------------------------
# STEP 2: STRUCTURE-style bar plots for K = k_min..k_max
# sort_mode = "alphabetical" : populations sorted A-Z, individuals
#             within each population sorted by dominant cluster
# sort_mode = "by_cluster"   : populations grouped by dominant cluster
#             (cluster with highest mean ancestry), then A-Z within cluster
# -------------------------------------------------------------------
build_admixture_plot <- function(q_df, cluster_cols, k, sort_mode) {

  # Compute individual-level dominant cluster and max proportion
  q_df$dominant_cluster <- apply(q_df[, cluster_cols, drop = FALSE], 1, which.max)
  q_df$max_prop         <- apply(q_df[, cluster_cols, drop = FALSE], 1, max)

  if (sort_mode == "alphabetical") {
    # Sort: population A-Z, then within population by dominant cluster desc
    q_df <- q_df[order(q_df$population,
                        q_df$dominant_cluster,
                        -q_df$max_prop), ]

  } else {  # by_cluster
    # Compute dominant cluster for each population (cluster with highest mean Q)
    pop_cluster <- tapply(q_df$dominant_cluster, q_df$population,
                          function(x) as.integer(names(which.max(table(x)))))
    q_df$pop_cluster <- pop_cluster[q_df$population]
    # Sort: population dominant cluster, then population A-Z, then individual
    q_df <- q_df[order(q_df$pop_cluster,
                        q_df$population,
                        q_df$dominant_cluster,
                        -q_df$max_prop), ]
    q_df$pop_cluster <- NULL
  }

  q_df$ind_order <- seq_len(nrow(q_df))

  # Population boundary positions (for separator lines and labels)
  pop_order  <- unique(q_df$population)
  pop_bounds <- do.call(rbind, lapply(pop_order, function(p) {
    idx <- which(q_df$population == p)
    data.frame(population = p,
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

  sort_label <- if (sort_mode == "alphabetical") "alphabetical order"
                else "grouped by dominant cluster"

  p <- ggplot(q_long, aes(x = ind_order, y = proportion, fill = cluster)) +
    geom_col(width = 1, position = "stack") +
    geom_vline(data = pop_bounds[-1, ],
               aes(xintercept = xmin),
               color = "white", linewidth = 0.6) +
    annotate("text",
             x = pop_bounds$xmid, y = -0.05,
             label = pop_bounds$population,
             angle = 45, hjust = 1, vjust = 1, size = 2.4) +
    scale_fill_manual(values = CLUSTER_COLORS[seq_len(k)],
                      labels = paste0("K", seq_len(k)),
                      name   = "Cluster") +
    scale_x_continuous(breaks = NULL, expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-0.35, 1.02)) +
    labs(
      title    = sprintf("ADMIXTURE K = %d  (n = %d individuals)", k, nrow(q_df)),
      subtitle = sort_label,
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
  p
}

make_admixture_plot <- function(k) {
  q_file <- file.path(out_dir, sprintf("%s.%d.named.Q", bed_base, k))
  if (!file.exists(q_file)) {
    cat(sprintf("WARN: Q file missing for K=%d, skipping\n", k))
    return(invisible(NULL))
  }

  q_raw <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
  cluster_cols <- paste0("cluster_", seq_len(k))
  colnames(q_raw) <- c("sample_id", cluster_cols)

  q_df <- merge(q_raw, map, by = "sample_id", all.x = TRUE, sort = FALSE)
  q_df$population[is.na(q_df$population) | q_df$population == ""] <- "UNKNOWN"

  pw <- max(12, nrow(q_df) * 0.07)

  for (smode in c("alphabetical", "by_cluster")) {
    p       <- build_admixture_plot(q_df, cluster_cols, k, smode)
    suffix  <- if (smode == "alphabetical") "alpha" else "bycluster"
    fn_base <- file.path(plot_dir,
                         sprintf("admixture_barplot_K%02d_%s", k, suffix))
    ggsave(paste0(fn_base, ".png"), p, width = pw, height = 4, dpi = 300)
    ggsave(paste0(fn_base, ".pdf"), p, width = pw, height = 4)
  }
  cat(sprintf("K=%d bar plots written (alphabetical + by_cluster)\n", k))
}

for (k in seq(k_min, k_max)) {
  make_admixture_plot(k)
}

# -------------------------------------------------------------------
# STEP 3: Combined multi-K panel (best K ± 2) — two sort modes
# -------------------------------------------------------------------
k_panel <- sort(unique(c(max(k_min, best_k - 2) : min(k_max, best_k + 2))))
k_panel <- k_panel[k_panel >= k_min & k_panel <= k_max]

# Build one panel strip per K for a given sort_mode
build_panel_strip <- function(k, sort_mode) {
  q_file <- file.path(out_dir, sprintf("%s.%d.named.Q", bed_base, k))
  if (!file.exists(q_file)) return(NULL)

  q_raw <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
  cluster_cols <- paste0("cluster_", seq_len(k))
  colnames(q_raw) <- c("sample_id", cluster_cols)

  q_df <- merge(q_raw, map, by = "sample_id", all.x = TRUE, sort = FALSE)
  q_df$population[is.na(q_df$population) | q_df$population == ""] <- "UNKNOWN"

  # Apply sorting
  q_df$dominant_cluster <- apply(q_df[, cluster_cols, drop = FALSE], 1, which.max)
  q_df$max_prop         <- apply(q_df[, cluster_cols, drop = FALSE], 1, max)

  if (sort_mode == "alphabetical") {
    q_df <- q_df[order(q_df$population, q_df$dominant_cluster, -q_df$max_prop), ]
  } else {
    pop_cluster <- tapply(q_df$dominant_cluster, q_df$population,
                          function(x) as.integer(names(which.max(table(x)))))
    q_df$pop_cluster <- pop_cluster[q_df$population]
    q_df <- q_df[order(q_df$pop_cluster, q_df$population,
                        q_df$dominant_cluster, -q_df$max_prop), ]
    q_df$pop_cluster <- NULL
  }

  q_df$ind_order <- seq_len(nrow(q_df))
  show_labels    <- (k == max(k_panel))
  y_floor        <- if (show_labels) -0.35 else -0.02

  pop_order  <- unique(q_df$population)
  pop_bounds <- do.call(rbind, lapply(pop_order, function(p) {
    idx <- which(q_df$population == p)
    data.frame(population = p, xmin = min(idx) - 0.5, xmid = mean(idx),
               stringsAsFactors = FALSE)
  }))

  q_long <- pivot_longer(q_df, cols = all_of(cluster_cols),
                          names_to = "cluster", values_to = "proportion")
  q_long$cluster <- factor(q_long$cluster, levels = cluster_cols)

  p <- ggplot(q_long, aes(x = ind_order, y = proportion, fill = cluster)) +
    geom_col(width = 1, position = "stack") +
    geom_vline(data = pop_bounds[-1, ], aes(xintercept = xmin),
               color = "white", linewidth = 0.5)

  if (show_labels) {
    p <- p + annotate("text",
                      x = pop_bounds$xmid, y = -0.05,
                      label = pop_bounds$population,
                      angle = 45, hjust = 1, vjust = 1, size = 2.2)
  }

  p +
    scale_fill_manual(values = CLUSTER_COLORS[seq_len(k)]) +
    scale_x_continuous(breaks = NULL, expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(y_floor, 1.02)) +
    labs(title = NULL, x = NULL,
         y = sprintf("K = %d%s", k, if (k == best_k) " *" else "")) +
    coord_cartesian(clip = "off") +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      panel.grid      = element_blank(),
      legend.position = "none",
      axis.title.y    = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 1),
      plot.margin     = margin(t = 2, r = 5,
                               b = if (show_labels) 50 else 2, l = 5)
    )
}

# Write combined panel for a given sort_mode
write_panel <- function(sort_mode) {
  panel_plots <- Filter(Negate(is.null),
                        lapply(k_panel, build_panel_strip, sort_mode = sort_mode))
  if (length(panel_plots) == 0) return(invisible(NULL))

  suffix       <- if (sort_mode == "alphabetical") "alpha" else "bycluster"
  library(grid)
  n_panels     <- length(panel_plots)
  pw           <- max(12, nrow(map) * 0.07)
  ph_per_panel <- 2.5
  ph_total     <- ph_per_panel * n_panels + 0.6

  for (ext in c("png", "pdf")) {
    fn <- file.path(plot_dir,
                    sprintf("admixture_panel_bestK_%s.%s", suffix, ext))
    if (ext == "png") {
      png(fn, width = round(pw * 150), height = round(ph_total * 150), res = 150)
    } else {
      pdf(fn, width = pw, height = ph_total)
    }
    grid.newpage()
    pushViewport(viewport(x = 0.5, y = 1, width = 1, height = 0.6 / ph_total,
                          just = c("center", "top")))
    sort_label <- if (sort_mode == "alphabetical") "alphabetical order"
                  else "grouped by dominant cluster"
    grid.text(
      sprintf("ADMIXTURE - K = %s (best K = %d marked *) - %s",
              paste(k_panel, collapse = ", "), best_k, sort_label),
      gp = gpar(fontface = "bold", fontsize = 10)
    )
    popViewport()
    for (i in seq_along(panel_plots)) {
      vp <- viewport(
        x      = 0.5,
        y      = 1 - 0.6 / ph_total - (i - 0.5) * ph_per_panel / ph_total,
        width  = 1,
        height = ph_per_panel / ph_total
      )
      print(panel_plots[[i]], vp = vp)
    }
    dev.off()
    cat(sprintf("Combined panel written: %s\n", fn))
  }
}

write_panel("alphabetical")
write_panel("by_cluster")

cat(sprintf("DONE 07.2 ADMIXTURE plots written to: %s\n", plot_dir))

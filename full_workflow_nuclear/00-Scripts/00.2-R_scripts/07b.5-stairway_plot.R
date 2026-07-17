#!/usr/bin/env Rscript
# 07b.5-stairway_plot.R
# Visualize Stairway Plot 2 demographic inference results for K=3 genetic
# clusters (Pop_1/Pop_2/Pop_3) across purity thresholds (T70/T80/T90/T95).
#
# Called by 07b.5-stairway_plot.slurm with args: out_base  tag  mu_ref  year_per_gen
#
# When tag is a specific threshold (e.g. "T80"):
#   Figure 1 — Ne(t): 3 populations overlaid, 95% CI ribbon, log-log axes.
#   Figure 2 — Mu sensitivity: 3-panel grid, mu = 2..8 × 10^-9.
#   Figure 3 — Ne(t) with total uncertainty band (CI × mu range).
#
# When tag == "ALL":
#   Figure 4 — Cross-threshold Ne(t): color by population, linetype by threshold.
#
# Mu rescaling: when mu changes, Ne and time both scale by mu_ref / mu_alt.
#
# Args: out_base  tag  mu_ref  year_per_gen

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
args         <- commandArgs(trailingOnly = TRUE)
out_base     <- if (length(args) >= 1) args[1] else stop("out_base required")
tag          <- if (length(args) >= 2) args[2] else "ALL"
mu_ref       <- as.numeric(if (length(args) >= 3) args[3] else 5e-9)
year_per_gen <- as.numeric(if (length(args) >= 4) args[4] else 30)

POPS       <- c("Pop_1", "Pop_2", "Pop_3")
THRESHOLDS <- c("T70", "T80", "T90", "T95")
MU_VALUES  <- c(2, 3, 4, 5, 6, 7, 8) * 1e-9

POP_COLORS <- c(Pop_1 = "#EE7600", Pop_2 = "#458B00", Pop_3 = "#CD2626")

cat(sprintf("INFO: out_base = %s\n", out_base))
cat(sprintf("INFO: tag      = %s\n", tag))

# -------------------------------------------------------------------
# Load a Stairway Plot .final.summary file for one population × threshold.
# Returns a data frame with population and threshold columns appended,
# or NULL if the file is missing.
# -------------------------------------------------------------------
load_summary <- function(pop, tg, out_base) {
  f <- file.path(out_base, tg, pop, "stairway",
                 paste0(pop, " (", tg, ").final.summary"))
  if (!file.exists(f)) {
    cat(sprintf("WARN: summary not found for %s %s: %s\n", pop, tg, f))
    return(NULL)
  }
  df <- tryCatch(
    read.table(f, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      cat(sprintf("ERROR reading %s: %s\n", f, e$message))
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  colnames(df)[colnames(df) == "Ne_median"]  <- "Ne_med"
  colnames(df)[colnames(df) == "Ne_2.5%"]   <- "Ne_lo"
  colnames(df)[colnames(df) == "Ne_97.5%"]  <- "Ne_hi"
  colnames(df)[colnames(df) == "Ne_12.5%"]  <- "Ne_lo2"
  colnames(df)[colnames(df) == "Ne_87.5%"]  <- "Ne_hi2"

  # Remove artifact row (year ~ 0 due to floating-point underflow)
  df <- df[!is.na(df$year) & df$year >= 1, ]
  df$population <- pop
  df$threshold  <- tg
  df
}

# -------------------------------------------------------------------
# Load all populations for one or more thresholds
# -------------------------------------------------------------------
load_all <- function(tags) {
  dfs <- list()
  for (tg in tags) {
    for (pop in POPS) {
      d <- load_summary(pop, tg, out_base)
      if (!is.null(d)) dfs[[length(dfs) + 1]] <- d
    }
  }
  if (length(dfs) == 0) return(NULL)
  do.call(rbind, dfs)
}

# -------------------------------------------------------------------
# Mu sensitivity rescaling helper
# -------------------------------------------------------------------
expand_mu <- function(df, mu_values, mu_ref) {
  do.call(rbind, lapply(mu_values, function(mu_alt) {
    f           <- mu_ref / mu_alt
    df$year_s   <- df$year   * f
    df$Ne_s     <- df$Ne_med * f
    df$Ne_lo_s  <- df$Ne_lo  * f
    df$Ne_hi_s  <- df$Ne_hi  * f
    df$mu_alt   <- mu_alt
    df
  }))
}

# -------------------------------------------------------------------
# Figures 1-3: per-threshold plots (3 populations overlaid)
# -------------------------------------------------------------------
make_per_threshold_plots <- function(combined_df, tg, plot_dir) {

  combined_df$population <- factor(combined_df$population, levels = POPS)

  # Figure 1: combined Ne(t), 95% CI ribbon
  p1 <- ggplot(combined_df,
               aes(x = year, color = population, fill = population)) +
    geom_ribbon(aes(ymin = Ne_lo, ymax = Ne_hi),
                alpha = 0.15, color = NA) +
    geom_step(aes(y = Ne_med), linewidth = 0.8, direction = "hv") +
    scale_color_manual(values = POP_COLORS, name = "Cluster") +
    scale_fill_manual(values  = POP_COLORS, name = "Cluster") +
    scale_x_log10(name = "Time (years before present)") +
    scale_y_log10(name = expression(paste("Effective population size (", N[e], ")"))) +
    annotation_logticks(sides = "bl", size = 0.3) +
    labs(
      title    = sprintf("Demographic history — Stairway Plot 2 (K=%d, %s)", 3L, tg),
      subtitle = sprintf(
        "Folded SFS | L = 194–221 Mbp callable (per cluster) | μ = %.0e | %d yr/gen | 95%% CI ribbon",
        mu_ref, as.integer(year_per_gen)
      )
    ) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )

  ggsave(file.path(plot_dir, sprintf("Ne_trajectories_%s.png", tg)), p1,
         width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, sprintf("Ne_trajectories_%s.pdf", tg)), p1,
         width = 10, height = 6)
  cat(sprintf("Figure 1 written for %s\n", tg))

  # Figure 2: mu sensitivity, one panel per population
  mu_df <- expand_mu(combined_df, MU_VALUES, mu_ref)
  mu_df$population <- factor(mu_df$population, levels = POPS)
  mu_df$mu_label   <- sprintf("%.0e", mu_df$mu_alt)
  mu_df$is_ref     <- abs(mu_df$mu_alt - mu_ref) < 1e-15

  mu_line_colors <- setNames(rep("grey70", length(MU_VALUES)),
                             sprintf("%.0e", MU_VALUES))
  mu_line_colors[sprintf("%.0e", mu_ref)] <- "black"

  p2 <- ggplot(mu_df,
               aes(x = year_s, y = Ne_s,
                   color = mu_label, linewidth = is_ref, group = mu_label)) +
    geom_step(direction = "hv") +
    scale_color_manual(
      values = mu_line_colors,
      name   = expression(paste(mu, " (sub/site/gen)"))
    ) +
    scale_linewidth_manual(values = c("TRUE" = 1.1, "FALSE" = 0.45), guide = "none") +
    scale_x_log10(name = "Time (years before present)") +
    scale_y_log10(name = expression(N[e])) +
    annotation_logticks(sides = "bl", size = 0.25) +
    facet_wrap(~ population, ncol = 3, scales = "free") +
    labs(
      title    = sprintf("Mu sensitivity — Stairway Plot 2 (K=3, %s)", tg),
      subtitle = "Reference mu (5×10⁻⁹) in black | grey: mu = 2–8 × 10⁻⁹"
    ) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 8),
      legend.title     = element_text(size = 8),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )

  ggsave(file.path(plot_dir, sprintf("Ne_mu_sensitivity_%s.png", tg)), p2,
         width = 12, height = 5, dpi = 300)
  ggsave(file.path(plot_dir, sprintf("Ne_mu_sensitivity_%s.pdf", tg)), p2,
         width = 12, height = 5)
  cat(sprintf("Figure 2 written for %s\n", tg))

  # Figure 3: Ne(t) with total uncertainty band (CI + mu range)
  mu_min <- min(MU_VALUES)
  mu_max <- max(MU_VALUES)

  band_df <- do.call(rbind, lapply(POPS, function(pop) {
    d <- combined_df[combined_df$population == pop, ]
    if (nrow(d) == 0) return(NULL)
    data.frame(
      population = pop,
      year       = d$year,
      Ne_med     = d$Ne_med,
      Ne_lo_ci   = d$Ne_lo,
      Ne_hi_ci   = d$Ne_hi,
      Ne_lo_band = d$Ne_lo * (mu_ref / mu_max),
      Ne_hi_band = d$Ne_hi * (mu_ref / mu_min)
    )
  }))
  band_df$population <- factor(band_df$population, levels = POPS)

  p3 <- ggplot(band_df, aes(color = population, fill = population, group = population)) +
    geom_ribbon(aes(x = year, ymin = Ne_lo_band, ymax = Ne_hi_band),
                alpha = 0.08, color = NA) +
    geom_ribbon(aes(x = year, ymin = Ne_lo_ci, ymax = Ne_hi_ci),
                alpha = 0.20, color = NA) +
    geom_step(aes(x = year, y = Ne_med), linewidth = 0.9, direction = "hv") +
    scale_color_manual(values = POP_COLORS, name = "Cluster") +
    scale_fill_manual(values  = POP_COLORS, name = "Cluster") +
    scale_x_log10(name = "Time (years before present)") +
    scale_y_log10(name = expression(paste("Effective population size (", N[e], ")"))) +
    annotation_logticks(sides = "bl", size = 0.3) +
    labs(
      title    = sprintf("Demographic history with total uncertainty — K=3, %s", tg),
      subtitle = paste0(
        "Solid = median (μ=5×10⁻⁹) | dark ribbon = 95% CI | ",
        "light band = CI × μ range (2–8×10⁻⁹)"
      )
    ) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )

  ggsave(file.path(plot_dir, sprintf("Ne_trajectories_with_mu_band_%s.png", tg)), p3,
         width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, sprintf("Ne_trajectories_with_mu_band_%s.pdf", tg)), p3,
         width = 10, height = 6)
  cat(sprintf("Figure 3 written for %s\n", tg))
}

# -------------------------------------------------------------------
# Figure 4: Cross-threshold comparison (all 4 thresholds × 3 populations)
# Faceted by population, color and linetype encode threshold.
# -------------------------------------------------------------------
make_cross_threshold_plot <- function(combined_df, plot_dir) {
  if (is.null(combined_df) || nrow(combined_df) == 0) return(invisible(NULL))

  combined_df$population <- factor(combined_df$population, levels = POPS)
  combined_df$threshold  <- factor(combined_df$threshold,  levels = THRESHOLDS)

  THRESH_LINES <- c(T70 = "solid", T80 = "dashed", T90 = "dotdash", T95 = "dotted")
  THRESH_ALPHA <- c(T70 = 0.20, T80 = 0.18, T90 = 0.15, T95 = 0.12)

  p4 <- ggplot(combined_df,
               aes(x = year, y = Ne_med,
                   color = population, linetype = threshold,
                   fill  = population, group = interaction(population, threshold))) +
    geom_ribbon(aes(ymin = Ne_lo, ymax = Ne_hi, alpha = threshold),
                color = NA, show.legend = FALSE) +
    geom_step(linewidth = 0.7, direction = "hv") +
    scale_color_manual(values = POP_COLORS, name = "Cluster") +
    scale_fill_manual(values  = POP_COLORS, name = "Cluster") +
    scale_linetype_manual(values = THRESH_LINES, name = "Threshold") +
    scale_alpha_manual(values = THRESH_ALPHA) +
    scale_x_log10(name = "Time (years before present)") +
    scale_y_log10(name = expression(paste("Effective population size (", N[e], ")"))) +
    annotation_logticks(sides = "bl", size = 0.25) +
    facet_wrap(~ population, ncol = 3, scales = "free_y") +
    labs(
      title    = "Demographic history across purity thresholds — Stairway Plot 2 (K=3)",
      subtitle = sprintf(
        "Folded SFS | μ = %.0e | %d yr/gen | 95%% CI ribbon | linetype = purity threshold",
        mu_ref, as.integer(year_per_gen)
      )
    ) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    ) +
    guides(color     = guide_legend(order = 1),
           fill      = "none",
           linetype  = guide_legend(order = 2))

  ggsave(file.path(plot_dir, "Ne_trajectories_all_thresholds.png"), p4,
         width = 14, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "Ne_trajectories_all_thresholds.pdf"), p4,
         width = 14, height = 6)
  cat(sprintf("Figure 4 written: %s\n",
              file.path(plot_dir, "Ne_trajectories_all_thresholds.png")))

  # Additional: overlay all pops × thresholds without facets (color = pop, linetype = threshold)
  p4b <- ggplot(combined_df,
                aes(x = year, y = Ne_med,
                    color = population, linetype = threshold,
                    group = interaction(population, threshold))) +
    geom_step(linewidth = 0.65, direction = "hv") +
    scale_color_manual(values = POP_COLORS, name = "Cluster") +
    scale_linetype_manual(values = THRESH_LINES, name = "Threshold") +
    scale_x_log10(name = "Time (years before present)") +
    scale_y_log10(name = expression(paste(N[e]))) +
    annotation_logticks(sides = "bl", size = 0.25) +
    labs(
      title    = "Demographic history overlay — Stairway Plot 2 (K=3)",
      subtitle = sprintf("μ = %.0e | %d yr/gen | all thresholds and clusters", mu_ref,
                         as.integer(year_per_gen))
    ) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )

  ggsave(file.path(plot_dir, "Ne_trajectories_overlay.png"), p4b,
         width = 10, height = 6, dpi = 300)
  cat("Figure 4b (overlay) written\n")
}

# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
if (tag != "ALL") {
  # Single-threshold mode
  cat(sprintf("Mode: per-threshold (%s)\n", tag))
  combined_df <- load_all(tag)

  if (is.null(combined_df)) {
    cat(sprintf("ERROR: no summary files found for %s — exiting\n", tag))
    quit(status = 1)
  }

  n_loaded <- length(unique(combined_df$population))
  cat(sprintf("INFO: %d populations loaded for %s\n", n_loaded, tag))

  plot_dir <- file.path(out_base, tag, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  make_per_threshold_plots(combined_df, tag, plot_dir)

  cat(sprintf("DONE 07b.5-stairway_plot.R — %s — plots: %s\n", tag, plot_dir))

} else {
  # Cross-threshold mode
  cat("Mode: cross-threshold (ALL)\n")
  combined_df <- load_all(THRESHOLDS)

  if (is.null(combined_df)) {
    cat("ERROR: no summary files found for any threshold — exiting\n")
    quit(status = 1)
  }

  n_combos <- length(unique(paste(combined_df$threshold, combined_df$population)))
  cat(sprintf("INFO: %d population × threshold combinations loaded\n", n_combos))

  plot_dir_all <- file.path(out_base, "plots")
  dir.create(plot_dir_all, recursive = TRUE, showWarnings = FALSE)

  make_cross_threshold_plot(combined_df, plot_dir_all)

  cat(sprintf("DONE 07b.5-stairway_plot.R — ALL — plots: %s\n", plot_dir_all))
}

#!/usr/bin/env Rscript
# 07.11-stairway_plot.R
# Visualize Stairway Plot 2 demographic inference results for 4 ADMIXTURE
# clusters (K=4). Produces three figures:
#   1. Combined Ne(t): 4 populations overlaid (POP_COLORS), 95% CI ribbon,
#      reference mu, log-log axes.
#   2. Mu sensitivity: 4-panel grid (one per population), Ne(t) for
#      mu = 2..8 × 10^-9; reference mu highlighted.
#   3. Combined Ne(t) with total uncertainty band (CI × mu range).
#
# Mu rescaling: when mu changes, both Ne and time scale by mu_ref / mu_alt
# (both estimated from theta = 4*Ne*mu inferred from SFS).
#
# Args: out_base  k_global  mu_ref  year_per_gen

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
k_global     <- as.integer(if (length(args) >= 2) args[2] else 4)
mu_ref       <- as.numeric(if (length(args) >= 3) args[3] else 5e-9)
year_per_gen <- as.numeric(if (length(args) >= 4) args[4] else 30)

pops <- c("Pop_1", "Pop_2", "Pop_3", "Pop_4")

POP_COLORS <- c(
  Pop_1 = "#EE7600",
  Pop_2 = "#458B00",
  Pop_3 = "#CD2626",
  Pop_4 = "#9A32CD"
)

MU_VALUES <- c(2, 3, 4, 5, 6, 7, 8) * 1e-9

plot_dir <- file.path(out_base, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Load Stairway Plot .final.summary files
# Columns renamed to avoid special characters (%) in aes() calls.
# -------------------------------------------------------------------
load_summary <- function(pop, k_global, out_base) {
  f <- file.path(out_base,
                 paste0(pop, "_K", k_global),
                 "stairway",
                 paste0(pop, " (K=", k_global, ").final.summary"))
  if (!file.exists(f)) {
    cat(sprintf("WARN: summary not found for %s: %s\n", pop, f))
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

  # Rename columns with special characters for safe use in aes()
  colnames(df)[colnames(df) == "Ne_median"]      <- "Ne_med"
  colnames(df)[colnames(df) == "Ne_2.5%"]        <- "Ne_lo"
  colnames(df)[colnames(df) == "Ne_97.5%"]        <- "Ne_hi"
  colnames(df)[colnames(df) == "Ne_12.5%"]        <- "Ne_lo2"
  colnames(df)[colnames(df) == "Ne_87.5%"]        <- "Ne_hi2"

  # Remove artifact row (year=0 due to floating-point underflow in first row)
  df <- df[!is.na(df$year) & df$year >= 1, ]
  df$population <- pop
  df
}

all_df <- lapply(pops, load_summary, k_global = k_global, out_base = out_base)
all_df <- Filter(Negate(is.null), all_df)

if (length(all_df) == 0) {
  cat("ERROR: no summary files found — exiting\n")
  quit(status = 1)
}

combined_df <- do.call(rbind, all_df)
combined_df$population <- factor(combined_df$population, levels = pops)

cat(sprintf("INFO: %d populations loaded (%s)\n",
            length(all_df), paste(sapply(all_df, function(d) d$population[1]), collapse = ", ")))

# -------------------------------------------------------------------
# FIGURE 1: Combined Ne(t) at reference mu
# -------------------------------------------------------------------
cat("Building Figure 1: combined Ne(t)...\n")

p1 <- ggplot(combined_df,
             aes(x = year, color = population, fill = population)) +
  geom_ribbon(aes(ymin = Ne_lo, ymax = Ne_hi),
              alpha = 0.15, color = NA) +
  geom_step(aes(y = Ne_med), linewidth = 0.8, direction = "hv") +
  scale_color_manual(values = POP_COLORS, name = "Population") +
  scale_fill_manual(values  = POP_COLORS, name = "Population") +
  scale_x_log10(name = "Time (years before present)") +
  scale_y_log10(name = expression(paste("Effective population size (", N[e], ")"))) +
  annotation_logticks(sides = "bl", size = 0.3) +
  labs(
    title    = sprintf("Demographic history — Stairway Plot 2 (K=%d)", k_global),
    subtitle = sprintf("Folded SFS | L = 194–221 Mbp callable sites (per-population) | μ = %.0e | %d yr/gen | 95%% CI ribbon",
                       mu_ref, as.integer(year_per_gen))
  ) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(plot_dir, "Ne_trajectories.png"), p1,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "Ne_trajectories.pdf"), p1,
       width = 10, height = 6)
cat("Figure 1 written\n")

# -------------------------------------------------------------------
# FIGURE 2: Mu sensitivity — one panel per population
# -------------------------------------------------------------------
cat("Building Figure 2: mu sensitivity...\n")

expand_mu <- function(pop_df, mu_values, mu_ref) {
  do.call(rbind, lapply(mu_values, function(mu_alt) {
    f           <- mu_ref / mu_alt
    pop_df$year_s  <- pop_df$year   * f
    pop_df$Ne_s    <- pop_df$Ne_med * f
    pop_df$Ne_lo_s <- pop_df$Ne_lo  * f
    pop_df$Ne_hi_s <- pop_df$Ne_hi  * f
    pop_df$mu_alt  <- mu_alt
    pop_df
  }))
}

mu_df <- do.call(rbind, lapply(unique(combined_df$population), function(pop) {
  expand_mu(combined_df[combined_df$population == pop, ], MU_VALUES, mu_ref)
}))
mu_df$population <- factor(mu_df$population, levels = pops)
mu_df$mu_label   <- sprintf("%.0e", mu_df$mu_alt)
mu_df$is_ref     <- abs(mu_df$mu_alt - mu_ref) < 1e-15

mu_line_colors <- setNames(
  rep("grey70", length(MU_VALUES)),
  sprintf("%.0e", MU_VALUES)
)
mu_line_colors[sprintf("%.0e", mu_ref)] <- "black"

p2 <- ggplot(mu_df,
             aes(x = year_s, y = Ne_s,
                 color = mu_label, linewidth = is_ref,
                 group = mu_label)) +
  geom_step(direction = "hv") +
  scale_color_manual(
    values = mu_line_colors,
    name   = expression(paste(mu, " (sub/site/gen)"))
  ) +
  scale_linewidth_manual(values = c("TRUE" = 1.1, "FALSE" = 0.45),
                         guide  = "none") +
  scale_x_log10(name = "Time (years before present)") +
  scale_y_log10(name = expression(N[e])) +
  annotation_logticks(sides = "bl", size = 0.25) +
  facet_wrap(~ population, ncol = 2, scales = "free") +
  labs(
    title    = sprintf("Mu sensitivity — Stairway Plot 2 (K=%d)", k_global),
    subtitle = "Reference mu (5×10⁻⁹) in black | grey lines: mu = 2–8 × 10⁻⁹"
  ) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 8),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 9),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(plot_dir, "Ne_mu_sensitivity.png"), p2,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "Ne_mu_sensitivity.pdf"), p2,
       width = 10, height = 8)
cat("Figure 2 written\n")

# -------------------------------------------------------------------
# FIGURE 3: Combined Ne(t) with total uncertainty band
# Outer (light) band = 97.5% CI at mu_min to 2.5% CI at mu_max.
# Inner (darker) band = 95% CI at reference mu.
# Solid step = median at reference mu.
# -------------------------------------------------------------------
cat("Building Figure 3: combined Ne(t) with mu uncertainty band...\n")

mu_min <- min(MU_VALUES)
mu_max <- max(MU_VALUES)

band_df <- do.call(rbind, lapply(unique(combined_df$population), function(pop) {
  d <- combined_df[combined_df$population == pop, ]
  data.frame(
    population = pop,
    year       = d$year,
    Ne_med     = d$Ne_med,
    Ne_lo_ci   = d$Ne_lo,
    Ne_hi_ci   = d$Ne_hi,
    # Outer band: worst-case combinations of CI and mu
    Ne_lo_band = d$Ne_lo  * (mu_ref / mu_max),
    Ne_hi_band = d$Ne_hi  * (mu_ref / mu_min)
  )
}))
band_df$population <- factor(band_df$population, levels = pops)

p3 <- ggplot(band_df, aes(color = population, fill = population, group = population)) +
  geom_ribbon(aes(x = year, ymin = Ne_lo_band, ymax = Ne_hi_band),
              alpha = 0.08, color = NA) +
  geom_ribbon(aes(x = year, ymin = Ne_lo_ci, ymax = Ne_hi_ci),
              alpha = 0.20, color = NA) +
  geom_step(aes(x = year, y = Ne_med), linewidth = 0.9, direction = "hv") +
  scale_color_manual(values = POP_COLORS, name = "Population") +
  scale_fill_manual(values  = POP_COLORS, name = "Population") +
  scale_x_log10(name = "Time (years before present)") +
  scale_y_log10(name = expression(paste("Effective population size (", N[e], ")"))) +
  annotation_logticks(sides = "bl", size = 0.3) +
  labs(
    title    = sprintf("Demographic history with total uncertainty — Stairway Plot 2 (K=%d)",
                       k_global),
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

ggsave(file.path(plot_dir, "Ne_trajectories_with_mu_band.png"), p3,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "Ne_trajectories_with_mu_band.pdf"), p3,
       width = 10, height = 6)
cat("Figure 3 written\n")

cat(sprintf("DONE 07.11-stairway_plot.R — plots: %s\n", plot_dir))

#!/usr/bin/env Rscript
# 07.12a-smcpp_estimate.R
# Visualises SMC++ per-population Ne(t) estimates for the 4 ADMIXTURE
# clusters (K=4). Produces three figures:
#   1. Combined Ne(t) — 4 populations overlaid (POP_COLORS, log-log)
#   2. Mu sensitivity — 4-panel grid, mu = 2–8 × 10^-9
#   3. Pop_3 comparison — SMC++ vs Stairway Plot 2 on the same axes
#
# Input: CSV files from `smc++ plot -g YEAR_PER_GEN --csv`
#   Columns: label (population), x (years BP), y (Ne)
#
# Args: out_base  stairway_base  k_global  mu_ref  year_per_gen

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
args          <- commandArgs(trailingOnly = TRUE)
out_base      <- if (length(args) >= 1) args[1] else stop("out_base required")
stairway_base <- if (length(args) >= 2) args[2] else stop("stairway_base required")
k_global      <- as.integer(if (length(args) >= 3) args[3] else 4)
mu_ref        <- as.numeric(if (length(args) >= 4) args[4] else 5e-9)
year_per_gen  <- as.numeric(if (length(args) >= 5) args[5] else 30)

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
# Load SMC++ CSV files (one per population)
# smc++ plot --csv produces: label, x (years), y (Ne)
# -------------------------------------------------------------------
load_smcpp <- function(pop, out_base) {
  f <- file.path(out_base, "csv", paste0("Ne_", pop, ".csv"))
  if (!file.exists(f)) {
    cat(sprintf("WARN: SMC++ CSV not found for %s: %s\n", pop, f))
    return(NULL)
  }
  df <- tryCatch(
    read.csv(f, stringsAsFactors = FALSE),
    error = function(e) {
      cat(sprintf("ERROR reading %s: %s\n", f, e$message))
      NULL
    }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Normalise column names (smc++ may output 'label', 'x', 'y')
  names(df) <- tolower(names(df))
  if (!"x" %in% names(df) || !"y" %in% names(df)) {
    cat(sprintf("WARN: unexpected columns in %s: %s\n", f, paste(names(df), collapse = ", ")))
    return(NULL)
  }

  # Remove t=0 artifact and very distant past (>5 My poorly resolved)
  df <- df[!is.na(df$x) & df$x >= 10 & df$x <= 5e6, ]
  df$population <- pop
  df$year       <- df$x
  df$Ne         <- df$y
  df[, c("population", "year", "Ne")]
}

smcpp_list  <- lapply(pops, load_smcpp, out_base = out_base)
names(smcpp_list) <- pops
smcpp_valid <- Filter(Negate(is.null), smcpp_list)

if (length(smcpp_valid) == 0) {
  cat("ERROR: no SMC++ CSV files found — exiting\n")
  quit(status = 1)
}

smcpp_df <- do.call(rbind, smcpp_valid)
smcpp_df$population <- factor(smcpp_df$population, levels = pops)

cat(sprintf("INFO: %d populations loaded from SMC++ (%s)\n",
            length(smcpp_valid),
            paste(names(smcpp_valid), collapse = ", ")))

# -------------------------------------------------------------------
# FIGURE 1: Combined Ne(t) — 4 populations
# -------------------------------------------------------------------
cat("Building Figure 1: combined SMC++ Ne(t)...\n")

p1 <- ggplot(smcpp_df, aes(x = year, y = Ne,
                             color = population, group = population)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = POP_COLORS, name = "Population") +
  scale_x_log10(name = "Time (years before present)") +
  scale_y_log10(name = expression(paste("Effective population size (", N[e], ")"))) +
  annotation_logticks(sides = "bl", size = 0.3) +
  labs(
    title    = sprintf("Demographic history — SMC++ (K=%d)", k_global),
    subtitle = sprintf("Per-population cv | μ = %.0e | r = 3×10⁻⁸ | %d yr/gen | TE-masked",
                       mu_ref, as.integer(year_per_gen))
  ) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(plot_dir, "Ne_smcpp_combined.png"), p1,
       width = 10, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "Ne_smcpp_combined.pdf"), p1,
       width = 10, height = 6)
cat("Figure 1 written\n")

# -------------------------------------------------------------------
# FIGURE 2: Mu sensitivity — 4-panel grid
# Rescaling: when mu changes from mu_ref to mu_alt, both Ne and time
# scale by (mu_ref / mu_alt) because theta = 4 * Ne * mu.
# -------------------------------------------------------------------
cat("Building Figure 2: mu sensitivity...\n")

expand_mu <- function(pop_df, mu_values, mu_ref) {
  do.call(rbind, lapply(mu_values, function(mu_alt) {
    f           <- mu_ref / mu_alt
    d           <- pop_df
    d$year_s    <- d$year * f
    d$Ne_s      <- d$Ne   * f
    d$mu_alt    <- mu_alt
    d
  }))
}

mu_df <- do.call(rbind, lapply(unique(as.character(smcpp_df$population)),
                                function(pop) {
  expand_mu(smcpp_df[smcpp_df$population == pop, ], MU_VALUES, mu_ref)
}))
mu_df$population <- factor(mu_df$population, levels = pops)
mu_df$mu_label   <- sprintf("%.0e", mu_df$mu_alt)
mu_df$is_ref     <- abs(mu_df$mu_alt - mu_ref) < 1e-15

mu_colors <- setNames(rep("grey70", length(MU_VALUES)),
                      sprintf("%.0e", MU_VALUES))
mu_colors[sprintf("%.0e", mu_ref)] <- "black"

p2 <- ggplot(mu_df,
             aes(x = year_s, y = Ne_s,
                 color = mu_label, linewidth = is_ref, group = mu_label)) +
  geom_line() +
  scale_color_manual(values = mu_colors,
                     name   = expression(paste(mu, " (sub/site/gen)"))) +
  scale_linewidth_manual(values = c("TRUE" = 1.1, "FALSE" = 0.45),
                         guide  = "none") +
  scale_x_log10(name = "Time (years before present)") +
  scale_y_log10(name = expression(N[e])) +
  annotation_logticks(sides = "bl", size = 0.25) +
  facet_wrap(~ population, ncol = 2, scales = "free") +
  labs(
    title    = sprintf("Mu sensitivity — SMC++ (K=%d)", k_global),
    subtitle = "Reference mu (5×10⁻⁹) in black | grey: mu = 2–8 × 10⁻⁹"
  ) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 8),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 9),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(plot_dir, "Ne_smcpp_mu_sensitivity.png"), p2,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "Ne_smcpp_mu_sensitivity.pdf"), p2,
       width = 10, height = 8)
cat("Figure 2 written\n")

# -------------------------------------------------------------------
# FIGURE 3: Pop_3 comparison — SMC++ vs Stairway Plot 2
# Overlays the two independent Ne(t) inferences to assess concordance.
# Stairway Plot provides a 95% CI ribbon; SMC++ provides a point estimate.
# -------------------------------------------------------------------
cat("Building Figure 3: Pop_3 SMC++ vs Stairway Plot comparison...\n")

load_stairway_pop3 <- function(stairway_base, k_global) {
  f <- file.path(stairway_base,
                 paste0("Pop_3_K", k_global),
                 "stairway",
                 paste0("Pop_3 (K=", k_global, ").final.summary"))
  if (!file.exists(f)) {
    cat(sprintf("WARN: Stairway summary not found: %s\n", f))
    return(NULL)
  }
  df <- tryCatch(
    read.table(f, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) { cat(sprintf("ERROR: %s\n", e$message)); NULL }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  names(df)[names(df) == "Ne_median"] <- "Ne_med"
  names(df)[names(df) == "Ne_2.5%"]  <- "Ne_lo"
  names(df)[names(df) == "Ne_97.5%"] <- "Ne_hi"
  df <- df[!is.na(df$year) & df$year >= 10, ]
  df$method <- "Stairway Plot 2"
  df[, c("year", "Ne_med", "Ne_lo", "Ne_hi", "method")]
}

smcpp_p3 <- smcpp_valid[["Pop_3"]]
sw_p3    <- load_stairway_pop3(stairway_base, k_global)

if (!is.null(smcpp_p3) && !is.null(sw_p3)) {

  # SMC++ data frame for the comparison
  smcpp_p3$method <- "SMC++"
  smcpp_p3$Ne_med <- smcpp_p3$Ne
  smcpp_p3$Ne_lo  <- NA_real_
  smcpp_p3$Ne_hi  <- NA_real_

  p3 <- ggplot() +
    # Stairway Plot 95% CI ribbon
    geom_ribbon(data = sw_p3,
                aes(x = year, ymin = Ne_lo, ymax = Ne_hi),
                fill = "#CD2626", alpha = 0.15) +
    # Stairway Plot median
    geom_step(data = sw_p3,
              aes(x = year, y = Ne_med, linetype = "Stairway Plot 2"),
              color = "#CD2626", linewidth = 0.8, direction = "hv") +
    # SMC++ trajectory
    geom_line(data = smcpp_p3,
              aes(x = year, y = Ne_med, linetype = "SMC++"),
              color = "#CD2626", linewidth = 1.1) +
    scale_linetype_manual(
      values = c("Stairway Plot 2" = "dashed", "SMC++" = "solid"),
      name   = "Method"
    ) +
    scale_x_log10(name = "Time (years before present)") +
    scale_y_log10(name = expression(paste("Ne — Pop_3"))) +
    annotation_logticks(sides = "bl", size = 0.3) +
    labs(
      title    = sprintf("Pop_3 demographic history — SMC++ vs Stairway Plot 2 (K=%d)", k_global),
      subtitle = paste0(
        "Solid = SMC++ cv | Dashed = Stairway Plot 2 median | ",
        "Ribbon = Stairway Plot 95% CI | μ = 5×10⁻⁹ | 30 yr/gen"
      )
    ) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )

  ggsave(file.path(plot_dir, "Ne_Pop3_smcpp_vs_stairway.png"), p3,
         width = 10, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "Ne_Pop3_smcpp_vs_stairway.pdf"), p3,
         width = 10, height = 6)
  cat("Figure 3 written\n")
} else {
  cat("WARN: Figure 3 skipped (SMC++ or Stairway data missing for Pop_3)\n")
}

# -------------------------------------------------------------------
# FIGURE 4: All populations — SMC++ vs Stairway Plot 2 (4-panel)
# One panel per population; SMC++ solid, Stairway dashed + CI ribbon.
# -------------------------------------------------------------------
cat("Building Figure 4: all populations SMC++ vs Stairway comparison...\n")

load_stairway_pop <- function(pop, stairway_base, k_global) {
  f <- file.path(stairway_base,
                 paste0(pop, "_K", k_global),
                 "stairway",
                 paste0(pop, " (K=", k_global, ").final.summary"))
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(
    read.table(f, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  names(df)[names(df) == "Ne_median"] <- "Ne_med"
  names(df)[names(df) == "Ne_2.5%"]  <- "Ne_lo"
  names(df)[names(df) == "Ne_97.5%"] <- "Ne_hi"
  df <- df[!is.na(df$year) & df$year >= 10, ]
  df$population <- pop
  df[, c("population", "year", "Ne_med", "Ne_lo", "Ne_hi")]
}

sw_all <- do.call(rbind, Filter(Negate(is.null),
                                lapply(pops, load_stairway_pop,
                                       stairway_base = stairway_base,
                                       k_global      = k_global)))

if (!is.null(sw_all) && nrow(sw_all) > 0 && length(smcpp_valid) > 0) {
  sw_all$population <- factor(sw_all$population, levels = pops)

  p4 <- ggplot() +
    geom_ribbon(data = sw_all,
                aes(x = year, ymin = Ne_lo, ymax = Ne_hi,
                    fill = population),
                alpha = 0.18, color = NA) +
    geom_step(data = sw_all,
              aes(x = year, y = Ne_med, color = population,
                  linetype = "Stairway Plot 2"),
              linewidth = 0.7, direction = "hv") +
    geom_line(data = smcpp_df,
              aes(x = year, y = Ne, color = population,
                  linetype = "SMC++"),
              linewidth = 1.0) +
    scale_color_manual(values = POP_COLORS, guide = "none") +
    scale_fill_manual(values  = POP_COLORS, guide = "none") +
    scale_linetype_manual(
      values = c("Stairway Plot 2" = "dashed", "SMC++" = "solid"),
      name   = "Method"
    ) +
    scale_x_log10(name = "Time (years before present)") +
    scale_y_log10(name = expression(N[e])) +
    annotation_logticks(sides = "bl", size = 0.25) +
    facet_wrap(~ population, ncol = 2, scales = "free_y") +
    labs(
      title    = sprintf("SMC++ vs Stairway Plot 2 — all clusters (K=%d)", k_global),
      subtitle = "Solid = SMC++ | Dashed = Stairway median | Ribbon = 95% CI | μ = 5×10⁻⁹"
    ) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = 9, color = "black"),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )

  ggsave(file.path(plot_dir, "Ne_all_smcpp_vs_stairway.png"), p4,
         width = 11, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "Ne_all_smcpp_vs_stairway.pdf"), p4,
         width = 11, height = 8)
  cat("Figure 4 written\n")
} else {
  cat("WARN: Figure 4 skipped (Stairway data missing)\n")
}

cat(sprintf("DONE 07.12a-smcpp_estimate.R — plots: %s\n", plot_dir))

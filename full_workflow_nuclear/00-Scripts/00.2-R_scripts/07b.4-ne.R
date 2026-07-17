#!/usr/bin/env Rscript
# 07b.4-ne.R
# Compiles Ne estimation results per genetic cluster (Pop_1/Pop_2/Pop_3)
# across purity thresholds (T70/T80/T90/T95).
# Tools: currentNe2 (default + -x) and GONE2 (default + -x).
# Called by 07b.4-ne.slurm with args: <out_base> <TAG|"ALL">

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

POP_COLORS <- c(Pop_1 = "#EE7600", Pop_2 = "#458B00", Pop_3 = "#CD2626")

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
out_base <- if (length(args) >= 1) args[1] else stop("out_base required")
tag      <- if (length(args) >= 2) args[2] else "ALL"

POPS      <- c("Pop_1", "Pop_2", "Pop_3")
THRESHOLDS <- c("T70", "T80", "T90", "T95")

cat(sprintf("INFO: out_base = %s\n", out_base))
cat(sprintf("INFO: tag      = %s\n", tag))

# -------------------------------------------------------------------
# Parsers (reused from 07.10-Ne_global.R)
# -------------------------------------------------------------------
parse_currentNe2 <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    cat(sprintf("  WARN: not found or empty: %s\n", path))
    return(NULL)
  }
  lines <- readLines(path, warn = FALSE)

  extract_next_num <- function(...) {
    for (lbl in c(...)) {
      idx <- grep(lbl, lines, fixed = TRUE)
      if (length(idx) > 0) {
        val <- suppressWarnings(as.numeric(trimws(lines[idx[1] + 1])))
        if (!is.na(val)) return(val)
      }
    }
    NA_real_
  }

  Ne      <- extract_next_num("# Ne point estimate:",
                               "# Ne of the entire metapopulation")
  CI50_lo <- extract_next_num("# Lower limit of 50% CI:", "# Lower bound of 50% CI:",
                               "# Lower bound of the 50% CI:",
                               "# Lower 50% limit of the Ne estimate:")
  CI50_hi <- extract_next_num("# Upper limit of 50% CI:", "# Upper bound of 50% CI:",
                               "# Upper bound of the 50% CI:",
                               "# Upper 50% limit of the Ne estimate:")
  CI90_lo <- extract_next_num("# Lower limit of 90% CI:", "# Lower bound of 90% CI:",
                               "# Lower bound of the 90% CI:",
                               "# Lower 90% limit of the Ne estimate:")
  CI90_hi <- extract_next_num("# Upper limit of 90% CI:", "# Upper bound of 90% CI:",
                               "# Upper bound of the 90% CI:",
                               "# Upper 90% limit of the Ne estimate:")
  NT      <- extract_next_num("# NT estimate:", "# NT:", "# N_T of the metapopulation")
  FST     <- extract_next_num("# FST estimate:", "# Estimated FST:", "# FST:", "# Fst (")
  m_rate  <- extract_next_num("# Migration rate estimate:", "# Estimated migration rate:",
                               "# Migration rate:")
  s_sub   <- extract_next_num("# Number of subpopulations estimate:",
                               "# Estimated number of subpopulations:",
                               "# Number of subpopulations:")

  data.frame(Ne = Ne, CI50_lower = CI50_lo, CI50_upper = CI50_hi,
             CI90_lower = CI90_lo, CI90_upper = CI90_hi,
             NT = NT, FST = FST, migration = m_rate, n_subpop = s_sub,
             stringsAsFactors = FALSE)
}

parse_GONE2 <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    cat(sprintf("  WARN: not found or empty: %s\n", path))
    return(NULL)
  }
  lines      <- readLines(path, warn = FALSE)
  lines      <- lines[nzchar(trimws(lines))]
  data_start <- which(grepl("^\\s*[0-9]", lines))[1]
  if (is.na(data_start)) return(NULL)
  d <- tryCatch(
    read.table(text = lines[data_start:length(lines)],
               header = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || ncol(d) < 2) return(NULL)
  data.frame(generation = as.integer(d[[1]]), Ne = as.numeric(d[[2]]))
}

# -------------------------------------------------------------------
# Collect results across thresholds and populations
# -------------------------------------------------------------------
collect_all <- function(tags_to_process) {
  rows_cne2  <- list()
  rows_gone2 <- list()

  for (tg in tags_to_process) {
    for (pop in POPS) {
      pop_dir  <- file.path(out_base, tg, pop)
      label    <- paste0(pop, "_", tg)

      # currentNe2 default
      f_def <- file.path(pop_dir, "currentNe2", paste0(label, "_default.txt"))
      r_def <- parse_currentNe2(f_def)
      if (!is.null(r_def)) {
        rows_cne2[[length(rows_cne2) + 1]] <- data.frame(
          threshold = tg, population = pop, mode = "default",
          r_def, stringsAsFactors = FALSE
        )
      }

      # currentNe2 -x
      f_mix <- file.path(pop_dir, "currentNe2", paste0(label, "_mix.txt"))
      r_mix <- parse_currentNe2(f_mix)
      if (!is.null(r_mix)) {
        rows_cne2[[length(rows_cne2) + 1]] <- data.frame(
          threshold = tg, population = pop, mode = "structure",
          r_mix, stringsAsFactors = FALSE
        )
      }

      # GONE2 default
      f_g_def <- file.path(pop_dir, "GONE2", "default",
                            paste0(label, ".ped_GONE2_Ne"))
      r_g_def <- parse_GONE2(f_g_def)
      if (!is.null(r_g_def)) {
        rows_gone2[[length(rows_gone2) + 1]] <- data.frame(
          threshold = tg, population = pop, mode = "default",
          r_g_def, stringsAsFactors = FALSE
        )
      }

      # GONE2 -x
      f_g_mix <- file.path(pop_dir, "GONE2", "structure",
                            paste0(label, ".ped_GONE2_Ne_mix"))
      r_g_mix <- parse_GONE2(f_g_mix)
      if (!is.null(r_g_mix)) {
        rows_gone2[[length(rows_gone2) + 1]] <- data.frame(
          threshold = tg, population = pop, mode = "structure",
          r_g_mix, stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    cne2  = if (length(rows_cne2)  > 0) do.call(rbind, rows_cne2)  else NULL,
    gone2 = if (length(rows_gone2) > 0) do.call(rbind, rows_gone2) else NULL
  )
}

# -------------------------------------------------------------------
# Plot helpers
# -------------------------------------------------------------------
plot_currentNe2_bar <- function(cne2_df, tg, plots_dir) {
  df <- cne2_df[cne2_df$threshold == tg & cne2_df$mode == "default", ]
  if (nrow(df) == 0 || all(is.na(df$Ne))) return(invisible(NULL))

  df$population <- factor(df$population, levels = POPS)

  p <- ggplot(df, aes(x = population, y = Ne, fill = population)) +
    geom_col(width = 0.55, show.legend = FALSE) +
    geom_errorbar(aes(ymin = CI90_lower, ymax = CI90_upper),
                  width = 0.2, linewidth = 0.7, colour = "grey30", na.rm = TRUE) +
    geom_errorbar(aes(ymin = CI50_lower, ymax = CI50_upper),
                  width = 0.12, linewidth = 1.3, colour = "grey15", na.rm = TRUE) +
    scale_fill_manual(values = POP_COLORS) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = sprintf("currentNe2 — current Ne per cluster (%s)", tg),
      subtitle = "Bars: 50% CI (thick) and 90% CI (thin)",
      x        = NULL, y = "Effective population size (Ne)"
    )

  ggsave(file.path(plots_dir, sprintf("currentNe2_default_%s.png", tg)),
         p, width = 6, height = 5, dpi = 200)
}

plot_GONE2_trajectories <- function(gone2_df, tg, plots_dir, max_gen = 500) {
  df <- gone2_df[gone2_df$threshold == tg & gone2_df$mode == "default", ]
  df <- df[df$generation <= max_gen & df$Ne > 0, ]
  if (nrow(df) == 0) return(invisible(NULL))

  df$population <- factor(df$population, levels = POPS)

  p <- ggplot(df, aes(x = generation, y = Ne, colour = population)) +
    geom_line(linewidth = 0.9) +
    scale_x_reverse(labels = comma) +
    scale_y_log10(labels = comma) +
    scale_colour_manual(values = POP_COLORS, name = "Cluster") +
    labs(
      title    = sprintf("GONE2 — temporal Ne trajectory per cluster (%s)", tg),
      subtitle = "Panmictic assumption; x-axis reversed (recent on right)",
      x        = "Generations ago", y = "Ne (log scale)"
    ) +
    theme(legend.position = "right")

  ggsave(file.path(plots_dir, sprintf("GONE2_default_trajectories_%s.png", tg)),
         p, width = 9, height = 5, dpi = 200)

  # default vs structure-corrected per population
  df_all <- gone2_df[gone2_df$threshold == tg, ]
  df_all <- df_all[df_all$generation <= max_gen & df_all$Ne > 0, ]
  if (length(unique(df_all$mode)) > 1) {
    df_all$population <- factor(df_all$population, levels = POPS)
    p2 <- ggplot(df_all, aes(x = generation, y = Ne,
                              colour = population, linetype = mode)) +
      geom_line(linewidth = 0.8) +
      scale_x_reverse(labels = comma) +
      scale_y_log10(labels = comma) +
      scale_colour_manual(values = POP_COLORS, name = "Cluster") +
      scale_linetype_manual(
        values = c("default" = "solid", "structure" = "dashed"),
        name   = "Mode"
      ) +
      labs(
        title    = sprintf("GONE2 — default vs structure-corrected (%s)", tg),
        subtitle = "Dashed: structure-corrected (-x); solid: panmictic",
        x        = "Generations ago", y = "Ne (log scale)"
      ) +
      theme(legend.position = "right")

    ggsave(file.path(plots_dir, sprintf("GONE2_comparison_%s.png", tg)),
           p2, width = 9, height = 5, dpi = 200)
  }
}

plot_currentNe2_across_thresholds <- function(cne2_df, plots_dir) {
  df <- cne2_df[cne2_df$mode == "default" & !is.na(cne2_df$Ne), ]
  if (nrow(df) == 0) return(invisible(NULL))

  df$threshold  <- factor(df$threshold, levels = THRESHOLDS)
  df$population <- factor(df$population, levels = POPS)

  p <- ggplot(df, aes(x = threshold, y = Ne, colour = population, group = population)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = CI90_lower, ymax = CI90_upper),
                  width = 0.15, linewidth = 0.6, na.rm = TRUE) +
    scale_colour_manual(values = POP_COLORS, name = "Cluster") +
    scale_y_continuous(labels = comma) +
    labs(
      title    = "currentNe2 — current Ne across purity thresholds",
      subtitle = "Error bars: 90% CI; default (panmictic) mode",
      x        = "Purity threshold", y = "Effective population size (Ne)"
    ) +
    theme(legend.position = "right")

  ggsave(file.path(plots_dir, "currentNe2_Ne_across_thresholds.png"),
         p, width = 8, height = 5, dpi = 200)
}

plot_GONE2_across_thresholds <- function(gone2_df, plots_dir, max_gen = 200) {
  df <- gone2_df[gone2_df$mode == "default" & gone2_df$generation <= max_gen & gone2_df$Ne > 0, ]
  if (nrow(df) == 0) return(invisible(NULL))

  df$threshold  <- factor(df$threshold, levels = THRESHOLDS)
  df$population <- factor(df$population, levels = POPS)

  p <- ggplot(df, aes(x = generation, y = Ne,
                      colour = population, linetype = threshold)) +
    geom_line(linewidth = 0.7) +
    scale_x_reverse(labels = comma) +
    scale_y_log10(labels = comma) +
    scale_colour_manual(values = POP_COLORS, name = "Cluster") +
    scale_linetype_manual(
      values = c(T70 = "solid", T80 = "dashed", T90 = "dotdash", T95 = "dotted"),
      name   = "Threshold"
    ) +
    labs(
      title    = "GONE2 — temporal Ne trajectories across purity thresholds",
      subtitle = "Panmictic assumption; last 200 generations",
      x        = "Generations ago", y = "Ne (log scale)"
    ) +
    facet_wrap(~ population, ncol = 3) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

  ggsave(file.path(plots_dir, "GONE2_trajectories_across_thresholds.png"),
         p, width = 12, height = 5, dpi = 200)
}

# -------------------------------------------------------------------
# Main: per-threshold run
# -------------------------------------------------------------------
if (tag != "ALL") {
  tags_to_run <- tag
  out_dir     <- file.path(out_base, tag)
  plots_dir   <- file.path(out_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  res <- collect_all(tags_to_run)

  if (!is.null(res$cne2)) {
    plot_currentNe2_bar(res$cne2, tag, plots_dir)
    cat(sprintf("INFO: currentNe2 bar plot saved for %s\n", tag))

    # Per-threshold summary table
    write.table(res$cne2[res$cne2$threshold == tag, ],
                file.path(out_dir, sprintf("07b.4-Ne_summary_%s.tsv", tag)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if (!is.null(res$gone2)) {
    plot_GONE2_trajectories(res$gone2, tag, plots_dir)
    cat(sprintf("INFO: GONE2 trajectory plots saved for %s\n", tag))
  }

  cat(sprintf("DONE: 07b.4-ne.R %s\n", tag))

} else {
  # -------------------------------------------------------------------
  # Cross-threshold summary (tag == "ALL")
  # -------------------------------------------------------------------
  plots_dir_all <- file.path(out_base, "plots")
  dir.create(plots_dir_all, recursive = TRUE, showWarnings = FALSE)

  res <- collect_all(THRESHOLDS)

  if (!is.null(res$cne2)) {
    write.table(res$cne2,
                file.path(out_base, "07b.4-Ne_summary_all.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("INFO: combined Ne summary table written\n")
    plot_currentNe2_across_thresholds(res$cne2, plots_dir_all)
  }

  if (!is.null(res$gone2)) {
    plot_GONE2_across_thresholds(res$gone2, plots_dir_all)
  }

  cat("DONE: 07b.4-ne.R ALL thresholds\n")
  cat(sprintf("Plots : %s/\n", plots_dir_all))
  cat(sprintf("Table : %s/07b.4-Ne_summary_all.tsv\n", out_base))
}

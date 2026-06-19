#!/usr/bin/env Rscript
# 07.10-Ne_global.R
# Compiles Ne estimation results from the pooled French Guiana dataset:
#   - currentNe2 (default): current Ne, panmictic assumption, 50%/90% CI
#   - currentNe2 (-x)     : metapopulation Ne + FST + migration rate + s
#   - GONE2 (default)     : temporal Ne trajectory, panmictic assumption
#   - GONE2 (-x)          : structure-corrected temporal Ne trajectory
# Produces summary tables and plots; called by 07.10-Ne_global.slurm.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else stop("out_dir required")
label   <- if (length(args) >= 2) args[2] else "guiana_all"

plots_dir <- file.path(out_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

cne2_dir <- file.path(out_dir, "currentNe2")
gone2_dir <- file.path(out_dir, "GONE2")

cat(sprintf("INFO: out_dir = %s\n", out_dir))
cat(sprintf("INFO: label   = %s\n", label))

# -------------------------------------------------------------------
# Parser: currentNe2 (default and -x/mix output)
# Tries multiple label variants to be robust to minor version changes.
# -------------------------------------------------------------------
parse_currentNe2 <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    cat(sprintf("WARN: not found or empty: %s\n", path))
    return(NULL)
  }
  lines <- readLines(path, warn = FALSE)

  extract_next_num <- function(...) {
    for (label in c(...)) {
      idx <- grep(label, lines, fixed = TRUE)
      if (length(idx) > 0) {
        val <- suppressWarnings(as.numeric(trimws(lines[idx[1] + 1])))
        if (!is.na(val)) return(val)
      }
    }
    NA_real_
  }

  ne_idx <- grep("# Ne point estimate:", lines, fixed = TRUE)
  Ne <- if (length(ne_idx) > 0)
    suppressWarnings(as.numeric(trimws(lines[ne_idx[1] + 1]))) else NA_real_

  CI50_lo <- extract_next_num(
    "# Lower limit of 50% CI:", "# Lower bound of 50% CI:",
    "# Lower bound of the 50% CI:")
  CI50_hi <- extract_next_num(
    "# Upper limit of 50% CI:", "# Upper bound of 50% CI:",
    "# Upper bound of the 50% CI:")
  CI90_lo <- extract_next_num(
    "# Lower limit of 90% CI:", "# Lower bound of 90% CI:",
    "# Lower bound of the 90% CI:")
  CI90_hi <- extract_next_num(
    "# Upper limit of 90% CI:", "# Upper bound of 90% CI:",
    "# Upper bound of the 90% CI:")

  NT  <- extract_next_num("# NT estimate:", "# NT:")
  FST <- extract_next_num("# FST estimate:", "# Estimated FST:", "# FST:")
  m   <- extract_next_num(
    "# Migration rate estimate:", "# Estimated migration rate:",
    "# Migration rate:")
  s   <- extract_next_num(
    "# Number of subpopulations estimate:",
    "# Estimated number of subpopulations:",
    "# Number of subpopulations:")

  data.frame(
    Ne = Ne, CI50_lower = CI50_lo, CI50_upper = CI50_hi,
    CI90_lower = CI90_lo, CI90_upper = CI90_hi,
    NT = NT, FST = FST, migration = m, n_subpop = s,
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------
# Parser: GONE2 trajectory (tab-separated, generation + Ne)
# -------------------------------------------------------------------
parse_GONE2 <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    cat(sprintf("WARN: not found or empty: %s\n", path))
    return(NULL)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  data_start <- which(grepl("^\\s*[0-9]", lines))[1]
  if (is.na(data_start)) {
    cat(sprintf("WARN: no numeric data in %s\n", path))
    return(NULL)
  }
  d <- tryCatch(
    read.table(text = lines[data_start:length(lines)],
               header = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || ncol(d) < 2) return(NULL)
  data.frame(generation = as.integer(d[[1]]),
             Ne         = as.numeric(d[[2]]),
             stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------
# STEP 1: Parse all outputs
# -------------------------------------------------------------------
cat("INFO: parsing currentNe2 default output...\n")
cne2_def <- parse_currentNe2(
  file.path(cne2_dir, paste0(label, "_currentNe2_default.txt"))
)

cat("INFO: parsing currentNe2 -x output...\n")
cne2_mix <- parse_currentNe2(
  file.path(cne2_dir, paste0(label, "_currentNe2_mix.txt"))
)

cat("INFO: parsing GONE2 default trajectory...\n")
gone2_def <- parse_GONE2(
  file.path(gone2_dir, "default", paste0(label, ".ped_GONE2_Ne"))
)

cat("INFO: parsing GONE2 -x trajectory...\n")
gone2_str <- parse_GONE2(
  file.path(gone2_dir, "structure", paste0(label, ".ped_GONE2_Ne"))
)

# -------------------------------------------------------------------
# STEP 2: Ne point estimate summary table
# -------------------------------------------------------------------
cat("INFO: building summary table...\n")

summary_rows <- list()

if (!is.null(cne2_def) && !is.na(cne2_def$Ne)) {
  summary_rows[["currentNe2_default"]] <- data.frame(
    tool        = "currentNe2",
    mode        = "default (panmictic)",
    Ne_estimate = round(cne2_def$Ne),
    CI50_lower  = round(cne2_def$CI50_lower),
    CI50_upper  = round(cne2_def$CI50_upper),
    CI90_lower  = round(cne2_def$CI90_lower),
    CI90_upper  = round(cne2_def$CI90_upper),
    stringsAsFactors = FALSE
  )
}

if (!is.null(cne2_mix) && !is.na(cne2_mix$Ne)) {
  summary_rows[["currentNe2_struct"]] <- data.frame(
    tool        = "currentNe2",
    mode        = "-x (metapopulation Ne)",
    Ne_estimate = round(cne2_mix$Ne),
    CI50_lower  = NA_integer_,
    CI50_upper  = NA_integer_,
    CI90_lower  = NA_integer_,
    CI90_upper  = NA_integer_,
    stringsAsFactors = FALSE
  )
}

if (!is.null(gone2_def) && nrow(gone2_def) > 0) {
  ne_gen1 <- gone2_def$Ne[gone2_def$generation == 1]
  if (length(ne_gen1) == 0) ne_gen1 <- gone2_def$Ne[1]
  summary_rows[["GONE2_default"]] <- data.frame(
    tool        = "GONE2",
    mode        = "default (panmictic, gen 1)",
    Ne_estimate = round(ne_gen1),
    CI50_lower  = NA_integer_,
    CI50_upper  = NA_integer_,
    CI90_lower  = NA_integer_,
    CI90_upper  = NA_integer_,
    stringsAsFactors = FALSE
  )
}

if (!is.null(gone2_str) && nrow(gone2_str) > 0) {
  ne_gen1_x <- gone2_str$Ne[gone2_str$generation == 1]
  if (length(ne_gen1_x) == 0) ne_gen1_x <- gone2_str$Ne[1]
  summary_rows[["GONE2_struct"]] <- data.frame(
    tool        = "GONE2",
    mode        = "-x (structure-corrected, gen 1)",
    Ne_estimate = round(ne_gen1_x),
    CI50_lower  = NA_integer_,
    CI50_upper  = NA_integer_,
    CI90_lower  = NA_integer_,
    CI90_upper  = NA_integer_,
    stringsAsFactors = FALSE
  )
}

if (length(summary_rows) > 0) {
  ne_summary <- do.call(rbind, summary_rows)
  rownames(ne_summary) <- NULL
  write.table(ne_summary,
              file.path(out_dir, "07.10-Ne_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("INFO: Ne summary table written\n")
  print(ne_summary[, c("tool", "mode", "Ne_estimate", "CI90_lower", "CI90_upper")])
}

# -------------------------------------------------------------------
# STEP 3: Structure parameters table (currentNe2 -x)
# -------------------------------------------------------------------
if (!is.null(cne2_mix)) {
  struct_df <- data.frame(
    parameter = c("Ne_meta", "NT", "FST", "migration_rate", "n_subpop"),
    value     = c(
      cne2_mix$Ne, cne2_mix$NT, cne2_mix$FST,
      cne2_mix$migration, cne2_mix$n_subpop
    ),
    description = c(
      "Effective size per subpopulation",
      "Total effective size (metapopulation)",
      "Average pairwise FST between subpopulations",
      "Migration rate (proportion per generation)",
      "Estimated number of subpopulations"
    ),
    stringsAsFactors = FALSE
  )
  write.table(struct_df,
              file.path(out_dir, "07.10-structure_parameters.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("INFO: structure parameters table written\n")
  cat(sprintf(
    "  Ne_meta  = %.1f\n  NT       = %.1f\n  FST      = %.4f\n  m        = %.5f\n  s        = %.2f\n",
    ifelse(is.na(cne2_mix$Ne), NA, cne2_mix$Ne),
    ifelse(is.na(cne2_mix$NT), NA, cne2_mix$NT),
    ifelse(is.na(cne2_mix$FST), NA, cne2_mix$FST),
    ifelse(is.na(cne2_mix$migration), NA, cne2_mix$migration),
    ifelse(is.na(cne2_mix$n_subpop), NA, cne2_mix$n_subpop)
  ))
}

# -------------------------------------------------------------------
# STEP 4: Plot — GONE2 temporal trajectories (default vs -x)
# -------------------------------------------------------------------
if (!is.null(gone2_def) || !is.null(gone2_str)) {
  cat("INFO: plotting GONE2 trajectories...\n")

  traj_data <- rbind(
    if (!is.null(gone2_def)) {
      cbind(gone2_def, run = "default (panmictic)")
    },
    if (!is.null(gone2_str)) {
      cbind(gone2_str, run = "-x (structure-corrected)")
    }
  )

  max_gen <- 500
  traj_data <- traj_data[traj_data$generation <= max_gen & traj_data$Ne > 0, ]

  p_traj <- ggplot(traj_data, aes(x = generation, y = Ne,
                                   colour = run, linetype = run)) +
    geom_line(linewidth = 0.8) +
    scale_x_reverse(labels = comma) +
    scale_y_log10(labels = comma) +
    scale_colour_manual(values = c(
      "default (panmictic)"      = "#2166AC",
      "-x (structure-corrected)" = "#D6604D"
    )) +
    scale_linetype_manual(values = c(
      "default (panmictic)"      = "solid",
      "-x (structure-corrected)" = "dashed"
    )) +
    labs(
      title    = "GONE2 — temporal Ne trajectory, all French Guiana individuals",
      subtitle = "default: panmictic assumption; -x: structure-corrected",
      x        = "Generations ago",
      y        = "Effective population size (Ne, log scale)",
      colour   = "Mode", linetype = "Mode"
    ) +
    theme(
      legend.position  = "bottom",
      plot.title       = element_text(size = 12, face = "bold"),
      plot.subtitle    = element_text(size = 10),
      axis.title       = element_text(size = 10)
    )

  ggsave(file.path(plots_dir, "07.10-GONE2_trajectories.png"),
         p_traj, width = 10, height = 6, dpi = 200)
  cat("INFO: GONE2 trajectory plot saved\n")
}

# -------------------------------------------------------------------
# STEP 5: Plot — currentNe2 Ne comparison (default vs -x, with CI)
# -------------------------------------------------------------------
if (!is.null(cne2_def) || !is.null(cne2_mix)) {
  cat("INFO: plotting currentNe2 comparison...\n")

  cne2_plot_df <- rbind(
    if (!is.null(cne2_def) && !is.na(cne2_def$Ne)) {
      data.frame(
        mode       = "default\n(panmictic)",
        Ne         = cne2_def$Ne,
        CI90_lower = cne2_def$CI90_lower,
        CI90_upper = cne2_def$CI90_upper,
        CI50_lower = cne2_def$CI50_lower,
        CI50_upper = cne2_def$CI50_upper,
        stringsAsFactors = FALSE
      )
    },
    if (!is.null(cne2_mix) && !is.na(cne2_mix$Ne)) {
      data.frame(
        mode       = "-x\n(Ne per subpopulation)",
        Ne         = cne2_mix$Ne,
        CI90_lower = NA_real_,
        CI90_upper = NA_real_,
        CI50_lower = NA_real_,
        CI50_upper = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )

  if (!is.null(cne2_plot_df) && nrow(cne2_plot_df) > 0) {
    p_cne2 <- ggplot(cne2_plot_df, aes(x = mode, y = Ne, fill = mode)) +
      geom_col(width = 0.5, show.legend = FALSE) +
      geom_errorbar(aes(ymin = CI90_lower, ymax = CI90_upper),
                    width = 0.15, linewidth = 0.7, colour = "grey30",
                    na.rm = TRUE) +
      geom_errorbar(aes(ymin = CI50_lower, ymax = CI50_upper),
                    width = 0.1, linewidth = 1.2, colour = "grey20",
                    na.rm = TRUE) +
      scale_fill_manual(values = c(
        "default\n(panmictic)"            = "#4393C3",
        "-x\n(Ne per subpopulation)"      = "#D6604D"
      )) +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
      labs(
        title    = "currentNe2 — current Ne estimate, all French Guiana individuals",
        subtitle = "Bars show 50% CI (thick) and 90% CI (thin) where available",
        x        = "Mode",
        y        = "Effective population size (Ne)"
      ) +
      theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title    = element_text(size = 10)
      )

    ggsave(file.path(plots_dir, "07.10-currentNe2_comparison.png"),
           p_cne2, width = 7, height = 6, dpi = 200)
    cat("INFO: currentNe2 comparison plot saved\n")
  }
}

# -------------------------------------------------------------------
# STEP 6: Plot — structure parameters barplot (currentNe2 -x)
# Plots FST, migration rate, and inferred s alongside ADMIXTURE K.
# -------------------------------------------------------------------
if (!is.null(cne2_mix) && any(!is.na(c(cne2_mix$FST, cne2_mix$migration, cne2_mix$n_subpop)))) {
  cat("INFO: plotting structure parameters...\n")

  param_df <- data.frame(
    parameter = c("FST", "Migration rate (m)", "Inferred s"),
    value     = c(
      ifelse(is.na(cne2_mix$FST), NA, cne2_mix$FST),
      ifelse(is.na(cne2_mix$migration), NA, cne2_mix$migration),
      ifelse(is.na(cne2_mix$n_subpop), NA, cne2_mix$n_subpop)
    ),
    stringsAsFactors = FALSE
  )
  param_df <- param_df[!is.na(param_df$value), ]

  if (nrow(param_df) > 0) {
    p_struct <- ggplot(param_df, aes(x = parameter, y = value, fill = parameter)) +
      geom_col(width = 0.5, show.legend = FALSE) +
      geom_text(aes(label = signif(value, 4)),
                vjust = -0.5, size = 3.5) +
      scale_fill_brewer(palette = "Set2") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title    = "currentNe2 -x — inferred population structure parameters",
        subtitle = "All 155 French Guiana individuals; s = inferred number of subpopulations",
        x        = NULL,
        y        = "Estimated value"
      ) +
      theme(
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.title.y  = element_text(size = 10)
      )

    ggsave(file.path(plots_dir, "07.10-structure_parameters.png"),
           p_struct, width = 7, height = 5, dpi = 200)
    cat("INFO: structure parameters plot saved\n")
  }
}

cat("DONE: 07.10-Ne_global.R completed\n")
cat(sprintf("Plots  : %s/\n", plots_dir))
cat(sprintf("Tables : %s/07.10-*.tsv\n", out_dir))

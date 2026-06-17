#!/usr/bin/env Rscript
# 07.8-Ne_cluster.R
# Compiles Ne estimation results per genetic cluster:
#   - currentNe : current Ne with 50% / 90% bootstrap CI
#   - GONE      : temporal Ne trajectory over ~2000 generations
# Produces summary table + plots for both estimators.

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
args      <- commandArgs(trailingOnly = TRUE)
out_dir   <- if (length(args) >= 1) args[1] else stop("out_dir required")
pops_file <- if (length(args) >= 2) args[2] else stop("pops_file required")
k_global  <- if (length(args) >= 3) args[3] else "?"

pops     <- readLines(pops_file)
pops     <- pops[nzchar(pops)]
pop_dir  <- file.path(out_dir, "per_cluster")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Discrete color palette for clusters (up to 4 populations)
POP_COLORS <- c(
  Pop_1 = "#EE7600",
  Pop_2 = "#458B00",
  Pop_3 = "#CD2626",
  Pop_4 = "#9A32CD"
)

# -------------------------------------------------------------------
# Helper: parse currentNe output file
# Expected format (lines containing key values):
#   Ne = XXXX
#   50% CI : XXXX - XXXX
#   90% CI : XXXX - XXXX
# Returns a one-row data frame or NULL on failure.
# -------------------------------------------------------------------
parse_currentNe <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  lines <- readLines(path, warn = FALSE)

  # currentNe format: label on one line (e.g. "# Ne point estimate:"),
  # numeric value on the next line. idx[1] = first match = genome-wide estimate.
  extract_next_num <- function(label) {
    idx <- grep(label, lines, fixed = TRUE)
    if (length(idx) == 0) return(NA_real_)
    suppressWarnings(as.numeric(trimws(lines[idx[1] + 1])))
  }

  data.frame(
    Ne         = extract_next_num("# Ne point estimate:"),
    CI50_lower = extract_next_num("# Lower bound of the 50% CI:"),
    CI50_upper = extract_next_num("# Upper bound of the 50% CI:"),
    CI90_lower = extract_next_num("# Lower bound of the 90% CI:"),
    CI90_upper = extract_next_num("# Upper bound of the 90% CI:"),
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------
# Helper: parse GONE Output_Ne_<pop> file
# Expected format (tab/space-separated, header on first data line):
#   Generation  Geometric_mean  ...
#   1           XXXX
#   ...
# Returns a data frame with columns: generation, Ne
# -------------------------------------------------------------------
parse_GONE <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]

  # Find first all-numeric line (skip header / comment lines)
  data_start <- which(grepl("^\\s*[0-9]", lines))[1]
  if (is.na(data_start)) {
    cat(sprintf("WARN: no numeric data found in %s\n", path))
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
# STEP 1: Collect results for all populations
# -------------------------------------------------------------------
cat("INFO: parsing currentNe outputs...\n")
cne_list <- lapply(pops, function(pop) {
  f <- file.path(pop_dir, pop, "currentNe", paste0(pop, "_currentNe.txt"))
  # currentNe may also produce the file without extension; try both
  if (!file.exists(f)) f <- file.path(pop_dir, pop, "currentNe",
                                       paste0(pop, "_currentNe"))
  d <- parse_currentNe(f)
  if (!is.null(d)) {
    d$population <- pop
    d$n_individuals <- length(readLines(
      file.path(pop_dir, paste0(pop, ".samples.txt")), warn = FALSE))
  }
  d
})
names(cne_list) <- pops
cne_valid <- Filter(Negate(is.null), cne_list)

cat("INFO: parsing GONE outputs...\n")
gone_list <- lapply(pops, function(pop) {
  f <- file.path(pop_dir, pop, "GONE", paste0("Output_Ne_", pop))
  d <- parse_GONE(f)
  if (!is.null(d)) d$population <- pop
  d
})
names(gone_list) <- pops
gone_valid <- Filter(Negate(is.null), gone_list)

# -------------------------------------------------------------------
# STEP 2: currentNe summary table
# -------------------------------------------------------------------
if (length(cne_valid) > 0) {
  cne_df <- do.call(rbind, cne_valid)
  cne_df <- cne_df[, c("population", "n_individuals",
                        "Ne", "CI50_lower", "CI50_upper",
                        "CI90_lower", "CI90_upper")]
  cne_df[sapply(cne_df, is.numeric)] <- lapply(
    cne_df[sapply(cne_df, is.numeric)], round, 0)
  cne_tsv <- file.path(out_dir, sprintf("07.8-currentNe_summary_K%s.tsv", k_global))
  write.table(cne_df, cne_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("currentNe summary written:", cne_tsv, "\n")
  print(cne_df)
} else {
  cat("WARN: no valid currentNe outputs found\n")
  cne_df <- NULL
}

# -------------------------------------------------------------------
# STEP 3: GONE summary table (Ne at generation 1 = most recent)
# -------------------------------------------------------------------
if (length(gone_valid) > 0) {
  gone_df <- do.call(rbind, gone_valid)
  gone_tsv <- file.path(out_dir, sprintf("07.8-GONE_Ne_trajectory_K%s.tsv", k_global))
  write.table(gone_df, gone_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("GONE trajectory written:", gone_tsv, "\n")

  # Most recent Ne per population (generation 1)
  gone_current <- do.call(rbind, lapply(gone_valid, function(d) {
    d_sorted <- d[order(d$generation), ]
    data.frame(population = d_sorted$population[1],
               Ne_gen1    = d_sorted$Ne[1],
               stringsAsFactors = FALSE)
  }))
} else {
  cat("WARN: no valid GONE outputs found\n")
  gone_df    <- NULL
  gone_current <- NULL
}

# -------------------------------------------------------------------
# STEP 4: Plot — currentNe barplot with CI
# -------------------------------------------------------------------
if (!is.null(cne_df) && any(!is.na(cne_df$Ne))) {
  df_p <- cne_df[!is.na(cne_df$Ne), ]
  df_p$population <- factor(df_p$population, levels = sort(unique(df_p$population)))

  p_cne <- ggplot(df_p, aes(x = population, y = Ne, fill = population)) +
    geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
    geom_errorbar(aes(ymin = CI90_lower, ymax = CI90_upper),
                  width = 0.18, color = "grey20", linewidth = 0.7) +
    geom_errorbar(aes(ymin = CI50_lower, ymax = CI50_upper),
                  width = 0.10, color = "black", linewidth = 1.1) +
    geom_text(aes(label = format(round(Ne, 0), big.mark = ",")),
              vjust = -0.5, size = 3.2) +
    scale_fill_manual(values = POP_COLORS, na.value = "grey50") +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
    labs(title    = sprintf("Current effective population size — K=%s", k_global),
         subtitle = "Error bars: thick = 50% CI, thin = 90% CI (neural-network bootstrap)",
         x = NULL, y = "Ne (currentNe)") +
    theme(legend.position = "none",
          axis.text.x     = element_text(size = 10),
          plot.subtitle   = element_text(size = 8, color = "grey40"))

  fn_cne <- file.path(plot_dir, sprintf("currentNe_K%s.png", k_global))
  ggsave(fn_cne, p_cne, width = max(6, length(unique(df_p$population)) * 1.8),
         height = 5, dpi = 300)
  cat("Plot written:", fn_cne, "\n")
}

# -------------------------------------------------------------------
# STEP 5: Plot — GONE temporal Ne trajectory
# -------------------------------------------------------------------
if (!is.null(gone_df) && nrow(gone_df) > 0) {
  gone_df$population <- factor(gone_df$population,
                                levels = sort(unique(gone_df$population)))
  # Restrict to meaningful generations (1–500, beyond that LD signal is weak)
  gone_plot <- gone_df[gone_df$generation <= 500 & !is.na(gone_df$Ne) & gone_df$Ne > 0, ]

  if (nrow(gone_plot) > 0) {
    p_gone <- ggplot(gone_plot, aes(x = generation, y = Ne,
                                    color = population, group = population)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = POP_COLORS, na.value = "grey50") +
      scale_x_continuous(breaks = c(1, 10, 50, 100, 200, 500)) +
      scale_y_log10(labels = scales::comma) +
      labs(title    = sprintf("Temporal Ne trajectory — GONE — K=%s", k_global),
           subtitle = "Generations 1–500 — log10(Ne) scale",
           x = "Generations ago", y = "Ne (GONE, log scale)",
           color = "Cluster") +
      theme(plot.subtitle = element_text(size = 8, color = "grey40"))

    fn_gone <- file.path(plot_dir, sprintf("GONE_trajectory_K%s.png", k_global))
    ggsave(fn_gone, p_gone, width = 8, height = 5, dpi = 300)
    cat("Plot written:", fn_gone, "\n")

    # Full trajectory (all generations)
    gone_all <- gone_df[!is.na(gone_df$Ne) & gone_df$Ne > 0, ]
    if (max(gone_all$generation) > 500) {
      p_gone_full <- p_gone %+% gone_all +
        scale_x_continuous() +
        labs(subtitle = "All generations — log10(Ne) scale")
      fn_full <- file.path(plot_dir, sprintf("GONE_trajectory_full_K%s.png", k_global))
      ggsave(fn_full, p_gone_full, width = 9, height = 5, dpi = 300)
      cat("Plot written:", fn_full, "\n")
    }
  }
}

# -------------------------------------------------------------------
# STEP 6: Combined comparison plot (current Ne from both tools)
# -------------------------------------------------------------------
if (!is.null(cne_df) && !is.null(gone_current)) {
  cne_current <- cne_df[!is.na(cne_df$Ne),
                         c("population", "Ne", "CI90_lower", "CI90_upper")]
  colnames(cne_current)[2:4] <- c("Ne_currentNe", "CI90_lo", "CI90_hi")

  comp <- merge(cne_current, gone_current, by = "population", all = TRUE)

  if (nrow(comp) > 0) {
    comp_long <- rbind(
      data.frame(population = comp$population,
                 Ne = comp$Ne_currentNe, tool = "currentNe",
                 stringsAsFactors = FALSE),
      data.frame(population = comp$population,
                 Ne = comp$Ne_gen1, tool = "GONE",
                 stringsAsFactors = FALSE)
    )
    comp_long <- comp_long[!is.na(comp_long$Ne), ]
    comp_long$population <- factor(comp_long$population,
                                   levels = sort(unique(comp_long$population)))

    p_comp <- ggplot(comp_long,
                     aes(x = population, y = Ne, fill = population,
                         alpha = tool)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6,
               color = "grey30", linewidth = 0.3) +
      scale_fill_manual(values = POP_COLORS, na.value = "grey50") +
      scale_alpha_manual(values = c(currentNe = 1.0, GONE = 0.55)) +
      scale_y_continuous(labels = scales::comma,
                         expand = expansion(mult = c(0, 0.12))) +
      labs(title    = sprintf("Current Ne — currentNe vs GONE (generation 1) — K=%s", k_global),
           x = NULL, y = "Ne", fill = "Cluster", alpha = "Tool") +
      theme(axis.text.x = element_text(size = 10))

    fn_comp <- file.path(plot_dir, sprintf("Ne_comparison_K%s.png", k_global))
    ggsave(fn_comp, p_comp,
           width = max(6, length(unique(comp_long$population)) * 2.5),
           height = 5, dpi = 300)
    cat("Plot written:", fn_comp, "\n")
  }
}

cat(sprintf("DONE 07.8 Ne R script completed — K=%s\n", k_global))

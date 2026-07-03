#!/usr/bin/env Rscript
# 07.8-Ne_cluster.R
# Compiles Ne estimation results per genetic cluster for four tools:
#   - currentNe  v1: current Ne  (neural-network LD, 50%/90% CI)
#   - currentNe2 v2: current Ne  (improved LD, default + structure -x)
#   - GONE   v1    : temporal Ne trajectory (~2000 generations)
#   - GONE2  v2    : temporal Ne trajectory (default + structure -x)
# Produces summary tables + plots organised in tool-specific subdirectories.

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

cne_out_dir  <- file.path(out_dir, "currentNe")
cne2_out_dir <- file.path(out_dir, "currentNe2")
g1_out_dir   <- file.path(out_dir, "GONE_1")
g2_out_dir   <- file.path(out_dir, "GONE_2")
for (d in c(cne_out_dir, cne2_out_dir, g1_out_dir, g2_out_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

POP_COLORS <- c(
  Pop_1 = "#EE7600",
  Pop_2 = "#458B00",
  Pop_3 = "#CD2626",
  Pop_4 = "#9A32CD"
)

# -------------------------------------------------------------------
# Helper: count individuals per population
# -------------------------------------------------------------------
n_ind_for <- function(pop) {
  f <- file.path(pop_dir, paste0(pop, ".samples.txt"))
  if (file.exists(f)) length(readLines(f, warn = FALSE)) else NA_integer_
}

# -------------------------------------------------------------------
# Parser: currentNe v1
# Label on one line, numeric value on the next line.
# idx[1] = first match = genome-wide estimate (two estimates present).
# -------------------------------------------------------------------
parse_currentNe <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  lines <- readLines(path, warn = FALSE)
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
# Parser: currentNe2 (default output and -x/mix output)
# Labels differ slightly from v1; mix output additionally contains
# structure parameters (FST, migration rate, subpopulation number).
# Tries multiple label variants to be robust to version differences.
# -------------------------------------------------------------------
parse_currentNe2 <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
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

  # Ne — panmictic (default) or metapopulation estimate (-x)
  # Default output: "# Ne point estimate:"
  # Mix (-x) output: "# Ne of the entire metapopulation (effective..."
  Ne <- extract_next_num(
    "# Ne point estimate:",
    "# Ne of the entire metapopulation")

  # CI — labels differ between default and -x output
  CI50_lo <- extract_next_num(
    "# Lower limit of 50% CI:",
    "# Lower bound of 50% CI:",
    "# Lower bound of the 50% CI:",
    "# Lower 50% limit of the Ne estimate:")
  CI50_hi <- extract_next_num(
    "# Upper bound of 50% CI:",
    "# Upper limit of 50% CI:",
    "# Upper bound of the 50% CI:",
    "# Upper 50% limit of the Ne estimate:")
  CI90_lo <- extract_next_num(
    "# Lower limit of 90% CI:",
    "# Lower bound of 90% CI:",
    "# Lower bound of the 90% CI:",
    "# Lower 90% limit of the Ne estimate:")
  CI90_hi <- extract_next_num(
    "# Upper limit of 90% CI:",
    "# Upper bound of 90% CI:",
    "# Upper bound of the 90% CI:",
    "# Upper 90% limit of the Ne estimate:")

  # Structure parameters (present only in mix/-x output)
  # -x labels: "# N_T of the metapopulation...:", "# Fst (subpopulation...):", "# Migration rate:"
  NT  <- extract_next_num("# NT estimate:", "# NT:", "# N_T of the metapopulation")
  FST <- extract_next_num("# FST estimate:", "# Estimated FST:", "# FST:", "# Fst (")
  m   <- extract_next_num(
    "# Migration rate estimate:",
    "# Estimated migration rate:",
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
# Parser: GONE v1 and GONE2 default (2-column tab-separated format)
# Skips comment lines (#), reads generation (col1) + Ne (col2).
# -------------------------------------------------------------------
parse_GONE_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
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
# Parser: GONE2 -x structure output (5-column format)
# Header: Rec_rate_bin / generation / N_T_metapop / Ne_metapop / d²
# Returns generation (col 2) + Ne_metapop (col 4).
# -------------------------------------------------------------------
parse_GONE2_mix <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  lines <- readLines(path, warn = FALSE)
  data_start <- which(grepl("Rec_rate_bin", lines, fixed = TRUE))[1]
  if (is.na(data_start)) {
    cat(sprintf("WARN: GONE2 -x header not found in %s\n", path))
    return(NULL)
  }
  d <- tryCatch(
    read.table(text  = lines[(data_start + 1):length(lines)],
               header = FALSE, stringsAsFactors = FALSE,
               col.names = c("rec_rate", "generation", "N_T", "Ne", "d2")),
    error = function(e) {
      cat(sprintf("ERROR reading GONE2 -x %s: %s\n", path, e$message))
      NULL
    }
  )
  if (is.null(d) || nrow(d) == 0) return(NULL)
  data.frame(generation = as.integer(d$generation),
             Ne         = as.numeric(d$Ne),
             stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------
# STEP 1: Collect per-population results for all four tools
# -------------------------------------------------------------------
cat("INFO: parsing currentNe (v1) outputs...\n")
cne1_list <- lapply(pops, function(pop) {
  f <- file.path(pop_dir, pop, "currentNe", paste0(pop, "_currentNe.txt"))
  if (!file.exists(f))
    f <- file.path(pop_dir, pop, "currentNe", paste0(pop, "_currentNe"))
  d <- parse_currentNe(f)
  if (!is.null(d)) { d$population <- pop; d$n_ind <- n_ind_for(pop) }
  d
})
names(cne1_list) <- pops
cne1_valid <- Filter(Negate(is.null), cne1_list)

cat("INFO: parsing currentNe2 (v2) outputs...\n")
cne2_default_list <- lapply(pops, function(pop) {
  f <- file.path(pop_dir, pop, "currentNe2", paste0(pop, "_currentNe2_default.txt"))
  d <- parse_currentNe2(f)
  if (!is.null(d)) { d$population <- pop; d$n_ind <- n_ind_for(pop) }
  d
})
names(cne2_default_list) <- pops
cne2_default_valid <- Filter(Negate(is.null), cne2_default_list)

cne2_mix_list <- lapply(pops, function(pop) {
  f <- file.path(pop_dir, pop, "currentNe2", paste0(pop, "_currentNe2_mix.txt"))
  d <- parse_currentNe2(f)
  if (!is.null(d)) { d$population <- pop; d$n_ind <- n_ind_for(pop) }
  d
})
names(cne2_mix_list) <- pops
cne2_mix_valid <- Filter(Negate(is.null), cne2_mix_list)

cat("INFO: parsing GONE (v1) outputs...\n")
gone1_list <- lapply(pops, function(pop) {
  f <- file.path(pop_dir, pop, "GONE", paste0("Output_Ne_", pop))
  d <- parse_GONE_file(f)
  if (!is.null(d)) d$population <- pop
  d
})
names(gone1_list) <- pops
gone1_valid <- Filter(Negate(is.null), gone1_list)

cat("INFO: parsing GONE2 default outputs...\n")
gone2_def_list <- lapply(pops, function(pop) {
  f <- file.path(pop_dir, pop, "GONE_2", "default",
                 paste0(pop, ".ped_GONE2_Ne"))
  d <- parse_GONE_file(f)
  if (!is.null(d)) d$population <- pop
  d
})
names(gone2_def_list) <- pops
gone2_def_valid <- Filter(Negate(is.null), gone2_def_list)

cat("INFO: parsing GONE2 -x (structure) outputs...\n")
gone2_str_list <- lapply(pops, function(pop) {
  f <- file.path(pop_dir, pop, "GONE_2", "structure",
                 paste0(pop, ".ped_GONE2_Ne_mix"))
  d <- parse_GONE2_mix(f)
  if (!is.null(d)) d$population <- pop
  d
})
names(gone2_str_list) <- pops
gone2_str_valid <- Filter(Negate(is.null), gone2_str_list)

# -------------------------------------------------------------------
# STEP 2: currentNe v1 — summary table + barplot
# -------------------------------------------------------------------
if (length(cne1_valid) > 0) {
  cne1_df <- do.call(rbind, cne1_valid)
  cne1_df <- cne1_df[, c("population", "n_ind", "Ne",
                           "CI50_lower", "CI50_upper",
                           "CI90_lower", "CI90_upper")]
  cne1_df[sapply(cne1_df, is.numeric)] <-
    lapply(cne1_df[sapply(cne1_df, is.numeric)], round, 0)
  tsv <- file.path(cne_out_dir,
                   sprintf("07.8-currentNe_summary_K%s.tsv", k_global))
  write.table(cne1_df, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("currentNe v1 table written:", tsv, "\n")
  print(cne1_df)

  df_p <- cne1_df[!is.na(cne1_df$Ne), ]
  df_p$population <- factor(df_p$population, levels = sort(unique(df_p$population)))
  p <- ggplot(df_p, aes(x = population, y = Ne, fill = population)) +
    geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
    geom_errorbar(aes(ymin = CI90_lower, ymax = CI90_upper),
                  width = 0.18, color = "grey20", linewidth = 0.7) +
    geom_errorbar(aes(ymin = CI50_lower, ymax = CI50_upper),
                  width = 0.10, color = "black", linewidth = 1.1) +
    geom_text(aes(label = format(round(Ne, 0), big.mark = ",")),
              vjust = -0.5, size = 3.2) +
    scale_fill_manual(values = POP_COLORS, na.value = "grey50") +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(title    = sprintf("Current Ne — currentNe v1 — K=%s", k_global),
         subtitle = "Error bars: thick = 50% CI, thin = 90% CI",
         x = NULL, y = "Ne (currentNe v1)") +
    theme(legend.position = "none", axis.text.x = element_text(size = 10),
          plot.subtitle = element_text(size = 8, color = "grey40"))
  fn <- file.path(cne_out_dir,
                  sprintf("currentNe_v1_K%s.png", k_global))
  ggsave(fn, p, width = max(6, length(unique(df_p$population)) * 1.8),
         height = 5, dpi = 300)
  cat("Plot written:", fn, "\n")
} else {
  cat("WARN: no valid currentNe v1 outputs found\n")
  cne1_df <- NULL
}

# -------------------------------------------------------------------
# STEP 3: currentNe2 — summary table + barplots (default + structure)
# -------------------------------------------------------------------
build_cne2_table <- function(valid_list, label) {
  if (length(valid_list) == 0) return(NULL)
  df <- do.call(rbind, valid_list)
  cols_base <- c("population", "n_ind", "Ne",
                 "CI50_lower", "CI50_upper", "CI90_lower", "CI90_upper")
  cols_struct <- c("NT", "FST", "migration", "n_subpop")
  cols_use <- c(cols_base,
                cols_struct[cols_struct %in% colnames(df) &
                              !all(is.na(df[, cols_struct[cols_struct %in% colnames(df)]]))])
  df <- df[, intersect(cols_use, colnames(df))]
  df[sapply(df, is.numeric)] <- lapply(df[sapply(df, is.numeric)], round, 4)
  tsv <- file.path(cne2_out_dir,
                   sprintf("07.8-currentNe2_%s_K%s.tsv", label, k_global))
  write.table(df, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("currentNe2 %s table written: %s\n", label, tsv))
  print(df)
  df
}

cne2_def_df <- build_cne2_table(cne2_default_valid, "default")
cne2_mix_df <- build_cne2_table(cne2_mix_valid,     "structure")

plot_cne2_barplot <- function(df, subtitle_txt, fn_suffix) {
  if (is.null(df) || all(is.na(df$Ne))) return(invisible(NULL))
  df <- df[!is.na(df$Ne), ]
  df$population <- factor(df$population, levels = sort(unique(df$population)))
  p <- ggplot(df, aes(x = population, y = Ne, fill = population)) +
    geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
    geom_errorbar(aes(ymin = CI90_lower, ymax = CI90_upper),
                  width = 0.18, color = "grey20", linewidth = 0.7) +
    geom_errorbar(aes(ymin = CI50_lower, ymax = CI50_upper),
                  width = 0.10, color = "black", linewidth = 1.1) +
    geom_text(aes(label = format(round(Ne, 0), big.mark = ",")),
              vjust = -0.5, size = 3.2) +
    scale_fill_manual(values = POP_COLORS, na.value = "grey50") +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(title    = sprintf("Current Ne — currentNe2 %s — K=%s",
                            fn_suffix, k_global),
         subtitle = subtitle_txt,
         x = NULL, y = "Ne (currentNe2)") +
    theme(legend.position = "none", axis.text.x = element_text(size = 10),
          plot.subtitle = element_text(size = 8, color = "grey40"))
  fn <- file.path(cne2_out_dir,
                  sprintf("currentNe2_%s_K%s.png", fn_suffix, k_global))
  ggsave(fn, p, width = max(6, nrow(df) * 1.8), height = 5, dpi = 300)
  cat("Plot written:", fn, "\n")

  # If mix output has FST, add a secondary barplot for structure params
  if (!is.null(df$FST) && any(!is.na(df$FST))) {
    struct_cols <- intersect(c("FST", "migration", "n_subpop"), colnames(df))
    struct_long <- reshape(
      df[, c("population", struct_cols)],
      varying = struct_cols, times = struct_cols,
      v.names = "value", timevar = "param",
      idvar = "population", direction = "long"
    )
    struct_long <- struct_long[!is.na(struct_long$value), ]
    if (nrow(struct_long) > 0) {
      ps <- ggplot(struct_long,
                   aes(x = population, y = value, fill = population)) +
        geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
        facet_wrap(~param, scales = "free_y") +
        scale_fill_manual(values = POP_COLORS, na.value = "grey50") +
        labs(title    = sprintf("Structure parameters — currentNe2 -x — K=%s",
                                k_global),
             x = NULL, y = "Estimate") +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 9))
      fn_s <- file.path(cne2_out_dir,
                        sprintf("currentNe2_structure_params_K%s.png", k_global))
      ggsave(fn_s, ps, width = max(8, nrow(df) * 2), height = 5, dpi = 300)
      cat("Plot written:", fn_s, "\n")
    }
  }
}

plot_cne2_barplot(cne2_def_df,
                  "Error bars: thick = 50% CI, thin = 90% CI | panmictic assumption",
                  "default")
plot_cne2_barplot(cne2_mix_df,
                  "Error bars: thick = 50% CI, thin = 90% CI | population structure assumed (-x)",
                  "structure")

# -------------------------------------------------------------------
# STEP 4: GONE v1 — summary table + trajectory plots
# -------------------------------------------------------------------
if (length(gone1_valid) > 0) {
  gone1_df <- do.call(rbind, gone1_valid)
  tsv <- file.path(g1_out_dir,
                   sprintf("07.8-GONE_v1_trajectory_K%s.tsv", k_global))
  write.table(gone1_df, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("GONE v1 trajectory table written:", tsv, "\n")

  gone1_df$population <- factor(gone1_df$population,
                                 levels = sort(unique(gone1_df$population)))
  gone1_plot <- gone1_df[gone1_df$generation <= 500 &
                           !is.na(gone1_df$Ne) & gone1_df$Ne > 0, ]
  if (nrow(gone1_plot) > 0) {
    p <- ggplot(gone1_plot, aes(x = generation, y = Ne,
                                color = population, group = population)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = POP_COLORS, na.value = "grey50") +
      scale_x_continuous(breaks = c(1, 10, 50, 100, 200, 500)) +
      scale_y_log10(labels = scales::comma) +
      labs(title    = sprintf("Temporal Ne trajectory — GONE v1 — K=%s", k_global),
           subtitle = "Generations 1–500 — log10(Ne) scale",
           x = "Generations ago", y = "Ne (GONE v1, log scale)",
           color = "Cluster") +
      theme(plot.subtitle = element_text(size = 8, color = "grey40"))
    fn <- file.path(g1_out_dir,
                    sprintf("GONE_v1_trajectory_K%s.png", k_global))
    ggsave(fn, p, width = 8, height = 5, dpi = 300)
    cat("Plot written:", fn, "\n")

    gone1_all <- gone1_df[!is.na(gone1_df$Ne) & gone1_df$Ne > 0, ]
    if (max(gone1_all$generation) > 500) {
      pfull <- p %+% gone1_all +
        scale_x_continuous() +
        labs(subtitle = "All generations — log10(Ne) scale")
      fn2 <- file.path(g1_out_dir,
                       sprintf("GONE_v1_trajectory_full_K%s.png", k_global))
      ggsave(fn2, pfull, width = 9, height = 5, dpi = 300)
      cat("Plot written:", fn2, "\n")
    }
  }
} else {
  cat("WARN: no valid GONE v1 outputs found\n")
  gone1_df <- NULL
}

# -------------------------------------------------------------------
# STEP 5: GONE2 — summary tables + trajectory plots (default + -x)
# -------------------------------------------------------------------
plot_gone2_trajectory <- function(def_valid, str_valid, label) {
  has_def <- length(def_valid) > 0
  has_str <- length(str_valid) > 0
  if (!has_def && !has_str) return(invisible(NULL))

  build_df <- function(vlist, run_type) {
    df <- do.call(rbind, vlist)
    df$run <- run_type
    df
  }

  if (has_def) {
    tsv <- file.path(g2_out_dir,
                     sprintf("07.8-GONE2_default_trajectory_K%s.tsv", k_global))
    write.table(do.call(rbind, def_valid), tsv,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("GONE2 default trajectory table written:", tsv, "\n")
  }
  if (has_str) {
    tsv <- file.path(g2_out_dir,
                     sprintf("07.8-GONE2_structure_trajectory_K%s.tsv", k_global))
    write.table(do.call(rbind, str_valid), tsv,
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("GONE2 structure trajectory table written:", tsv, "\n")
  }

  # Overlay default vs -x per population (500-generation window)
  if (has_def && has_str) {
    combined <- rbind(build_df(def_valid, "default"),
                      build_df(str_valid, "structure (-x)"))
    combined$population <- factor(combined$population,
                                  levels = sort(unique(combined$population)))
    combined_plot <- combined[combined$generation <= 500 &
                                !is.na(combined$Ne) & combined$Ne > 0, ]
    if (nrow(combined_plot) > 0) {
      p <- ggplot(combined_plot,
                  aes(x = generation, y = Ne,
                      color = population, linetype = run, group = interaction(population, run))) +
        geom_line(linewidth = 0.8) +
        scale_color_manual(values = POP_COLORS, na.value = "grey50") +
        scale_linetype_manual(values = c("default" = "solid",
                                         "structure (-x)" = "dashed")) +
        scale_x_continuous(breaks = c(1, 10, 50, 100, 200, 500)) +
        scale_y_log10(labels = scales::comma) +
        labs(title    = sprintf("Temporal Ne trajectory — GONE2 — K=%s", k_global),
             subtitle = "Solid = panmictic (default) | Dashed = structure-corrected (-x) | Generations 1–500",
             x = "Generations ago", y = "Ne (GONE2, log scale)",
             color = "Cluster", linetype = "Model") +
        theme(plot.subtitle = element_text(size = 8, color = "grey40"))
      fn <- file.path(g2_out_dir,
                      sprintf("GONE2_trajectory_overlay_K%s.png", k_global))
      ggsave(fn, p, width = 9, height = 5, dpi = 300)
      cat("Plot written:", fn, "\n")

      # Full trajectory
      combined_all <- combined[!is.na(combined$Ne) & combined$Ne > 0, ]
      if (max(combined_all$generation) > 500) {
        pfull <- p %+% combined_all +
          scale_x_continuous() +
          labs(subtitle = "Solid = panmictic | Dashed = structure-corrected | All generations")
        fn2 <- file.path(g2_out_dir,
                         sprintf("GONE2_trajectory_overlay_full_K%s.png", k_global))
        ggsave(fn2, pfull, width = 10, height = 5, dpi = 300)
        cat("Plot written:", fn2, "\n")
      }
    }
  } else {
    # Only one mode available — single trajectory
    vlist <- if (has_def) def_valid else str_valid
    run_lbl <- if (has_def) "default" else "structure (-x)"
    df_all <- do.call(rbind, vlist)
    df_all$population <- factor(df_all$population,
                                levels = sort(unique(df_all$population)))
    df_plot <- df_all[df_all$generation <= 500 & !is.na(df_all$Ne) & df_all$Ne > 0, ]
    if (nrow(df_plot) > 0) {
      p <- ggplot(df_plot, aes(x = generation, y = Ne,
                               color = population, group = population)) +
        geom_line(linewidth = 0.9) +
        scale_color_manual(values = POP_COLORS, na.value = "grey50") +
        scale_x_continuous(breaks = c(1, 10, 50, 100, 200, 500)) +
        scale_y_log10(labels = scales::comma) +
        labs(title    = sprintf("Temporal Ne trajectory — GONE2 %s — K=%s",
                                run_lbl, k_global),
             x = "Generations ago", y = "Ne (GONE2, log scale)",
             color = "Cluster")
      fn <- file.path(g2_out_dir,
                      sprintf("GONE2_trajectory_%s_K%s.png",
                              gsub(" |\\(|\\)", "_", run_lbl), k_global))
      ggsave(fn, p, width = 8, height = 5, dpi = 300)
      cat("Plot written:", fn, "\n")
    }
  }
}

plot_gone2_trajectory(gone2_def_valid, gone2_str_valid, "GONE2")

# -------------------------------------------------------------------
# STEP 6: All-tools current Ne comparison barplot
# Assembles gen-1 Ne from GONE/GONE2 and point-estimate Ne from
# currentNe/currentNe2 into a single grouped barplot.
# -------------------------------------------------------------------
get_gen1_Ne <- function(valid_list) {
  if (length(valid_list) == 0) return(NULL)
  do.call(rbind, lapply(valid_list, function(d) {
    d_sorted <- d[order(d$generation), ]
    data.frame(population = d_sorted$population[1],
               Ne = d_sorted$Ne[1],
               stringsAsFactors = FALSE)
  }))
}

comp_rows <- list()

if (!is.null(cne1_df)) {
  r <- cne1_df[!is.na(cne1_df$Ne), c("population", "Ne")]
  r$tool <- "currentNe v1"
  comp_rows[[length(comp_rows) + 1]] <- r
}
if (!is.null(cne2_def_df) && any(!is.na(cne2_def_df$Ne))) {
  r <- cne2_def_df[!is.na(cne2_def_df$Ne), c("population", "Ne")]
  r$tool <- "currentNe2 default"
  comp_rows[[length(comp_rows) + 1]] <- r
}
if (!is.null(cne2_mix_df) && any(!is.na(cne2_mix_df$Ne))) {
  r <- cne2_mix_df[!is.na(cne2_mix_df$Ne), c("population", "Ne")]
  r$tool <- "currentNe2 -x"
  comp_rows[[length(comp_rows) + 1]] <- r
}
g1_curr <- get_gen1_Ne(gone1_valid)
if (!is.null(g1_curr)) {
  g1_curr$tool <- "GONE v1 (gen 1)"
  comp_rows[[length(comp_rows) + 1]] <- g1_curr
}
g2_curr <- get_gen1_Ne(gone2_def_valid)
if (!is.null(g2_curr)) {
  g2_curr$tool <- "GONE2 default (gen 1)"
  comp_rows[[length(comp_rows) + 1]] <- g2_curr
}
g2s_curr <- get_gen1_Ne(gone2_str_valid)
if (!is.null(g2s_curr)) {
  g2s_curr$tool <- "GONE2 -x (gen 1)"
  comp_rows[[length(comp_rows) + 1]] <- g2s_curr
}

if (length(comp_rows) > 0) {
  comp_df <- do.call(rbind, comp_rows)
  comp_df <- comp_df[!is.na(comp_df$Ne), ]
  comp_df$population <- factor(comp_df$population,
                                levels = sort(unique(comp_df$population)))
  comp_df$tool <- factor(comp_df$tool, levels = unique(comp_df$tool))

  tool_alphas <- setNames(
    seq(1.0, 0.3, length.out = nlevels(comp_df$tool)),
    levels(comp_df$tool)
  )

  p_comp <- ggplot(comp_df,
                   aes(x = population, y = Ne,
                       fill = population, alpha = tool)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65,
             color = "grey30", linewidth = 0.25) +
    scale_fill_manual(values = POP_COLORS, na.value = "grey50") +
    scale_alpha_manual(values = tool_alphas) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.12))) +
    labs(title = sprintf("Current Ne — all tools comparison — K=%s", k_global),
         subtitle = "Grouped by tool (transparency). GONE/GONE2 values = generation 1.",
         x = NULL, y = "Ne", fill = "Cluster", alpha = "Tool") +
    theme(axis.text.x = element_text(size = 10),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          legend.key.size = unit(0.5, "cm"))

  fn_comp <- file.path(out_dir,
                       sprintf("07.8-Ne_all_tools_K%s.png", k_global))
  ggsave(fn_comp, p_comp,
         width = max(8, nlevels(comp_df$population) * 3),
         height = 5, dpi = 300)
  cat("All-tools comparison plot written:", fn_comp, "\n")

  tsv_comp <- file.path(out_dir,
                        sprintf("07.8-Ne_all_tools_K%s.tsv", k_global))
  write.table(comp_df, tsv_comp, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("All-tools table written:", tsv_comp, "\n")
}

cat(sprintf("DONE 07.8 Ne R script completed — K=%s\n", k_global))

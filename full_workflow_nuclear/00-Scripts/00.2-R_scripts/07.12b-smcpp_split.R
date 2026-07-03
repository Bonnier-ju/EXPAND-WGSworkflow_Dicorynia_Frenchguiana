#!/usr/bin/env Rscript
# 07.12b-smcpp_split.R
# Visualises SMC++ pairwise split analysis for the 4 ADMIXTURE clusters.
# Reads Ne(t) CSV files from smc++ split and produces:
#   1. Joint Ne(t) panel grid — 6 pairs, each showing both populations
#      and their shared ancestral Ne trajectory
#   2. Divergence-time heatmap — T_split (years BP) for all C(4,2)=6 pairs
#   3. Pop_3 focus — 3 panels showing Pop_3 vs each other population
#
# Note: smc++ split CSV format (per pair):
#   label, x (years), y (Ne)
#   label = pop name (Pop_A or Pop_B) with two trajectories per pair
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

pairs <- list(
  c("Pop_1", "Pop_2"), c("Pop_1", "Pop_3"), c("Pop_1", "Pop_4"),
  c("Pop_2", "Pop_3"), c("Pop_2", "Pop_4"), c("Pop_3", "Pop_4")
)

plot_dir <- file.path(out_base, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Load split CSV files
# smc++ split outputs two trajectories per pair (one per population).
# The CSV from `smc++ plot --csv` has: label, x (years), y (Ne).
# The split time can be inferred as the year at which the two
# trajectories converge (merge into the same Ne value).
# -------------------------------------------------------------------
load_split <- function(pop_a, pop_b, out_base) {
  f <- file.path(out_base, "csv", "split",
                 paste0("Ne_", pop_a, "_", pop_b, ".csv"))
  if (!file.exists(f)) {
    cat(sprintf("WARN: split CSV not found: %s\n", f))
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

  names(df) <- tolower(names(df))
  if (!"x" %in% names(df) || !"y" %in% names(df)) {
    cat(sprintf("WARN: unexpected columns in %s: %s\n",
                f, paste(names(df), collapse = ", ")))
    return(NULL)
  }

  df <- df[!is.na(df$x) & df$x >= 10 & df$x <= 5e6, ]
  df$pair    <- paste0(pop_a, "/", pop_b)
  df$pop_a   <- pop_a
  df$pop_b   <- pop_b
  df$year    <- df$x
  df$Ne      <- df$y
  # label column from smc++ identifies the population within the pair
  df$population <- if ("label" %in% names(df)) df$label else NA_character_
  df[, c("pair", "pop_a", "pop_b", "population", "year", "Ne")]
}

split_list <- lapply(pairs, function(p) load_split(p[1], p[2], out_base))
split_valid <- Filter(Negate(is.null), split_list)

if (length(split_valid) == 0) {
  cat("ERROR: no split CSV files found — exiting\n")
  quit(status = 1)
}

cat(sprintf("INFO: %d/%d pairs loaded\n", length(split_valid), length(pairs)))

# -------------------------------------------------------------------
# Infer T_split: the oldest year at which the two trajectories DIVERGE
# (going forward in time). We define T_split as the year where the
# relative difference between the two population Ne values drops below
# a threshold (5%), i.e., they are effectively merged.
# -------------------------------------------------------------------
infer_tsplit <- function(df) {
  pops_in <- unique(df$population)
  if (length(pops_in) < 2 || any(is.na(pops_in))) return(NA_real_)
  d1 <- df[df$population == pops_in[1], ]
  d2 <- df[df$population == pops_in[2], ]
  # Interpolate on common time grid
  tgrid <- sort(unique(c(d1$year, d2$year)))
  ne1 <- approx(d1$year, d1$Ne, xout = tgrid, rule = 2)$y
  ne2 <- approx(d2$year, d2$Ne, xout = tgrid, rule = 2)$y
  rel_diff <- abs(ne1 - ne2) / pmax(ne1, ne2)
  # T_split = oldest time point where trajectories still diverge (>5%)
  diverged <- tgrid[rel_diff > 0.05]
  if (length(diverged) == 0) return(min(tgrid))
  max(diverged)
}

tsplit_df <- do.call(rbind, lapply(split_valid, function(df) {
  ts <- infer_tsplit(df)
  data.frame(pair  = df$pair[1],
             pop_a = df$pop_a[1],
             pop_b = df$pop_b[1],
             T_split_yr = ts,
             stringsAsFactors = FALSE)
}))
tsplit_df$T_split_kyr <- tsplit_df$T_split_yr / 1000

cat("Split times (years BP):\n")
print(tsplit_df[, c("pair", "T_split_yr")])

# -------------------------------------------------------------------
# FIGURE 1: Joint Ne(t) panels — 6 pairs, 2 trajectories each
# -------------------------------------------------------------------
cat("Building Figure 1: joint Ne(t) panels...\n")

split_df <- do.call(rbind, split_valid)
# Assign population colors: population column should contain pop names
split_df$pop_color_key <- ifelse(
  split_df$population %in% pops,
  split_df$population,
  split_df$pop_a
)
split_df$pop_color_key <- factor(split_df$pop_color_key, levels = pops)

p1 <- ggplot(split_df,
             aes(x = year, y = Ne,
                 color = pop_color_key,
                 group = interaction(pair, population))) +
  geom_line(linewidth = 0.8) +
  # Annotate T_split
  geom_vline(data = tsplit_df,
             aes(xintercept = T_split_yr),
             linetype = "dotted", color = "grey50", linewidth = 0.5) +
  scale_color_manual(values = POP_COLORS, name = "Population",
                     drop = FALSE) +
  scale_x_log10(name = "Time (years before present)") +
  scale_y_log10(name = expression(N[e])) +
  annotation_logticks(sides = "bl", size = 0.25) +
  facet_wrap(~ pair, ncol = 3, scales = "free_y") +
  labs(
    title    = sprintf("Pairwise joint Ne(t) — SMC++ split (K=%d)", k_global),
    subtitle = "Two trajectories per panel | Dotted = inferred T_split | μ = 5×10⁻⁹ | 30 yr/gen"
  ) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 8),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(plot_dir, "Ne_split_joint_panels.png"), p1,
       width = 14, height = 8, dpi = 300)
ggsave(file.path(plot_dir, "Ne_split_joint_panels.pdf"), p1,
       width = 14, height = 8)
cat("Figure 1 written\n")

# -------------------------------------------------------------------
# FIGURE 2: Divergence-time heatmap — T_split for all 6 pairs
# Symmetric matrix showing T_split (kyr BP) between each population pair.
# -------------------------------------------------------------------
cat("Building Figure 2: T_split heatmap...\n")

# Build symmetric matrix
all_pops_in_split <- unique(c(tsplit_df$pop_a, tsplit_df$pop_b))
all_pops_in_split <- intersect(pops, all_pops_in_split)

mat_rows <- do.call(rbind, lapply(all_pops_in_split, function(pa) {
  do.call(rbind, lapply(all_pops_in_split, function(pb) {
    if (pa == pb) return(data.frame(pop_a = pa, pop_b = pb,
                                    T_split_kyr = NA_real_,
                                    stringsAsFactors = FALSE))
    ts_row <- tsplit_df[
      (tsplit_df$pop_a == pa & tsplit_df$pop_b == pb) |
      (tsplit_df$pop_a == pb & tsplit_df$pop_b == pa), ]
    ts_val <- if (nrow(ts_row) > 0) ts_row$T_split_kyr[1] else NA_real_
    data.frame(pop_a = pa, pop_b = pb, T_split_kyr = ts_val,
               stringsAsFactors = FALSE)
  }))
}))
mat_rows$pop_a <- factor(mat_rows$pop_a, levels = pops)
mat_rows$pop_b <- factor(mat_rows$pop_b, levels = rev(pops))
mat_rows$label <- ifelse(is.na(mat_rows$T_split_kyr),
                         "—",
                         sprintf("%.0f kyr", mat_rows$T_split_kyr))

p2 <- ggplot(mat_rows, aes(x = pop_a, y = pop_b)) +
  geom_tile(aes(fill = T_split_kyr), color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3.5, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#1a3d6b", high = "#c8d8f0",
                      na.value = "grey85",
                      name = "T_split\n(kyr BP)") +
  labs(
    title    = sprintf("Divergence times between clusters — SMC++ split (K=%d)", k_global),
    subtitle = "T_split estimated as oldest year with > 5% relative Ne divergence between populations",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text  = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    panel.grid = element_blank()
  )

ggsave(file.path(plot_dir, "split_tsplit_heatmap.png"), p2,
       width = 7, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "split_tsplit_heatmap.pdf"), p2,
       width = 7, height = 6)
cat("Figure 2 written\n")

# -------------------------------------------------------------------
# FIGURE 3: Pop_3 focus — 3 panels (Pop_3 vs Pop_1, Pop_2, Pop_4)
# Highlights Pop_3 trajectory relative to each other population,
# directly addressing the GONE vs Stairway Plot discrepancy question.
# -------------------------------------------------------------------
cat("Building Figure 3: Pop_3 focus panels...\n")

pop3_pairs <- c("Pop_1/Pop_3", "Pop_2/Pop_3", "Pop_3/Pop_4")
pop3_split <- split_df[split_df$pair %in% pop3_pairs, ]

if (nrow(pop3_split) > 0) {
  pop3_tsplit <- tsplit_df[
    tsplit_df$pop_a == "Pop_3" | tsplit_df$pop_b == "Pop_3", ]

  p3 <- ggplot(pop3_split,
               aes(x = year, y = Ne,
                   color = pop_color_key,
                   group = interaction(pair, population))) +
    geom_line(linewidth = 1.0) +
    geom_vline(data = pop3_tsplit,
               aes(xintercept = T_split_yr),
               linetype = "dotted", color = "grey40", linewidth = 0.6) +
    scale_color_manual(values = POP_COLORS, name = "Population",
                       drop = FALSE) +
    scale_x_log10(name = "Time (years before present)") +
    scale_y_log10(name = expression(N[e])) +
    annotation_logticks(sides = "bl", size = 0.3) +
    facet_wrap(~ pair, ncol = 3, scales = "free_y") +
    labs(
      title    = sprintf("Pop_3 pairwise history — SMC++ split (K=%d)", k_global),
      subtitle = "Dotted = T_split | Red = Pop_3 trajectory | μ = 5×10⁻⁹ | 30 yr/gen"
    ) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )

  ggsave(file.path(plot_dir, "Ne_split_Pop3_focus.png"), p3,
         width = 13, height = 5, dpi = 300)
  ggsave(file.path(plot_dir, "Ne_split_Pop3_focus.pdf"), p3,
         width = 13, height = 5)
  cat("Figure 3 written\n")
} else {
  cat("WARN: Figure 3 skipped (no Pop_3 split pairs found)\n")
}

# -------------------------------------------------------------------
# Write T_split summary table
# -------------------------------------------------------------------
tsv_out <- file.path(out_base, "csv", "split", "07.12b-tsplit_summary.tsv")
write.table(tsplit_df[order(tsplit_df$T_split_yr), ],
            tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("T_split table written:", tsv_out, "\n")

cat(sprintf("DONE 07.12b-smcpp_split.R — plots: %s\n", plot_dir))

library(ggplot2)
library(dplyr)

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 06.6a-ld_decay.R <geno.ld file> <output plot dir>")
}
ld_file  <- args[1]
plot_dir <- args[2]

# -------------------------------------------------------------------
# Read LD table
# -------------------------------------------------------------------
# vcftools --geno-r2 output format:
#   CHR  POS1  POS2  N_INDV  R^2
# R^2 is the squared Pearson correlation between allele counts at two sites.
# N_INDV is the number of individuals with non-missing genotypes at both sites.
message("Reading LD table: ", ld_file)
ld <- read.table(ld_file, header = TRUE, sep = "\t",
                 col.names = c("CHR", "POS1", "POS2", "N_INDV", "R2"))

# Remove pairs with missing r² or fewer than 10 individuals
ld <- ld[!is.na(ld$R2) & ld$N_INDV >= 10, ]

# Compute pairwise distance (bp)
ld$DIST <- abs(ld$POS2 - ld$POS1)

message(sprintf("Total SNP pairs: %d", nrow(ld)))
message(sprintf("Distance range: %d - %d bp", min(ld$DIST), max(ld$DIST)))

# -------------------------------------------------------------------
# Bin r² by distance
# -------------------------------------------------------------------
# Bins of 10 kb up to 1 Mb to capture full LD decay curve.
# Mean r² per bin is the standard summary used in LD decay plots.
BREAKS <- c(0, 1e3, 5e3, 10e3, 25e3, 50e3, 75e3,
            100e3, 200e3, 300e3, 400e3, 500e3, 750e3, 1e6)

ld$DIST_BIN <- cut(ld$DIST, breaks = BREAKS, include.lowest = TRUE)

ld_bins <- ld %>%
  group_by(DIST_BIN) %>%
  summarise(
    MEAN_R2   = mean(R2, na.rm = TRUE),
    MEDIAN_R2 = median(R2, na.rm = TRUE),
    N_PAIRS   = n(),
    MID_DIST  = mean(DIST, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(DIST_BIN))

# -------------------------------------------------------------------
# Identify distance where mean r² drops below threshold
# -------------------------------------------------------------------
# r² = 0.2 is the conventional threshold used to define the boundary
# of effective LD blocks and to set the LD pruning window size.
# The pruning window should be >= this distance to ensure all
# correlated SNP pairs are captured.
R2_THRESHOLD <- 0.2

# Find the first bin where mean r² drops below threshold
below_threshold <- ld_bins[ld_bins$MEAN_R2 < R2_THRESHOLD, ]

if (nrow(below_threshold) > 0) {
  recommended_window_bp <- ceiling(min(below_threshold$MID_DIST))
  recommended_window_kb <- ceiling(recommended_window_bp / 1000)
  message(sprintf("Mean r² drops below %.1f at ~%d bp (%d kb)",
                  R2_THRESHOLD, recommended_window_bp, recommended_window_kb))
} else {
  recommended_window_kb <- 1000  # default fallback
  message("WARNING: mean r² never drops below threshold within 1 Mb window")
}

# -------------------------------------------------------------------
# LD decay plot
# -------------------------------------------------------------------
p <- ggplot(ld_bins, aes(x = MID_DIST / 1000, y = MEAN_R2)) +
  geom_point(aes(size = N_PAIRS), alpha = 0.7, color = "#2c7bb6") +
  geom_line(color = "#2c7bb6", linewidth = 0.8) +
  geom_smooth(method = "loess", span = 0.4, se = TRUE,
              color = "#d7191c", fill = "#d7191c", alpha = 0.15, linewidth = 1) +
  geom_hline(yintercept = R2_THRESHOLD, linetype = "dashed",
             color = "grey40", linewidth = 0.7) +
  annotate("text", x = max(ld_bins$MID_DIST / 1000) * 0.7,
           y = R2_THRESHOLD + 0.02,
           label = paste0("r² = ", R2_THRESHOLD, " threshold"),
           color = "grey40", size = 3.5) +
  {if (nrow(below_threshold) > 0)
    geom_vline(xintercept = recommended_window_kb, linetype = "dotted",
               color = "#1a9641", linewidth = 0.8)} +
  {if (nrow(below_threshold) > 0)
    annotate("text", x = recommended_window_kb + max(ld_bins$MID_DIST / 1000) * 0.02,
             y = max(ld_bins$MEAN_R2) * 0.9,
             label = paste0("~", recommended_window_kb, " kb"),
             color = "#1a9641", size = 3.5, hjust = 0)} +
  scale_size_continuous(name = "N pairs", range = c(1, 5)) +
  labs(
    title = "LD decay — Dicorynia guianensis nuclear SNPs",
    subtitle = sprintf("Based on %s randomly sampled SNP pairs | r² threshold = %.1f",
                       format(nrow(ld), big.mark = ","), R2_THRESHOLD),
    x = "Pairwise distance (kb)",
    y = expression(paste("Mean ", r^2))
  ) +
  scale_x_continuous(labels = scales::comma) +
  coord_cartesian(ylim = c(0, NA))

ggsave(file.path(plot_dir, "ld_decay_curve.png"),
       plot = p, width = 10, height = 6, dpi = 300)
message("Plot saved: ld_decay_curve.png")

# -------------------------------------------------------------------
# Write window recommendation
# -------------------------------------------------------------------
rec_file <- file.path(plot_dir, "ld_window_recommendation.txt")
writeLines(c(
  "# LD decay analysis — window size recommendation for 06.6 LD pruning",
  sprintf("# Based on %d SNP pairs from a random subsample of SNPs", nrow(ld)),
  sprintf("# r² threshold used: %.1f", R2_THRESHOLD),
  "",
  sprintf("Recommended LD_WINDOW_KB for 06.6: %d", recommended_window_kb),
  "",
  "# Summary of mean r² per distance bin:",
  paste(sprintf("%-20s  MEAN_R2=%.4f  MEDIAN_R2=%.4f  N=%d",
                as.character(ld_bins$DIST_BIN),
                ld_bins$MEAN_R2,
                ld_bins$MEDIAN_R2,
                ld_bins$N_PAIRS), collapse = "\n")
), rec_file)

message("Recommendation written: ", rec_file)
message(sprintf("Suggested LD_WINDOW_KB for 06.6: %d", recommended_window_kb))

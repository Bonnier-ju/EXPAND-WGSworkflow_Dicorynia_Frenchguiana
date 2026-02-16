# coverage_summary_plot.R
# Usage: Rscript 04.1-coverage_summary_plot.R input_tsv output_dir

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 04.1-coverage_summary_plot.R <input_tsv> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load required libraries
library(ggplot2)
library(readr)

# Read input coverage file
coverage <- read_tsv(input_file, col_types = cols())

# Current workflow format only:
# sample_id, total_reads, mapped_reads, mapped_percent, avg_depth_q20
required_cols <- c("sample_id", "mapped_percent", "avg_depth_q20")
if (!all(required_cols %in% colnames(coverage))) {
  stop("Unsupported input format. Missing columns: ",
       paste(setdiff(required_cols, colnames(coverage)), collapse = ", "))
}

coverage$Sample <- coverage$sample_id
coverage$Average_Coverage <- suppressWarnings(as.numeric(coverage$avg_depth_q20))
coverage$Mapped_Percent <- suppressWarnings(as.numeric(coverage$mapped_percent))

# Keep rows with valid coverage values for plotting
coverage_plot <- coverage[!is.na(coverage$Average_Coverage), , drop = FALSE]
if (nrow(coverage_plot) == 0) {
  stop("No valid coverage values found in input file.")
}

# Plot coverage per sample with threshold line
p_cov <- ggplot(coverage_plot, aes(x = reorder(Sample, -Average_Coverage), y = Average_Coverage)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(Average_Coverage, 1)), vjust = -0.5, size = 3) +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 1, y = 10, label = "10x threshold", vjust = -1.2, hjust = 0, color = "red", size = 3) +
  labs(title = "Average Coverage per Sample",
       x = "Sample",
       y = "Average Coverage (Q >= 20)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save coverage plot
coverage_png <- file.path(output_dir, "coverage_barplot.png")
ggsave(filename = coverage_png, plot = p_cov, width = 10, height = 6)

# Mapping percentage plot
mapped_png <- file.path(output_dir, "mapped_percent_barplot.png")
mapping_plot <- coverage[!is.na(coverage$Mapped_Percent), , drop = FALSE]
if (nrow(mapping_plot) == 0) {
  stop("No valid mapped_percent values found in input file.")
}
p_map <- ggplot(mapping_plot, aes(x = reorder(Sample, -Mapped_Percent), y = Mapped_Percent)) +
  geom_col(fill = "darkseagreen4") +
  geom_text(aes(label = round(Mapped_Percent, 1)), vjust = -0.5, size = 3) +
  geom_hline(yintercept = 90, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Mapped Reads Percentage per Sample",
       x = "Sample",
       y = "Mapped reads (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = mapped_png, plot = p_map, width = 10, height = 6)

# Histogram of average depth
depth_hist_png <- file.path(output_dir, "coverage_histogram.png")
p_hist <- ggplot(coverage_plot, aes(x = Average_Coverage)) +
  geom_histogram(binwidth = 2, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Average Coverage",
       x = "Average coverage (Q >= 20)",
       y = "Number of samples") +
  theme_minimal()
ggsave(filename = depth_hist_png, plot = p_hist, width = 8, height = 6)

# Scatter plot: mapped percentage vs average coverage
scatter_png <- file.path(output_dir, "mapped_vs_coverage_scatter.png")
p_scatter <- ggplot(mapping_plot, aes(x = Mapped_Percent, y = Average_Coverage, label = Sample)) +
  geom_point(size = 3, color = "darkorange3") +
  geom_text(vjust = -0.8, size = 3) +
  geom_vline(xintercept = 90, color = "red", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Mapped Percentage vs Average Coverage",
       x = "Mapped reads (%)",
       y = "Average coverage (Q >= 20)") +
  theme_minimal()
ggsave(filename = scatter_png, plot = p_scatter, width = 9, height = 6)

message("Coverage plot saved to:")
message("  - ", coverage_png)
message("Mapping plot saved to:")
message("  - ", mapped_png)
message("Coverage histogram saved to:")
message("  - ", depth_hist_png)
message("Mapped vs coverage scatter saved to:")
message("  - ", scatter_png)

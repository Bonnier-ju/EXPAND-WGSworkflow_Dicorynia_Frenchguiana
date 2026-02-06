# coverage_summary_plot.R
# Usage: Rscript coverage_summary_plot.R input_tsv output_dir

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
parent_output_dir <- args[2]

# Create R_results subfolder if it doesn't exist
output_dir <- file.path(parent_output_dir, "R_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load required libraries
library(ggplot2)
library(readr)

# Read input coverage file
coverage <- read_tsv(input_file, col_types = cols())

# Clean column names (if needed)
colnames(coverage) <- c("Sample", "Average_Coverage")

# Plot coverage per sample with threshold line
p <- ggplot(coverage, aes(x = reorder(Sample, -Average_Coverage), y = Average_Coverage)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(Average_Coverage, 1)), vjust = -0.5, size = 3) +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 1, y = 10, label = "10x threshold", vjust = -1.2, hjust = 0, color = "red", size = 3) +
  labs(title = "Average Coverage per Sample",
       x = "Sample",
       y = "Average Coverage (filtered Q ≥ 20)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
output_png <- file.path(output_dir, "coverage_barplot.png")
ggsave(filename = output_png, plot = p, width = 10, height = 6)

message("✅ Coverage plot saved to:")
message("  - ", output_png)

#!/usr/bin/env Rscript
# 08.0.1-chelsa_extraction.R
# Extract CHELSA BIO1-BIO19 values at site coordinates.
# Requires: terra (for raster extraction)
# Input:  geoloc_site.csv, chelsa_rasters/ directory
# Output: chelsa_env_per_site.csv, chelsa_env_scaled.csv, correlation plots

suppressPackageStartupMessages({
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required. Install with: install.packages('terra')")
  library(terra)
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
geoloc_file <- if (length(args) >= 1) args[1] else stop("geoloc_file required")
chelsa_dir  <- if (length(args) >= 2) args[2] else stop("chelsa_dir required")
out_dir     <- if (length(args) >= 3) args[3] else stop("out_dir required")

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# STEP 1: Load site coordinates
# -------------------------------------------------------------------
geoloc <- read.csv(geoloc_file, stringsAsFactors = FALSE)
colnames(geoloc) <- tolower(colnames(geoloc))  # normalize to lowercase

# Remove sites with missing coordinates
n_before <- nrow(geoloc)
geoloc   <- geoloc[!is.na(geoloc$lat) & !is.na(geoloc$long), ]
n_after  <- nrow(geoloc)
if (n_before > n_after) {
  cat(sprintf("INFO: %d site(s) with missing coordinates excluded: %s\n",
              n_before - n_after,
              paste(geoloc$sites[is.na(geoloc$lat) | is.na(geoloc$long)],
                    collapse = ", ")))
}
cat(sprintf("INFO: %d sites with valid coordinates\n", n_after))

# SpatVector for terra extraction
pts <- vect(geoloc, geom = c("long", "lat"), crs = "EPSG:4326")

# -------------------------------------------------------------------
# STEP 2: Load CHELSA rasters and extract values
# -------------------------------------------------------------------
bio_files <- file.path(chelsa_dir,
  sprintf("CHELSA_bio%d_1981-2010_V.2.1.tif", 1:19))

missing_files <- bio_files[!file.exists(bio_files)]
if (length(missing_files) > 0) {
  stop("Missing CHELSA raster files:\n", paste(missing_files, collapse = "\n"))
}

cat("Loading CHELSA rasters and extracting values...\n")
bio_stack <- rast(bio_files)
names(bio_stack) <- paste0("BIO", 1:19)

# Extract values at site coordinates (bilinear interpolation)
extracted <- extract(bio_stack, pts, method = "bilinear")
extracted$ID <- NULL  # remove terra row index

# CHELSA BIO variables are stored as integers with scale factors:
# BIO1, BIO2, BIO5, BIO6, BIO7, BIO8, BIO9, BIO10, BIO11 : temp * 10 (°C * 10)
# BIO3 : dimensionless (isothermality, %)
# BIO4 : temperature seasonality * 100
# BIO12-BIO19 : precipitation (mm)
# Apply CHELSA v2.1 scale factors to convert to standard units
temp_vars <- c("BIO1","BIO2","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11")
for (v in temp_vars) {
  extracted[[v]] <- extracted[[v]] / 10
}
# BIO4: divide by 100 to get standard CV units
extracted[["BIO4"]] <- extracted[["BIO4"]] / 100

env_df <- cbind(geoloc[, c("sites", "lat", "long")], extracted)
colnames(env_df)[1] <- "site"

# Save raw extraction
out_raw <- file.path(out_dir, "chelsa_env_per_site.csv")
write.csv(env_df, file = out_raw, row.names = FALSE, quote = FALSE)
cat("Raw CHELSA values written:", out_raw, "\n")

# -------------------------------------------------------------------
# STEP 3: Standardize (z-score) for use in LFMM/RDA
# -------------------------------------------------------------------
bio_cols   <- paste0("BIO", 1:19)
env_scaled <- env_df
env_scaled[bio_cols] <- scale(env_df[bio_cols])

out_scaled <- file.path(out_dir, "chelsa_env_scaled.csv")
write.csv(env_scaled, file = out_scaled, row.names = FALSE, quote = FALSE)
cat("Scaled CHELSA values written:", out_scaled, "\n")

# -------------------------------------------------------------------
# STEP 4: Correlation matrix among BIO variables
# -------------------------------------------------------------------
cor_mat <- cor(env_df[bio_cols], use = "pairwise.complete.obs")

# Long format for ggplot heatmap
cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Var1", "Var2", "Correlation")
cor_long$Var1 <- factor(cor_long$Var1, levels = paste0("BIO", 1:19))
cor_long$Var2 <- factor(cor_long$Var2, levels = paste0("BIO", 1:19))

p_cor <- ggplot(cor_long, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2.0, color = "black") +
  scale_fill_gradientn(
    colours = c("#d73027", "#ffffff", "#4575b4"),
    limits  = c(-1, 1),
    name    = "Pearson r"
  ) +
  labs(title = "CHELSA BIO1-BIO19 — Pearson correlation matrix", x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title  = element_text(face = "bold")
  )
ggsave(file.path(plot_dir, "chelsa_bio_correlation.png"), p_cor,
       width = 10, height = 9, dpi = 300)
ggsave(file.path(plot_dir, "chelsa_bio_correlation.pdf"), p_cor,
       width = 10, height = 9)

# -------------------------------------------------------------------
# STEP 5: Identify highly correlated variable pairs (|r| >= 0.8)
# -------------------------------------------------------------------
high_cor <- which(abs(cor_mat) >= 0.8 & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  hc_df <- data.frame(
    Var1 = rownames(cor_mat)[high_cor[, 1]],
    Var2 = colnames(cor_mat)[high_cor[, 2]],
    r    = round(cor_mat[high_cor], 3)
  )
  hc_df <- hc_df[order(-abs(hc_df$r)), ]
  out_hc <- file.path(out_dir, "chelsa_high_correlations.tsv")
  write.table(hc_df, file = out_hc, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("INFO: %d highly correlated pairs (|r| >= 0.8) written to: %s\n",
              nrow(hc_df), out_hc))
} else {
  cat("INFO: no highly correlated pairs (|r| >= 0.8) detected\n")
}

# -------------------------------------------------------------------
# STEP 6: Summary table
# -------------------------------------------------------------------
summary_df <- data.frame(
  variable    = bio_cols,
  mean        = round(colMeans(env_df[bio_cols], na.rm = TRUE), 3),
  sd          = round(apply(env_df[bio_cols], 2, sd, na.rm = TRUE), 3),
  min         = round(apply(env_df[bio_cols], 2, min, na.rm = TRUE), 3),
  max         = round(apply(env_df[bio_cols], 2, max, na.rm = TRUE), 3),
  n_sites     = colSums(!is.na(env_df[bio_cols])),
  stringsAsFactors = FALSE
)
out_summ <- file.path(out_dir, "chelsa_env_summary.tsv")
write.table(summary_df, file = out_summ, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary statistics written:", out_summ, "\n")

cat("DONE 08.0.1 CHELSA extraction completed\n")
cat("Outputs:\n")
cat("  Raw values  :", out_raw, "\n")
cat("  Scaled      :", out_scaled, "\n")
cat("  Correlation :", file.path(plot_dir, "chelsa_bio_correlation.png"), "\n")
cat("  High corr.  :", file.path(out_dir, "chelsa_high_correlations.tsv"), "\n")

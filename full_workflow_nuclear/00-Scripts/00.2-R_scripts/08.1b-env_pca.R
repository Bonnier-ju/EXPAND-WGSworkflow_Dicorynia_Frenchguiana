#!/usr/bin/env Rscript
# 08.1b-env_pca.R
# PCA on user-selected environmental variables at site level.
#
# Inputs:
#   chelsa_env_per_site.csv   — CHELSA BIO variables (temp in K after /10, precip raw)
#   envirem_env_per_site.csv  — ENVIREM variables (physical units)
#   manual_variables_per_site.csv — elevation + habitat dummies
#   variables_to_keep_template.txt — user selection (uncommented = kept)
#
# Outputs (plots/):
#   pca_screeplot.png/pdf
#   pca_scores_pc12.png/pdf  — sites on PC1 × PC2
#   pca_scores_pc13.png/pdf  — sites on PC1 × PC3
#   pca_biplot_pc12.png/pdf  — biplot: sites + variable loadings
#   pca_loadings_circle.png/pdf — correlation circle (variables only)
#
# Scale factor corrections applied here (chelsa_env_per_site.csv):
#   Absolute temperatures (BIO1,5,6,8,9,10,11): stored in K after /10 → -273.15
#   Precipitation (BIO12-BIO19): stored raw (×10 of mm) → /10
#   ENVIREM variables: already in physical units
#   Elevation: already in metres

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) library(ggrepel)

SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else {
  warning("Site color file not found: ", SITE_COLORS_FILE)
  c()
}

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args            <- commandArgs(trailingOnly = TRUE)
chelsa_csv      <- if (length(args) >= 1) args[1] else stop("chelsa_csv required")
envirem_csv     <- if (length(args) >= 2) args[2] else stop("envirem_csv required")
manual_csv      <- if (length(args) >= 3) args[3] else stop("manual_csv required")
vars_file       <- if (length(args) >= 4) args[4] else stop("vars_file required")
out_dir         <- if (length(args) >= 5) args[5] else stop("out_dir required")

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# STEP 1: Parse selected variables from template file
# -------------------------------------------------------------------
cat("Parsing selected variables from:", vars_file, "\n")
lines <- readLines(vars_file)
# Keep lines that are NOT pure comments (don't start with optional whitespace + #)
# and that contain a variable name after trimming
selected_vars <- character(0)
for (ln in lines) {
  stripped <- trimws(ln)
  # Skip blank or comment-only lines
  if (nchar(stripped) == 0 || grepl("^#", stripped)) next
  # Extract variable name: first token before any # or whitespace
  var_name <- sub("\\s.*", "", sub("#.*", "", stripped))
  var_name <- trimws(var_name)
  if (nchar(var_name) > 0) selected_vars <- c(selected_vars, var_name)
}
cat(sprintf("INFO: %d variables selected: %s\n", length(selected_vars),
            paste(selected_vars, collapse = ", ")))

# -------------------------------------------------------------------
# STEP 2: Load site-level data
# -------------------------------------------------------------------
cat("Loading CHELSA site data:", chelsa_csv, "\n")
chelsa <- read.csv(chelsa_csv, stringsAsFactors = FALSE, check.names = FALSE)

cat("Loading ENVIREM site data:", envirem_csv, "\n")
envirem <- read.csv(envirem_csv, stringsAsFactors = FALSE, check.names = FALSE)

cat("Loading manual site data:", manual_csv, "\n")
manual <- read.csv(manual_csv, stringsAsFactors = FALSE, check.names = FALSE)

# -------------------------------------------------------------------
# STEP 3: Apply CHELSA scale factor corrections
# The chelsa_env_per_site.csv stores:
#   Absolute temp vars (BIO1,5,6,8,9,10,11): Kelvin = raw_int16 / 10 → subtract 273.15
#   Precip vars (BIO12-BIO19): raw_int16 values → divide by 10 for mm
# -------------------------------------------------------------------
temp_abs_vars  <- intersect(c("BIO1","BIO5","BIO6","BIO8","BIO9","BIO10","BIO11"),
                            names(chelsa))
precip_vars    <- intersect(paste0("BIO", 12:19), names(chelsa))

for (v in temp_abs_vars)  chelsa[[v]] <- chelsa[[v]] - 273.15
for (v in precip_vars)    chelsa[[v]] <- chelsa[[v]] / 10

cat(sprintf("INFO: BIO5 range after correction: %.1f - %.1f °C\n",
            min(chelsa$BIO5, na.rm = TRUE), max(chelsa$BIO5, na.rm = TRUE)))
cat(sprintf("INFO: BIO12 range after correction: %.0f - %.0f mm\n",
            min(chelsa$BIO12, na.rm = TRUE), max(chelsa$BIO12, na.rm = TRUE)))

# -------------------------------------------------------------------
# STEP 4: Merge data sources and filter to French Guiana sites
# -------------------------------------------------------------------
env_all <- merge(chelsa[, c("site", names(chelsa)[names(chelsa) != "lat" & names(chelsa) != "long"])],
                 envirem[, c("site", names(envirem)[names(envirem) != "lat" & names(envirem) != "long"])],
                 by = "site", all.x = TRUE)
env_all <- merge(env_all,
                 manual[, c("site", "elevation")],
                 by = "site", all.x = TRUE)

# Exclude outgroup site
env_all <- env_all[env_all$site != "Cameroun_Benin", ]
cat(sprintf("INFO: %d French Guiana sites after filtering\n", nrow(env_all)))

# -------------------------------------------------------------------
# STEP 5: Select chosen variables and handle NAs
# -------------------------------------------------------------------
missing_vars <- setdiff(selected_vars, names(env_all))
if (length(missing_vars) > 0) {
  warning("Variables not found in merged data: ", paste(missing_vars, collapse = ", "))
  selected_vars <- intersect(selected_vars, names(env_all))
}

pca_data <- env_all[, c("site", selected_vars)]

# Report missingness
n_complete <- sum(complete.cases(pca_data[, selected_vars]))
cat(sprintf("INFO: %d / %d sites have complete data for all selected variables\n",
            n_complete, nrow(pca_data)))

# Keep only sites with complete data
pca_complete <- pca_data[complete.cases(pca_data[, selected_vars]), ]
cat(sprintf("INFO: Proceeding with %d sites\n", nrow(pca_complete)))

# -------------------------------------------------------------------
# STEP 6: Run PCA
# -------------------------------------------------------------------
cat("Running PCA (prcomp, scale = TRUE)...\n")
pca_mat <- as.matrix(pca_complete[, selected_vars])
rownames(pca_mat) <- pca_complete$site

pca_res <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

# Variance explained
var_exp   <- pca_res$sdev^2
prop_var  <- var_exp / sum(var_exp)
cum_var   <- cumsum(prop_var)
n_pc      <- min(length(selected_vars), nrow(pca_complete) - 1)

cat("Variance explained per PC:\n")
for (i in seq_len(min(n_pc, 6))) {
  cat(sprintf("  PC%d: %.1f%% (cumulative: %.1f%%)\n",
              i, prop_var[i] * 100, cum_var[i] * 100))
}

# Scores (site positions) and loadings (variable contributions)
scores   <- as.data.frame(pca_res$x)
scores$site <- rownames(scores)
loadings <- as.data.frame(pca_res$rotation)
loadings$variable <- rownames(loadings)

# Axis labels with variance explained
lab_pc <- function(i) sprintf("PC%d (%.1f%%)", i, prop_var[i] * 100)

# -------------------------------------------------------------------
# STEP 7: Scree plot
# -------------------------------------------------------------------
scree_df <- data.frame(
  PC      = factor(paste0("PC", seq_len(n_pc)), levels = paste0("PC", seq_len(n_pc))),
  var_pct = prop_var[seq_len(n_pc)] * 100,
  cum_pct = cum_var[seq_len(n_pc)] * 100
)

p_scree <- ggplot(scree_df, aes(x = PC, y = var_pct)) +
  geom_col(fill = "#4575b4", width = 0.65) +
  geom_line(aes(y = cum_pct, group = 1), color = "#d73027", linewidth = 0.8) +
  geom_point(aes(y = cum_pct), color = "#d73027", size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  labs(title = "PCA - Variance explained per axis",
       subtitle = "Bars: individual variance; red line: cumulative; dashed: 80%",
       x = NULL, y = "Variance explained (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "pca_screeplot.png"), p_scree,
       width = max(6, n_pc * 0.6), height = 4, dpi = 300)
ggsave(file.path(plot_dir, "pca_screeplot.pdf"), p_scree,
       width = max(6, n_pc * 0.6), height = 4)
cat("Scree plot written\n")

# -------------------------------------------------------------------
# Helper: build site score plot
# -------------------------------------------------------------------
make_score_plot <- function(scores, pc_x = 1, pc_y = 2, site_colors = c()) {
  xcol <- paste0("PC", pc_x)
  ycol <- paste0("PC", pc_y)

  p <- ggplot(scores, aes(x = .data[[xcol]], y = .data[[ycol]],
                           color = site, label = site)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_point(size = 3.5) +
    labs(title = sprintf("PCA — Site scores (PC%d x PC%d)", pc_x, pc_y),
         x = lab_pc(pc_x), y = lab_pc(pc_y),
         color = "Site") +
    theme(legend.position = "right",
          legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 7))

  if (length(site_colors) > 0)
    p <- p + scale_color_manual(values = site_colors, na.value = "grey50")

  if (has_ggrepel) {
    p <- p + ggrepel::geom_text_repel(size = 2.8, show.legend = FALSE,
                                      max.overlaps = 20, seed = 42)
  } else {
    p <- p + geom_text(vjust = -0.7, size = 2.8, show.legend = FALSE)
  }
  p
}

# -------------------------------------------------------------------
# STEP 8: Site score plots (PC1×PC2 and PC1×PC3)
# -------------------------------------------------------------------
p_sc12 <- make_score_plot(scores, 1, 2, site_colors)
ggsave(file.path(plot_dir, "pca_scores_pc12.png"), p_sc12,
       width = 9, height = 6.5, dpi = 300)
ggsave(file.path(plot_dir, "pca_scores_pc12.pdf"), p_sc12,
       width = 9, height = 6.5)
cat("Score plot PC1×PC2 written\n")

if (n_pc >= 3) {
  p_sc13 <- make_score_plot(scores, 1, 3, site_colors)
  ggsave(file.path(plot_dir, "pca_scores_pc13.png"), p_sc13,
         width = 9, height = 6.5, dpi = 300)
  ggsave(file.path(plot_dir, "pca_scores_pc13.pdf"), p_sc13,
         width = 9, height = 6.5)
  cat("Score plot PC1×PC3 written\n")
}

# -------------------------------------------------------------------
# STEP 9: Variable correlation circle (loadings × sqrt(eigenvalue))
# Correlation between variable and PC: loading × sqrt(eigenvalue)
# -------------------------------------------------------------------
# Scale loadings to correlation with PCs
eig     <- pca_res$sdev^2
cor_mat <- sweep(pca_res$rotation, 2, pca_res$sdev, `*`)
cor_df  <- as.data.frame(cor_mat)
cor_df$variable <- rownames(cor_mat)

make_circle_plot <- function(cor_df, pc_x = 1, pc_y = 2) {
  xcol <- paste0("PC", pc_x)
  ycol <- paste0("PC", pc_y)

  # Unit circle for reference
  theta  <- seq(0, 2 * pi, length.out = 200)
  circle <- data.frame(x = cos(theta), y = sin(theta))

  p <- ggplot(cor_df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_path(data = circle, aes(x = x, y = y),
              color = "grey60", linewidth = 0.4, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
    geom_segment(aes(x = 0, y = 0,
                     xend = .data[[xcol]], yend = .data[[ycol]]),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 color = "#4575b4", linewidth = 0.7) +
    coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
    labs(title = sprintf("PCA — Correlation circle (PC%d x PC%d)", pc_x, pc_y),
         x = lab_pc(pc_x), y = lab_pc(pc_y))

  if (has_ggrepel) {
    p <- p + ggrepel::geom_text_repel(aes(label = variable),
                                      size = 2.8, color = "#d73027",
                                      max.overlaps = 30, seed = 42)
  } else {
    p <- p + geom_text(aes(label = variable),
                       vjust = -0.5, size = 2.8, color = "#d73027")
  }
  p
}

p_circ12 <- make_circle_plot(cor_df, 1, 2)
ggsave(file.path(plot_dir, "pca_loadings_circle_pc12.png"), p_circ12,
       width = 7, height = 7, dpi = 300)
ggsave(file.path(plot_dir, "pca_loadings_circle_pc12.pdf"), p_circ12,
       width = 7, height = 7)
cat("Correlation circle PC1×PC2 written\n")

if (n_pc >= 3) {
  p_circ13 <- make_circle_plot(cor_df, 1, 3)
  ggsave(file.path(plot_dir, "pca_loadings_circle_pc13.png"), p_circ13,
         width = 7, height = 7, dpi = 300)
  ggsave(file.path(plot_dir, "pca_loadings_circle_pc13.pdf"), p_circ13,
         width = 7, height = 7)
  cat("Correlation circle PC1×PC3 written\n")
}

# -------------------------------------------------------------------
# STEP 10: Combined biplot (sites + variable arrows, PC1×PC2)
# Scale factor: fit loadings to same range as scores
# -------------------------------------------------------------------
score_range <- max(abs(c(scores$PC1, scores$PC2)))
load_range  <- max(abs(c(loadings$PC1, loadings$PC2)))
scale_fac   <- score_range / load_range * 0.75

load_scaled       <- loadings
load_scaled$PC1   <- loadings$PC1 * scale_fac
load_scaled$PC2   <- loadings$PC2 * scale_fac

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_segment(data = load_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "grey40", linewidth = 0.6) +
  geom_point(data = scores,
             aes(x = PC1, y = PC2, color = site),
             size = 3.5) +
  labs(title = "PCA — Biplot (PC1 x PC2)",
       x = lab_pc(1), y = lab_pc(2), color = "Site") +
  theme(legend.position = "right",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 7))

if (length(site_colors) > 0)
  p_biplot <- p_biplot +
    scale_color_manual(values = site_colors, na.value = "grey50")

if (has_ggrepel) {
  p_biplot <- p_biplot +
    ggrepel::geom_text_repel(data = scores,
                             aes(x = PC1, y = PC2, label = site, color = site),
                             size = 2.6, show.legend = FALSE,
                             max.overlaps = 20, seed = 42) +
    ggrepel::geom_text_repel(data = load_scaled,
                             aes(x = PC1, y = PC2, label = variable),
                             size = 2.8, color = "grey30",
                             fontface = "italic",
                             max.overlaps = 30, seed = 42)
} else {
  p_biplot <- p_biplot +
    geom_text(data = scores,
              aes(x = PC1, y = PC2, label = site, color = site),
              vjust = -0.7, size = 2.6, show.legend = FALSE) +
    geom_text(data = load_scaled,
              aes(x = PC1, y = PC2, label = variable),
              vjust = -0.5, size = 2.8, color = "grey30", fontface = "italic")
}

ggsave(file.path(plot_dir, "pca_biplot_pc12.png"), p_biplot,
       width = 10, height = 7, dpi = 300)
ggsave(file.path(plot_dir, "pca_biplot_pc12.pdf"), p_biplot,
       width = 10, height = 7)
cat("Biplot PC1×PC2 written\n")

# -------------------------------------------------------------------
# STEP 11: Save PCA scores and loadings as TSVs
# -------------------------------------------------------------------
pca_out_dir <- file.path(out_dir, "pca")
dir.create(pca_out_dir, recursive = TRUE, showWarnings = FALSE)

write.table(scores, file = file.path(pca_out_dir, "pca_site_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(loadings, file = file.path(pca_out_dir, "pca_variable_loadings.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(PC = paste0("PC", seq_len(n_pc)),
             eigenvalue = var_exp[seq_len(n_pc)],
             prop_var   = round(prop_var[seq_len(n_pc)] * 100, 2),
             cum_var    = round(cum_var[seq_len(n_pc)] * 100, 2)),
  file = file.path(pca_out_dir, "pca_variance_explained.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("PCA tables written to:", pca_out_dir, "\n")

cat("\nDONE 08.1b PCA on selected environmental variables completed\n")
cat("Plots:\n")
for (fn in list.files(plot_dir, pattern = "^pca_", full.names = TRUE)) {
  cat("  ", fn, "\n")
}

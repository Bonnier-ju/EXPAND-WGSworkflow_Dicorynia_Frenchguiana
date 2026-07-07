#!/usr/bin/env Rscript
# 07b.0-purity_filtering.R
# Filter individuals by ADMIXTURE ancestry proportion purity at K=3.
#
# For each threshold T in {0.70, 0.80, 0.90, 0.95}: an individual is
# considered "pure" for a cluster when max(Q1,Q2,Q3) >= T.
# Produces one combined PLINK --keep file per threshold (all pure individuals,
# no group label) for downstream PCA and ADMIXTURE re-runs (07b.1, 07b.2).
# Per-cluster keep files are produced by 07b.2 from new ADMIXTURE results.
#
# Inputs  (9 positional args):
#   1  q_file          admixture_input.3.named.Q  (space-sep, no header: IND_ID Q1 Q2 Q3)
#   2  fam_file        admixture_input.fam         (PLINK FAM)
#   3  meta_file       sample_sheet_complete_filtered.csv
#   4  geoloc_file     geoloc_site.csv             (Sites, lat, long)
#   5  sites_clusters  sites_by_clusters.csv       (Sites, K=2, K=3, K=4)
#   6  sites_colors    sites_couleurs.csv          (site, couleur_hex)
#   7  out_lists       output dir for keep files
#   8  out_summary     output dir for summary TSV
#   9  out_plots       output dir for plots

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

POP_COLORS <- c(Pop_1 = "#EE7600", Pop_2 = "#458B00", Pop_3 = "#CD2626")
THRESHOLDS <- c(0.70, 0.80, 0.90, 0.95)

# -------------------------------------------------------------------
# STEP 1: Parse arguments
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: 07b.0-purity_filtering.R q_file fam_file meta_file geoloc_file sites_clusters sites_colors out_lists out_summary out_plots")
}
q_file         <- args[1]
fam_file       <- args[2]
meta_file      <- args[3]
geoloc_file    <- args[4]
sites_clusters <- args[5]
sites_colors   <- args[6]
out_lists      <- args[7]
out_summary    <- args[8]
out_plots      <- args[9]

for (d in c(out_lists, out_summary, out_plots))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# STEP 2: Load input files
# -------------------------------------------------------------------
cat("Loading input files...\n")

# Q file: no header, space-separated
q_raw <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
if (ncol(q_raw) != 4)
  stop("Q file must have 4 columns: IND_ID Q1 Q2 Q3")
colnames(q_raw) <- c("IND_ID", "Q1", "Q2", "Q3")
cat(sprintf("  Q file: %d individuals\n", nrow(q_raw)))

# FAM file: FID IID 0 0 sex pheno
fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
colnames(fam)[1:2] <- c("FID", "IID")
cat(sprintf("  FAM file: %d individuals\n", nrow(fam)))

# Sample metadata: sample_id, site (+ other columns)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
cat(sprintf("  Sample metadata: %d rows, columns: %s\n",
            nrow(meta), paste(names(meta), collapse = ", ")))

# Geoloc: Sites, lat, long (capital S)
geoloc <- read.csv(geoloc_file, stringsAsFactors = FALSE)
names(geoloc)[1] <- "site"  # normalize to lowercase
cat(sprintf("  Geoloc: %d sites\n", nrow(geoloc)))

# Sites × clusters: Sites, K=2, K=3, K=4
sc <- read.csv(sites_clusters, stringsAsFactors = FALSE, check.names = FALSE)
names(sc)[1] <- "site"
cat(sprintf("  Sites-clusters: %d sites\n", nrow(sc)))

# Site colors
sc_colors <- read.csv(sites_colors, stringsAsFactors = FALSE)
site_colors <- setNames(sc_colors$couleur_hex, sc_colors$site)
cat(sprintf("  Site colors: %d sites\n", length(site_colors)))

# -------------------------------------------------------------------
# STEP 3: Merge Q with individual metadata
# -------------------------------------------------------------------
cat("Merging Q with metadata...\n")

# Merge Q with site info (sample_id = IND_ID in Q file)
merged <- merge(q_raw,
                meta[, c("sample_id", "site")],
                by.x = "IND_ID", by.y = "sample_id",
                all.x = FALSE)
cat(sprintf("  After merge with metadata: %d individuals\n", nrow(merged)))

# Attach longitude for west→east ordering
merged <- merge(merged, geoloc[, c("site", "long")],
                by = "site", all.x = TRUE)

# Keep only French Guiana sites (those in sites_by_clusters)
fg_sites <- unique(sc$site)
merged_fg <- merged[merged$site %in% fg_sites, ]
cat(sprintf("  After filtering to FG sites: %d individuals (%d sites)\n",
            nrow(merged_fg), length(unique(merged_fg$site))))

# Exclude Angela/treemutation individuals if present (defensive check —
# guiana_only Q file already excludes them, but kept for robustness)
EXCLUDE_IDS <- c("DBR_AGJ", "DBR_AGK", "DBR_AGL", "DBR_AGM",
                  "DBR_ADT")  # sampling error
excl_found  <- intersect(merged_fg$IND_ID, EXCLUDE_IDS)
if (length(excl_found) > 0) {
  cat(sprintf("  Excluding Angela/treemutation individuals: %s\n",
              paste(excl_found, collapse = ", ")))
  merged_fg <- merged_fg[!merged_fg$IND_ID %in% EXCLUDE_IDS, ]
  cat(sprintf("  After exclusion: %d individuals (%d sites)\n",
              nrow(merged_fg), length(unique(merged_fg$site))))
}

if (nrow(merged_fg) == 0)
  stop("No French Guiana individuals found after filtering — check site name matching.")

# -------------------------------------------------------------------
# STEP 4: Determine Q-column → cluster name mapping
# -------------------------------------------------------------------
# Uses individual-level dominant Q (argmax per individual) rather than
# site-level means. This is robust to admixed sites where the site-level
# mean may be low even though the cluster's column is correct.
# -------------------------------------------------------------------
cat("Determining Q-column to cluster mapping (individual-dominant approach)...\n")

k3_col <- "K=3"
sc_fg  <- sc[sc$site %in% fg_sites, c("site", k3_col)]
names(sc_fg)[2] <- "cluster_k3"

# For each Q column, find individuals for which that column is the highest,
# then determine the modal site-cluster assignment for those individuals
q_mat_all    <- as.matrix(merged_fg[, c("Q1", "Q2", "Q3")])
dominant_idx <- apply(q_mat_all, 1, which.max)

col_to_cluster <- character(3)
names(col_to_cluster) <- c("Q1", "Q2", "Q3")

for (i in 1:3) {
  qcol     <- c("Q1", "Q2", "Q3")[i]
  dom_inds <- merged_fg[dominant_idx == i, , drop = FALSE]
  n_dom    <- nrow(dom_inds)

  if (n_dom == 0) {
    warning(sprintf("No individuals with %s as dominant column — cannot map.", qcol))
    next
  }

  dom_sites <- merge(dom_inds[, "site", drop = FALSE],
                     sc_fg, by = "site", all.x = TRUE)
  dom_sites <- dom_sites[!is.na(dom_sites$cluster_k3), ]

  if (nrow(dom_sites) == 0) {
    warning(sprintf("Cannot determine cluster for %s — no site match in sites_by_clusters.", qcol))
    next
  }

  modal_cluster <- names(sort(table(dom_sites$cluster_k3), decreasing = TRUE))[1]
  mean_q        <- mean(dom_inds[[qcol]])
  col_to_cluster[qcol] <- modal_cluster
  cat(sprintf("  %s → %s  (n_dominant=%d, mean Q=%.3f)\n",
              qcol, modal_cluster, n_dom, mean_q))
}

# Conflict check: each cluster must be assigned to exactly one Q column
assigned_clusters <- col_to_cluster[nchar(col_to_cluster) > 0]
if (any(duplicated(assigned_clusters))) {
  stop(sprintf(
    "Mapping conflict: two Q columns assigned to the same cluster (%s).\n  Check sites_by_clusters.csv for K=3 against the Q file.",
    paste(assigned_clusters[duplicated(assigned_clusters)], collapse = ", ")
  ))
}

cluster_to_col <- setNames(names(col_to_cluster), col_to_cluster)
cat("  Final mapping:",
    paste(sprintf("%s=%s", names(col_to_cluster), col_to_cluster), collapse = " | "), "\n")
cluster_names <- sort(unique(col_to_cluster[nchar(col_to_cluster) > 0]))

# -------------------------------------------------------------------
# STEP 5: Compute purity and cluster assignment per individual
# -------------------------------------------------------------------
cat("Computing individual purity scores...\n")

q_vals        <- as.matrix(merged_fg[, c("Q1", "Q2", "Q3")])
purity        <- apply(q_vals, 1, max)
best_col_idx  <- apply(q_vals, 1, which.max)
best_col_name <- c("Q1", "Q2", "Q3")[best_col_idx]
assigned_cluster <- unname(col_to_cluster[best_col_name])

merged_fg$purity           <- purity
merged_fg$assigned_cluster <- assigned_cluster

# West→east site ordering by longitude
site_lon_order <- geoloc[geoloc$site %in% fg_sites, ]
site_lon_order <- site_lon_order[order(site_lon_order$long), ]
site_order     <- site_lon_order$site

merged_fg$site <- factor(merged_fg$site, levels = site_order)
merged_fg      <- merged_fg[order(merged_fg$site, merged_fg$IND_ID), ]

cat(sprintf("  Purity range: [%.3f, %.3f]  median: %.3f\n",
            min(purity), max(purity), median(purity)))

# -------------------------------------------------------------------
# STEP 6: Per-threshold filtering and keep file output
# -------------------------------------------------------------------
# Only COMBINED keep files are produced here (all pure individuals,
# no group label). These feed directly into:
#   - PCA re-run (07b.1)       → structure discovered de novo
#   - ADMIXTURE re-run (07b.2) → NEW Q values + NEW cluster assignments
#
# Per-cluster keep files for diversity (07b.3), FST (07b.4), and Ne
# (07b.5) are generated by 07b.2 from the NEW ADMIXTURE results, not
# from the original K=3 assignments used here for purity filtering.
# -------------------------------------------------------------------
cat("\nSTEP 6: Per-threshold filtering...\n")

summary_rows <- list()

for (T in THRESHOLDS) {
  tag <- sprintf("T%02d", as.integer(T * 100))
  cat(sprintf("  Threshold %s (%.0f%%):\n", tag, T * 100))

  kept <- merged_fg[merged_fg$purity >= T, ]
  n_total <- nrow(kept)

  # Combined keep file (FID IID) — used for PCA and ADMIXTURE re-run
  fam_all <- merge(
    data.frame(IND_ID = kept$IND_ID, stringsAsFactors = FALSE),
    fam[, 1:2], by.x = "IND_ID", by.y = "IID", all.x = TRUE
  )
  fam_all$FID[is.na(fam_all$FID)] <- 0

  out_all <- file.path(out_lists, sprintf("inds_%s_all.txt", tag))
  write.table(fam_all[, c("FID", "IND_ID")], out_all,
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  cat(sprintf("    -> %s: %d individuals total (no group label)\n",
              basename(out_all), n_total))

  # Per-cluster counts for summary (informational only — based on original Q)
  for (cl in cluster_names) {
    kept_cl <- kept[kept$assigned_cluster == cl, ]
    n <- nrow(kept_cl)
    flag <- if (n < 15) " [!n<15]" else if (n < 20) " [warn n<20]" else ""
    cat(sprintf("       %s: %d individuals%s\n", cl, n, flag))

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      threshold      = T,
      threshold_tag  = tag,
      cluster        = cl,
      n_individuals  = n,
      flag_n_lt_15   = n < 15,
      flag_n_lt_20   = n < 20,
      stringsAsFactors = FALSE
    )
  }
}

# -------------------------------------------------------------------
# STEP 7: Write summary table
# -------------------------------------------------------------------
cat("\nSTEP 7: Writing summary...\n")

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df,
            file.path(out_summary, "purity_filtering_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== PURITY FILTERING SUMMARY ===\n")
print(summary_df[, c("threshold_tag", "cluster", "n_individuals",
                      "flag_n_lt_15", "flag_n_lt_20")])

# -------------------------------------------------------------------
# STEP 8: Purity distribution plot (faceted by cluster, per site)
# -------------------------------------------------------------------
cat("\nSTEP 8: Purity distribution plot...\n")

plot_df <- merged_fg[, c("IND_ID", "site", "purity", "assigned_cluster")]

p_dist <- ggplot(plot_df,
                 aes(x = site, y = purity, fill = assigned_cluster)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.6) +
  geom_jitter(aes(color = assigned_cluster),
              width = 0.15, size = 0.7, alpha = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = THRESHOLDS[1], linetype = "solid",   color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = THRESHOLDS[2], linetype = "dashed",  color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = THRESHOLDS[3], linetype = "dotdash", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = THRESHOLDS[4], linetype = "dotted",  color = "grey30", linewidth = 0.5) +
  geom_text(
    data = data.frame(
      y    = THRESHOLDS + 0.013,
      label = paste0(THRESHOLDS * 100, "%"),
      assigned_cluster = sort(cluster_names)[1]   # draw labels in first facet only
    ),
    aes(x = 1, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0, size = 2.8, color = "grey30"
  ) +
  scale_fill_manual(values  = POP_COLORS, name = "Cluster") +
  scale_color_manual(values = POP_COLORS) +
  scale_y_continuous(limits = c(0.3, 1.02), breaks = seq(0.3, 1.0, 0.1)) +
  facet_grid(assigned_cluster ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Individual purity scores by site and cluster (K=3)",
    subtitle = "Purity = max(Q1, Q2, Q3) — horizontal lines: T70/T80/T90/T95",
    x = "Site (west → east)",
    y = "Purity"
  ) +
  theme(
    axis.text.x  = element_text(angle = 50, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    strip.text       = element_text(size = 9, face = "bold")
  )

ggsave(file.path(out_plots, "purity_distribution_by_site.png"),
       p_dist, width = 14, height = 8, dpi = 150)
ggsave(file.path(out_plots, "purity_distribution_by_site.pdf"),
       p_dist, width = 14, height = 8)
cat("  Saved: purity_distribution_by_site.png\n")

# -------------------------------------------------------------------
# STEP 9: ADMIXTURE barplot — full FG dataset, west→east
# -------------------------------------------------------------------
cat("STEP 9: ADMIXTURE barplot (full dataset)...\n")

# Build a wide data frame with columns renamed Q1/Q2/Q3 → Pop_X cluster names
long_df <- merged_fg[, c("IND_ID", "site")]
for (qcol in c("Q1", "Q2", "Q3")) {
  cl_name <- col_to_cluster[qcol]
  long_df[[cl_name]] <- merged_fg[[qcol]]
}
long_df$ind_order <- seq_len(nrow(long_df))

long_melt <- pivot_longer(long_df,
                           cols      = all_of(cluster_names),
                           names_to  = "cluster",
                           values_to = "Q")
long_melt$cluster <- factor(long_melt$cluster, levels = cluster_names)

# Site midpoints for x-axis labels
site_mid <- long_df %>%
  group_by(site) %>%
  summarise(mid = mean(ind_order), .groups = "drop")
site_mid <- site_mid[match(levels(merged_fg$site), site_mid$site), ]
site_mid <- site_mid[!is.na(site_mid$site), ]

p_bar <- ggplot(long_melt, aes(x = ind_order, y = Q, fill = cluster)) +
  geom_col(width = 1, position = "stack") +
  scale_fill_manual(values = POP_COLORS, name = "Cluster") +
  scale_x_continuous(
    breaks = site_mid$mid,
    labels = site_mid$site,
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "ADMIXTURE K=3 — French Guiana individuals (west → east)",
    x = "Site",
    y = "Ancestry proportion"
  ) +
  theme(
    axis.text.x  = element_text(angle = 50, hjust = 1, size = 7),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(out_plots, "admixture_barplot_full.png"),
       p_bar, width = 16, height = 4, dpi = 150)
ggsave(file.path(out_plots, "admixture_barplot_full.pdf"),
       p_bar, width = 16, height = 4)
cat("  Saved: admixture_barplot_full.png\n")

# -------------------------------------------------------------------
# STEP 10: ADMIXTURE barplot per threshold (kept individuals only)
# -------------------------------------------------------------------
cat("STEP 10: ADMIXTURE barplots per threshold...\n")

for (T in THRESHOLDS) {
  tag  <- sprintf("T%02d", as.integer(T * 100))
  kept <- merged_fg[merged_fg$purity >= T, ]

  if (nrow(kept) == 0) {
    cat(sprintf("  %s: 0 individuals, skipping barplot\n", tag))
    next
  }

  kept$site     <- droplevels(kept$site)   # remove absent site levels → no NA breaks
  kept$ind_order <- seq_len(nrow(kept))

  # Rename Q columns to cluster names for long format
  kept_long_df <- kept[, c("IND_ID", "site", "ind_order")]
  for (qcol in c("Q1", "Q2", "Q3")) {
    cl_name <- col_to_cluster[qcol]
    kept_long_df[[cl_name]] <- kept[[qcol]]
  }

  long_k <- pivot_longer(kept_long_df,
                          cols      = all_of(cluster_names),
                          names_to  = "cluster",
                          values_to = "Q")
  long_k$cluster <- factor(long_k$cluster, levels = cluster_names)

  site_mid_k <- kept %>%
    group_by(site) %>%
    summarise(mid = mean(ind_order), .groups = "drop")

  p_k <- ggplot(long_k, aes(x = ind_order, y = Q, fill = cluster)) +
    geom_col(width = 1, position = "stack") +
    scale_fill_manual(values = POP_COLORS, name = "Cluster") +
    scale_x_continuous(
      breaks = site_mid_k$mid,
      labels = site_mid_k$site,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = sprintf("ADMIXTURE K=3 — pure individuals (purity ≥ %.0f%%) — %d individuals",
                      T * 100, nrow(kept)),
      x = "Site (west → east)",
      y = "Ancestry proportion"
    ) +
    theme(
      axis.text.x  = element_text(angle = 50, hjust = 1, size = 7),
      panel.grid   = element_blank(),
      legend.position = "right"
    )

  ggsave(file.path(out_plots, sprintf("admixture_barplot_%s.png", tag)),
         p_k, width = 16, height = 4, dpi = 150)
  cat(sprintf("  Saved: admixture_barplot_%s.png (%d individuals)\n", tag, nrow(kept)))
}

cat("\n=== 07b.0 purity filtering COMPLETE ===\n")
cat(sprintf("Keep files:   %s/\n", out_lists))
cat(sprintf("Summary:      %s/purity_filtering_summary.tsv\n", out_summary))
cat(sprintf("Plots:        %s/\n", out_plots))

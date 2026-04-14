#!/usr/bin/env Rscript
# 07.4-ibd_mantel.R
# Isolation by Distance (IBD) analysis — guiana_only group
#
# Outputs:
#   1. Mantel test results table
#   2. IBD scatter plot: FST/(1-FST) vs ln(geographic distance)
#   3. Heatmap of pairwise FST matrix
#   4. Heatmap of pairwise geographic distances (km)
#   5. Mantel statistic permutation distribution plot

suppressPackageStartupMessages({
  library(ggplot2)
  library(vegan)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# Site color palette
SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else c()

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
out_dir     <- if (length(args) >= 1) args[1] else stop("out_dir required")
map_file    <- if (length(args) >= 2) args[2] else stop("map_file required")
geoloc_file <- if (length(args) >= 3) args[3] else stop("geoloc_file required")
n_perm      <- if (length(args) >= 4) as.integer(args[4]) else 999L

group    <- "guiana_only"
fst_file <- file.path(out_dir, "fst_vcftools", group, "pairwise_fst_summary.tsv")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------
pop_map <- read.table(map_file, header = FALSE, sep = "\t",
                      col.names = c("sample_id", "population"),
                      stringsAsFactors = FALSE)

geoloc <- read.csv(geoloc_file, stringsAsFactors = FALSE)
colnames(geoloc) <- tolower(colnames(geoloc))
if (!"sites" %in% colnames(geoloc) && "site" %in% colnames(geoloc))
  colnames(geoloc)[colnames(geoloc) == "site"] <- "sites"

if (!file.exists(fst_file)) stop("FST summary file not found: ", fst_file)
fst_df <- read.table(fst_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# If summary file has NA FST values, recompute from .weir.fst files directly
weir_dir <- file.path(out_dir, "fst_vcftools", group)
n_na <- sum(fst_df$weighted_fst == "NA" | is.na(fst_df$weighted_fst))
if (n_na > 0) {
  cat(sprintf("INFO: %d NA FST values in summary — reading from .weir.fst files\n", n_na))
  for (k in seq_len(nrow(fst_df))) {
    if (!is.na(fst_df$weighted_fst[k]) && fst_df$weighted_fst[k] != "NA") next
    s1 <- gsub("[^A-Za-z0-9._-]", "_", fst_df$pop1[k])
    s2 <- gsub("[^A-Za-z0-9._-]", "_", fst_df$pop2[k])
    wf <- file.path(weir_dir, sprintf("%s_vs_%s.weir.fst", s1, s2))
    if (!file.exists(wf)) next
    wdat <- tryCatch(
      read.table(wf, header = TRUE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(wdat) || !"WEIR_AND_COCKERHAM_FST" %in% colnames(wdat)) next
    vals <- suppressWarnings(as.numeric(wdat$WEIR_AND_COCKERHAM_FST))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) next
    fst_df$mean_fst[k]     <- round(mean(vals), 6)
    fst_df$weighted_fst[k] <- round(mean(vals), 6)
  }
  n_still_na <- sum(fst_df$weighted_fst == "NA" | is.na(fst_df$weighted_fst))
  cat(sprintf("INFO: %d NA FST values remaining after .weir.fst fallback\n", n_still_na))
}

# -------------------------------------------------------------------
# Outgroups excluded from guiana_only (same as SLURM)
# -------------------------------------------------------------------
OUTGROUP_IDS     <- c("DBR_AGN", "DBR_AGO", "DBR_AGP", "DBR_AGQ", "DBR_AGR")
TREEMUTATION_IDS <- c("DBR_AGJ", "DBR_AGK", "DBR_AGL", "DBR_AGM")
excluded         <- c(OUTGROUP_IDS, TREEMUTATION_IDS)

guiana_map <- pop_map[!pop_map$sample_id %in% excluded, ]
pops_all   <- sort(unique(guiana_map$population))
pops       <- pops_all[pops_all %in% geoloc$sites]

missing_geo <- setdiff(pops_all, geoloc$sites)
if (length(missing_geo) > 0)
  cat(sprintf("WARN: populations without coordinates excluded: %s\n",
              paste(missing_geo, collapse = ", ")))

cat(sprintf("Populations included: %d\n", length(pops)))

# -------------------------------------------------------------------
# Haversine distance (km)
# -------------------------------------------------------------------
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R     <- 6371
  dlat  <- (lat2 - lat1) * pi / 180
  dlon  <- (lon2 - lon1) * pi / 180
  a     <- sin(dlat / 2)^2 +
           cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dlon / 2)^2
  2 * R * asin(sqrt(a))
}

# -------------------------------------------------------------------
# Build geographic distance matrix
# -------------------------------------------------------------------
geo_sub <- geoloc[match(pops, geoloc$sites), ]
n       <- length(pops)
geo_mat <- matrix(0, n, n, dimnames = list(pops, pops))
for (i in seq_len(n))
  for (j in seq_len(n))
    if (i != j)
      geo_mat[i, j] <- haversine_km(geo_sub$lat[i], geo_sub$long[i],
                                     geo_sub$lat[j], geo_sub$long[j])

# -------------------------------------------------------------------
# Build FST matrix
# -------------------------------------------------------------------
fst_mat <- matrix(NA_real_, n, n, dimnames = list(pops, pops))
diag(fst_mat) <- 0

for (k in seq_len(nrow(fst_df))) {
  p1  <- fst_df$pop1[k]; p2 <- fst_df$pop2[k]
  fst <- suppressWarnings(as.numeric(fst_df$weighted_fst[k]))
  if (is.na(fst)) fst <- suppressWarnings(as.numeric(fst_df$mean_fst[k]))
  if (!is.na(fst) && p1 %in% pops && p2 %in% pops) {
    fst <- max(0, fst)
    fst_mat[p1, p2] <- fst
    fst_mat[p2, p1] <- fst
  }
}

# Remove populations with missing FST
complete_pops <- pops[!apply(fst_mat, 1, function(x) any(is.na(x)))]
if (length(complete_pops) < 3) stop("Too few populations with complete FST data")

fst_mat <- fst_mat[complete_pops, complete_pops]
geo_mat <- geo_mat[complete_pops, complete_pops]
n_pops  <- length(complete_pops)
cat(sprintf("Populations with complete FST: %d\n", n_pops))

# -------------------------------------------------------------------
# Bootstrap significance of pairwise FST
# -------------------------------------------------------------------
# For each pair, read per-site FST values from the .weir.fst file,
# resample sites with replacement (nboot=1000), compute bootstrap mean.
# p-value = proportion of bootstrap means <= 0 (one-sided: H0: FST <= 0).
# Also report 95% bootstrap CI (lower, upper).
# -------------------------------------------------------------------
NBOOT <- 1000L
cat(sprintf("Bootstrap pairwise FST significance (nboot=%d)...\n", NBOOT))

pval_mat  <- matrix(NA_real_, n_pops, n_pops, dimnames = list(complete_pops, complete_pops))
lower_mat <- matrix(NA_real_, n_pops, n_pops, dimnames = list(complete_pops, complete_pops))
upper_mat <- matrix(NA_real_, n_pops, n_pops, dimnames = list(complete_pops, complete_pops))
diag(pval_mat) <- diag(lower_mat) <- diag(upper_mat) <- 0

set.seed(42L)
for (i in seq_len(n_pops)) {
  for (j in seq_len(n_pops)) {
    if (j <= i) next
    p1 <- complete_pops[i]; p2 <- complete_pops[j]
    s1 <- gsub("[^A-Za-z0-9._-]", "_", p1)
    s2 <- gsub("[^A-Za-z0-9._-]", "_", p2)
    wf <- file.path(weir_dir, sprintf("%s_vs_%s.weir.fst", s1, s2))
    if (!file.exists(wf)) next
    wdat <- tryCatch(read.table(wf, header = TRUE, stringsAsFactors = FALSE),
                     error = function(e) NULL)
    if (is.null(wdat) || !"WEIR_AND_COCKERHAM_FST" %in% colnames(wdat)) next
    vals <- suppressWarnings(as.numeric(wdat$WEIR_AND_COCKERHAM_FST))
    vals <- vals[!is.na(vals)]
    if (length(vals) < 2) next

    # Bootstrap: resample sites with replacement
    boot_means <- replicate(NBOOT, mean(sample(vals, length(vals), replace = TRUE)))
    pval  <- mean(boot_means <= 0)   # proportion of bootstraps with mean FST <= 0
    lower <- quantile(boot_means, 0.025)
    upper <- quantile(boot_means, 0.975)

    pval_mat[p1, p2]  <- pval_mat[p2, p1]  <- pval
    lower_mat[p1, p2] <- lower_mat[p2, p1] <- round(lower, 6)
    upper_mat[p1, p2] <- upper_mat[p2, p1] <- round(upper, 6)
  }
}

# --- Export 1: FST matrix (weighted FST, rounded) ---
fst_out  <- round(fst_mat, 6)
diag(fst_out) <- 0
out_fst_mat <- file.path(out_dir, "07.4-pairwise_fst_matrix.tsv")
write.table(
  cbind(population = rownames(fst_out), as.data.frame(fst_out)),
  file = out_fst_mat, sep = "\t", quote = FALSE, row.names = FALSE
)
cat("FST matrix written:", out_fst_mat, "\n")

# --- Export 2: p-value matrix ---
out_pval_mat <- file.path(out_dir, "07.4-pairwise_fst_pvalue_matrix.tsv")
write.table(
  cbind(population = rownames(pval_mat), as.data.frame(round(pval_mat, 4))),
  file = out_pval_mat, sep = "\t", quote = FALSE, row.names = FALSE
)
cat("FST p-value matrix written:", out_pval_mat, "\n")

# --- Export 3: long-format table with FST + CI + p-value + significance ---
sig_pairs <- do.call(rbind, lapply(seq_len(n_pops), function(i) {
  do.call(rbind, lapply(seq_len(n_pops), function(j) {
    if (j <= i) return(NULL)
    data.frame(
      pop1        = complete_pops[i],
      pop2        = complete_pops[j],
      weighted_fst = round(fst_mat[i, j], 6),
      FST_lower95 = lower_mat[i, j],
      FST_upper95 = upper_mat[i, j],
      pval        = pval_mat[i, j],
      significant = !is.na(pval_mat[i, j]) & pval_mat[i, j] < 0.05,
      stringsAsFactors = FALSE
    )
  }))
}))
sig_pairs <- sig_pairs[order(sig_pairs$pval, -sig_pairs$weighted_fst), ]
out_sig <- file.path(out_dir, "07.4-pairwise_fst_significance.tsv")
write.table(sig_pairs, file = out_sig, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("FST significance table written (%d pairs, %d significant at p<0.05): %s\n",
            nrow(sig_pairs),
            sum(sig_pairs$significant, na.rm = TRUE),
            out_sig))
cat(sprintf("INFO: bootstrap method — p-value = proportion of %d bootstrap means <= 0\n", NBOOT))

# Linearised FST and log distance
fst_lin <- fst_mat / (1 - fst_mat)
fst_lin[fst_mat == 0] <- 0
geo_ln  <- log(geo_mat)
diag(geo_ln) <- 0

# -------------------------------------------------------------------
# Mantel test
# -------------------------------------------------------------------
cat(sprintf("Running Mantel test (%d permutations)...\n", n_perm))
mantel_res <- vegan::mantel(as.dist(fst_lin), as.dist(geo_ln),
                             method = "pearson", permutations = n_perm)
cat(sprintf("Mantel r = %.4f | p = %.4f\n", mantel_res$statistic, mantel_res$signif))

# Save results table
res_df <- data.frame(
  group          = group,
  n_populations  = n_pops,
  mantel_r       = round(mantel_res$statistic, 6),
  p_value        = mantel_res$signif,
  n_permutations = n_perm,
  method         = "Pearson",
  x_variable     = "ln(geographic distance km)",
  y_variable     = "FST / (1 - FST)",
  stringsAsFactors = FALSE
)
out_tsv <- file.path(out_dir, "07.4-mantel_results_guiana_only.tsv")
write.table(res_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Mantel results written:", out_tsv, "\n")

# -------------------------------------------------------------------
# Upper triangle data for scatter plots
# -------------------------------------------------------------------
idx <- which(upper.tri(fst_lin), arr.ind = TRUE)
plot_df <- data.frame(
  pop1    = complete_pops[idx[, 1]],
  pop2    = complete_pops[idx[, 2]],
  geo_km  = geo_mat[idx],
  geo_ln  = geo_ln[idx],
  fst     = fst_mat[idx],
  fst_lin = fst_lin[idx],
  stringsAsFactors = FALSE
)

lm_fit    <- lm(fst_lin ~ geo_ln, data = plot_df)
r2        <- summary(lm_fit)$r.squared
pval_lbl  <- if (mantel_res$signif < 0.001) "p < 0.001" else
             sprintf("p = %.3f", mantel_res$signif)
annot_lbl <- sprintf("Mantel r = %.3f\n%s\nR\u00b2 = %.3f",
                     mantel_res$statistic, pval_lbl, r2)

# -------------------------------------------------------------------
# PLOT 1: IBD scatter — FST/(1-FST) vs ln(distance)
# -------------------------------------------------------------------
x_rng <- range(plot_df$geo_ln)
y_rng <- range(plot_df$fst_lin)

p1 <- ggplot(plot_df, aes(x = geo_ln, y = fst_lin, color = pop1)) +
  geom_point(size = 2.8, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              linewidth = 0.9, linetype = "solid", inherit.aes = FALSE,
              mapping = aes(x = geo_ln, y = fst_lin)) +
  annotate("label", fill = "white", label.size = 0,
           x = x_rng[1] + 0.04 * diff(x_rng),
           y = y_rng[2] - 0.02 * diff(y_rng),
           label = annot_lbl,
           hjust = 0, vjust = 1, size = 3.3, family = "mono") +
  labs(
    title    = "Isolation by Distance — Guiana populations",
    subtitle = expression(F[ST]/(1-F[ST]) ~ "vs." ~ ln(geographic~distance~(km))),
    x = "ln(geographic distance, km)",
    y = expression(F[ST] / (1 - F[ST])),
    color = "Population (pop1)"
  ) +
  theme(legend.key.size = unit(0.35, "cm"),
        legend.text     = element_text(size = 7))

if (length(site_colors) > 0)
  p1 <- p1 + scale_color_manual(values = site_colors, na.value = "grey50")

ggsave(file.path(plot_dir, "07.4-ibd_scatter_fst_lin_vs_ln_dist.png"),
       p1, width = 9, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "07.4-ibd_scatter_fst_lin_vs_ln_dist.pdf"),
       p1, width = 9, height = 6)
cat("Plot 1 written: IBD scatter FST/(1-FST) vs ln(dist)\n")

# -------------------------------------------------------------------
# PLOT 2: IBD scatter — raw FST vs geographic distance (km)
# -------------------------------------------------------------------
lm_raw   <- lm(fst ~ geo_km, data = plot_df)
r2_raw   <- summary(lm_raw)$r.squared
annot_raw <- sprintf("R\u00b2 = %.3f", r2_raw)

p2 <- ggplot(plot_df, aes(x = geo_km, y = fst, color = pop1)) +
  geom_point(size = 2.8, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              linewidth = 0.9, inherit.aes = FALSE,
              mapping = aes(x = geo_km, y = fst)) +
  annotate("label", fill = "white", label.size = 0,
           x = min(plot_df$geo_km) + 0.04 * diff(range(plot_df$geo_km)),
           y = max(plot_df$fst)    - 0.02 * diff(range(plot_df$fst)),
           label = annot_raw, hjust = 0, vjust = 1, size = 3.3) +
  labs(
    title    = "Isolation by Distance — Guiana populations (raw FST)",
    subtitle = expression(F[ST] ~ "vs. geographic distance (km)"),
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    color = "Population (pop1)"
  ) +
  theme(legend.key.size = unit(0.35, "cm"),
        legend.text     = element_text(size = 7))

if (length(site_colors) > 0)
  p2 <- p2 + scale_color_manual(values = site_colors, na.value = "grey50")

ggsave(file.path(plot_dir, "07.4-ibd_scatter_fst_vs_dist_km.png"),
       p2, width = 9, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "07.4-ibd_scatter_fst_vs_dist_km.pdf"),
       p2, width = 9, height = 6)
cat("Plot 2 written: IBD scatter raw FST vs km\n")

# -------------------------------------------------------------------
# PLOT 3: Heatmap — pairwise FST matrix
# -------------------------------------------------------------------
fst_long <- do.call(rbind, lapply(seq_len(n_pops), function(i) {
  do.call(rbind, lapply(seq_len(n_pops), function(j) {
    data.frame(pop1 = complete_pops[i], pop2 = complete_pops[j],
               fst  = fst_mat[i, j], stringsAsFactors = FALSE)
  }))
}))
fst_long$pop1 <- factor(fst_long$pop1, levels = complete_pops)
fst_long$pop2 <- factor(fst_long$pop2, levels = rev(complete_pops))

p3 <- ggplot(fst_long, aes(x = pop1, y = pop2, fill = fst)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(!is.na(fst) & pop1 != pop2,
                               sprintf("%.3f", fst), "")),
            size = 2.2, color = "black") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#d6604d",
                       midpoint = median(fst_mat[upper.tri(fst_mat)], na.rm = TRUE),
                       name = expression(F[ST]), na.value = "grey90") +
  labs(
    title = expression("Pairwise" ~ F[ST] ~ "matrix — Guiana populations"),
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y  = element_text(size = 7),
    panel.grid   = element_blank(),
    legend.key.height = unit(0.8, "cm")
  )

hw <- max(7, n_pops * 0.45)
ggsave(file.path(plot_dir, "07.4-ibd_heatmap_fst.png"),
       p3, width = hw, height = hw * 0.85, dpi = 300)
ggsave(file.path(plot_dir, "07.4-ibd_heatmap_fst.pdf"),
       p3, width = hw, height = hw * 0.85)
cat("Plot 3 written: FST heatmap\n")

# -------------------------------------------------------------------
# PLOT 4: Heatmap — pairwise geographic distances (km)
# -------------------------------------------------------------------
geo_long <- do.call(rbind, lapply(seq_len(n_pops), function(i) {
  do.call(rbind, lapply(seq_len(n_pops), function(j) {
    data.frame(pop1 = complete_pops[i], pop2 = complete_pops[j],
               dist_km = geo_mat[i, j], stringsAsFactors = FALSE)
  }))
}))
geo_long$pop1 <- factor(geo_long$pop1, levels = complete_pops)
geo_long$pop2 <- factor(geo_long$pop2, levels = rev(complete_pops))

p4 <- ggplot(geo_long, aes(x = pop1, y = pop2, fill = dist_km)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(dist_km > 0, sprintf("%.0f", dist_km), "")),
            size = 2.0, color = "black") +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                      name = "Distance\n(km)") +
  labs(
    title = "Pairwise geographic distance matrix — Guiana populations (km)",
    x = NULL, y = NULL
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y  = element_text(size = 7),
    panel.grid   = element_blank(),
    legend.key.height = unit(0.8, "cm")
  )

ggsave(file.path(plot_dir, "07.4-ibd_heatmap_geo_dist.png"),
       p4, width = hw, height = hw * 0.85, dpi = 300)
ggsave(file.path(plot_dir, "07.4-ibd_heatmap_geo_dist.pdf"),
       p4, width = hw, height = hw * 0.85)
cat("Plot 4 written: geographic distance heatmap\n")

# -------------------------------------------------------------------
# PLOT 5: Mantel permutation distribution
# -------------------------------------------------------------------
perm_r    <- as.numeric(mantel_res$perm)
obs_r     <- mantel_res$statistic
perm_df   <- data.frame(r = perm_r)
pval_lbl2 <- if (mantel_res$signif < 0.001) "p < 0.001" else
             sprintf("p = %.3f", mantel_res$signif)

p5 <- ggplot(perm_df, aes(x = r)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_density(color = "steelblue4", linewidth = 0.7) +
  geom_vline(xintercept = obs_r, color = "red", linewidth = 1.1, linetype = "dashed") +
  annotate("label", fill = "white", label.size = 0,
           x = obs_r, y = Inf,
           label = sprintf("Observed r = %.3f\n%s", obs_r, pval_lbl2),
           hjust = -0.1, vjust = 1.2, color = "red", size = 3.3) +
  labs(
    title    = "Mantel test — permutation distribution",
    subtitle = sprintf("%d permutations | Pearson correlation", n_perm),
    x = "Mantel statistic r (permuted)",
    y = "Density"
  )

ggsave(file.path(plot_dir, "07.4-ibd_mantel_permutation_dist.png"),
       p5, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(plot_dir, "07.4-ibd_mantel_permutation_dist.pdf"),
       p5, width = 7, height = 4.5)
cat("Plot 5 written: Mantel permutation distribution\n")

cat(sprintf("\nDONE 07.4 IBD analysis — %d plots written to: %s\n",
            5L, plot_dir))
print(res_df)

#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.3-genetic_diversity_cp.R
# Description : Chloroplast genetic diversity statistics for Dicorynia guianensis.
#               Computes Hd, pi, Watterson's theta, Tajima's D, Fu's Fs,
#               rarefied haplotype richness, and pairwise FST per sampling site.
# Author  : Julien Bonnier
# Usage   : Rscript --vanilla 07.3-genetic_diversity_cp.R \
#               <alignment_fasta> <haplotype_table> <metadata_csv> <output_dir>
# =============================================================================

.libPaths(c(path.expand("~/work/R"), .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

theme_set(theme_minimal(base_size = 11) + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

set.seed(42)  # reproducibility for Monte Carlo (Fu's Fs)

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: script.R <aln> <hap_table> <metadata> <outdir>")

aln_file   <- args[1]
hap_file   <- args[2]
meta_file  <- args[3]
output_dir <- args[4]

cat("INFO: alignment    =", aln_file, "\n")
cat("INFO: hap_table    =", hap_file, "\n")
cat("INFO: metadata     =", meta_file, "\n")
cat("INFO: output_dir   =", output_dir, "\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Read data ─────────────────────────────────────────────────────────────────
cat("INFO: reading alignment...\n")
aln <- read.dna(aln_file, format = "fasta")
cat("INFO: alignment loaded —", nrow(aln), "sequences x", ncol(aln), "bp\n")

# Subset to segregating sites only (speeds up computations)
seg_idx <- seg.sites(aln)
cat("INFO: segregating sites:", length(seg_idx), "\n")
aln_seg <- aln[, seg_idx]

meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(sample_id, site) %>%
  distinct() %>%
  mutate(
    # Recode Herbier to D.paranensis for display consistency
    site   = ifelse(site == "Herbier", "D.paranensis", site),
    region = case_when(
      site == "Cameroun_Benin" ~ "Outgroup",
      site == "D.paranensis"   ~ "Outgroup",
      TRUE                      ~ "French Guiana"
    )
  )

# French Guiana samples only — outgroups excluded from all analyses
meta_fg <- meta %>% filter(region == "French Guiana")

hap_table <- read.table(hap_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  select(sample_id, haplotype_id)

sample_info <- hap_table %>% left_join(meta_fg, by = "sample_id") %>% filter(!is.na(site))
cat("INFO: metadata loaded —", nrow(meta_fg), "FG samples,", n_distinct(meta_fg$site), "sites (outgroups excluded)\n")

# ── Helper: Fu's Fs via Monte Carlo Ewens sampling ────────────────────────────
# Simulates the distribution of K (number of haplotypes) under Ewens sampling
# formula given theta estimated from nucleotide diversity (theta = pi * L).
# Fs = ln(S / (1-S)) where S = P(K_sim <= k_obs).
# Significantly negative Fs → excess of haplotypes → population expansion signal.
compute_fu_fs <- function(n, k_obs, theta, n_sim = 5000) {
  if (n < 3 || is.na(theta) || theta <= 0) return(list(Fs = NA_real_, pval = NA_real_))
  sim_k <- replicate(n_sim, {
    k <- 1L
    for (i in 2:n) if (runif(1) < theta / (theta + i - 1L)) k <- k + 1L
    k
  })
  S   <- mean(sim_k <= k_obs)
  S   <- max(0.5 / n_sim, min(1 - 0.5 / n_sim, S))
  Fs  <- log(S / (1 - S))
  list(Fs = round(Fs, 4), pval = round(S, 4))
}

# ── Helper: compute diversity stats for a DNAbin subset ───────────────────────
compute_diversity <- function(dna_subset, label) {
  n  <- nrow(dna_subset)
  if (n < 2) {
    return(data.frame(
      group        = label,
      n_samples    = n,
      n_haplotypes = NA_integer_,
      Hd           = NA_real_,
      pi           = NA_real_,
      theta_w      = NA_real_,
      tajima_D     = NA_real_,
      tajima_pval  = NA_real_,
      fu_fs        = NA_real_,
      fu_fs_pval   = NA_real_
    ))
  }

  haps   <- haplotype(dna_subset)
  hd     <- hap.div(dna_subset)
  pi_val <- nuc.div(dna_subset)

  # Watterson's theta (per site, same normalization as pi)
  # theta_W = S / (a1 * L) where S = segregating sites in subset, L = total columns
  s_local <- length(seg.sites(dna_subset))
  a1      <- sum(1 / seq_len(n - 1))
  theta_w <- if (s_local > 0) round(s_local / a1 / ncol(dna_subset), 6) else 0

  # Tajima's D
  taj <- tryCatch(
    tajima.test(dna_subset),
    error = function(e) list(D = NA_real_, Pval.normal = NA_real_)
  )

  # Fu's Fs (Monte Carlo Ewens): theta estimated as total mean pairwise diffs
  theta_seq <- pi_val * ncol(dna_subset)
  fu <- compute_fu_fs(n, nrow(haps), theta_seq)

  data.frame(
    group        = label,
    n_samples    = n,
    n_haplotypes = nrow(haps),
    Hd           = round(hd, 4),
    pi           = round(pi_val, 6),
    theta_w      = theta_w,
    tajima_D     = round(taj$D, 4),
    tajima_pval  = round(taj$Pval.normal, 4),
    fu_fs        = fu$Fs,
    fu_fs_pval   = fu$pval
  )
}

# ── Global statistics (French Guiana only) ────────────────────────────────────
cat("INFO: computing global diversity (FG only)...\n")
fg_samps_global <- meta_fg$sample_id[meta_fg$sample_id %in% rownames(aln_seg)]
stats_global <- compute_diversity(aln_seg[fg_samps_global, , drop = FALSE], "Global (FG)")

# ── Per-site statistics (French Guiana only) ──────────────────────────────────
cat("INFO: computing per-site diversity...\n")

sites <- sort(unique(meta_fg$site))
stats_sites <- lapply(sites, function(s) {
  samples_s <- meta_fg$sample_id[meta_fg$site == s]
  samples_s <- samples_s[samples_s %in% rownames(aln_seg)]
  if (length(samples_s) < 1) return(NULL)
  cat("  site:", s, "—", length(samples_s), "samples\n")
  compute_diversity(aln_seg[samples_s, , drop = FALSE], s)
})
stats_sites <- bind_rows(Filter(Negate(is.null), stats_sites))

# ── Private haplotypes per site ───────────────────────────────────────────────
cat("INFO: computing private haplotypes...\n")

hap_site <- sample_info %>%
  group_by(haplotype_id) %>%
  summarise(sites_present = n_distinct(site), .groups = "drop")

private_haps <- sample_info %>%
  left_join(hap_site, by = "haplotype_id") %>%
  filter(sites_present == 1) %>%
  group_by(site) %>%
  summarise(n_private_haplotypes = n_distinct(haplotype_id), .groups = "drop")

# ── Rarefied haplotype richness ───────────────────────────────────────────────
cat("INFO: computing rarefied haplotype richness...\n")

# Rarefaction to the smallest FG site with >= 2 samples
fg_sizes <- meta_fg %>%
  filter(sample_id %in% rownames(aln_seg)) %>%
  group_by(site) %>% summarise(n = n(), .groups = "drop")
min_n_rarefy <- max(2L, min(fg_sizes$n))
cat("INFO: rarefaction target n =", min_n_rarefy, "\n")

compute_rarefied_richness <- function(dna_subset, min_n) {
  n <- nrow(dna_subset)
  if (n < min_n) return(NA_real_)
  haps <- haplotype(dna_subset)
  freq <- sapply(attr(haps, "index"), length)
  # Exact formula: E[K_min_n] = sum_h [ 1 - C(n-f_h, min_n) / C(n, min_n) ]
  expected_k <- sum(vapply(freq, function(f) {
    if (n - f < min_n) return(1)
    1 - exp(lchoose(n - f, min_n) - lchoose(n, min_n))
  }, numeric(1)))
  round(expected_k, 2)
}

rarefy_df <- do.call(rbind, lapply(sort(unique(meta_fg$site)), function(s) {
  samp <- meta_fg$sample_id[meta_fg$site == s]
  samp <- samp[samp %in% rownames(aln_seg)]
  data.frame(
    site              = s,
    rarefied_richness = compute_rarefied_richness(aln_seg[samp, , drop = FALSE], min_n_rarefy)
  )
}))

# ── Pairwise FST between FG sites (Hudson 1992 estimator) ─────────────────────
cat("INFO: computing pairwise FST...\n")

fg_sites  <- sort(unique(meta_fg$site))
fg_samps  <- meta_fg$sample_id[meta_fg$sample_id %in% rownames(aln_seg)]

# Precompute full pairwise distance matrix for FG samples (raw = proportion differences)
dist_mat <- as.matrix(dist.dna(aln_seg[fg_samps, , drop = FALSE],
                                model = "raw", pairwise.deletion = TRUE))

hudson_fst <- function(mat, idx_a, idx_b) {
  if (length(idx_a) < 2 || length(idx_b) < 2) return(NA_real_)
  m_aa <- mat[idx_a, idx_a]; dAA <- mean(m_aa[upper.tri(m_aa)])
  m_bb <- mat[idx_b, idx_b]; dBB <- mean(m_bb[upper.tri(m_bb)])
  dAB  <- mean(mat[idx_a, idx_b])
  if (is.na(dAB) || dAB == 0) return(NA_real_)
  round(max(0, min(1, (dAB - (dAA + dBB) / 2) / dAB)), 4)
}

n_sites  <- length(fg_sites)
fst_mat  <- matrix(NA_real_, n_sites, n_sites, dimnames = list(fg_sites, fg_sites))
diag(fst_mat) <- 0

for (i in seq_len(n_sites - 1)) {
  for (j in (i + 1):n_sites) {
    sa <- fg_samps[meta_fg$site[match(fg_samps, meta_fg$sample_id)] == fg_sites[i]]
    sb <- fg_samps[meta_fg$site[match(fg_samps, meta_fg$sample_id)] == fg_sites[j]]
    fst_mat[i, j] <- fst_mat[j, i] <- hudson_fst(dist_mat, sa, sb)
  }
}

# Save FST matrix as TSV
fst_df <- as.data.frame(fst_mat)
fst_df <- cbind(site = rownames(fst_df), fst_df)
write.table(fst_df, file.path(output_dir, "pairwise_fst.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("INFO: saved pairwise_fst.tsv\n")

# ── Merge all site stats ───────────────────────────────────────────────────────
stats_sites_full <- stats_sites %>%
  left_join(private_haps, by = c("group" = "site")) %>%
  mutate(n_private_haplotypes = replace_na(n_private_haplotypes, 0)) %>%
  left_join(rarefy_df, by = c("group" = "site")) %>%
  left_join(meta_fg %>% select(site, region) %>% distinct(), by = c("group" = "site"))

# ── Write stats table ─────────────────────────────────────────────────────────
out_stats <- file.path(output_dir, "diversity_stats_per_site.tsv")
write.table(
  bind_rows(
    stats_global %>% mutate(region = "French Guiana", n_private_haplotypes = NA, rarefied_richness = NA),
    stats_sites_full
  ) %>% select(group, region, n_samples, n_haplotypes, Hd, pi, theta_w,
               tajima_D, tajima_pval, fu_fs, fu_fs_pval,
               n_private_haplotypes, rarefied_richness),
  file = out_stats, sep = "\t", row.names = FALSE, quote = FALSE
)
cat("INFO: saved", out_stats, "\n")

# ── Color palette — French Guiana sites only ──────────────────────────────────
sites_fg  <- sort(unique(meta_fg$site))
pal_fg    <- colorRampPalette(brewer.pal(12, "Paired"))(length(sites_fg))
site_cols <- setNames(pal_fg, sites_fg)

site_order <- sites_fg  # FG sites alphabetically

stats_fg <- stats_sites_full %>%
  filter(group %in% site_order) %>%
  mutate(group = factor(group, levels = site_order))

# ── Plot 1: Haplotype diversity (Hd) per site ─────────────────────────────────
cat("INFO: generating Hd plot...\n")

p_hd <- ggplot(stats_fg, aes(x = group, y = Hd, fill = group)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0("n=", n_samples)), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = site_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Haplotype diversity per site"),
    subtitle = paste0("Chloroplast | n=", nrow(meta), " individuals | ", nrow(hap_site), " haplotypes"),
    x        = "Sampling site",
    y        = "Haplotype diversity (Hd)"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "cp_diversity_hd_per_site.png"),
       p_hd, width = 14, height = 7, dpi = 300)
cat("INFO: saved cp_diversity_hd_per_site.png\n")

# ── Plot 2: Nucleotide diversity (pi) per site ────────────────────────────────
cat("INFO: generating pi plot...\n")

p_pi <- ggplot(stats_fg, aes(x = group, y = pi, fill = group)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0("n=", n_samples)), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = site_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Nucleotide diversity per site"),
    subtitle = paste0("Chloroplast | n=", nrow(meta), " individuals"),
    x        = "Sampling site",
    y        = expression(paste("Nucleotide diversity (", pi, ")"))
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "cp_diversity_pi_per_site.png"),
       p_pi, width = 14, height = 7, dpi = 300)
cat("INFO: saved cp_diversity_pi_per_site.png\n")

# ── Plot 3: Hd vs pi scatter per site ────────────────────────────────────────
cat("INFO: generating Hd vs pi scatter...\n")

stats_fg_clean <- stats_fg %>% filter(!is.na(Hd) & !is.na(pi))

p_scatter <- ggplot(stats_fg_clean, aes(x = pi, y = Hd, color = group, size = n_samples)) +
  geom_point(alpha = 0.85) +
  geom_text_repel(aes(label = group), size = 2.8, max.overlaps = 20) +
  scale_color_manual(values = site_cols, guide = "none") +
  scale_size_continuous(name = "n individuals", range = c(3, 9)) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Hd vs. pi per site"),
    subtitle = "Chloroplast | point size = sample size",
    x        = expression(paste("Nucleotide diversity (", pi, ")")),
    y        = "Haplotype diversity (Hd)"
  )

ggsave(file.path(output_dir, "cp_diversity_hd_vs_pi.png"),
       p_scatter, width = 10, height = 8, dpi = 300)
cat("INFO: saved cp_diversity_hd_vs_pi.png\n")

# ── Plot 4: Number of haplotypes + private haplotypes per site ────────────────
cat("INFO: generating haplotype count plot...\n")

hap_counts <- stats_fg %>%
  filter(!is.na(n_haplotypes)) %>%
  select(group, n_haplotypes, n_private_haplotypes) %>%
  pivot_longer(cols = c(n_haplotypes, n_private_haplotypes),
               names_to = "type", values_to = "count") %>%
  mutate(type = recode(type,
    "n_haplotypes"         = "Total haplotypes",
    "n_private_haplotypes" = "Private haplotypes"
  ))

p_hap_count <- ggplot(hap_counts, aes(x = group, y = count, fill = type)) +
  geom_col(position = "dodge", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Total haplotypes" = "#4DAF4A", "Private haplotypes" = "#984EA3"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Haplotype richness per site"),
    subtitle = "Total vs. private haplotypes (chloroplast)",
    x        = "Sampling site",
    y        = "Number of haplotypes"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(output_dir, "cp_diversity_haplotype_counts.png"),
       p_hap_count, width = 14, height = 7, dpi = 300)
cat("INFO: saved cp_diversity_haplotype_counts.png\n")

# ── Plot 5: Rarefied haplotype richness per site ──────────────────────────────
cat("INFO: generating rarefied richness plot...\n")

stats_fg_rarefy <- stats_fg %>%
  filter(!is.na(rarefied_richness))

p_rarefy <- ggplot(stats_fg_rarefy, aes(x = group, y = rarefied_richness, fill = group)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0("n=", n_samples)), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = site_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Rarefied haplotype richness per site"),
    subtitle = paste0("Chloroplast | rarefied to n=", min_n_rarefy, " individuals per site"),
    caption  = "Exact rarefaction formula (Coleman 1981) — sites with n < target excluded",
    x        = "Sampling site",
    y        = paste0("Expected haplotypes (n=", min_n_rarefy, ")")
  ) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "cp_diversity_rarefied_richness.png"),
       p_rarefy, width = 14, height = 7, dpi = 300)
cat("INFO: saved cp_diversity_rarefied_richness.png\n")

# ── Plot 6: Pairwise FST heatmap ───────────────────────────────────────────────
cat("INFO: generating pairwise FST heatmap...\n")

fst_long <- as.data.frame(as.table(fst_mat)) %>%
  rename(site1 = Var1, site2 = Var2, fst = Freq) %>%
  mutate(site1 = factor(site1, levels = fg_sites),
         site2 = factor(site2, levels = rev(fg_sites)))

p_fst <- ggplot(fst_long, aes(x = site1, y = site2, fill = fst)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(fst), "", sprintf("%.2f", fst))),
            size = 2.2, color = "grey20") +
  scale_fill_gradientn(
    colours  = c("#FFFFFF", "#FEE391", "#FEC44F", "#FE9929", "#CC4C02"),
    na.value = "grey90",
    name     = expression(F[ST]),
    limits   = c(0, 1)
  ) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Pairwise F"[ST]~"between sites"),
    subtitle = "Chloroplast | Hudson (1992) estimator | French Guiana sites only",
    x        = NULL, y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.key.height = unit(1.2, "cm")
  ) +
  coord_fixed()

ggsave(file.path(output_dir, "cp_diversity_pairwise_fst.png"),
       p_fst, width = 12, height = 11, dpi = 300)
cat("INFO: saved cp_diversity_pairwise_fst.png\n")

# ── Print summary ─────────────────────────────────────────────────────────────
cat("\n── Global diversity summary ──────────────────────────────────────────\n")
cat(sprintf("  Samples     : %d\n", stats_global$n_samples))
cat(sprintf("  Haplotypes  : %d\n", stats_global$n_haplotypes))
cat(sprintf("  Hd (global) : %.4f\n", stats_global$Hd))
cat(sprintf("  pi (global) : %.6f\n", stats_global$pi))
cat(sprintf("  Tajima's D  : %.4f (p=%.4f)\n", stats_global$tajima_D, stats_global$tajima_pval))
cat("──────────────────────────────────────────────────────────────────────\n\n")

cat("DONE genetic diversity statistics complete\n")

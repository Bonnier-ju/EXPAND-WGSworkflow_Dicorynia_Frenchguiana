#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.3-genetic_diversity_cp.R
# Description : Chloroplast genetic diversity statistics for Dicorynia guianensis.
#               Computes haplotype diversity (Hd), nucleotide diversity (pi),
#               and Tajima's D globally and per sampling site.
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
  mutate(region = ifelse(site == "Cameroun_Benin", "Outgroup", "French Guiana"))

hap_table <- read.table(hap_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  select(sample_id, haplotype_id)

sample_info <- hap_table %>% left_join(meta, by = "sample_id")
cat("INFO: metadata loaded —", nrow(meta), "samples,", n_distinct(meta$site), "sites\n")

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
      tajima_D     = NA_real_,
      tajima_pval  = NA_real_
    ))
  }

  # Haplotype diversity
  haps <- haplotype(dna_subset)
  hd   <- hap.div(dna_subset)

  # Nucleotide diversity (raw pairwise differences / sites)
  pi_val <- nuc.div(dna_subset)

  # Tajima's D (requires >= 3 sequences and segregating sites)
  taj <- tryCatch(
    tajima.test(dna_subset),
    error = function(e) list(D = NA_real_, Pval.normal = NA_real_)
  )

  data.frame(
    group        = label,
    n_samples    = n,
    n_haplotypes = nrow(haps),
    Hd           = round(hd, 4),
    pi           = round(pi_val, 6),
    tajima_D     = round(taj$D, 4),
    tajima_pval  = round(taj$Pval.normal, 4)
  )
}

# ── Global statistics ─────────────────────────────────────────────────────────
cat("INFO: computing global diversity...\n")
stats_global <- compute_diversity(aln_seg, "Global")

# ── Per-site statistics ───────────────────────────────────────────────────────
cat("INFO: computing per-site diversity...\n")

sites <- sort(unique(meta$site))
stats_sites <- lapply(sites, function(s) {
  samples_s <- meta$sample_id[meta$site == s]
  samples_s <- samples_s[samples_s %in% rownames(aln_seg)]
  if (length(samples_s) < 1) return(NULL)
  cat("  site:", s, "—", length(samples_s), "samples\n")
  compute_diversity(aln_seg[samples_s, , drop = FALSE], s)
})
stats_sites <- bind_rows(Filter(Negate(is.null), stats_sites))

# ── Per-region statistics ─────────────────────────────────────────────────────
cat("INFO: computing per-region diversity...\n")

regions <- sort(unique(meta$region))
stats_regions <- lapply(regions, function(r) {
  samples_r <- meta$sample_id[meta$region == r]
  samples_r <- samples_r[samples_r %in% rownames(aln_seg)]
  compute_diversity(aln_seg[samples_r, , drop = FALSE], paste0("Region: ", r))
})
stats_regions <- bind_rows(Filter(Negate(is.null), stats_regions))

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

# ── Merge all site stats ───────────────────────────────────────────────────────
stats_sites_full <- stats_sites %>%
  left_join(private_haps, by = c("group" = "site")) %>%
  mutate(n_private_haplotypes = replace_na(n_private_haplotypes, 0)) %>%
  left_join(meta %>% select(site, region) %>% distinct(), by = c("group" = "site"))

# ── Write stats table ─────────────────────────────────────────────────────────
out_stats <- file.path(output_dir, "diversity_stats_per_site.tsv")
write.table(
  bind_rows(
    stats_global %>% mutate(region = "All", n_private_haplotypes = NA),
    stats_regions %>% mutate(region = NA, n_private_haplotypes = NA),
    stats_sites_full %>% rename(group = group)
  ) %>% select(group, region, n_samples, n_haplotypes, Hd, pi, tajima_D, tajima_pval, n_private_haplotypes),
  file = out_stats, sep = "\t", row.names = FALSE, quote = FALSE
)
cat("INFO: saved", out_stats, "\n")

# ── Color palette (consistent with 07.2) ─────────────────────────────────────
sites_fg  <- sort(unique(meta$site[meta$region == "French Guiana"]))
pal_fg    <- colorRampPalette(brewer.pal(12, "Paired"))(length(sites_fg))
site_cols <- setNames(c(pal_fg, "#E41A1C"), c(sites_fg, "Cameroun_Benin"))

# Ordered for plots: FG sites alphabetically, outgroup last
site_order <- c(sites_fg, "Cameroun_Benin")

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

# ── Print summary ─────────────────────────────────────────────────────────────
cat("\n── Global diversity summary ──────────────────────────────────────────\n")
cat(sprintf("  Samples     : %d\n", stats_global$n_samples))
cat(sprintf("  Haplotypes  : %d\n", stats_global$n_haplotypes))
cat(sprintf("  Hd (global) : %.4f\n", stats_global$Hd))
cat(sprintf("  pi (global) : %.6f\n", stats_global$pi))
cat(sprintf("  Tajima's D  : %.4f (p=%.4f)\n", stats_global$tajima_D, stats_global$tajima_pval))
cat("──────────────────────────────────────────────────────────────────────\n\n")

cat("DONE genetic diversity statistics complete\n")

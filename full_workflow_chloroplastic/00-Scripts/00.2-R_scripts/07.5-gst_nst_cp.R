#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.5-gst_nst_cp.R
# Description : Comparison of GST (allelic differentiation) vs NST
#               (phylogeographic differentiation) following Pons & Petit (1996).
#               Tests whether NST > GST, indicating significant phylogeographic
#               structure (geographically proximate populations share more
#               closely related chloroplast haplotypes).
#               Outputs:
#                 - cp_gst_nst_summary.tsv        : GST, NST, p-value
#                 - cp_gst_nst_barplot.png         : GST vs NST observed values
#                 - cp_gst_nst_permutation.png     : permutation distribution
#                 - cp_gst_nst_site_diversity.png  : per-site Hd and pi
# Author  : Julien Bonnier
# Usage   : Rscript --vanilla 07.5-gst_nst_cp.R \
#               <alignment_fasta> <metadata_csv> <output_dir>
# =============================================================================

.libPaths(c(path.expand("~/work/R"), .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ── Global white background ───────────────────────────────────────────────────
theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 07.5-gst_nst_cp.R <aln_fasta> <metadata_csv> <output_dir>")
}
aln_file   <- args[1]
meta_file  <- args[2]
output_dir <- args[3]

cat("INFO: alignment  =", aln_file,   "\n")
cat("INFO: metadata   =", meta_file,  "\n")
cat("INFO: output_dir =", output_dir, "\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
N_PERM <- 5000

# ── Read alignment ─────────────────────────────────────────────────────────────
cat("INFO: reading alignment...\n")
aln <- read.dna(aln_file, format = "fasta")
cat("INFO:", nrow(aln), "sequences x", ncol(aln), "bp\n")

seg_idx <- seg.sites(aln)
cat("INFO: subsetting to", length(seg_idx), "segregating sites\n")
aln_seg <- aln[, seg_idx]

# ── Read metadata ──────────────────────────────────────────────────────────────
meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(sample_id, site, project) %>%
  distinct() %>%
  mutate(
    site = case_when(
      project == "Treemutation" ~ "Angela",
      site    == "Herbier"      ~ "D.paranensis",
      TRUE                      ~ site
    )
  )

# Exclude Angela and outgroups — FG Dicorynia only
EXCLUDE <- c("Angela", "D.paranensis", "Cameroun_Benin")
meta_fg <- meta %>% filter(!site %in% EXCLUDE)

cat("INFO: FG samples:", nrow(meta_fg), "individuals,",
    n_distinct(meta_fg$site), "sites\n")

# ── Subset alignment to FG samples ────────────────────────────────────────────
fg_samps <- meta_fg$sample_id[meta_fg$sample_id %in% rownames(aln_seg)]
aln_fg   <- aln_seg[fg_samps, ]
cat("INFO: FG alignment:", nrow(aln_fg), "sequences x", ncol(aln_fg), "sites\n")

# Site assignment vector (same order as aln_fg rows)
site_vec <- meta_fg$site[match(rownames(aln_fg), meta_fg$sample_id)]
sites    <- sort(unique(site_vec))
n_sites  <- length(sites)
site_n   <- as.integer(table(site_vec)[sites])
cat("INFO:", n_sites, "sites | min n =", min(site_n), "| max n =", max(site_n), "\n")

# ── Precompute haplotype assignments (for fast GST permutation) ────────────────
cat("INFO: computing haplotype assignments...\n")
haps_all   <- haplotype(aln_fg)
hap_index  <- attr(haps_all, "index")   # list: haplotype i -> row indices in aln_fg
hap_assign <- integer(nrow(aln_fg))     # sample -> haplotype number
for (h in seq_along(hap_index)) hap_assign[hap_index[[h]]] <- h

# ── Precompute pairwise raw distances (for fast NST permutation) ───────────────
cat("INFO: computing pairwise distances...\n")
D_mat <- as.matrix(dist.dna(aln_fg, model = "raw", pairwise.deletion = TRUE))

# ── Fast GST computation (haplotype-frequency-based) ──────────────────────────
# Uses precomputed hap_assign — no haplotype() call per permutation
fast_gst <- function(hap_asgn, sv) {
  n_total <- length(sv)
  sv_u    <- unique(sv)

  # HT: total haplotype diversity
  cnt_total <- table(hap_asgn)
  p_total   <- cnt_total / n_total
  HT <- n_total / (n_total - 1) * (1 - sum(p_total^2))

  # HS: weighted mean within-site Hd
  hs_vals <- vapply(sv_u, function(s) {
    idx <- which(sv == s)
    n_s <- length(idx)
    if (n_s < 2L) return(NA_real_)
    cnt <- table(hap_asgn[idx])
    p_s <- cnt / n_s
    n_s / (n_s - 1) * (1 - sum(p_s^2))
  }, numeric(1))

  n_s_vec <- vapply(sv_u, function(s) sum(sv == s), integer(1))
  ok      <- !is.na(hs_vals)
  HS      <- sum(hs_vals[ok] * n_s_vec[ok]) / sum(n_s_vec[ok])

  if (HT <= 0) NA_real_ else (HT - HS) / HT
}

# ── Fast NST computation (distance-matrix-based) ──────────────────────────────
# piT is constant across permutations; only piS changes
fast_nst <- function(D, sv, piT) {
  sv_u    <- unique(sv)

  piS_vals <- vapply(sv_u, function(s) {
    idx <- which(sv == s)
    n_s <- length(idx)
    if (n_s < 2L) return(NA_real_)
    d_sub <- D[idx, idx]
    mean(d_sub[upper.tri(d_sub)])
  }, numeric(1))

  n_s_vec <- vapply(sv_u, function(s) sum(sv == s), integer(1))
  ok      <- !is.na(piS_vals)
  piS     <- sum(piS_vals[ok] * n_s_vec[ok]) / sum(n_s_vec[ok])

  if (piT <= 0) NA_real_ else (piT - piS) / piT
}

# ── Observed values ────────────────────────────────────────────────────────────
cat("INFO: computing observed GST and NST...\n")

obs_gst <- fast_gst(hap_assign, site_vec)

piT_obs <- mean(D_mat[upper.tri(D_mat)])
obs_nst <- fast_nst(D_mat, site_vec, piT_obs)

obs_diff <- obs_nst - obs_gst

cat(sprintf("INFO: GST      = %.4f\n", obs_gst))
cat(sprintf("INFO: NST      = %.4f\n", obs_nst))
cat(sprintf("INFO: NST-GST  = %.4f\n", obs_diff))

# ── Permutation test ───────────────────────────────────────────────────────────
cat("INFO: running", N_PERM, "permutations...\n")

perm_gst  <- numeric(N_PERM)
perm_nst  <- numeric(N_PERM)
perm_diff <- numeric(N_PERM)

for (i in seq_len(N_PERM)) {
  sv_p         <- sample(site_vec)
  perm_gst[i]  <- fast_gst(hap_assign, sv_p)
  perm_nst[i]  <- fast_nst(D_mat, sv_p, piT_obs)
  perm_diff[i] <- perm_nst[i] - perm_gst[i]
}

p_value <- mean(perm_diff >= obs_diff, na.rm = TRUE)
cat(sprintf("INFO: permutation p-value (NST > GST) = %.4f\n", p_value))

# ── Per-site diversity (for plot and table) ────────────────────────────────────
site_stats <- lapply(sites, function(s) {
  idx <- which(site_vec == s)
  n_s <- length(idx)
  if (n_s < 2) return(NULL)
  hd  <- hap.div(haplotype(aln_fg[idx, ]))
  pi  <- nuc.div(aln_fg[idx, ])
  data.frame(site = s, n = n_s, Hd = hd, pi = pi, stringsAsFactors = FALSE)
}) %>% bind_rows()

# Weighted HS and piS (for reference lines in plot)
HS_obs  <- sum(site_stats$Hd * site_stats$n) / sum(site_stats$n)
piS_obs <- sum(site_stats$pi * site_stats$n) / sum(site_stats$n)

# ── Summary table ──────────────────────────────────────────────────────────────
sig_label <- dplyr::case_when(
  p_value < 0.001 ~ "***",
  p_value < 0.01  ~ "**",
  p_value < 0.05  ~ "*",
  TRUE            ~ "ns"
)

summary_tbl <- data.frame(
  metric        = c("GST",   "NST",   "NST - GST"),
  observed      = c(obs_gst, obs_nst, obs_diff),
  perm_mean     = c(NA,      NA,      mean(perm_diff, na.rm = TRUE)),
  perm_ci_low   = c(NA,      NA,      quantile(perm_diff, 0.025, na.rm = TRUE)),
  perm_ci_high  = c(NA,      NA,      quantile(perm_diff, 0.975, na.rm = TRUE)),
  p_value       = c(NA,      NA,      p_value),
  significance  = c(NA,      NA,      sig_label),
  n_permutations = c(NA,     NA,      N_PERM)
)

write.table(summary_tbl,
  file = file.path(output_dir, "cp_gst_nst_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
cat("INFO: saved cp_gst_nst_summary.tsv\n")

# ── Plot 1: GST vs NST barplot ─────────────────────────────────────────────────
cat("INFO: generating Plot 1 — GST vs NST barplot...\n")

p1_data <- data.frame(
  index = factor(c("GST", "NST"), levels = c("GST", "NST")),
  value = c(obs_gst, obs_nst)
)

sig_annot <- sprintf("NST − GST = %.4f\n%s (p = %.4f, n = %d perm.)",
                     obs_diff, sig_label, p_value, N_PERM)

p1 <- ggplot(p1_data, aes(x = index, y = value, fill = index)) +
  geom_col(width = 0.45, alpha = 0.88) +
  geom_text(aes(label = sprintf("%.4f", value)),
            vjust = -0.5, size = 4.5, fontface = "bold") +
  annotate("text", x = 1.5, y = max(obs_gst, obs_nst) * 1.1,
           label = sig_annot, size = 3.5, color = "grey30", hjust = 0.5) +
  scale_fill_manual(values = c("GST" = "#4E9AF1", "NST" = "#E05C5C")) +
  scale_y_continuous(limits = c(0, max(obs_gst, obs_nst) * 1.35), expand = c(0, 0)) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast phylogeographic structure"),
    subtitle = "GST = allelic differentiation | NST = phylogeographic differentiation",
    x = NULL, y = "Differentiation index",
    caption  = "Pons & Petit (1996) | French Guiana populations only | Angela & outgroups excluded"
  ) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 11),
    plot.subtitle   = element_text(size = 8, color = "grey40"),
    plot.caption    = element_text(size = 7, color = "grey50"),
    axis.text.x     = element_text(size = 14, face = "bold")
  )

ggsave(file.path(output_dir, "cp_gst_nst_barplot.png"),
       p1, width = 6, height = 6, dpi = 300)
cat("INFO: saved cp_gst_nst_barplot.png\n")

# ── Plot 2: Permutation distribution ──────────────────────────────────────────
cat("INFO: generating Plot 2 — permutation distribution...\n")

perm_df <- data.frame(diff = perm_diff)

p2 <- ggplot(perm_df, aes(x = diff)) +
  geom_histogram(bins = 70, fill = "grey70", color = "white", alpha = 0.85) +
  geom_vline(xintercept = obs_diff,
             color = "#E05C5C", linewidth = 1.3, linetype = "dashed") +
  annotate("text",
           x = obs_diff, y = Inf, vjust = 2, hjust = -0.08,
           label = sprintf("Observed\nNST − GST = %.4f\n%s (p = %.4f)",
                           obs_diff, sig_label, p_value),
           color = "#E05C5C", size = 3.5, fontface = "bold") +
  labs(
    title    = "Permutation test: NST > GST ?",
    subtitle = sprintf("%d permutations of site assignments | Seed = 42", N_PERM),
    x        = "NST − GST (permuted)",
    y        = "Count",
    caption  = "Red dashed line = observed NST − GST. Proportion to its right = p-value."
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    plot.caption  = element_text(size = 7, color = "grey50")
  )

ggsave(file.path(output_dir, "cp_gst_nst_permutation.png"),
       p2, width = 8, height = 5, dpi = 300)
cat("INFO: saved cp_gst_nst_permutation.png\n")

# ── Plot 3: Per-site Hd and pi ────────────────────────────────────────────────
cat("INFO: generating Plot 3 — per-site diversity...\n")

site_long <- site_stats %>%
  pivot_longer(cols = c(Hd, pi), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("Hd", "pi"),
                    labels = c("Haplotype diversity (Hd)",
                               "Nucleotide diversity (π)")),
    site = reorder(site, value, FUN = mean)
  )

ref_lines <- data.frame(
  metric = factor(c("Haplotype diversity (Hd)", "Nucleotide diversity (π)"),
                  levels = c("Haplotype diversity (Hd)", "Nucleotide diversity (π)")),
  ref    = c(HS_obs, piS_obs),
  lab    = c(sprintf("HS = %.3f", HS_obs), sprintf("πS = %.5f", piS_obs))
)

p3 <- ggplot(site_long, aes(x = site, y = value, fill = metric)) +
  geom_col(alpha = 0.85) +
  geom_hline(data = ref_lines,
             aes(yintercept = ref), color = "grey20",
             linetype = "dashed", linewidth = 0.7) +
  geom_text(data = ref_lines,
            aes(x = -Inf, y = ref, label = lab),
            hjust = -0.1, vjust = -0.4, size = 3, color = "grey20") +
  facet_wrap(~ metric, scales = "free_x") +
  scale_fill_manual(values = c(
    "Haplotype diversity (Hd)" = "#4E9AF1",
    "Nucleotide diversity (π)" = "#E05C5C"
  )) +
  coord_flip() +
  labs(
    title    = "Within-site diversity contributing to GST and NST",
    subtitle = "Dashed line = weighted population mean (HS / πS)",
    x = NULL, y = "Diversity index",
    caption  = "French Guiana only | Angela & outgroups excluded"
  ) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 11),
    plot.subtitle   = element_text(size = 8, color = "grey40"),
    plot.caption    = element_text(size = 7, color = "grey50"),
    strip.text      = element_text(face = "bold", size = 9)
  )

ggsave(file.path(output_dir, "cp_gst_nst_site_diversity.png"),
       p3, width = 11, height = 7, dpi = 300)
cat("INFO: saved cp_gst_nst_site_diversity.png\n")

cat("DONE GST/NST analysis complete — results in", output_dir, "\n")

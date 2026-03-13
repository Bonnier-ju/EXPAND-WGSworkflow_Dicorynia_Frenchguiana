#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.2-haplotype_networks_cp.R
# Description : Chloroplast haplotype network for Dicorynia guianensis.
#               Network topology from pegas::haploNet (MSN).
#               Layout via MDS on pairwise genetic distances.
#               Nodes = pie charts (scatterpie), fixed uniform size.
#               Two plots:
#                 1. All samples (outgroup + Dicorynia paranensis included)
#                 2. French Guiana only (outgroup + D.paranensis excluded)
# Author  : Julien Bonnier
# Usage   : Rscript --vanilla 07.2-haplotype_networks_cp.R \
#               <alignment_fasta> <haplotype_table> <metadata_csv> <output_dir>
# =============================================================================

.libPaths(c(path.expand("~/work/R"), .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(dplyr)
  library(ggplot2)
  library(scatterpie)
  library(RColorBrewer)
  library(ggrepel)
})

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

# ── Read alignment & subset to variable sites ──────────────────────────────────
cat("INFO: reading alignment...\n")
aln <- read.dna(aln_file, format = "fasta")
cat("INFO: alignment loaded —", nrow(aln), "sequences x", ncol(aln), "bp\n")

seg_idx <- seg.sites(aln)
cat("INFO: subsetting to", length(seg_idx), "segregating sites\n")
aln <- aln[, seg_idx]

# ── Read metadata ──────────────────────────────────────────────────────────────
hap_table <- read.table(hap_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  select(sample_id, haplotype_id)

meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(sample_id, site, project) %>%
  distinct() %>%
  mutate(
    site   = case_when(
      project == "Treemutation" ~ "Angela",
      site    == "Herbier"      ~ "Dicorynia paranensis",
      TRUE                      ~ site
    ),
    region = ifelse(site == "Cameroun_Benin", "Outgroup (Cameroun/Benin)", "French Guiana")
  )

cat("INFO: metadata loaded —", nrow(meta), "samples,", n_distinct(meta$site), "sites\n")

sample_info <- hap_table %>% left_join(meta, by = "sample_id")

# ── Compute haplotypes ─────────────────────────────────────────────────────────
cat("INFO: computing haplotypes...\n")
haps <- haplotype(aln)
cat("INFO:", nrow(haps), "unique haplotypes\n")

# ── Build MSN ─────────────────────────────────────────────────────────────────
cat("INFO: building minimum spanning network...\n")
net <- haploNet(haps, getProb = FALSE)
n_haps <- nrow(haps)
cat("INFO: network built —", n_haps, "nodes\n")

# ── Extract edge list from network ────────────────────────────────────────────
# pegas haploNet object: matrix columns = [from, to, steps, ...]
net_mat   <- unclass(net)
edge_from <- net_mat[, 1]
edge_to   <- net_mat[, 2]
edge_step <- net_mat[, 3]

# ── MDS layout on pairwise haplotype distances ────────────────────────────────
# Positions haplotypes by genetic distance: close haplotypes → close nodes
# pos_all : original MDS (no scaling) — used for all-samples plot
# pos_fg  : scaled x10 to spread FG nodes — used for FG-only plot
cat("INFO: computing MDS layout...\n")
d_mat <- as.matrix(dist.dna(haps, model = "N", pairwise.deletion = TRUE))
mds   <- cmdscale(d_mat, k = 2)

pos_all <- data.frame(hap_id = seq_len(n_haps),
                      x = mds[, 1],
                      y = mds[, 2])

pos_fg <- data.frame(hap_id = seq_len(n_haps),
                     x = mds[, 1] * 20,
                     y = mds[, 2] * 20)

# ── Site composition per haplotype ───────────────────────────────────────────
hap_membership <- attr(haps, "index")
seq_names      <- rownames(aln)

all_sites <- sort(unique(meta$site))

pie_sites <- matrix(0, nrow = n_haps, ncol = length(all_sites),
                    dimnames = list(NULL, all_sites))

for (i in seq_along(hap_membership)) {
  samples <- seq_names[hap_membership[[i]]]
  info    <- sample_info %>% filter(sample_id %in% samples)
  for (s in info$site) pie_sites[i, s] <- pie_sites[i, s] + 1
}

node_sizes <- sapply(hap_membership, length)

# ── Color palettes ─────────────────────────────────────────────────────────────
# FG sites (includes Angela, excludes Cameroun_Benin and Dicorynia paranensis)
sites_fg <- sort(unique(meta$site[meta$region == "French Guiana" &
                                   meta$site != "Dicorynia paranensis"]))
pal_fg   <- colorRampPalette(brewer.pal(12, "Paired"))(length(sites_fg))

site_colors <- setNames(
  c(pal_fg, "#E41A1C", "#888888"),
  c(sites_fg, "Cameroun_Benin", "Dicorynia paranensis")
)

# ── Node radii ────────────────────────────────────────────────────────────────
# r_pie      : radius for haplotypes with n >= 2 (shown as pie charts)
# r_singleton: radius for singletons (n == 1, shown as small grey dots)
r_pie       <- diff(range(mds[, 1])) * 0.012
r_singleton <- r_pie * 0.4

# ── Build edge data frames (one per layout) ───────────────────────────────────
make_edges <- function(pos_df) {
  data.frame(from = edge_from, to = edge_to, steps = edge_step) %>%
    left_join(pos_df, by = c("from" = "hap_id")) %>%
    rename(x_from = x, y_from = y) %>%
    left_join(pos_df, by = c("to" = "hap_id")) %>%
    rename(x_to = x, y_to = y)
}

edges_all <- make_edges(pos_all)
edges_fg  <- make_edges(pos_fg)

# ── Plot helper function ───────────────────────────────────────────────────────
make_network_plot <- function(node_df, edge_df, pie_cols, fill_colors,
                              legend_title, title_text, subtitle_text) {
  ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
                 color = "grey50", linewidth = 0.4, alpha = 0.7) +
    geom_text(data = edge_df %>% filter(steps <= 5),
              aes(x = (x_from + x_to) / 2, y = (y_from + y_to) / 2, label = steps),
              size = 1.8, color = "grey30") +
    geom_scatterpie(data = node_df,
                    aes(x = x, y = y, r = r),
                    cols = pie_cols,
                    color = "white", linewidth = 0.2, alpha = 0.95) +
    scale_fill_manual(values = fill_colors, name = legend_title) +
    coord_equal() +
    theme_void(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position  = "right",
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.4, "cm"),
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle    = element_text(size = 8, hjust = 0.5, color = "grey40"),
      plot.caption     = element_text(size = 7, color = "grey50"),
      plot.margin      = margin(10, 10, 10, 10)
    ) +
    labs(
      title    = title_text,
      subtitle = subtitle_text,
      caption  = "Uniform node size | Edge labels = mutational steps (≤5 shown)"
    ) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 4)))
}

# ── Plot 1: All samples (outgroup + Dicorynia paranensis included) ─────────────
cat("INFO: generating all-samples network plot...\n")

node_all_base <- pos_all %>%
  bind_cols(as.data.frame(pie_sites)) %>%
  mutate(n = node_sizes,
         is_singleton = n == 1,
         r = ifelse(is_singleton, r_singleton, r_pie))

node_all_pies <- node_all_base %>% filter(!is_singleton)
node_all_sing <- node_all_base %>% filter(is_singleton)

p_all <- make_network_plot(
  node_df      = node_all_pies,
  edge_df      = edges_all,
  pie_cols     = all_sites,
  fill_colors  = site_colors,
  legend_title = "Sampling site",
  title_text   = expression(italic("Dicorynia guianensis") ~ "— Chloroplast haplotype network"),
  subtitle_text = paste0("All samples | n=", nrow(meta), " individuals | ", n_haps, " haplotypes | ",
                         nrow(node_all_sing), " singletons as grey dots")
) +
  geom_point(data = node_all_sing,
             aes(x = x, y = y), shape = 21,
             fill = "grey60", color = "white", size = r_singleton * 80, stroke = 0.3)

out_all <- file.path(output_dir, "cp_hapnet_all_samples.png")
ggsave(out_all, p_all, width = 16, height = 14, dpi = 300, device = "png")
cat("INFO: saved", out_all, "\n")

# ── Repulsion algorithm ────────────────────────────────────────────────────────
# Iteratively pushes overlapping nodes apart while preserving relative layout.
# Returns repelled x/y; original positions kept for connector lines.
repel_positions <- function(x, y, r, max_iter = 2000, padding = 0.05) {
  n  <- length(x)
  px <- x
  py <- y
  for (iter in seq_len(max_iter)) {
    moved <- FALSE
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1, n)) {
        dx     <- px[j] - px[i]
        dy     <- py[j] - py[i]
        dist_ij <- sqrt(dx^2 + dy^2)
        min_d  <- r[i] + r[j] + padding
        if (dist_ij < min_d && dist_ij > 1e-10) {
          push   <- (min_d - dist_ij) / 2
          angle  <- atan2(dy, dx)
          px[i]  <- px[i] - push * cos(angle)
          py[i]  <- py[i] - push * sin(angle)
          px[j]  <- px[j] + push * cos(angle)
          py[j]  <- py[j] + push * sin(angle)
          moved  <- TRUE
        }
      }
    }
    if (!moved) break
  }
  cat("INFO: repulsion converged after", iter, "iterations\n")
  list(x = px, y = py)
}

# ── Plot 2: French Guiana only (no outgroup, no Dicorynia paranensis) ──────────
cat("INFO: generating French Guiana-only network plot...\n")

# Sites for FG-only plot: exclude Cameroun_Benin and Dicorynia paranensis
sites_fg_only  <- sort(unique(meta$site[meta$region == "French Guiana" &
                                          meta$site != "Dicorynia paranensis"]))
site_colors_fg <- site_colors[sites_fg_only]

# Node data: keep only FG site columns, exclude haplotypes with 0 FG individuals
node_fg <- pos_fg %>%
  bind_cols(as.data.frame(pie_sites[, sites_fg_only, drop = FALSE])) %>%
  mutate(
    n_fg        = rowSums(across(all_of(sites_fg_only))),
    n           = node_sizes,
    is_singleton = n_fg == 1,
    r           = ifelse(is_singleton, r_singleton, r_pie)
  ) %>%
  filter(n_fg > 0)

n_fg_haps <- nrow(node_fg)
n_fg_ind  <- sum(node_fg$n_fg)
cat("INFO: French Guiana —", n_fg_haps, "haplotypes,", n_fg_ind, "individuals\n")

# Apply repulsion to FG nodes — store original positions for connector lines
cat("INFO: applying repulsion to FG nodes...\n")
repelled <- repel_positions(node_fg$x, node_fg$y, node_fg$r)
node_fg <- node_fg %>%
  mutate(
    x_orig = x,
    y_orig = y,
    x      = repelled$x,
    y      = repelled$y,
    moved  = sqrt((x - x_orig)^2 + (y - y_orig)^2) > (r_fixed * 0.1)
  )

# Connector lines: only for nodes that were actually moved
connectors_fg <- node_fg %>% filter(moved) %>% select(x_orig, y_orig, x, y)

# Edges: keep only edges where both endpoints have FG individuals
# Built from scratch using repelled node positions
fg_hap_ids <- node_fg$hap_id
edges_fg_plot <- data.frame(from = edge_from, to = edge_to, steps = edge_step) %>%
  filter(from %in% fg_hap_ids & to %in% fg_hap_ids) %>%
  left_join(node_fg %>% select(hap_id, x, y), by = c("from" = "hap_id")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(node_fg %>% select(hap_id, x, y), by = c("to" = "hap_id")) %>%
  rename(x_to = x, y_to = y)

p_fg <- ggplot() +
  # Connector lines from original to repelled position
  geom_segment(data = connectors_fg,
               aes(x = x_orig, y = y_orig, xend = x, yend = y),
               color = "grey70", linewidth = 0.3, linetype = "dashed") +
  # Network edges (at repelled positions)
  geom_segment(data = edges_fg_plot,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
               color = "grey50", linewidth = 0.4, alpha = 0.7) +
  geom_text(data = edges_fg_plot %>% filter(steps <= 5),
            aes(x = (x_from + x_to) / 2, y = (y_from + y_to) / 2, label = steps),
            size = 1.8, color = "grey30") +
  # Singleton nodes: small grey dots
  geom_point(data = node_fg %>% filter(is_singleton),
             aes(x = x, y = y), shape = 21,
             fill = "grey60", color = "white", size = r_singleton * 80, stroke = 0.3) +
  # Pie chart nodes for haplotypes with n >= 2
  geom_scatterpie(data = node_fg %>% filter(!is_singleton),
                  aes(x = x, y = y, r = r),
                  cols = sites_fg_only,
                  color = "white", linewidth = 0.2, alpha = 0.95) +
  scale_fill_manual(values = site_colors_fg, name = "Sampling site") +
  coord_equal() +
  theme_void(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "right",
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.4, "cm"),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle    = element_text(size = 8, hjust = 0.5, color = "grey40"),
    plot.caption     = element_text(size = 7, color = "grey50"),
    plot.margin      = margin(10, 10, 10, 10)
  ) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast haplotype network (French Guiana)"),
    subtitle = paste0("Outgroup & Dicorynia paranensis excluded | ", n_fg_haps,
                      " haplotypes | ", n_fg_ind, " individuals | singletons = grey dots"),
    caption  = "Pie charts: n≥2 | Grey dots: singletons (n=1) | Dashed lines = repelled nodes | Edge labels = mutational steps (≤5)"
  ) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 4)))

out_fg <- file.path(output_dir, "cp_hapnet_fg_only.png")
ggsave(out_fg, p_fg, width = 16, height = 14, dpi = 300, device = "png")
cat("INFO: saved", out_fg, "\n")

# ── Export summary table ───────────────────────────────────────────────────────
hap_summary <- data.frame(
  haplotype_id  = seq_len(n_haps),
  n_individuals = node_sizes,
  x_mds = pos$x,
  y_mds = pos$y,
  sites = sapply(hap_membership, function(idx) {
    samples <- seq_names[idx]
    paste(sort(unique(sample_info$site[sample_info$sample_id %in% samples])), collapse = ";")
  }),
  stringsAsFactors = FALSE
)

out_table <- file.path(output_dir, "cp_hapnet_summary.tsv")
write.table(hap_summary, file = out_table, sep = "\t", row.names = FALSE, quote = FALSE)
cat("INFO: saved", out_table, "\n")

cat("DONE haplotype network complete\n")

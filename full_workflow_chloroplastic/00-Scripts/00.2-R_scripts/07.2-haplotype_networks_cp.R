#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.2-haplotype_networks_cp.R
# Description : Chloroplast haplotype network for Dicorynia guianensis.
#               Network topology from pegas::haploNet (MSN).
#               Layout via MDS on pairwise genetic distances.
#               Nodes = pie charts (scatterpie) sized by haplotype frequency.
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
  mutate(region = ifelse(site == "Cameroun_Benin", "Outgroup (Cameroun/Benin)", "French Guiana"))

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
cat("INFO: computing MDS layout...\n")
d_mat <- as.matrix(dist.dna(haps, model = "N", pairwise.deletion = TRUE))
mds   <- cmdscale(d_mat, k = 2)
pos   <- data.frame(hap_id = seq_len(n_haps), x = mds[, 1], y = mds[, 2])

# ── Site/region composition per haplotype ────────────────────────────────────
hap_membership <- attr(haps, "index")
seq_names      <- rownames(aln)

all_sites   <- sort(unique(meta$site))
all_regions <- c("French Guiana", "Outgroup (Cameroun/Benin)")

pie_sites   <- matrix(0, nrow = n_haps, ncol = length(all_sites),
                      dimnames = list(NULL, all_sites))
pie_regions <- matrix(0, nrow = n_haps, ncol = length(all_regions),
                      dimnames = list(NULL, all_regions))

for (i in seq_along(hap_membership)) {
  samples <- seq_names[hap_membership[[i]]]
  info    <- sample_info %>% filter(sample_id %in% samples)
  for (s in info$site)   pie_sites[i, s]   <- pie_sites[i, s]   + 1
  for (r in info$region) pie_regions[i, r] <- pie_regions[i, r] + 1
}

node_sizes <- sapply(hap_membership, length)

# ── Color palettes ─────────────────────────────────────────────────────────────
sites_fg    <- sort(unique(meta$site[meta$region == "French Guiana"]))
pal_fg      <- colorRampPalette(brewer.pal(12, "Paired"))(length(sites_fg))
site_colors <- setNames(c(pal_fg, "#E41A1C"), c(sites_fg, "Cameroun_Benin"))

region_colors <- c("French Guiana" = "#2166AC", "Outgroup (Cameroun/Benin)" = "#E41A1C")

# ── Build node data frames ─────────────────────────────────────────────────────
# Pie radius: sqrt(n) scaled to plot coordinates
r_scale <- diff(range(pos$x)) * 0.03

node_sites <- pos %>%
  bind_cols(as.data.frame(pie_sites)) %>%
  mutate(n = node_sizes, r = pmax(sqrt(n) * r_scale, r_scale * 0.6))

node_regions <- pos %>%
  bind_cols(as.data.frame(pie_regions)) %>%
  mutate(n = node_sizes, r = pmax(sqrt(n) * r_scale, r_scale * 0.6))

# ── Build edge data frame ─────────────────────────────────────────────────────
edges <- data.frame(from = edge_from, to = edge_to, steps = edge_step) %>%
  left_join(pos, by = c("from" = "hap_id")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(pos, by = c("to" = "hap_id")) %>%
  rename(x_to = x, y_to = y)

# ── Plot helper function ───────────────────────────────────────────────────────
make_network_plot <- function(node_df, pie_cols, fill_colors, legend_title, subtitle_text) {
  ggplot() +
    # Edges
    geom_segment(data = edges,
                 aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
                 color = "grey50", linewidth = 0.4, alpha = 0.7) +
    # Edge step labels (only for edges with few steps, to avoid clutter)
    geom_text(data = edges %>% filter(steps <= 5),
              aes(x = (x_from + x_to) / 2, y = (y_from + y_to) / 2, label = steps),
              size = 1.8, color = "grey30") +
    # Pie chart nodes
    geom_scatterpie(data = node_df,
                    aes(x = x, y = y, r = r),
                    cols = pie_cols,
                    color = "white", linewidth = 0.2,
                    alpha = 0.95) +
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
      title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast haplotype network"),
      subtitle = subtitle_text,
      caption  = "Node size ∝ haplotype frequency | Edge labels = mutational steps (≤5 shown)"
    ) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 4)))
}

# ── Plot 1: by sampling site ──────────────────────────────────────────────────
cat("INFO: generating site-level network plot...\n")
p_sites <- make_network_plot(
  node_df      = node_sites,
  pie_cols     = all_sites,
  fill_colors  = site_colors,
  legend_title = "Sampling site",
  subtitle_text = paste0("Minimum spanning network | MDS layout | n=", nrow(meta), " individuals | ",
                         n_haps, " haplotypes")
)

out_sites <- file.path(output_dir, "cp_hapnet_sites.png")
ggsave(out_sites, p_sites, width = 16, height = 14, dpi = 300, device = "png")
cat("INFO: saved", out_sites, "\n")

# ── Plot 2: by region ─────────────────────────────────────────────────────────
cat("INFO: generating region-level network plot...\n")
p_regions <- make_network_plot(
  node_df      = node_regions,
  pie_cols     = all_regions,
  fill_colors  = region_colors,
  legend_title = "Region",
  subtitle_text = "French Guiana vs. Cameroun/Benin outgroup"
)

out_regions <- file.path(output_dir, "cp_hapnet_regions.png")
ggsave(out_regions, p_regions, width = 14, height = 12, dpi = 300, device = "png")
cat("INFO: saved", out_regions, "\n")

# ── Plot 3: French Guiana only — by site (outgroup excluded) ──────────────────
cat("INFO: generating French Guiana-only site network plot...\n")

# Sites excluding outgroup
sites_fg_only <- sort(unique(meta$site[meta$region == "French Guiana"]))
site_colors_fg <- site_colors[sites_fg_only]

# Node data: keep only FG site columns, exclude haplotypes with 0 FG individuals
node_fg <- pos %>%
  bind_cols(as.data.frame(pie_sites[, sites_fg_only, drop = FALSE])) %>%
  mutate(
    n_fg = rowSums(across(all_of(sites_fg_only))),
    n    = node_sizes,
    r    = pmax(sqrt(n_fg) * r_scale, r_scale * 0.6)
  ) %>%
  filter(n_fg > 0)  # remove outgroup-only haplotypes

n_fg_haps <- nrow(node_fg)
n_fg_ind  <- sum(node_fg$n_fg)
cat("INFO: French Guiana —", n_fg_haps, "haplotypes,", n_fg_ind, "individuals\n")

# Edges: keep only edges where both endpoints have FG individuals
fg_hap_ids <- node_fg$hap_id
edges_fg <- edges %>%
  filter(from %in% fg_hap_ids & to %in% fg_hap_ids)

p_fg <- ggplot() +
  geom_segment(data = edges_fg,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
               color = "grey50", linewidth = 0.4, alpha = 0.7) +
  geom_text(data = edges_fg %>% filter(steps <= 5),
            aes(x = (x_from + x_to) / 2, y = (y_from + y_to) / 2, label = steps),
            size = 1.8, color = "grey30") +
  geom_scatterpie(data = node_fg,
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
    subtitle = paste0("Outgroup excluded | ", n_fg_haps, " haplotypes | ", n_fg_ind, " individuals | 21 sites"),
    caption  = "Node size ∝ haplotype frequency | Edge labels = mutational steps (≤5 shown)"
  ) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 4)))

out_fg <- file.path(output_dir, "cp_hapnet_fg_only.png")
ggsave(out_fg, p_fg, width = 16, height = 14, dpi = 300, device = "png")
cat("INFO: saved", out_fg, "\n")

# ── Plot 4: French Guiana only, Treemutation samples highlighted ───────────────
cat("INFO: generating Treemutation-highlighted network plot...\n")

treemut_samples <- meta %>%
  filter(project == "Treemutation") %>%
  pull(sample_id)

cat("INFO: Treemutation samples —", length(treemut_samples), ":",
    paste(treemut_samples, collapse = ", "), "\n")

# For each haplotype, count how many Treemutation individuals it contains
n_treemut_per_hap <- sapply(hap_membership, function(idx) {
  samples <- seq_names[idx]
  sum(samples %in% treemut_samples)
})

# Build node data: FG only + Treemutation flag
node_treemut <- node_fg %>%
  mutate(
    n_treemut   = n_treemut_per_hap[hap_id],
    has_treemut = n_treemut > 0,
    border_col  = ifelse(has_treemut, "#FF6600", "white"),
    border_size = ifelse(has_treemut, 1.2, 0.2)
  )

# Haplotypes carrying only Treemutation individuals (private)
treemut_only_ids <- node_treemut %>% filter(has_treemut & n_fg == n_treemut) %>% pull(hap_id)

p_treemut <- ggplot() +
  geom_segment(data = edges_fg,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
               color = "grey50", linewidth = 0.4, alpha = 0.7) +
  geom_text(data = edges_fg %>% filter(steps <= 5),
            aes(x = (x_from + x_to) / 2, y = (y_from + y_to) / 2, label = steps),
            size = 1.8, color = "grey30") +
  # All FG nodes (no Treemutation border) — drawn first
  geom_scatterpie(data = node_treemut %>% filter(!has_treemut),
                  aes(x = x, y = y, r = r),
                  cols = sites_fg_only,
                  color = "white", linewidth = 0.2, alpha = 0.95) +
  # Treemutation nodes — drawn on top with orange border
  geom_scatterpie(data = node_treemut %>% filter(has_treemut),
                  aes(x = x, y = y, r = r),
                  cols = sites_fg_only,
                  color = "#FF6600", linewidth = 1.2, alpha = 0.95) +
  scale_fill_manual(values = site_colors_fg, name = "Sampling site") +
  # Dummy layer for legend entry for Treemutation border
  geom_point(data = data.frame(x = NA_real_, y = NA_real_),
             aes(x = x, y = y, shape = "Treemutation individuals"),
             color = "#FF6600", size = 4, stroke = 1.5) +
  scale_shape_manual(values = c("Treemutation individuals" = 1), name = NULL) +
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
    subtitle = paste0("Outgroup excluded | Treemutation individuals highlighted (orange border) | n=",
                      length(treemut_samples), " Treemutation samples"),
    caption  = "Node size ∝ haplotype frequency | Edge labels = mutational steps (≤5 shown)"
  ) +
  guides(
    fill  = guide_legend(ncol = 1, override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4, stroke = 1.5))
  )

out_treemut <- file.path(output_dir, "cp_hapnet_treemutation.png")
ggsave(out_treemut, p_treemut, width = 16, height = 14, dpi = 300, device = "png")
cat("INFO: saved", out_treemut, "\n")

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

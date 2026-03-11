#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.1-iqtree_viz_cp.R
# Description : Visualisation of the chloroplast ML phylogenetic tree
#               produced by IQ-TREE3. Tips colored by sampling site.
#               Node support displayed as SH-aLRT/UFBoot values.
# Author  : Julien Bonnier
# Usage   : Rscript --vanilla 07.1-iqtree_viz_cp.R <treefile> <metadata_csv> <output_dir>
# =============================================================================

# Ensure personal library is accessible when run with --vanilla
.libPaths(c(path.expand("~/work/R"), .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 07.1-iqtree_viz_cp.R <treefile> <metadata_csv> <output_dir>")
}
treefile   <- args[1]
meta_file  <- args[2]
output_dir <- args[3]

cat("INFO: treefile  =", treefile, "\n")
cat("INFO: metadata  =", meta_file, "\n")
cat("INFO: output_dir=", output_dir, "\n")

if (!file.exists(treefile))  stop("treefile not found: ", treefile)
if (!file.exists(meta_file)) stop("metadata not found: ", meta_file)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Read tree ─────────────────────────────────────────────────────────────────
tree <- read.tree(treefile)
cat("INFO: tree loaded —", Ntip(tree), "tips,", Nnode(tree), "internal nodes\n")

# ── Parse node support labels (SH-aLRT/UFBoot format) ─────────────────────────
# IQ-TREE with -alrt + -B writes: "85.3/96" on internal nodes
# Split into two separate columns for flexible plotting
node_labels <- tree$node.label
support_df <- data.frame(
  node    = (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)),
  label   = node_labels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    alrt   = as.numeric(sub("/.*", "", label)),
    ufboot = as.numeric(sub(".*/", "", label)),
    # Well-supported node: SH-aLRT >= 80 AND UFBoot >= 95
    supported = !is.na(alrt) & !is.na(ufboot) & alrt >= 80 & ufboot >= 95
  )

# ── Read metadata ──────────────────────────────────────────────────────────────
meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(sample_id, site) %>%
  distinct() %>%
  mutate(
    region = ifelse(site == "Cameroun_Benin", "Outgroup (Cameroun/Benin)", "French Guiana")
  )

cat("INFO: metadata loaded —", nrow(meta), "samples,", n_distinct(meta$site), "sites\n")

# ── Color palette ──────────────────────────────────────────────────────────────
# 21 French Guiana sites + 1 outgroup = 22 distinct colors
sites_fg  <- sort(unique(meta$site[meta$region == "French Guiana"]))
sites_out <- "Cameroun_Benin"
all_sites <- c(sites_fg, sites_out)

# Use a mix of color palettes for 22 distinct colors
pal_fg  <- colorRampPalette(brewer.pal(12, "Paired"))(length(sites_fg))
pal_out <- "#E41A1C"  # Red for outgroup — visually distinctive
site_colors <- setNames(c(pal_fg, pal_out), all_sites)

# ── Annotate tip metadata ──────────────────────────────────────────────────────
tip_data <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(meta, by = c("label" = "sample_id"))

# ── Join support data to tree via %<+% ────────────────────────────────────────
# This ensures isTip column is available internally when geom_nodelab is called

# ── Plot 1: Rectangular tree ───────────────────────────────────────────────────
p_rect <- ggtree(tree, layout = "rectangular", linewidth = 0.3) %<+% tip_data %<+% support_df +

  # Tip points colored by site
  geom_tippoint(aes(color = site), size = 1.5, alpha = 0.9) +

  # Node support: show UFBoot only at well-supported internal nodes
  geom_nodelab(
    aes(label = ifelse(!is.na(supported) & supported, as.character(ufboot), "")),
    size  = 1.8,
    hjust = 1.2,
    vjust = -0.4,
    color = "grey30"
  ) +

  # Outgroup annotation bar
  geom_strip(
    taxa1 = "DBR_AGN", taxa2 = "DBR_AGR",
    label = "Outgroup\n(Cameroun/Benin)",
    color = pal_out, fontsize = 2.5, offset = 0.0002, offset.text = 0.0001
  ) +

  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "right") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree"),
    subtitle = paste0("IQ-TREE3 | ModelFinder | UFBoot 1000 | SH-aLRT 1000 | n=",
                      Ntip(tree), " individuals"),
    caption  = "Node labels: UFBoot support (shown if SH-aLRT ≥ 80 & UFBoot ≥ 95)"
  ) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    plot.caption  = element_text(size = 7, color = "grey50"),
    legend.text   = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

# ── Plot 2: Circular (fan) tree ────────────────────────────────────────────────
p_circ <- ggtree(tree, layout = "circular", linewidth = 0.3) %<+% tip_data %<+% support_df +

  geom_tippoint(aes(color = site), size = 1.2, alpha = 0.9) +

  geom_nodelab(
    aes(label = ifelse(!is.na(supported) & supported, as.character(ufboot), "")),
    size  = 1.5,
    color = "grey30"
  ) +

  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "right") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (circular)"),
    subtitle = paste0("IQ-TREE3 | ModelFinder | UFBoot 1000 | SH-aLRT 1000 | n=",
                      Ntip(tree), " individuals")
  ) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    legend.text   = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

# ── Plot 3: Region-level tree (FG vs Outgroup) ─────────────────────────────────
region_colors <- c("French Guiana" = "#2166AC", "Outgroup (Cameroun/Benin)" = "#E41A1C")

p_region <- ggtree(tree, layout = "rectangular", linewidth = 0.3) %<+% tip_data %<+% support_df +

  geom_tippoint(aes(color = region), size = 1.5, alpha = 0.9) +

  geom_nodelab(
    aes(label = ifelse(!is.na(supported) & supported, as.character(ufboot), "")),
    size  = 1.8,
    hjust = 1.2,
    vjust = -0.4,
    color = "grey30"
  ) +

  scale_color_manual(values = region_colors, name = "Region") +
  theme_tree2(legend.position = "right") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (by region)"),
    subtitle = "French Guiana vs. Cameroun/Benin outgroup"
  ) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title  = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9)
  )

# ── Save outputs ───────────────────────────────────────────────────────────────
out_rect   <- file.path(output_dir, "cp_phylo_tree_rectangular.png")
out_circ   <- file.path(output_dir, "cp_phylo_tree_circular.png")
out_region <- file.path(output_dir, "cp_phylo_tree_region.png")

ggsave(out_rect,   p_rect,   width = 14, height = 28, units = "in", dpi = 300, device = "png")
ggsave(out_circ,   p_circ,   width = 14, height = 14, units = "in", dpi = 300, device = "png")
ggsave(out_region, p_region, width = 14, height = 28, units = "in", dpi = 300, device = "png")

cat("DONE tree plots saved:\n")
cat("  ", out_rect,   "\n")
cat("  ", out_circ,   "\n")
cat("  ", out_region, "\n")

# ── Export tree summary table ──────────────────────────────────────────────────
summary_file <- file.path(output_dir, "cp_phylo_node_support.tsv")
write.table(
  support_df %>% select(node, alrt, ufboot, supported),
  file      = summary_file,
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
cat("  ", summary_file, "\n")

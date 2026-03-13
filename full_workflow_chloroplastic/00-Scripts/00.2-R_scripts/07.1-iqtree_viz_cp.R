#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.1-iqtree_viz_cp.R
# Description : Visualisation of the chloroplast ML phylogenetic tree
#               produced by IQ-TREE3. Site colors from external CSV.
#               Seven plots:
#                 1. All samples, rectangular, colored by site
#                 2. All samples, rectangular, sqrt branch lengths
#                 3. All samples, circular, colored by site
#                 4. All samples, circular, sqrt branch lengths
#                 5. French Guiana only (no outgroups), rectangular
#                 6. Site-level tree, rectangular (true branch lengths)
#                 7. Site-level tree, rectangular, cladogram (topology only)
#                 8. Site-level tree, circular (true branch lengths)
#                 9. Site-level tree, circular, cladogram (topology only)
# Author  : Julien Bonnier
# Usage   : Rscript --vanilla 07.1-iqtree_viz_cp.R \
#               <treefile> <metadata_csv> <colors_csv> <output_dir>
# =============================================================================

.libPaths(c(path.expand("~/work/R"), .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 07.1-iqtree_viz_cp.R <treefile> <metadata_csv> <colors_csv> <output_dir>")
}
treefile    <- args[1]
meta_file   <- args[2]
colors_file <- args[3]
output_dir  <- args[4]

cat("INFO: treefile   =", treefile,    "\n")
cat("INFO: metadata   =", meta_file,   "\n")
cat("INFO: colors     =", colors_file, "\n")
cat("INFO: output_dir =", output_dir,  "\n")

if (!file.exists(treefile))    stop("treefile not found: ",    treefile)
if (!file.exists(meta_file))   stop("metadata not found: ",   meta_file)
if (!file.exists(colors_file)) stop("colors CSV not found: ", colors_file)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Read tree ─────────────────────────────────────────────────────────────────
tree <- read.tree(treefile)
cat("INFO: tree loaded —", Ntip(tree), "tips,", Nnode(tree), "internal nodes\n")

# ── Parse node support labels (SH-aLRT/UFBoot format) ─────────────────────────
node_labels <- tree$node.label
support_df <- data.frame(
  node  = (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)),
  label = node_labels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    alrt   = as.numeric(sub("/.*", "", label)),
    ufboot = as.numeric(sub(".*/", "", label)),
    supported = !is.na(alrt) & !is.na(ufboot) & alrt >= 80 & ufboot >= 95
  )

# ── Read metadata ──────────────────────────────────────────────────────────────
meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(sample_id, site, project) %>%
  distinct() %>%
  mutate(
    # Recode Herbier to D.paranensis; Treemutation already Angela in CSV
    site = case_when(
      site == "Herbier" ~ "D.paranensis",
      TRUE              ~ site
    ),
    region = case_when(
      site == "Cameroun_Benin" ~ "Outgroup",
      site == "D.paranensis"   ~ "Outgroup",
      TRUE                     ~ "French Guiana"
    )
  )

cat("INFO: metadata loaded —", nrow(meta), "samples,", n_distinct(meta$site), "sites\n")

# ── Read site colors from CSV ─────────────────────────────────────────────────
# Normalize CSV site keys to match metadata site names
color_raw <- read.csv(colors_file, stringsAsFactors = FALSE)

color_df <- color_raw %>%
  mutate(site_meta = case_when(
    site == "Trinite"        ~ "Trinité",
    site == "Petit_Croissant" ~ "Petit_croissant",
    site == "D.benthamianus"  ~ "Cameroun_Benin",
    TRUE                      ~ site
  ))

site_colors <- setNames(color_df$couleur_hex, color_df$site_meta)

# Verify all metadata sites have a color
missing_colors <- setdiff(unique(meta$site), names(site_colors))
if (length(missing_colors) > 0) {
  cat("WARNING: no color defined for sites:", paste(missing_colors, collapse = ", "), "\n")
  # Assign grey to missing sites
  site_colors[missing_colors] <- "#AAAAAA"
}

cat("INFO: site colors loaded —", nrow(color_df), "sites\n")

# ── Annotate tip metadata ──────────────────────────────────────────────────────
tip_data <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(meta, by = c("label" = "sample_id"))

# ── Common theme ──────────────────────────────────────────────────────────────
tree_theme <- theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA),
  plot.title       = element_text(face = "bold", size = 11),
  plot.subtitle    = element_text(size = 8, color = "grey40"),
  plot.caption     = element_text(size = 7, color = "grey50"),
  legend.text      = element_text(size = 7),
  legend.key.size  = unit(0.4, "cm")
)

node_label_layer <- function() {
  geom_nodelab(
    aes(label = ifelse(!is.na(supported) & supported, as.character(ufboot), "")),
    size  = 1.8, hjust = 1.2, vjust = -0.4, color = "grey30"
  )
}

# ── Plot 1: Rectangular — all samples ─────────────────────────────────────────
cat("INFO: generating plot 1 — all samples, rectangular...\n")

p_all_rect <- ggtree(tree, layout = "rectangular", linewidth = 0.3) %<+%
  tip_data %<+% support_df +
  geom_tippoint(aes(color = site), size = 1.2, alpha = 0.9) +
  node_label_layer() +
  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "right") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree"),
    subtitle = paste0("IQ-TREE3 | K3Pu+F+R3 | UFBoot 1000 | SH-aLRT 1000 | n=",
                      Ntip(tree), " individuals"),
    caption  = "Node labels: UFBoot (shown if SH-aLRT ≥ 80 & UFBoot ≥ 95)"
  ) +
  tree_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(output_dir, "cp_phylo_all_rect.png"),
       p_all_rect, width = 14, height = 28, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_all_rect.png\n")

# ── Plot 2: Circular — all samples ────────────────────────────────────────────
cat("INFO: generating plot 2 — all samples, circular...\n")

p_all_circ <- ggtree(tree, layout = "circular", linewidth = 0.3) %<+%
  tip_data %<+% support_df +
  geom_tippoint(aes(color = site), size = 1.2, alpha = 0.9) +
  geom_nodelab(
    aes(label = ifelse(!is.na(supported) & supported, as.character(ufboot), "")),
    size = 1.5, color = "grey30"
  ) +
  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "right") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (circular)"),
    subtitle = paste0("IQ-TREE3 | K3Pu+F+R3 | UFBoot 1000 | n=", Ntip(tree), " individuals")
  ) +
  tree_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(output_dir, "cp_phylo_all_circ.png"),
       p_all_circ, width = 14, height = 14, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_all_circ.png\n")

# ── Plots 2b & 3b: sqrt-transformed branch lengths — all samples ──────────────
cat("INFO: generating plots 2b & 3b — all samples, sqrt branch lengths...\n")

tree_sqrt <- tree
tree_sqrt$edge.length <- sqrt(tree$edge.length)

p_all_rect_sqrt <- ggtree(tree_sqrt, layout = "rectangular", linewidth = 0.3) %<+%
  tip_data %<+% support_df +
  geom_tippoint(aes(color = site), size = 1.2, alpha = 0.9) +
  node_label_layer() +
  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "right") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (sqrt branch lengths)"),
    subtitle = paste0("IQ-TREE3 | sqrt-transformed branch lengths | n=", Ntip(tree), " individuals"),
    caption  = "Branch lengths = sqrt(ML distances) — short branches expanded for readability"
  ) +
  tree_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(output_dir, "cp_phylo_all_rect_sqrt.png"),
       p_all_rect_sqrt, width = 14, height = 28, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_all_rect_sqrt.png\n")

p_all_circ_sqrt <- ggtree(tree_sqrt, layout = "circular", linewidth = 0.3) %<+%
  tip_data %<+% support_df +
  geom_tippoint(aes(color = site), size = 1.2, alpha = 0.9) +
  geom_nodelab(
    aes(label = ifelse(!is.na(supported) & supported, as.character(ufboot), "")),
    size = 1.5, color = "grey30"
  ) +
  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "right") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (circular, sqrt branch lengths)"),
    subtitle = paste0("IQ-TREE3 | sqrt-transformed branch lengths | n=", Ntip(tree), " individuals"),
    caption  = "Branch lengths = sqrt(ML distances) — short branches expanded for readability"
  ) +
  tree_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(output_dir, "cp_phylo_all_circ_sqrt.png"),
       p_all_circ_sqrt, width = 14, height = 14, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_all_circ_sqrt.png\n")

# ── Plot 3: Rectangular — French Guiana only (no outgroups) ───────────────────
cat("INFO: generating plot 3 — French Guiana only...\n")

outgroup_tips <- meta$sample_id[meta$region == "Outgroup"]
outgroup_tips <- outgroup_tips[outgroup_tips %in% tree$tip.label]
cat("INFO: dropping", length(outgroup_tips), "outgroup tips\n")

tree_fg <- drop.tip(tree, outgroup_tips)
tip_data_fg <- data.frame(label = tree_fg$tip.label, stringsAsFactors = FALSE) %>%
  left_join(meta, by = c("label" = "sample_id"))

# Support df for FG tree (same node indices after drop.tip may shift)
support_df_fg <- data.frame(
  node  = (Ntip(tree_fg) + 1):(Ntip(tree_fg) + Nnode(tree_fg)),
  label = tree_fg$node.label,
  stringsAsFactors = FALSE
) %>%
  mutate(
    alrt      = as.numeric(sub("/.*", "", label)),
    ufboot    = as.numeric(sub(".*/", "", label)),
    supported = !is.na(alrt) & !is.na(ufboot) & alrt >= 80 & ufboot >= 95
  )

site_colors_fg <- site_colors[names(site_colors) != "Cameroun_Benin" &
                                names(site_colors) != "D.paranensis"]

p_fg_rect <- ggtree(tree_fg, layout = "rectangular", linewidth = 0.3) %<+%
  tip_data_fg %<+% support_df_fg +
  geom_tippoint(aes(color = site), size = 1.2, alpha = 0.9) +
  node_label_layer() +
  scale_color_manual(values = site_colors_fg, name = "Sampling site") +
  theme_tree2(legend.position = "right") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (French Guiana)"),
    subtitle = paste0("Outgroups excluded | n=", Ntip(tree_fg), " individuals"),
    caption  = "Node labels: UFBoot (shown if SH-aLRT ≥ 80 & UFBoot ≥ 95)"
  ) +
  tree_theme +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(output_dir, "cp_phylo_fg_rect.png"),
       p_fg_rect, width = 14, height = 26, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_fg_rect.png\n")

# ── Plot 4 & 5: Site-level tree — French Guiana only (one tip per FG site) ─────
cat("INFO: generating plots 4 & 5 — site-level FG only (rect + circ)...\n")

# Keep one representative per site (all sites including outgroups)
site_reps_fg <- meta %>%
  filter(sample_id %in% tree$tip.label,
         site != "Cameroun_Benin") %>%
  group_by(site) %>%
  slice(1) %>%
  ungroup()

tree_site_fg <- keep.tip(tree, site_reps_fg$sample_id)
tip_map_fg   <- setNames(site_reps_fg$site, site_reps_fg$sample_id)
tree_site_fg$tip.label <- tip_map_fg[tree_site_fg$tip.label]

# Support for site-level tree
support_df_site_fg <- data.frame(
  node  = (Ntip(tree_site_fg) + 1):(Ntip(tree_site_fg) + Nnode(tree_site_fg)),
  label = tree_site_fg$node.label,
  stringsAsFactors = FALSE
) %>%
  mutate(
    alrt      = as.numeric(sub("/.*", "", label)),
    ufboot    = as.numeric(sub(".*/", "", label)),
    supported = !is.na(alrt) & !is.na(ufboot) & alrt >= 80 & ufboot >= 95
  )

tip_data_site_fg <- data.frame(
  label = tree_site_fg$tip.label,
  site  = tree_site_fg$tip.label,
  stringsAsFactors = FALSE
)

# Plot 4: Rectangular
p_sites_rect <- ggtree(tree_site_fg, layout = "rectangular", linewidth = 0.5) %<+%
  tip_data_site_fg %<+% support_df_site_fg +
  geom_tippoint(aes(color = site), size = 3, alpha = 0.95) +
  geom_tiplab(aes(color = site), size = 3.5, offset = 0.00005, fontface = "bold") +
  node_label_layer() +
  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "none") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (site-level)"),
    subtitle = paste0("One representative per site | ", Ntip(tree_site_fg), " sites"),
    caption  = "Node labels: UFBoot (shown if SH-aLRT ≥ 80 & UFBoot ≥ 95)"
  ) +
  tree_theme

ggsave(file.path(output_dir, "cp_phylo_sites_rect.png"),
       p_sites_rect, width = 12, height = 10, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_sites_rect.png\n")

# Plot 5: Circular
p_sites_circ <- ggtree(tree_site_fg, layout = "circular", linewidth = 0.5) %<+%
  tip_data_site_fg %<+% support_df_site_fg +
  geom_tippoint(aes(color = site), size = 3, alpha = 0.95) +
  geom_tiplab(aes(color = site), size = 3.5, offset = 0.00005, fontface = "bold") +
  geom_nodelab(
    aes(label = ifelse(!is.na(supported) & supported, as.character(ufboot), "")),
    size = 2, color = "grey30"
  ) +
  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "none") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (site-level, circular)"),
    subtitle = paste0("One representative per site | ", Ntip(tree_site_fg), " sites")
  ) +
  tree_theme

ggsave(file.path(output_dir, "cp_phylo_sites_circ.png"),
       p_sites_circ, width = 12, height = 12, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_sites_circ.png\n")

# ── Plots 6b & 7b: Site-level cladograms (topology only) ──────────────────────
cat("INFO: generating plots 6b & 7b — site-level cladograms...\n")

# Plot 6b: Rectangular cladogram
p_sites_rect_clado <- ggtree(tree_site_fg, branch.length = "none",
                              layout = "rectangular", linewidth = 0.5) %<+%
  tip_data_site_fg %<+% support_df_site_fg +
  geom_tippoint(aes(color = site), size = 3, alpha = 0.95) +
  geom_tiplab(aes(color = site), size = 3.5, offset = 0.3, fontface = "bold") +
  node_label_layer() +
  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "none") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (site-level, cladogram)"),
    subtitle = paste0("Topology only — branch lengths not shown | ", Ntip(tree_site_fg), " sites"),
    caption  = "Node labels: UFBoot (shown if SH-aLRT ≥ 80 & UFBoot ≥ 95)"
  ) +
  tree_theme

ggsave(file.path(output_dir, "cp_phylo_sites_rect_clado.png"),
       p_sites_rect_clado, width = 12, height = 10, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_sites_rect_clado.png\n")

# Plot 7b: Circular cladogram
p_sites_circ_clado <- ggtree(tree_site_fg, branch.length = "none",
                              layout = "circular", linewidth = 0.5) %<+%
  tip_data_site_fg %<+% support_df_site_fg +
  geom_tippoint(aes(color = site), size = 3, alpha = 0.95) +
  geom_tiplab(aes(color = site), size = 3.5, offset = 0.3, fontface = "bold") +
  geom_nodelab(
    aes(label = ifelse(!is.na(supported) & supported, as.character(ufboot), "")),
    size = 2, color = "grey30"
  ) +
  scale_color_manual(values = site_colors, name = "Sampling site") +
  theme_tree2(legend.position = "none") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast ML tree (site-level, circular cladogram)"),
    subtitle = paste0("Topology only — branch lengths not shown | ", Ntip(tree_site_fg), " sites")
  ) +
  tree_theme

ggsave(file.path(output_dir, "cp_phylo_sites_circ_clado.png"),
       p_sites_circ_clado, width = 12, height = 12, units = "in", dpi = 300)
cat("INFO: saved cp_phylo_sites_circ_clado.png\n")

# ── Export node support table ──────────────────────────────────────────────────
write.table(
  support_df %>% select(node, alrt, ufboot, supported),
  file      = file.path(output_dir, "cp_phylo_node_support.tsv"),
  sep       = "\t", row.names = FALSE, quote = FALSE
)

cat("DONE tree plots saved to", output_dir, "\n")

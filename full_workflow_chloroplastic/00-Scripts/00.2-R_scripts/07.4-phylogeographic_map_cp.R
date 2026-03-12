#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.4-phylogeographic_map_cp.R
# Description : Phylogeographic map of Dicorynia guianensis chloroplast
#               haplotype distribution across sampling sites.
#               Pie charts at each site show haplotype composition.
#               Outgroup (Cameroun/Benin) excluded from maps.
#               Uses ggplot2::map_data() — no sf/GDAL dependency.
# Author  : Julien Bonnier
# Usage   : Rscript --vanilla 07.4-phylogeographic_map_cp.R \
#               <haplotype_table> <metadata_csv> <geoloc_csv> <output_dir>
# =============================================================================

.libPaths(c(path.expand("~/work/R"), .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scatterpie)
  library(RColorBrewer)
  library(ggrepel)
})

theme_set(theme_minimal(base_size = 11) + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: script.R <hap_table> <metadata> <geoloc> <outdir>")

hap_file    <- args[1]
meta_file   <- args[2]
geoloc_file <- args[3]
output_dir  <- args[4]

cat("INFO: hap_table    =", hap_file, "\n")
cat("INFO: metadata     =", meta_file, "\n")
cat("INFO: geoloc       =", geoloc_file, "\n")
cat("INFO: output_dir   =", output_dir, "\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Read data ─────────────────────────────────────────────────────────────────
hap_table <- read.table(hap_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  select(sample_id, haplotype_id)

meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(sample_id, site) %>%
  distinct() %>%
  mutate(region = ifelse(site == "Cameroun_Benin", "Outgroup", "French Guiana"))

geoloc <- read.csv(geoloc_file, stringsAsFactors = FALSE) %>%
  rename(site = Sites) %>%
  filter(!is.na(lat) & !is.na(long))

cat("INFO: sites with coordinates:", nrow(geoloc), "\n")

sample_info <- hap_table %>% left_join(meta, by = "sample_id")

# ── Base map from ggplot2 (uses maps package — no GDAL needed) ────────────────
world_map <- map_data("world")
# Focus on French Guiana bounding box
xlim_fg <- c(-55.2, -51.2)
ylim_fg <- c(1.8, 6.4)

# ── Haplotype frequencies per site ───────────────────────────────────────────
N_TOP_HAPS <- 20

hap_counts_global <- sample_info %>%
  count(haplotype_id, name = "n_global") %>%
  arrange(desc(n_global))

top_haps <- hap_counts_global %>%
  slice_head(n = N_TOP_HAPS) %>%
  pull(haplotype_id)

hap_site <- sample_info %>%
  filter(site != "Cameroun_Benin") %>%
  mutate(hap_group = ifelse(haplotype_id %in% top_haps, haplotype_id, "Other")) %>%
  count(site, hap_group) %>%
  group_by(site) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

all_hap_groups <- c(sort(top_haps), "Other")

hap_wide <- hap_site %>%
  select(site, hap_group, freq) %>%
  pivot_wider(names_from = hap_group, values_from = freq, values_fill = 0)

# Ensure all haplotype columns are present
for (h in setdiff(all_hap_groups, names(hap_wide))) hap_wide[[h]] <- 0
hap_wide <- hap_wide %>% select(site, all_of(all_hap_groups))

n_per_site <- sample_info %>%
  filter(site != "Cameroun_Benin") %>%
  count(site, name = "n_samples")

map_data_sites <- hap_wide %>%
  left_join(geoloc, by = "site") %>%
  left_join(n_per_site, by = "site") %>%
  filter(!is.na(lat))

cat("INFO: FG sites on map:", nrow(map_data_sites), "\n")

# ── Color palette ─────────────────────────────────────────────────────────────
n_top <- length(top_haps)
pal_haps <- colorRampPalette(c(
  brewer.pal(8, "Set1"), brewer.pal(8, "Set2"), brewer.pal(8, "Dark2")
))(n_top)
hap_colors <- setNames(c(pal_haps, "grey70"), c(sort(top_haps), "Other"))

# ── Fixed pie radius (avoids large pies in dense areas) ───────────────────────
r_fixed <- diff(xlim_fg) * 0.022   # uniform radius for all sites
map_data_sites <- map_data_sites %>% mutate(r = r_fixed)

# ── Repulsion layout: push overlapping pies apart ─────────────────────────────
# Iteratively separates pies that overlap; connector lines drawn to true coords
repel_pies <- function(x, y, r, max_iter = 800, step_frac = 0.5) {
  px <- x; py <- y
  n  <- length(x)
  for (iter in seq_len(max_iter)) {
    moved <- FALSE
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1, n)) {
        dx   <- px[j] - px[i]
        dy   <- py[j] - py[i]
        dist <- sqrt(dx^2 + dy^2)
        gap  <- r[i] + r[j] + 0.04    # minimum separation (degrees)
        if (dist < gap && dist > 1e-9) {
          push   <- (gap - dist) * step_frac
          px[i]  <- px[i] - push * dx / dist
          py[i]  <- py[i] - push * dy / dist
          px[j]  <- px[j] + push * dx / dist
          py[j]  <- py[j] + push * dy / dist
          moved  <- TRUE
        }
      }
    }
    if (!moved) break
  }
  list(x = px, y = py)
}

repelled <- repel_pies(map_data_sites$long, map_data_sites$lat,
                       rep(r_fixed, nrow(map_data_sites)))

map_data_sites <- map_data_sites %>%
  mutate(px = repelled$x, py = repelled$y,
         moved = sqrt((px - long)^2 + (py - lat)^2) > 0.01)

# Connector lines: true site coords → repelled pie center
connectors <- map_data_sites %>%
  filter(moved) %>%
  select(long, lat, px, py)

cat("INFO: repulsion displaced", sum(map_data_sites$moved), "sites\n")

# ── Plot 1: Pie chart map — haplotype composition per site ───────────────────
cat("INFO: generating haplotype composition map...\n")

p_fg <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "#F0EFE7", color = "grey55", linewidth = 0.3) +
  coord_fixed(ratio = 1, xlim = xlim_fg, ylim = ylim_fg) +
  # Connector lines from true site to repelled pie
  geom_segment(data = connectors,
               aes(x = long, y = lat, xend = px, yend = py),
               color = "grey45", linewidth = 0.35, linetype = "solid") +
  # True site anchor points
  geom_point(data = map_data_sites %>% filter(moved),
             aes(x = long, y = lat),
             shape = 3, size = 1.2, color = "grey30", stroke = 0.6) +
  geom_scatterpie(
    data = map_data_sites,
    aes(x = px, y = py, r = r),
    cols = all_hap_groups,
    color = "white", linewidth = 0.2, alpha = 0.95
  ) +
  geom_text_repel(
    data = map_data_sites,
    aes(x = px, y = py, label = paste0(site, "\n(n=", n_samples, ")")),
    size = 2.0, color = "grey15", max.overlaps = 40,
    box.padding = 0.25, segment.size = 0.25, segment.color = "grey60",
    min.segment.length = 0.1
  ) +
  scale_fill_manual(values = hap_colors, name = "Haplotype") +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast haplotype distribution"),
    subtitle = paste0("French Guiana | n=", sum(map_data_sites$n_samples),
                      " individuals | ", nrow(hap_counts_global),
                      " unique haplotypes | top ", N_TOP_HAPS, " shown individually"),
    x = "Longitude", y = "Latitude",
    caption = "Uniform pie size | lines connect overlapping pies to true site location (\u002b)"
  ) +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 6.5),
    legend.key.size  = unit(0.35, "cm"),
    panel.grid       = element_line(color = "grey88", linewidth = 0.2),
    panel.border     = element_rect(color = "grey60", fill = NA, linewidth = 0.4),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  ) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(output_dir, "cp_map_fg_haplotypes.png"),
       p_fg, width = 14, height = 12, dpi = 300)
cat("INFO: saved cp_map_fg_haplotypes.png\n")

# ── Per-site diversity for bubble map ────────────────────────────────────────
hap_site_stats <- hap_site %>%
  group_by(site) %>%
  summarise(
    Hd = {
      p  <- freq
      ns <- sum(n)
      if (ns < 2) NA_real_ else (1 - sum(p^2)) * ns / (ns - 1)
    },
    n_haplotypes = n_distinct(hap_group[n > 0]),
    n_samples    = sum(n),
    .groups = "drop"
  )

bubble_data <- hap_site_stats %>%
  left_join(geoloc, by = "site") %>%
  filter(!is.na(lat) & !is.na(Hd))

# ── Plot 2: Bubble map — Hd and N haplotypes per site ────────────────────────
cat("INFO: generating diversity bubble map...\n")

p_bubble <- ggplot() +
  geom_polygon(data = world_map,
               aes(x = long, y = lat, group = group),
               fill = "#F0EFE7", color = "grey55", linewidth = 0.3) +
  coord_fixed(ratio = 1, xlim = xlim_fg, ylim = ylim_fg) +
  geom_point(data = bubble_data,
             aes(x = long, y = lat, size = n_haplotypes, fill = Hd),
             shape = 21, color = "white", stroke = 0.6, alpha = 0.9) +
  geom_text_repel(
    data = bubble_data,
    aes(x = long, y = lat, label = site),
    size = 2.2, color = "grey15", max.overlaps = 25,
    box.padding = 0.3, segment.size = 0.3, segment.color = "grey50"
  ) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Haplotype\ndiversity (Hd)") +
  scale_size_continuous(name = "N haplotypes", range = c(3, 12)) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast diversity map"),
    subtitle = "French Guiana | bubble size = N haplotypes | color = Hd",
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    legend.position = "right",
    panel.grid      = element_line(color = "grey88", linewidth = 0.2),
    panel.border    = element_rect(color = "grey60", fill = NA, linewidth = 0.4),
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(output_dir, "cp_map_fg_diversity.png"),
       p_bubble, width = 12, height = 10, dpi = 300)
cat("INFO: saved cp_map_fg_diversity.png\n")

# ── Export site summary table ─────────────────────────────────────────────────
out_table <- file.path(output_dir, "cp_map_site_summary.tsv")
write.table(
  hap_site_stats %>%
    left_join(geoloc, by = "site") %>%
    arrange(desc(Hd)),
  file = out_table, sep = "\t", row.names = FALSE, quote = FALSE
)
cat("INFO: saved", out_table, "\n")

cat("DONE phylogeographic map complete\n")

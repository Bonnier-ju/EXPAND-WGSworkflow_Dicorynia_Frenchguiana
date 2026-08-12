#!/usr/bin/env Rscript
# 08.1b.3-topo_extraction.R
#
# Extract 2 topographic variables at the 19 French Guiana study sites from
# the Copernicus DEM GLO-30 tiles already downloaded for the project
# (08.1-env_variable_selection/data/elevation/copernicus_dem30/*.tif).
#
# Simplified to 2 variables (confirmed 2026-08-11): elevation + TRI (Terrain
# Ruggedness Index). TWI (topographic wetness index) dropped - no flow
# accumulation tool available on Genotoul without a new binary/package
# install (whitebox_tools/SAGA/GRASS all absent), and TRI/TWI are often
# correlated anyway so the information loss is limited.
#
# Input:
#   geoloc_site.csv          (metadata/, columns: Sites, lat, long)
#   copernicus_dem30/*.tif   (7 tiles, ~30 m resolution, EPSG:4326)
#
# Output:
#   data/topo_variables_19sites.csv   (site, elevation, tri)
#   data/topo_variables_summary.tsv
#   plots/topo_{elevation,tri}_per_site.png
#
# Args: dem_dir  geoloc_csv  out_dir

suppressPackageStartupMessages({
  ok_terra <- requireNamespace("terra", quietly = TRUE)
  library(ggplot2)
})
if (!ok_terra) {
  stop("Package 'terra' not available in R_LIBS_USER. Check with: search_R_package terra")
}

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else {
  warning("sites_couleurs.csv not found - site labels will be grey")
  c()
}

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
dem_dir   <- if (length(args) >= 1) args[1] else stop("dem_dir required")
geoloc_csv <- if (length(args) >= 2) args[2] else stop("geoloc_csv required")
out_dir   <- if (length(args) >= 3) args[3] else stop("out_dir required")

data_dir <- file.path(out_dir, "data")
plot_dir <- file.path(out_dir, "plots")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Guard clause: verify inputs before running ----
tiles <- list.files(dem_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(tiles) == 0) stop(sprintf("ERROR: no .tif tiles found in %s", dem_dir))
if (!file.exists(geoloc_csv) || file.info(geoloc_csv)$size == 0) {
  stop(sprintf("ERROR: required input missing or empty: %s", geoloc_csv))
}

cat(sprintf("Found %d DEM tiles in %s\n", length(tiles), dem_dir))

# -------------------------------------------------------------------
# STEP 1: Load sites and mosaic DEM tiles (virtual raster, lazy)
# -------------------------------------------------------------------
cat("STEP 1: loading sites and building DEM mosaic ...\n")

sites_df <- read.csv(geoloc_csv, stringsAsFactors = FALSE)
stopifnot(all(c("Sites", "lat", "long") %in% colnames(sites_df)))
cat(sprintf("  %d sites in %s\n", nrow(sites_df), geoloc_csv))

dem <- terra::vrt(tiles)
pts <- terra::vect(sites_df, geom = c("long", "lat"), crs = terra::crs(dem))

# -------------------------------------------------------------------
# STEP 2: Extract elevation
# -------------------------------------------------------------------
cat("STEP 2: extracting elevation ...\n")

elev_vals <- terra::extract(dem, pts)[, 2]

# -------------------------------------------------------------------
# STEP 3: Compute TRI over the mosaic and extract at sites
# -------------------------------------------------------------------
cat("STEP 3: computing TRI and extracting ...\n")

tri_rast <- terra::terrain(dem, v = "TRI")
tri_vals <- terra::extract(tri_rast, pts)[, 2]

topo_df <- data.frame(
  site      = sites_df$Sites,
  elevation = elev_vals,
  tri       = tri_vals,
  stringsAsFactors = FALSE
)
topo_df <- topo_df[order(topo_df$site), ]

n_na <- sum(is.na(topo_df$elevation) | is.na(topo_df$tri))
if (n_na > 0) {
  warning(sprintf("%d site(s) have NA elevation/tri - check DEM tile coverage:\n%s",
                  n_na, paste(topo_df$site[is.na(topo_df$elevation) | is.na(topo_df$tri)],
                              collapse = ", ")))
}

# -------------------------------------------------------------------
# STEP 4: Write outputs
# -------------------------------------------------------------------
cat("STEP 4: writing outputs ...\n")

out_csv <- file.path(data_dir, "topo_variables_19sites.csv")
write.csv(topo_df, out_csv, row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_csv))

topo_vars <- c("elevation", "tri")
summ_rows <- lapply(topo_vars, function(v) {
  vals <- topo_df[[v]]
  m <- mean(vals, na.rm = TRUE)
  s <- sd(vals, na.rm = TRUE)
  data.frame(
    variable = v,
    n_sites  = sum(!is.na(vals)),
    mean     = round(m, 4),
    sd       = round(s, 4),
    min      = round(min(vals, na.rm = TRUE), 4),
    max      = round(max(vals, na.rm = TRUE), 4),
    cv_pct   = round(100 * abs(s / m), 2),
    stringsAsFactors = FALSE
  )
})
summ_df <- do.call(rbind, summ_rows)
out_summ <- file.path(data_dir, "topo_variables_summary.tsv")
write.table(summ_df, out_summ, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("  Wrote: %s\n", out_summ))

# -------------------------------------------------------------------
# STEP 5: Per-site plots
# -------------------------------------------------------------------
cat("STEP 5: per-site plots ...\n")

for (v in topo_vars) {
  plot_df <- topo_df[, c("site", v)]
  colnames(plot_df) <- c("site", "value")
  plot_df <- plot_df[!is.na(plot_df$value), ]
  plot_df <- plot_df[order(plot_df$value, decreasing = TRUE), ]
  plot_df$site <- factor(plot_df$site, levels = plot_df$site)

  p <- ggplot(plot_df, aes(x = site, y = value, fill = site)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = site_colors, na.value = "grey60", guide = "none") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(title = sprintf("Topography - %s per site", v),
         subtitle = "Copernicus DEM GLO-30 (~30 m)", x = NULL, y = v) +
    theme(axis.text.x = element_text(size = 7, hjust = 1),
          plot.subtitle = element_text(size = 8, color = "grey40"))

  ggsave(file.path(plot_dir, sprintf("topo_%s_per_site.png", v)), p,
         width = 10, height = 5, dpi = 300)
}
cat(sprintf("  %d per-site plots written\n", length(topo_vars)))

cat(sprintf("\nDONE 08.1b.3-topo_extraction.R\n"))
cat(sprintf("  %d sites | elevation + tri\n", nrow(topo_df)))

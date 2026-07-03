#!/usr/bin/env Rscript
# 08.1a-env_variable_correlation.R
# Landscape-level pairwise correlation between all continuous environmental
# variables (CHELSA BIO1-19, ENVIREM 18 vars, elevation Copernicus DEM).
#
# Approach: terra samples N random points within French Guiana (~3°N–6°N,
# 55°W–51°W), extracts raster values, applies CHELSA/ENVIREM scale factors,
# and computes Pearson correlations across the landscape.
# Scale factors are validated against the already-computed site-level CSVs.
#
# Inputs:
#   chelsa_dir            — directory with CHELSA_bio{N}_1981-2010_V.2.1.tif
#   envirem_root          — ENVIREM raster root directory
#   dem_dir               — Copernicus DEM GLO-30 tile directory
#   chelsa_site_csv       — chelsa_env_per_site.csv (physical units, 08.0.1 output)
#   envirem_site_csv      — envirem_env_per_site.csv (physical units, 08.0.2 output)
#   manual_site_csv       — manual_variables_per_site.csv (08.0.3 output)
#   out_dir               — output root (08.1-env_variable_selection/)
#   n_points              — number of random landscape points (default 5000)
#   cor_threshold         — |r| threshold for flagging high correlations (default 0.7)
#   terraclimate_site_csv — terraclimate_env_per_site.csv (08.0.4 output)
#
# Outputs:
#   correlation_analysis/landscape_correlation_matrix.csv
#   correlation_analysis/landscape_high_correlations.tsv
#   correlation_analysis/site_level_correlation_matrix.csv   [all 4 sources, n=19 sites]
#   correlation_analysis/site_level_high_correlations.tsv    [all 4 sources, n=19 sites]
#   correlation_analysis/variable_summary.tsv
#   correlation_analysis/variables_to_keep_template.txt
#   plots/landscape_correlation_heatmap.png
#   plots/landscape_correlation_heatmap.pdf
#   plots/site_level_correlation_heatmap.png                 [all 4 sources]
#   plots/site_level_correlation_heatmap.pdf                 [all 4 sources]

suppressPackageStartupMessages({
  if (!requireNamespace("terra", quietly = TRUE))
    stop("Package 'terra' is required. Install with: install.packages('terra')")
  library(terra)
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args           <- commandArgs(trailingOnly = TRUE)
chelsa_dir     <- args[1]
envirem_root   <- args[2]
dem_dir        <- args[3]
chelsa_site_csv  <- args[4]
envirem_site_csv <- args[5]
manual_site_csv  <- args[6]
out_dir        <- args[7]
n_points       <- if (length(args) >= 8) as.integer(args[8]) else 5000L
cor_threshold         <- if (length(args) >= 9) as.numeric(args[9]) else 0.7
terraclimate_site_csv <- if (length(args) >= 10) args[10] else NULL

cor_dir  <- file.path(out_dir, "correlation_analysis")
plot_dir <- file.path(out_dir, "plots")
dir.create(cor_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("INFO: n_points = %d, cor_threshold = %.2f\n", n_points, cor_threshold))

# -------------------------------------------------------------------
# STEP 1: Define raster paths
# -------------------------------------------------------------------
set1     <- file.path(envirem_root, "SAmerica_current_30arcsec_geotiff_set1")
set2     <- file.path(envirem_root, "SAmerica_current_30arcsec_geotiff_set2")
set3     <- file.path(envirem_root, "SAmerica_current_30arcsec_geotiff_set3")
set4     <- file.path(envirem_root, "SAmerica_current_30arcsec_geotiff_set4")
elev_set <- file.path(envirem_root, "elev_SAmerica_current_30arcsec_geotiff")

chelsa_paths <- setNames(
  file.path(chelsa_dir, paste0("CHELSA_bio", 1:19, "_1981-2010_V.2.1.tif")),
  paste0("BIO", 1:19)
)

envirem_paths <- c(
  annualPET                = file.path(set1, "current_30arcsec_annualPET.tif"),
  aridityIndexThornthwaite = file.path(set1, "current_30arcsec_aridityIndexThornthwaite.tif"),
  climaticMoistureIndex    = file.path(set1, "current_30arcsec_climaticMoistureIndex.tif"),
  continentality           = file.path(set1, "current_30arcsec_continentality.tif"),
  embergerQ                = file.path(set2, "current_30arcsec_embergerQ.tif"),
  growingDegDays0          = file.path(set2, "current_30arcsec_growingDegDays0.tif"),
  growingDegDays5          = file.path(set2, "current_30arcsec_growingDegDays5.tif"),
  maxTempColdest           = file.path(set2, "current_30arcsec_maxTempColdest.tif"),
  minTempWarmest           = file.path(set3, "current_30arcsec_minTempWarmest.tif"),
  monthCountByTemp10       = file.path(set3, "current_30arcsec_monthCountByTemp10.tif"),
  PETColdestQuarter        = file.path(set3, "current_30arcsec_PETColdestQuarter.tif"),
  PETDriestQuarter         = file.path(set3, "current_30arcsec_PETDriestQuarter.tif"),
  PETseasonality           = file.path(set4, "current_30arcsec_PETseasonality.tif"),
  PETWarmestQuarter        = file.path(set4, "current_30arcsec_PETWarmestQuarter.tif"),
  PETWettestQuarter        = file.path(set4, "current_30arcsec_PETWettestQuarter.tif"),
  thermicityIndex          = file.path(set4, "current_30arcsec_thermicityIndex.tif"),
  topoWet                  = file.path(elev_set, "current_30arcsec_topoWet.tif"),
  tri                      = file.path(elev_set, "current_30arcsec_tri.tif")
)

# Verify all rasters exist
missing <- c(chelsa_paths, envirem_paths)[!file.exists(c(chelsa_paths, envirem_paths))]
if (length(missing) > 0)
  stop("Missing rasters:\n", paste(" ", names(missing), collapse = "\n"))

# -------------------------------------------------------------------
# STEP 2: Sample random landscape points within French Guiana
# -------------------------------------------------------------------
# Extent covers mainland French Guiana (land only, avoids Atlantic coast noise)
guyane_ext <- terra::ext(-55, -51, 3, 6)

cat("STEP 2: sampling landscape points\n")
# Use BIO1 as reference raster for sampling (define valid land pixels)
ref_rast <- terra::rast(chelsa_paths["BIO1"])
ref_crop  <- terra::crop(ref_rast, guyane_ext)

set.seed(42)
pts    <- terra::spatSample(ref_crop, size = n_points, method = "random",
                             na.rm = TRUE, xy = TRUE)
pts_sv <- terra::vect(pts[, c("x", "y")], geom = c("x", "y"), crs = "EPSG:4326")
cat(sprintf("INFO: %d valid landscape points sampled\n", nrow(pts)))

# -------------------------------------------------------------------
# STEP 3: Extract CHELSA values
# -------------------------------------------------------------------
cat("STEP 3: extracting CHELSA values\n")
chelsa_vals <- data.frame(matrix(NA_real_, nrow = nrow(pts), ncol = 19))
colnames(chelsa_vals) <- paste0("BIO", 1:19)

for (v in paste0("BIO", 1:19)) {
  r   <- terra::rast(chelsa_paths[v])
  ex  <- terra::extract(r, pts_sv, ID = FALSE)
  chelsa_vals[[v]] <- ex[[1]]
}
cat("  CHELSA extraction done\n")

# -------------------------------------------------------------------
# STEP 4: Apply CHELSA v2.1 scale factors
# -------------------------------------------------------------------
# terra applies the internal Int16 Scale=0.1 automatically:
# → absolute temp vars returned in K; precip in mm; differences in °C.
# BIO3 (Float32, Scale=0.1 in metadata): terra returns value×0.1 → multiply ×10.
# Validation against site-level CSV is performed below.
cat("STEP 4: applying CHELSA scale factors\n")

temp_abs_vars  <- c("BIO1","BIO5","BIO6","BIO8","BIO9","BIO10","BIO11")
temp_diff_vars <- c("BIO2","BIO7")

# Detect whether terra applied the Int16 scale (returned in K) or not (raw integers)
bio1_median <- median(chelsa_vals[["BIO1"]], na.rm = TRUE)
cat(sprintf("  BIO1 median raw value from terra: %.2f\n", bio1_median))

if (bio1_median < 50) {
  # terra returned values already in physical units (°C for temps, mm for precip)
  cat("  Detected: terra returned physical units directly — no scaling needed\n")
} else if (bio1_median > 200 & bio1_median < 500) {
  # terra applied Scale=0.1 but not offset → values in K (~298)
  cat("  Detected: values in K → subtracting 273.15\n")
  for (v in temp_abs_vars) chelsa_vals[[v]] <- chelsa_vals[[v]] - 273.15
} else {
  # terra returned raw Int16 integers (~2986)
  cat("  Detected: raw Int16 integers → applying /10 corrections\n")
  for (v in temp_abs_vars)  chelsa_vals[[v]] <- chelsa_vals[[v]] / 10 - 273.15
  for (v in temp_diff_vars) chelsa_vals[[v]] <- chelsa_vals[[v]] / 10
  chelsa_vals[["BIO4"]] <- chelsa_vals[["BIO4"]] / 10
  for (v in paste0("BIO", 12:19)) chelsa_vals[[v]] <- chelsa_vals[[v]] / 10
}

# BIO3 (isothermality %): only multiply ×10 if terra returned val×0.1 (i.e. not already in %)
if (bio1_median < 50) {
  # terra returned physical units directly — BIO3 already in %
} else {
  chelsa_vals[["BIO3"]] <- chelsa_vals[["BIO3"]] * 10
}

# -------------------------------------------------------------------
# STEP 5: Validate CHELSA scale factors against site-level CSV
# -------------------------------------------------------------------
cat("STEP 5: validating CHELSA scale factors\n")
site_chelsa <- read.csv(chelsa_site_csv, stringsAsFactors = FALSE)
# Extract BIO1 at site coordinates using terra
site_sv <- terra::vect(site_chelsa[, c("long", "lat")],
                       geom = c("long", "lat"), crs = "EPSG:4326")
r_bio1  <- terra::rast(chelsa_paths["BIO1"])
bio1_terra_sites <- terra::extract(r_bio1, site_sv, ID = FALSE)[[1]]

# Apply same scaling as above
if (bio1_median > 200) {
  bio1_terra_sites <- bio1_terra_sites - 273.15
} else {
  bio1_terra_sites <- bio1_terra_sites / 10 - 273.15
}

bio1_ref  <- site_chelsa[["BIO1"]]
max_delta <- max(abs(bio1_terra_sites - bio1_ref), na.rm = TRUE)
cat(sprintf("  BIO1 max abs. diff (terra vs site CSV): %.4f°C\n", max_delta))
if (max_delta > 1) {
  warning(sprintf(
    "BIO1 scale validation: max delta = %.3f°C (>1°C threshold).\n",
    max_delta), "Check CHELSA scale factors in the script.")
} else {
  cat("  Scale validation passed (max delta <= 1°C)\n")
}

# -------------------------------------------------------------------
# STEP 6: Extract ENVIREM values
# -------------------------------------------------------------------
cat("STEP 6: extracting ENVIREM values\n")
envirem_vals <- data.frame(matrix(NA_real_, nrow = nrow(pts),
                                  ncol = length(envirem_paths)))
colnames(envirem_vals) <- names(envirem_paths)

for (v in names(envirem_paths)) {
  r  <- terra::rast(envirem_paths[v])
  ex <- terra::extract(r, pts_sv, ID = FALSE)
  envirem_vals[[v]] <- ex[[1]]
}

# ENVIREM scale factors: maxTempColdest and minTempWarmest are Int16 ÷10 → °C
# terra applies Scale=0.1 automatically if present; same detection logic:
maxT_median <- median(envirem_vals[["maxTempColdest"]], na.rm = TRUE)
if (maxT_median > 100) {
  # terra applied scale → already in °C
  cat("  ENVIREM maxTempColdest: terra applied scale (already °C)\n")
} else {
  # raw Int16 → apply ÷10
  cat("  ENVIREM maxTempColdest: applying /10\n")
  envirem_vals[["maxTempColdest"]] <- envirem_vals[["maxTempColdest"]] / 10
  envirem_vals[["minTempWarmest"]] <- envirem_vals[["minTempWarmest"]] / 10
}
cat("  ENVIREM extraction done\n")

# -------------------------------------------------------------------
# STEP 7: Extract elevation (Copernicus DEM, mosaic of tiles)
# -------------------------------------------------------------------
cat("STEP 7: extracting elevation\n")
dem_files <- list.files(dem_dir, pattern = "\\.tif$", full.names = TRUE)
elev_vals <- rep(NA_real_, nrow(pts))

if (length(dem_files) > 0) {
  dem_mosaic <- tryCatch({
    if (length(dem_files) == 1) {
      terra::rast(dem_files)
    } else {
      terra::mosaic(terra::sprc(lapply(dem_files, terra::rast)))
    }
  }, error = function(e) {
    cat(sprintf("  WARN: could not mosaic DEM tiles: %s\n", conditionMessage(e)))
    NULL
  })
  if (!is.null(dem_mosaic)) {
    ex_elev  <- terra::extract(dem_mosaic, pts_sv, ID = FALSE)
    elev_vals <- ex_elev[[1]]
    cat(sprintf("  Elevation range: %.0f – %.0f m\n",
                min(elev_vals, na.rm = TRUE), max(elev_vals, na.rm = TRUE)))
  }
} else {
  cat("  WARN: no DEM tiles found — elevation will be NA\n")
}

# -------------------------------------------------------------------
# STEP 8: Assemble landscape data frame
# -------------------------------------------------------------------
land_df <- cbind(
  pts[, c("x", "y")],
  chelsa_vals,
  envirem_vals,
  elevation = elev_vals
)
colnames(land_df)[1:2] <- c("lon", "lat")

# Remove rows with all-NA environmental values
env_cols_all <- c(paste0("BIO", 1:19), names(envirem_paths), "elevation")
n_before <- nrow(land_df)
land_df   <- land_df[rowSums(!is.na(land_df[, env_cols_all])) > 0, ]
cat(sprintf("INFO: %d/%d points retained after NA filtering\n", nrow(land_df), n_before))

# Drop variables with zero variance (constant across all landscape points)
variances <- sapply(land_df[, env_cols_all], var, na.rm = TRUE)
zero_var  <- names(variances)[variances == 0 | is.na(variances)]
if (length(zero_var) > 0) {
  cat(sprintf("  Dropping zero-variance variables: %s\n", paste(zero_var, collapse = ", ")))
  env_cols_all <- setdiff(env_cols_all, zero_var)
}

# -------------------------------------------------------------------
# STEP 9: Correlation matrix
# -------------------------------------------------------------------
cat("STEP 9: computing correlation matrix\n")
cor_mat  <- cor(land_df[, env_cols_all], use = "pairwise.complete.obs")

out_cor_csv <- file.path(cor_dir, "landscape_correlation_matrix.csv")
write.csv(round(cor_mat, 4), file = out_cor_csv, quote = FALSE)
cat("Correlation matrix written:", out_cor_csv, "\n")

# -------------------------------------------------------------------
# STEP 10: Identify high-correlation pairs
# -------------------------------------------------------------------
high_idx <- which(abs(cor_mat) >= cor_threshold & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(high_idx) > 0) {
  hc_df <- data.frame(
    Var1   = rownames(cor_mat)[high_idx[, 1]],
    Var2   = colnames(cor_mat)[high_idx[, 2]],
    r      = round(cor_mat[high_idx], 4),
    source1 = ifelse(grepl("^BIO", rownames(cor_mat)[high_idx[, 1]]), "CHELSA", "ENVIREM/elev"),
    source2 = ifelse(grepl("^BIO", colnames(cor_mat)[high_idx[, 2]]), "CHELSA", "ENVIREM/elev"),
    stringsAsFactors = FALSE
  )
  hc_df <- hc_df[order(-abs(hc_df$r)), ]
  out_hc <- file.path(cor_dir, "landscape_high_correlations.tsv")
  write.table(hc_df, file = out_hc, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("INFO: %d pairs with |r| >= %.2f written to %s\n",
              nrow(hc_df), cor_threshold, out_hc))
  print(hc_df)
} else {
  cat(sprintf("INFO: no pairs with |r| >= %.2f detected\n", cor_threshold))
}

# -------------------------------------------------------------------
# STEP 11: Variable summary table
# -------------------------------------------------------------------
source_map <- c(
  setNames(rep("CHELSA", 19), paste0("BIO", 1:19)),
  setNames(rep("ENVIREM", 18), names(envirem_paths)),
  elevation = "manual (Copernicus DEM GLO-30)"
)
units_map <- c(
  BIO1="°C", BIO2="°C", BIO3="%", BIO4="SD×100", BIO5="°C", BIO6="°C",
  BIO7="°C", BIO8="°C", BIO9="°C", BIO10="°C", BIO11="°C",
  BIO12="mm", BIO13="mm", BIO14="mm", BIO15="CV", BIO16="mm",
  BIO17="mm", BIO18="mm", BIO19="mm",
  annualPET="mm/yr", aridityIndexThornthwaite="index",
  climaticMoistureIndex="index", continentality="index",
  embergerQ="index", growingDegDays0="°C·days", growingDegDays5="°C·days",
  maxTempColdest="°C", minTempWarmest="°C", monthCountByTemp10="months",
  PETColdestQuarter="mm", PETDriestQuarter="mm", PETseasonality="SD×100",
  PETWarmestQuarter="mm", PETWettestQuarter="mm", thermicityIndex="index",
  topoWet="index", tri="index", elevation="m"
)

summ_df <- data.frame(
  variable = env_cols_all,
  source   = source_map[env_cols_all],
  units    = units_map[env_cols_all],
  n_valid  = colSums(!is.na(land_df[, env_cols_all])),
  mean     = round(colMeans(land_df[, env_cols_all], na.rm = TRUE), 3),
  sd       = round(apply(land_df[, env_cols_all], 2, sd, na.rm = TRUE), 3),
  stringsAsFactors = FALSE,
  row.names = NULL
)
out_summ <- file.path(cor_dir, "variable_summary.tsv")
write.table(summ_df, file = out_summ, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Variable summary written:", out_summ, "\n")
print(summ_df)

# -------------------------------------------------------------------
# STEP 12: Correlation heatmap
# -------------------------------------------------------------------
cat("STEP 12: generating correlation heatmap\n")
cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Var1", "Var2", "r")
cor_long$Var1 <- factor(cor_long$Var1, levels = env_cols_all)
cor_long$Var2 <- factor(cor_long$Var2, levels = env_cols_all)

# Source annotation for axis labels (CHELSA vs ENVIREM)
label_colors <- ifelse(env_cols_all %in% paste0("BIO", 1:19), "#1a6696",
                ifelse(env_cols_all == "elevation", "#8B4513", "#2e8b57"))

p_cor <- ggplot(cor_long, aes(x = Var2, y = Var1, fill = r)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 1.4, color = "black") +
  scale_fill_gradientn(
    colours = c("#d73027", "#f7f7f7", "#4575b4"),
    limits  = c(-1, 1), name = "Pearson r"
  ) +
  labs(title = "Landscape-level Pearson correlation — CHELSA + ENVIREM + elevation",
       subtitle = sprintf("n = %d random points, French Guiana (3-6°N, 55-51°W)",
                          nrow(land_df)),
       x = NULL, y = NULL) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 6,
                                     colour = label_colors),
    axis.text.y       = element_text(size = 6, colour = label_colors),
    plot.title        = element_text(face = "bold", size = 9),
    plot.subtitle     = element_text(size = 7),
    legend.key.height = unit(0.8, "cm")
  )

ggsave(file.path(plot_dir, "landscape_correlation_heatmap.png"),
       p_cor, width = 16, height = 15, dpi = 300)
ggsave(file.path(plot_dir, "landscape_correlation_heatmap.pdf"),
       p_cor, width = 16, height = 15)
cat("Heatmap written\n")

# -------------------------------------------------------------------
# STEP 12b: Site-level joint correlation — all 4 sources + TerraClimate
#
# The landscape raster analysis (CHELSA + ENVIREM, n=5000 pts) captures
# spatial redundancy across French Guiana but excludes TerraClimate.
# The site-level analysis (n=19 sites) covers all 4 sources and is the
# operative correlation reference for GEA variable selection, since
# the GEA models operate on per-site environmental values.
# -------------------------------------------------------------------
cat("STEP 12b: site-level joint correlation (all 4 sources + TerraClimate)\n")

site_hc_df    <- data.frame(Var1 = character(), Var2 = character(),
                             r = numeric(), stringsAsFactors = FALSE)
terra_vars     <- character(0)
site_env_cols  <- NULL

if (!is.null(terraclimate_site_csv) && file.exists(terraclimate_site_csv)) {

  site_chelsa  <- read.csv(chelsa_site_csv,         stringsAsFactors = FALSE)
  site_envirem <- read.csv(envirem_site_csv,         stringsAsFactors = FALSE)
  site_manual  <- read.csv(manual_site_csv,          stringsAsFactors = FALSE)
  site_terra   <- read.csv(terraclimate_site_csv,    stringsAsFactors = FALSE)

  terra_vars <- c("aet_mean", "aet_sd", "def_mean", "def_sd", "ppt_sd",
                  "PDSI_mean", "PDSI_sd", "soil_mean", "soil_sd",
                  "srad_mean", "vpd_mean")
  terra_vars <- terra_vars[terra_vars %in% names(site_terra)]

  site_all <- Reduce(
    function(a, b) merge(a, b, by = "site"),
    list(
      site_chelsa [, c("site", paste0("BIO", 1:19))],
      site_envirem[, c("site", names(envirem_paths)[names(envirem_paths) %in% names(site_envirem)])],
      site_manual [, c("site", "elevation")],
      site_terra  [, c("site", terra_vars)]
    )
  )

  site_env_cols <- c(
    paste0("BIO", 1:19),
    names(envirem_paths)[names(envirem_paths) %in% names(site_all)],
    "elevation",
    terra_vars
  )
  site_env_cols <- site_env_cols[site_env_cols %in% names(site_all)]

  cat(sprintf("  %d sites x %d variables\n", nrow(site_all), length(site_env_cols)))

  # Remove columns with zero variance (constant among sites)
  site_var <- sapply(site_all[, site_env_cols], var, na.rm = TRUE)
  zero_site <- names(site_var)[site_var == 0 | is.na(site_var)]
  if (length(zero_site) > 0) {
    cat(sprintf("  Dropping zero-variance (site-level): %s\n",
                paste(zero_site, collapse = ", ")))
    site_env_cols <- setdiff(site_env_cols, zero_site)
  }

  site_cor_mat <- cor(site_all[, site_env_cols], use = "pairwise.complete.obs")

  out_site_cor <- file.path(cor_dir, "site_level_correlation_matrix.csv")
  write.csv(round(site_cor_mat, 4), file = out_site_cor, quote = FALSE)
  cat("Site-level correlation matrix:", out_site_cor, "\n")

  site_hc_idx <- which(abs(site_cor_mat) >= cor_threshold &
                         upper.tri(site_cor_mat), arr.ind = TRUE)

  var_source <- function(v) {
    if (grepl("^BIO", v)) "CHELSA"
    else if (v %in% terra_vars) "TerraClimate"
    else if (v == "elevation") "manual"
    else "ENVIREM"
  }

  if (nrow(site_hc_idx) > 0) {
    site_hc_df <- data.frame(
      Var1    = rownames(site_cor_mat)[site_hc_idx[, 1]],
      Var2    = colnames(site_cor_mat)[site_hc_idx[, 2]],
      r       = round(site_cor_mat[site_hc_idx], 4),
      stringsAsFactors = FALSE
    )
    site_hc_df$source1 <- sapply(site_hc_df$Var1, var_source)
    site_hc_df$source2 <- sapply(site_hc_df$Var2, var_source)
    site_hc_df <- site_hc_df[order(-abs(site_hc_df$r)), ]
    out_site_hc <- file.path(cor_dir, "site_level_high_correlations.tsv")
    write.table(site_hc_df, file = out_site_hc, sep = "\t",
                quote = FALSE, row.names = FALSE)
    cat(sprintf("INFO: %d site-level pairs |r| >= %.2f -> %s\n",
                nrow(site_hc_df), cor_threshold, out_site_hc))
    print(site_hc_df)
  } else {
    cat(sprintf("INFO: no site-level pairs with |r| >= %.2f\n", cor_threshold))
  }

  # Site-level heatmap — color-coded by source
  site_col_long <- as.data.frame(as.table(site_cor_mat))
  colnames(site_col_long) <- c("Var1", "Var2", "r")
  site_col_long$Var1 <- factor(site_col_long$Var1, levels = site_env_cols)
  site_col_long$Var2 <- factor(site_col_long$Var2, levels = site_env_cols)

  site_label_cols <- sapply(site_env_cols, function(v) {
    if (grepl("^BIO", v)) "#1a6696"
    else if (v %in% terra_vars) "#8B0000"
    else if (v == "elevation") "#8B4513"
    else "#2e8b57"
  })

  p_site_cor <- ggplot(site_col_long, aes(x = Var2, y = Var1, fill = r)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 1.1, color = "black") +
    scale_fill_gradientn(
      colours = c("#d73027", "#f7f7f7", "#4575b4"),
      limits  = c(-1, 1), name = "Pearson r"
    ) +
    labs(
      title    = "Site-level Pearson r — CHELSA + ENVIREM + elevation + TerraClimate",
      subtitle = sprintf(
        "n = %d sites | blue=CHELSA | green=ENVIREM | brown=elevation | red=TerraClimate",
        nrow(site_all)),
      x = NULL, y = NULL
    ) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 5.5,
                                       colour = site_label_cols),
      axis.text.y       = element_text(size = 5.5, colour = site_label_cols),
      plot.title        = element_text(face = "bold", size = 9),
      plot.subtitle     = element_text(size = 7),
      legend.key.height = unit(0.8, "cm")
    )

  ggsave(file.path(plot_dir, "site_level_correlation_heatmap.png"),
         p_site_cor, width = 18, height = 17, dpi = 300)
  ggsave(file.path(plot_dir, "site_level_correlation_heatmap.pdf"),
         p_site_cor, width = 18, height = 17)
  cat("Site-level heatmap written\n")

} else {
  cat("WARN: terraclimate_site_csv not provided or not found — skipping site-level joint correlation\n")
}

# -------------------------------------------------------------------
# STEP 13: Generate variable selection template for user
#
# Correlation counts in the template use the SITE-LEVEL matrix (all 4
# sources, n=19 sites) because GEA models operate on per-site values.
# Existing user selections are preserved if a template already exists.
# -------------------------------------------------------------------
cat("STEP 13: generating variable selection template\n")

tmpl_path <- file.path(cor_dir, "variables_to_keep_template.txt")

# Preserve selections from an existing template (user may have uncommented vars)
old_kept <- character(0)
if (file.exists(tmpl_path)) {
  old_lines <- readLines(tmpl_path)
  for (line in old_lines) {
    lt <- trimws(line)
    if (!startsWith(lt, "#") && nchar(lt) > 0) {
      var_name <- strsplit(lt, "[[:space:]#]")[[1]][1]
      if (nchar(var_name) > 0) old_kept <- c(old_kept, var_name)
    }
  }
  cat(sprintf("INFO: preserving %d selections from existing template\n", length(old_kept)))
}

# Helper: count high-correlation partners from the site-level matrix
n_hc_site <- function(v) {
  if (nrow(site_hc_df) == 0) return(0L)
  sum(site_hc_df$Var1 == v | site_hc_df$Var2 == v)
}

tmpl_lines <- c(
  "# 08.1a -- Variable selection template",
  "# ======================================",
  "# Instructions:",
  "#   1. Site-level heatmap (all 4 sources): plots/site_level_correlation_heatmap.png",
  "#   2. Site-level high-corr pairs:         correlation_analysis/site_level_high_correlations.tsv",
  "#   3. Landscape heatmap (CHELSA+ENVIREM): plots/landscape_correlation_heatmap.png",
  "#   4. For each cluster of correlated variables (|r| >= threshold),",
  "#      keep ONE representative based on ecological relevance.",
  "#   5. Edit this file: uncomment (remove #) the variables you want to KEEP.",
  "#   6. Save as 'variables_to_keep.txt' (same directory).",
  "#   7. This file feeds 08.1b (RDA forward selection + VIF).",
  "#",
  "# Note: elevation, forest habitat dummies, and pedology dummies are treated",
  "#       separately in 08.1b -- do not list them here.",
  "#",
  sprintf("# Generated: %s", Sys.time()),
  sprintf("# Site-level n = 19 sites | Landscape n = %d pts | |r| >= %.2f",
          nrow(land_df), cor_threshold),
  "# Correlation counts refer to site-level Pearson r (all 4 sources combined).",
  "#",
  "# --- CHELSA BIO variables ---"
)

for (v in paste0("BIO", 1:19)) {
  if (!v %in% env_cols_all) next
  n_hc  <- n_hc_site(v)
  flag  <- if (n_hc > 0) sprintf(" # correlated with %d other variable(s)", n_hc) else ""
  pfx   <- if (v %in% old_kept) " " else "# "
  tmpl_lines <- c(tmpl_lines, sprintf("%s%s%s", pfx, v, flag))
}

tmpl_lines <- c(tmpl_lines, "#", "# --- ENVIREM variables ---")
for (v in names(envirem_paths)) {
  if (!v %in% env_cols_all) next
  n_hc  <- n_hc_site(v)
  flag  <- if (n_hc > 0) sprintf(" # correlated with %d other variable(s)", n_hc) else ""
  pfx   <- if (v %in% old_kept) " " else "# "
  tmpl_lines <- c(tmpl_lines, sprintf("%s%s%s", pfx, v, flag))
}

tmpl_lines <- c(tmpl_lines, "#", "# --- Continuous manual variable ---")
{
  v    <- "elevation"
  n_hc <- n_hc_site(v)
  flag <- if (n_hc > 0) sprintf(" # correlated with %d other variable(s)", n_hc) else ""
  pfx  <- if (v %in% old_kept) " " else "# "
  tmpl_lines <- c(tmpl_lines, sprintf("%s%s%s", pfx, v, flag))
}

# TerraClimate section — all commented by default (new source, user must evaluate)
terra_tmpl_vars <- if (length(terra_vars) > 0 && !is.null(site_env_cols))
  terra_vars[terra_vars %in% site_env_cols] else character(0)

if (length(terra_tmpl_vars) > 0) {
  tmpl_lines <- c(
    tmpl_lines,
    "#",
    "# --- TerraClimate variables (1991-2020 30-year normals, ~4 km resolution) ---",
    "# [NEW] Review site_level_correlation_heatmap.png before selecting."
  )
  for (v in terra_tmpl_vars) {
    n_hc <- n_hc_site(v)
    flag <- if (n_hc > 0) sprintf(" # correlated with %d other variable(s)", n_hc) else ""
    pfx  <- if (v %in% old_kept) " " else "# "
    tmpl_lines <- c(tmpl_lines, sprintf("%s%s%s", pfx, v, flag))
  }
}

writeLines(tmpl_lines, tmpl_path)
cat("Variable selection template written:", tmpl_path, "\n")

cat(sprintf(
  "\nSUMMARY\n-------\n  Landscape variables  : %d\n  Landscape points     : %d\n",
  length(env_cols_all), nrow(land_df)))
cat(sprintf("  Landscape high-corr (|r|>=%.2f): %d pairs\n",
            cor_threshold, if (exists("hc_df") && nrow(hc_df) > 0) nrow(hc_df) else 0))
cat(sprintf("  Site-level variables : %d\n",
            if (!is.null(site_env_cols)) length(site_env_cols) else 0))
cat(sprintf("  Site-level high-corr : %d pairs\n", nrow(site_hc_df)))
cat(sprintf("  Zero-variance dropped: %s\n",
            if (length(zero_var) > 0) paste(zero_var, collapse = ", ") else "none"))

cat("\nDONE 08.1a landscape + site-level correlation analysis completed\n")

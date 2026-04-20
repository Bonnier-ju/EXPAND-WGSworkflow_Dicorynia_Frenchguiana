#!/usr/bin/env Rscript
# 08.0.3-manual_variables_extraction.R
# Dummy coding for forest habitat and pedology variables.
#
# The pedology spatial join (point-in-polygon) is performed upstream in the
# SLURM script using Python/OGR, so no spatial libraries are needed here.
#
# Inputs:
#   combined_raw.tsv  — site, lat, long, forest_habitat_code, elevation, pedology_label
#                       (produced by STEP 4 of the .slurm script)
#
# Outputs:
#   data/manual_variables/manual_variables_per_site.csv
#     — one row per site, columns:
#         site, lat, long,
#         forest_habitat_code       (raw nominal code, for reference)
#         elevation                 (continuous, metres, Copernicus DEM GLO-30)
#         pedology_label            (legende_si label, for reference)
#         hab_<CODE>                (0/1 dummy, one per habitat code with >= MIN_SITES)
#         pedol_<LABEL>             (0/1 dummy, one per pedology type with >= MIN_SITES)
#   data/manual_variables/manual_variables_category_summary.tsv
#     — distribution of categories across sites + dummy retention decision
#   plots/manual_variables_habitat_barplot.png
#   plots/manual_variables_pedology_barplot.png
#
# Note: dummy variables with fewer than MIN_SITES sites coded as 1 are
# excluded — insufficient variation for GEA. These variables are listed
# in the summary table with retained = FALSE.
# Categorical variables are NOT included in the landscape correlation step
# (08.1a) but ARE available for RDA forward selection (08.1b) as factors.

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# Site color palette
SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else { c() }

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
input_tsv <- if (length(args) >= 1) args[1] else stop("combined_raw.tsv required")
out_dir   <- if (length(args) >= 2) args[2] else stop("out_dir required")
min_sites <- if (length(args) >= 3) as.integer(args[3]) else 3L

man_out_dir <- file.path(out_dir, "data/manual_variables")
plot_dir    <- file.path(out_dir, "plots")
dir.create(man_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,    recursive = TRUE, showWarnings = FALSE)

cat(sprintf("INFO: dummy coding threshold: >= %d sites per category\n", min_sites))

# -------------------------------------------------------------------
# STEP 1: Load combined table
# -------------------------------------------------------------------
cat("Loading combined table:", input_tsv, "\n")
df <- read.table(input_tsv, header = TRUE, sep = "\t",
                 quote = "", stringsAsFactors = FALSE, check.names = FALSE)
df$forest_habitat_code <- as.character(as.integer(df$forest_habitat_code))
df$elevation           <- suppressWarnings(as.numeric(df$elevation))
cat(sprintf("INFO: %d sites loaded\n", nrow(df)))
cat(sprintf("INFO: elevation range: %.0f – %.0f m\n",
            min(df$elevation, na.rm = TRUE), max(df$elevation, na.rm = TRUE)))

cat("INFO: unique pedology types at sites:\n")
print(table(df$pedology_label))

# -------------------------------------------------------------------
# STEP 2: Dummy coding helper
# -------------------------------------------------------------------
make_dummies <- function(df, col, prefix, min_n) {
  vals     <- df[[col]]
  tab      <- sort(table(vals), decreasing = TRUE)
  retained <- names(tab)[tab >= min_n]
  dropped  <- names(tab)[tab < min_n]

  if (length(dropped) > 0)
    cat(sprintf("  Dropped (%s, < %d sites): %s\n",
                prefix, min_n, paste(dropped, collapse = ", ")))
  if (length(retained) == 0) {
    cat(sprintf("  WARNING: no %s category passes the threshold — no dummies created\n", prefix))
    return(list(df = df, summary = data.frame(
      variable   = names(tab),
      n_sites    = as.integer(tab),
      retained   = names(tab) %in% retained,
      dummy_col  = NA_character_
    )))
  }

  dummy_cols <- character(0)
  for (v in retained) {
    # Build clean column name: replace non-alphanumeric with underscore, truncate
    clean <- gsub("[^A-Za-z0-9]", "_", v)
    clean <- gsub("_+", "_", clean)
    clean <- sub("^_|_$", "", clean)
    clean <- substr(clean, 1, 40)
    cname <- paste0(prefix, "_", clean)
    df[[cname]] <- as.integer(vals == v)
    dummy_cols  <- c(dummy_cols, cname)
    cat(sprintf("  Dummy: %-50s n=%d\n", cname, sum(df[[cname]])))
  }

  summ <- data.frame(
    variable  = names(tab),
    n_sites   = as.integer(tab),
    retained  = names(tab) %in% retained,
    dummy_col = ifelse(names(tab) %in% retained,
                       paste0(prefix, "_", substr(
                         gsub("_+", "_", sub("^_|_$", "",
                           gsub("[^A-Za-z0-9]", "_", names(tab)))), 1, 40)),
                       NA_character_),
    stringsAsFactors = FALSE
  )
  list(df = df, summary = summ, dummy_cols = dummy_cols)
}

# -------------------------------------------------------------------
# STEP 3: Dummy code forest habitat
# -------------------------------------------------------------------
cat("\nSTEP 3: dummy coding forest habitat (prefix: hab)\n")
res_hab  <- make_dummies(df, "forest_habitat_code", "hab", min_sites)
df       <- res_hab$df
summ_hab <- cbind(source = "forest_habitat", res_hab$summary)

# -------------------------------------------------------------------
# STEP 4: Dummy code pedology
# -------------------------------------------------------------------
cat("\nSTEP 4: dummy coding pedology (prefix: pedol)\n")
res_pedol  <- make_dummies(df, "pedology_label", "pedol", min_sites)
df         <- res_pedol$df
summ_pedol <- cbind(source = "pedology", res_pedol$summary)

# -------------------------------------------------------------------
# STEP 5: Write outputs
# -------------------------------------------------------------------
out_csv <- file.path(man_out_dir, "manual_variables_per_site.csv")
write.csv(df, file = out_csv, row.names = FALSE, quote = FALSE)
cat("\nPer-site table written:", out_csv, "\n")
print(df[, c("site", "forest_habitat_code", "elevation", "pedology_label")])

summ_all <- rbind(summ_hab, summ_pedol)
out_summ <- file.path(man_out_dir, "manual_variables_category_summary.tsv")
write.table(summ_all, file = out_summ, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Category summary written:", out_summ, "\n")
cat("\n--- Category summary ---\n")
print(summ_all)

# -------------------------------------------------------------------
# STEP 6: Barplots — category distribution across sites
# -------------------------------------------------------------------
# Forest habitat
hab_tab <- as.data.frame(table(site = df$site,
                               code = df$forest_habitat_code))
hab_tab <- hab_tab[hab_tab$Freq > 0, ]
p_hab <- ggplot(hab_tab, aes(x = site, fill = code)) +
  geom_bar(stat = "identity", aes(y = Freq), width = 0.7) +
  labs(title = "Forest habitat type per site",
       x = NULL, y = "Count", fill = "Habitat code") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title  = element_text(face = "bold"))
ggsave(file.path(plot_dir, "manual_variables_habitat_barplot.png"),
       p_hab, width = 10, height = 4, dpi = 300)

# Pedology
pedol_tab <- as.data.frame(table(site  = df$site,
                                 label = df$pedology_label))
pedol_tab <- pedol_tab[pedol_tab$Freq > 0, ]
pedol_tab$label_short <- substr(pedol_tab$label, 1, 35)
p_pedol <- ggplot(pedol_tab, aes(x = site, fill = label_short)) +
  geom_bar(stat = "identity", aes(y = Freq), width = 0.7) +
  labs(title = "Pedology type per site",
       x = NULL, y = "Count", fill = "Pedology (legende_si)") +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
        plot.title   = element_text(face = "bold"),
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 8))
ggsave(file.path(plot_dir, "manual_variables_pedology_barplot.png"),
       p_pedol, width = 12, height = 5, dpi = 300)

cat("Barplots written\n")

# -------------------------------------------------------------------
# Summary of retained dummy variables
# -------------------------------------------------------------------
all_dummy_cols <- c(res_hab$dummy_cols, res_pedol$dummy_cols)
cat(sprintf("\n--- %d dummy variable(s) retained for GEA (threshold: >= %d sites) ---\n",
            length(all_dummy_cols), min_sites))
for (col in all_dummy_cols) {
  cat(sprintf("  %-55s  n_sites_1 = %d\n", col, sum(df[[col]] == 1, na.rm = TRUE)))
}

cat("\nNOTE: these dummy variables are excluded from landscape-level correlation (08.1a)")
cat("\n      and from LFMM/BayPass/GEMMA. They are available for RDA (08.1b) as factors.\n")

cat("\nDONE 08.0.3 manual variables extraction completed\n")
cat("Output:", out_csv, "\n")

#!/usr/bin/env Rscript
# 07.0-prep_cluster_assignments.R
# Builds cluster_assignments.tsv from sites_by_clusters.csv + complete_metadata.csv.
# Also creates the subdirectory skeleton for PCA, Admixture and Diversity analyses.
#
# Run: Rscript --vanilla 07.0-prep_cluster_assignments.R

BASE            <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
SITES_CLUSTERS  <- file.path(BASE, "metadata/sites_by_clusters.csv")
METADATA        <- file.path(BASE, "metadata/complete_metadata.csv")
OUT_FILE        <- file.path(BASE, "metadata/cluster_assignments.tsv")

EXCLUDE_SITES   <- c("Cameroun_Benin", "Herbier")
EXCLUDE_IDS     <- c("DBR_AGJ", "DBR_AGK", "DBR_AGL", "DBR_AGM")
K_VALUES        <- c(2L, 3L, 4L)

ANALYSIS_DIRS <- list(
  list(dir = file.path(BASE, "07-population_structure/workflow_full/07.1-pca_ibd_fst"),
       prefix = "07.1"),
  list(dir = file.path(BASE, "07-population_structure/workflow_full/07.2-admixture"),
       prefix = "07.2"),
  list(dir = file.path(BASE, "07-population_structure/workflow_full/07.3-diversity_metrics"),
       prefix = "07.3")
)

# -------------------------------------------------------------------
# 1. Read inputs (all files are now pure ASCII — no encoding issues)
# -------------------------------------------------------------------
sites_clust <- read.csv(SITES_CLUSTERS, stringsAsFactors = FALSE)
colnames(sites_clust) <- c("site", "K2", "K3", "K4")
sites_clust[] <- lapply(sites_clust, trimws)

meta <- read.csv(METADATA, stringsAsFactors = FALSE)
meta <- meta[, c("id_genoscope", "sites")]
colnames(meta) <- c("id_genoscope", "site")
meta[] <- lapply(meta, trimws)

# -------------------------------------------------------------------
# 2. Filter: remove outgroups and treemutation individuals
# -------------------------------------------------------------------
meta_filt <- meta[!meta$site %in% EXCLUDE_SITES & !meta$id_genoscope %in% EXCLUDE_IDS, ]
cat(sprintf("Individuals retained (guiana_only): %d\n", nrow(meta_filt)))

# -------------------------------------------------------------------
# 3. Join on site name
# -------------------------------------------------------------------
assignments <- merge(meta_filt, sites_clust, by = "site", all.x = TRUE)
assignments <- assignments[order(assignments$id_genoscope),
                           c("id_genoscope", "site", "K2", "K3", "K4")]

unmatched <- assignments[is.na(assignments$K2), ]
if (nrow(unmatched) > 0) {
  cat("WARNING — unmatched sites (NA in assignments):\n")
  print(unique(unmatched$site))
} else {
  cat("All sites matched successfully.\n")
}

# -------------------------------------------------------------------
# 4. Write cluster_assignments.tsv
# -------------------------------------------------------------------
write.table(assignments, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Written: %s (%d individuals)\n", OUT_FILE, nrow(assignments)))

for (k in K_VALUES) {
  col  <- paste0("K", k)
  pops <- sort(unique(assignments[[col]]))
  cat(sprintf("K=%d — %d populations: %s\n", k, length(pops), paste(pops, collapse = ", ")))
}

# -------------------------------------------------------------------
# 5. Create subdirectory skeleton
# -------------------------------------------------------------------
cat("\nCreating directory structure...\n")

for (item in ANALYSIS_DIRS) {
  base_dir <- item$dir
  prefix   <- item$prefix

  # K=1 subdir — placeholder for existing global analysis
  dir.create(file.path(base_dir, sprintf("%s.1-K=1", prefix)),
             recursive = TRUE, showWarnings = FALSE)

  for (k in K_VALUES) {
    sub_dir <- file.path(base_dir, sprintf("%s.%d-K=%d", prefix, k, k))
    pops    <- sort(unique(assignments[[paste0("K", k)]]))
    for (pop in pops) {
      for (subf in c("plots", "logs", "tmp")) {
        dir.create(file.path(sub_dir, pop, subf), recursive = TRUE, showWarnings = FALSE)
      }
    }
    cat(sprintf("  %s  [%s]\n", sub_dir, paste(pops, collapse = ", ")))
  }
}

cat("\n--- NOTE ---\n")
cat("Existing global outputs are NOT moved automatically.\n")
cat("Copy them manually into the K=1 placeholders if needed:\n")
cat(sprintf("  cp %s/07.1-pca_ibd_fst/pca_guiana_only.* ",  BASE))
cat(sprintf("%s/07.1-pca_ibd_fst/07.1.1-K=1/\n", BASE))
cat(sprintf("  cp -r %s/07.2-admixture/guiana_only/ ",       BASE))
cat(sprintf("%s/07.2-admixture/07.2.1-K=1/\n", BASE))
cat("---\n")
cat("DONE.\n")

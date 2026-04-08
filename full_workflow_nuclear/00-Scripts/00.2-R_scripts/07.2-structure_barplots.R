#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# Site-specific color palette (consistent across all population graphics)
SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else {
  warning("Site color file not found: ", SITE_COLORS_FILE, " — using ggplot2 defaults")
  c()
}

args <- commandArgs(trailingOnly = TRUE)

base_dir_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_full/07.2-structure"
fam_file_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_full/07.1-pca_ibd_fst/genotypes.fam"
pop_map_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_full/07.1-pca_ibd_fst/sample_population_map.tsv"

base_dir <- if (length(args) >= 1) args[1] else base_dir_default
fam_file <- if (length(args) >= 2) args[2] else fam_file_default
pop_map_file <- if (length(args) >= 3) args[3] else pop_map_default

if (!file.exists(base_dir)) stop("Missing base dir: ", base_dir)
if (!file.exists(fam_file)) stop("Missing fam file: ", fam_file)
if (!file.exists(pop_map_file)) stop("Missing population map: ", pop_map_file)

out_dir <- file.path(base_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
sample_ids <- sample_fam$V2
n_samples <- length(sample_ids)

pop_map <- read.table(pop_map_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(pop_map) <- c("sample_id", "population")
pop_by_sample <- setNames(pop_map$population, pop_map$sample_id)

parse_structure_file <- function(fpath, n_expected) {
  lines <- readLines(fpath, warn = FALSE)

  fname <- basename(fpath)
  k <- as.integer(sub(".*_k([0-9]+)_r[0-9]+_f$", "\\1", fname))
  rep <- as.integer(sub(".*_k[0-9]+_r([0-9]+)_f$", "\\1", fname))

  lnp_line <- grep("Estimated Ln Prob of Data", lines, value = TRUE)
  lnp <- NA_real_
  if (length(lnp_line) > 0) {
    v <- sub(".*=\\s*", "", lnp_line[1])
    lnp <- suppressWarnings(as.numeric(v))
  }

  start <- grep("^\\s*Inferred ancestry of individuals", lines)
  if (length(start) == 0) return(NULL)

  qrows <- list()
  for (i in seq.int(start[1] + 1, length(lines))) {
    ln <- lines[i]
    if (grepl("^\\s*$", ln)) next
    if (grepl("^\\s*Estimated Allele Frequencies", ln)) break
    if (!grepl(":", ln, fixed = TRUE)) next

    rhs <- strsplit(ln, ":", fixed = TRUE)[[1]]
    if (length(rhs) < 2) next
    q <- suppressWarnings(as.numeric(strsplit(trimws(rhs[2]), "\\s+")[[1]]))
    q <- q[!is.na(q)]
    if (length(q) > 0) qrows[[length(qrows) + 1]] <- q
  }

  if (length(qrows) != n_expected) return(NULL)

  mat <- do.call(rbind, qrows)
  if (!is.matrix(mat)) return(NULL)

  list(k = k, rep = rep, lnp = lnp, qmat = mat, file = fpath)
}

files <- list.files(base_dir, pattern = "^str_k[0-9]+_r[0-9]+_f$", full.names = TRUE)
if (length(files) == 0) stop("No STRUCTURE output files found in: ", base_dir)

parsed <- lapply(files, parse_structure_file, n_expected = n_samples)
parsed <- parsed[!vapply(parsed, is.null, logical(1))]
if (length(parsed) == 0) stop("No parseable STRUCTURE outputs found.")

meta <- data.frame(
  k = vapply(parsed, function(x) x$k, integer(1)),
  rep = vapply(parsed, function(x) x$rep, integer(1)),
  lnp = vapply(parsed, function(x) x$lnp, numeric(1)),
  file = vapply(parsed, function(x) x$file, character(1)),
  stringsAsFactors = FALSE
)

# choose best replicate per K (highest LnP(D), fallback first if NA)
best_idx <- c()
for (k in sort(unique(meta$k))) {
  idx <- which(meta$k == k)
  if (all(is.na(meta$lnp[idx]))) {
    best_idx <- c(best_idx, idx[1])
  } else {
    best_idx <- c(best_idx, idx[which.max(meta$lnp[idx])])
  }
}
best <- parsed[best_idx]

records <- list()
for (obj in best) {
  q <- obj$qmat
  k <- obj$k
  for (cidx in seq_len(ncol(q))) {
    records[[length(records) + 1]] <- data.frame(
      sample_id = sample_ids,
      K = paste0("K", k),
      cluster = paste0("C", cidx),
      ancestry = q[, cidx],
      stringsAsFactors = FALSE
    )
  }
}
plot_df <- do.call(rbind, records)
plot_df$population <- pop_by_sample[plot_df$sample_id]
plot_df$population[is.na(plot_df$population) | plot_df$population == ""] <- "UNKNOWN"

sample_order <- unique(plot_df[order(plot_df$population, plot_df$sample_id), "sample_id"])
plot_df$sample_id <- factor(plot_df$sample_id, levels = sample_order)
plot_df$K <- factor(plot_df$K, levels = paste0("K", sort(unique(meta$k))))

# Population annotation strip: one colored tile per individual showing site membership.
# Built as a separate thin ggplot and combined with the ancestry barplot via gridExtra.
# Falls back to saving ancestry plot only if gridExtra is unavailable.
make_pop_strip <- function(sample_lvls, pop_by_sample, color_values) {
  strip_df <- data.frame(
    sample_id = factor(sample_lvls, levels = sample_lvls),
    population = pop_by_sample[sample_lvls],
    y = 1,
    stringsAsFactors = FALSE
  )
  strip_df$population[is.na(strip_df$population) | strip_df$population == ""] <- "UNKNOWN"
  p <- ggplot(strip_df, aes(x = sample_id, y = y, fill = population)) +
    geom_col(width = 1) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "none"
    )
  if (length(color_values) > 0) {
    p <- p + scale_fill_manual(values = color_values, na.value = "grey80")
  }
  p
}

save_with_strip <- function(p_ancestry, p_strip, filepath_base, width, height) {
  # Try gridExtra for combining ancestry barplot + population strip
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    g <- gridExtra::arrangeGrob(p_ancestry, p_strip, nrow = 2, heights = c(10, 0.6))
    ggplot2::ggsave(paste0(filepath_base, ".png"), g, width = width, height = height + 0.5, dpi = 300)
    ggplot2::ggsave(paste0(filepath_base, ".pdf"), g, width = width, height = height + 0.5)
  } else {
    # gridExtra not available: save ancestry barplot only
    ggsave(paste0(filepath_base, ".png"), p_ancestry, width = width, height = height, dpi = 300)
    ggsave(paste0(filepath_base, ".pdf"), p_ancestry, width = width, height = height)
    ggsave(paste0(filepath_base, "_pop_strip.png"), p_strip, width = width, height = 0.6, dpi = 300)
  }
}

# Global multi-K plot
p_all <- ggplot(plot_df, aes(x = sample_id, y = ancestry, fill = cluster)) +
  geom_col(width = 1) +
  facet_wrap(~K, ncol = 3, scales = "free_x") +
  labs(
    title = "STRUCTURE Admixture Barplots (best replicate per K)",
    x = "Individuals",
    y = "Ancestry proportion"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "structure_barplot_allK_best_rep.png"), p_all, width = 14, height = 10, dpi = 300)
ggsave(file.path(out_dir, "structure_barplot_allK_best_rep.pdf"), p_all, width = 14, height = 10)

# Per-K plots with population annotation strip
for (k in levels(plot_df$K)) {
  d <- plot_df[plot_df$K == k, , drop = FALSE]
  k_sample_order <- levels(droplevels(d$sample_id))

  p_k <- ggplot(d, aes(x = sample_id, y = ancestry, fill = cluster)) +
    geom_col(width = 1) +
    labs(
      title = paste0("STRUCTURE Barplot - ", k, " (best replicate)"),
      x = NULL,
      y = "Ancestry proportion"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  p_strip_k <- make_pop_strip(k_sample_order, pop_by_sample, site_colors)
  fn <- file.path(out_dir, paste0("structure_barplot_", k, "_best_rep"))
  save_with_strip(p_k, p_strip_k, fn, width = 12, height = 5)
}

write.table(
  meta[order(meta$k, meta$rep), ],
  file = file.path(out_dir, "structure_runs_lnpd.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

best_meta <- meta[best_idx, ]
best_meta <- best_meta[order(best_meta$k), ]
write.table(
  best_meta,
  file = file.path(out_dir, "structure_best_replicate_per_k.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("STRUCTURE barplots written to:", out_dir, "\n")

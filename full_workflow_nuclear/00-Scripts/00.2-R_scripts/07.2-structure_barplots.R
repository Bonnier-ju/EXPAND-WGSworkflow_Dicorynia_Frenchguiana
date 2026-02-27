#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

base_dir_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_test/07.2-structure"
fam_file_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_test/07.1-pca_ibd_fst/genotypes.fam"
pop_map_default <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/07-population_structure/workflow_test/07.1-pca_ibd_fst/sample_population_map.tsv"

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

# Per-K plots
for (k in levels(plot_df$K)) {
  d <- plot_df[plot_df$K == k, , drop = FALSE]
  p_k <- ggplot(d, aes(x = sample_id, y = ancestry, fill = cluster)) +
    geom_col(width = 1) +
    labs(
      title = paste0("STRUCTURE Barplot - ", k, " (best replicate)"),
      x = "Individuals",
      y = "Ancestry proportion"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  fn <- paste0("structure_barplot_", k, "_best_rep")
  ggsave(file.path(out_dir, paste0(fn, ".png")), p_k, width = 12, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0(fn, ".pdf")), p_k, width = 12, height = 5)
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

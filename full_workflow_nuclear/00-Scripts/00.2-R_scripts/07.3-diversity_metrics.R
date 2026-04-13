#!/usr/bin/env Rscript
# 07.3-diversity_metrics.R
# Per-population diversity metrics:
#   - Ho, He, FIS, Ar (rarefied) : hierfstat on LD-pruned VCF (06.6b)
#   - pi, Tajima's D             : vcftools outputs (06.5, all SNPs)
#   - Private alleles            : vcftools --freq2 outputs (06.5)

suppressPackageStartupMessages({
  library(ggplot2)
  library(hierfstat)
  library(gaston)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

# Site-specific color palette
SITE_COLORS_FILE <- "/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sites_couleurs.csv"
site_colors <- if (file.exists(SITE_COLORS_FILE)) {
  pal <- read.csv(SITE_COLORS_FILE, stringsAsFactors = FALSE)
  setNames(pal$couleur_hex, pal$site)
} else {
  warning("Site color file not found: ", SITE_COLORS_FILE)
  c()
}

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
out_dir     <- if (length(args) >= 1) args[1] else stop("out_dir required")
map_file    <- if (length(args) >= 2) args[2] else stop("map_file required")
pops_file   <- if (length(args) >= 3) args[3] else stop("pops_file required")
min_n_ar    <- if (length(args) >= 4) as.integer(args[4]) else 4L
vcf_hf_file <- if (length(args) >= 5) args[5] else stop("vcf_hf_file required (LD-pruned VCF)")

pop_dir  <- file.path(out_dir, "per_population")
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Load population map and population list
# -------------------------------------------------------------------
map  <- read.table(map_file, header = FALSE, sep = "\t",
                   col.names = c("sample_id", "population"),
                   stringsAsFactors = FALSE)
pops <- readLines(pops_file)
pops <- pops[nzchar(pops)]

# -------------------------------------------------------------------
# Helper: safe read of vcftools output
# -------------------------------------------------------------------
read_vcftools <- function(path, ...) {
  if (file.exists(path) && file.info(path)$size > 0) {
    tryCatch(read.table(path, header = TRUE, stringsAsFactors = FALSE, ...),
             error = function(e) NULL)
  } else NULL
}

# -------------------------------------------------------------------
# STEP 1: hierfstat — Ho, He, FIS, Ar from LD-pruned VCF (06.6b)
# -------------------------------------------------------------------
cat(sprintf("INFO: loading LD-pruned VCF: %s\n", vcf_hf_file))
if (!file.exists(vcf_hf_file)) stop("LD-pruned VCF not found: ", vcf_hf_file)

x <- gaston::read.vcf(vcf_hf_file, max.snps = 1e6, verbose = FALSE)
cat(sprintf("INFO: %d individuals x %d SNPs loaded\n", nrow(x), ncol(x)))

# Match VCF samples to population map
pop_vector <- map$population[match(x@ped$id, map$sample_id)]
n_unknown  <- sum(is.na(pop_vector))
if (n_unknown > 0) {
  cat(sprintf("WARN: %d sample(s) not in population map — excluded from hierfstat\n", n_unknown))
}

# Keep only mapped samples
keep_mask  <- !is.na(pop_vector)
x_hf       <- x[keep_mask, ]
pop_vector <- pop_vector[keep_mask]

# Population factor with sorted levels (ensures consistent integer codes)
pop_levels <- sort(unique(pop_vector))
pop_factor <- factor(pop_vector, levels = pop_levels)

# Convert gaston dosage (0/1/2) to hierfstat encoding (11/12/22)
cat("INFO: converting genotype matrix to hierfstat format...\n")
geno_dos <- as.matrix(x_hf)  # individuals x SNPs, values 0/1/2/NA
geno_hf  <- matrix(NA_integer_, nrow = nrow(geno_dos), ncol = ncol(geno_dos))
geno_hf[geno_dos == 0] <- 11L
geno_hf[geno_dos == 1] <- 12L
geno_hf[geno_dos == 2] <- 22L

hf_df <- cbind(data.frame(pop = pop_factor), as.data.frame(geno_hf))
rm(geno_dos, geno_hf, x, x_hf); gc()

# --- basic.stats: Ho, He (Hs = within-population gene diversity), FIS ---
cat("INFO: running hierfstat::basic.stats()...\n")
bs <- hierfstat::basic.stats(hf_df)

# bs$Ho, bs$Hs, bs$Fis are loci x populations matrices
# Column names match the population factor levels
ho_per_pop  <- colMeans(bs$Ho,  na.rm = TRUE)
he_per_pop  <- colMeans(bs$Hs,  na.rm = TRUE)
fis_per_pop <- colMeans(bs$Fis, na.rm = TRUE)

# --- allelic.richness ---
cat(sprintf("INFO: running hierfstat::allelic.richness(min.n=%d)...\n", min_n_ar))
ar_res     <- hierfstat::allelic.richness(hf_df, min.n = min_n_ar)
ar_per_pop <- colMeans(ar_res$Ar, na.rm = TRUE)

# Assemble hierfstat results aligned on pop_levels
hf_metrics <- data.frame(
  population = pop_levels,
  Ho         = as.numeric(ho_per_pop[as.character(pop_levels)]),
  He         = as.numeric(he_per_pop[as.character(pop_levels)]),
  FIS        = as.numeric(fis_per_pop[as.character(pop_levels)]),
  Ar         = as.numeric(ar_per_pop[as.character(pop_levels)]),
  stringsAsFactors = FALSE
)

rm(hf_df, bs, ar_res); gc()
cat("INFO: hierfstat metrics done\n")

# -------------------------------------------------------------------
# STEP 2: pi and Tajima's D from vcftools outputs (06.5, all SNPs)
# -------------------------------------------------------------------
vcftools_metrics <- lapply(pops, function(pop) {
  safe  <- gsub("[^A-Za-z0-9._-]", "_", pop)
  pop_d <- file.path(pop_dir, safe)
  n_ind <- sum(map$population == pop)

  pi_dat   <- read_vcftools(file.path(pop_d, "pi.sites.pi"))
  pi_mean  <- if (!is.null(pi_dat) && "PI" %in% colnames(pi_dat))
    mean(pi_dat$PI, na.rm = TRUE) else NA_real_

  taj_dat  <- read_vcftools(file.path(pop_d, "tajima.Tajima.D"))
  taj_mean <- if (!is.null(taj_dat) && "TajimaD" %in% colnames(taj_dat))
    mean(taj_dat$TajimaD, na.rm = TRUE) else NA_real_

  data.frame(population = pop, n_individuals = n_ind,
             pi_mean = pi_mean, TajimaD_mean = taj_mean,
             stringsAsFactors = FALSE)
})
vcftools_df <- do.call(rbind, vcftools_metrics)

# -------------------------------------------------------------------
# STEP 3: Private alleles from vcftools --freq2 outputs (06.5)
# -------------------------------------------------------------------
cat("INFO: computing private alleles from vcftools freq tables...\n")

freq_tables <- lapply(pops, function(pop) {
  safe   <- gsub("[^A-Za-z0-9._-]", "_", pop)
  freq_f <- file.path(pop_dir, safe, "freq.frq")
  d <- read_vcftools(freq_f, comment.char = "", sep = "\t")
  if (is.null(d) || ncol(d) < 6) return(NULL)
  d$locus_id   <- paste0(d[[1]], ":", d[[2]])
  d$freq_minor <- as.numeric(d[[6]])
  d[, c("locus_id", "freq_minor")]
})
names(freq_tables) <- pops

valid_pops <- pops[!sapply(freq_tables, is.null)]
if (length(valid_pops) >= 2) {
  merged <- Reduce(function(a, b) merge(a, b, by = "locus_id"),
    lapply(valid_pops, function(p) {
      d <- freq_tables[[p]]; colnames(d)[2] <- p; d
    }))

  private_counts <- sapply(valid_pops, function(pop) {
    others    <- setdiff(valid_pops, pop)
    own_freq  <- merged[[pop]]
    other_mat <- as.matrix(merged[, others, drop = FALSE])
    sum(own_freq > 0 & rowSums(other_mat > 0, na.rm = TRUE) == 0, na.rm = TRUE)
  })

  private_df <- data.frame(population      = valid_pops,
                            private_alleles = as.integer(private_counts),
                            stringsAsFactors = FALSE)
} else {
  private_df <- data.frame(population      = pops,
                            private_alleles = NA_integer_,
                            stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------
# STEP 4: Merge all metrics into final summary table
# -------------------------------------------------------------------
final_df <- Reduce(function(a, b) merge(a, b, by = "population", all = TRUE),
                   list(vcftools_df, hf_metrics, private_df))

num_cols <- c("pi_mean", "Ho", "He", "FIS", "TajimaD_mean", "Ar")
final_df[num_cols] <- lapply(final_df[num_cols], function(x) round(x, 6))
final_df <- final_df[order(final_df$population), ]

out_tsv <- file.path(out_dir, "07.3-diversity_metrics_summary.tsv")
write.table(final_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary table written:", out_tsv, "\n")
print(final_df)

# -------------------------------------------------------------------
# STEP 5: Barplots per metric
# -------------------------------------------------------------------
metrics_to_plot <- list(
  list(col = "pi_mean",         label = "Mean nucleotide diversity (\u03c0)"),
  list(col = "Ho",              label = "Observed heterozygosity (Ho)"),
  list(col = "He",              label = "Expected heterozygosity (He)"),
  list(col = "FIS",             label = "Inbreeding coefficient (FIS)"),
  list(col = "TajimaD_mean",    label = "Mean Tajima's D"),
  list(col = "Ar",              label = sprintf("Allelic richness rarefied (n=%d ind.)", min_n_ar)),
  list(col = "private_alleles", label = "Number of private alleles")
)

for (m in metrics_to_plot) {
  col   <- m$col
  label <- m$label
  df_p  <- final_df[!is.na(final_df[[col]]), ]
  if (nrow(df_p) == 0) next
  df_p  <- df_p[order(-df_p[[col]]), ]
  df_p$population <- factor(df_p$population, levels = df_p$population)

  p <- ggplot(df_p, aes_string(x = "population", y = col, fill = "population")) +
    geom_col(width = 0.7) +
    geom_text(aes_string(label = sprintf('sprintf("%%.4g", %s)', col)),
              vjust = -0.4, size = 2.6) +
    labs(title = label, x = NULL, y = label) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      plot.title      = element_text(face = "bold"),
      legend.position = "none"
    )

  if (length(site_colors) > 0)
    p <- p + scale_fill_manual(values = site_colors, na.value = "grey50")

  fn <- file.path(plot_dir, paste0(gsub("[^A-Za-z0-9]", "_", col), ".png"))
  ggsave(fn, p, width = max(10, nrow(df_p) * 0.5), height = 5, dpi = 300)
  cat("Plot written:", fn, "\n")
}

cat("DONE 07.3 diversity metrics R script completed\n")

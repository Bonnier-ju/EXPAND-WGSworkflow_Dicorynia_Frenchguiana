#!/usr/bin/env Rscript
# 07b.3-diversity_metrics.R
# Diversity metrics for purity-filtered individuals, grouped by K=3 cluster.
# Metrics:
#   - Ho, He, FIS, Ar (rarefied, hierfstat)   : LD-pruned VCF (06.6b)
#   - pi, Tajima's D (vcftools)                : all-SNP VCF (06.5)
#   - Private alleles (count + rarefied PAR)   : vcftools --freq2 outputs
#
# Arguments:
#   1  out_dir        threshold output directory (e.g. 07b.3-diversity_metrics/T70)
#   2  map_file       TSV: sample_id <tab> cluster (Pop_1/2/3)
#   3  min_n_ar       rarefaction sample size (individuals; default 10)
#   4  vcf_ld_file    LD-pruned VCF (06.6b) for hierfstat
#   5  cluster_dir    directory containing per-cluster vcftools outputs
#   6  threshold_tag  label for filenames (T70|T80|T90|T95)

suppressPackageStartupMessages({
  library(ggplot2)
  library(hierfstat)
  library(gaston)
})

theme_set(theme_minimal() + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

POP_COLORS <- c(Pop_1 = "#EE7600", Pop_2 = "#458B00", Pop_3 = "#CD2626")

# -------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------
args          <- commandArgs(trailingOnly = TRUE)
out_dir       <- if (length(args) >= 1) args[1] else stop("out_dir required")
map_file      <- if (length(args) >= 2) args[2] else stop("map_file required")
min_n_ar      <- if (length(args) >= 3) as.integer(args[3]) else 10L
vcf_ld_file   <- if (length(args) >= 4) args[4] else stop("vcf_ld_file required")
cluster_dir   <- if (length(args) >= 5) args[5] else file.path(out_dir, "per_cluster")
threshold_tag <- if (length(args) >= 6) args[6] else "T??"

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# Load population map; derive population list
# -------------------------------------------------------------------
map  <- read.table(map_file, header = FALSE, sep = "\t",
                   col.names = c("sample_id", "population"),
                   stringsAsFactors = FALSE)
pops <- sort(unique(map$population))
cat(sprintf("INFO: %s | %d individuals | %d clusters: %s\n",
            threshold_tag, nrow(map), length(pops), paste(pops, collapse = ", ")))

# -------------------------------------------------------------------
# Helper: safe vcftools output reader
# -------------------------------------------------------------------
read_vcftools <- function(path, ...) {
  if (file.exists(path) && file.info(path)$size > 0)
    tryCatch(read.table(path, header = TRUE, stringsAsFactors = FALSE, ...),
             error = function(e) NULL)
  else NULL
}

# -------------------------------------------------------------------
# STEP 1: hierfstat — Ho, He, FIS, Ar from LD-pruned VCF
# -------------------------------------------------------------------
cat(sprintf("INFO: loading LD-pruned VCF: %s\n", vcf_ld_file))
if (!file.exists(vcf_ld_file)) stop("LD-pruned VCF not found: ", vcf_ld_file)

x <- gaston::read.vcf(vcf_ld_file, max.snps = 1e6, verbose = FALSE)
cat(sprintf("INFO: VCF loaded: %d individuals x %d SNPs\n", nrow(x), ncol(x)))

pop_vector <- map$population[match(x@ped$id, map$sample_id)]
keep_mask  <- !is.na(pop_vector)
n_unknown  <- sum(!keep_mask)
if (n_unknown > 0)
  cat(sprintf("INFO: %d individuals not in cluster map — excluded from hierfstat\n", n_unknown))

x_hf       <- x[keep_mask, ]
pop_vector <- pop_vector[keep_mask]
pop_levels <- sort(unique(pop_vector))
pop_factor <- factor(pop_vector, levels = pop_levels)

cat("INFO: converting genotypes to hierfstat format...\n")
geno_dos <- as.matrix(x_hf)
geno_hf  <- matrix(NA_integer_, nrow = nrow(geno_dos), ncol = ncol(geno_dos))
geno_hf[geno_dos == 0] <- 11L
geno_hf[geno_dos == 1] <- 12L
geno_hf[geno_dos == 2] <- 22L
hf_df <- cbind(data.frame(pop = pop_factor), as.data.frame(geno_hf))
rm(geno_dos, x, x_hf); gc()

cat("INFO: running hierfstat::basic.stats()...\n")
bs <- hierfstat::basic.stats(hf_df)
ho_per_pop  <- colMeans(bs$Ho,  na.rm = TRUE)
he_per_pop  <- colMeans(bs$Hs,  na.rm = TRUE)
fis_per_pop <- colMeans(bs$Fis, na.rm = TRUE)

cat(sprintf("INFO: running hierfstat::boot.ppfis() (nboot=500)...\n"))
fis_boot <- hierfstat::boot.ppfis(hf_df, nboot = 500)
fis_ll   <- as.numeric(unlist(fis_boot$ll))
fis_ul   <- as.numeric(unlist(fis_boot$ul))
if (length(fis_ll) != length(pop_levels) || length(fis_ul) != length(pop_levels)) {
  cat(sprintf("WARN: boot.ppfis CI extraction failed — attempting fallback\n"))
  if (!is.null(fis_boot$fis)) {
    fis_mat <- as.matrix(fis_boot$fis)
    if (ncol(fis_mat) == 2 && nrow(fis_mat) == length(pop_levels)) {
      fis_ll <- fis_mat[, 1]; fis_ul <- fis_mat[, 2]
    } else if (ncol(fis_mat) == length(pop_levels)) {
      fis_ll <- apply(fis_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
      fis_ul <- apply(fis_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
    } else {
      fis_ll <- rep(NA_real_, length(pop_levels))
      fis_ul <- rep(NA_real_, length(pop_levels))
    }
  } else {
    fis_ll <- rep(NA_real_, length(pop_levels))
    fis_ul <- rep(NA_real_, length(pop_levels))
  }
}
fis_sig <- !(fis_ll <= 0 & fis_ul >= 0)
fis_sig[is.na(fis_ll) | is.na(fis_ul)] <- NA

cat(sprintf("INFO: running hierfstat::allelic.richness(min.n=%d)...\n", min_n_ar))
ar_res     <- hierfstat::allelic.richness(hf_df, min.n = min_n_ar)
ar_per_pop <- colMeans(ar_res$Ar, na.rm = TRUE)

hf_metrics <- data.frame(
  population  = pop_levels,
  Ho          = as.numeric(ho_per_pop),
  He          = as.numeric(he_per_pop),
  FIS         = as.numeric(fis_per_pop),
  FIS_lower95 = round(fis_ll, 6),
  FIS_upper95 = round(fis_ul, 6),
  FIS_sig     = fis_sig,
  Ar          = as.numeric(ar_per_pop),
  stringsAsFactors = FALSE
)
rm(hf_df, bs, ar_res, fis_boot); gc()
cat("INFO: hierfstat metrics done\n")

# -------------------------------------------------------------------
# STEP 2: pi and Tajima's D from vcftools outputs
# -------------------------------------------------------------------
vcftools_metrics <- lapply(pops, function(pop) {
  pop_d  <- file.path(cluster_dir, pop)
  n_ind  <- sum(map$population == pop)

  pi_dat  <- read_vcftools(file.path(pop_d, "pi.sites.pi"))
  pi_mean <- if (!is.null(pi_dat) && "PI" %in% colnames(pi_dat))
    mean(pi_dat$PI, na.rm = TRUE) else NA_real_

  taj_dat  <- read_vcftools(file.path(pop_d, "tajima.Tajima.D"))
  taj_vals <- if (!is.null(taj_dat) && "TajimaD" %in% colnames(taj_dat))
    taj_dat$TajimaD[!is.na(taj_dat$TajimaD)] else numeric(0)
  taj_mean <- if (length(taj_vals) > 0) mean(taj_vals) else NA_real_
  taj_pval <- if (length(taj_vals) >= 5)
    tryCatch(wilcox.test(taj_vals, mu = 0, exact = FALSE)$p.value,
             error = function(e) NA_real_) else NA_real_

  data.frame(population = pop, n_individuals = n_ind,
             pi_mean = pi_mean, TajimaD_mean = taj_mean,
             TajimaD_pval = taj_pval, stringsAsFactors = FALSE)
})
vcftools_df <- do.call(rbind, vcftools_metrics)

# -------------------------------------------------------------------
# STEP 3: Private alleles (count) and rarefied PAR
# -------------------------------------------------------------------
cat("INFO: computing private alleles and rarefied PAR...\n")
g_par <- 2L * min_n_ar   # gene copies for rarefaction

# Load frequency tables: locus_id, N_CHR, n_allele (alt allele count)
freq_tables <- lapply(pops, function(pop) {
  freq_f <- file.path(cluster_dir, pop, "freq.frq")
  if (!file.exists(freq_f) || file.info(freq_f)$size == 0) return(NULL)
  d <- tryCatch(
    read.table(freq_f, header = FALSE, skip = 1, sep = "\t",
               col.names = c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ1", "FREQ2"),
               stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(d)) return(NULL)
  d$locus_id <- paste0(d$CHROM, ":", d$POS)
  d$N_CHR    <- as.integer(d$N_CHR)
  d$n_allele <- as.integer(round(as.numeric(d$FREQ2) * d$N_CHR))
  d[, c("locus_id", "N_CHR", "n_allele")]
})
names(freq_tables) <- pops

valid_pops <- pops[!sapply(freq_tables, is.null)]

if (length(valid_pops) >= 2) {
  # Merge tables across populations
  parts <- lapply(valid_pops, function(p) {
    d <- freq_tables[[p]]
    colnames(d)[colnames(d) == "N_CHR"]    <- paste0("N_", p)
    colnames(d)[colnames(d) == "n_allele"] <- paste0("n_", p)
    d
  })
  mg <- Reduce(function(a, b) merge(a, b, by = "locus_id", all = TRUE), parts)

  private_counts <- integer(length(valid_pops))
  par_rarefied   <- numeric(length(valid_pops))

  for (pi_i in seq_along(valid_pops)) {
    pop        <- valid_pops[pi_i]
    n_col      <- paste0("n_", pop)
    N_col      <- paste0("N_", pop)
    other_n    <- paste0("n_", setdiff(valid_pops, pop))

    own_n      <- mg[[n_col]]
    other_mat  <- mg[, other_n, drop = FALSE]

    # Private allele: present in this pop (n > 0), absent in all others (n = 0 or NA)
    is_private <- !is.na(own_n) & own_n > 0 &
                  rowSums(!is.na(other_mat) & other_mat > 0, na.rm = FALSE) == 0
    is_private[is.na(is_private)] <- FALSE

    private_counts[pi_i] <- sum(is_private)

    # Rarefied PAR — Kalinowski (2004) rarefaction formula
    # Expected number of private alleles in a random draw of g gene copies
    n_a <- own_n[is_private]
    N_p <- mg[[N_col]][is_private]

    # Keep only loci where total gene copies >= g (rarefaction possible)
    ok  <- !is.na(N_p) & N_p >= g_par & !is.na(n_a) & n_a > 0
    n_a <- n_a[ok]
    N_p <- N_p[ok]

    if (length(n_a) == 0) {
      par_rarefied[pi_i] <- 0
    } else {
      contrib <- vapply(seq_along(n_a), function(j) {
        ni <- n_a[j]; Ni <- N_p[j]
        if (Ni - ni < g_par) return(1)   # all non-private copies < g → guaranteed in sample
        log_r <- lchoose(Ni - ni, g_par) - lchoose(Ni, g_par)
        if (is.nan(log_r) || is.infinite(log_r)) return(1)
        1 - exp(log_r)
      }, numeric(1))
      par_rarefied[pi_i] <- sum(contrib, na.rm = TRUE)
    }

    cat(sprintf("  %s: %d private alleles | PAR(g=%d) = %.2f\n",
                pop, private_counts[pi_i], g_par, par_rarefied[pi_i]))
  }

  private_df <- data.frame(
    population      = valid_pops,
    private_alleles = private_counts,
    PAR             = round(par_rarefied, 4),
    stringsAsFactors = FALSE
  )
} else {
  private_df <- data.frame(
    population      = pops,
    private_alleles = NA_integer_,
    PAR             = NA_real_,
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------
# STEP 4: Merge all metrics into final summary table
# -------------------------------------------------------------------
final_df <- Reduce(function(a, b) merge(a, b, by = "population", all = TRUE),
                   list(vcftools_df, hf_metrics, private_df))

num_cols <- c("pi_mean", "Ho", "He", "FIS", "FIS_lower95", "FIS_upper95",
              "TajimaD_mean", "TajimaD_pval", "Ar", "PAR")
num_cols <- num_cols[num_cols %in% colnames(final_df)]
final_df[num_cols] <- lapply(final_df[num_cols], function(x) round(as.numeric(x), 6))

if ("TajimaD_pval" %in% colnames(final_df))
  final_df$TajimaD_sig <- !is.na(final_df$TajimaD_pval) & final_df$TajimaD_pval < 0.05
final_df <- final_df[order(final_df$population), ]

out_tsv <- file.path(out_dir,
                     sprintf("07b.3-diversity_metrics_%s_summary.tsv", threshold_tag))
write.table(final_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary table written:", out_tsv, "\n")
print(final_df)

# -------------------------------------------------------------------
# STEP 5: Barplots per metric
# -------------------------------------------------------------------
metrics_to_plot <- list(
  list(col = "pi_mean",         label = "Mean nucleotide diversity (π)",
       sig_col = NULL),
  list(col = "Ho",              label = "Observed heterozygosity (Ho)",
       sig_col = NULL),
  list(col = "He",              label = "Expected heterozygosity (He)",
       sig_col = NULL),
  list(col = "FIS",             label = "Inbreeding coefficient (FIS)",
       sig_col = "FIS_sig"),
  list(col = "TajimaD_mean",    label = "Mean Tajima's D",
       sig_col = "TajimaD_sig"),
  list(col = "Ar",              label = sprintf("Allelic richness rarefied (n=%d ind.)", min_n_ar),
       sig_col = NULL),
  list(col = "private_alleles", label = "Private alleles (count)",
       sig_col = NULL),
  list(col = "PAR",             label = sprintf("Private allelic richness rarefied (g=%d gene copies)", g_par),
       sig_col = NULL)
)

for (m in metrics_to_plot) {
  col     <- m$col
  label   <- m$label
  sig_col <- m$sig_col
  df_p    <- final_df[!is.na(final_df[[col]]), ]
  if (nrow(df_p) == 0) next

  df_p    <- df_p[order(-df_p[[col]]), ]
  df_p$population <- factor(df_p$population, levels = df_p$population)
  df_p$.val_label <- sprintf("%.4g", df_p[[col]])
  df_p$.vjust_val <- ifelse(df_p[[col]] >= 0, -0.4, 1.4)

  p <- ggplot(df_p, aes(x = population, y = .data[[col]], fill = population)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = .val_label, vjust = .vjust_val), size = 3) +
    scale_fill_manual(values = POP_COLORS, na.value = "grey50") +
    labs(title    = sprintf("%s — %s", label, threshold_tag),
         subtitle = sprintf("purity-filtered individuals (n=%d)", nrow(map)),
         x = NULL, y = label) +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1, size = 9),
      plot.title      = element_text(face = "bold"),
      legend.position = "none"
    )

  if (!is.null(sig_col) && sig_col %in% colnames(df_p)) {
    df_sig <- df_p[!is.na(df_p[[sig_col]]) & df_p[[sig_col]], ]
    if (nrow(df_sig) > 0) {
      sig_offset <- diff(range(df_p[[col]], na.rm = TRUE)) * 0.06
      df_sig$.star_y     <- ifelse(df_sig[[col]] >= 0,
                                   df_sig[[col]] + sig_offset,
                                   df_sig[[col]] - sig_offset)
      df_sig$.star_vjust <- ifelse(df_sig[[col]] >= 0, 0, 1)
      p <- p + geom_text(data = df_sig,
                         aes(x = population, y = .star_y, vjust = .star_vjust),
                         label = "*", size = 5, color = "red", inherit.aes = FALSE)
    }
    sig_subtitle <- if (sig_col == "FIS_sig")
      "* p < 0.05 (bootstrap 95% CI excludes 0, nboot = 500)"
    else if (sig_col == "TajimaD_sig")
      "* p < 0.05 (Wilcoxon signed-rank test vs. 0)"
    else NULL
    if (!is.null(sig_subtitle))
      p <- p + labs(subtitle = sig_subtitle) +
               theme(plot.subtitle = element_text(size = 7.5, color = "red"))
  }

  fn <- file.path(plot_dir,
                  sprintf("%s_%s.png", gsub("[^A-Za-z0-9]", "_", col), threshold_tag))
  ggsave(fn, p, width = 6, height = 5, dpi = 300)
  cat("Plot written:", fn, "\n")
}

cat(sprintf("DONE 07b.3 diversity metrics — %s\n", threshold_tag))

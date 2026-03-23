#!/usr/bin/env Rscript
# =============================================================================
# Script : 07.6-amova_cp.R
# Description : AMOVA (Analysis of Molecular Variance) for Dicorynia guianensis
#               chloroplast data. Partitions genetic variance at two levels:
#               - 2-level : Among sites / Within sites  → PhiST
#               - 3-level : Among haplogroups / Among sites within groups /
#                           Within sites                → PhiCT, PhiSC, PhiST
#               Haplogroups defined from haplotype network, ML tree + pairwise FST:
#               West        = Acarouany, Apatou, Maripasoula (refugial core)
#               Center-East = all remaining FG sites (incl. Chutes_Voltaires, Saul)
# Author  : Julien Bonnier
# Usage   : Rscript --vanilla 07.6-amova_cp.R \
#               <alignment_fasta> <metadata_csv> <output_dir>
# =============================================================================

.libPaths(c(path.expand("~/work/R"), .libPaths()))

suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
})

theme_set(theme_minimal(base_size = 11) + theme(
  plot.background  = element_rect(fill = "white", color = NA),
  panel.background = element_rect(fill = "white", color = NA)
))

set.seed(42)

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: script.R <aln> <metadata> <outdir>")

aln_file   <- args[1]
meta_file  <- args[2]
output_dir <- args[3]

cat("INFO: alignment   =", aln_file, "\n")
cat("INFO: metadata    =", meta_file, "\n")
cat("INFO: output_dir  =", output_dir, "\n")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

N_PERM <- 9999

# ── Haplogroup assignments ─────────────────────────────────────────────────────
# Defined from haplotype network topology, ML tree, and pairwise FST matrix.
# Option A (2 groups — corrected):
#   West        = Acarouany, Apatou, Maripasoula (refugial core, positive D/Fs,
#                 all private haplotypes, Suriname corridor)
#   Center_East = all remaining FG sites, including:
#                 - Chutes_Voltaires (FST/Trinité = 0.053 vs FST/Acarouany = 0.370
#                   → belongs to CE lineage despite geographic proximity to West)
#                 - Saul (FST pattern closer to CE than West)
#                 - Trinité (central node of star network, CE ancestral haplotypes)
HAPLOGROUPS <- list(
  West        = c("Acarouany", "Apatou", "Maripasoula"),
  Center_East = c("Angela", "Cacao", "Chutes_Voltaires", "Crique_Tigre",
                  "Foret_Regina_St_Georges", "MC_87", "MC_88", "Montagne_Fer",
                  "Nouragues_Inselberg", "Petit_croissant", "Piste_St_Elie",
                  "Regina", "Saut_Lavilette", "Saut_Takari", "St_georges",
                  "Saul", "Trinité")
)

site_to_group <- stack(HAPLOGROUPS) %>%
  rename(site = values, haplogroup = ind) %>%
  mutate(haplogroup = as.character(haplogroup))

# ── Read alignment ─────────────────────────────────────────────────────────────
cat("INFO: reading alignment...\n")
aln <- read.dna(aln_file, format = "fasta")
cat("INFO:", nrow(aln), "sequences x", ncol(aln), "bp\n")

seg_idx <- seg.sites(aln)
cat("INFO:", length(seg_idx), "segregating sites retained\n")
aln_seg <- aln[, seg_idx]

# ── Read metadata ──────────────────────────────────────────────────────────────
meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(sample_id, site) %>%
  distinct() %>%
  mutate(site = ifelse(site == "Herbier", "D.paranensis", site)) %>%
  filter(!site %in% c("Cameroun_Benin", "D.paranensis")) %>%
  filter(sample_id %in% rownames(aln_seg))

meta <- meta %>% left_join(site_to_group, by = "site")

n_unassigned <- sum(is.na(meta$haplogroup))
if (n_unassigned > 0) {
  cat("WARNING:", n_unassigned, "samples without haplogroup — will be excluded from 3-level AMOVA:\n")
  print(meta %>% filter(is.na(haplogroup)) %>% select(sample_id, site))
}

cat("INFO:", nrow(meta), "FG samples |",
    n_distinct(meta$site), "sites |",
    n_distinct(na.omit(meta$haplogroup)), "haplogroups\n")

# ── Helper: extract Phi statistics and p-values from pegas amova object ────────
extract_amova <- function(amova_obj, model_label) {
  vc    <- amova_obj$varcomp
  rows  <- rownames(vc)
  total <- sum(vc[, "sigma2"])
  # In 3-level AMOVA, pegas names the within-group site level "site" — disambiguate
  if (model_label == "3-level") rows[rows == "site"] <- "site within haplogroup"

  data.frame(
    model        = model_label,
    source       = rows,
    sigma2       = round(vc[, "sigma2"], 8),
    percent_var  = round(vc[, "sigma2"] / total * 100, 3),
    p_value      = if ("P.value" %in% colnames(vc)) round(vc[, "P.value"], 4) else NA_real_,
    stringsAsFactors = FALSE
  )
}

compute_phi <- function(amova_obj, model) {
  # pegas does not store Phi in the object — compute from varcomp sigma2 values.
  # trimws() removes padding spaces that pegas adds to row names; direct [row,col]
  # indexing avoids named-vector issues with round().
  vc     <- amova_obj$varcomp
  rnames <- trimws(rownames(vc))
  rownames(vc) <- rnames

  total <- sum(vc[, "sigma2"])

  get_s <- function(name) as.numeric(vc[name, "sigma2"])
  get_p <- function(name) {
    if (!"P.value" %in% colnames(vc)) return(NA_real_)
    as.numeric(vc[name, "P.value"])
  }

  if (model == "2-level") {
    phi_st <- get_s("site") / total
    data.frame(
      statistic      = "PhiST",
      value          = round(phi_st, 5),
      p_value        = round(get_p("site"), 4),
      model          = "2-level",
      interpretation = "Total differentiation among sites",
      stringsAsFactors = FALSE
    )
  } else {
    s_grp <- get_s("haplogroup")
    s_pop <- get_s("site")
    s_err <- get_s("Error")
    phi_ct <- s_grp / total
    phi_sc <- s_pop / (s_pop + s_err)
    phi_st <- (s_grp + s_pop) / total
    data.frame(
      statistic = c("PhiCT", "PhiSC", "PhiST"),
      value     = round(c(phi_ct, phi_sc, phi_st), 5),
      p_value   = round(c(get_p("haplogroup"), get_p("site"), NA_real_), 4),
      model     = "3-level",
      interpretation = c(
        "Differentiation among haplogroups",
        "Differentiation among sites within haplogroups",
        "Total differentiation"
      ),
      stringsAsFactors = FALSE
    )
  }
}

# ── 2-level AMOVA: Among sites / Within sites ──────────────────────────────────
cat("INFO: running 2-level AMOVA (", N_PERM, "permutations)...\n")

samps_2  <- meta$sample_id
dist_2   <- dist.dna(aln_seg[samps_2, , drop = FALSE],
                     model = "raw", pairwise.deletion = TRUE)
strata_2 <- data.frame(
  site = factor(meta$site[match(samps_2, meta$sample_id)]),
  row.names = samps_2
)

amova_2lvl <- amova(dist_2 ~ site, data = strata_2, nperm = N_PERM)

cat("\n── 2-level AMOVA result ──────────────────────────────────────────────\n")
print(amova_2lvl)

# ── 3-level AMOVA: Among haplogroups / Among sites / Within ───────────────────
cat("\nINFO: running 3-level AMOVA (", N_PERM, "permutations)...\n")

meta_3  <- meta %>% filter(!is.na(haplogroup))
samps_3 <- meta_3$sample_id
dist_3  <- dist.dna(aln_seg[samps_3, , drop = FALSE],
                    model = "raw", pairwise.deletion = TRUE)
strata_3 <- data.frame(
  haplogroup = factor(meta_3$haplogroup[match(samps_3, meta_3$sample_id)]),
  site       = factor(meta_3$site[match(samps_3, meta_3$sample_id)]),
  row.names  = samps_3
)

amova_3lvl <- amova(dist_3 ~ haplogroup/site, data = strata_3, nperm = N_PERM)

cat("\n── 3-level AMOVA result ──────────────────────────────────────────────\n")
print(amova_3lvl)

# ── Compile results ────────────────────────────────────────────────────────────
vc_table <- bind_rows(
  extract_amova(amova_2lvl, "2-level"),
  extract_amova(amova_3lvl, "3-level")
)

phi_table <- bind_rows(
  compute_phi(amova_2lvl, "2-level"),
  compute_phi(amova_3lvl, "3-level")
) %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    p_value < 0.1   ~ ".",
    TRUE            ~ "ns"
  ))

write.table(vc_table,  file.path(output_dir, "cp_amova_variance_components.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(phi_table, file.path(output_dir, "cp_amova_phi_statistics.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("INFO: saved cp_amova_variance_components.tsv\n")
cat("INFO: saved cp_amova_phi_statistics.tsv\n")

# ── Plot 1: Variance partition stacked barplot ─────────────────────────────────
cat("INFO: generating variance partition plot...\n")

source_levels <- c("Among haplogroups", "Among sites within haplogroups",
                   "Among sites", "Within sites")
source_colors <- c(
  "Among haplogroups"              = "#D73027",
  "Among sites within haplogroups" = "#FC8D59",
  "Among sites"                    = "#4575B4",
  "Within sites"                   = "#91BFDB"
)

plot_vc <- vc_table %>%
  mutate(
    source = recode(source,
      "haplogroup"              = "Among haplogroups",
      "site within haplogroup"  = "Among sites within haplogroups",
      "site"                    = "Among sites",
      "Error"              = "Within sites"
    ),
    source = factor(source, levels = rev(source_levels)),
    model  = factor(model, levels = c("2-level", "3-level"),
                    labels = c("2-level\n(sites only)", "3-level\n(haplogroups / sites)"))
  ) %>%
  filter(percent_var >= 0)

p_var <- ggplot(plot_vc, aes(x = model, y = percent_var, fill = source)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(sprintf("%.1f", percent_var), "%")),
            position = position_stack(vjust = 0.5),
            size = 3.4, color = "white", fontface = "bold") +
  scale_fill_manual(values = source_colors, name = "Variance source",
                    guide  = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast AMOVA"),
    subtitle = paste0("Partitioning of molecular variance | ", N_PERM, " permutations"),
    x        = "AMOVA model",
    y        = "% Total variance"
  ) +
  theme(
    legend.position  = "right",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "cp_amova_variance_partition.png"),
       p_var, width = 9, height = 6, dpi = 300)
cat("INFO: saved cp_amova_variance_partition.png\n")

# ── Plot 2: Phi statistics barplot ────────────────────────────────────────────
cat("INFO: generating Phi statistics plot...\n")

phi_plot <- phi_table %>%
  mutate(
    statistic = factor(statistic, levels = c("PhiST", "PhiCT", "PhiSC")),
    model     = factor(model, levels = c("2-level", "3-level"))
  )

p_phi <- ggplot(phi_plot, aes(x = statistic, y = value, fill = model)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(aes(label = significance), vjust = -0.5, size = 5.5, fontface = "bold") +
  geom_text(aes(label = sprintf("%.4f", value)),
            vjust = 1.8, size = 3.4, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("2-level" = "#4575B4", "3-level" = "#D73027"),
                    name = "Model") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, 1)) +
  labs(
    title    = expression(italic("Dicorynia guianensis") ~ "— Chloroplast" ~ Phi ~ "statistics"),
    subtitle = paste0("2-level and 3-level AMOVA | ", N_PERM, " permutations"),
    caption  = "Significance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns = not significant",
    x        = expression(Phi ~ "statistic"),
    y        = "Value"
  ) +
  theme(
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir, "cp_amova_phi_stats.png"),
       p_phi, width = 8, height = 6, dpi = 300)
cat("INFO: saved cp_amova_phi_stats.png\n")

# ── Print summary ──────────────────────────────────────────────────────────────
cat("\n══ AMOVA SUMMARY ═════════════════════════════════════════════════════\n")
cat(sprintf("  Samples : %d | Sites : %d | Haplogroups : %d\n",
            nrow(meta), n_distinct(meta$site), n_distinct(na.omit(meta$haplogroup))))
cat(sprintf("  Permutations : %d\n\n", N_PERM))
cat("  2-level model:\n")
cat(sprintf("    PhiST = %.5f  (p = %.4f)\n",
            phi_table$value[phi_table$statistic == "PhiST" & phi_table$model == "2-level"],
            phi_table$p_value[phi_table$statistic == "PhiST" & phi_table$model == "2-level"]))
cat("  3-level model:\n")
for (stat in c("PhiCT", "PhiSC", "PhiST")) {
  row <- phi_table[phi_table$statistic == stat & phi_table$model == "3-level", ]
  cat(sprintf("    %-6s = %.5f  (p = %.4f)  %s\n",
              stat, row$value, row$p_value, row$significance))
}
cat("══════════════════════════════════════════════════════════════════════\n\n")

cat("DONE AMOVA complete — outputs saved to", output_dir, "\n")

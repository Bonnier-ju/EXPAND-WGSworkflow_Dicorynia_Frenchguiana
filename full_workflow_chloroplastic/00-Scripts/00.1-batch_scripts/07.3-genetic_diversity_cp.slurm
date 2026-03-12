#!/bin/bash
#SBATCH -p workq
#SBATCH --time=01:00:00
#SBATCH -J cp-07.3-diversity
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load statistics/R/4.3.0

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
WORKFLOW="workflow_full"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
INPUT_ALN="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/alignments_phylo/alignment_trimmed.fasta"
HAP_TABLE="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/haplotype_tables/haplotype_table.tsv"
METADATA="${BASE_DIR}/metadata/sample_sheet_chloroplastic_complete_filtered.csv"
OUTPUT_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/07.3-genetic_diversity"
R_SCRIPT="${BASE_DIR}/00-Scripts/00.2-R_scripts/07.3-genetic_diversity_cp.R"

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$INPUT_ALN" ];  then echo "ERROR: alignment not found: $INPUT_ALN" >&2; exit 1; fi
if [ ! -s "$HAP_TABLE" ];  then echo "ERROR: haplotype table not found: $HAP_TABLE" >&2; exit 1; fi
if [ ! -s "$METADATA" ];   then echo "ERROR: metadata not found: $METADATA" >&2; exit 1; fi

echo "START $(date -Iseconds) genetic_diversity cp"

# -------------------------------------------------------------------
# STEP 2: Run R script
# -------------------------------------------------------------------
Rscript --vanilla "$R_SCRIPT" \
    "$INPUT_ALN" \
    "$HAP_TABLE" \
    "$METADATA" \
    "$OUTPUT_DIR"

if [ -s "${OUTPUT_DIR}/diversity_stats_per_site.tsv" ]; then
    n_sites=$(tail -n +2 "${OUTPUT_DIR}/diversity_stats_per_site.tsv" | wc -l)
    echo "DONE $(date -Iseconds) n_groups=${n_sites} outputs saved to ${OUTPUT_DIR}"
else
    echo "ERROR: R script produced no output" >&2; exit 1
fi

# -------------------------------------------------------------------
# Description
# Genetic diversity analysis of the chloroplast dataset.
# Inputs:
#   - alignment_trimmed.fasta : trimmed consensus alignment (from 06.6)
#   - haplotype_table.tsv     : per-sample haplotype assignment (from 06.4)
#   - sample_sheet_chloroplastic_complete_filtered.csv : metadata (164 samples)
# Statistics computed globally, per region, and per sampling site:
#   - n_samples       : number of individuals
#   - n_haplotypes    : number of unique haplotypes
#   - Hd              : haplotype diversity (Nei 1987)
#   - pi              : nucleotide diversity (mean pairwise differences / site)
#   - Tajima's D      : test of selective neutrality / demographic history
#   - n_private_haps  : haplotypes unique to one site
# Outputs (03-genetic_diversity/):
#   diversity_stats_per_site.tsv      : full stats table
#   cp_diversity_hd_per_site.png      : barplot of Hd per site
#   cp_diversity_pi_per_site.png      : barplot of pi per site
#   cp_diversity_hd_vs_pi.png         : Hd vs. pi scatter per site
#   cp_diversity_haplotype_counts.png : total vs. private haplotypes per site
# -------------------------------------------------------------------

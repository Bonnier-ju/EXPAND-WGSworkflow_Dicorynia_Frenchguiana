#!/bin/bash
#SBATCH -p workq
#SBATCH --time=01:00:00
#SBATCH -J cp-07.2-hapnet
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
METADATA="${BASE_DIR}/metadata/sample_sheet_chloroplastic_complete.csv"
OUTPUT_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/02-haplotype_networks"
R_SCRIPT="${BASE_DIR}/00-Scripts/00.2-R_scripts/07.2-haplotype_networks_cp.R"

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$INPUT_ALN" ];  then echo "ERROR: alignment not found: $INPUT_ALN" >&2; exit 1; fi
if [ ! -s "$HAP_TABLE" ];  then echo "ERROR: haplotype table not found: $HAP_TABLE" >&2; exit 1; fi
if [ ! -s "$METADATA" ];   then echo "ERROR: metadata not found: $METADATA" >&2; exit 1; fi

echo "START $(date -Iseconds) haplotype_networks cp"

# -------------------------------------------------------------------
# STEP 2: Run R script
# -------------------------------------------------------------------
Rscript --vanilla "$R_SCRIPT" \
    "$INPUT_ALN" \
    "$HAP_TABLE" \
    "$METADATA" \
    "$OUTPUT_DIR"

if ls "${OUTPUT_DIR}/"*.png &>/dev/null; then
    echo "DONE $(date -Iseconds) outputs saved to ${OUTPUT_DIR}"
else
    echo "ERROR: R script produced no output" >&2; exit 1
fi

# -------------------------------------------------------------------
# Description
# Haplotype network inference for the chloroplast dataset.
# Input: coordinate-aligned consensus sequences (160,511 bp, 180 samples)
# Method: minimum spanning network (pegas::haploNet, getProb=FALSE)
# Nodes: sized by haplotype frequency, pie charts showing site composition
# Outputs:
#   cp_hapnet_sites.png    : network colored by sampling site (all individuals)
#   cp_hapnet_regions.png  : network colored by region (FG vs outgroup)
#   cp_hapnet_fg_only.png  : network colored by site, French Guiana only (outgroup excluded)
#   cp_hapnet_treemutation.png : FG only, Treemutation samples highlighted with orange border
# -------------------------------------------------------------------

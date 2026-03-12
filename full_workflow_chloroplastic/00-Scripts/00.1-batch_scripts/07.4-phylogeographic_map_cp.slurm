#!/bin/bash
#SBATCH -p workq
#SBATCH --time=01:00:00
#SBATCH -J cp-07.4-map
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
HAP_TABLE="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/haplotype_tables/haplotype_table.tsv"
METADATA="${BASE_DIR}/metadata/sample_sheet_chloroplastic_complete_filtered.csv"
GEOLOC="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/geoloc_site.csv"
OUTPUT_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/07.4-phylogeographic_map"
R_SCRIPT="${BASE_DIR}/00-Scripts/00.2-R_scripts/07.4-phylogeographic_map_cp.R"

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$HAP_TABLE" ];  then echo "ERROR: haplotype table not found: $HAP_TABLE" >&2; exit 1; fi
if [ ! -s "$METADATA" ];   then echo "ERROR: metadata not found: $METADATA" >&2; exit 1; fi
if [ ! -s "$GEOLOC" ];     then echo "ERROR: geoloc file not found: $GEOLOC" >&2; exit 1; fi

echo "START $(date -Iseconds) phylogeographic_map cp"

# -------------------------------------------------------------------
# STEP 2: Run R script
# -------------------------------------------------------------------
Rscript --vanilla "$R_SCRIPT" \
    "$HAP_TABLE" \
    "$METADATA" \
    "$GEOLOC" \
    "$OUTPUT_DIR"

if [ -s "${OUTPUT_DIR}/cp_map_fg_haplotypes.png" ]; then
    echo "DONE $(date -Iseconds) outputs saved to ${OUTPUT_DIR}"
else
    echo "ERROR: R script produced no output" >&2; exit 1
fi

# -------------------------------------------------------------------
# Description
# Phylogeographic map of chloroplast haplotype distribution across
# French Guiana sampling sites for Dicorynia guianensis.
# Outgroup (Cameroun/Benin) excluded from the map.
# Inputs:
#   - haplotype_table.tsv    : per-sample haplotype assignment (from 06.4)
#   - sample_sheet_chloroplastic_complete_filtered.csv : metadata (164 samples)
#   - geoloc_site.csv        : GPS coordinates per sampling site
# Outputs (04-phylogeographic_map/):
#   cp_map_fg_haplotypes.png  : pie chart map — haplotype composition per site
#   cp_map_fg_diversity.png   : bubble map — Hd and N haplotypes per site
#   cp_map_site_summary.tsv   : per-site summary table (Hd, N haplotypes, coords)
# -------------------------------------------------------------------

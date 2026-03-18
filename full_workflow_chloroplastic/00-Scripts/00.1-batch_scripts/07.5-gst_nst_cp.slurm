#!/bin/bash
#SBATCH -p workq
#SBATCH --time=02:00:00
#SBATCH -J cp-07.5-gst-nst
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =============================================================================
# STEP 07.5 — GST vs NST phylogeographic structure test (Pons & Petit 1996)
# Inputs:
#   - alignment_trimmed.fasta  : trimmed chloroplast alignment (all samples)
#   - sample_sheet_*.csv       : metadata with site assignments
# Outputs (in OUTPUT_DIR):
#   - cp_gst_nst_summary.tsv       : GST, NST, permutation p-value
#   - cp_gst_nst_barplot.png       : GST vs NST observed values
#   - cp_gst_nst_permutation.png   : permutation distribution of NST-GST
#   - cp_gst_nst_site_diversity.png: per-site Hd and pi
# =============================================================================

WORKFLOW="workflow_full"
BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"

INPUT_ALN="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/alignments_phylo/alignment_trimmed.fasta"
META="${BASE_DIR}/metadata/sample_sheet_chloroplastic_complete_filtered.csv"
RSCRIPT="${BASE_DIR}/00-Scripts/00.2-R_scripts/07.5-gst_nst_cp.R"
OUTPUT_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/07.5-gst_nst"

mkdir -p "$OUTPUT_DIR"

echo "START $(date -Iseconds) GST/NST analysis"
echo "  alignment : $INPUT_ALN"
echo "  metadata  : $META"
echo "  output    : $OUTPUT_DIR"

# Check inputs
[ -s "$INPUT_ALN" ] || { echo "ERROR: alignment not found: $INPUT_ALN"; exit 1; }
[ -s "$META"      ] || { echo "ERROR: metadata not found: $META";        exit 1; }
[ -s "$RSCRIPT"   ] || { echo "ERROR: R script not found: $RSCRIPT";    exit 1; }

module load statistics/R/4.3.0

Rscript --vanilla "$RSCRIPT" \
    "$INPUT_ALN" \
    "$META" \
    "$OUTPUT_DIR"

# Check outputs
for f in cp_gst_nst_summary.tsv \
          cp_gst_nst_barplot.png \
          cp_gst_nst_permutation.png \
          cp_gst_nst_site_diversity.png; do
    if [ -s "${OUTPUT_DIR}/${f}" ]; then
        echo "  [OK] $f"
    else
        echo "  [MISSING] $f"
    fi
done

echo "DONE $(date -Iseconds) outputs saved to ${OUTPUT_DIR}"

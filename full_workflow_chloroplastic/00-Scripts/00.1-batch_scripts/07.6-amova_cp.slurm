#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J cp-07.6-amova
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =============================================================================
# STEP 07.6 — AMOVA (Analysis of Molecular Variance) — chloroplast cpDNA
# Inputs:
#   - alignment_trimmed.fasta : trimmed chloroplast alignment (all samples)
#   - sample_sheet_*.csv      : metadata with site assignments
# Outputs (in OUTPUT_DIR):
#   - cp_amova_variance_components.tsv : sigma2 and % variance per level
#   - cp_amova_phi_statistics.tsv      : PhiST, PhiCT, PhiSC with p-values
#   - cp_amova_variance_partition.png  : stacked barplot of variance partition
#   - cp_amova_phi_stats.png           : barplot of Phi statistics
#
# AMOVA models:
#   2-level : Among sites / Within sites               → PhiST
#   3-level : Among haplogroups / Among sites within   → PhiCT, PhiSC, PhiST
#             haplogroups / Within sites
#
# Haplogroups (from haplotype network + ML tree):
#   West        : Acarouany, Apatou, Chutes_Voltaires, Maripasoula, Saul
#   Center-East : all remaining FG sites
# =============================================================================

WORKFLOW="workflow_full"
BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"

INPUT_ALN="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/alignments_phylo/alignment_trimmed.fasta"
META="${BASE_DIR}/metadata/sample_sheet_chloroplastic_complete_filtered.csv"
RSCRIPT="${BASE_DIR}/00-Scripts/00.2-R_scripts/07.6-amova_cp.R"
OUTPUT_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/07.6-amova"

mkdir -p "$OUTPUT_DIR"

echo "START $(date -Iseconds) AMOVA"
echo "  alignment : $INPUT_ALN"
echo "  metadata  : $META"
echo "  output    : $OUTPUT_DIR"

# Check inputs
[ -s "$INPUT_ALN" ] || { echo "ERROR: alignment not found: $INPUT_ALN"; exit 1; }
[ -s "$META"      ] || { echo "ERROR: metadata not found: $META";        exit 1; }
[ -s "$RSCRIPT"   ] || { echo "ERROR: R script not found: $RSCRIPT";     exit 1; }

module load statistics/R/4.3.0

Rscript --vanilla "$RSCRIPT" \
    "$INPUT_ALN" \
    "$META" \
    "$OUTPUT_DIR"

# Check outputs
echo ""
echo "Output verification:"
for f in cp_amova_variance_components.tsv \
          cp_amova_phi_statistics.tsv \
          cp_amova_variance_partition.png \
          cp_amova_phi_stats.png; do
    if [ -s "${OUTPUT_DIR}/${f}" ]; then
        echo "  [OK] $f"
    else
        echo "  [MISSING] $f"
    fi
done

echo ""
echo "DONE $(date -Iseconds) outputs saved to ${OUTPUT_DIR}"

#!/bin/bash
#SBATCH -p workq
#SBATCH --time=12:00:00
#SBATCH -J 08.2-lfmm2
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# -------------------------------------------------------------------
# STEP 1: Define paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"

# Genotype matrix (155 ind x 433,798 SNPs, 0/1/2 imputed)
GENO_MATRIX="${BASE}/08-gea/workflow_full/08.2-genotype_matrix/genotype_012_imputed.tsv"
SNP_IDS="${BASE}/08-gea/workflow_full/08.2-genotype_matrix/snp_ids.txt"
SAMPLE_META="${BASE}/08-gea/workflow_full/08.2-genotype_matrix/sample_metadata.tsv"

# Final environmental variables selected by RDA forward selection + VIF (08.1d)
VARS_FINAL="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/rda_vif/variables_final.txt"

# Site-level environmental CSVs (08.0.x)
CHELSA_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/chelsa/chelsa_env_per_site.csv"
ENVIREM_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/envirem/envirem_env_per_site.csv"
MANUAL_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/manual_variables/manual_variables_per_site.csv"
TERRA_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/terraclimate/terraclimate_env_per_site.csv"

# Number of latent factors (K from fastSTRUCTURE)
K=4

# FDR threshold for candidate SNP identification
FDR=0.05

OUT_DIR="${BASE}/08-gea/workflow_full/08.2-lfmm2"
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.2-lfmm2.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.2-lfmm2.R"

mkdir -p "${OUT_DIR}/output" "${OUT_DIR}/plots" "${OUT_DIR}/input" "$LOG_DIR"
echo "START $(date -Iseconds) 08.2 LFMM2 GEA" > "$LOG_FILE"

# Check required inputs
for f in "$GENO_MATRIX" "$SNP_IDS" "$SAMPLE_META" "$VARS_FINAL" \
          "$CHELSA_CSV" "$ENVIREM_CSV" "$MANUAL_CSV" "$TERRA_CSV"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: required input missing or empty: $f" | tee -a "$LOG_FILE"
        exit 1
    fi
done

# -------------------------------------------------------------------
# STEP 2: R — LFMM2 analysis
# -------------------------------------------------------------------
echo "STEP 2: LFMM2 GEA (K=${K}, FDR=${FDR})" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$GENO_MATRIX" \
    "$SNP_IDS" \
    "$SAMPLE_META" \
    "$VARS_FINAL" \
    "$CHELSA_CSV" \
    "$ENVIREM_CSV" \
    "$MANUAL_CSV" \
    "$TERRA_CSV" \
    "${OUT_DIR}/output" \
    "${OUT_DIR}/plots" \
    "$K" \
    "$FDR" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.2 LFMM2 completed" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}"

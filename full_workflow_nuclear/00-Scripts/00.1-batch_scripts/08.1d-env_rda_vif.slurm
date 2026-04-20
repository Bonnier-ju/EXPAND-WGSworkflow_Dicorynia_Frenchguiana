#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J 08.1c-env_rda_vif
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=64G
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

# Genotype matrix from 08.2 (155 individuals x N SNPs, 0/1/2 imputed)
GENO_MATRIX="${BASE}/08-gea/workflow_full/08.2-genotype_matrix/genotype_012_imputed.tsv"
SAMPLE_META="${BASE}/08-gea/workflow_full/08.2-genotype_matrix/sample_metadata.tsv"

# Environmental data (site-level CSVs from 08.0.x)
CHELSA_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/chelsa/chelsa_env_per_site.csv"
ENVIREM_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/envirem/envirem_env_per_site.csv"
MANUAL_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/manual_variables/manual_variables_per_site.csv"

# Variables selected in 08.1a/08.1b (uncommented lines = kept)
VARS_FILE="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/correlation_analysis/variables_to_keep_template.txt"

OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/rda_vif"
LOG_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/logs"
LOG_FILE="${LOG_DIR}/08.1c-env_rda_vif.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1c-env_rda_vif.R"

# Number of permutations for forward selection significance test
N_PERM=999

mkdir -p "$OUT_DIR" "$LOG_DIR"
echo "START $(date -Iseconds) 08.1c RDA forward selection + VIF" > "$LOG_FILE"

# Check inputs
for f in "$GENO_MATRIX" "$SAMPLE_META" "$CHELSA_CSV" "$ENVIREM_CSV" "$MANUAL_CSV" "$VARS_FILE"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: required input missing: $f" | tee -a "$LOG_FILE"
        exit 1
    fi
done

# -------------------------------------------------------------------
# STEP 2: R — RDA forward selection + VIF
# -------------------------------------------------------------------
echo "STEP 2: RDA forward selection + VIF" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$GENO_MATRIX" \
    "$SAMPLE_META" \
    "$CHELSA_CSV" \
    "$ENVIREM_CSV" \
    "$MANUAL_CSV" \
    "$VARS_FILE" \
    "$OUT_DIR" \
    "$N_PERM" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.1c RDA + VIF completed" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}"

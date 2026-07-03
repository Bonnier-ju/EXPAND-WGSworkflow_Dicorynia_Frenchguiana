#!/bin/bash
#SBATCH -p workq
#SBATCH --time=00:30:00
#SBATCH -J 08.1b-env_pca
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
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

CHELSA_SITE_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/chelsa/chelsa_env_per_site.csv"
ENVIREM_SITE_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/envirem/envirem_env_per_site.csv"
MANUAL_SITE_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/manual_variables/manual_variables_per_site.csv"
TERRA_SITE_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/terraclimate/terraclimate_env_per_site.csv"

VARS_TO_KEEP="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/correlation_analysis/variables_to_keep.txt"

OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection"
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.1b-env_pca.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b-env_pca.R"

mkdir -p "$LOG_DIR"
echo "START $(date -Iseconds) 08.1b environmental PCA" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: R — PCA on selected variables, site-level plots
# -------------------------------------------------------------------
echo "STEP 2: R PCA on selected environmental variables" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$CHELSA_SITE_CSV" \
    "$ENVIREM_SITE_CSV" \
    "$MANUAL_SITE_CSV" \
    "$VARS_TO_KEEP" \
    "$OUT_DIR" \
    "$TERRA_SITE_CSV" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.1b env PCA completed" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}/plots/"

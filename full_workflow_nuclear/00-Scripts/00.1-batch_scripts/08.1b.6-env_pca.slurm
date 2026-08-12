#!/bin/bash
# 08.1b.6-env_pca.slurm
#
# PCA on the environmental variables Julien selected by editing
# correlation_analysis/variable_selection_template.txt (08.1b.5 output).
# Site scores coloured by the K=3 GEA cluster (cf. CLAUDE.md §4).
#
# Do NOT submit before the template has been edited - the R script will
# fail loudly (not silently) if no variable is uncommented.

#SBATCH -p workq
#SBATCH --time=00:15:00
#SBATCH -J 08.1b.6-env_pca
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ---- STEP 0: Environment ----
module purge
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# ---- STEP 1: Paths ----
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/08.1b-TerraClimate_var_selection"
CLUSTERS_CSV="${BASE}/metadata/sites_by_clusters.csv"
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.1b.6-env_pca.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b.6-env_pca.R"

mkdir -p "${LOG_DIR}"

# ---- Guard clause ----
for f in "${OUT_DIR}/correlation_analysis/env_variables_merged_31vars.csv" \
         "${OUT_DIR}/correlation_analysis/variable_selection_template.txt" \
         "${CLUSTERS_CSV}" "${R_SCRIPT}"; do
  [ -s "$f" ] || { echo "ERROR: required input missing or empty: $f"; exit 1; }
done

echo "START $(date -Iseconds) 08.1b.6 env PCA" | tee "${LOG_FILE}"

# ---- STEP 2: Run ----
Rscript --vanilla "${R_SCRIPT}" "${OUT_DIR}" "${CLUSTERS_CSV}" \
  >> "${LOG_FILE}" 2>&1

echo "DONE $(date -Iseconds) 08.1b.6 env PCA" | tee -a "${LOG_FILE}"
echo "Output: ${OUT_DIR}/pca/"

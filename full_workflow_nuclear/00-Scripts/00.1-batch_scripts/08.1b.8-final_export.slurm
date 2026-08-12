#!/bin/bash
# 08.1b.8-final_export.slurm
#
# Final export of the refonte 08.1 pipeline: individual-level environmental
# matrix (154 ind x final variables) + markdown summary report built from
# the real outputs of 08.1b.1-08.1b.7. This is the last script of the 08.1b
# sub-pipeline - its output (env_matrix_final_154ind.tsv) is the official
# input for 08.2 (LFMM2) onward.
#
# Lightweight: just merges small tables already computed.

#SBATCH -p workq
#SBATCH --time=00:10:00
#SBATCH -J 08.1b.8-final_export
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=2G
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
META_FILE="${BASE}/08-gea/workflow_full/08.2-genotype_matrix/sample_metadata.tsv"
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.1b.8-final_export.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b.8-final_export.R"

mkdir -p "${LOG_DIR}"

# ---- Guard clause ----
for f in "${OUT_DIR}/rda_vif/variables_final.txt" \
         "${OUT_DIR}/correlation_analysis/env_variables_merged_31vars.csv" \
         "${META_FILE}" "${R_SCRIPT}"; do
  [ -s "$f" ] || { echo "ERROR: required input missing or empty: $f"; exit 1; }
done

echo "START $(date -Iseconds) 08.1b.8 final export" | tee "${LOG_FILE}"

# ---- STEP 2: Run ----
Rscript --vanilla "${R_SCRIPT}" "${OUT_DIR}" "${META_FILE}" >> "${LOG_FILE}" 2>&1

echo "DONE $(date -Iseconds) 08.1b.8 final export" | tee -a "${LOG_FILE}"
echo "Output: ${OUT_DIR}/final/"
echo ""
echo ">>> This is the official env-variable input for 08.2 (LFMM2) onward. <<<"

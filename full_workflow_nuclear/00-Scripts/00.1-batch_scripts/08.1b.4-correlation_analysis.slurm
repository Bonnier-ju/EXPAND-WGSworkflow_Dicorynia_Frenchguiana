#!/bin/bash
# 08.1b.4-correlation_analysis.slurm
#
# Merge the 10 custom + 19 BIO + 2 topo candidate variables (31 total) and
# compute pairwise Pearson correlations, flagging |r|>=0.70 pairs ahead of
# the manual selection step (08.1b.5).
#
# Lightweight: 19 sites x 31 variables, base R + ggplot2 only.

#SBATCH -p workq
#SBATCH --time=00:20:00
#SBATCH -J 08.1b.4-correlation_analysis
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
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.1b.4-correlation_analysis.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b.4-correlation_analysis.R"

CORR_THRESHOLD=0.70

mkdir -p "${LOG_DIR}"

# ---- Guard clause: verify inputs before launching ----
for f in "${OUT_DIR}/data/terraclimate_derived_variables_19sites.csv" \
         "${OUT_DIR}/data/bioclim_variables_19sites.csv" \
         "${OUT_DIR}/data/topo_variables_19sites.csv" \
         "${R_SCRIPT}"; do
  [ -s "$f" ] || { echo "ERROR: required input missing or empty: $f"; exit 1; }
done

echo "START $(date -Iseconds) 08.1b.4 correlation analysis" | tee "${LOG_FILE}"

# ---- STEP 2: Run ----
Rscript --vanilla "${R_SCRIPT}" "${OUT_DIR}" "${CORR_THRESHOLD}" \
  >> "${LOG_FILE}" 2>&1

echo "DONE $(date -Iseconds) 08.1b.4 correlation analysis" | tee -a "${LOG_FILE}"
echo "Output: ${OUT_DIR}/correlation_analysis/"

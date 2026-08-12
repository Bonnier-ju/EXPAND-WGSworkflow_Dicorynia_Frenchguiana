#!/bin/bash
# 08.1b.5-manual_selection_prep.slurm
#
# Semi-manual step (cf. CLAUDE.md §0): prepares correlation groups, CV, and
# descriptions for the 31 candidate variables so Julien can pick ~8-10 final
# variables by editing the generated template. Does NOT select variables
# itself. STOP after this job - do not submit 08.1b.6 until the edited
# variable_selection_template.txt is ready.
#
# Trivial computation, no heavy resources needed.

#SBATCH -p workq
#SBATCH --time=00:10:00
#SBATCH -J 08.1b.5-manual_selection_prep
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
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.1b.5-manual_selection_prep.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b.5-manual_selection_prep.R"

mkdir -p "${LOG_DIR}"

# ---- Guard clause ----
for f in "${OUT_DIR}/correlation_analysis/variable_summary.tsv" \
         "${OUT_DIR}/correlation_analysis/correlated_pairs_table.tsv" \
         "${R_SCRIPT}"; do
  [ -s "$f" ] || { echo "ERROR: required input missing or empty: $f"; exit 1; }
done

echo "START $(date -Iseconds) 08.1b.5 manual selection prep" | tee "${LOG_FILE}"

# ---- STEP 2: Run ----
Rscript --vanilla "${R_SCRIPT}" "${OUT_DIR}" >> "${LOG_FILE}" 2>&1

echo "DONE $(date -Iseconds) 08.1b.5 manual selection prep" | tee -a "${LOG_FILE}"
echo "Output: ${OUT_DIR}/correlation_analysis/variable_selection_template.txt"
echo ""
echo ">>> STOP: edit variable_selection_template.txt before running 08.1b.6 <<<"

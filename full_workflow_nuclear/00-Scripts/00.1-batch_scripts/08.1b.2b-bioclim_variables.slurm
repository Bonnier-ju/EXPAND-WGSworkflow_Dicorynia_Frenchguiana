#!/bin/bash
# 08.1b.2b-bioclim_variables.slurm
#
# Compute BIO1-BIO19 from the TerraClimate 1970-2000 monthly climatology,
# using dismo::biovars() (advised by Sylvain Schmitt, 2026-08-11).
# Recomputes CHELSA-equivalent bioclim variables from TerraClimate alone -
# avoids mixing CHELSA + TerraClimate as separate sources for the same
# variable type. See 08.1b.2b-bioclim_variables.R header for details.
#
# Requires the R package 'dismo' - NOT yet confirmed available in
# R_LIBS_USER on Genotoul. Check with `search_R_package dismo` before
# submitting; if absent, install into the project renv first (CLAUDE.md §5).
#
# Lightweight: 19 sites x 12 months, base R + one function call.

#SBATCH -p workq
#SBATCH --time=00:30:00
#SBATCH -J 08.1b.2b-bioclim_variables
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
# Combined: packages are split across two personal libraries on this account
# (~/R/R and the R-default ~/R/x86_64-pc-linux-gnu-library/4.3) - see CLAUDE.md §5.
export R_LIBS_USER="${HOME}/R/R:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# ---- STEP 1: Paths ----
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/08.1b-TerraClimate_var_selection"
RAW_TSV="${OUT_DIR}/data/terraclimate_raw_monthly_19sites.tsv"
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.1b.2b-bioclim_variables.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b.2b-bioclim_variables.R"

YEAR_START=1970
YEAR_END=2000

mkdir -p "${LOG_DIR}"

# ---- Guard clause: verify each input before launching ----
for f in "${RAW_TSV}" "${R_SCRIPT}"; do
  [ -s "$f" ] || { echo "ERROR: required input missing or empty: $f"; exit 1; }
done

echo "START $(date -Iseconds) 08.1b.2b bioclim variables" | tee "${LOG_FILE}"

# ---- STEP 2: Run ----
Rscript --vanilla "${R_SCRIPT}" "${RAW_TSV}" "${OUT_DIR}" "${YEAR_START}" "${YEAR_END}" \
  >> "${LOG_FILE}" 2>&1

echo "DONE $(date -Iseconds) 08.1b.2b bioclim variables" | tee -a "${LOG_FILE}"
echo "Output: ${OUT_DIR}/data/bioclim_variables_19sites.csv"

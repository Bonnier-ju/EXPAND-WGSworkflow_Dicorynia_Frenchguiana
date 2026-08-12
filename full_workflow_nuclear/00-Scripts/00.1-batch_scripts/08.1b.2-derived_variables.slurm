#!/bin/bash
# 08.1b.2-derived_variables.slurm
#
# Compute the 10 custom TerraClimate-derived environmental variables
# (dry/wet season duration + intensity, precip_contrast, cwd_annual,
# tmax_dry_season, soil_moisture_min, vpd_max_annual) from the raw monthly
# data produced by 08.1b.1-terraclimate_monthly_extraction.slurm.
#
# See 08.1b.2-derived_variables.R header for the exact formulas and the
# 2026-08-11 decisions (dry season = aet>ppt, dry_season_intensity = deficit
# magnitude, vpd_max_annual = mean of annual maxima).
#
# Lightweight: ~19 sites x 31 years x 12 months, base R only. No heavy
# compute - small resource request.

#SBATCH -p workq
#SBATCH --time=00:30:00
#SBATCH -J 08.1b.2-derived_variables
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
LOG_FILE="${LOG_DIR}/08.1b.2-derived_variables.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b.2-derived_variables.R"

YEAR_START=1970
YEAR_END=2000

mkdir -p "${LOG_DIR}"

# ---- Guard clause: verify each input before launching ----
for f in "${RAW_TSV}" "${R_SCRIPT}"; do
  [ -s "$f" ] || { echo "ERROR: required input missing or empty: $f"; exit 1; }
done

echo "START $(date -Iseconds) 08.1b.2 derived variables" | tee "${LOG_FILE}"

# ---- STEP 2: Run ----
Rscript --vanilla "${R_SCRIPT}" "${RAW_TSV}" "${OUT_DIR}" "${YEAR_START}" "${YEAR_END}" \
  >> "${LOG_FILE}" 2>&1

echo "DONE $(date -Iseconds) 08.1b.2 derived variables" | tee -a "${LOG_FILE}"
echo "Output: ${OUT_DIR}/data/terraclimate_derived_variables_19sites.csv"

#!/bin/bash
# 08.1b.7-rda_forward_selection.slurm
#
# RDA forward selection + VIF screening on the Hellinger-transformed
# genotype matrix (154 ind x 433,798 SNPs) against the environmental
# variables Julien selected in 08.1b.5/08.1b.6.
#
# Heaviest step of the 08.1b pipeline - sized like the previous validated
# run of this same analysis (08.1c, July 2026: 137 ind x 433,798 SNPs,
# completed within these resources). Run `seff <jobid>` afterwards and
# adjust for future GEA methods (cf. CLAUDE.md §9).
#
# Do NOT submit before 08.1b.6 (PCA) has been reviewed and the same edited
# variable_selection_template.txt is final - this reads the same file.

#SBATCH -p workq
#SBATCH --time=12:00:00
#SBATCH -J 08.1b.7-rda_forward_selection
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ---- STEP 0: Environment ----
module purge
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# ---- STEP 1: Paths ----
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GEA_DIR="${BASE}/08-gea/workflow_full"
OUT_DIR="${GEA_DIR}/08.1-env_variable_selection/08.1b-TerraClimate_var_selection"

GENO_FILE="${GEA_DIR}/08.2-genotype_matrix/genotype_012_imputed.tsv"
META_FILE="${GEA_DIR}/08.2-genotype_matrix/sample_metadata.tsv"
ENV_MERGED="${OUT_DIR}/correlation_analysis/env_variables_merged_31vars.csv"
TEMPLATE="${OUT_DIR}/correlation_analysis/variable_selection_template.txt"
RDA_OUT="${OUT_DIR}/rda_vif"

LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.1b.7-rda_forward_selection.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b.7-rda_forward_selection.R"

N_PERM=999
SEED=$(head -n 1 "${BASE}/00-Scripts/00.1-batch_scripts/seed.txt")

mkdir -p "${LOG_DIR}" "${RDA_OUT}"

# ---- Guard clause: verify each input before launching ----
for f in "${GENO_FILE}" "${META_FILE}" "${ENV_MERGED}" "${TEMPLATE}" "${R_SCRIPT}"; do
  [ -s "$f" ] || { echo "ERROR: required input missing or empty: $f"; exit 1; }
done

echo "START $(date -Iseconds) 08.1b.7 RDA forward selection" | tee "${LOG_FILE}"
echo "  Seed: ${SEED} | Permutations: ${N_PERM}" | tee -a "${LOG_FILE}"

# ---- STEP 2: Run ----
Rscript --vanilla "${R_SCRIPT}" \
  "${GENO_FILE}" "${META_FILE}" "${ENV_MERGED}" "${TEMPLATE}" "${RDA_OUT}" "${N_PERM}" "${SEED}" \
  >> "${LOG_FILE}" 2>&1

echo "DONE $(date -Iseconds) 08.1b.7 RDA forward selection" | tee -a "${LOG_FILE}"
echo "Output: ${RDA_OUT}/variables_final.txt"
echo "Remember: seff \$SLURM_JOB_ID to check actual resource usage and adjust future GEA jobs."

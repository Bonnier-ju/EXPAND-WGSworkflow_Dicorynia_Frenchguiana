#!/bin/bash
# 08.1b.3-topo_extraction.slurm
#
# Extract elevation + TRI (Terrain Ruggedness Index) at the 19 study sites
# from the Copernicus DEM GLO-30 tiles already downloaded for the project.
# TWI dropped (confirmed 2026-08-11) - no flow accumulation tool available
# on Genotoul without a new install; simplified to 2 topographic variables.
#
# No download needed - tiles already present in data/elevation/copernicus_dem30/.

#SBATCH -p workq
#SBATCH --time=01:00:00
#SBATCH -J 08.1b.3-topo_extraction
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ---- STEP 0: Environment ----
module purge
module load statistics/R/4.3.0
# Combined: packages are split across two personal libraries on this account
# (~/R/R and the R-default ~/R/x86_64-pc-linux-gnu-library/4.3) - terra lives
# in the latter - see CLAUDE.md §5.
export R_LIBS_USER="${HOME}/R/R:${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# ---- STEP 1: Paths ----
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/08.1b-TerraClimate_var_selection"
DEM_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/elevation/copernicus_dem30"
GEOLOC="${BASE}/metadata/geoloc_site.csv"
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.1b.3-topo_extraction.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1b.3-topo_extraction.R"

mkdir -p "${LOG_DIR}"

# ---- Guard clause: verify each input before launching ----
for f in "${GEOLOC}" "${R_SCRIPT}"; do
  [ -s "$f" ] || { echo "ERROR: required input missing or empty: $f"; exit 1; }
done
if [ ! -d "${DEM_DIR}" ] || [ -z "$(ls -A "${DEM_DIR}"/*.tif 2>/dev/null)" ]; then
  echo "ERROR: no DEM tiles found in ${DEM_DIR}"; exit 1
fi

echo "START $(date -Iseconds) 08.1b.3 topo extraction" | tee "${LOG_FILE}"

# ---- STEP 2: Run ----
Rscript --vanilla "${R_SCRIPT}" "${DEM_DIR}" "${GEOLOC}" "${OUT_DIR}" \
  >> "${LOG_FILE}" 2>&1

echo "DONE $(date -Iseconds) 08.1b.3 topo extraction" | tee -a "${LOG_FILE}"
echo "Output: ${OUT_DIR}/data/topo_variables_19sites.csv"

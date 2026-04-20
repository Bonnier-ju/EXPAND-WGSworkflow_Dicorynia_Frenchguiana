#!/bin/bash
#SBATCH -p workq
#SBATCH --time=02:00:00
#SBATCH -J 08.1a-env_variable_correlation
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Define paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"

# CHELSA rasters (19 BIO variables, 30 arcsec, ~110 MB each)
CHELSA_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/chelsa/chelsa_rasters"

# ENVIREM rasters (18 variables, 30 arcsec)
ENVIREM_ROOT="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/envirem/elev_SAmerica_current_30arcsec_geotiff"

# Copernicus DEM GLO-30 tiles (elevation, ~30 m)
DEM_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/elevation/copernicus_dem30"

# Site-level CSVs (already in physical units) — used for scale-factor validation
CHELSA_SITE_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/chelsa/chelsa_env_per_site.csv"
ENVIREM_SITE_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/envirem/envirem_env_per_site.csv"
MANUAL_SITE_CSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/manual_variables/manual_variables_per_site.csv"

R_OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection"
LOG_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/logs"
LOG_FILE="${LOG_DIR}/08.1a-env_variable_correlation.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.1a-env_variable_correlation.R"

# Number of random landscape points for correlation estimation.
# 5000 points gives stable Pearson r estimates while keeping extraction fast.
N_POINTS=5000

# Correlation threshold for reporting high-correlation pairs.
# Pairs with |r| >= threshold are flagged for potential redundancy.
COR_THRESHOLD=0.7

mkdir -p "$LOG_DIR"
echo "START $(date -Iseconds) 08.1a landscape variable correlation" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: R — landscape sampling, extraction, correlation, outputs
# -------------------------------------------------------------------
echo "STEP 2: R landscape correlation analysis" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$CHELSA_DIR" \
    "$ENVIREM_ROOT" \
    "$DEM_DIR" \
    "$CHELSA_SITE_CSV" \
    "$ENVIREM_SITE_CSV" \
    "$MANUAL_SITE_CSV" \
    "$R_OUT_DIR" \
    "$N_POINTS" \
    "$COR_THRESHOLD" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.1a landscape correlation completed" | tee -a "$LOG_FILE"
echo "Outputs: ${R_OUT_DIR}/correlation_analysis/"

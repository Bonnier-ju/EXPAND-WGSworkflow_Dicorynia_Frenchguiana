#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J 08.0.1-chelsa_extraction
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
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
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GEOLOC="${BASE}/metadata/geoloc_site.csv"

OUT_DIR="${BASE}/08-gea/workflow_full/08.0-env_variables"
CHELSA_DIR="${OUT_DIR}/chelsa_rasters"
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/08.0.1-chelsa_extraction.log"

R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.0.1-chelsa_extraction.R"

mkdir -p "$OUT_DIR" "$CHELSA_DIR" "$LOG_DIR"
echo "START $(date -Iseconds) 08.0.1 CHELSA extraction" > "$LOG_FILE"

if [ ! -s "$GEOLOC" ]; then
    echo "ERROR: geoloc file not found: $GEOLOC" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Download CHELSA BIO1-BIO19 rasters (v2.1, 30 arcsec ~1 km)
# -------------------------------------------------------------------
# CHELSA v2.1 bioclimatic variables (1981-2010 climatology).
# Files are GeoTIFF (~200 MB each), downloaded once and reused.
# URL pattern: https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio{N}_1981-2010_V.2.1.tif
echo "STEP 2: downloading CHELSA BIO1-BIO19 rasters" | tee -a "$LOG_FILE"

CHELSA_BASE="https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio"

for N in $(seq 1 19); do
    fname="CHELSA_bio${N}_1981-2010_V.2.1.tif"
    fpath="${CHELSA_DIR}/${fname}"
    if [ -s "$fpath" ]; then
        echo "  BIO${N}: already present, skipping" | tee -a "$LOG_FILE"
    else
        echo "  BIO${N}: downloading..." | tee -a "$LOG_FILE"
        curl -fsSL --retry 3 --retry-delay 5 \
            "${CHELSA_BASE}/${fname}" \
            -o "$fpath" >> "$LOG_FILE" 2>&1
        echo "  BIO${N}: done ($(du -sh "$fpath" | cut -f1))" | tee -a "$LOG_FILE"
    fi
done

echo "STEP 2 completed: all BIO rasters present" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Extract CHELSA values at site coordinates (R)
# -------------------------------------------------------------------
echo "STEP 3: extracting CHELSA values at site coordinates" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$GEOLOC" \
    "$CHELSA_DIR" \
    "$OUT_DIR" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.0.1 CHELSA extraction completed" | tee -a "$LOG_FILE"
echo "Outputs: $OUT_DIR"

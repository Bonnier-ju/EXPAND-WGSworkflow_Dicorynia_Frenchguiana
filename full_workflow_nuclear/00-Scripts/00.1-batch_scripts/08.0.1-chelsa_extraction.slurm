#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J 08.0.1-chelsa_extraction
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
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
# gdallocationinfo is available system-wide in /usr/bin — no module needed

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

# Raw GDAL extraction output (tab-separated: site lat long BIO1 ... BIO19)
EXTRACTED_TSV="${OUT_DIR}/chelsa_raw_extracted.tsv"

mkdir -p "$OUT_DIR" "$CHELSA_DIR" "$LOG_DIR"
echo "START $(date -Iseconds) 08.0.1 CHELSA extraction" > "$LOG_FILE"

if [ ! -s "$GEOLOC" ]; then
    echo "ERROR: geoloc file not found: $GEOLOC" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Download CHELSA BIO1-BIO19 rasters (v2.1, 30 arcsec ~1 km)
# -------------------------------------------------------------------
# Files are ~110 MB each (GeoTIFF), downloaded once and reused.
CHELSA_BASE="https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio"
echo "STEP 2: downloading CHELSA BIO1-BIO19 rasters" | tee -a "$LOG_FILE"

for N in $(seq 1 19); do
    fname="CHELSA_bio${N}_1981-2010_V.2.1.tif"
    fpath="${CHELSA_DIR}/${fname}"
    if [ -s "$fpath" ]; then
        echo "  BIO${N}: already present, skipping" | tee -a "$LOG_FILE"
    else
        echo "  BIO${N}: downloading..." | tee -a "$LOG_FILE"
        curl -fsSL --retry 3 --retry-delay 5 \
            "${CHELSA_BASE}/${fname}" -o "$fpath" >> "$LOG_FILE" 2>&1
        echo "  BIO${N}: done ($(du -sh "$fpath" | cut -f1))" | tee -a "$LOG_FILE"
    fi
done
echo "STEP 2 completed" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Extract CHELSA values at site coordinates using gdallocationinfo
# -------------------------------------------------------------------
# gdallocationinfo -valonly -wgs84 <raster> <lon> <lat>
# Returns the pixel value at the given WGS84 coordinate (bilinear not
# available via CLI — nearest neighbour; for ~1km resolution and site-level
# extraction this is appropriate).
echo "STEP 3: extracting values with gdallocationinfo" | tee -a "$LOG_FILE"

# Build header line
HEADER="site\tlat\tlong"
for N in $(seq 1 19); do HEADER="${HEADER}\tBIO${N}"; done
echo -e "$HEADER" > "$EXTRACTED_TSV"

# Read geoloc CSV (skip header), extract for each site
tail -n +2 "$GEOLOC" | while IFS=',' read -r site lat lon; do
    # Remove carriage returns / trailing spaces
    site=$(echo "$site" | tr -d '\r' | xargs)
    lat=$(echo "$lat"   | tr -d '\r' | xargs)
    lon=$(echo "$lon"   | tr -d '\r' | xargs)

    row="${site}\t${lat}\t${lon}"
    for N in $(seq 1 19); do
        fname="${CHELSA_DIR}/CHELSA_bio${N}_1981-2010_V.2.1.tif"
        val=$(gdallocationinfo -valonly -wgs84 "$fname" "$lon" "$lat" 2>/dev/null \
              | tr -d '[:space:]')
        if [ -z "$val" ]; then val="NA"; fi
        row="${row}\t${val}"
    done
    echo -e "$row" >> "$EXTRACTED_TSV"
    echo "  Extracted: ${site}" | tee -a "$LOG_FILE"
done

echo "STEP 3 completed: $EXTRACTED_TSV" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: R post-processing — scale factors, standardisation, plots
# -------------------------------------------------------------------
echo "STEP 4: R post-processing" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$EXTRACTED_TSV" \
    "$OUT_DIR" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.0.1 CHELSA extraction completed" | tee -a "$LOG_FILE"
echo "Outputs: $OUT_DIR"

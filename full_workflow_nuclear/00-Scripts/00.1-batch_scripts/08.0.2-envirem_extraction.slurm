#!/bin/bash
#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 08.0.2-envirem_extraction
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
# gdallocationinfo is available system-wide in /usr/bin — no module needed

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GEOLOC="${BASE}/metadata/geoloc_site.csv"

# Sites to exclude from extraction:
#   - Outgroup sites outside French Guiana (Cameroun_Benin, Venezuela)
#   - Sites with NA coordinates (Herbier — herbarium specimens without precise geoloc)
# Note: St_georges is kept — it is a French Guiana site.
#   The 4 Angela/treemutation individuals (DBR_AGJ/K/L/M) assigned to St_georges
#   are excluded at the individual level in all downstream GEA analyses.
EXCLUDE_SITES="Cameroun_Benin"

OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/envirem"
RASTER_DIR="${OUT_DIR}/envirem_rasters"
LOG_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/logs"
LOG_FILE="${LOG_DIR}/08.0.2-envirem_extraction.log"

R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.0.2-envirem_extraction.R"

# Raw GDAL extraction output (tab-separated: site lat long + 36 monthly variables)
EXTRACTED_TSV="${OUT_DIR}/envirem_raw_extracted.tsv"

# BIO variables from CHELSA step — needed for embergerQ, thermicityIndex, continentality
CHELSA_PER_SITE="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/chelsa/chelsa_env_per_site.csv"

# R post-processing outputs go to the 08.1 root for downstream variable selection
R_OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection"

mkdir -p "$OUT_DIR" "$RASTER_DIR" "$LOG_DIR"
echo "START $(date -Iseconds) 08.0.2 ENVIREM extraction" > "$LOG_FILE"

if [ ! -s "$GEOLOC" ]; then
    echo "ERROR: geoloc file not found: $GEOLOC" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$CHELSA_PER_SITE" ]; then
    echo "ERROR: CHELSA per-site file not found (run 08.0.1 first): $CHELSA_PER_SITE" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Download CHELSA v2.1 monthly rasters (1981-2010 climatologies)
# -------------------------------------------------------------------
# Three variable types × 12 months = 36 rasters (~4 GB total).
# Files are downloaded once and reused on subsequent runs.
#
# CHELSA v2.1 monthly encoding (confirmed via gdalinfo):
#   tas (monthly mean air temp)  : UInt16, Scale=0.1, Offset=-273.15  → °C = raw×0.1 - 273.15
#   pr  (monthly precipitation)  : UInt16, Scale=0.1, Offset=0        → mm/month = raw×0.1
#   pet_penman_monthly (PET Penman-Monteith) : UInt16, Scale=0.1, Offset=0 → mm/month = raw×0.1
#
# NOTE: verify pet encoding on first download with:
#   gdalinfo ${RASTER_DIR}/CHELSA_pet_penman_monthly_01_1981-2010_V.2.1.tif | grep -E 'Scale|Offset|Type'
#
CHELSA_BASE="https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010"
echo "STEP 2: downloading CHELSA monthly rasters (tas, pr, pet_penman)" | tee -a "$LOG_FILE"

# Monthly mean air temperature (tas)
for MM in $(seq -w 1 12); do
    fname="CHELSA_tas_${MM}_1981-2010_V.2.1.tif"
    fpath="${RASTER_DIR}/${fname}"
    if [ -s "$fpath" ]; then
        echo "  tas_${MM}: already present, skipping" | tee -a "$LOG_FILE"
    else
        echo "  tas_${MM}: downloading..." | tee -a "$LOG_FILE"
        curl -fsSL --retry 3 --retry-delay 5 \
            "${CHELSA_BASE}/tas/${fname}" -o "$fpath" >> "$LOG_FILE" 2>&1
        echo "  tas_${MM}: done ($(du -sh "$fpath" | cut -f1))" | tee -a "$LOG_FILE"
    fi
done

# Monthly precipitation (pr)
for MM in $(seq -w 1 12); do
    fname="CHELSA_pr_${MM}_1981-2010_V.2.1.tif"
    fpath="${RASTER_DIR}/${fname}"
    if [ -s "$fpath" ]; then
        echo "  pr_${MM}: already present, skipping" | tee -a "$LOG_FILE"
    else
        echo "  pr_${MM}: downloading..." | tee -a "$LOG_FILE"
        curl -fsSL --retry 3 --retry-delay 5 \
            "${CHELSA_BASE}/pr/${fname}" -o "$fpath" >> "$LOG_FILE" 2>&1
        echo "  pr_${MM}: done ($(du -sh "$fpath" | cut -f1))" | tee -a "$LOG_FILE"
    fi
done

# Monthly Penman-Monteith PET (pet_penman_monthly)
for MM in $(seq -w 1 12); do
    fname="CHELSA_pet_penman_monthly_${MM}_1981-2010_V.2.1.tif"
    fpath="${RASTER_DIR}/${fname}"
    if [ -s "$fpath" ]; then
        echo "  pet_${MM}: already present, skipping" | tee -a "$LOG_FILE"
    else
        echo "  pet_${MM}: downloading..." | tee -a "$LOG_FILE"
        curl -fsSL --retry 3 --retry-delay 5 \
            "${CHELSA_BASE}/pet_penman_monthly/${fname}" -o "$fpath" >> "$LOG_FILE" 2>&1
        echo "  pet_${MM}: done ($(du -sh "$fpath" | cut -f1))" | tee -a "$LOG_FILE"
    fi
done

echo "STEP 2 completed (36 rasters)" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2b: Print encoding info for first raster of each type
# -------------------------------------------------------------------
# This block verifies scale/offset so the R script uses correct factors.
echo "STEP 2b: verifying raster encoding (gdalinfo)" | tee -a "$LOG_FILE"
for var in tas pr pet_penman_monthly; do
    case "$var" in
        tas)                fname="${RASTER_DIR}/CHELSA_tas_01_1981-2010_V.2.1.tif" ;;
        pr)                 fname="${RASTER_DIR}/CHELSA_pr_01_1981-2010_V.2.1.tif" ;;
        pet_penman_monthly) fname="${RASTER_DIR}/CHELSA_pet_penman_monthly_01_1981-2010_V.2.1.tif" ;;
    esac
    echo "  ${var}:" | tee -a "$LOG_FILE"
    gdalinfo "$fname" 2>/dev/null \
        | grep -E 'Type=|Scale=|Offset=|NoData' \
        | head -6 | tee -a "$LOG_FILE"
done

# -------------------------------------------------------------------
# STEP 3: Extract monthly values at site coordinates using gdallocationinfo
# -------------------------------------------------------------------
# Output columns:
#   site  lat  long  tas_01..tas_12  pr_01..pr_12  pet_01..pet_12
echo "STEP 3: extracting monthly values with gdallocationinfo" | tee -a "$LOG_FILE"

# Build header
HEADER="site\tlat\tlong"
for MM in $(seq -w 1 12); do HEADER="${HEADER}\ttas_${MM}"; done
for MM in $(seq -w 1 12); do HEADER="${HEADER}\tpr_${MM}";  done
for MM in $(seq -w 1 12); do HEADER="${HEADER}\tpet_${MM}"; done
echo -e "$HEADER" > "$EXTRACTED_TSV"

tail -n +2 "$GEOLOC" | while IFS=',' read -r site lat lon; do
    site=$(echo "$site" | tr -d '\r' | xargs)
    lat=$(echo  "$lat"  | tr -d '\r' | xargs)
    lon=$(echo  "$lon"  | tr -d '\r' | xargs)

    # Skip outgroup sites and sites with missing coordinates
    if echo "$EXCLUDE_SITES" | grep -qw "$site"; then
        echo "  Skipping: ${site} (outgroup site)" | tee -a "$LOG_FILE"
        continue
    fi
    if [ "$lat" = "NA" ] || [ "$lon" = "NA" ] || [ -z "$lat" ] || [ -z "$lon" ]; then
        echo "  Skipping: ${site} (missing coordinates)" | tee -a "$LOG_FILE"
        continue
    fi

    row="${site}\t${lat}\t${lon}"

    # Monthly mean air temperature
    for MM in $(seq -w 1 12); do
        fname="${RASTER_DIR}/CHELSA_tas_${MM}_1981-2010_V.2.1.tif"
        val=$(gdallocationinfo -valonly -wgs84 "$fname" "$lon" "$lat" 2>/dev/null \
              | tr -d '[:space:]')
        [ -z "$val" ] && val="NA"
        row="${row}\t${val}"
    done

    # Monthly precipitation
    for MM in $(seq -w 1 12); do
        fname="${RASTER_DIR}/CHELSA_pr_${MM}_1981-2010_V.2.1.tif"
        val=$(gdallocationinfo -valonly -wgs84 "$fname" "$lon" "$lat" 2>/dev/null \
              | tr -d '[:space:]')
        [ -z "$val" ] && val="NA"
        row="${row}\t${val}"
    done

    # Monthly PET (Penman-Monteith)
    for MM in $(seq -w 1 12); do
        fname="${RASTER_DIR}/CHELSA_pet_penman_monthly_${MM}_1981-2010_V.2.1.tif"
        val=$(gdallocationinfo -valonly -wgs84 "$fname" "$lon" "$lat" 2>/dev/null \
              | tr -d '[:space:]')
        [ -z "$val" ] && val="NA"
        row="${row}\t${val}"
    done

    echo -e "$row" >> "$EXTRACTED_TSV"
    echo "  Extracted: ${site}" | tee -a "$LOG_FILE"
done

echo "STEP 3 completed: $EXTRACTED_TSV" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: R post-processing — compute ENVIREM variables, scale, plot
# -------------------------------------------------------------------
echo "STEP 4: R post-processing" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$EXTRACTED_TSV" \
    "$CHELSA_PER_SITE" \
    "$R_OUT_DIR" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.0.2 ENVIREM extraction completed" | tee -a "$LOG_FILE"
echo "Outputs: $OUT_DIR"

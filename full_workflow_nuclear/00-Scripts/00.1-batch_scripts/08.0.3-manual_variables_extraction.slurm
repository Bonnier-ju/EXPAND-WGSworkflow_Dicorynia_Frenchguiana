#!/bin/bash
#SBATCH -p workq
#SBATCH --time=01:00:00
#SBATCH -J 08.0.3-manual_variables_extraction
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=4G
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
# gdallocationinfo available system-wide

# -------------------------------------------------------------------
# STEP 1: Define paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GEOLOC="${BASE}/metadata/geoloc_site.csv"

MANUAL_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/manual_variables"
HABITAT_RASTER="${MANUAL_DIR}/Habitats forestiers.tif"
PEDOL_SHP="${MANUAL_DIR}/Pedologie/pedologie_guyane.shp"

LOG_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/logs"
LOG_FILE="${LOG_DIR}/08.0.3-manual_variables_extraction.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.0.3-manual_variables_extraction.R"

RAW_TSV="${MANUAL_DIR}/forest_habitat_raw.tsv"
R_OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection"

# Sites to exclude (outgroup and NA coordinates)
EXCLUDE_SITES="Cameroun_Benin"

# Minimum number of sites a category must appear in to generate a dummy variable.
# Categories present at fewer sites provide insufficient power for GEA.
MIN_SITES_DUMMY=2

mkdir -p "$LOG_DIR"
echo "START $(date -Iseconds) 08.0.3 manual variables extraction" > "$LOG_FILE"

if [ ! -s "$GEOLOC" ]; then
    echo "ERROR: geoloc file not found: $GEOLOC" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$HABITAT_RASTER" ]; then
    echo "ERROR: forest habitat raster not found: $HABITAT_RASTER" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$PEDOL_SHP" ]; then
    echo "ERROR: pedology shapefile not found: $PEDOL_SHP" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Extract forest habitat code at site coordinates
# -------------------------------------------------------------------
# "Habitats forestiers.tif": Float32 raster, nominal categorical codes
# (e.g. 41220, 41420, 41610). Values define forest habitat types based
# on pedology and vegetation. No scale factor — values are integer codes.
echo "STEP 2: extracting forest habitat codes" | tee -a "$LOG_FILE"

echo -e "site\tlat\tlong\tforest_habitat_code" > "$RAW_TSV"

tail -n +2 "$GEOLOC" | while IFS=',' read -r site lat lon; do
    site=$(echo "$site" | tr -d '\r' | xargs)
    lat=$(echo  "$lat"  | tr -d '\r' | xargs)
    lon=$(echo  "$lon"  | tr -d '\r' | xargs)

    if echo "$EXCLUDE_SITES" | grep -qw "$site"; then
        echo "  Skipping: ${site} (outgroup site)" | tee -a "$LOG_FILE"
        continue
    fi
    if [ "$lat" = "NA" ] || [ "$lon" = "NA" ] || [ -z "$lat" ] || [ -z "$lon" ]; then
        echo "  Skipping: ${site} (missing coordinates)" | tee -a "$LOG_FILE"
        continue
    fi

    val=$(gdallocationinfo -valonly -wgs84 "$HABITAT_RASTER" "$lon" "$lat" 2>/dev/null \
          | tr -d '[:space:]')
    [ -z "$val" ] && val="NA"

    echo -e "${site}\t${lat}\t${lon}\t${val}" >> "$RAW_TSV"
    echo "  ${site}: forest_habitat = ${val}" | tee -a "$LOG_FILE"
done

echo "STEP 2 completed: $RAW_TSV" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: R — pedology spatial join + dummy coding for both variables
# -------------------------------------------------------------------
echo "STEP 3: R spatial join and dummy coding" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$RAW_TSV" \
    "$PEDOL_SHP" \
    "$R_OUT_DIR" \
    "$MIN_SITES_DUMMY" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.0.3 manual variables extraction completed" | tee -a "$LOG_FILE"
echo "Outputs: ${R_OUT_DIR}/data/manual_variables/"

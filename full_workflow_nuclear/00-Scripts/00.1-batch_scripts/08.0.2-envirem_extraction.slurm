#!/bin/bash
#SBATCH -p workq
#SBATCH --time=03:00:00
#SBATCH -J 08.0.2-envirem_extraction
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
# gdallocationinfo and python3 are available system-wide

# -------------------------------------------------------------------
# STEP 1: Define paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GEOLOC="${BASE}/metadata/geoloc_site.csv"

# ENVIREM rasters (pre-downloaded, 18 variables from WorldClim/ENVIREM)
ENVIREM_ROOT="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/envirem/elev_SAmerica_current_30arcsec_geotiff"
SET1="${ENVIREM_ROOT}/SAmerica_current_30arcsec_geotiff_set1"
SET2="${ENVIREM_ROOT}/SAmerica_current_30arcsec_geotiff_set2"
SET3="${ENVIREM_ROOT}/SAmerica_current_30arcsec_geotiff_set3"
SET4="${ENVIREM_ROOT}/SAmerica_current_30arcsec_geotiff_set4"
ELEV_SET="${ENVIREM_ROOT}/elev_SAmerica_current_30arcsec_geotiff"

# Copernicus DEM GLO-30 tiles (downloaded here, 30 m resolution, Float32, metres)
DEM_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/elevation/copernicus_dem30"

LOG_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/logs"
LOG_FILE="${LOG_DIR}/08.0.2-envirem_extraction.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.0.2-envirem_extraction.R"
EXTRACTED_TSV="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/envirem/envirem_raw_extracted.tsv"
R_OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection"

# Sites to exclude (outgroup and NA coordinates — see 08.0.1 for rationale)
EXCLUDE_SITES="Cameroun_Benin"

mkdir -p "$DEM_DIR" "$LOG_DIR"
echo "START $(date -Iseconds) 08.0.2 ENVIREM + elevation extraction" > "$LOG_FILE"

if [ ! -s "$GEOLOC" ]; then
    echo "ERROR: geoloc file not found: $GEOLOC" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Verify ENVIREM rasters
# -------------------------------------------------------------------
# Encoding summary:
#   Float32, physical units, no scaling: annualPET (mm/yr), aridityIndexThornthwaite,
#     climaticMoistureIndex, continentality, embergerQ, growingDegDays0 (°C·days),
#     growingDegDays5, monthCountByTemp10, PETColdestQuarter (mm),
#     PETDriestQuarter (mm), PETseasonality (SD×100), PETWarmestQuarter (mm),
#     PETWettestQuarter (mm), thermicityIndex, topoWet, tri
#   Int16, ÷10 in R to get °C: maxTempColdest, minTempWarmest
echo "STEP 2: verifying ENVIREM rasters" | tee -a "$LOG_FILE"

declare -A RASTERS
RASTERS["annualPET"]="${SET1}/current_30arcsec_annualPET.tif"
RASTERS["aridityIndexThornthwaite"]="${SET1}/current_30arcsec_aridityIndexThornthwaite.tif"
RASTERS["climaticMoistureIndex"]="${SET1}/current_30arcsec_climaticMoistureIndex.tif"
RASTERS["continentality"]="${SET1}/current_30arcsec_continentality.tif"
RASTERS["embergerQ"]="${SET2}/current_30arcsec_embergerQ.tif"
RASTERS["growingDegDays0"]="${SET2}/current_30arcsec_growingDegDays0.tif"
RASTERS["growingDegDays5"]="${SET2}/current_30arcsec_growingDegDays5.tif"
RASTERS["maxTempColdest"]="${SET2}/current_30arcsec_maxTempColdest.tif"
RASTERS["minTempWarmest"]="${SET3}/current_30arcsec_minTempWarmest.tif"
RASTERS["monthCountByTemp10"]="${SET3}/current_30arcsec_monthCountByTemp10.tif"
RASTERS["PETColdestQuarter"]="${SET3}/current_30arcsec_PETColdestQuarter.tif"
RASTERS["PETDriestQuarter"]="${SET3}/current_30arcsec_PETDriestQuarter.tif"
RASTERS["PETseasonality"]="${SET4}/current_30arcsec_PETseasonality.tif"
RASTERS["PETWarmestQuarter"]="${SET4}/current_30arcsec_PETWarmestQuarter.tif"
RASTERS["PETWettestQuarter"]="${SET4}/current_30arcsec_PETWettestQuarter.tif"
RASTERS["thermicityIndex"]="${SET4}/current_30arcsec_thermicityIndex.tif"
RASTERS["topoWet"]="${ELEV_SET}/current_30arcsec_topoWet.tif"
RASTERS["tri"]="${ELEV_SET}/current_30arcsec_tri.tif"

ENVIREM_VARS="annualPET aridityIndexThornthwaite climaticMoistureIndex continentality \
              embergerQ growingDegDays0 growingDegDays5 maxTempColdest minTempWarmest \
              monthCountByTemp10 PETColdestQuarter PETDriestQuarter PETseasonality \
              PETWarmestQuarter PETWettestQuarter thermicityIndex topoWet tri"

MISSING=0
for var in $ENVIREM_VARS; do
    f="${RASTERS[$var]}"
    if [ ! -s "$f" ]; then
        echo "  ERROR: missing or empty raster: $f" | tee -a "$LOG_FILE"
        MISSING=$((MISSING + 1))
    fi
done
if [ "$MISSING" -gt 0 ]; then
    echo "ERROR: $MISSING ENVIREM raster(s) missing — run downloads first" | tee -a "$LOG_FILE"
    exit 1
fi
echo "  18/18 ENVIREM rasters verified" | tee -a "$LOG_FILE"
echo "STEP 2 completed" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Download Copernicus DEM GLO-30 tiles for French Guiana sites
# -------------------------------------------------------------------
# Resolution: ~30 m (1 arcsec), Float32, values in metres (EGM2008 geoid).
# Tiles needed for French Guiana sites (lat 3–6°N, lon 52–55°W): 7 tiles.
# Source: AWS Open Data — https://registry.opendata.aws/copernicus-dem/
# URL: https://copernicus-dem-30m.s3.amazonaws.com/<TILENAME>/<TILENAME>.tif
# Tile naming: Copernicus_DSM_COG_10_N{lat:02d}_00_W{lon:03d}_00_DEM
echo "STEP 3: downloading Copernicus DEM GLO-30 tiles" | tee -a "$LOG_FILE"

COPERNICUS_BASE="https://copernicus-dem-30m.s3.amazonaws.com"

for TILE in N03_W053 N03_W054 N04_W052 N04_W053 N04_W054 N05_W054 N05_W055; do
    TILENAME="Copernicus_DSM_COG_10_${TILE}_00_DEM"
    TILEPATH="${DEM_DIR}/${TILENAME}.tif"
    if [ -s "$TILEPATH" ]; then
        echo "  ${TILE}: already present, skipping" | tee -a "$LOG_FILE"
    else
        echo "  ${TILE}: downloading..." | tee -a "$LOG_FILE"
        curl -fsSL --retry 3 --retry-delay 5 \
            "${COPERNICUS_BASE}/${TILENAME}/${TILENAME}.tif" \
            -o "$TILEPATH" >> "$LOG_FILE" 2>&1
        echo "  ${TILE}: done ($(du -sh "$TILEPATH" | cut -f1))" | tee -a "$LOG_FILE"
    fi
done
echo "STEP 3 completed" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: Extract ENVIREM + elevation values at site coordinates
# -------------------------------------------------------------------
# Output columns: site  lat  long  <18 ENVIREM vars>  elevation
echo "STEP 4: extracting values at site coordinates" | tee -a "$LOG_FILE"

# Build TSV header
HEADER="site\tlat\tlong"
for var in $ENVIREM_VARS; do HEADER="${HEADER}\t${var}"; done
HEADER="${HEADER}\televation"
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

    # Extract 18 ENVIREM variables
    for var in $ENVIREM_VARS; do
        val=$(gdallocationinfo -valonly -wgs84 "${RASTERS[$var]}" "$lon" "$lat" 2>/dev/null \
              | tr -d '[:space:]')
        [ -z "$val" ] && val="NA"
        row="${row}\t${val}"
    done

    # Determine Copernicus DEM tile for this site, then extract elevation
    # Tile lat = floor(lat), tile lon = abs(floor(lon)) for Western hemisphere
    TILE_LAT=$(python3 -c "import math; print(int(math.floor($lat)))")
    TILE_LON=$(python3 -c "import math; print(int(abs(math.floor($lon))))")
    TILENAME=$(printf "Copernicus_DSM_COG_10_N%02d_00_W%03d_00_DEM" "$TILE_LAT" "$TILE_LON")
    TILEPATH="${DEM_DIR}/${TILENAME}.tif"

    if [ -s "$TILEPATH" ]; then
        elev=$(gdallocationinfo -valonly -wgs84 "$TILEPATH" "$lon" "$lat" 2>/dev/null \
               | tr -d '[:space:]')
        [ -z "$elev" ] && elev="NA"
    else
        echo "  WARN: DEM tile not found for ${site}: ${TILENAME}" | tee -a "$LOG_FILE"
        elev="NA"
    fi
    row="${row}\t${elev}"

    echo -e "$row" >> "$EXTRACTED_TSV"
    echo "  Extracted: ${site} (elev=${elev} m)" | tee -a "$LOG_FILE"
done

echo "STEP 4 completed: $EXTRACTED_TSV" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 5: R post-processing — scale factors, standardisation, plots
# -------------------------------------------------------------------
echo "STEP 5: R post-processing" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$EXTRACTED_TSV" \
    "$R_OUT_DIR" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.0.2 ENVIREM + elevation extraction completed" | tee -a "$LOG_FILE"
echo "Outputs: ${R_OUT_DIR}/data/envirem/"

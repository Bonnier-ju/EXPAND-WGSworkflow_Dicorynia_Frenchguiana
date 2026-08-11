#!/bin/bash
# 08.1b.1-terraclimate_monthly_extraction.slurm
#
# Download TerraClimate monthly NetCDF files (1970-2000) and extract
# values at the 19 French Guiana study sites.
#
# Refonte du plan 08.1 (cf. CLAUDE.md §4) : source unique TerraClimate,
# fenetre 1970-2000, methode d'extraction adaptee de 08.0.4-terraclimate_extraction.slurm
# (meme mecanique curl + gdallocationinfo, deja validee sur Genotoul) mais avec
# de nouvelles variables (tmax, tmin, pet) et une nouvelle fenetre temporelle.
#
# Variables telechargees (8) - celles requises par le plan 08.1b (etape 1) :
#   aet  : Actual Evapotranspiration (mm/month)
#   def  : Climatic Water Deficit    (mm/month)
#   pet  : Potential Evapotranspiration (mm/month)
#   ppt  : Precipitation (mm/month)
#   soil : Soil moisture (mm)
#   tmax : Maximum temperature (degC)
#   tmin : Minimum temperature (degC)
#   vpd  : Vapor pressure deficit (kPa)
#
# All variables stored as Int16 x scale_factor in NetCDF; GDAL applies
# scale_factor/add_offset automatically -> gdallocationinfo returns physical units.
#
# Source: https://www.climatologylab.org/terraclimate.html
#   Files: https://climate.northwestknowledge.net/TERRACLIMATE-DATA/
#   Format: NetCDF4, 1/24 deg (~4 km), monthly, global
#
# Strategy: download one file at a time, extract immediately, delete -> max ~300 MB
# disk usage at any one point regardless of file size.
#
# This script only extracts raw monthly values. Annual aggregation is done
# downstream by 08.1b.2 (10 custom derived variables) and 08.1b.2b (BIO1-19
# via dismo::biovars() on the monthly climatology) - kept separate since the
# two consume the raw monthly data differently (per-year vs climatological).
#
# Output:
#   08.1b-TerraClimate_var_selection/data/terraclimate_raw_monthly_19sites.tsv
#     (site lat long variable year m01 m02 ... m12)

#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 08.1b.1-terraclimate_extraction
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment
# -------------------------------------------------------------------
module purge
# gdallocationinfo available system-wide on Genotoul (no module load required
# - confirmed working via 08.0.4-terraclimate_extraction.slurm, same mechanism).

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
YEAR_START=1970
YEAR_END=2000

# TerraClimate variable names (case-sensitive - match the filename and the
# NetCDF internal variable name; all lowercase for this set, unlike PDSI).
VARS="aet def pet ppt soil tmax tmin vpd"

# Sites to exclude (outgroup - no French Guiana coordinates). geoloc_site.csv
# currently holds only the 19 Guyane sites, but this guard is kept for safety
# in case the file is ever extended.
EXCLUDE_SITES="Cameroun_Benin"

TC_BASE="https://climate.northwestknowledge.net/TERRACLIMATE-DATA"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GEOLOC="${BASE}/metadata/geoloc_site.csv"
OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/08.1b-TerraClimate_var_selection"
DATA_DIR="${OUT_DIR}/data"
LOG_DIR="${OUT_DIR}/logs"
LOG="${LOG_DIR}/08.1b.1-terraclimate_monthly_extraction.log"

RAW_TSV="${DATA_DIR}/terraclimate_raw_monthly_19sites.tsv"
TMP_NC="${DATA_DIR}/.tmp_terraclimate.nc"   # single temp file, reused per download

mkdir -p "${DATA_DIR}" "${LOG_DIR}"
echo "START $(date -Iseconds) 08.1b.1 TerraClimate monthly extraction" | tee "${LOG}"
echo "  Period: ${YEAR_START}-${YEAR_END} | Variables: ${VARS}" | tee -a "${LOG}"

# ---- Guard clause: verify required input before launching ----
if [ ! -s "${GEOLOC}" ]; then
    echo "ERROR: required input missing or empty: ${GEOLOC}" | tee -a "${LOG}"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Read site coordinates into arrays (skip header + excluded sites)
# -------------------------------------------------------------------
echo "STEP 2: loading site coordinates ..." | tee -a "${LOG}"

declare -a SITE_NAMES SITE_LATS SITE_LONS

while IFS=',' read -r site lat lon; do
    site=$(echo "$site" | tr -d '\r' | xargs)
    lat=$(echo  "$lat"  | tr -d '\r' | xargs)
    lon=$(echo  "$lon"  | tr -d '\r' | xargs)
    if echo "${EXCLUDE_SITES}" | grep -qw "$site"; then continue; fi
    if [ "$lat" = "NA" ] || [ -z "$lat" ]; then continue; fi
    SITE_NAMES+=("$site")
    SITE_LATS+=("$lat")
    SITE_LONS+=("$lon")
done < <(tail -n +2 "${GEOLOC}")

N_SITES=${#SITE_NAMES[@]}
echo "  ${N_SITES} sites loaded" | tee -a "${LOG}"

if [ "${N_SITES}" -ne 19 ]; then
    echo "WARN: expected 19 sites, got ${N_SITES} - check ${GEOLOC}" | tee -a "${LOG}"
fi

# -------------------------------------------------------------------
# STEP 3: Write raw TSV header
# Columns: site  lat  long  variable  year  m01 m02 ... m12
# One row per site x variable x year.
# -------------------------------------------------------------------
printf "site\tlat\tlong\tvariable\tyear\tm01\tm02\tm03\tm04\tm05\tm06\tm07\tm08\tm09\tm10\tm11\tm12\n" \
    > "${RAW_TSV}"

# -------------------------------------------------------------------
# STEP 4: Download, extract, delete - one file at a time
# For each variable x year: download NetCDF -> gdallocationinfo -> delete.
# gdallocationinfo without -b returns ALL bands -> 12 lines (months).
# paste -s converts 12 lines to one tab-separated row.
# -------------------------------------------------------------------
echo "STEP 4: downloading and extracting monthly values ..." | tee -a "${LOG}"

N_VARS=$(echo "${VARS}" | wc -w)
var_idx=0
for VAR in ${VARS}; do
    var_idx=$((var_idx + 1))
    echo "  [${var_idx}/${N_VARS}] ${VAR}" | tee -a "${LOG}"
    n_ok=0
    n_fail=0

    for YEAR in $(seq "${YEAR_START}" "${YEAR_END}"); do
        NC_URL="${TC_BASE}/TerraClimate_${VAR}_${YEAR}.nc"
        SUBDS="NETCDF:${TMP_NC}:${VAR}"

        # Download with retry (overwrite temp file each time)
        if ! curl -fsSL --retry 5 --retry-delay 10 \
                --connect-timeout 30 --max-time 1800 \
                "${NC_URL}" -o "${TMP_NC}" 2>> "${LOG}"; then
            echo "    WARN: download failed for ${VAR} ${YEAR}" | tee -a "${LOG}"
            n_fail=$((n_fail + 1))
            continue
        fi

        # Extract all 12 monthly bands at each site
        for i in $(seq 0 $((N_SITES - 1))); do
            SNAME="${SITE_NAMES[$i]}"
            SLAT="${SITE_LATS[$i]}"
            SLON="${SITE_LONS[$i]}"

            # gdallocationinfo returns one value per line (band = month)
            monthly=$(gdallocationinfo -valonly -wgs84 \
                          "${SUBDS}" "${SLON}" "${SLAT}" 2>/dev/null \
                      | paste -s -d '\t')

            if [ -z "${monthly}" ]; then
                monthly=$(printf 'NA\t%.0s' {1..12} | sed 's/\t$//')
            fi

            printf "%s\t%s\t%s\t%s\t%d\t%s\n" \
                "${SNAME}" "${SLAT}" "${SLON}" "${VAR}" "${YEAR}" "${monthly}" \
                >> "${RAW_TSV}"
        done

        n_ok=$((n_ok + 1))
    done

    # Delete temp NetCDF immediately to free disk
    rm -f "${TMP_NC}"
    echo "    done: ${n_ok} years OK, ${n_fail} failed" | tee -a "${LOG}"
done

N_ROWS=$(wc -l < "${RAW_TSV}")
echo "STEP 4 completed: ${N_ROWS} data rows in ${RAW_TSV}" | tee -a "${LOG}"

if [ "${N_ROWS}" -lt 2 ]; then
    echo "ERROR: raw TSV is empty - extraction failed" | tee -a "${LOG}"; exit 1
fi

EXPECTED_ROWS=$(( N_SITES * N_VARS * (YEAR_END - YEAR_START + 1) + 1 ))
echo "  Expected ~${EXPECTED_ROWS} rows (19 sites x 8 vars x 31 years + header)" | tee -a "${LOG}"

echo "DONE $(date -Iseconds) 08.1b.1 TerraClimate monthly extraction" | tee -a "${LOG}"
echo "Output: ${RAW_TSV}"
echo "Next: 08.1b.2-derived_variables.slurm and 08.1b.2b-bioclim_variables.slurm"

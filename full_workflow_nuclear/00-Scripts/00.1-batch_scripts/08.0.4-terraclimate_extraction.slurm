#!/bin/bash
# 08.0.4-terraclimate_extraction.slurm
#
# Download TerraClimate monthly NetCDF files (1991-2020) and extract
# values at French Guiana study sites.
#
# Variables downloaded (7) and derived outputs (11).
# Retained variables are unique to TerraClimate — not covered by CHELSA BIO1-19
# or ENVIREM. Dropped: ppt (≈BIO12), tmax/tmin (≈BIO1/5/6), pet (≈annualPET),
# vap/q/ws (low ecological relevance or covered by vpd/soil).
#
# Downloaded → derived output columns:
#   aet   : Actual Evapotranspiration (mm/month) → annual sum → aet_mean + aet_sd
#   def   : Climatic Water Deficit   (mm/month) → annual sum → def_mean + def_sd
#   PDSI  : Palmer Drought Severity Index        → annual mean → PDSI_mean + PDSI_sd
#   ppt   : Precipitation (mm/month) → annual sum → ppt_sd only (mean ≈ BIO12)
#   soil  : Soil moisture (mm)       → annual mean → soil_mean + soil_sd
#   srad  : Shortwave radiation (W/m²/day)       → annual mean → srad_mean only
#   vpd   : Vapor pressure deficit (kPa)         → annual mean → vpd_mean only
#
# All variables stored as Int16 × scale_factor in NetCDF; GDAL applies
# scale_factor/add_offset automatically → gdallocationinfo returns physical units.
#
# Source: https://www.climatologylab.org/terraclimate.html
#   Files: https://climate.northwestknowledge.net/TERRACLIMATE-DATA/
#   Format: NetCDF4, 1/24° (~4 km), monthly, global
#
# Strategy: download one file at a time, extract immediately, delete → max ~300 MB
# disk usage at any one point regardless of file size.
#
# Output:
#   data/terraclimate/terraclimate_raw_monthly.tsv
#     (site lat long variable year m01 m02 ... m12)
#   Then calls 08.0.4-terraclimate_extraction.R for:
#     data/terraclimate/terraclimate_env_per_site.csv
#     data/terraclimate/terraclimate_env_scaled.csv
#     plots/terraclimate_*

#SBATCH -p workq
#SBATCH --time=06:00:00
#SBATCH -J 08.0.4-terraclimate
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
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"
# gdallocationinfo available system-wide

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
YEAR_START=1991
YEAR_END=2020

# TerraClimate variable names (case-sensitive — match the filename and the
# NetCDF internal variable name; PDSI uses uppercase).
VARS="aet def PDSI ppt soil srad vpd"

# Sites to exclude (outgroup — no French Guiana coordinates)
EXCLUDE_SITES="Cameroun_Benin"

TC_BASE="https://climate.northwestknowledge.net/TERRACLIMATE-DATA"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GEOLOC="${BASE}/metadata/geoloc_site.csv"
OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection"
DATA_DIR="${OUT_DIR}/data/terraclimate"
LOG_DIR="${OUT_DIR}/logs"
LOG="${LOG_DIR}/08.0.4-terraclimate_extraction.log"

RAW_TSV="${DATA_DIR}/terraclimate_raw_monthly.tsv"
TMP_NC="${DATA_DIR}/.tmp_terraclimate.nc"   # single temp file, reused per download

R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.0.4-terraclimate_extraction.R"

mkdir -p "${DATA_DIR}" "${LOG_DIR}" "${OUT_DIR}/plots"
echo "START $(date -Iseconds) 08.0.4 TerraClimate extraction" | tee "${LOG}"
echo "  Period: ${YEAR_START}–${YEAR_END} | Variables: ${VARS}" | tee -a "${LOG}"

if [ ! -s "${GEOLOC}" ]; then
    echo "ERROR: geoloc file not found: ${GEOLOC}" | tee -a "${LOG}"; exit 1
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

# -------------------------------------------------------------------
# STEP 3: Write raw TSV header
# Columns: site  lat  long  variable  year  m01 m02 ... m12
# One row per site × variable × year.
# -------------------------------------------------------------------
printf "site\tlat\tlong\tvariable\tyear\tm01\tm02\tm03\tm04\tm05\tm06\tm07\tm08\tm09\tm10\tm11\tm12\n" \
    > "${RAW_TSV}"

# -------------------------------------------------------------------
# STEP 4: Download, extract, delete — one file at a time
# For each variable × year: download NetCDF → gdallocationinfo → delete.
# gdallocationinfo without -b returns ALL bands → 12 lines (months).
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
    echo "ERROR: raw TSV is empty — extraction failed" | tee -a "${LOG}"; exit 1
fi

# -------------------------------------------------------------------
# STEP 5: R post-processing
# Annual aggregation (sum/mean), 30-year climatological statistics,
# z-score standardisation, correlation plots, per-site variable plots.
# -------------------------------------------------------------------
echo "STEP 5: R post-processing ..." | tee -a "${LOG}"
Rscript --vanilla "${R_SCRIPT}" \
    "${RAW_TSV}" \
    "${OUT_DIR}" \
    "${YEAR_START}" \
    "${YEAR_END}" >> "${LOG}" 2>&1

echo "DONE $(date -Iseconds) 08.0.4 TerraClimate extraction" | tee -a "${LOG}"
echo "Outputs: ${DATA_DIR}/"

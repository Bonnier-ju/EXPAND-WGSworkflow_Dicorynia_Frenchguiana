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
# gdallocationinfo and python3 available system-wide

# -------------------------------------------------------------------
# STEP 1: Define paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
GEOLOC="${BASE}/metadata/geoloc_site.csv"

MANUAL_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/manual_variables"
HABITAT_RASTER="${MANUAL_DIR}/Habitats forestiers.tif"
PEDOL_SHP="${MANUAL_DIR}/Pedologie/pedologie_guyane.shp"

# Copernicus DEM GLO-30 tiles (~30 m, Float32, metres EGM2008 geoid)
DEM_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/data/elevation/copernicus_dem30"

LOG_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection/logs"
LOG_FILE="${LOG_DIR}/08.0.3-manual_variables_extraction.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.0.3-manual_variables_extraction.R"

RAW_TSV="${MANUAL_DIR}/forest_habitat_elev_raw.tsv"
COMBINED_TSV="${MANUAL_DIR}/combined_raw.tsv"
R_OUT_DIR="${BASE}/08-gea/workflow_full/08.1-env_variable_selection"

# Sites to exclude (outgroup and NA coordinates)
EXCLUDE_SITES="Cameroun_Benin"

# Minimum number of sites a category must appear in to generate a dummy variable.
MIN_SITES_DUMMY=2

mkdir -p "$LOG_DIR" "$DEM_DIR"
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
# STEP 2: Download Copernicus DEM GLO-30 tiles for French Guiana sites
# -------------------------------------------------------------------
# Resolution: ~30 m (1 arcsec), Float32, values in metres (EGM2008 geoid).
# Tiles covering French Guiana sites (lat 3-6°N, lon 52-55°W).
# Source: AWS Open Data — https://registry.opendata.aws/copernicus-dem/
# Tile naming: Copernicus_DSM_COG_10_N{lat:02d}_00_W{lon:03d}_00_DEM
# Tiles that do not exist in the dataset (e.g. ocean-only tiles) are skipped
# gracefully — sites falling in missing tiles will receive elevation = NA.
echo "STEP 2: downloading Copernicus DEM GLO-30 tiles" | tee -a "$LOG_FILE"

# Regional endpoint required — global s3.amazonaws.com returns 512-byte XML errors
COPERNICUS_BASE="https://copernicus-dem-30m.s3.eu-central-1.amazonaws.com"

for TILE in N03_W053 N03_W054 N04_W052 N04_W053 N04_W054 N05_W054 N05_W055; do
    TILENAME="Copernicus_DSM_COG_10_${TILE}_00_DEM"
    TILEPATH="${DEM_DIR}/${TILENAME}.tif"
    if [ -s "$TILEPATH" ] && [ "$(stat -c%s "$TILEPATH")" -gt 1048576 ]; then
        echo "  ${TILE}: already present, skipping" | tee -a "$LOG_FILE"
    else
        rm -f "$TILEPATH"
        echo "  ${TILE}: downloading..." | tee -a "$LOG_FILE"
        curl -sSL --retry 3 --retry-delay 5 \
            "${COPERNICUS_BASE}/${TILENAME}/${TILENAME}.tif" \
            -o "$TILEPATH" >> "$LOG_FILE" 2>&1
        # Validate: real TIF tiles are several MB; <1 MB indicates an error response
        if [ -s "$TILEPATH" ] && [ "$(stat -c%s "$TILEPATH")" -gt 1048576 ]; then
            echo "  ${TILE}: done ($(du -sh "$TILEPATH" | cut -f1))" | tee -a "$LOG_FILE"
        else
            rm -f "$TILEPATH"
            echo "  ${TILE}: WARN — tile not available, elevation will be NA for sites in this tile" | tee -a "$LOG_FILE"
        fi
    fi
done
echo "STEP 2 completed" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Extract forest habitat code and elevation at site coordinates
# -------------------------------------------------------------------
# "Habitats forestiers.tif": Float32 raster, nominal categorical codes
# (e.g. 41220, 41420, 41610). No scale factor — values are integer codes.
# Elevation: Copernicus DEM GLO-30, Float32, metres. Tile determined
# dynamically from site coordinates (tile_lat = floor(lat), tile_lon = abs(floor(lon))).
echo "STEP 3: extracting forest habitat codes and elevation" | tee -a "$LOG_FILE"

echo -e "site\tlat\tlong\tforest_habitat_code\televation" > "$RAW_TSV"

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

    # Forest habitat
    hab=$(gdallocationinfo -valonly -wgs84 "$HABITAT_RASTER" "$lon" "$lat" 2>/dev/null \
          | tr -d '[:space:]')
    [ -z "$hab" ] && hab="NA"

    # Elevation from Copernicus DEM GLO-30
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

    echo -e "${site}\t${lat}\t${lon}\t${hab}\t${elev}" >> "$RAW_TSV"
    echo "  ${site}: forest_habitat = ${hab}, elevation = ${elev} m" | tee -a "$LOG_FILE"
done

echo "STEP 3 completed: $RAW_TSV" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: Extract pedology label at site coordinates using Python
# -------------------------------------------------------------------
# Point-in-polygon query on pedologie_guyane.shp (WGS84, field: legende_si).
# Pure Python — reads .shp/.dbf directly, no external packages required.
# For sites outside all polygons the nearest polygon centroid is used.
echo "STEP 4: extracting pedology labels (pure-Python spatial join)" | tee -a "$LOG_FILE"

python3 - "$RAW_TSV" "$PEDOL_SHP" "$COMBINED_TSV" << 'PYEOF'
import sys, csv, struct

raw_tsv, shp_path, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
dbf_path = shp_path.replace('.shp', '.dbf')

# ---- read DBF -------------------------------------------------------
def read_dbf(path):
    records = []
    with open(path, 'rb') as f:
        f.read(4)
        n_records   = struct.unpack('<I', f.read(4))[0]
        header_size = struct.unpack('<H', f.read(2))[0]
        f.read(22)
        fields = []
        while True:
            d = f.read(32)
            if not d or d[0] == 0x0D:
                break
            name   = d[:11].split(b'\x00')[0].decode('latin-1').strip()
            length = d[16]
            fields.append((name, length))
        f.seek(header_size)
        for _ in range(n_records):
            flag = f.read(1)
            if not flag:
                break
            row = {}
            for name, length in fields:
                row[name] = f.read(length).decode('latin-1').strip()
            if flag != b'*':
                records.append(row)
    return records

# ---- read SHP polygons ----------------------------------------------
def read_polygons(path):
    polys = []
    with open(path, 'rb') as f:
        f.seek(100)
        while True:
            hdr = f.read(8)
            if len(hdr) < 8:
                break
            cl   = struct.unpack('>I', hdr[4:8])[0] * 2
            data = f.read(cl)
            if len(data) < 4:
                break
            stype = struct.unpack('<I', data[:4])[0]
            if stype not in (5, 15, 25):
                polys.append(None)
                continue
            nparts = struct.unpack('<I', data[36:40])[0]
            npts   = struct.unpack('<I', data[40:44])[0]
            parts  = [struct.unpack('<I', data[44+i*4:48+i*4])[0] for i in range(nparts)]
            off    = 44 + nparts * 4
            pts    = [struct.unpack('<dd', data[off+i*16:off+i*16+16]) for i in range(npts)]
            rings  = [pts[parts[i]: parts[i+1] if i+1 < nparts else npts] for i in range(nparts)]
            polys.append(rings)
    return polys

# ---- geometry helpers -----------------------------------------------
def pip(x, y, ring):
    inside, j = False, len(ring) - 1
    for i in range(len(ring)):
        xi, yi = ring[i]; xj, yj = ring[j]
        if ((yi > y) != (yj > y)) and x < (xj - xi) * (y - yi) / (yj - yi) + xi:
            inside = not inside
        j = i
    return inside

def in_poly(x, y, rings):
    if rings is None or not pip(x, y, rings[0]):
        return False
    return all(not pip(x, y, h) for h in rings[1:])

def centroid(rings):
    pts = rings[0]
    return sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)

# ---- main -----------------------------------------------------------
records = read_dbf(dbf_path)
polys   = read_polygons(shp_path)

with open(raw_tsv) as fin, open(out_tsv, 'w', newline='') as fout:
    reader = csv.DictReader(fin, delimiter='\t')
    writer = csv.writer(fout, delimiter='\t')
    writer.writerow(['site', 'lat', 'long', 'forest_habitat_code', 'elevation', 'pedology_label'])
    for row in reader:
        x, y  = float(row['long']), float(row['lat'])
        label = None
        for rec, poly in zip(records, polys):
            if in_poly(x, y, poly):
                label = rec.get('legende_si', '')
                break
        if not label:
            best_d, best_lab = float('inf'), 'NA'
            for rec, poly in zip(records, polys):
                if poly:
                    cx, cy = centroid(poly)
                    d = (x - cx)**2 + (y - cy)**2
                    if d < best_d:
                        best_d, best_lab = d, rec.get('legende_si', 'NA')
            label = best_lab
        writer.writerow([row['site'], row['lat'], row['long'],
                         row['forest_habitat_code'], row['elevation'], label or 'NA'])
        print(f"  {row['site']}: pedology = {label}", flush=True)

print(f"Pedology join completed: {out_tsv}", flush=True)
PYEOF

echo "STEP 4 completed: $COMBINED_TSV" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 5: R — dummy coding for forest habitat and pedology
# -------------------------------------------------------------------
# R reads the combined TSV. Elevation is passed through as a continuous
# variable — no dummy coding applied.
echo "STEP 5: R dummy coding" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$COMBINED_TSV" \
    "$R_OUT_DIR" \
    "$MIN_SITES_DUMMY" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.0.3 manual variables extraction completed" | tee -a "$LOG_FILE"
echo "Outputs: ${R_OUT_DIR}/data/manual_variables/"

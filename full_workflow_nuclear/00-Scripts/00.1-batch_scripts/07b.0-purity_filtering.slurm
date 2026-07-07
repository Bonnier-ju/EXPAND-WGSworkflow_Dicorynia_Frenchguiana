#!/bin/bash
#SBATCH -p workq
#SBATCH --time=00:30:00
#SBATCH -J 07b.0-purity_filtering
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
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# -------------------------------------------------------------------
# STEP 1: Define paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
WD="${BASE}/07b-population_structure_pur_inds"

Q_FILE="${WD}/input/admixture_input.3.named.Q"
FAM_FILE="${WD}/input/admixture_input.fam"
META_FILE="${BASE}/metadata/sample_sheet_complete_filtered.csv"
GEOLOC_FILE="${BASE}/metadata/geoloc_site.csv"
SITES_CLUSTERS="${BASE}/metadata/sites_by_clusters.csv"
SITES_COLORS="${BASE}/metadata/sites_couleurs.csv"

OUT_LISTS="${WD}/07b.0-purity_filtering/individual_lists"
OUT_SUMMARY="${WD}/07b.0-purity_filtering/summary"
OUT_PLOTS="${WD}/07b.0-purity_filtering/plots"
LOG_FILE="${WD}/logs/07b.0-purity_filtering.log"

R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07b.0-purity_filtering.R"

echo "START $(date -Iseconds) 07b.0 purity filtering" > "$LOG_FILE"

for f in "$Q_FILE" "$FAM_FILE" "$META_FILE" "$GEOLOC_FILE" "$SITES_CLUSTERS" "$SITES_COLORS"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: required input missing: $f" | tee -a "$LOG_FILE"
        exit 1
    fi
done

# -------------------------------------------------------------------
# STEP 2: R — purity filtering
# -------------------------------------------------------------------
echo "STEP 2: purity filtering (thresholds: 0.70 0.80 0.90 0.95)" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$Q_FILE" \
    "$FAM_FILE" \
    "$META_FILE" \
    "$GEOLOC_FILE" \
    "$SITES_CLUSTERS" \
    "$SITES_COLORS" \
    "$OUT_LISTS" \
    "$OUT_SUMMARY" \
    "$OUT_PLOTS" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 07b.0 purity filtering completed" | tee -a "$LOG_FILE"
echo "Outputs: ${WD}/00-purity_filtering/"

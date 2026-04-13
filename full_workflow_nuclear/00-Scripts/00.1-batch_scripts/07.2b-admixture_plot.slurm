#!/bin/bash
#SBATCH -p workq
#SBATCH --time=01:00:00
#SBATCH -J 07.2b-admixture_plot
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
OUT_DIR="${BASE}/07-population_structure/workflow_full/07.2-admixture"
MAP_FILE="${OUT_DIR}/sample_population_map.tsv"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.2-admixture_plot.R"

K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-10}"

if [ ! -s "$MAP_FILE" ]; then
    echo "ERROR: population map not found: $MAP_FILE"
    exit 1
fi

echo "START $(date -Iseconds) 07.2b ADMIXTURE plotting (3 groups)"

for group in all_samples no_outgroup guiana_only; do
    group_dir="${OUT_DIR}/${group}"
    if [ ! -d "$group_dir" ]; then
        echo "WARN: group directory not found, skipping: $group_dir"
        continue
    fi
    echo "Plotting group: ${group}"
    Rscript --vanilla "$R_SCRIPT" \
        "$group_dir" \
        "$MAP_FILE" \
        "$K_MIN" \
        "$K_MAX" \
        "admixture_input" \
        || echo "WARN: R plotting failed for group ${group}"
done

echo "DONE $(date -Iseconds) — plots written to ${OUT_DIR}/<group>/plots/"

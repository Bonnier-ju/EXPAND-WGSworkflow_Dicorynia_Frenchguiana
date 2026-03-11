#!/bin/bash
#SBATCH -p workq
#SBATCH --time=02:00:00
#SBATCH -J 04-coverage_summary_plot
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.out
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u
set -o pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup (Genotoul best practice)
# -------------------------------------------------------------------
module purge
module load statistics/R/4.3.0

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# To run on the test dataset, change WORKFLOW:
#   WORKFLOW="workflow_test"
# -------------------------------------------------------------------
WORKFLOW="workflow_full"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
INPUT_TSV="${BASE_DIR}/04-coverage_mapping/${WORKFLOW}/mapping_coverage_summary.tsv"
OUTPUT_DIR="${BASE_DIR}/04-coverage_mapping/${WORKFLOW}/plots"
R_SCRIPT="${BASE_DIR}/00-Scripts/00.2-R_scripts/04.1-coverage_summary_plot.R"

if [ ! -f "$INPUT_TSV" ]; then
    echo "ERROR: input TSV not found: $INPUT_TSV" >&2
    exit 1
fi

if [ ! -f "$R_SCRIPT" ]; then
    echo "ERROR: R script not found: $R_SCRIPT" >&2
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Run R script in batch mode
# -------------------------------------------------------------------
Rscript --vanilla "$R_SCRIPT" "$INPUT_TSV" "$OUTPUT_DIR"

echo "Coverage plotting completed."

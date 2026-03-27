#!/bin/bash
#SBATCH -p workq
#SBATCH --time=00:30:00
#SBATCH -J cp-07.7-samova-input
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
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

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
WORKFLOW="workflow_full"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
INPUT_ALN="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/alignments_phylo/alignment_trimmed.fasta"
METADATA="${BASE_DIR}/metadata/sample_sheet_chloroplastic_complete_filtered.csv"
GEOLOC="${BASE_DIR}/metadata/geoloc_site.csv"
OUTPUT_DIR="${BASE_DIR}/07-phylogeny_phylogeography_structure/${WORKFLOW}/07.7-samova_input"
R_SCRIPT="${BASE_DIR}/00-Scripts/00.2-R_scripts/07.7-samova_input_prep_cp.R"

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$INPUT_ALN" ];  then echo "ERROR: alignment not found: $INPUT_ALN" >&2; exit 1; fi
if [ ! -s "$METADATA" ];   then echo "ERROR: metadata not found: $METADATA" >&2; exit 1; fi
if [ ! -s "$GEOLOC" ];     then echo "ERROR: geoloc file not found: $GEOLOC" >&2; exit 1; fi

echo "START $(date -Iseconds) samova_input_prep cp"

# -------------------------------------------------------------------
# STEP 2: Run R script
# -------------------------------------------------------------------
Rscript --vanilla "$R_SCRIPT" \
    "$INPUT_ALN" \
    "$METADATA" \
    "$GEOLOC" \
    "$OUTPUT_DIR"

# Output verification
ok=0; fail=0
for f in dicorynia_cpDNA.geo dicorynia_cpDNA.arp \
         samova_site_summary.tsv session_info.txt; do
    if [ -s "${OUTPUT_DIR}/${f}" ]; then echo "  [OK] $f"; ok=$((ok + 1))
    else echo "  [MISSING] $f"; fail=$((fail + 1)); fi
done

if [ "$fail" -gt 0 ]; then
    echo "ERROR: $fail output(s) missing" >&2; exit 1
fi
echo "DONE $(date -Iseconds) ${ok} outputs saved to ${OUTPUT_DIR}"

# -------------------------------------------------------------------
# Description
# Generates SAMOVA 2.0 input files from the cpDNA trimmed alignment.
# Excludes: Angela site (uncertain cluster assignment) and outgroups
#   (Cameroun_Benin, Herbier).
# Inputs:
#   - alignment_trimmed.fasta : trimmed consensus alignment (from 06.6)
#   - sample_sheet_chloroplastic_complete_filtered.csv : metadata
#   - geoloc_site.csv : GPS coordinates (19 FG sites)
# Outputs (07.7-samova_input/):
#   dicorynia_cpDNA.geo    : SAMOVA geographic coordinates file
#                            format: N<TAB>"Site"<TAB>lon<TAB>lat<TAB>1
#   dicorynia_cpDNA.arp    : Arlequin 3.5 DNA haplotype data
#                            HaplotypeDefinition + SampleData per site
#   samova_site_summary.tsv : QC table of n_samples and n_haplotypes per site
#   session_info.txt        : R session info for reproducibility
# -------------------------------------------------------------------

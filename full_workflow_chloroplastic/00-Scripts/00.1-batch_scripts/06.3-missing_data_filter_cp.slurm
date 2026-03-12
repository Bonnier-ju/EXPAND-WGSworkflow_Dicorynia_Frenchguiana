#!/bin/bash
#SBATCH -p workq
#SBATCH --time=02:00:00
#SBATCH -J cp-06.3-missing_data_filter
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# To run on the test dataset, change WORKFLOW and SAMPLE_SHEET:
#   WORKFLOW="workflow_test_9ind"
#   SAMPLE_SHEET="...sample_sheet_chloroplastic_test.csv"
# Thresholds (haploid chloroplast):
#   MAX_SITE_MISSING  : max missing rate per site (0.2 = 20%)
#   MAX_SAMPLE_MISSING: max missing rate per sample (0.5 = 50%)
#   MIN_MAC           : minimum minor allele count (1 = at least 1 copy)
# -------------------------------------------------------------------
WORKFLOW="workflow_full"
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete_filtered.csv"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
INPUT_VCF="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/normalization/variants.normalized.vcf.gz"
OUTPUT_DIR="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/missing_data_filter"
LOG_FILE="${OUTPUT_DIR}/missing_data_filter_cp.log"

MAX_SITE_MISSING=0.2
MAX_SAMPLE_MISSING=0.5
MIN_MAC=1

# Samples explicitly excluded regardless of missing rate.
# All low-quality samples already removed via filtered sample sheet (164 samples).
FORCE_EXCLUDE=""

SAMPLES_TO_EXCLUDE="${OUTPUT_DIR}/samples_high_missing.txt"
OUTPUT_VCF="${OUTPUT_DIR}/variants.filtered.vcf.gz"

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found or empty: $INPUT_VCF" >&2; exit 1
fi
if [ -s "$OUTPUT_VCF" ]; then
    echo "SKIP: output already exists: $OUTPUT_VCF"; exit 0
fi

echo "START $(date -Iseconds) missing_data_filter cp" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Identify samples with too much missing data
#   bcftools stats -s- : per-sample statistics
#   Extract missing genotype rate per sample
# -------------------------------------------------------------------
echo "INFO: computing per-sample missing rate" >> "$LOG_FILE"
bcftools stats -s - "$INPUT_VCF" 2>> "$LOG_FILE" | \
    awk -v thresh="$MAX_SAMPLE_MISSING" '
        /^PSC/ {
            sample=$3; total=$4+$5+$6; missing=$6
            if (total > 0 && (missing/total) > thresh) print sample
        }
    ' > "$SAMPLES_TO_EXCLUDE"

# Append forced exclusions to the exclusion list
for s in $FORCE_EXCLUDE; do
    echo "$s" >> "$SAMPLES_TO_EXCLUDE"
done
# Remove duplicates
sort -u "$SAMPLES_TO_EXCLUDE" -o "$SAMPLES_TO_EXCLUDE"

n_excluded=$(wc -l < "$SAMPLES_TO_EXCLUDE")
echo "INFO: ${n_excluded} samples excluded total (missing > ${MAX_SAMPLE_MISSING} or forced)" >> "$LOG_FILE"
echo "INFO: forced exclusions (low cp coverage): ${FORCE_EXCLUDE}" >> "$LOG_FILE"
if [ "$n_excluded" -gt 0 ]; then
    echo "INFO: excluded samples:" >> "$LOG_FILE"
    cat "$SAMPLES_TO_EXCLUDE" >> "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 3: Filter sites and samples
#   -e 'F_MISSING > MAX_SITE_MISSING' : exclude sites with too much missing data
#   -c MIN_MAC:minor                  : retain sites with ≥1 minor allele
#   --samples-file with ^ prefix      : exclude high-missing samples
# -------------------------------------------------------------------
echo "INFO: filtering sites (missing > ${MAX_SITE_MISSING}) and MAC < ${MIN_MAC}" >> "$LOG_FILE"

EXCLUDE_OPT=""
if [ "$n_excluded" -gt 0 ]; then
    EXCLUDE_OPT="--samples-file ^${SAMPLES_TO_EXCLUDE}"
fi

bcftools view \
    -e "F_MISSING > ${MAX_SITE_MISSING}" \
    -c "${MIN_MAC}:minor" \
    ${EXCLUDE_OPT} \
    -O z -o "$OUTPUT_VCF" \
    "$INPUT_VCF" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 4: Index the filtered VCF
# -------------------------------------------------------------------
echo "INFO: indexing filtered VCF" >> "$LOG_FILE"
bcftools index --tbi "$OUTPUT_VCF" >> "$LOG_FILE" 2>&1

if [ -s "$OUTPUT_VCF" ]; then
    n_var=$(bcftools stats "$OUTPUT_VCF" 2>/dev/null | grep "^SN" | grep "number of records" | awk '{print $NF}')
    n_samples=$(bcftools query -l "$OUTPUT_VCF" 2>/dev/null | wc -l)
    echo "DONE $(date -Iseconds) variants=${n_var} samples=${n_samples} output=${OUTPUT_VCF}" | tee -a "$LOG_FILE"
else
    echo "ERROR: missing data filter output missing" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# Description
# Step 1: bcftools stats -s- → per-sample missing genotype rate
# Step 2: Exclude samples with missing rate > MAX_SAMPLE_MISSING (default 50%)
# Step 3: bcftools view -e 'F_MISSING > 0.2' → remove incomplete sites
#         -c 1:minor → retain polymorphic sites (MAC ≥ 1)
# Step 4: TBI indexing
# Haploid chloroplast: no HWE filter, no strict minimum MAF filter
# -------------------------------------------------------------------

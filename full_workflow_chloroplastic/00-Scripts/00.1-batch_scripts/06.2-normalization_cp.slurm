#!/bin/bash
#SBATCH -p workq
#SBATCH --time=02:00:00
#SBATCH -J cp-06.2-normalization
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
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# To run on the test dataset, change WORKFLOW and SAMPLE_SHEET:
#   WORKFLOW="workflow_test_9ind"
#   SAMPLE_SHEET="...sample_sheet_chloroplastic_test.csv"
# -------------------------------------------------------------------
WORKFLOW="workflow_full"
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete.csv"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
REF_GENOME="$(awk -F, 'NR==2{print $7}' "$SAMPLE_SHEET")"
INPUT_VCF="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/hard_filter/variants.hard_filtered.pass.vcf.gz"
OUTPUT_DIR="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/normalization"
LOG_FILE="${OUTPUT_DIR}/normalization_cp.log"

OUTPUT_VCF="${OUTPUT_DIR}/variants.normalized.vcf.gz"

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found or empty: $INPUT_VCF" >&2; exit 1
fi
if [ -z "$REF_GENOME" ] || [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found: $REF_GENOME" >&2; exit 1
fi
if [ -s "$OUTPUT_VCF" ]; then
    echo "SKIP: output already exists: $OUTPUT_VCF"; exit 0
fi

echo "START $(date -Iseconds) normalization cp" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Split multiallelic variants then normalize
#   -m -any  : split each alternate allele onto a separate line
#   -f REF   : left-align indels against the reference
# -------------------------------------------------------------------
echo "INFO: decomposing multiallelic sites and left-aligning indels" >> "$LOG_FILE"
bcftools norm -m -any "$INPUT_VCF" 2>> "$LOG_FILE" | \
    bcftools norm -f "$REF_GENOME" -O z -o "$OUTPUT_VCF" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 3: Index the normalized VCF
# -------------------------------------------------------------------
echo "INFO: indexing normalized VCF" >> "$LOG_FILE"
bcftools index --tbi "$OUTPUT_VCF" >> "$LOG_FILE" 2>&1

if [ -s "$OUTPUT_VCF" ]; then
    n_var=$(bcftools stats "$OUTPUT_VCF" 2>/dev/null | grep "^SN" | grep "number of records" | awk '{print $NF}')
    echo "DONE $(date -Iseconds) variants_normalized=${n_var} output=${OUTPUT_VCF}" | tee -a "$LOG_FILE"
else
    echo "ERROR: normalization output missing" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# Description
# Step 1: bcftools norm -m -any  → split multiallelic sites into biallelic records
# Step 2: bcftools norm -f REF   → left-align indels against the reference
# Step 3: bcftools index --tbi   → TBI index of the compressed VCF
# Output: normalized VCF ready for missing data filtering (06.3)
# -------------------------------------------------------------------

#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=50:00:00
#SBATCH -J cp-05.3-genotypegvcfs
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load devel/python/Python-3.11.1
module load devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs (chloroplastic workflow full)
# -------------------------------------------------------------------
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete.csv"
REF_GENOME="$(awk -F, 'NR==2{print $7}' "$SAMPLE_SHEET")"
DB_PATH="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/05-variant_calling/workflow_full/genomicsdbimport/GATK_database_dicorynia_cp"
OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/05-variant_calling/workflow_full/genotypegvcfs"
TMP_DIR="${OUTPUT_DIR}/tmp"
LOG_FILE="${OUTPUT_DIR}/genotypegvcfs_cp.log"

OUTPUT_VCF="${OUTPUT_DIR}/jointgenotypes.raw.vcf.gz"
OUTPUT_TBI="${OUTPUT_VCF}.tbi"

JAVA_MEM_GB=28
SAMPLE_PLOIDY=1

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found or empty: $SAMPLE_SHEET" >&2
    exit 1
fi
if [ -z "$REF_GENOME" ] || [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found from sample sheet: $REF_GENOME" >&2
    exit 1
fi
if [ ! -d "$DB_PATH" ]; then
    echo "ERROR: GenomicsDB workspace not found: $DB_PATH" >&2
    exit 1
fi

echo "START $(date -Iseconds) GenotypeGVCFs cp" > "$LOG_FILE"

# Ensure reference indexes exist
if [ ! -s "${REF_GENOME}.fai" ]; then
    echo "INFO: creating reference fai index" >> "$LOG_FILE"
    samtools faidx "$REF_GENOME" >> "$LOG_FILE" 2>&1
fi

REF_DICT="${REF_GENOME%.fa}.dict"
if [ "$REF_DICT" = "$REF_GENOME.dict" ]; then
    REF_DICT="${REF_GENOME%.fasta}.dict"
fi
if [ ! -s "$REF_DICT" ]; then
    echo "INFO: creating reference sequence dictionary" >> "$LOG_FILE"
    gatk CreateSequenceDictionary -R "$REF_GENOME" -O "$REF_DICT" >> "$LOG_FILE" 2>&1
fi

# Skip if output already present
if [ -s "$OUTPUT_VCF" ] && [ -s "$OUTPUT_TBI" ]; then
    echo "SKIP: output already exists: $OUTPUT_VCF" | tee -a "$LOG_FILE"
    exit 0
fi

# -------------------------------------------------------------------
# STEP 2: Joint genotyping (chloroplastic adapted)
# -------------------------------------------------------------------
gatk --java-options "-Xms${JAVA_MEM_GB}G -Xmx${JAVA_MEM_GB}G -XX:ParallelGCThreads=4" GenotypeGVCFs \
    -R "$REF_GENOME" \
    -V "gendb://${DB_PATH}" \
    -O "$OUTPUT_VCF" \
    --sample-ploidy "$SAMPLE_PLOIDY" \
    --genomicsdb-use-bcf-codec true \
    --create-output-variant-index true \
    --tmp-dir "$TMP_DIR" >> "$LOG_FILE" 2>&1

if [ -s "$OUTPUT_VCF" ] && [ -s "$OUTPUT_TBI" ]; then
    echo "DONE $(date -Iseconds) output=${OUTPUT_VCF}" | tee -a "$LOG_FILE"
else
    echo "ERROR: GenotypeGVCFs output missing" | tee -a "$LOG_FILE"
    exit 1
fi

rm -rf "$TMP_DIR"

# Step 1: Validate chloroplastic GenomicsDB workspace and reference files
# Step 2: Run GATK GenotypeGVCFs on gendb:// workspace
# Step 3: Use chloroplastic ploidy setting (--sample-ploidy 1)
# Step 4: Validate output VCF and index

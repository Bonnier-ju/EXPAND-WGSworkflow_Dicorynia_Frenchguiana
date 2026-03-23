#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=10-00:00:00
#SBATCH -J 05.3-genotypegvcfs
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL


set -euo pipefail

echo "======================================================"
echo "START $(date -Iseconds) GenotypeGVCFs"
echo "  SLURM job ID : ${SLURM_JOB_ID}"
echo "  CPUs         : ${SLURM_CPUS_PER_TASK}"
echo "  Memory       : 32G"
echo "======================================================"

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load devel/python/Python-3.11.1 devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
REF_GENOME="/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa"
DB_PATH="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/05-variant_calling/workflow_full/genomicsdbimport/GATK_database_dicorynia"
OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/05-variant_calling/workflow_full/genotypegvcfs"
TMP_DIR="/tmp/${USER}/genotypegvcfs_${SLURM_JOB_ID}"
LOG_FILE="${OUTPUT_DIR}/genotypegvcfs.log"

OUTPUT_VCF="${OUTPUT_DIR}/jointgenotypes.raw.vcf.gz"
OUTPUT_TBI="${OUTPUT_VCF}.tbi"

JAVA_MEM_GB=28

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

if [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found: $REF_GENOME" >&2
    exit 1
fi
if [ ! -d "$DB_PATH" ]; then
    echo "ERROR: GenomicsDB workspace not found: $DB_PATH" >&2
    exit 1
fi

echo "INFO: modules loaded"
echo "START $(date -Iseconds) GenotypeGVCFs" > "$LOG_FILE"

# Ensure reference indexes exist
if [ ! -s "${REF_GENOME}.fai" ]; then
    echo "INFO: creating reference fai index" >> "$LOG_FILE"
    samtools faidx "$REF_GENOME" >> "$LOG_FILE" 2>&1
fi
REF_DICT="${REF_GENOME%.fa}.dict"
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
# STEP 2: Joint genotyping
# -------------------------------------------------------------------
echo "INFO: running GenotypeGVCFs..." | tee -a "$LOG_FILE"
echo "  DB path    : ${DB_PATH}"     | tee -a "$LOG_FILE"
echo "  Output VCF : ${OUTPUT_VCF}"  | tee -a "$LOG_FILE"
echo "  Java heap  : ${JAVA_MEM_GB}G" | tee -a "$LOG_FILE"

gatk --java-options "-Xms${JAVA_MEM_GB}G -Xmx${JAVA_MEM_GB}G -XX:ParallelGCThreads=2" GenotypeGVCFs \
    -R "$REF_GENOME" \
    -V "gendb://${DB_PATH}" \
    -O "$OUTPUT_VCF" \
    --create-output-variant-index true \
    --tmp-dir "$TMP_DIR" \
    --verbosity WARNING \
    >> "$LOG_FILE" 2>&1

if [ -s "$OUTPUT_VCF" ] && [ -s "$OUTPUT_TBI" ]; then
    VCF_SIZE=$(du -sh "$OUTPUT_VCF" | cut -f1)
    echo "INFO: output verified — size: ${VCF_SIZE}" | tee -a "$LOG_FILE"
else
    echo "ERROR: joint genotyping output missing" | tee -a "$LOG_FILE"
    exit 1
fi

# --- Cleanup ------------------------------------------------------------------
rm -rf "$TMP_DIR"
echo "INFO: tmp cleaned up" | tee -a "$LOG_FILE"

echo "======================================================"
echo "DONE $(date -Iseconds) GenotypeGVCFs"
echo "  VCF output : ${OUTPUT_VCF} (${VCF_SIZE})"
echo "  Log        : ${LOG_FILE}"
echo "  Next step  : sbatch 06.x-variant_filtering.slurm"
echo "======================================================"
echo "DONE $(date -Iseconds)" >> "$LOG_FILE"

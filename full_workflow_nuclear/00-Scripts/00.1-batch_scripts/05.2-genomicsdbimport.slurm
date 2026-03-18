#!/bin/bash
#SBATCH -p workq
#SBATCH --time=80:00:00
#SBATCH -J 05.2-genomicsdbimport
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =============================================================================
# STEP 05.2 — GenomicsDBImport
# Aggregates 164 per-sample gVCFs into a single TileDB workspace
# for joint genotyping in step 05.3 (GenotypeGVCFs).
# =============================================================================

echo "======================================================"
echo "START $(date -Iseconds) GenomicsDBImport"
echo "  SLURM job ID : ${SLURM_JOB_ID}"
echo "  CPUs         : ${SLURM_CPUS_PER_TASK}"
echo "  Memory       : 64G"
echo "======================================================"

# --- Modules ------------------------------------------------------------------
module purge
module load devel/python/Python-3.11.1 devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/samtools/1.21
echo "INFO: modules loaded"

# --- Paths --------------------------------------------------------------------
GVCF_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/05-variant_calling/workflow_full/haplotypecaller_gvcf"
AGGREGATION_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/05-variant_calling/workflow_full/genomicsdbimport"
REF_GENOME="/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa"

SAMPLE_MAP="${AGGREGATION_DIR}/sample_map.txt"
DB_PATH="${AGGREGATION_DIR}/GATK_database_dicorynia"
SEQ_LIST="${AGGREGATION_DIR}/reference.sequences.list"
LOG_FILE="${AGGREGATION_DIR}/genomicsdbimport.log"
TMP_DIR="/tmp/${USER}/genomicsdbimport_${SLURM_JOB_ID}"

# GATK tuning — optimized for 164 samples / 64G RAM / 8 CPUs
JAVA_MEM_GB=56          # Java heap (leave ~8G for OS + JVM overhead)
READER_THREADS=8        # Parallel gVCF readers (= nb CPUs)
BATCH_SIZE=50           # Samples per DB commit (164 samples → ~4 batches)

mkdir -p "$AGGREGATION_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) GenomicsDBImport" > "$LOG_FILE"

# --- Input checks -------------------------------------------------------------
[ -d "$GVCF_DIR" ]   || { echo "ERROR: GVCF dir not found: $GVCF_DIR";    exit 1; }
[ -s "$REF_GENOME" ] || { echo "ERROR: reference not found: $REF_GENOME"; exit 1; }

if [ -e "$DB_PATH" ]; then
    echo "ERROR: GenomicsDB workspace already exists: $DB_PATH"
    echo "Remove it manually if you want a clean rebuild."
    exit 1
fi

# --- Step 1: Reference sequence list -----------------------------------------
echo "INFO: preparing reference sequence list..." | tee -a "$LOG_FILE"
if [ ! -s "${REF_GENOME}.fai" ]; then
    echo "INFO: creating .fai index..." | tee -a "$LOG_FILE"
    samtools faidx "$REF_GENOME" 2>> "$LOG_FILE"
fi
cut -f1 "${REF_GENOME}.fai" > "$SEQ_LIST"
N_CONTIGS=$(wc -l < "$SEQ_LIST")
echo "INFO: ${N_CONTIGS} contigs in sequence list" | tee -a "$LOG_FILE"

# --- Step 2: Build sample map -------------------------------------------------
echo "INFO: building sample map..." | tee -a "$LOG_FILE"
find "$GVCF_DIR" -maxdepth 1 -type f -name "*.g.vcf.gz" | sort \
    > "${AGGREGATION_DIR}/gvcf_files.list"

N_GVCF=$(wc -l < "${AGGREGATION_DIR}/gvcf_files.list")
[ "$N_GVCF" -gt 0 ] || { echo "ERROR: no gVCF files found in $GVCF_DIR"; exit 1; }
echo "INFO: found ${N_GVCF} gVCF files" | tee -a "$LOG_FILE"

: > "$SAMPLE_MAP"
while read -r gvcf; do
    sample_name="$(basename "$gvcf" .g.vcf.gz)"
    if [ ! -s "${gvcf}.tbi" ]; then
        echo "INFO: indexing ${sample_name}..." | tee -a "$LOG_FILE"
        gatk IndexFeatureFile -I "$gvcf" >> "$LOG_FILE" 2>&1
    fi
    printf "%s\t%s\n" "$sample_name" "$gvcf" >> "$SAMPLE_MAP"
done < "${AGGREGATION_DIR}/gvcf_files.list"

echo "INFO: sample map ready — $(wc -l < "$SAMPLE_MAP") samples" | tee -a "$LOG_FILE"

# --- Step 3: GenomicsDBImport -------------------------------------------------
echo "INFO: running GenomicsDBImport..." | tee -a "$LOG_FILE"
echo "  batch-size     : ${BATCH_SIZE}"      | tee -a "$LOG_FILE"
echo "  reader-threads : ${READER_THREADS}"  | tee -a "$LOG_FILE"
echo "  java heap      : ${JAVA_MEM_GB}G"    | tee -a "$LOG_FILE"
echo "  DB path        : ${DB_PATH}"         | tee -a "$LOG_FILE"

gatk --java-options "-Xms${JAVA_MEM_GB}G -Xmx${JAVA_MEM_GB}G -XX:ParallelGCThreads=4" \
    GenomicsDBImport \
    --genomicsdb-workspace-path "$DB_PATH" \
    --sample-name-map "$SAMPLE_MAP" \
    -L "$SEQ_LIST" \
    --batch-size "$BATCH_SIZE" \
    --reader-threads "$READER_THREADS" \
    --genomicsdb-shared-posixfs-optimizations true \
    --tmp-dir "$TMP_DIR" \
    --verbosity WARNING \
    2>> "$LOG_FILE"

# --- Verification -------------------------------------------------------------
if [ -d "$DB_PATH" ]; then
    DB_SIZE=$(du -sh "$DB_PATH" | cut -f1)
    echo "INFO: GenomicsDB created — size: ${DB_SIZE}" | tee -a "$LOG_FILE"
else
    echo "ERROR: DB_PATH not found after import — check log: $LOG_FILE"
    exit 1
fi

# --- Cleanup ------------------------------------------------------------------
rm -rf "$TMP_DIR"
echo "INFO: tmp directory cleaned up" | tee -a "$LOG_FILE"

echo "======================================================"
echo "DONE $(date -Iseconds) GenomicsDB = ${DB_PATH}"
echo "  DB size    : ${DB_SIZE}"
echo "  Log        : ${LOG_FILE}"
echo "  Next step  : sbatch 05.3-genotypegvcfs.slurm"
echo "======================================================"
echo "DONE $(date -Iseconds)" >> "$LOG_FILE"

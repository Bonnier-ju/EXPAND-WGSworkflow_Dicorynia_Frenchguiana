#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J cp-03-bwa_mem2_alignment
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.err
#SBATCH --array=1-4
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# Chloroplastic alignment from already-trimmed reads
# - Default mode: test set (7 individuals)
# - Reuses trimmed reads generated in nuclear workflow
# -------------------------------------------------------------------

set -euo pipefail

module purge
module load bioinfo/bwa-mem2/2.2.1
module load bioinfo/samtools/1.21

TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-4}"
SORT_THREADS=$((TOTAL_THREADS / 3))
[ "$SORT_THREADS" -lt 1 ] && SORT_THREADS=1
VIEW_THREADS=1
BWA_THREADS=$((TOTAL_THREADS - SORT_THREADS - VIEW_THREADS))
[ "$BWA_THREADS" -lt 1 ] && BWA_THREADS=1
SORT_MEM_PER_THREAD="1G"

SAMPLE_SHEET="${SAMPLE_SHEET:-/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_test.csv}"
READS_DIR="${READS_DIR:-/home/jbonnier/work/chapitre_5/full_workflow_nuclear/02-filtering_indexing/workflow_full/cleaned_reads}"
ALIGN_DIR="${ALIGN_DIR:-/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/03-alignment/workflow_test_9ind/alignments}"
LOG_DIR="${ALIGN_DIR}/logs"
TMP_DIR="${ALIGN_DIR}/tmp"

mkdir -p "$ALIGN_DIR" "$LOG_DIR" "$TMP_DIR"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found or empty: $SAMPLE_SHEET" >&2
    exit 1
fi

FILES_LIST="${ALIGN_DIR}/alignment_files.list"
# Columns: sample_id (1), reference_genome (7), read_group (9)
tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1 "\t" $7 "\t" $9}' > "$FILES_LIST"

if [ ! -s "$FILES_LIST" ]; then
    echo "ERROR: empty file list generated from $SAMPLE_SHEET" >&2
    exit 1
fi

# Ensure reference index exists (single chloroplast reference expected)
REF_GENOME=$(awk -F, 'NR==2{print $7}' "$SAMPLE_SHEET")
if [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found: $REF_GENOME" >&2
    exit 1
fi

if [ ! -s "${REF_GENOME}.bwt.2bit.64" ]; then
    echo "Building bwa-mem2 index for reference: $REF_GENOME"
    bwa-mem2 index "$REF_GENOME"
fi

if [ ! -s "${REF_GENOME}.fai" ]; then
    echo "Building samtools faidx: ${REF_GENOME}.fai"
    samtools faidx "$REF_GENOME"
fi

SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-2}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

sed -n "${START_LINE},${END_LINE}p" "$FILES_LIST" | while IFS=$'\t' read -r sample_id reference read_group; do
    [ -z "$sample_id" ] && continue

    R1="${READS_DIR}/${sample_id}_R1.trimmed.fastq.gz"
    R2="${READS_DIR}/${sample_id}_R2.trimmed.fastq.gz"
    BAM="${ALIGN_DIR}/${sample_id}.sorted.bam"
    LOG="${LOG_DIR}/${sample_id}.alignment.log"

    if [ ! -s "$reference" ]; then
        echo "ERROR: reference genome missing for $sample_id: $reference" | tee -a "$LOG"
        continue
    fi

    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
        echo "ERROR: trimmed reads missing for $sample_id: $R1 $R2" | tee -a "$LOG"
        continue
    fi

    if [ -s "$BAM" ]; then
        echo "SKIP: BAM already exists for $sample_id: $BAM" | tee -a "$LOG"
        continue
    fi

    if [ -n "$read_group" ]; then
        RG="$read_group"
    else
        RG="@RG\tID:${sample_id}\tSM:${sample_id}\tLB:chloroplastic\tPL:ILLUMINA"
    fi

    echo "START $(date -Iseconds) sample=$sample_id" >> "$LOG"
    echo "THREADS sample=${sample_id} bwa=${BWA_THREADS} view=${VIEW_THREADS} sort=${SORT_THREADS} sort_mem=${SORT_MEM_PER_THREAD}" >> "$LOG"

    if { bwa-mem2 mem -M -R "$RG" -t "$BWA_THREADS" "$reference" "$R1" "$R2" |
         samtools view -b -@ "$VIEW_THREADS" - |
         samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM_PER_THREAD" -T "${TMP_DIR}/${sample_id}" -o "$BAM" - &&
         samtools index "$BAM"; } >> "$LOG" 2>&1; then
        if [ -s "$BAM" ] && [ -s "${BAM}.bai" ]; then
            echo "DONE $(date -Iseconds) sample=$sample_id bam=$BAM" >> "$LOG"
        else
            echo "ERROR: BAM/BAI missing after alignment for $sample_id" | tee -a "$LOG"
        fi
    else
        echo "ERROR: alignment pipeline failed for $sample_id" | tee -a "$LOG"
    fi
done

echo "All requested samples processed for task ${TASK_ID}."

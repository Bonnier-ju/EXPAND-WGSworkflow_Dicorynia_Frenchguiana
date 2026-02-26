#!/bin/bash
#SBATCH -p workq
#SBATCH --time=48:00:00
#SBATCH -J 03-bwa_mem2_alignment
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.out
#SBATCH --array=1-36
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# Clean environment
module purge

# Load bwa-mem2 and samtools
module load bioinfo/bwa-mem2/2.2.1
module load bioinfo/samtools/1.21

set -u
set -o pipefail

# -------------------------------------------------------------------
# STEP 0.1: Thread/memory tuning for alignment pipeline
# -------------------------------------------------------------------
TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-4}"
SORT_THREADS=$((TOTAL_THREADS / 3))
[ "$SORT_THREADS" -lt 1 ] && SORT_THREADS=1
VIEW_THREADS=1
BWA_THREADS=$((TOTAL_THREADS - SORT_THREADS - VIEW_THREADS))
[ "$BWA_THREADS" -lt 1 ] && BWA_THREADS=1
SORT_MEM_PER_THREAD="1G"

# File paths
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_complete.csv"
READS_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/02-filtering_indexing/workflow_full/cleaned_reads"
ALIGN_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/03-reading_sorting/workflow_full/alignments"
LOG_DIR="${ALIGN_DIR}/logs"
TMP_DIR="${ALIGN_DIR}/tmp"
mkdir -p "$ALIGN_DIR" "$LOG_DIR" "$TMP_DIR"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found or empty: $SAMPLE_SHEET" >&2
    exit 1
fi

FILES_LIST="${ALIGN_DIR}/alignment_files.list"
# Columns: sample_id (1), reference_genome (7), read_group (9), read1 (10), read2 (11)
tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1 "\t" $7 "\t" $9 "\t" $10 "\t" $11}' > "$FILES_LIST"

SAMPLES_PER_TASK=5
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

# Loop over sample sheet (skip header)
# Expect read1/read2 as the last two columns in the CSV, and reference genome as the last column
sed -n "${START_LINE},${END_LINE}p" "$FILES_LIST" | while IFS=$'\t' read -r sample_id reference read_group read1 read2; do
    [ -z "$sample_id" ] && continue
    echo "Aligning sample: $sample_id"

    # Define input files
    R1="${READS_DIR}/${sample_id}_R1.trimmed.fastq.gz"
    R2="${READS_DIR}/${sample_id}_R2.trimmed.fastq.gz"
    BAM="${ALIGN_DIR}/${sample_id}.sorted.bam"
    LOG="${LOG_DIR}/${sample_id}.alignment.log"
    # Alignment pipeline for paired-end reads using bwa-mem2 and samtools
    # --------------------------------------------------------------------
    if [ ! -f "$reference" ]; then
        echo "ERROR: reference genome not found: $reference" | tee -a "$LOG"
        continue
    fi
    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
        echo "ERROR: missing reads for $sample_id: $R1 $R2" | tee -a "$LOG"
        continue
    fi
    if [ -s "$BAM" ]; then
        echo "SKIP: BAM already exists for $sample_id: $BAM" | tee -a "$LOG"
        continue
    fi

    # Use read group from sample sheet if provided, otherwise build one
    if [ -n "$read_group" ]; then
        RG="$read_group"
    else
        filename=$(basename "$read1")
        LB=$(echo "$filename" | cut -d'_' -f1,2)
        LANE=$(echo "$filename" | cut -d'_' -f3)
        FLOWCELL=$(echo "$filename" | cut -d'_' -f5 | cut -d'.' -f1)
        PU="${FLOWCELL}.${LANE}"
        PLATFORM="ILLUMINA"
        RG="@RG\tID:${sample_id}\tSM:${sample_id}\tLB:${LB}\tPU:${PU}\tPL:${PLATFORM}"
    fi

    echo "START $(date -Iseconds) sample=$sample_id" >> "$LOG"

    echo "THREADS sample=${sample_id} bwa=${BWA_THREADS} view=${VIEW_THREADS} sort=${SORT_THREADS} sort_mem=${SORT_MEM_PER_THREAD}" >> "$LOG"

    if { bwa-mem2 mem -M -R "$RG" -t "$BWA_THREADS" "$reference" "$R1" "$R2" |
         samtools view -b -@ "$VIEW_THREADS" - |
         samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM_PER_THREAD" -T "${TMP_DIR}/${sample_id}" -o "$BAM" - &&
         samtools index "$BAM"; } >> "$LOG" 2>&1; then
        if [ -s "$BAM" ]; then
            echo "DONE $(date -Iseconds) sample=$sample_id bam=$BAM" >> "$LOG"
        else
            echo "ERROR: BAM not created for $sample_id" | tee -a "$LOG"
        fi
    else
        echo "ERROR: alignment pipeline failed for $sample_id" | tee -a "$LOG"
    fi
done

echo "All samples aligned, sorted, and indexed."

# Step 1: Construct read group header
    # - Required by tools like GATK for downstream processing
# Step 2: Align reads with bwa-mem2
    # - -t ${SLURM_CPUS_PER_TASK} : use all allocated CPUs
    # - -M : mark shorter split hits as secondary (GATK compatibility)
    # - -R : include read group information in the header
# Step 3: Convert SAM to BAM
    # - -b : output BAM format
    # - -@ ${SLURM_CPUS_PER_TASK} : use allocated CPUs for conversion
# Step 4: Sort BAM by coordinates
    # - -@ ${SLURM_CPUS_PER_TASK} : use allocated CPUs for sorting
    # - -o : output sorted BAM file
    # - && : ensure that the next command runs only if the previous one succeeds
# Step 5: Index sorted BAM 
    # - Indexing is necessary for random access to the BAM file

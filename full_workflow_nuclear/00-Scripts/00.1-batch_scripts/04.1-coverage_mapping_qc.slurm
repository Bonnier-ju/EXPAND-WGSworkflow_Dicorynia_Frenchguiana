#!/bin/bash
#SBATCH -p workq
#SBATCH --time=50:00:00
#SBATCH -J 04-coverage_mapping_qc
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.out
#SBATCH --array=1-36
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u
set -o pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/samtools/1.21
module load bioinfo/Qualimap/31-08-20

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs and initialize directories
# -------------------------------------------------------------------
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_complete.csv"
ALIGN_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/03-reading_sorting/workflow_full/alignments"
OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/04-coverage_mapping/workflow_full"

FLAGSTAT_DIR="${OUTPUT_DIR}/flagstat"
DEPTH_DIR="${OUTPUT_DIR}/depth"
QUALIMAP_DIR="${OUTPUT_DIR}/qualimap"
LOG_DIR="${OUTPUT_DIR}/logs"
SUMMARY_TSV="${OUTPUT_DIR}/mapping_coverage_summary.tsv"
SUMMARY_ROWS_DIR="${OUTPUT_DIR}/summary_rows"
SAMPLES_LIST="${OUTPUT_DIR}/samples.list"
SAMPLES_PER_TASK=5
SAMTOOLS_THREADS="${SLURM_CPUS_PER_TASK:-4}"
QUALIMAP_THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$FLAGSTAT_DIR" "$DEPTH_DIR" "$QUALIMAP_DIR" "$LOG_DIR" "$SUMMARY_ROWS_DIR"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found: $SAMPLE_SHEET" >&2
    exit 1
fi

# Build sample list once; reuse it across array tasks.
if [ ! -s "$SAMPLES_LIST" ]; then
    tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1}' > "$SAMPLES_LIST"
fi

# -------------------------------------------------------------------
# STEP 2: Per-sample mapping and coverage QC
# -------------------------------------------------------------------
process_sample() {
    sample_id="$1"

    BAMFILE="${ALIGN_DIR}/${sample_id}.sorted.bam"
    FLAGSTAT_OUT="${FLAGSTAT_DIR}/${sample_id}.flagstat.txt"
    DEPTH_SUMMARY="${DEPTH_DIR}/${sample_id}.avg_depth_q20.txt"
    SAMPLE_QUALIMAP_DIR="${QUALIMAP_DIR}/${sample_id}"
    SAMPLE_LOG="${LOG_DIR}/${sample_id}.coverage_qc.log"
    SAMPLE_ROW="${SUMMARY_ROWS_DIR}/${sample_id}.tsv"

    echo "START $(date -Iseconds) sample=${sample_id}" > "$SAMPLE_LOG"

    # 2.1 Validate BAM input
    if [ ! -s "$BAMFILE" ]; then
        echo "ERROR: missing or empty BAM: $BAMFILE" | tee -a "$SAMPLE_LOG"
        echo -e "${sample_id}\tNA\tNA\tNA\tNA" > "$SAMPLE_ROW"
        return 0
    fi

    # 2.2 Ensure BAM index is available
    if [ ! -s "${BAMFILE}.bai" ]; then
        echo "INFO: BAM index missing, creating ${BAMFILE}.bai" >> "$SAMPLE_LOG"
        samtools index "$BAMFILE" >> "$SAMPLE_LOG" 2>&1
    fi

    # 2.3 Compute mapping statistics (flagstat)
    samtools flagstat -@ "$SAMTOOLS_THREADS" "$BAMFILE" > "$FLAGSTAT_OUT" 2>> "$SAMPLE_LOG"
    total_reads="$(awk '/in total/ {print $1; exit}' "$FLAGSTAT_OUT")"
    mapped_reads="$(awk '/ mapped \(/ {print $1; exit}' "$FLAGSTAT_OUT")"
    mapped_percent="$(awk -F'[()%]' '/ mapped \(/ {gsub(/ /, "", $2); print $2; exit}' "$FLAGSTAT_OUT")"

    # 2.4 Compute mean depth with mapping quality filter
    avg_depth_q20="$(samtools depth -@ "$SAMTOOLS_THREADS" -a -Q 20 "$BAMFILE" 2>> "$SAMPLE_LOG" | awk '{sum+=$3; n++} END {if (n>0) printf "%.4f", sum/n; else print "0"}')"
    echo "$avg_depth_q20" > "$DEPTH_SUMMARY"

    # 2.5 Run Qualimap bamqc report
    mkdir -p "$SAMPLE_QUALIMAP_DIR"
    qualimap bamqc \
        -bam "$BAMFILE" \
        -outdir "$SAMPLE_QUALIMAP_DIR" \
        -nt "$QUALIMAP_THREADS" \
        --java-mem-size=8G >> "$SAMPLE_LOG" 2>&1

    echo -e "${sample_id}\t${total_reads:-NA}\t${mapped_reads:-NA}\t${mapped_percent:-NA}\t${avg_depth_q20}" > "$SAMPLE_ROW"
    echo "DONE $(date -Iseconds) sample=${sample_id}" >> "$SAMPLE_LOG"
}

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    START_LINE=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
    END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

    sed -n "${START_LINE},${END_LINE}p" "$SAMPLES_LIST" | while read -r sample_id; do
        [ -z "$sample_id" ] && continue
        process_sample "$sample_id"
    done
else
    while read -r sample_id; do
        [ -z "$sample_id" ] && continue
        process_sample "$sample_id"
    done < "$SAMPLES_LIST"
fi

# -------------------------------------------------------------------
# STEP 3: Final status and summary merge
# -------------------------------------------------------------------
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo -e "sample_id\ttotal_reads\tmapped_reads\tmapped_percent\tavg_depth_q20" > "$SUMMARY_TSV"
    cat "${SUMMARY_ROWS_DIR}"/*.tsv >> "$SUMMARY_TSV"
    echo "Coverage and mapping QC completed."
    echo "Summary: $SUMMARY_TSV"
else
    echo "Array task ${SLURM_ARRAY_TASK_ID} done."
    echo "Per-sample summary rows in: ${SUMMARY_ROWS_DIR}/"
    echo "After all array tasks, merge rows with:"
    echo "echo -e \"sample_id\\ttotal_reads\\tmapped_reads\\tmapped_percent\\tavg_depth_q20\" > ${SUMMARY_TSV}"
    echo "cat ${SUMMARY_ROWS_DIR}/*.tsv >> ${SUMMARY_TSV}"

    # Submit summary merge once (from array task 1)
    if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
        MERGE_JOB_ID="$(sbatch --dependency=afterok:${SLURM_ARRAY_JOB_ID} --wrap "echo -e 'sample_id\ttotal_reads\tmapped_reads\tmapped_percent\tavg_depth_q20' > '${SUMMARY_TSV}'; cat '${SUMMARY_ROWS_DIR}'/*.tsv >> '${SUMMARY_TSV}'" | awk '{print $4}')"
        echo "Submitted merge job ${MERGE_JOB_ID}."
    fi
fi

# Step 1: Load sample IDs from sample sheet
# - Uses first column (sample_id) and ignores header
# - Builds paths to BAM files produced in previous alignment step
#
# Step 2: Validate BAM inputs and index if needed
# - Checks BAM existence and non-empty files
# - Runs samtools index when .bai is missing
#
# Step 3: Compute mapping statistics with samtools flagstat
# - Produces per-sample flagstat reports
# - Extracts total reads, mapped reads, and mapped percentage
#
# Step 4: Estimate average depth with samtools depth
# - Uses -Q 20 to focus on reads with mapping quality >= 20
# - Stores per-sample mean depth and appends global summary table
#
# Step 5: Run Qualimap bamqc for detailed QC
# - Generates per-sample QC reports (coverage, GC bias, mapping profiles)
# - Writes logs per sample for troubleshooting
#
# Step 6: Write final summary outputs
# - In array mode, writes one summary row per sample in summary_rows/
# - In sequential mode, also builds mapping_coverage_summary.tsv directly
#
# Parallelization note:
# - Current setting is 5 samples per array task (SAMPLES_PER_TASK=5)
# - To scale to full dataset, change SAMPLES_PER_TASK and --array range together

#!/bin/bash
#SBATCH -p workq
#SBATCH --time=95:00:00
#SBATCH -J 02-fastp_cleaning
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.out
#SBATCH --array=1-18
#SBATCH --mem=32G 
#SBATCH --cpus-per-task=4  
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u
set -o pipefail

# Clean environment
module purge

# Load fastp
module load bioinfo/fastp/0.23.4

# Define file paths
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_complete.csv"
OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/02-filtering_indexing/workflow_full/cleaned_reads"
FAILED_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/02-filtering_indexing/workflow_full/failed_reads"
LOG_DIR="${OUTPUT_DIR}/logs"
SAMPLES_PER_TASK=10

# Create output directories
mkdir -p "$OUTPUT_DIR" "$FAILED_DIR" "$LOG_DIR"

echo "Starting fastp read filtering based on sample sheet: $SAMPLE_SHEET"

# Read the sample sheet (skip header line)
# Expect read1/read2 as the last two columns in the CSV
FILES_LIST="${OUTPUT_DIR}/fastp_files.list"
tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1 "\t" $(NF-1) "\t" $NF}' > "$FILES_LIST"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

sed -n "${START_LINE},${END_LINE}p" "$FILES_LIST" | while IFS=$'\t' read -r sample_id read1 read2; do
    [ -z "$sample_id" ] && continue
    echo "Processing sample: $sample_id"

    if [ ! -s "$read1" ] || [ ! -s "$read2" ]; then
        echo "WARNING: missing input FASTQ for ${sample_id}" >&2
        continue
    fi

    # Define output file names
    out_r1="${OUTPUT_DIR}/${sample_id}_R1.trimmed.fastq.gz"
    out_r2="${OUTPUT_DIR}/${sample_id}_R2.trimmed.fastq.gz"
    failed_out="${FAILED_DIR}/${sample_id}.failed.fastq.gz"
    html_report="${LOG_DIR}/${sample_id}.fastp.html"
    json_report="${LOG_DIR}/${sample_id}.fastp.json"
    log_file="${LOG_DIR}/${sample_id}.fastp.log"

    if [ -s "$out_r1" ] && [ -s "$out_r2" ] && [ -s "$json_report" ]; then
        echo "SKIP: outputs already exist for $sample_id"
        continue
    fi

    # Run fastp (see previous detailed block for full options and comments)
    fastp --in1 "$read1" --in2 "$read2" \
          --out1 "$out_r1" --out2 "$out_r2" \
          --failed_out "$failed_out" \
          --detect_adapter_for_pe \
          --thread "${SLURM_CPUS_PER_TASK:-4}" \
          --trim_poly_g \
          --cut_front --cut_tail \
          --cut_window_size 4 \
          --cut_mean_quality 20 \
          --length_required 100 \
          --qualified_quality_phred 20 \
          --unqualified_percent_limit 30 \
          --dup_calc_accuracy 2 \
          --html "$html_report" \
          --json "$json_report" &> "$log_file"

    echo "Finished processing $sample_id"
done

echo "Task ${TASK_ID} processed sample lines ${START_LINE}-${END_LINE}."


   # Run fastp with filtering and quality trimming
    # ----------------------------------------------
    # --in1 / --in2           : Input paired-end FASTQ files
    # --out1 / --out2         : Output cleaned FASTQ files
    # --failed_out            : Reads that fail quality filters
    # --detect_adapter_for_pe: Automatically detect adapter sequences (useful for paired-end)
    # --thread 4              : Use 4 CPU threads
    # --trim_poly_g           : Remove long poly-G tails (common NovaSeq artifacts)
    # --length_required 100   : Keep only reads ≥ 100 bp after trimming
    # --dup_calc_accuracy 2   : Increase accuracy of duplication level estimation
    # --html / --json         : Generate quality reports (interactive + machine-readable)
    # &> log_file             : Redirect all output and errors to a log file

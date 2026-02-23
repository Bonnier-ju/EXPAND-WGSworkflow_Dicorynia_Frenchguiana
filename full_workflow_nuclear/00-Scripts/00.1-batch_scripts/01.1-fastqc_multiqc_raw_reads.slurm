#!/bin/bash
#SBATCH -p workq                             
#SBATCH --time=24:00:00                      # Maximum wall time (hh:mm:ss)
#SBATCH -J 01-fastqc_multiqc_raw             # Job name
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out   # Standard output log
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err   # Standard error log
#SBATCH --mem=24G                            
#SBATCH --cpus-per-task=16                  # Number of CPUs to use
#SBATCH --mail-user=julien.bonnier@inrae.fr # Email address for notifications
#SBATCH --mail-type=END,FAIL                # When to send email notifications

# -------------------------------------------------------------------
#  DESCRIPTION
# -------------------------------------------------------------------
# FastQC evaluates raw reads quality (per-base quality, adapters, GC content,
# duplication levels, etc.) on each FASTQ file. MultiQC aggregates all FastQC
# reports into a single summary to quickly spot outliers across samples.

# -------------------------------------------------------------------
#  STEP 0: Environment setup
# -------------------------------------------------------------------
# Remove all previously loaded modules
module purge

# Load required bioinformatics tools
module load bioinfo/FastQC/0.12.1
module load bioinfo/MultiQC/1.14

# -------------------------------------------------------------------
#  STEP 1: Define inputs and output directories
# -------------------------------------------------------------------
FASTQC_OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/01-quality_control/workflow_full/fastqc"
MULTIQC_OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/01-quality_control/workflow_full/multiqc"
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_complete.csv"
FASTQC_THREADS_PER_FILE=1

# Create output directories if they do not exist
mkdir -p "$FASTQC_OUTPUT_DIR" "$MULTIQC_OUTPUT_DIR"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found or empty: $SAMPLE_SHEET" >&2
    exit 1
fi

echo "Starting FastQC using sample sheet: $SAMPLE_SHEET"
echo "Results will be saved to: $FASTQC_OUTPUT_DIR"

# -------------------------------------------------------------------
#  STEP 2: Run FastQC on files listed in the sample sheet
# -------------------------------------------------------------------
# The sample sheet is a CSV with headers; read1/read2 are the last two columns
# 2) Check file existence
# 3) Avoid duplicates
# 4) Log processed files
FILES_LIST="${FASTQC_OUTPUT_DIR}/fastqc_files.list"
tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $(NF-1); print $NF}' | sort -u > "$FILES_LIST"

PARALLEL_JOBS=$(( ${SLURM_CPUS_PER_TASK:-8} / FASTQC_THREADS_PER_FILE ))
[ "$PARALLEL_JOBS" -lt 1 ] && PARALLEL_JOBS=1

while read -r fq_file; do
    [ -z "$fq_file" ] && continue
    if [ ! -f "$fq_file" ]; then
        echo "WARNING: missing file: $fq_file" >&2
        continue
    fi
    echo "$fq_file"
done < "$FILES_LIST" | xargs -r -P "$PARALLEL_JOBS" -I {} bash -c 'echo "Processing file: {}"; fastqc -t '"$FASTQC_THREADS_PER_FILE"' -o "'"$FASTQC_OUTPUT_DIR"'" "{}"'

echo "FastQC completed on all samples."

# -------------------------------------------------------------------
#  STEP 3: Aggregate FastQC results using MultiQC
# -------------------------------------------------------------------
echo "Generating MultiQC summary report..."
multiqc "$FASTQC_OUTPUT_DIR" -o "$MULTIQC_OUTPUT_DIR"

echo "MultiQC report created in: $MULTIQC_OUTPUT_DIR"

# -------------------------------------------------------------------
#  DONE
# -------------------------------------------------------------------
echo "Quality control workflow complete"

#!/bin/bash
#SBATCH -p workq                             
#SBATCH --time=04:00:00                      # Maximum wall time (hh:mm:ss)
#SBATCH -J 01-fastqc_multiqc_raw             # Job name
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out   # Standard output log
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.out   # Standard error log
#SBATCH --mem=16G                            
#SBATCH --cpus-per-task=4                    # Number of CPUs to use
#SBATCH --mail-user=julien.bonnier@inrae.fr # Email address for notifications
#SBATCH --mail-type=END,FAIL                 # When to send email notifications

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
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_test.csv"
FASTQC_OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/01-quality_control/workflow_test/fastqc_results"
MULTIQC_OUTPUT_DIR="${FASTQC_OUTPUT_DIR}/multiqc"

# Create output directories if they do not exist
mkdir -p "$MULTIQC_OUTPUT_DIR"

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

while read fq_file; do
    if [ -z "$fq_file" ]; then
        continue
    fi
    if [ ! -f "$fq_file" ]; then
        echo "WARNING: missing file: $fq_file" >&2
        continue
    fi
    echo "Processing file: $fq_file"
    fastqc -t "${SLURM_CPUS_PER_TASK:-4}" -o "$FASTQC_OUTPUT_DIR" "$fq_file"
done < "$FILES_LIST"

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

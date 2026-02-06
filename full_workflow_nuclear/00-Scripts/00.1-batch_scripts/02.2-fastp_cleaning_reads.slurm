#!/bin/bash
#SBATCH -p workq
#SBATCH --time=06:00:00
#SBATCH -J 02-fastp_cleaning
#SBATCH -o /work/user/jbonnier/local_adaptation_dicorynia/full_workflow/00-batch_scripts/output/output-%x-%j.out
#SBATCH -e /work/user/jbonnier/local_adaptation_dicorynia/full_workflow/00-batch_scripts/error/error-%x-%j.out
#SBATCH --mem=32G 
#SBATCH --cpus-per-task=4  
#SBATCH --mail-user=julien.bonnier@ecofog.gf
#SBATCH --mail-type=END,FAIL

# Clean environment
module purge

# Load fastp
module load bioinfo/fastp/0.23.4

# Define file paths
SAMPLE_SHEET="/home/jbonnier/work/local_adaptation_dicorynia/full_workflow/metadata/sample_sheet.csv"
OUTPUT_DIR="/work/user/jbonnier/local_adaptation_dicorynia/full_workflow/03-cleaned_reads"
FAILED_DIR="${OUTPUT_DIR}/failed_reads"
LOG_DIR="${OUTPUT_DIR}/logs"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$FAILED_DIR" "$LOG_DIR"

echo "Starting fastp read filtering based on sample sheet: $SAMPLE_SHEET"

# Read the sample sheet (skip header line)
tail -n +2 "$SAMPLE_SHEET" | while IFS=',' read -r sample_id species read1 read2 reference_genome; do
    echo "🔍 Processing sample: $sample_id ($species)"

    # Define output file names
    out_r1="${OUTPUT_DIR}/${sample_id}_R1.trimmed.fastq.gz"
    out_r2="${OUTPUT_DIR}/${sample_id}_R2.trimmed.fastq.gz"
    failed_out="${FAILED_DIR}/${sample_id}.failed.fastq.gz"
    html_report="${LOG_DIR}/${sample_id}.fastp.html"
    json_report="${LOG_DIR}/${sample_id}.fastp.json"
    log_file="${LOG_DIR}/${sample_id}.fastp.log"

    # Run fastp (see previous detailed block for full options and comments)
    fastp --in1 "$read1" --in2 "$read2" \
          --out1 "$out_r1" --out2 "$out_r2" \
          --failed_out "$failed_out" \
          --detect_adapter_for_pe \
          --thread 4 \
          --trim_poly_g \
          --length_required 100 \
          --dup_calc_accuracy 2 \
          --html "$html_report" \
          --json "$json_report" &> "$log_file"

    echo "✅ Finished processing $sample_id"
done

echo "🎉 All samples processed with fastp."


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
#!/bin/bash
#SBATCH -p workq
#SBATCH --time=30:00:00
#SBATCH -J 05-coverage_analysis
#SBATCH -o /work/user/jbonnier/local_adaptation_dicorynia/full_workflow/00-batch_scripts/output/output-%x-%j.out
#SBATCH -e /work/user/jbonnier/local_adaptation_dicorynia/full_workflow/00-batch_scripts/error/error-%x-%j.out
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@ecofog.gf
#SBATCH --mail-type=END,FAIL

# Clean environment
module purge

# Load required modules
module load bioinfo/samtools/1.19
module load bioinfo/Qualimap/31-08-20
module load statistics/R/4.2.2

# Set directories
SAMPLE_SHEET="/home/jbonnier/work/local_adaptation_dicorynia/full_workflow/metadata/sample_sheet.csv"
ALIGN_DIR="/work/user/jbonnier/local_adaptation_dicorynia/full_workflow/04-alignments"
OUTPUT_DIR="/work/user/jbonnier/local_adaptation_dicorynia/full_workflow/05-coverage"
mkdir -p "$OUTPUT_DIR"

# Define summary file
COVERAGE_SUMMARY_FILE="${OUTPUT_DIR}/coverage_summary.tsv"
echo -e "Sample\tAverage_Coverage" > "$COVERAGE_SUMMARY_FILE"

# Initialize counters
TOTAL_COVERAGE=0
SAMPLE_COUNT=0

# Loop through sample sheet
# Columns: sample_id,species,read1,read2,reference
tail -n +2 "$SAMPLE_SHEET" | while IFS=',' read -r sample_id species read1 read2 reference; do
    BAMFILE="${ALIGN_DIR}/${sample_id}.sorted.bam"

    if [[ -f "$BAMFILE" ]]; then
        echo "Processing $sample_id..."

        # Compute mean coverage using quality filter
        AVG_COV=$(samtools depth -a -Q 20 "$BAMFILE" | awk '{sum+=$3} END { if (NR>0) print sum/NR; else print 0}')
        echo -e "${sample_id}\t${AVG_COV}" >> "$COVERAGE_SUMMARY_FILE"

        # Accumulate values
        TOTAL_COVERAGE=$(echo "$TOTAL_COVERAGE + $AVG_COV" | bc)
        SAMPLE_COUNT=$((SAMPLE_COUNT + 1))

        # Run Qualimap
        SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${sample_id}"
        mkdir -p "$SAMPLE_OUTPUT_DIR"
        qualimap bamqc -bam "$BAMFILE" -outdir "${SAMPLE_OUTPUT_DIR}/${sample_id}_qualimap_report"

        echo "Analysis complete for $sample_id"
    else
        echo "[Warning] BAM file not found for $sample_id"
    fi

done

# Compute overall average coverage
if [ $SAMPLE_COUNT -gt 0 ]; then
    OVERALL_AVG_COV=$(echo "scale=2; $TOTAL_COVERAGE / $SAMPLE_COUNT" | bc)
    echo -e "Overall Average Coverage:\t$OVERALL_AVG_COV" >> "$COVERAGE_SUMMARY_FILE"
    echo "Overall average coverage calculated and saved to $COVERAGE_SUMMARY_FILE."
else
    echo "No BAM files processed. No average calculated."
fi

# Call R script to visualize results
Rscript /home/jbonnier/work/local_adaptation_dicorynia/full_workflow/00-R_scripts/05.1-coverage_summary_plot.R "$COVERAGE_SUMMARY_FILE" "$OUTPUT_DIR"
echo "Coverage analysis and visualization completed."
echo "All tasks completed successfully."

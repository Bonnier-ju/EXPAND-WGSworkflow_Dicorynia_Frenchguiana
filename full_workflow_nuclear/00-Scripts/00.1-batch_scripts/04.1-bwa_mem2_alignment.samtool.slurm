#!/bin/bash
#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 03-bwa_mem2_alignment
#SBATCH -o /work/user/jbonnier/local_adaptation_dicorynia/full_workflow/00-batch_scripts/output/output-%x-%j.out
#SBATCH -e /work/user/jbonnier/local_adaptation_dicorynia/full_workflow/00-batch_scripts/error/error-%x-%j.out
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@ecofog.gf
#SBATCH --mail-type=END,FAIL

# Clean environment
module purge

# Load bwa-mem2 and samtools
module load bioinfo/bwa-mem2/2.2.1
module load bioinfo/samtools/1.21

# File paths
SAMPLE_SHEET="/home/jbonnier/work/local_adaptation_dicorynia/full_workflow/metadata/sample_sheet.csv"
READS_DIR="/work/user/jbonnier/local_adaptation_dicorynia/full_workflow/03-cleaned_reads"
ALIGN_DIR="/work/user/jbonnier/local_adaptation_dicorynia/full_workflow/04-alignments"
LOG_DIR="${ALIGN_DIR}/logs"
FLAGSTAT_DIR="${ALIGN_DIR}/flagstat"
SUMMARY_TSV="${FLAGSTAT_DIR}/flagstat_summary.tsv"

mkdir -p "$ALIGN_DIR" "$LOG_DIR" "$FLAGSTAT_DIR"

# Write header to summary file
echo -e "sample_id\ttotal_reads\tmapped_reads\tpercent_mapped" > "$SUMMARY_TSV"

# Loop over sample sheet (skip header)
tail -n +2 "$SAMPLE_SHEET" | while IFS=',' read -r sample_id species read1 read2 reference; do
    echo "Aligning sample: $sample_id $species"

    # Define input files
    R1="${READS_DIR}/${sample_id}_R1.trimmed.fastq.gz"
    R2="${READS_DIR}/${sample_id}_R2.trimmed.fastq.gz"
    BAM="${ALIGN_DIR}/${sample_id}.sorted.bam"
    LOG="${LOG_DIR}/${sample_id}.alignment.log"
    FLAGSTAT_OUT="${FLAGSTAT_DIR}/${sample_id}.flagstat.txt"

    # Alignment pipeline for paired-end reads using bwa-mem2 and samtools
    # --------------------------------------------------------------------

    # Extract library (LB) and platform unit (PU) from FASTQ filename
    filename=$(basename "$read1")
    LB=$(echo "$filename" | cut -d'_' -f1,2)
    LANE=$(echo "$filename" | cut -d'_' -f3)
    FLOWCELL=$(echo "$filename" | cut -d'_' -f5 | cut -d'.' -f1)
    PU="${FLOWCELL}.${LANE}"
    PLATFORM="ILLUMINA"

    # Construct read group header
    RG="@RG\tID:${sample_id}\tSM:${sample_id}\tLB:${LB}\tPU:${PU}\tPL:${PLATFORM}"

    #Align reads and convert to sorted BAM
    bwa-mem2 mem -M -R "$RG" -t 8 "$reference" "$R1" "$R2" |
    samtools view -b -@ 4 - |
    samtools sort -@ 4 -o "$BAM" - &&
    samtools index "$BAM"

    #Run samtools flagstat to summarize mapping statistics
    samtools flagstat "$BAM" > "$FLAGSTAT_OUT"

    # Extract stats for summary
    total_reads=$(grep "in total" "$FLAGSTAT_OUT" | awk '{print $1}')
    mapped_reads=$(grep " mapped (" "$FLAGSTAT_OUT" | head -n1 | awk '{print $1}')
    percent_mapped=$(grep " mapped (" "$FLAGSTAT_OUT" | head -n1 | awk -F '[()]' '{print $2}' | cut -d'%' -f1)

    echo -e "$sample_id\t$total_reads\t$mapped_reads\t$percent_mapped" >> "$SUMMARY_TSV"

    echo "Alignment and flagstat completed for $sample_id" &>> "$LOG"
done

echo "All samples aligned, indexed, and checked with samtools flagstat. Summary saved in $SUMMARY_TSV."

# Step 1: Construct read group header
    # - Required by tools like GATK for downstream processing
# Step 2: Align reads with bwa-mem2
    # - -t ${SLURM_CPUS_PER_TASK} : use all allocated CPUs
    # - -M : mark shorter split hits as secondary (GATK compatibility)
    # - -R : include read group information in the header
# Step 3: Convert SAM to BAM
    # - -b : output BAM format
    # - -@ 4 : use 4 threads for conversion
# Step 4: Sort BAM by coordinates
    # - -@ 4 : use 4 threads for sorting
    # - -o : output sorted BAM file
    # - && : ensure that the next command runs only if the previous one succeeds
# Step 5: Index sorted BAM 
    # - Indexing is necessary for random access to the BAM file
# Step 6: Generate flagstat summary
    # - Provides alignment quality metrics (e.g., % mapped, duplicates, total reads)
    # - Summary file with key metrics created for all samples

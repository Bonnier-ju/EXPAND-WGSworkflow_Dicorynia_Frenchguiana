#!/bin/bash
#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 04-bwa_mem2_alignment
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.out
#SBATCH --array=1-2
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@ecofog.gf
#SBATCH --mail-type=END,FAIL

# Clean environment
module purge

# Load bwa-mem2 and samtools
module load bioinfo/bwa-mem2/2.2.1
module load bioinfo/samtools/1.21

# File paths
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_test.csv"
READS_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/02-filtering_indexing/workflow_test/cleaned_reads"
ALIGN_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/03-reading_sorting/workflow_test/alignments"
LOG_DIR="${ALIGN_DIR}/logs"
mkdir -p "$ALIGN_DIR" "$LOG_DIR"

FILES_LIST="${ALIGN_DIR}/alignment_files.list"
tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1 "\t" $(NF-2) "\t" $(NF-1) "\t" $NF}' > "$FILES_LIST"

SAMPLES_PER_TASK=3
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

# Loop over sample sheet (skip header)
# Expect read1/read2 as the last two columns in the CSV, and reference genome as the last column
sed -n "${START_LINE},${END_LINE}p" "$FILES_LIST" | while IFS=$'\t' read -r sample_id reference read1 read2; do
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
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "ERROR: missing reads for $sample_id: $R1 $R2" | tee -a "$LOG"
        continue
    fi
    if [ -s "$BAM" ]; then
        echo "SKIP: BAM already exists for $sample_id: $BAM" | tee -a "$LOG"
        continue
    fi

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
    bwa-mem2 mem -M -R "$RG" -t "${SLURM_CPUS_PER_TASK:-8}" "$reference" "$R1" "$R2" |
    samtools view -b -@ "${SLURM_CPUS_PER_TASK:-8}" - |
    samtools sort -@ "${SLURM_CPUS_PER_TASK:-8}" -o "$BAM" - &&
    samtools index "$BAM"

    if [ -s "$BAM" ]; then
        echo "Alignment completed for $sample_id" &>> "$LOG"
    else
        echo "ERROR: BAM not created for $sample_id" | tee -a "$LOG"
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

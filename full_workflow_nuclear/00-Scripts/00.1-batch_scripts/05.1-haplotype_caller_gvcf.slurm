#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=48:00:00
#SBATCH -J 05-haplotype_caller
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.out
#SBATCH --array=1-6
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u
set -o pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_test.csv"
ALIGN_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/03-reading_sorting/workflow_test/alignments"
OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/05-variant_calling/workflow_test/haplotypecaller_gvcf"
LOG_DIR="${OUTPUT_DIR}/logs"
SAMPLES_LIST="${OUTPUT_DIR}/samples.list"

SAMPLES_PER_TASK=1
JAVA_MEM_GB=20
TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-4}"
PAIRHMM_THREADS=$((TOTAL_THREADS - 1))
[ "$PAIRHMM_THREADS" -lt 1 ] && PAIRHMM_THREADS=1
TMP_ROOT="/tmp/${USER}/gatk_hc_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
mkdir -p "$TMP_ROOT"

if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found: $SAMPLE_SHEET" >&2
    exit 1
fi

# Build sample list once: sample_id + reference genome path (col 1 and 7)
if [ ! -s "$SAMPLES_LIST" ]; then
    tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1 "\t" $7}' > "$SAMPLES_LIST"
fi

# -------------------------------------------------------------------
# STEP 2: Select samples for current array task
# -------------------------------------------------------------------
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

sed -n "${START_LINE},${END_LINE}p" "$SAMPLES_LIST" | while IFS=$'\t' read -r sample_id reference; do
    [ -z "$sample_id" ] && continue

    BAM="${ALIGN_DIR}/${sample_id}.sorted.bam"
    BAI="${BAM}.bai"
    GVCF="${OUTPUT_DIR}/${sample_id}.g.vcf.gz"
    GVCF_INDEX="${GVCF}.tbi"
    SAMPLE_LOG="${LOG_DIR}/${sample_id}.haplotypecaller.log"

    echo "START $(date -Iseconds) sample=${sample_id}" > "$SAMPLE_LOG"

    # 2.1 Validate core inputs
    if [ ! -s "$reference" ]; then
        echo "ERROR: missing reference genome: $reference" | tee -a "$SAMPLE_LOG"
        continue
    fi
    if [ ! -s "$BAM" ]; then
        echo "ERROR: missing BAM: $BAM" | tee -a "$SAMPLE_LOG"
        continue
    fi

    # 2.2 Ensure reference indexes exist
    if [ ! -s "${reference}.fai" ]; then
        echo "INFO: creating reference .fai index" >> "$SAMPLE_LOG"
        samtools faidx "$reference" >> "$SAMPLE_LOG" 2>&1
    fi
    REF_DICT="${reference%.fa}.dict"
    if [ ! -s "$REF_DICT" ]; then
        echo "INFO: creating reference sequence dictionary" >> "$SAMPLE_LOG"
        gatk CreateSequenceDictionary -R "$reference" -O "$REF_DICT" >> "$SAMPLE_LOG" 2>&1
    fi

    # 2.3 Ensure BAM index exists
    if [ ! -s "$BAI" ]; then
        echo "INFO: creating BAM index" >> "$SAMPLE_LOG"
        samtools index "$BAM" >> "$SAMPLE_LOG" 2>&1
    fi

    # 2.4 Skip if output already present
    if [ -s "$GVCF" ] && [ -s "$GVCF_INDEX" ]; then
        echo "SKIP: output already exists: $GVCF" | tee -a "$SAMPLE_LOG"
        continue
    fi

    # 2.5 Run HaplotypeCaller in GVCF mode
    gatk --java-options "-Xms${JAVA_MEM_GB}G -Xmx${JAVA_MEM_GB}G -XX:ParallelGCThreads=2" HaplotypeCaller \
        -R "$reference" \
        -I "$BAM" \
        -O "$GVCF" \
        -ERC GVCF \
        --native-pair-hmm-threads "$PAIRHMM_THREADS" \
        --create-output-variant-index true \
        --tmp-dir "$TMP_ROOT" >> "$SAMPLE_LOG" 2>&1

    if [ -s "$GVCF" ] && [ -s "$GVCF_INDEX" ]; then
        echo "DONE $(date -Iseconds) sample=${sample_id} gvcf=${GVCF}" >> "$SAMPLE_LOG"
    else
        echo "ERROR: GVCF output missing for ${sample_id}" | tee -a "$SAMPLE_LOG"
    fi
done

rm -rf "$TMP_ROOT"

echo "HaplotypeCaller array task ${TASK_ID} completed."

# Step 1: Build sample list from sample sheet
# - Uses sample_id and reference_genome columns
#
# Step 2: Validate BAM/reference inputs and indexes
# - Checks BAM and reference existence
# - Creates missing .bai/.fai/.dict indexes
#
# Step 3: Run GATK HaplotypeCaller
# - One sample per array task
# - Emits per-sample GVCF output
#
# Step 4: Validate outputs
# - Ensures both .g.vcf.gz and .g.vcf.gz.tbi exist

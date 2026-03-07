#!/bin/bash
#SBATCH -p workq
#SBATCH --time=48:00:00
#SBATCH -J cp-05-haplotype_caller
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.err
#SBATCH --array=1-18
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# Chloroplastic variant calling (per sample) with GATK HaplotypeCaller
# Input : cp-03 BAM + BAI
# Output: sample.g.vcf.gz + sample.g.vcf.gz.tbi
# -------------------------------------------------------------------

module purge
module load devel/python/Python-3.11.1
module load devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/samtools/1.21

SAMPLE_SHEET="${SAMPLE_SHEET:-/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete.csv}"
ALIGN_DIR="${ALIGN_DIR:-/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/03-alignment/workflow_full/alignments}"
OUTPUT_DIR="${OUTPUT_DIR:-/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/05-variant_calling/workflow_full/haplotypecaller_gvcf}"
LOG_DIR="${OUTPUT_DIR}/logs"

SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-10}"
JAVA_MEM_GB=20
SAMPLE_PLOIDY=1
MAX_READS_PER_ALIGNMENT_START=0
TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-4}"
PAIRHMM_THREADS=$((TOTAL_THREADS - 1))
[ "$PAIRHMM_THREADS" -lt 1 ] && PAIRHMM_THREADS=1
TMP_ROOT="/tmp/${USER}/cp_gatk_hc_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$TMP_ROOT"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found or empty: $SAMPLE_SHEET" >&2
    exit 1
fi

# Columns: sample_id (1), reference_genome (7)
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

# Build per-task sample list to avoid race conditions between array tasks
TASK_SAMPLES_LIST="${TMPDIR:-/tmp}/cp_hc_samples_${SLURM_JOB_ID:-manual}_${TASK_ID}.list"
awk -F, 'NR>1{print $1 "\t" $7}' "$SAMPLE_SHEET" | sed -n "${START_LINE},${END_LINE}p" > "$TASK_SAMPLES_LIST"
if [ ! -s "$TASK_SAMPLES_LIST" ]; then
    echo "WARNING: no samples selected for task ${TASK_ID} (lines ${START_LINE}-${END_LINE})."
fi

while IFS=$'\t' read -r sample_id reference; do
    [ -z "$sample_id" ] && continue

    BAM="${ALIGN_DIR}/${sample_id}.sorted.bam"
    BAI="${BAM}.bai"
    GVCF="${OUTPUT_DIR}/${sample_id}.g.vcf.gz"
    GVCF_INDEX="${GVCF}.tbi"
    SAMPLE_LOG="${LOG_DIR}/${sample_id}.haplotypecaller.log"

    echo "START $(date -Iseconds) sample=${sample_id}" > "$SAMPLE_LOG"

    if [ ! -s "$reference" ]; then
        echo "ERROR: missing reference genome: $reference" | tee -a "$SAMPLE_LOG"
        continue
    fi
    if [ ! -s "$BAM" ]; then
        echo "ERROR: missing BAM: $BAM" | tee -a "$SAMPLE_LOG"
        continue
    fi

    if [ ! -s "${reference}.fai" ]; then
        echo "INFO: creating reference fai index" >> "$SAMPLE_LOG"
        samtools faidx "$reference" >> "$SAMPLE_LOG" 2>&1
    fi
    if [[ "$reference" == *.fasta ]]; then
        REF_DICT="${reference%.fasta}.dict"
    elif [[ "$reference" == *.fa ]]; then
        REF_DICT="${reference%.fa}.dict"
    else
        REF_DICT="${reference}.dict"
    fi
    if [ ! -s "$REF_DICT" ]; then
        echo "INFO: creating reference sequence dictionary" >> "$SAMPLE_LOG"
        # Array tasks can race here; if another task creates the dict first, continue.
        gatk CreateSequenceDictionary -R "$reference" -O "$REF_DICT" >> "$SAMPLE_LOG" 2>&1 || {
            if [ ! -s "$REF_DICT" ]; then
                echo "ERROR: failed to create reference dictionary: $REF_DICT" | tee -a "$SAMPLE_LOG"
                continue
            fi
        }
    fi

    if [ ! -s "$BAI" ]; then
        echo "INFO: creating BAM index" >> "$SAMPLE_LOG"
        samtools index "$BAM" >> "$SAMPLE_LOG" 2>&1
    fi

    if [ -s "$GVCF" ] && [ -s "$GVCF_INDEX" ]; then
        echo "SKIP: output already exists: $GVCF" | tee -a "$SAMPLE_LOG"
        continue
    fi

    gatk --java-options "-Xms${JAVA_MEM_GB}G -Xmx${JAVA_MEM_GB}G -XX:ParallelGCThreads=2" HaplotypeCaller \
        -R "$reference" \
        -I "$BAM" \
        -O "$GVCF" \
        -ERC GVCF \
        --sample-ploidy "$SAMPLE_PLOIDY" \
        --max-reads-per-alignment-start "$MAX_READS_PER_ALIGNMENT_START" \
        --native-pair-hmm-threads "$PAIRHMM_THREADS" \
        --create-output-variant-index true \
        --tmp-dir "$TMP_ROOT" >> "$SAMPLE_LOG" 2>&1

    if [ -s "$GVCF" ] && [ -s "$GVCF_INDEX" ]; then
        echo "DONE $(date -Iseconds) sample=${sample_id} gvcf=${GVCF}" >> "$SAMPLE_LOG"
    else
        echo "ERROR: GVCF output missing for ${sample_id}" | tee -a "$SAMPLE_LOG"
    fi
done < "$TASK_SAMPLES_LIST"

rm -f "$TASK_SAMPLES_LIST"
rm -rf "$TMP_ROOT"

echo "HaplotypeCaller cp array task ${TASK_ID} completed."

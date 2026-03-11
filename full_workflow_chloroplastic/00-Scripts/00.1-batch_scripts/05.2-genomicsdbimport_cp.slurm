#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J cp-05.2-genomicsdbimport
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load devel/python/Python-3.11.1
module load devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/samtools/1.21

SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete_filtered.csv"
GVCF_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/05-variant_calling/workflow_full/haplotypecaller_gvcf"
AGGREGATION_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/05-variant_calling/workflow_full/genomicsdbimport"

REF_GENOME="$(awk -F, 'NR==2{print $7}' "$SAMPLE_SHEET")"
SAMPLE_MAP="${AGGREGATION_DIR}/sample_map.txt"
DB_PATH="${AGGREGATION_DIR}/GATK_database_dicorynia_cp"
SEQ_LIST="${AGGREGATION_DIR}/reference.sequences.list"
GVCF_LIST="${AGGREGATION_DIR}/gvcf_files.list"
LOG_FILE="${AGGREGATION_DIR}/genomicsdbimport_cp.log"

JAVA_MEM_GB=24
READER_THREADS=8
BATCH_SIZE=20
TMP_DIR="/tmp/${USER}/cp_genomicsdbimport_${SLURM_JOB_ID}"

mkdir -p "$AGGREGATION_DIR" "$TMP_DIR"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found or empty: $SAMPLE_SHEET" >&2
    exit 1
fi
if [ -z "$REF_GENOME" ] || [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found from sample sheet: $REF_GENOME" >&2
    exit 1
fi
if [ ! -d "$GVCF_DIR" ]; then
    echo "ERROR: GVCF directory not found: $GVCF_DIR" >&2
    exit 1
fi
if [ -e "$DB_PATH" ]; then
    echo "ERROR: GenomicsDB workspace already exists: $DB_PATH" >&2
    echo "Remove it manually for a clean rebuild." >&2
    exit 1
fi

echo "START $(date -Iseconds) GenomicsDBImport cp" > "$LOG_FILE"

if [ ! -s "${REF_GENOME}.fai" ]; then
    echo "INFO: creating reference fai index" | tee -a "$LOG_FILE"
    samtools faidx "$REF_GENOME" >> "$LOG_FILE" 2>&1
fi
cut -f1 "${REF_GENOME}.fai" > "$SEQ_LIST"
echo "INFO: sequence list created: $SEQ_LIST" >> "$LOG_FILE"

# Build GVCF list from filtered sample sheet only (164 samples)
: > "$GVCF_LIST"
while IFS=',' read -r sample_id _rest; do
    gvcf="${GVCF_DIR}/${sample_id}.g.vcf.gz"
    if [ -s "$gvcf" ]; then
        echo "$gvcf" >> "$GVCF_LIST"
    else
        echo "WARNING: GVCF not found for sample ${sample_id}: ${gvcf}" | tee -a "$LOG_FILE"
    fi
done < <(tail -n +2 "$SAMPLE_SHEET")

if [ ! -s "$GVCF_LIST" ]; then
    echo "ERROR: no .g.vcf.gz files found from sample sheet samples" | tee -a "$LOG_FILE"
    exit 1
fi

: > "$SAMPLE_MAP"
while read -r gvcf; do
    sample_name="$(basename "$gvcf" .g.vcf.gz)"

    if ! gzip -t "$gvcf" 2>/dev/null; then
        echo "WARNING: corrupted gvcf skipped: $gvcf" | tee -a "$LOG_FILE"
        continue
    fi

    if [ ! -s "${gvcf}.tbi" ]; then
        echo "INFO: missing index for $gvcf, creating .tbi" >> "$LOG_FILE"
        gatk IndexFeatureFile -I "$gvcf" >> "$LOG_FILE" 2>&1
    fi

    echo -e "${sample_name}\t${gvcf}" >> "$SAMPLE_MAP"
done < "$GVCF_LIST"

if [ ! -s "$SAMPLE_MAP" ]; then
    echo "ERROR: no valid gvcf available after integrity checks" | tee -a "$LOG_FILE"
    exit 1
fi

echo "INFO: sample map created: $SAMPLE_MAP" >> "$LOG_FILE"
echo "INFO: number of samples: $(wc -l < "$SAMPLE_MAP")" >> "$LOG_FILE"

gatk --java-options "-Xms${JAVA_MEM_GB}G -Xmx${JAVA_MEM_GB}G -XX:ParallelGCThreads=2" GenomicsDBImport \
    --genomicsdb-workspace-path "$DB_PATH" \
    --batch-size "$BATCH_SIZE" \
    -L "$SEQ_LIST" \
    --sample-name-map "$SAMPLE_MAP" \
    --tmp-dir "$TMP_DIR" \
    --reader-threads "$READER_THREADS" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) GenomicsDBImport cp workspace=${DB_PATH}" | tee -a "$LOG_FILE"

rm -rf "$TMP_DIR"

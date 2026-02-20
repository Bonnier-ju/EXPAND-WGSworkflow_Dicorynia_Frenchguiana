#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J 05.2-genomicsdbimport
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.out
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
GVCF_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/05-variant_calling/workflow_test/haplotypecaller_gvcf"
AGGREGATION_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/05-variant_calling/workflow_test/genomicsdbimport"
REF_GENOME="/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa"

SAMPLE_MAP="${AGGREGATION_DIR}/sample_map.txt"
DB_PATH="${AGGREGATION_DIR}/GATK_database_dicorynia"
TMP_DIR="${AGGREGATION_DIR}/tmp"
SEQ_LIST="${AGGREGATION_DIR}/reference.sequences.list"
LOG_FILE="${AGGREGATION_DIR}/genomicsdbimport.log"

# GATK tuning
JAVA_MEM_GB=16
READER_THREADS=4
BATCH_SIZE=10

TMP_DIR="/tmp/${USER}/genomicsdbimport_${SLURM_JOB_ID}"
mkdir -p "$AGGREGATION_DIR" "$TMP_DIR"

if [ ! -d "$GVCF_DIR" ]; then
    echo "ERROR: GVCF directory not found: $GVCF_DIR" >&2
    exit 1
fi
if [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: Reference genome not found: $REF_GENOME" >&2
    exit 1
fi
if [ -e "$DB_PATH" ]; then
    echo "ERROR: GenomicsDB workspace already exists: $DB_PATH" >&2
    echo "Remove it manually if you want a clean rebuild." >&2
    exit 1
fi

echo "START $(date -Iseconds) GenomicsDBImport" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Prepare reference sequence list
# -------------------------------------------------------------------
if [ ! -s "${REF_GENOME}.fai" ]; then
    echo "INFO: creating reference fai index" | tee -a "$LOG_FILE"
    samtools faidx "$REF_GENOME" >> "$LOG_FILE" 2>&1
fi
cut -f1 "${REF_GENOME}.fai" > "$SEQ_LIST"
echo "INFO: sequence list created: $SEQ_LIST" >> "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Build sample map from *.g.vcf.gz files
# -------------------------------------------------------------------
find "$GVCF_DIR" -maxdepth 1 -type f -name "*.g.vcf.gz" | sort > "${AGGREGATION_DIR}/gvcf_files.list"

if [ ! -s "${AGGREGATION_DIR}/gvcf_files.list" ]; then
    echo "ERROR: no .g.vcf.gz files found in $GVCF_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

: > "$SAMPLE_MAP"
while read -r gvcf; do
    sample_name="$(basename "$gvcf" .g.vcf.gz)"
    if [ ! -s "${gvcf}.tbi" ]; then
        echo "INFO: missing index for $gvcf, creating .tbi" >> "$LOG_FILE"
        gatk IndexFeatureFile -I "$gvcf" >> "$LOG_FILE" 2>&1
    fi
    echo -e "${sample_name}\t${gvcf}" >> "$SAMPLE_MAP"
done < "${AGGREGATION_DIR}/gvcf_files.list"

echo "INFO: sample map created: $SAMPLE_MAP" >> "$LOG_FILE"
echo "INFO: number of samples: $(wc -l < "$SAMPLE_MAP")" >> "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: Run GenomicsDBImport
# -------------------------------------------------------------------
gatk --java-options "-Xms${JAVA_MEM_GB}G -Xmx${JAVA_MEM_GB}G -XX:ParallelGCThreads=2" GenomicsDBImport \
    --genomicsdb-workspace-path "$DB_PATH" \
    --batch-size "$BATCH_SIZE" \
    -L "$SEQ_LIST" \
    --sample-name-map "$SAMPLE_MAP" \
    --tmp-dir "$TMP_DIR" \
    --reader-threads "$READER_THREADS" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) GenomicsDBImport workspace=${DB_PATH}" | tee -a "$LOG_FILE"

# Step 1: Build sample map from all per-sample GVCFs
# - sample_name <tab> /path/to/sample.g.vcf.gz
#
# Step 2: Build interval list from reference .fai
# - One contig name per line
#
# Step 3: Aggregate all GVCFs with GenomicsDBImport
# - Outputs intermediate database at GATK_database_dicorynia

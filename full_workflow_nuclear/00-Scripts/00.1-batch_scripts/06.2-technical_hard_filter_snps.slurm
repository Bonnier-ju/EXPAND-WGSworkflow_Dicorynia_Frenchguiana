#!/bin/bash
#SBATCH -p workq
#SBATCH --time=12:00:00
#SBATCH -J 06.2-hard_filter_snps
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u
set -o pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load devel/python/Python-3.11.1
module load devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/bcftools/1.21
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
REF_GENOME="/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa"
INPUT_VCF="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.1-normalization_selection/joint.norm.biallelic.snps.vcf.gz"
INPUT_TBI="${INPUT_VCF}.tbi"

OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.2-technical_hard_filter"
LOG_DIR="${OUTPUT_DIR}/logs"
STATS_DIR="${OUTPUT_DIR}/stats"
TMP_DIR="${OUTPUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/06.2-technical_hard_filter_snps.log"

TAGGED_VCF="${OUTPUT_DIR}/snps.filtered.tagged.vcf.gz"
HARDFILTERED_VCF="${OUTPUT_DIR}/snps.hardfiltered.vcf.gz"
HARDFILTERED_TBI="${HARDFILTERED_VCF}.tbi"

RAW_STATS="${STATS_DIR}/joint.norm.biallelic.snps.stats.txt"
HARD_STATS="${STATS_DIR}/snps.hardfiltered.stats.txt"

MIN_DP="${MIN_DP:-5}"
MAX_DP="${MAX_DP:-100}"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$STATS_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 06.2 hard filter SNPs" > "$LOG_FILE"

if [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found: $REF_GENOME" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$INPUT_VCF" ] || [ ! -s "$INPUT_TBI" ]; then
    echo "ERROR: input SNP VCF missing: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "${REF_GENOME}.fai" ]; then
    samtools faidx "$REF_GENOME" >> "$LOG_FILE" 2>&1
fi
REF_DICT="${REF_GENOME%.fa}.dict"
if [ ! -s "$REF_DICT" ]; then
    gatk CreateSequenceDictionary -R "$REF_GENOME" -O "$REF_DICT" >> "$LOG_FILE" 2>&1
fi

# -------------------------------------------------------------------
# STEP 2: Tag low-quality variants with GATK VariantFiltration
# -------------------------------------------------------------------
gatk --java-options "-Xms8G -Xmx8G -XX:ParallelGCThreads=2" VariantFiltration \
    -R "$REF_GENOME" \
    -V "$INPUT_VCF" \
    -O "$TAGGED_VCF" \
    --filter-name "QD_lt_2" --filter-expression "QD < 2.0" \
    --filter-name "QUAL_lt_30" --filter-expression "QUAL < 30.0" \
    --filter-name "FS_gt_60" --filter-expression "FS > 60.0" \
    --filter-name "SOR_gt_3" --filter-expression "SOR > 3.0" \
    --filter-name "MQ_lt_40" --filter-expression "MQ < 40.0" \
    --filter-name "MQRankSum_lt_-12.5" --filter-expression "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum_lt_-8" --filter-expression "ReadPosRankSum < -8.0" \
    --filter-name "DP_out_of_range" --filter-expression "DP < ${MIN_DP} || DP > ${MAX_DP}" \
    --create-output-variant-index true \
    --tmp-dir "$TMP_DIR" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 3: Keep only PASS variants
# -------------------------------------------------------------------
gatk --java-options "-Xms8G -Xmx8G -XX:ParallelGCThreads=2" SelectVariants \
    -R "$REF_GENOME" \
    -V "$TAGGED_VCF" \
    --exclude-filtered true \
    -O "$HARDFILTERED_VCF" \
    --create-output-variant-index true \
    --tmp-dir "$TMP_DIR" >> "$LOG_FILE" 2>&1

if [ ! -s "$HARDFILTERED_VCF" ] || [ ! -s "$HARDFILTERED_TBI" ]; then
    echo "ERROR: hard-filtered SNP VCF missing" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 4: Build technical filtering statistics
# -------------------------------------------------------------------
bcftools stats "$INPUT_VCF" > "$RAW_STATS"
bcftools stats "$HARDFILTERED_VCF" > "$HARD_STATS"

echo "DONE $(date -Iseconds) 06.2 completed" | tee -a "$LOG_FILE"
echo "Output: $HARDFILTERED_VCF"
echo "Tagged VCF: $TAGGED_VCF"
echo "Stats: $STATS_DIR"

# Step 1: Input SNP VCF from 06.1 normalization/selection
# - joint.norm.biallelic.snps.vcf.gz
#
# Step 2: Apply technical hard filters (QD, QUAL, FS, SOR, MQ, MQRankSum, ReadPosRankSum, DP)
# - VariantFiltration tags failing variants in FILTER column
#
# Step 3: Keep only PASS variants
# - SelectVariants --exclude-filtered true
# - Output: snps.hardfiltered.vcf.gz
#
# Step 4: Generate bcftools stats before/after filtering

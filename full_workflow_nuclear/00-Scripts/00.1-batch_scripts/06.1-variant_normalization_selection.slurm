#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J 06.1-variant_norm
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.out
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
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
module load bioinfo/Bcftools/1.9
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs and validate files
# -------------------------------------------------------------------
REF_GENOME="/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa"
INPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/05-variant_calling/workflow_test/genotypegvcfs"
INPUT_VCF="${INPUT_DIR}/jointgenotypes.raw.vcf.gz"
INPUT_TBI="${INPUT_VCF}.tbi"

OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.1-normalization_selection"
LOG_DIR="${OUTPUT_DIR}/logs"
STATS_DIR="${OUTPUT_DIR}/stats"
LOG_FILE="${LOG_DIR}/06.1-variant_normalization_selection.log"

NORM_VCF="${OUTPUT_DIR}/joint.norm.vcf.gz"
NORM_TBI="${NORM_VCF}.tbi"
SNP_VCF="${OUTPUT_DIR}/joint.norm.biallelic.snps.vcf.gz"
SNP_TBI="${SNP_VCF}.tbi"
INDEL_VCF="${OUTPUT_DIR}/joint.norm.biallelic.indels.vcf.gz"
INDEL_TBI="${INDEL_VCF}.tbi"

RAW_STATS="${STATS_DIR}/jointgenotypes.raw.stats.txt"
NORM_STATS="${STATS_DIR}/joint.norm.stats.txt"
SNP_STATS="${STATS_DIR}/joint.norm.biallelic.snps.stats.txt"
INDEL_STATS="${STATS_DIR}/joint.norm.biallelic.indels.stats.txt"

BCFTOOLS_THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$STATS_DIR"
echo "START $(date -Iseconds) 06.1 variant normalization and selection" > "$LOG_FILE"

if [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found: $REF_GENOME" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$INPUT_VCF" ] || [ ! -s "$INPUT_TBI" ]; then
    echo "ERROR: input VCF or index missing: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -s "${REF_GENOME}.fai" ]; then
    echo "INFO: creating reference fai index" | tee -a "$LOG_FILE"
    samtools faidx "$REF_GENOME" >> "$LOG_FILE" 2>&1
fi
REF_DICT="${REF_GENOME%.fa}.dict"
if [ ! -s "$REF_DICT" ]; then
    echo "INFO: creating reference sequence dictionary" | tee -a "$LOG_FILE"
    gatk CreateSequenceDictionary -R "$REF_GENOME" -O "$REF_DICT" >> "$LOG_FILE" 2>&1
fi

# -------------------------------------------------------------------
# STEP 2: Normalize variants with bcftools norm
# -------------------------------------------------------------------
echo "STEP 2: bcftools norm on raw joint VCF" | tee -a "$LOG_FILE"
bcftools norm \
    --threads "$BCFTOOLS_THREADS" \
    -f "$REF_GENOME" \
    -m-any \
    -Oz \
    -o "$NORM_VCF" \
    "$INPUT_VCF" >> "$LOG_FILE" 2>&1

bcftools index --threads "$BCFTOOLS_THREADS" -t "$NORM_VCF" >> "$LOG_FILE" 2>&1

if [ ! -s "$NORM_VCF" ] || [ ! -s "$NORM_TBI" ]; then
    echo "ERROR: normalized VCF output missing" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 3: Select biallelic SNPs and INDELs with GATK SelectVariants
# -------------------------------------------------------------------
echo "STEP 3: GATK SelectVariants (biallelic SNPs)" | tee -a "$LOG_FILE"
gatk --java-options "-Xms8G -Xmx8G -XX:ParallelGCThreads=2" SelectVariants \
    -R "$REF_GENOME" \
    -V "$NORM_VCF" \
    --select-type-to-include SNP \
    --restrict-alleles-to BIALLELIC \
    -O "$SNP_VCF" >> "$LOG_FILE" 2>&1

echo "STEP 3: GATK SelectVariants (biallelic INDELs)" | tee -a "$LOG_FILE"
gatk --java-options "-Xms8G -Xmx8G -XX:ParallelGCThreads=2" SelectVariants \
    -R "$REF_GENOME" \
    -V "$NORM_VCF" \
    --select-type-to-include INDEL \
    --restrict-alleles-to BIALLELIC \
    -O "$INDEL_VCF" >> "$LOG_FILE" 2>&1

if [ ! -s "$SNP_VCF" ] || [ ! -s "$SNP_TBI" ]; then
    echo "ERROR: SNP VCF output missing" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$INDEL_VCF" ] || [ ! -s "$INDEL_TBI" ]; then
    echo "ERROR: INDEL VCF output missing" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 4: Compute bcftools stats on raw/normalized/selected VCFs
# -------------------------------------------------------------------
echo "STEP 4: bcftools stats" | tee -a "$LOG_FILE"
bcftools stats --threads "$BCFTOOLS_THREADS" "$INPUT_VCF" > "$RAW_STATS"
bcftools stats --threads "$BCFTOOLS_THREADS" "$NORM_VCF" > "$NORM_STATS"
bcftools stats --threads "$BCFTOOLS_THREADS" "$SNP_VCF" > "$SNP_STATS"
bcftools stats --threads "$BCFTOOLS_THREADS" "$INDEL_VCF" > "$INDEL_STATS"

echo "DONE $(date -Iseconds) 06.1 completed" | tee -a "$LOG_FILE"
echo "Output files:"
echo "  - $NORM_VCF"
echo "  - $SNP_VCF"
echo "  - $INDEL_VCF"
echo "Stats files in: $STATS_DIR"

# Step 1: Validate raw joint VCF from 05.3 and reference genome
# - Input: jointgenotypes.raw.vcf.gz + .tbi
# - Ensures .fai/.dict for the reference
#
# Step 2: Normalize variants with bcftools norm
# - Left-align and split multiallelic records (-m-any)
# - Output: joint.norm.vcf.gz + .tbi
#
# Step 3: Select biallelic SNPs and biallelic INDELs
# - Uses GATK SelectVariants with --restrict-alleles-to BIALLELIC
# - Outputs:
#   joint.norm.biallelic.snps.vcf.gz
#   joint.norm.biallelic.indels.vcf.gz
#
# Step 4: Generate bcftools stats reports
# - Creates stats for raw, normalized, SNP, and INDEL VCFs

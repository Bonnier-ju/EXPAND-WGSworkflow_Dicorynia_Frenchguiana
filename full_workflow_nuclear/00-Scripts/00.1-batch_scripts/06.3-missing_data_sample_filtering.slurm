#!/bin/bash
#SBATCH -p workq
#SBATCH --time=48:00:00
#SBATCH -J 06.3-missing_filter
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/PLINK/2.00a4
module load bioinfo/Bcftools/1.9
module load bioinfo/samtools/1.23
module load bioinfo/VCFtools/0.1.16

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs and thresholds
# -------------------------------------------------------------------
INPUT_VCF="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_full/06.2-technical_hard_filter/snps.hardfiltered.vcf.gz"
INPUT_TBI="${INPUT_VCF}.tbi"

OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_full/06.3-missing_data_sample_filtering"
LOG_DIR="${OUTPUT_DIR}/logs"
TMP_DIR="${OUTPUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/06.3-missing_data_sample_filtering.log"

MAX_SITE_MISSING="${MAX_SITE_MISSING:-0.8}"      # keep SNP if >=80% genotypes are present
MAX_SAMPLE_MISSING="${MAX_SAMPLE_MISSING:-0.2}"  # remove sample if >20% missing genotypes

SITES_FILTERED_VCF="${OUTPUT_DIR}/snps.site_missing_filtered.vcf.gz"
FINAL_VCF="${OUTPUT_DIR}/snps.missing_filtered.vcf.gz"
FINAL_TBI="${FINAL_VCF}.tbi"

MISS_PREFIX="${OUTPUT_DIR}/snps.site_missing_filtered"
REMOVE_SAMPLES_TXT="${OUTPUT_DIR}/samples_to_remove_missing.txt"
KEPT_SAMPLES_TXT="${OUTPUT_DIR}/samples_kept_missing.txt"

PLINK_PREFIX="${OUTPUT_DIR}/snps.missing_filtered"
STATS_DIR="${OUTPUT_DIR}/stats"
FINAL_STATS="${STATS_DIR}/snps.missing_filtered.stats.txt"

BCFTOOLS_THREADS="${SLURM_CPUS_PER_TASK:-2}"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$TMP_DIR" "$STATS_DIR"
echo "START $(date -Iseconds) 06.3 missing data and sample filtering" > "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ] || [ ! -s "$INPUT_TBI" ]; then
    echo "ERROR: input VCF missing: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Site-level missing data filtering (vcftools --max-missing)
# -------------------------------------------------------------------
vcftools \
    --gzvcf "$INPUT_VCF" \
    --max-missing "$MAX_SITE_MISSING" \
    --recode --recode-INFO-all \
    --stdout 2>> "$LOG_FILE" | bgzip -c > "$SITES_FILTERED_VCF"

bcftools index --threads "$BCFTOOLS_THREADS" -t "$SITES_FILTERED_VCF" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 3: Sample-level missingness report and sample removal
# -------------------------------------------------------------------
vcftools \
    --gzvcf "$SITES_FILTERED_VCF" \
    --missing-indv \
    --out "$MISS_PREFIX" >> "$LOG_FILE" 2>&1

awk -v thr="$MAX_SAMPLE_MISSING" 'NR>1 && $5>thr {print $1}' "${MISS_PREFIX}.imiss" > "$REMOVE_SAMPLES_TXT"
awk -v thr="$MAX_SAMPLE_MISSING" 'NR>1 && $5<=thr {print $1}' "${MISS_PREFIX}.imiss" > "$KEPT_SAMPLES_TXT"

if [ -s "$REMOVE_SAMPLES_TXT" ]; then
    vcftools \
        --gzvcf "$SITES_FILTERED_VCF" \
        --remove "$REMOVE_SAMPLES_TXT" \
        --recode --recode-INFO-all \
        --stdout 2>> "$LOG_FILE" | bgzip -c > "$FINAL_VCF"
else
    cp "$SITES_FILTERED_VCF" "$FINAL_VCF"
fi

bcftools index --threads "$BCFTOOLS_THREADS" -t "$FINAL_VCF" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 4: Export PLINK2 files (.bed .bim .fam)
# -------------------------------------------------------------------
plink2 \
    --vcf "$FINAL_VCF" \
    --allow-extra-chr \
    --make-bed \
    --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "$FINAL_VCF" ] || [ ! -s "$FINAL_TBI" ]; then
    echo "ERROR: final missing-filtered VCF missing" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "${PLINK_PREFIX}.bed" ] || [ ! -s "${PLINK_PREFIX}.bim" ] || [ ! -s "${PLINK_PREFIX}.fam" ]; then
    echo "ERROR: PLINK outputs missing" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 5: Statistics
# -------------------------------------------------------------------
bcftools stats --threads "$BCFTOOLS_THREADS" "$FINAL_VCF" > "$FINAL_STATS"

echo "DONE $(date -Iseconds) 06.3 completed" | tee -a "$LOG_FILE"
echo "Output VCF: $FINAL_VCF"
echo "PLINK files: ${PLINK_PREFIX}.bed/.bim/.fam"
echo "Removed samples: $(wc -l < "$REMOVE_SAMPLES_TXT")"
echo "Kept samples: $(wc -l < "$KEPT_SAMPLES_TXT")"
echo "Stats: $FINAL_STATS"

# Step 1: Input SNP VCF from 06.2 hard filtering
# - snps.hardfiltered.vcf.gz
#
# Step 2: Site-level missing data filter (vcftools --max-missing 0.8)
# - Remove sites where >20% of genotypes are missing
#
# Step 3: Sample-level missingness report and removal
# - vcftools --missing-indv: compute F_MISS per sample
# - Remove samples with >20% missing genotypes
#
# Step 4: Export PLINK2 binary files
# - snps.missing_filtered.bed/.bim/.fam
#
# Step 5: bcftools stats on final VCF (Ts/Tv, SNP count)

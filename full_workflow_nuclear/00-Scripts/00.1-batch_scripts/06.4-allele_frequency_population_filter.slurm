#!/bin/bash
#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 06.4-allele_freq_filter
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u
set -o pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/plink2
module load bioinfo/bcftools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
IN_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.3-missing_data_sample_filtering"
IN_VCF="${IN_DIR}/snps.missing_filtered.vcf.gz"
IN_TBI="${IN_VCF}.tbi"
IN_BFILE="${IN_DIR}/snps.missing_filtered"

OUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.4-allele_frequency_population_filter"
LOG_DIR="${OUT_DIR}/logs"
STATS_DIR="${OUT_DIR}/stats"
TMP_DIR="${OUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/06.4-allele_frequency_population_filter.log"

OUT_BFILE="${OUT_DIR}/snps_non_informative_filtered"
OUT_VCF_PREFIX="${OUT_DIR}/snps_non_informative_filtered"
OUT_VCF="${OUT_VCF_PREFIX}.vcf.gz"
OUT_TBI="${OUT_VCF}.tbi"
OUT_STATS="${STATS_DIR}/snps_non_informative_filtered.stats.txt"

# Remove singleton/doubleton by MAC threshold (MAC <=2 removed)
MIN_MAC="${MIN_MAC:-3}"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$STATS_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 06.4 allele frequency/population filtering" > "$LOG_FILE"

if [ ! -s "$IN_VCF" ] || [ ! -s "$IN_TBI" ]; then
    echo "ERROR: input VCF missing: $IN_VCF" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "${IN_BFILE}.bed" ] || [ ! -s "${IN_BFILE}.bim" ] || [ ! -s "${IN_BFILE}.fam" ]; then
    echo "ERROR: input PLINK files missing: ${IN_BFILE}.bed/.bim/.fam" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Filter non-informative variants and rare MAC <=2
# -------------------------------------------------------------------
# Using --mac MIN_MAC removes monomorphic variants and singleton/doubleton sites
# when MIN_MAC=3.
plink2 \
    --bfile "$IN_BFILE" \
    --allow-extra-chr \
    --mac "$MIN_MAC" \
    --make-bed \
    --out "$OUT_BFILE" >> "$LOG_FILE" 2>&1

if [ ! -s "${OUT_BFILE}.bed" ] || [ ! -s "${OUT_BFILE}.bim" ] || [ ! -s "${OUT_BFILE}.fam" ]; then
    echo "ERROR: filtered PLINK files missing" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 3: Export filtered VCF
# -------------------------------------------------------------------
plink2 \
    --bfile "$OUT_BFILE" \
    --allow-extra-chr \
    --recode vcf bgz \
    --out "$OUT_VCF_PREFIX" >> "$LOG_FILE" 2>&1

bcftools index -t "$OUT_VCF" >> "$LOG_FILE" 2>&1

if [ ! -s "$OUT_VCF" ] || [ ! -s "$OUT_TBI" ]; then
    echo "ERROR: filtered VCF output missing: $OUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 4: Variant statistics after filtering
# -------------------------------------------------------------------
bcftools stats "$OUT_VCF" > "$OUT_STATS"

echo "DONE $(date -Iseconds) 06.4 completed" | tee -a "$LOG_FILE"
echo "Output VCF: $OUT_VCF"
echo "Output PLINK: ${OUT_BFILE}.bed/.bim/.fam"
echo "Stats: $OUT_STATS"

# Step 1: Input from 06.3
# - snps.missing_filtered.vcf.gz
# - snps.missing_filtered.bed/.bim/.fam
#
# Step 2: Population-level allele filtering
# - Removes non-informative and ultra-rare variants with --mac 3
# - Removes singleton/doubleton sites (MAC 1 and 2)
#
# Step 3: Export filtered VCF
# - snps_non_informative_filtered.vcf.gz + .tbi
#
# Step 4: Compute bcftools stats

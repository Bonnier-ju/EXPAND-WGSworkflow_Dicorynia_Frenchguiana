#!/bin/bash
#SBATCH -p workq
#SBATCH --time=02:00:00
#SBATCH -J 08.2-genotype_matrix
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/PLINK/2.00a4
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# -------------------------------------------------------------------
# STEP 1: Define paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"

# Input: LD-pruned GEA VCF (155 individuals, MAF>=5%) from 06.9
GEA_VCF="${BASE}/06-variant_filtering/workflow_full/06.9-gea_vcf_preparation/snps_gea_ldpruned.vcf.gz"

# Metadata for individual-to-site mapping
SAMPLE_SHEET="${BASE}/metadata/sample_sheet_complete_filtered.csv"

OUT_DIR="${BASE}/08-gea/workflow_full/08.2-genotype_matrix"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/08.2-genotype_matrix.log"

PLINK_PREFIX="${TMP_DIR}/snps_gea_ldpruned"
RAW_012="${TMP_DIR}/snps_gea_ldpruned.raw"   # plink2 --export A output
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/08.2-genotype_matrix.R"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 08.2 genotype matrix preparation" > "$LOG_FILE"

if [ ! -s "$GEA_VCF" ]; then
    echo "ERROR: GEA VCF not found: $GEA_VCF" | tee -a "$LOG_FILE"
    echo "       Run 06.9-gea_vcf_preparation.slurm first." | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Convert VCF to 0/1/2 dosage matrix via plink2 --export A
# Dosage is counted for the ALT allele (0 = hom REF, 1 = het, 2 = hom ALT).
# -------------------------------------------------------------------
echo "STEP 2: converting VCF to 0/1/2 dosage matrix (plink2 --export A)" | tee -a "$LOG_FILE"

plink2 \
    --vcf "$GEA_VCF" \
    --allow-extra-chr \
    --threads "$THREADS" \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 50 \
    --export A \
    --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "$RAW_012" ]; then
    echo "ERROR: plink2 --export A output missing: $RAW_012" | tee -a "$LOG_FILE"
    exit 1
fi

N_ROWS=$(awk 'END{print NR-1}' "$RAW_012")
N_COLS=$(head -1 "$RAW_012" | awk '{print NF-6}')
echo "INFO: raw matrix: ${N_ROWS} individuals x ${N_COLS} SNPs" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: R — clean matrix, add sample metadata, write final outputs
# -------------------------------------------------------------------
echo "STEP 3: R post-processing (clean matrix, sample metadata)" | tee -a "$LOG_FILE"

Rscript --vanilla "$R_SCRIPT" \
    "$RAW_012" \
    "$SAMPLE_SHEET" \
    "$OUT_DIR" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 08.2 genotype matrix completed" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}"

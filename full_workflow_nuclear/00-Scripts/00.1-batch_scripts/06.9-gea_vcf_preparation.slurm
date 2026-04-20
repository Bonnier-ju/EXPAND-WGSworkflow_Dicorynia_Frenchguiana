#!/bin/bash
#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 06.9-gea_vcf_preparation
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
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

# -------------------------------------------------------------------
# STEP 1: Define paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"

# Input: TE-filtered VCF (164 individuals, 17.4 M SNPs)
IN_VCF="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"

# Metadata: sample sheet with site assignments
SAMPLE_SHEET="${BASE}/metadata/sample_sheet_complete_filtered.csv"

OUT_DIR="${BASE}/06-variant_filtering/workflow_full/06.9-gea_vcf_preparation"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/06.9-gea_vcf_preparation.log"

# Output files
KEEP_LIST="${OUT_DIR}/gea_155_individuals.txt"
MAF_PREFIX="${TMP_DIR}/snps_gea_maf05"
LD_PREFIX="${OUT_DIR}/snps_gea_ldpruned"
PRUNE_IN="${LD_PREFIX}.prune.in"
GEA_VCF="${LD_PREFIX}.vcf.gz"
GEA_TBI="${GEA_VCF}.tbi"

# LD pruning parameters (same as 06.6b, tuned on 06.6a LD decay results)
# Physical LD decays within 1-5 kb; 50 kb window gives ample coverage.
LD_WINDOW_KB=50
LD_STEP_SNP=1
LD_R2=0.4
MAF_THRESHOLD=0.05

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 06.9 GEA VCF preparation" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Build the 155 FG individual keep-list
# Exclude: Cameroun_Benin (outgroup), Venezuela/Herbier (DBR_AGO),
#          DBR_AGJ/K/L/M (treemutation/Angela, St_georges)
# -------------------------------------------------------------------
echo "STEP 2: building 155-individual keep-list" | tee -a "$LOG_FILE"

awk -F',' 'NR>1 \
    && $4 != "Cameroun_Benin" \
    && $4 != "Herbier" \
    && $4 != "Venezuela" \
    && $1 !~ /^DBR_AG[JKLM]$/ \
    { print $1 }' "$SAMPLE_SHEET" | sort > "$KEEP_LIST"

N_KEEP=$(wc -l < "$KEEP_LIST")
echo "INFO: ${N_KEEP} individuals in keep-list: $KEEP_LIST" | tee -a "$LOG_FILE"

if [ "$N_KEEP" -ne 155 ]; then
    echo "ERROR: expected 155 individuals, got ${N_KEEP}" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 3: Subset to 155 individuals + MAF >= 5% filter
# Variant IDs set to CHR:POS:REF:ALT for unambiguous downstream use.
# -------------------------------------------------------------------
echo "STEP 3: subsetting to 155 individuals, MAF >= ${MAF_THRESHOLD}" | tee -a "$LOG_FILE"

plink2 \
    --vcf "$IN_VCF" \
    --keep "$KEEP_LIST" \
    --maf "$MAF_THRESHOLD" \
    --allow-extra-chr \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 50 \
    --threads "$THREADS" \
    --make-bed \
    --out "$MAF_PREFIX" >> "$LOG_FILE" 2>&1

N_MAF=$(wc -l < "${MAF_PREFIX}.bim")
echo "INFO: ${N_MAF} SNPs retained after MAF filter" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: LD pruning
# Same parameters as 06.6b (50 kb window, step=1, r2=0.4).
# -------------------------------------------------------------------
echo "STEP 4: LD pruning (${LD_WINDOW_KB}kb window, r2 <= ${LD_R2})" | tee -a "$LOG_FILE"

plink2 \
    --bfile "$MAF_PREFIX" \
    --allow-extra-chr \
    --threads "$THREADS" \
    --indep-pairwise "${LD_WINDOW_KB}kb" "$LD_STEP_SNP" "$LD_R2" \
    --out "$LD_PREFIX" >> "$LOG_FILE" 2>&1

N_PRUNE_IN=$(wc -l < "$PRUNE_IN")
echo "INFO: ${N_PRUNE_IN} SNPs retained after LD pruning" | tee -a "$LOG_FILE"

plink2 \
    --bfile "$MAF_PREFIX" \
    --allow-extra-chr \
    --threads "$THREADS" \
    --extract "$PRUNE_IN" \
    --make-bed \
    --out "$LD_PREFIX" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 5: Export pruned VCF
# -------------------------------------------------------------------
echo "STEP 5: exporting pruned VCF" | tee -a "$LOG_FILE"

plink2 \
    --bfile "$LD_PREFIX" \
    --allow-extra-chr \
    --threads "$THREADS" \
    --export vcf bgz \
    --out "$LD_PREFIX" >> "$LOG_FILE" 2>&1

bcftools index --threads "$THREADS" -t "$GEA_VCF" >> "$LOG_FILE" 2>&1

if [ ! -s "$GEA_VCF" ] || [ ! -s "$GEA_TBI" ]; then
    echo "ERROR: pruned GEA VCF missing: $GEA_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 6: Summary
# -------------------------------------------------------------------
N_SAMPLES=$(bcftools query -l "$GEA_VCF" | wc -l)
N_SNPS=$(bcftools stats "$GEA_VCF" 2>/dev/null | grep "^SN.*number of SNPs:" | cut -f4)

echo "DONE $(date -Iseconds) 06.9 GEA VCF preparation completed" | tee -a "$LOG_FILE"
echo "Individuals : ${N_SAMPLES}"                   | tee -a "$LOG_FILE"
echo "SNPs        : ${N_SNPS}"                      | tee -a "$LOG_FILE"
echo "Keep-list   : ${KEEP_LIST}"                   | tee -a "$LOG_FILE"
echo "GEA VCF     : ${GEA_VCF}"                     | tee -a "$LOG_FILE"

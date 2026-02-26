#!/bin/bash
#SBATCH -p workq
#SBATCH --time=20:00:00
#SBATCH -J 06.6-ld_pruning
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
module load bioinfo/bcftools/1.21
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
IN_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.5-transposable_element_filtering"
IN_VCF="${IN_DIR}/snps_TE_filtered.vcf.gz"
IN_TBI="${IN_VCF}.tbi"

OUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.6-ld_pruning"
LOG_DIR="${OUT_DIR}/logs"
STATS_DIR="${OUT_DIR}/stats"
TMP_DIR="${OUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/06.6-ld_pruning.log"

WORK_PREFIX="${OUT_DIR}/snps_TE_filtered_for_ld"
PRUNE_PREFIX="${OUT_DIR}/snps_TE_filtered.ld_pruned"
PRUNE_IN="${PRUNE_PREFIX}.prune.in"
PRUNE_OUT="${PRUNE_PREFIX}.prune.out"
PRUNED_VCF="${PRUNE_PREFIX}.vcf.gz"
PRUNED_TBI="${PRUNED_VCF}.tbi"

# Optional sample exclusion file (one sample ID per line).
EXCLUDE_SAMPLES_FILE="${EXCLUDE_SAMPLES_FILE:-}"

# LD pruning parameters (same idea as your historical script).
LD_WINDOW_KB="${LD_WINDOW_KB:-10000}"
LD_STEP_SNP="${LD_STEP_SNP:-5000}"
LD_R2="${LD_R2:-0.4}"

# Optional full LD computation output (0: disabled, 1: enabled).
COMPUTE_LD_R2="${COMPUTE_LD_R2:-0}"

RAW_STATS="${STATS_DIR}/snps_TE_filtered.stats.txt"
PRUNED_STATS="${STATS_DIR}/snps_TE_filtered.ld_pruned.stats.txt"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$STATS_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 06.6 LD pruning" > "$LOG_FILE"

if [ ! -s "$IN_VCF" ] || [ ! -s "$IN_TBI" ]; then
    echo "ERROR: input VCF missing: $IN_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

REMOVE_ARGS=()
if [ -n "$EXCLUDE_SAMPLES_FILE" ] && [ -s "$EXCLUDE_SAMPLES_FILE" ]; then
    REMOVE_ARGS=(--remove "$EXCLUDE_SAMPLES_FILE")
    echo "INFO: excluding samples listed in $EXCLUDE_SAMPLES_FILE" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 2: Convert VCF to PLINK binary
# -------------------------------------------------------------------
plink2 \
    --vcf "$IN_VCF" \
    --allow-extra-chr \
    "${REMOVE_ARGS[@]}" \
    --make-bed \
    --out "$WORK_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "${WORK_PREFIX}.bed" ] || [ ! -s "${WORK_PREFIX}.bim" ] || [ ! -s "${WORK_PREFIX}.fam" ]; then
    echo "ERROR: intermediate PLINK files missing" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 3: LD pruning
# -------------------------------------------------------------------
plink2 \
    --bfile "$WORK_PREFIX" \
    --allow-extra-chr \
    --indep-pairwise "$LD_WINDOW_KB" "$LD_STEP_SNP" "$LD_R2" \
    --out "$PRUNE_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "$PRUNE_IN" ]; then
    echo "ERROR: prune.in file is missing/empty: $PRUNE_IN" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 4: Keep only pruned SNP set and export
# -------------------------------------------------------------------
plink2 \
    --bfile "$WORK_PREFIX" \
    --allow-extra-chr \
    --extract "$PRUNE_IN" \
    --make-bed \
    --out "$PRUNE_PREFIX" >> "$LOG_FILE" 2>&1

plink2 \
    --bfile "$PRUNE_PREFIX" \
    --allow-extra-chr \
    --export vcf bgz \
    --out "$PRUNE_PREFIX" >> "$LOG_FILE" 2>&1

bcftools index -t "$PRUNED_VCF" >> "$LOG_FILE" 2>&1

if [ ! -s "$PRUNED_VCF" ] || [ ! -s "$PRUNED_TBI" ]; then
    echo "ERROR: pruned VCF output missing: $PRUNED_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 5: Optional LD pairwise output (for diagnostics/plot)
# -------------------------------------------------------------------
if [ "$COMPUTE_LD_R2" -eq 1 ]; then
    LD_OUT="${OUT_DIR}/ld_results"
    plink2 \
        --bfile "$WORK_PREFIX" \
        --allow-extra-chr \
        --r2 \
        --ld-window-kb 100000 \
        --out "$LD_OUT" >> "$LOG_FILE" 2>&1 || true

    if [ -s "${LD_OUT}.vcor" ]; then
        awk 'NR==1 || rand() < 0.1' "${LD_OUT}.vcor" > "${LD_OUT}_10percent.vcor"
    fi
fi

# -------------------------------------------------------------------
# STEP 6: Stats before/after pruning
# -------------------------------------------------------------------
bcftools stats "$IN_VCF" > "$RAW_STATS"
bcftools stats "$PRUNED_VCF" > "$PRUNED_STATS"

RAW_N=$(bcftools view -H "$IN_VCF" | wc -l)
PRUNED_N=$(bcftools view -H "$PRUNED_VCF" | wc -l)
REMOVED_N=$((RAW_N - PRUNED_N))

echo "DONE $(date -Iseconds) 06.6 completed" | tee -a "$LOG_FILE"
echo "Input VCF         : $IN_VCF"
echo "Pruned VCF        : $PRUNED_VCF"
echo "Prune.in          : $PRUNE_IN"
echo "Prune.out         : $PRUNE_OUT"
echo "Variants raw      : $RAW_N"
echo "Variants pruned   : $PRUNED_N"
echo "Variants removed  : $REMOVED_N"
echo "Stats             : $STATS_DIR"

# Step 1: Input from 06.5 TE-filtered VCF
# Step 2: Convert to PLINK binary format
# Step 3: Run LD pruning (--indep-pairwise)
# Step 4: Keep pruned SNP set and export VCF
# Step 5: Optional pairwise LD table for diagnostics
# Step 6: Generate before/after stats

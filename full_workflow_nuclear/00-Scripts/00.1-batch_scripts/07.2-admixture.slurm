#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=120:00:00
#SBATCH -J 07.2-admixture
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/PLINK/1.90b7
module load bioinfo/ADMIXTURE/1.3.0

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"

VCF_66="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
VCF_65="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"

if [ -s "$VCF_66" ]; then
    INPUT_VCF="$VCF_66"
    INPUT_STEP="06.6b"
else
    INPUT_VCF="$VCF_65"
    INPUT_STEP="06.5"
fi

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.2-admixture"
LOG_DIR="${OUT_DIR}/logs"
LOG_FILE="${LOG_DIR}/07.2-admixture.log"

# K range for ADMIXTURE runs
K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-10}"

# Cross-validation fold count (10-fold is standard)
CV_FOLDS="${CV_FOLDS:-10}"

PLINK_THREADS="${SLURM_CPUS_PER_TASK:-16}"
ADMIXTURE_THREADS="${SLURM_CPUS_PER_TASK:-16}"

PLINK_PREFIX="${OUT_DIR}/admixture_input"

mkdir -p "$OUT_DIR" "$LOG_DIR"
echo "START $(date -Iseconds) 07.2 ADMIXTURE" > "$LOG_FILE"
echo "Input VCF: $INPUT_VCF (from $INPUT_STEP)" | tee -a "$LOG_FILE"
echo "K range: ${K_MIN}..${K_MAX} | CV folds: ${CV_FOLDS} | threads: ${ADMIXTURE_THREADS}" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Convert VCF to PLINK binary format
# -------------------------------------------------------------------
echo "STEP 2: VCF → PLINK binary" | tee -a "$LOG_FILE"
plink \
    --vcf "$INPUT_VCF" \
    --make-bed \
    --allow-extra-chr \
    --const-fid \
    --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "${PLINK_PREFIX}.bed" ]; then
    echo "ERROR: PLINK conversion failed" | tee -a "$LOG_FILE"
    exit 1
fi

N_INDS=$(wc -l < "${PLINK_PREFIX}.fam")
N_SNPS=$(wc -l < "${PLINK_PREFIX}.bim")
echo "INFO: ${N_INDS} individuals, ${N_SNPS} SNPs" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Rename chromosomes in .bim to numeric (required by ADMIXTURE)
# -------------------------------------------------------------------
# ADMIXTURE requires numeric chromosome IDs (1-22/X/Y).
# Dicorynia scaffold names are non-numeric — set all to 0 (unknown/unplaced).
echo "STEP 3: renaming chromosomes in .bim → 0" | tee -a "$LOG_FILE"
awk '{$1=0; print $0}' "${PLINK_PREFIX}.bim" > "${PLINK_PREFIX}_temp.bim"
mv "${PLINK_PREFIX}_temp.bim" "${PLINK_PREFIX}.bim"

# -------------------------------------------------------------------
# STEP 4: Run ADMIXTURE for K = K_MIN..K_MAX
# -------------------------------------------------------------------
# ADMIXTURE writes output files (prefix.K.Q, prefix.K.P) in the
# current working directory — must cd to OUT_DIR before running.
echo "STEP 4: running ADMIXTURE K=${K_MIN}..${K_MAX}" | tee -a "$LOG_FILE"
cd "$OUT_DIR"

for K in $(seq "$K_MIN" "$K_MAX"); do
    echo "  Running K=${K}..." | tee -a "$LOG_FILE"
    admixture \
        --cv="$CV_FOLDS" \
        -j"$ADMIXTURE_THREADS" \
        "${PLINK_PREFIX}.bed" "$K" \
        | tee "${OUT_DIR}/log_K${K}.out" >> "$LOG_FILE"
    echo "  K=${K} completed." | tee -a "$LOG_FILE"
done

# -------------------------------------------------------------------
# STEP 5: Collect cross-validation errors
# -------------------------------------------------------------------
echo "STEP 5: collecting cross-validation errors" | tee -a "$LOG_FILE"
CV_FILE="${OUT_DIR}/cross_validation.txt"
grep "CV error" "${OUT_DIR}"/log_K*.out > "$CV_FILE" || true
echo "CV summary: $CV_FILE" | tee -a "$LOG_FILE"
cat "$CV_FILE" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 6: Associate individual names with Q matrices
# -------------------------------------------------------------------
# Extract individual IDs (IID column 2 of .fam) and prepend to each
# .Q file to produce named Q files readable by the R plotting script.
echo "STEP 6: associating individual names with Q matrices" | tee -a "$LOG_FILE"
IDS_FILE="${OUT_DIR}/individual_ids.txt"
awk '{print $2}' "${PLINK_PREFIX}.fam" > "$IDS_FILE"

BED_BASE=$(basename "$PLINK_PREFIX")
for K in $(seq "$K_MIN" "$K_MAX"); do
    Q_FILE="${OUT_DIR}/${BED_BASE}.${K}.Q"
    NAMED_Q="${OUT_DIR}/${BED_BASE}.${K}.named.Q"
    if [ -s "$Q_FILE" ]; then
        paste "$IDS_FILE" "$Q_FILE" > "$NAMED_Q"
        echo "  Created: $(basename "$NAMED_Q")" | tee -a "$LOG_FILE"
    else
        echo "  WARN: Q file missing for K=${K}: $Q_FILE" | tee -a "$LOG_FILE"
    fi
done

echo "DONE $(date -Iseconds) 07.2 ADMIXTURE completed" | tee -a "$LOG_FILE"
echo "Outputs: $OUT_DIR"
echo "CV errors: $CV_FILE"

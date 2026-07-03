#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=24:00:00
#SBATCH -J 07.12b-smcpp_split
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   sbatch 07.12b-smcpp_split.slurm
#   Depends on: 07.12a-smcpp_estimate.slurm (marginal models required)
#
# Joint demographic inference for all 6 pairwise combinations of the
# 4 ADMIXTURE clusters using smc++ split.  Estimates:
#   - Per-population Ne(t) conditioned on the split event
#   - Divergence time (T_split) between each pair
#
# Pipeline per pair (Pop_A × Pop_B):
#   STEP 3: vcf2smc cross-population
#           NDIST_SPLIT=10 distinguished individuals from each population
#           → 10 (from A) + 10 (from B) = 20 files per scaffold × 13 = 260
#           → 6 pairs × 260 = 1,560 files total (parallel)
#   STEP 4: smc++ split
#           Initialises from marginal models (07.12a), estimates T_split
#           Pairs run sequentially, --cores CPUS each
#   STEP 5: smc++ plot --csv export per pair
#   STEP 6: R — pairwise split-time heatmap + joint Ne(t) panels
#
# 6 pairs: Pop_1/2, Pop_1/3, Pop_1/4, Pop_2/3, Pop_2/4, Pop_3/4
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment
# -------------------------------------------------------------------
module purge
module load devel/Miniconda/Miniconda3
module load bioinfo/SMC++/1.15.5
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

CPUS="${SLURM_CPUS_PER_TASK:-16}"

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL="${K_GLOBAL:-4}"
MU_REF="5e-9"
YEAR_PER_GEN=30
NDIST_SPLIT=10    # distinguished individuals per population per pair

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
TE_BED="${BASE}/07-population_structure/workflow_full/07.11-stairway_plot/07.11b-compute_L/te_merged.bed"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.12b-smcpp_split.R"
OUT_BASE="${BASE}/07-population_structure/workflow_full/07.12-smcpp"
LOG_DIR="${OUT_BASE}/logs"
FG_VCF="${LOG_DIR}/fg_155_chr13.vcf.gz"

SCAFFOLDS=(
    Super-Scaffold_1  Super-Scaffold_2  Super-Scaffold_4
    Super-Scaffold_5  Super-Scaffold_6  Super-Scaffold_7
    Super-Scaffold_8  Super-Scaffold_12 Super-Scaffold_13
    Super-Scaffold_14 Super-Scaffold_16 Super-Scaffold_17
    Super-Scaffold_21
)

# All 6 pairwise combinations of 4 populations
declare -a PAIRS=("Pop_1:Pop_2" "Pop_1:Pop_3" "Pop_1:Pop_4"
                  "Pop_2:Pop_3" "Pop_2:Pop_4" "Pop_3:Pop_4")

mkdir -p "${OUT_BASE}/split" "${OUT_BASE}/csv/split" "${OUT_BASE}/plots"
LOG_FILE="${LOG_DIR}/07.12b-smcpp_split.log"
echo "START $(date -Iseconds) 07.12b SMC++ split K=${K_GLOBAL}" > "$LOG_FILE"
echo "Pairs : ${PAIRS[*]}" | tee -a "$LOG_FILE"
echo "NDIST_SPLIT = ${NDIST_SPLIT} | CPUs = ${CPUS}" | tee -a "$LOG_FILE"

for f in "$TE_BED" "$FG_VCF"; do
    [ -s "$f" ] || { echo "ERROR: not found: $f" | tee -a "$LOG_FILE"; exit 1; }
done

# Check that all marginal models from 07.12a exist
for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    MODEL="${OUT_BASE}/${POP}/estimate/model.final.json"
    [ -s "$MODEL" ] || {
        echo "ERROR: marginal model missing for ${POP}: ${MODEL}" | tee -a "$LOG_FILE"
        echo "       Run 07.12a-smcpp_estimate.slurm first." | tee -a "$LOG_FILE"
        exit 1
    }
done
echo "INFO: all marginal models verified" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Load sample lists (written by 07.12a)
# -------------------------------------------------------------------
echo "STEP 2: loading sample lists" | tee -a "$LOG_FILE"
for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    SFILE="${LOG_DIR}/samples_${POP}.txt"
    [ -s "$SFILE" ] || {
        echo "ERROR: sample list missing: $SFILE — run 07.12a first" | tee -a "$LOG_FILE"
        exit 1
    }
    N=$(wc -l < "$SFILE")
    echo "  ${POP}: ${N} individuals" | tee -a "$LOG_FILE"
done

# -------------------------------------------------------------------
# STEP 3: vcf2smc cross-population
#
# For each pair (Pop_A, Pop_B) and each scaffold, generate:
#   - NDIST_SPLIT files with distinguished individuals from Pop_A
#   - NDIST_SPLIT files with distinguished individuals from Pop_B
# Each file contains genotypes for ALL individuals of BOTH populations
# so SMC++ can model cross-population coalescence patterns.
#
# Command:
#   smc++ vcf2smc --mask TE_BED -d DIST DIST VCF OUT SCAFFOLD \
#                 POP_A:csv_A POP_B:csv_B
# -------------------------------------------------------------------
echo "STEP 3: vcf2smc cross-population (parallel)" | tee -a "$LOG_FILE"

SPLIT_CMDS="${LOG_DIR}/vcf2smc_split_commands.txt"
> "$SPLIT_CMDS"
N_SKIP=0
N_GEN=0

for PAIR in "${PAIRS[@]}"; do
    POP_A="${PAIR%%:*}"
    POP_B="${PAIR##*:}"
    PAIR_DIR="${OUT_BASE}/split/${POP_A}_${POP_B}"
    mkdir -p "${PAIR_DIR}/vcf2smc"

    SAMPLES_A_FILE="${LOG_DIR}/samples_${POP_A}.txt"
    SAMPLES_B_FILE="${LOG_DIR}/samples_${POP_B}.txt"
    N_A=$(wc -l < "$SAMPLES_A_FILE")
    N_B=$(wc -l < "$SAMPLES_B_FILE")
    N_DIST_A=$(( N_A < NDIST_SPLIT ? N_A : NDIST_SPLIT ))
    N_DIST_B=$(( N_B < NDIST_SPLIT ? N_B : NDIST_SPLIT ))
    SAMPLES_A_CSV=$(paste -sd',' "$SAMPLES_A_FILE")
    SAMPLES_B_CSV=$(paste -sd',' "$SAMPLES_B_FILE")

    echo "  ${POP_A}×${POP_B}: ${N_DIST_A} + ${N_DIST_B} dist, 13 scaffolds" | tee -a "$LOG_FILE"

    # Distinguished individuals from Pop_A
    while IFS= read -r DIST; do
        for CHR in "${SCAFFOLDS[@]}"; do
            OUT_SMC="${PAIR_DIR}/vcf2smc/${POP_A}.${CHR}.${DIST}.smc.gz"
            if [ -s "$OUT_SMC" ]; then N_SKIP=$((N_SKIP + 1)); continue; fi
            printf 'smc++ vcf2smc --mask "%s" -d "%s" "%s" "%s" "%s" "%s" "%s:%s" "%s:%s"\n' \
                "$TE_BED" "$DIST" "$DIST" "$FG_VCF" "$OUT_SMC" "$CHR" \
                "$POP_A" "$SAMPLES_A_CSV" "$POP_B" "$SAMPLES_B_CSV" \
                >> "$SPLIT_CMDS"
            N_GEN=$((N_GEN + 1))
        done
    done < <(head -n "${N_DIST_A}" "$SAMPLES_A_FILE")

    # Distinguished individuals from Pop_B
    while IFS= read -r DIST; do
        for CHR in "${SCAFFOLDS[@]}"; do
            OUT_SMC="${PAIR_DIR}/vcf2smc/${POP_B}.${CHR}.${DIST}.smc.gz"
            if [ -s "$OUT_SMC" ]; then N_SKIP=$((N_SKIP + 1)); continue; fi
            printf 'smc++ vcf2smc --mask "%s" -d "%s" "%s" "%s" "%s" "%s" "%s:%s" "%s:%s"\n' \
                "$TE_BED" "$DIST" "$DIST" "$FG_VCF" "$OUT_SMC" "$CHR" \
                "$POP_A" "$SAMPLES_A_CSV" "$POP_B" "$SAMPLES_B_CSV" \
                >> "$SPLIT_CMDS"
            N_GEN=$((N_GEN + 1))
        done
    done < <(head -n "${N_DIST_B}" "$SAMPLES_B_FILE")
done

echo "  INFO: ${N_GEN} split vcf2smc commands, ${N_SKIP} already done" | tee -a "$LOG_FILE"

if [ "$N_GEN" -gt 0 ]; then
    echo "  running in parallel (-j ${CPUS})..." | tee -a "$LOG_FILE"
    parallel --jobs "${CPUS}" < "$SPLIT_CMDS" >> "$LOG_FILE" 2>&1 \
        || echo "WARN: some vcf2smc split jobs failed" | tee -a "$LOG_FILE"
    echo "  vcf2smc split done $(date -Iseconds)" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 4: smc++ split — joint Ne(t) + divergence time per pair
#
# Uses marginal models (07.12a model.final.json) as initialisation.
# smc++ split optimises a joint model that includes a split parameter:
# Ne of each population follows its own trajectory; they merge into
# an ancestral population at T_split generations in the past.
#
# Pairs run SEQUENTIALLY (memory and CPU constraints).
# -------------------------------------------------------------------
echo "STEP 4: smc++ split (pairs sequential)" | tee -a "$LOG_FILE"

for PAIR in "${PAIRS[@]}"; do
    POP_A="${PAIR%%:*}"
    POP_B="${PAIR##*:}"
    PAIR_DIR="${OUT_BASE}/split/${POP_A}_${POP_B}"
    MODEL_A="${OUT_BASE}/${POP_A}/estimate/model.final.json"
    MODEL_B="${OUT_BASE}/${POP_B}/estimate/model.final.json"
    SPLIT_MODEL="${PAIR_DIR}/estimate/model.final.json"

    if [ -s "$SPLIT_MODEL" ]; then
        echo "  [${POP_A}×${POP_B}] split model already exists — skipping" | tee -a "$LOG_FILE"
        continue
    fi

    N_SMC=$(ls "${PAIR_DIR}/vcf2smc/"*.smc.gz 2>/dev/null | wc -l)
    if [ "$N_SMC" -eq 0 ]; then
        echo "  WARN: [${POP_A}×${POP_B}] no vcf2smc files — skipping split" | tee -a "$LOG_FILE"
        continue
    fi

    mkdir -p "${PAIR_DIR}/estimate"
    echo "  [${POP_A}×${POP_B}] smc++ split: ${N_SMC} files $(date -Iseconds)" | tee -a "$LOG_FILE"

    smc++ split \
        --cores "${CPUS}" \
        -o "${PAIR_DIR}/estimate" \
        "${MODEL_A}" \
        "${MODEL_B}" \
        "${PAIR_DIR}/vcf2smc/"*.smc.gz >> "$LOG_FILE" 2>&1 \
        || echo "  WARN: [${POP_A}×${POP_B}] split returned non-zero" | tee -a "$LOG_FILE"

    if [ -s "$SPLIT_MODEL" ]; then
        echo "  [${POP_A}×${POP_B}] split done $(date -Iseconds)" | tee -a "$LOG_FILE"
    else
        echo "  WARN: [${POP_A}×${POP_B}] split model not found" | tee -a "$LOG_FILE"
    fi
done

# -------------------------------------------------------------------
# STEP 5: smc++ plot --csv — export split Ne(t) trajectories for R
# -------------------------------------------------------------------
echo "STEP 5: smc++ plot --csv export (split models)" | tee -a "$LOG_FILE"

for PAIR in "${PAIRS[@]}"; do
    POP_A="${PAIR%%:*}"
    POP_B="${PAIR##*:}"
    PAIR_DIR="${OUT_BASE}/split/${POP_A}_${POP_B}"
    SPLIT_MODEL="${PAIR_DIR}/estimate/model.final.json"
    CSV_PREFIX="${OUT_BASE}/csv/split/Ne_${POP_A}_${POP_B}"
    CSV_FILE="${CSV_PREFIX}.csv"

    if [ ! -s "$SPLIT_MODEL" ]; then
        echo "  WARN: [${POP_A}×${POP_B}] model not found — skipping" | tee -a "$LOG_FILE"
        continue
    fi
    if [ -s "$CSV_FILE" ]; then
        echo "  [${POP_A}×${POP_B}] CSV exists — skipping" | tee -a "$LOG_FILE"
        continue
    fi

    smc++ plot \
        -g "${YEAR_PER_GEN}" \
        --csv \
        "${CSV_PREFIX}" \
        "${SPLIT_MODEL}" >> "$LOG_FILE" 2>&1 \
        || echo "  WARN: [${POP_A}×${POP_B}] smc++ plot failed" | tee -a "$LOG_FILE"

    [ -s "$CSV_FILE" ] && echo "  [${POP_A}×${POP_B}] CSV: ${CSV_FILE}" | tee -a "$LOG_FILE"
done

# -------------------------------------------------------------------
# STEP 6: R — split-time visualisation
# -------------------------------------------------------------------
echo "STEP 6: R visualization" | tee -a "$LOG_FILE"

if [ -s "$R_SCRIPT" ]; then
    Rscript --vanilla "$R_SCRIPT" \
        "$OUT_BASE" \
        "$K_GLOBAL" \
        "$MU_REF" \
        "$YEAR_PER_GEN" >> "$LOG_FILE" 2>&1 \
        || echo "WARN: R script returned non-zero" | tee -a "$LOG_FILE"
else
    echo "WARN: R script not found: $R_SCRIPT" | tee -a "$LOG_FILE"
fi

echo "DONE $(date -Iseconds) 07.12b SMC++ split" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_BASE}/split/"

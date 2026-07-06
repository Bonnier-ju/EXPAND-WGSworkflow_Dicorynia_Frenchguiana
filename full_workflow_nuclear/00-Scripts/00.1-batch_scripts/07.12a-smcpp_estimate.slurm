#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=160:00:00
#SBATCH -J 07.12a-smcpp_estimate
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   sbatch 07.12a-smcpp_estimate.slurm
#   Optional: K_GLOBAL=4 (default)
#
# Per-population demographic Ne(t) inference with SMC++ v1.15.5.
# Complements Stairway Plot (SFS-based) and GONE (LD-based) to resolve
# discordances in recent Ne trajectories across the 4 ADMIXTURE clusters.
#
# smc++ cv is used instead of estimate: k-fold cross-validation
# automatically selects the regularization parameter (rho), avoiding
# both overfitting (too many breakpoints) and underfitting.
#
# Pipeline per population:
#   STEP 3 : Shared VCF extraction (155 FG individuals, 13 Super-Scaffolds)
#   STEP 4 : vcf2smc — one file per (distinguished_ind × scaffold)
#             NDIST=15 individuals used as distinguished (min(n, 15))
#             → 15 × 13 = 195 composite-likelihood input files per population
#             → Parallelised with GNU parallel (-j CPUS)
#   STEP 5 : smc++ cv — cross-validated Ne(t) estimation per population
#             --folds 3 | --knots 32 | --em-iterations 30
#             Populations run sequentially, --cores CPUS each (~6h/pop)
#   STEP 6 : smc++ plot --csv — export Ne(t) as CSV for R
#   STEP 7 : R — combined Ne(t) plot, mu sensitivity, Pop_3 vs Stairway
#
# Mask: TE regions (te_merged.bed) treated as missing in vcf2smc so
# that SMC++ does not count invariant TE sites as observable.
#
# Parameters:
#   mu  = 5e-9  per site per generation (consistent with Stairway Plot)
#   r   = 3e-8  per bp per generation (angiosperm average; no Dicorynia map)
#   Tg  = 30    years per generation
#   NDIST = 15  distinguished individuals (maximises sensitivity to recent Ne)
#   CV_FOLDS = 3 (good balance: 3 fits + 1 final per population)
#
# Output: 07-population_structure/workflow_full/07.12-smcpp/
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment
# -------------------------------------------------------------------
module purge
module load devel/Miniconda/Miniconda3
module load bioinfo/SMC++/1.15.5
module load bioinfo/Bcftools/1.9
module load bioinfo/samtools/1.21
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

CPUS="${SLURM_CPUS_PER_TASK:-16}"

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL="${K_GLOBAL:-4}"
K_COL="K${K_GLOBAL}"
MU_REF="5e-9"
RECOMB="3e-8"
YEAR_PER_GEN=30
NDIST=15
CV_FOLDS=3

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"
TE_BED="${BASE}/07-population_structure/workflow_full/07.11-stairway_plot/07.11b-compute_L/te_merged.bed"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.12a-smcpp_estimate.R"
STAIRWAY_BASE="${BASE}/07-population_structure/workflow_full/07.11-stairway_plot"
OUT_BASE="${BASE}/07-population_structure/workflow_full/07.12-smcpp"

# 13 chromosome-level Super-Scaffolds (≥30 Mb); same set as 07.8
SCAFFOLDS=(
    Super-Scaffold_1  Super-Scaffold_2  Super-Scaffold_4
    Super-Scaffold_5  Super-Scaffold_6  Super-Scaffold_7
    Super-Scaffold_8  Super-Scaffold_12 Super-Scaffold_13
    Super-Scaffold_14 Super-Scaffold_16 Super-Scaffold_17
    Super-Scaffold_21
)
SCAFFOLDS_CSV=$(IFS=,; echo "${SCAFFOLDS[*]}")

# Scaffold lengths from Dgu_HS1_HYBRID_SCAFFOLD.fa.fai
# Required by smc++ vcf2smc --length: subset VCF headers lack ##contig lines
declare -A SCAFFOLD_LENGTHS=(
    [Super-Scaffold_1]=41062991
    [Super-Scaffold_2]=43610879
    [Super-Scaffold_4]=38450403
    [Super-Scaffold_5]=40225412
    [Super-Scaffold_6]=37865487
    [Super-Scaffold_7]=37763175
    [Super-Scaffold_8]=39463349
    [Super-Scaffold_12]=40647968
    [Super-Scaffold_13]=41045391
    [Super-Scaffold_14]=34837072
    [Super-Scaffold_16]=32669108
    [Super-Scaffold_17]=36622932
    [Super-Scaffold_21]=33998326
)

LOG_DIR="${OUT_BASE}/logs"
mkdir -p "$LOG_DIR" "${OUT_BASE}/csv" "${OUT_BASE}/plots"
LOG_FILE="${LOG_DIR}/07.12a-smcpp_estimate.log"
echo "START $(date -Iseconds) 07.12a SMC++ estimate K=${K_GLOBAL}" > "$LOG_FILE"
echo "Input VCF  : $INPUT_VCF"     | tee -a "$LOG_FILE"
echo "TE mask    : $TE_BED"         | tee -a "$LOG_FILE"
echo "mu=${MU_REF} | r=${RECOMB} | Tg=${YEAR_PER_GEN} yr/gen" | tee -a "$LOG_FILE"
echo "NDIST=${NDIST} | CV_FOLDS=${CV_FOLDS} | CPUs=${CPUS}" | tee -a "$LOG_FILE"

for f in "$INPUT_VCF" "$TE_BED" "$CLUSTER_ASSIGNMENTS"; do
    [ -s "$f" ] || { echo "ERROR: not found: $f" | tee -a "$LOG_FILE"; exit 1; }
done

# -------------------------------------------------------------------
# STEP 2: Per-population sample lists from cluster_assignments.tsv
# -------------------------------------------------------------------
echo "STEP 2: building sample lists" | tee -a "$LOG_FILE"

awk -v kcol="$K_COL" '
NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == kcol) col = i
    if (!col) { print "ERROR: column " kcol " not found" > "/dev/stderr"; exit 2 }
    next
}
$col != "" && $col != "NA" { print $1 "\t" $col }
' "$CLUSTER_ASSIGNMENTS" > "${LOG_DIR}/sample_pop_map.tsv"

awk '{print $1}' "${LOG_DIR}/sample_pop_map.tsv" \
    > "${LOG_DIR}/fg_samples.txt"
N_FG=$(wc -l < "${LOG_DIR}/fg_samples.txt")
echo "  FG individuals: ${N_FG}" | tee -a "$LOG_FILE"

for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    awk -v p="$POP" '$2 == p { print $1 }' "${LOG_DIR}/sample_pop_map.tsv" \
        > "${LOG_DIR}/samples_${POP}.txt"
    N=$(wc -l < "${LOG_DIR}/samples_${POP}.txt")
    echo "  ${POP}: ${N} individuals" | tee -a "$LOG_FILE"
done

# -------------------------------------------------------------------
# STEP 2.5: bgzip + tabix-index TE mask for smc++ vcf2smc --mask
#
# smc++ vcf2smc passes the mask to pysam.TabixFile, which requires a
# bgzipped BED (.bed.gz) with a tabix index (.bed.gz.tbi).
# A plain .bed file causes: OSError: index .tbi not found
# -------------------------------------------------------------------
TE_BED_GZ="${TE_BED}.gz"
if [ ! -s "${TE_BED_GZ}.tbi" ]; then
    echo "STEP 2.5: bgzipping and indexing TE mask" | tee -a "$LOG_FILE"
    bgzip -c "$TE_BED" > "$TE_BED_GZ"
    tabix -p bed "$TE_BED_GZ"
    echo "  TE mask indexed: $TE_BED_GZ" | tee -a "$LOG_FILE"
else
    echo "STEP 2.5: TE mask index exists — reusing $TE_BED_GZ" | tee -a "$LOG_FILE"
fi
TE_BED="$TE_BED_GZ"

# -------------------------------------------------------------------
# STEP 3: Shared VCF — 155 FG individuals, 13 Super-Scaffolds
#         SMC++ uses Super-Scaffold names directly (no integer renaming).
#         All populations share this VCF; the per-population split
#         is handled at vcf2smc time via the POP:samples_csv argument.
# -------------------------------------------------------------------
echo "STEP 3: shared VCF extraction" | tee -a "$LOG_FILE"

FG_VCF="${LOG_DIR}/fg_155_chr13.vcf.gz"
if [ -s "$FG_VCF" ]; then
    echo "  reusing existing shared VCF" | tee -a "$LOG_FILE"
else
    echo "  extracting 155 individuals, 13 scaffolds..." | tee -a "$LOG_FILE"
    bcftools view \
        --samples-file "${LOG_DIR}/fg_samples.txt" \
        -r "${SCAFFOLDS_CSV}" \
        -O z -o "$FG_VCF" \
        "$INPUT_VCF" >> "$LOG_FILE" 2>&1
    bcftools index -t "$FG_VCF" >> "$LOG_FILE" 2>&1
    N_SNPS=$(bcftools index -n "$FG_VCF")
    echo "  INFO: ${N_SNPS} SNPs retained" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 4: vcf2smc — per (population × distinguished_individual × scaffold)
#
# For each population, the first min(n, NDIST) individuals are used
# as distinguished (one at a time). Each distinguished individual
# contributes one file per scaffold.
#
# Command:
#   smc++ vcf2smc --mask TE_BED -d DIST DIST VCF OUT SCAFFOLD POP:csv
#   -d DIST DIST  : same individual for both haplotypes (unphased mode)
#   POP:csv       : population label + comma-separated sample list
#
# All commands are written to a file and executed in parallel.
# -------------------------------------------------------------------
echo "STEP 4: generating vcf2smc commands" | tee -a "$LOG_FILE"

# Remove empty/truncated .smc.gz files from previous failed runs.
# Legitimate files are several MB; files from a run without --length are ~450 bytes.
SMC_MIN_BYTES=1000
N_PURGED=0
for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    OUT_POP="${OUT_BASE}/${POP}"
    while IFS= read -r f; do
        rm -f "$f"
        N_PURGED=$((N_PURGED + 1))
    done < <(find "${OUT_POP}/vcf2smc" -name "*.smc.gz" -size -${SMC_MIN_BYTES}c 2>/dev/null)
done
[ "$N_PURGED" -gt 0 ] && echo "  INFO: purged ${N_PURGED} truncated .smc.gz files" | tee -a "$LOG_FILE"

VCF2SMC_CMDS="${LOG_DIR}/vcf2smc_commands.txt"
> "$VCF2SMC_CMDS"
N_SKIP=0
N_GEN=0

for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    OUT_POP="${OUT_BASE}/${POP}"
    mkdir -p "${OUT_POP}/vcf2smc"

    SAMPLES_FILE="${LOG_DIR}/samples_${POP}.txt"
    N_IND=$(wc -l < "$SAMPLES_FILE")
    N_DIST=$(( N_IND < NDIST ? N_IND : NDIST ))
    SAMPLES_CSV=$(paste -sd',' "$SAMPLES_FILE")

    echo "  ${POP}: ${N_IND} ind, ${N_DIST} distinguished, ${SCAFFOLDS[*]}" | tee -a "$LOG_FILE"

    while IFS= read -r DIST; do
        for CHR in "${SCAFFOLDS[@]}"; do
            OUT_SMC="${OUT_POP}/vcf2smc/${POP}.${CHR}.${DIST}.smc.gz"
            if [ -s "$OUT_SMC" ]; then
                N_SKIP=$((N_SKIP + 1))
                continue
            fi
            printf 'smc++ vcf2smc --mask "%s" --length "%s" -d "%s" "%s" "%s" "%s" "%s" "%s:%s"\n' \
                "$TE_BED" "${SCAFFOLD_LENGTHS[$CHR]}" "$DIST" "$DIST" "$FG_VCF" \
                "$OUT_SMC" "$CHR" "$POP" "$SAMPLES_CSV" \
                >> "$VCF2SMC_CMDS"
            N_GEN=$((N_GEN + 1))
        done
    done < <(head -n "${N_DIST}" "$SAMPLES_FILE")
done

echo "  INFO: ${N_GEN} commands to run, ${N_SKIP} already done" | tee -a "$LOG_FILE"

if [ "$N_GEN" -gt 0 ]; then
    echo "  running vcf2smc in parallel (-j ${CPUS})..." | tee -a "$LOG_FILE"
    parallel --jobs "${CPUS}" < "$VCF2SMC_CMDS" >> "$LOG_FILE" 2>&1 \
        || echo "WARN: some vcf2smc jobs failed — check log" | tee -a "$LOG_FILE"
    echo "  vcf2smc done $(date -Iseconds)" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 5: smc++ cv — cross-validated Ne(t) estimation
#
# Populations run SEQUENTIALLY; each uses --cores CPUS for internal
# parallelism. The cv command:
#   1. Splits data files into CV_FOLDS folds
#   2. Fits CV_FOLDS models (train k-1, test 1) to select rho (regularisation)
#   3. Fits a final model on all data using the best rho
# Output: estimate/model.final.json (same format as smc++ estimate)
# -------------------------------------------------------------------
echo "STEP 5: smc++ cv (populations sequential)" | tee -a "$LOG_FILE"

for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    OUT_POP="${OUT_BASE}/${POP}"
    MODEL_FILE="${OUT_POP}/estimate/model.final.json"

    if [ -s "$MODEL_FILE" ]; then
        echo "  [${POP}] model.final.json already exists — skipping cv" | tee -a "$LOG_FILE"
        continue
    fi

    N_SMC=$(ls "${OUT_POP}/vcf2smc/"*.smc.gz 2>/dev/null | wc -l)
    if [ "$N_SMC" -eq 0 ]; then
        echo "  WARN: [${POP}] no vcf2smc files — skipping cv" | tee -a "$LOG_FILE"
        continue
    fi

    mkdir -p "${OUT_POP}/estimate"
    echo "  [${POP}] smc++ cv: ${N_SMC} files, folds=${CV_FOLDS}, cores=${CPUS} — $(date -Iseconds)" \
        | tee -a "$LOG_FILE"

    smc++ cv \
        --cores "${CPUS}" \
        --knots 32 \
        --em-iterations 30 \
        --folds "${CV_FOLDS}" \
        -o "${OUT_POP}/estimate" \
        "${MU_REF}" \
        "${OUT_POP}/vcf2smc/"*.smc.gz >> "$LOG_FILE" 2>&1 \
        || echo "  WARN: [${POP}] smc++ cv returned non-zero — check log" | tee -a "$LOG_FILE"

    if [ -s "$MODEL_FILE" ]; then
        echo "  [${POP}] cv done — model.final.json written $(date -Iseconds)" | tee -a "$LOG_FILE"
    else
        echo "  WARN: [${POP}] model.final.json not found after cv" | tee -a "$LOG_FILE"
    fi
done

# -------------------------------------------------------------------
# STEP 6: smc++ plot --csv — export Ne(t) trajectories for R
#         -g YEAR_PER_GEN: converts time axis from generations to years
#         --csv: writes <output>.csv instead of a PDF/PNG
#         One CSV per population; combined plot produced in R.
# -------------------------------------------------------------------
echo "STEP 6: smc++ plot --csv export" | tee -a "$LOG_FILE"

for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    MODEL_FILE="${OUT_BASE}/${POP}/estimate/model.final.json"
    CSV_PREFIX="${OUT_BASE}/csv/Ne_${POP}"
    CSV_FILE="${CSV_PREFIX}.csv"

    if [ ! -s "$MODEL_FILE" ]; then
        echo "  WARN: [${POP}] model not found — skipping plot export" | tee -a "$LOG_FILE"
        continue
    fi
    if [ -s "$CSV_FILE" ]; then
        echo "  [${POP}] CSV already exists — skipping" | tee -a "$LOG_FILE"
        continue
    fi

    smc++ plot \
        -g "${YEAR_PER_GEN}" \
        --csv \
        "${CSV_PREFIX}" \
        "${MODEL_FILE}" >> "$LOG_FILE" 2>&1 \
        || echo "  WARN: [${POP}] smc++ plot failed" | tee -a "$LOG_FILE"

    [ -s "$CSV_FILE" ] && echo "  [${POP}] CSV: ${CSV_FILE}" | tee -a "$LOG_FILE"
done

# -------------------------------------------------------------------
# STEP 7: R — Ne(t) visualisation
#   Figure 1: Combined Ne(t) for 4 populations (POP_COLORS, log-log)
#   Figure 2: Mu sensitivity per population (mu = 2–8 × 10^-9)
#   Figure 3: Pop_3 comparison SMC++ vs Stairway Plot 2 (log-log)
# -------------------------------------------------------------------
echo "STEP 7: R visualization" | tee -a "$LOG_FILE"

if [ -s "$R_SCRIPT" ]; then
    Rscript --vanilla "$R_SCRIPT" \
        "$OUT_BASE" \
        "$STAIRWAY_BASE" \
        "$K_GLOBAL" \
        "$MU_REF" \
        "$YEAR_PER_GEN" >> "$LOG_FILE" 2>&1 \
        || echo "WARN: R script returned non-zero — check log" | tee -a "$LOG_FILE"
else
    echo "WARN: R script not found: $R_SCRIPT" | tee -a "$LOG_FILE"
fi

echo "DONE $(date -Iseconds) 07.12a SMC++ estimate" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_BASE}/"

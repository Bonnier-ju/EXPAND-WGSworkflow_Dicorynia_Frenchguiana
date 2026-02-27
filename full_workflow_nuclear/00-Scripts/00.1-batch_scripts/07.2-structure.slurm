#!/bin/bash
#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 07.2-structure
#SBATCH --array=1-12
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "Starting 07.2 STRUCTURE runs"

module purge
module load bioinfo/PLINK/1.90b7
module load bioinfo/Structure/v2.3.4_console
module load statistics/R/4.3.0 >/dev/null 2>&1 || module load devel/R >/dev/null 2>&1 || true

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
IN_071_DIR="${BASE}/07-population_structure/workflow_test/07.1-pca_ibd_fst"
OUT_DIR="${BASE}/07-population_structure/workflow_test/07.2-structure"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
LOG_FILE="${LOG_DIR}/07.2-structure_task${TASK_ID}.log"
SUMMARY_FILE="${OUT_DIR}/07.2-structure_runs_task${TASK_ID}.tsv"

BFILE_PREFIX="${IN_071_DIR}/genotypes"
STRUCT_PREFIX="${OUT_DIR}/structure_input"
STRUCT_IN_DEFAULT="${STRUCT_PREFIX}.recode.strct_in"
STRUCT_IN="${STRUCT_IN:-$STRUCT_IN_DEFAULT}"

MAIN_PARA="${OUT_DIR}/main_param.file"
EXTRA_PARA="${OUT_DIR}/extra_param.file"
R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.2-structure_barplots.R"
FAM_FILE="${IN_071_DIR}/genotypes.fam"
POP_MAP_FILE="${IN_071_DIR}/sample_population_map.tsv"

# STRUCTURE run design
K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-12}"
N_REPS="${N_REPS:-3}"
BASE_SEED="${BASE_SEED:-1000}"
BURNIN="${BURNIN:-50000}"
NUMREPS="${NUMREPS:-500000}"

echo "START $(date -Iseconds) 07.2 STRUCTURE" > "$LOG_FILE"
echo "Input 07.1 bed prefix: $BFILE_PREFIX" | tee -a "$LOG_FILE"
echo "Output dir: $OUT_DIR" | tee -a "$LOG_FILE"
echo "K range: ${K_MIN}-${K_MAX}; replicates: ${N_REPS}" | tee -a "$LOG_FILE"
echo "BURNIN=$BURNIN; NUMREPS=$NUMREPS" | tee -a "$LOG_FILE"
echo "Array task: $TASK_ID" | tee -a "$LOG_FILE"

if [ ! -s "${BFILE_PREFIX}.bed" ] || [ ! -s "${BFILE_PREFIX}.bim" ] || [ ! -s "${BFILE_PREFIX}.fam" ]; then
    echo "ERROR: missing 07.1 PLINK files (${BFILE_PREFIX}.bed/.bim/.fam)" | tee -a "$LOG_FILE"
    exit 1
fi

# Build STRUCTURE input from 07.1 PLINK files if not provided
if [ ! -s "$STRUCT_IN" ]; then
    echo "Creating STRUCTURE input with PLINK --recode structure" | tee -a "$LOG_FILE"
    plink \
        --bfile "$BFILE_PREFIX" \
        --allow-extra-chr \
        --recode structure \
        --out "$STRUCT_PREFIX" >> "$LOG_FILE" 2>&1
fi

if [ ! -s "$STRUCT_IN" ]; then
    echo "ERROR: STRUCTURE input not found: $STRUCT_IN" | tee -a "$LOG_FILE"
    exit 1
fi

NUMINDS=$(wc -l < "${BFILE_PREFIX}.fam")
NUMLOCI=$(wc -l < "${BFILE_PREFIX}.bim")

# Auto-generate parameter files (overwrite by default to avoid stale mismatches)
if [ "${USE_EXISTING_PARAMS:-0}" != "1" ] || [ ! -s "$MAIN_PARA" ]; then
    cat > "$MAIN_PARA" <<EOF
#define INFILE $STRUCT_IN
#define OUTFILE ${OUT_DIR}/str
#define NUMINDS $NUMINDS
#define NUMLOCI $NUMLOCI
#define PLOIDY 2
#define MISSING -9
#define ONEROWPERIND 1
#define LABEL 1
#define POPDATA 1
#define POPFLAG 0
#define PHENOTYPE 0
#define EXTRACOLS 0
#define MARKERNAMES 1
#define RECESSIVEALLELES 0
#define MAPDISTANCES 1
#define BURNIN $BURNIN
#define NUMREPS $NUMREPS
EOF
fi

if [ "${USE_EXISTING_PARAMS:-0}" != "1" ] || [ ! -s "$EXTRA_PARA" ]; then
    cat > "$EXTRA_PARA" <<'EOF'
#define NOADMIX 0
#define LINKAGE 0
#define LOCPRIOR 0
#define FREQSCORR 1
#define ONEFST 0
#define INFERALPHA 1
#define ALPHA 1.0
#define INFERLAMBDA 0
#define LAMBDA 1.0
#define POPALPHAS 0
#define COMPUTEPROB 1
#define PRINTKLD 1
#define PRINTLIKES 1
#define RANDOMIZE 0
EOF
fi

echo -e "k\treplicate\tseed\toutput_prefix" > "$SUMMARY_FILE"

K=$((K_MIN + TASK_ID - 1))
if [ "$K" -lt "$K_MIN" ] || [ "$K" -gt "$K_MAX" ]; then
    echo "INFO: task ${TASK_ID} mapped to K=${K}, out of range ${K_MIN}-${K_MAX}. Exiting." | tee -a "$LOG_FILE"
    exit 0
fi

for REP in $(seq 1 "$N_REPS"); do
    SEED=$((BASE_SEED + K * 100 + REP))
    OUT_PREFIX="${OUT_DIR}/str_k${K}_r${REP}"
    echo "Running STRUCTURE K=${K} REP=${REP} SEED=${SEED}" | tee -a "$LOG_FILE"
    structure \
        -i "$STRUCT_IN" \
        -m "$MAIN_PARA" \
        -e "$EXTRA_PARA" \
        -K "$K" \
        -D "$SEED" \
        -o "$OUT_PREFIX" >> "$LOG_FILE" 2>&1
    echo -e "${K}\t${REP}\t${SEED}\t${OUT_PREFIX}" >> "$SUMMARY_FILE"
done

echo "DONE $(date -Iseconds) 07.2 STRUCTURE completed for K=${K}" | tee -a "$LOG_FILE"
echo "Run summary: $SUMMARY_FILE"

# Optional: auto-generate STRUCTURE barplots when all expected runs are present
if [ "${RUN_R_PLOT:-1}" = "1" ]; then
    expected_runs=$(( (K_MAX - K_MIN + 1) * N_REPS ))
    found_runs=$(find "$OUT_DIR" -maxdepth 1 -type f -name 'str_k*_r*_f' | wc -l)
    echo "R plot check: found_runs=${found_runs}; expected_runs=${expected_runs}" | tee -a "$LOG_FILE"
    if [ "$found_runs" -ge "$expected_runs" ]; then
        if [ -s "$R_PLOT_SCRIPT" ] && command -v Rscript >/dev/null 2>&1; then
            LOCK_DIR="${OUT_DIR}/.plot_lock"
            if mkdir "$LOCK_DIR" 2>/dev/null; then
                echo "STEP: STRUCTURE barplot generation in R" | tee -a "$LOG_FILE"
                Rscript --vanilla "$R_PLOT_SCRIPT" "$OUT_DIR" "$FAM_FILE" "$POP_MAP_FILE" >> "$LOG_FILE" 2>&1 || \
                    echo "WARN: STRUCTURE R plotting failed" | tee -a "$LOG_FILE"
                rmdir "$LOCK_DIR" || true
            else
                echo "INFO: plot step already handled by another task" | tee -a "$LOG_FILE"
            fi
        else
            echo "WARN: R plotting skipped (Rscript or 07.2 R script missing)" | tee -a "$LOG_FILE"
        fi
    else
        echo "INFO: plotting deferred (not all K/rep outputs available yet)" | tee -a "$LOG_FILE"
    fi
fi

#!/bin/bash
#SBATCH -p workq
#SBATCH --time=48:00:00
#SBATCH -J 07.6-admixture_cluster
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   K_GLOBAL=2 POP_ID=1 sbatch 07.6-admixture_cluster.slurm
#   K_GLOBAL=3 POP_ID=2 sbatch 07.6-admixture_cluster.slurm
#
# Parameters:
#   K_GLOBAL  — global ADMIXTURE K used to define populations (2, 3 or 4)
#   POP_ID    — population index to analyse (1 .. n_pop for this K)
#   K_MIN     — minimum K for intra-cluster run (default 1)
#   K_MAX     — maximum K for intra-cluster run (default 5)
#   CV_FOLDS  — cross-validation folds (default 10)
#   N_REP     — replicates per K (default 20)
#
# Runs ADMIXTURE K=K_MIN..K_MAX on the subset of guiana_only individuals
# belonging to Pop_${POP_ID} under the K=${K_GLOBAL} grouping.
# For each K, N_REP replicates are run in parallel (1 thread each, seed
# = K*1000+rep for reproducibility). The best replicate by log-likelihood
# is retained; all replicate logs are kept under replicates/.
# Output: 07.2-admixture/07.2.${K_GLOBAL}-K=${K_GLOBAL}/Pop_${POP_ID}/
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/1.90b7
module load bioinfo/ADMIXTURE/1.3.0
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL="${K_GLOBAL:-4}"
POP_ID="${POP_ID:-4}"
POP_LABEL="Pop_${POP_ID}"
K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-5}"
CV_FOLDS="${CV_FOLDS:-10}"
N_REP="${N_REP:-20}"
MAX_PARALLEL="${SLURM_CPUS_PER_TASK:-8}"  # run MAX_PARALLEL replicates concurrently, 1 thread each

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
VCF_66="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
VCF_65="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"
SITE_MAP="${BASE}/07-population_structure/workflow_full/07.1-pca_ibd_fst/07.1.1-K=1/sample_population_map.tsv"

R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.2-admixture_plot.R"

if [ -s "$VCF_66" ]; then
    INPUT_VCF="$VCF_66"; INPUT_STEP="06.6b"
else
    INPUT_VCF="$VCF_65"; INPUT_STEP="06.5"
fi

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.2-admixture/07.2.${K_GLOBAL}-K=${K_GLOBAL}/${POP_LABEL}"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
REP_DIR="${OUT_DIR}/replicates"
LOG_FILE="${LOG_DIR}/07.6-admixture_K${K_GLOBAL}_${POP_LABEL}.log"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR" "$REP_DIR"
echo "START $(date -Iseconds) 07.6 ADMIXTURE cluster K=${K_GLOBAL} ${POP_LABEL}" > "$LOG_FILE"
echo "Input VCF : $INPUT_VCF ($INPUT_STEP)" | tee -a "$LOG_FILE"
echo "Intra-cluster K range: ${K_MIN}..${K_MAX} | CV folds: ${CV_FOLDS} | reps per K: ${N_REP} | parallel: ${MAX_PARALLEL}" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"; exit 1
fi
if [ ! -s "$CLUSTER_ASSIGNMENTS" ]; then
    echo "ERROR: cluster_assignments.tsv not found — run 07.0-prep_cluster_assignments.R first" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build keep file for this K_GLOBAL / POP_LABEL
# Keep file uses FID=0 (--const-fid convention for ADMIXTURE)
# -------------------------------------------------------------------
K_COL="K${K_GLOBAL}"
KEEP_FILE="${TMP_DIR}/keep_${POP_LABEL}.txt"

awk -v kcol="$K_COL" -v pop="$POP_LABEL" '
NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == kcol) col = i
    if (!col) { print "ERROR: column " kcol " not found" > "/dev/stderr"; exit 2 }
    next
}
$col == pop { print "0\t" $1 }
' "$CLUSTER_ASSIGNMENTS" > "$KEEP_FILE"

N_IND=$(wc -l < "$KEEP_FILE")
echo "Individuals in K=${K_GLOBAL} ${POP_LABEL}: ${N_IND}" | tee -a "$LOG_FILE"

if [ "$N_IND" -lt 3 ]; then
    echo "ERROR: fewer than 3 individuals — cannot run ADMIXTURE" | tee -a "$LOG_FILE"; exit 1
fi

# Cap K_MAX to avoid degenerate Q matrices: ADMIXTURE is unreliable when K > N/3
K_MAX_EFF=$(( N_IND / 3 ))
[ "$K_MAX_EFF" -lt 1 ] && K_MAX_EFF=1
[ "$K_MAX_EFF" -gt "$K_MAX" ] && K_MAX_EFF="$K_MAX"
echo "INFO: effective K_MAX = ${K_MAX_EFF} (N_IND=${N_IND}, requested K_MAX=${K_MAX})" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: PLINK VCF -> BED (subset to this population, done once)
# -------------------------------------------------------------------
PLINK_PREFIX="${OUT_DIR}/admixture_input"
echo "STEP 3: PLINK conversion" | tee -a "$LOG_FILE"
plink \
    --vcf "$INPUT_VCF" \
    --make-bed \
    --allow-extra-chr \
    --const-fid \
    --keep "$KEEP_FILE" \
    --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "${PLINK_PREFIX}.bed" ]; then
    echo "ERROR: PLINK conversion failed" | tee -a "$LOG_FILE"; exit 1
fi
N_SNP=$(wc -l < "${PLINK_PREFIX}.bim")
echo "INFO: ${N_IND} individuals, ${N_SNP} SNPs" | tee -a "$LOG_FILE"

# ADMIXTURE requires numeric chromosome IDs — remap all to 0
awk '{$1=0; print $0}' "${PLINK_PREFIX}.bim" > "${PLINK_PREFIX}_tmp.bim"
mv "${PLINK_PREFIX}_tmp.bim" "${PLINK_PREFIX}.bim"

# -------------------------------------------------------------------
# STEP 4: ADMIXTURE — N_REP replicates per K, parallel, best-by-loglik
#
# Each replicate runs in its own subdirectory (replicates/KX/repY/) so
# ADMIXTURE output files do not collide when multiple instances run
# concurrently. Seed = K*1000+rep for full reproducibility.
# MAX_PARALLEL replicates run at once, each using 1 thread, keeping
# total CPU usage equal to the SLURM allocation while maximising
# throughput via parallelism over sequential single-threaded runs.
# -------------------------------------------------------------------
ABS_BED="$(realpath "${PLINK_PREFIX}.bed")"
BED_BASE="$(basename "${PLINK_PREFIX}")"

> "${OUT_DIR}/cross_validation.txt"   # reset CV summary

echo "STEP 4: ADMIXTURE K=${K_MIN}..${K_MAX_EFF} x ${N_REP} replicates" | tee -a "$LOG_FILE"

for K in $(seq "$K_MIN" "$K_MAX_EFF"); do
    echo "  K=${K}: launching ${N_REP} replicates..." | tee -a "$LOG_FILE"

    for REP in $(seq 1 "$N_REP"); do
        SEED=$(( K * 1000 + REP ))
        REP_OUT="${REP_DIR}/K${K}/rep${REP}"
        mkdir -p "$REP_OUT"

        (
            cd "$REP_OUT"
            admixture \
                --cv="${CV_FOLDS}" \
                -j1 \
                --seed="${SEED}" \
                "${ABS_BED}" "${K}" \
                > "${REP_OUT}/log.out" 2>&1
        ) &

        # Throttle: hold when MAX_PARALLEL jobs are already running
        while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]; do
            sleep 2
        done
    done
    wait   # ensure all replicates for this K finish before selecting best

    echo "    K=${K}: all replicates finished, selecting best by log-likelihood" | tee -a "$LOG_FILE"

    # --- Select best replicate by highest log-likelihood ---
    LOGLIK_SUMMARY="${OUT_DIR}/loglik_K${K}.txt"
    echo "# K=${K} log-likelihoods (seed = K*1000+rep)" > "$LOGLIK_SUMMARY"

    best_rep=1
    best_ll="-999999999"

    for REP in $(seq 1 "$N_REP"); do
        REP_OUT="${REP_DIR}/K${K}/rep${REP}"
        ll=$(grep "Loglikelihood:" "${REP_OUT}/log.out" 2>/dev/null \
             | awk '{print $2}' | tail -1)
        [ -z "$ll" ] && ll="-999999999"
        printf "  rep%02d  seed=%-6d  LL=%s\n" "$REP" "$(( K * 1000 + REP ))" "$ll" \
            | tee -a "$LOGLIK_SUMMARY"

        is_better=$(awk -v a="$ll" -v b="$best_ll" 'BEGIN { print (a+0 > b+0) ? 1 : 0 }')
        if [ "$is_better" = "1" ]; then
            best_ll="$ll"
            best_rep="$REP"
        fi
    done

    echo "  best: rep${best_rep}  LL=${best_ll}" | tee -a "$LOGLIK_SUMMARY"
    echo "    K=${K} → best=rep${best_rep}, LL=${best_ll}" | tee -a "$LOG_FILE"
    cat "$LOGLIK_SUMMARY" >> "$LOG_FILE"

    # Copy best Q and P to OUT_DIR
    BEST_Q="${REP_DIR}/K${K}/rep${best_rep}/${BED_BASE}.${K}.Q"
    BEST_P="${REP_DIR}/K${K}/rep${best_rep}/${BED_BASE}.${K}.P"

    if [ -s "$BEST_Q" ]; then
        cp "$BEST_Q" "${OUT_DIR}/${BED_BASE}.${K}.Q"
        cp "$BEST_P" "${OUT_DIR}/${BED_BASE}.${K}.P" 2>/dev/null || true
    else
        echo "    WARN: Q file missing for K=${K} rep${best_rep}" | tee -a "$LOG_FILE"
    fi

    # Append CV error from best replicate to summary
    grep "CV error" "${REP_DIR}/K${K}/rep${best_rep}/log.out" 2>/dev/null \
        >> "${OUT_DIR}/cross_validation.txt" || true

done

echo "CV errors (best replicate per K):" | tee -a "$LOG_FILE"
cat "${OUT_DIR}/cross_validation.txt" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 5: Attach individual names to Q matrices
# -------------------------------------------------------------------
IDS_FILE="${OUT_DIR}/individual_ids.txt"
awk '{print $2}' "${PLINK_PREFIX}.fam" > "$IDS_FILE"

for K in $(seq "$K_MIN" "$K_MAX_EFF"); do
    Q_FILE="${OUT_DIR}/${BED_BASE}.${K}.Q"
    NAMED_Q="${OUT_DIR}/${BED_BASE}.${K}.named.Q"
    if [ -s "$Q_FILE" ]; then
        paste "$IDS_FILE" "$Q_FILE" > "$NAMED_Q"
    else
        echo "WARN: Q file missing for K=${K}" | tee -a "$LOG_FILE"
    fi
done

# -------------------------------------------------------------------
# STEP 6: R plotting (reuses existing 07.2-admixture_plot.R)
#   arg1 = out_dir   (directory with named.Q and cross_validation.txt)
#   arg2 = map_file  (sample_id -> site, for labels in barplot)
#   arg3 = k_min
#   arg4 = k_max
#   arg5 = bed_base  (prefix of Q files)
# -------------------------------------------------------------------
if [ -s "$R_PLOT_SCRIPT" ]; then
    echo "STEP 6: R plotting" | tee -a "$LOG_FILE"
    Rscript --vanilla "$R_PLOT_SCRIPT" \
        "$OUT_DIR" \
        "$SITE_MAP" \
        "$K_MIN" \
        "$K_MAX_EFF" \
        "$BED_BASE" >> "$LOG_FILE" 2>&1 \
        || echo "WARN: R plotting failed" | tee -a "$LOG_FILE"
else
    echo "WARN: R plot script not found: $R_PLOT_SCRIPT" | tee -a "$LOG_FILE"
fi

echo "DONE $(date -Iseconds) 07.6 ADMIXTURE K=${K_GLOBAL} ${POP_LABEL}" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}/"

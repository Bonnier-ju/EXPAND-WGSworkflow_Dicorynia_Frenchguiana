#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=180:00:00
#SBATCH -J 07.9-fstruct_cluster
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   K_GLOBAL=2 POP_ID=1 sbatch 07.9-fastStructure_cluster.slurm
#   K_GLOBAL=3 POP_ID=2 sbatch 07.9-fastStructure_cluster.slurm
#   K_GLOBAL=4 POP_ID=3 sbatch 07.9-fastStructure_cluster.slurm
#
# Parameters:
#   K_GLOBAL   — global ADMIXTURE K used to define populations (2, 3 or 4)
#   POP_ID     — population index to analyse (1 .. n_pop for this K)
#   K_MIN      — minimum K for intra-cluster run (default 1)
#   K_MAX      — maximum K for intra-cluster run (default 5)
#   SEEDS      — space-separated list of random seeds (default "42 123 999")
#   N_SNP_THIN — number of SNPs to retain after subsampling (default 50000)
#   THIN_SEED  — PLINK random seed for reproducible thinning (default 2024)
#
# Runs fastStructure K=K_MIN..K_MAX × SEEDS replicates on the subset of
# guiana_only individuals belonging to Pop_${POP_ID} under K=${K_GLOBAL}.
# Input SNPs are subsampled to N_SNP_THIN (from the LD-pruned VCF) for
# computational tractability; 50k LD-pruned SNPs provide sufficient power
# for structure inference in n=16-90 individuals (Raj et al. 2014).
# Seeds within each K run in parallel (one per CPU, up to cpus-per-task).
# Each K × seed is an independent run (--prior=logistic, --cv=10, --tol=1e-6).
# The replicate with the highest marginal likelihood is retained per K.
# K selection uses chooseK.py on the best-replicate set.
#
# Output: 07.6-fastStructure/07.6.${K_GLOBAL}-K=${K_GLOBAL}/Pop_${POP_ID}/
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/1.90b7
module load bioinfo/fastStructure/1.0
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL="${K_GLOBAL:-4}"
POP_ID="${POP_ID:-1}"
POP_LABEL="Pop_${POP_ID}"
K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-5}"
SEEDS="${SEEDS:-42 123 999}"
N_SNP_THIN="${N_SNP_THIN:-50000}"
THIN_SEED="${THIN_SEED:-2024}"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"
SITE_MAP="${BASE}/07-population_structure/workflow_full/07.1-pca_ibd_fst/07.1.1-K=1/sample_population_map.tsv"
R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.6-fastStructure_plot.R"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.6-fastStructure/07.6.${K_GLOBAL}-K=${K_GLOBAL}/${POP_LABEL}"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/07.9-fstruct_K${K_GLOBAL}_${POP_LABEL}.log"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 07.9 fastStructure cluster K=${K_GLOBAL} ${POP_LABEL}" > "$LOG_FILE"
echo "Input VCF : $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Intra-cluster K range: ${K_MIN}..${K_MAX} | seeds: ${SEEDS} | SNPs retained: ${N_SNP_THIN} (thin seed: ${THIN_SEED})" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"; exit 1
fi
if [ ! -s "$CLUSTER_ASSIGNMENTS" ]; then
    echo "ERROR: cluster_assignments.tsv not found — run 07.0-prep_cluster_assignments.R first" \
        | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build keep file for this K_GLOBAL / POP_LABEL
# Keep file uses FID=0 (--const-fid convention)
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
    echo "ERROR: fewer than 3 individuals — cannot run fastStructure" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 3: PLINK VCF -> BED (subset to this population)
# fastStructure accepts PLINK BED; no chromosome remapping needed
# (unlike ADMIXTURE, fastStructure tolerates non-numeric chr IDs).
# -------------------------------------------------------------------
PLINK_PREFIX="${OUT_DIR}/fastStructure_input"
echo "STEP 3: PLINK conversion" | tee -a "$LOG_FILE"

plink \
    --vcf "$INPUT_VCF" \
    --make-bed \
    --allow-extra-chr \
    --const-fid \
    --keep "$KEEP_FILE" \
    --thin-count "$N_SNP_THIN" \
    --seed "$THIN_SEED" \
    --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "${PLINK_PREFIX}.bed" ]; then
    echo "ERROR: PLINK conversion failed" | tee -a "$LOG_FILE"; exit 1
fi
N_SNP=$(wc -l < "${PLINK_PREFIX}.bim")
echo "INFO: ${N_IND} individuals, ${N_SNP} SNPs (subsampled from LD-pruned VCF, thin seed=${THIN_SEED})" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: fastStructure K=K_MIN..K_MAX × SEEDS replicates
# Each replicate uses a unique output prefix: run_K${K}_s${SEED}
# Parameters: logistic prior, 10-fold CV, convergence tol=1e-6
# -------------------------------------------------------------------
echo "STEP 4: fastStructure K=${K_MIN}..${K_MAX} x seeds [${SEEDS}]" | tee -a "$LOG_FILE"

BEST_SEEDS_FILE="${OUT_DIR}/best_seeds.txt"
> "$BEST_SEEDS_FILE"

for K in $(seq "$K_MIN" "$K_MAX"); do
    echo "  --- K=${K} ---" | tee -a "$LOG_FILE"

    # Launch all seeds in parallel (one background process per seed)
    pids=()
    for SEED in $SEEDS; do
        echo "    K=${K} seed=${SEED}: launching..." | tee -a "$LOG_FILE"
        seed_log="${TMP_DIR}/run_K${K}_s${SEED}.log"
        (
            python "${HOME_FAST}/structure.py" \
                -K "${K}" \
                --input="${PLINK_PREFIX}" \
                --output="${OUT_DIR}/run_K${K}_s${SEED}" \
                --prior=logistic \
                --cv=10 \
                --tol=1e-6 \
                --format=bed \
                --seed="${SEED}" \
                > "$seed_log" 2>&1
        ) &
        pids+=($!)
    done

    # Wait for all seeds; tally failures
    n_failed=0
    for pid in "${pids[@]}"; do
        wait "$pid" || n_failed=$((n_failed + 1))
    done

    # Append individual seed logs to main log in seed order
    for SEED in $SEEDS; do
        seed_log="${TMP_DIR}/run_K${K}_s${SEED}.log"
        [ -s "$seed_log" ] && cat "$seed_log" >> "$LOG_FILE"
    done
    [ "$n_failed" -gt 0 ] && \
        echo "  WARN: ${n_failed} seed(s) failed for K=${K}" | tee -a "$LOG_FILE"

    echo "  K=${K}: all seeds done" | tee -a "$LOG_FILE"

    # Select best replicate: highest marginal likelihood from per-run .log files
    best_seed=""
    best_marg=""
    for SEED in $SEEDS; do
        run_log="${OUT_DIR}/run_K${K}_s${SEED}.${K}.log"
        if [ -s "$run_log" ]; then
            marg=$(awk '/Marginal likelihood/ {print $NF}' "$run_log" | tail -1)
            if [ -n "$marg" ]; then
                if [ -z "$best_marg" ] || \
                   awk -v a="$marg" -v b="$best_marg" 'BEGIN {exit !(a > b)}'; then
                    best_marg="$marg"
                    best_seed="$SEED"
                fi
            fi
        fi
    done

    if [ -z "$best_seed" ]; then
        best_seed=$(echo "$SEEDS" | awk '{print $1}')
        echo "  WARN: could not compare marginal likelihoods for K=${K} — using first seed ${best_seed}" \
            | tee -a "$LOG_FILE"
    fi

    echo "  K=${K}: best replicate = seed ${best_seed} (marginal=${best_marg})" | tee -a "$LOG_FILE"
    echo "${K} ${best_seed} ${best_marg}" >> "$BEST_SEEDS_FILE"

    # Symlink best meanQ for chooseK.py (expects prefix.K.meanQ pattern)
    ln -sf "${OUT_DIR}/run_K${K}_s${best_seed}.${K}.meanQ" \
           "${OUT_DIR}/fastStructure_best.${K}.meanQ"
done

# K model selection on best replicates
echo "STEP 4b: chooseK.py on best replicates" | tee -a "$LOG_FILE"
python "${HOME_FAST}/chooseK.py" \
    --input="${OUT_DIR}/fastStructure_best" \
    >> "$LOG_FILE" 2>&1 \
    || echo "WARN: chooseK.py returned non-zero" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 5: Attach individual names to .meanQ matrices
# The resulting .named.Q files follow the same format as ADMIXTURE
# output, enabling direct reuse of 07.2-admixture_plot.R.
# -------------------------------------------------------------------
IDS_FILE="${OUT_DIR}/individual_ids.txt"
awk '{print $2}' "${PLINK_PREFIX}.fam" > "$IDS_FILE"

echo "STEP 5: building named.Q files from best-replicate .meanQ" | tee -a "$LOG_FILE"
for K in $(seq "$K_MIN" "$K_MAX"); do
    best_seed=$(awk -v k="$K" '$1 == k {print $2}' "$BEST_SEEDS_FILE")
    MEANQ_FILE="${OUT_DIR}/run_K${K}_s${best_seed}.${K}.meanQ"
    NAMED_Q="${OUT_DIR}/fastStructure_input.${K}.named.Q"
    if [ -s "$MEANQ_FILE" ]; then
        paste "$IDS_FILE" "$MEANQ_FILE" > "$NAMED_Q"
        echo "  K=${K}: named.Q from seed ${best_seed} ($(wc -l < "$NAMED_Q") rows)" \
            | tee -a "$LOG_FILE"
    else
        echo "WARN: meanQ missing for K=${K} seed=${best_seed}" | tee -a "$LOG_FILE"
    fi
done

# -------------------------------------------------------------------
# STEP 6: R barplot (reuses 07.2-admixture_plot.R)
#   arg1 = out_dir   (directory with named.Q files)
#   arg2 = map_file  (sample_id -> site, for individual labels)
#   arg3 = k_min
#   arg4 = k_max
#   arg5 = bed_base  (prefix shared by all named.Q files)
# -------------------------------------------------------------------
if [ -s "$R_PLOT_SCRIPT" ]; then
    echo "STEP 6: R barplot" | tee -a "$LOG_FILE"
    Rscript --vanilla "$R_PLOT_SCRIPT" \
        "$OUT_DIR" \
        "$SITE_MAP" \
        "$K_MIN" \
        "$K_MAX" \
        "fastStructure_input" \
        "$LOG_FILE" >> "$LOG_FILE" 2>&1 \
        || echo "WARN: R plotting failed" | tee -a "$LOG_FILE"
else
    echo "WARN: R plot script not found: $R_PLOT_SCRIPT" | tee -a "$LOG_FILE"
fi

echo "DONE $(date -Iseconds) 07.9 fastStructure K=${K_GLOBAL} ${POP_LABEL}" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}/plots/"

#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=160:00:00
#SBATCH -J 07b.2-admixture
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=48G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment
# -------------------------------------------------------------------
module purge
module load bioinfo/PLINK/1.90b7
module load bioinfo/ADMIXTURE/1.3.0
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Paths and parameters
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
WD="${BASE}/07b-population_structure_pur_inds"

VCF_LD="${WD}/input/snps_TE_filtered.ld_pruned.vcf.gz"
SAMPLE_SHEET="${BASE}/metadata/sample_sheet_complete_filtered.csv"
LISTS_DIR="${WD}/07b.0-purity_filtering/individual_lists"
OUT_BASE="${WD}/07b.2-admixture"
LOG_FILE="${WD}/logs/07b.2-admixture.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.2-admixture_plot.R"

K_MIN=1
K_MAX=6
N_REP=10
CV_FOLDS=100
THREADS_PER_RUN=3   # 10 reps × 3 threads = 30/32 CPUs in parallel

for f in "$VCF_LD" "$SAMPLE_SHEET" "$R_SCRIPT"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: missing required input: $f"
        exit 1
    fi
done

mkdir -p "$OUT_BASE" "${WD}/logs"
echo "START $(date -Iseconds) 07b.2 ADMIXTURE on purity-filtered individuals" > "$LOG_FILE"
echo "Parameters: K=${K_MIN}..${K_MAX} | reps=${N_REP} | CV=${CV_FOLDS} | threads/run=${THREADS_PER_RUN} | CPUs=32" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 1b: Build site-exclusion list (St_georges + Petit_Croissant)
# -------------------------------------------------------------------
EXCLUDE_FILE="${OUT_BASE}/exclude_sites.txt"
awk -F',' 'NR>1 && ($4 == "St_georges" || $4 == "Petit_Croissant") {print $1}' \
    "$SAMPLE_SHEET" > "$EXCLUDE_FILE"
N_EXCL=$(wc -l < "$EXCLUDE_FILE")
echo "Site exclusion: $N_EXCL individuals from St_georges + Petit_Croissant" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Loop over purity thresholds
# -------------------------------------------------------------------
for THRESH in 70 80 90 95; do
    TAG="T${THRESH}"
    KEEP_RAW="${LISTS_DIR}/inds_${TAG}_all.txt"

    if [ ! -s "$KEEP_RAW" ]; then
        echo "WARN: keep file not found for $TAG — skipping" | tee -a "$LOG_FILE"
        continue
    fi

    OUT_T="${OUT_BASE}/${TAG}"
    mkdir -p "${OUT_T}/reps"

    # Skip threshold if R plots already exist (threshold fully complete)
    if ls "${OUT_T}/plots/"*.png 2>/dev/null | grep -q .; then
        echo "=== $TAG: already completed — skipping ===" | tee -a "$LOG_FILE"
        continue
    fi

    # STEP 2a: Apply site exclusion — keep FID=0 format for --const-fid
    KEEP_FILE="${OUT_T}/keep_${TAG}.txt"
    awk 'NR==FNR {excl[$1]=1; next} !($2 in excl)' \
        "$EXCLUDE_FILE" "$KEEP_RAW" > "$KEEP_FILE"

    N_IND=$(wc -l < "$KEEP_FILE")
    echo "=== $TAG: $N_IND individuals (after site exclusion) ===" | tee -a "$LOG_FILE"

    # STEP 2b: PLINK VCF → BED (skip if already done)
    PLINK_PREFIX="${OUT_T}/admixture_input"
    if [ ! -s "${PLINK_PREFIX}.bed" ]; then
        plink \
            --vcf "$VCF_LD" \
            --make-bed \
            --allow-extra-chr \
            --const-fid \
            --keep "$KEEP_FILE" \
            --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

        if [ ! -s "${PLINK_PREFIX}.bed" ]; then
            echo "ERROR: PLINK conversion failed for $TAG" | tee -a "$LOG_FILE"
            continue
        fi

        # STEP 2c: Rename chromosomes to 0 (ADMIXTURE requires numeric chr IDs)
        awk '{$1=0; print $0}' "${PLINK_PREFIX}.bim" > "${PLINK_PREFIX}.bim.tmp"
        mv "${PLINK_PREFIX}.bim.tmp" "${PLINK_PREFIX}.bim"
    fi

    N_SNP=$(wc -l < "${PLINK_PREFIX}.bim")
    echo "  $N_IND individuals | $N_SNP SNPs" | tee -a "$LOG_FILE"

    # Initialize output files only if they don't exist (preserve entries from prior runs)
    [ -s "${OUT_T}/cross_validation.txt" ] || > "${OUT_T}/cross_validation.txt"
    [ -s "${OUT_T}/loglik_summary.txt" ]   || printf "K\trep\tloglikelihood\tcv_error\n" > "${OUT_T}/loglik_summary.txt"

    # -------------------------------------------------------------------
    # STEP 3: ADMIXTURE K=K_MIN..K_MAX — N_REP replicates in parallel
    # -------------------------------------------------------------------
    for K in $(seq "$K_MIN" "$K_MAX"); do
        K_DIR="${OUT_T}/reps/K${K}"
        mkdir -p "$K_DIR"

        # Skip K if best Q already copied to threshold directory
        if [ -s "${OUT_T}/admixture_input.${K}.Q" ]; then
            echo "  K=${K}: best Q already exists — skipping" | tee -a "$LOG_FILE"
            continue
        fi

        echo "  K=${K}: launching $N_REP replicates (${THREADS_PER_RUN} threads each)..." | tee -a "$LOG_FILE"

        # Launch all replicates simultaneously
        pids=()
        for rep in $(seq -w 1 "$N_REP"); do
            REP_DIR="${K_DIR}/rep${rep}"
            mkdir -p "$REP_DIR"
            (
                cd "$REP_DIR"
                admixture \
                    --cv="${CV_FOLDS}" \
                    -j"${THREADS_PER_RUN}" \
                    --seed="${rep}" \
                    "${PLINK_PREFIX}.bed" "${K}" \
                    > "run.log" 2>&1
            ) &
            pids+=($!)
        done

        # Wait for all replicates; allow individual failures
        for pid in "${pids[@]}"; do
            wait "$pid" || true
        done

        # Select best replicate by highest log-likelihood
        BEST_REP=""
        BEST_LL="-999999999999"
        for rep in $(seq -w 1 "$N_REP"); do
            REP_DIR="${K_DIR}/rep${rep}"
            LOG_R="${REP_DIR}/run.log"

            LL=$(grep "^Loglikelihood:" "$LOG_R" 2>/dev/null | tail -1 | awk '{print $2}') || true
            CV_ERR=$(grep "CV error" "$LOG_R" 2>/dev/null | tail -1 | awk -F': ' '{print $2}') || true
            [ -z "$LL" ]     && LL="-999999999999"
            [ -z "$CV_ERR" ] && CV_ERR="NA"

            printf "%d\t%s\t%s\t%s\n" "$K" "$rep" "$LL" "$CV_ERR" \
                >> "${OUT_T}/loglik_summary.txt"

            if awk "BEGIN {exit !($LL > $BEST_LL)}"; then
                BEST_LL="$LL"
                BEST_REP="$rep"
            fi
        done

        if [ -z "$BEST_REP" ]; then
            echo "  ERROR: no valid replicate for K=${K} — skipping" | tee -a "$LOG_FILE"
            continue
        fi

        echo "  K=${K}: best replicate = rep${BEST_REP} (LL=${BEST_LL})" | tee -a "$LOG_FILE"

        # Copy best Q and P files to threshold directory
        Q_BEST="${K_DIR}/rep${BEST_REP}/admixture_input.${K}.Q"
        P_BEST="${K_DIR}/rep${BEST_REP}/admixture_input.${K}.P"
        [ -s "$Q_BEST" ] && cp "$Q_BEST" "${OUT_T}/admixture_input.${K}.Q"
        [ -s "$P_BEST" ] && cp "$P_BEST" "${OUT_T}/admixture_input.${K}.P"

        # Append CV error from best replicate to cross_validation.txt
        CV_BEST=$(grep "CV error" "${K_DIR}/rep${BEST_REP}/run.log" 2>/dev/null | tail -1)
        [ -n "$CV_BEST" ] && echo "$CV_BEST" >> "${OUT_T}/cross_validation.txt"
    done

    # -------------------------------------------------------------------
    # STEP 4: Generate named Q files (individual IDs + Q matrix)
    # -------------------------------------------------------------------
    IDS_FILE="${OUT_T}/individual_ids.txt"
    awk '{print $2}' "${PLINK_PREFIX}.fam" > "$IDS_FILE"

    for K in $(seq "$K_MIN" "$K_MAX"); do
        Q_FILE="${OUT_T}/admixture_input.${K}.Q"
        if [ -s "$Q_FILE" ]; then
            paste "$IDS_FILE" "$Q_FILE" > "${OUT_T}/admixture_input.${K}.named.Q"
        else
            echo "  WARN: Q file missing for K=${K} — named.Q not generated" | tee -a "$LOG_FILE"
        fi
    done

    # -------------------------------------------------------------------
    # STEP 5: Build population map for R plotting (sample_id → site)
    # -------------------------------------------------------------------
    SAMPLE_MAP="${OUT_T}/sample_map.tsv"
    awk -F',' 'NR==FNR {keep[$1]=1; next}
               FNR>1 && $1 in keep && $4!="" {print $1"\t"$4}' \
        "$IDS_FILE" "$SAMPLE_SHEET" > "$SAMPLE_MAP"

    # -------------------------------------------------------------------
    # STEP 6: R plots (CV error curve + STRUCTURE bar plots)
    # -------------------------------------------------------------------
    echo "  R plotting — $TAG" | tee -a "$LOG_FILE"
    Rscript --vanilla "$R_SCRIPT" \
        "$OUT_T" \
        "$SAMPLE_MAP" \
        "$K_MIN" \
        "$K_MAX" \
        "admixture_input" >> "$LOG_FILE" 2>&1 \
        || echo "  WARN: R plotting failed for $TAG — check log" | tee -a "$LOG_FILE"

    echo "  DONE $TAG — plots: ${OUT_T}/plots/" | tee -a "$LOG_FILE"
done

echo "DONE $(date -Iseconds) 07b.2 ADMIXTURE completed" | tee -a "$LOG_FILE"
echo "Outputs:"
for THRESH in 70 80 90 95; do
    echo "  T${THRESH}: ${OUT_BASE}/T${THRESH}/plots/"
done

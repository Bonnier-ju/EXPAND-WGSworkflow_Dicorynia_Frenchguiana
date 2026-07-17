#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=180:00:00
#SBATCH -J 07b.4-ne
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   sbatch 07b.4-ne.slurm
#
# Current Ne (currentNe2) and temporal Ne (GONE2) estimation on
# purity-filtered individuals (T70/T80/T90/T95), per genetic cluster
# (Pop_1/Pop_2/Pop_3, K=3 from sites_by_clusters.csv).
#
# Tools (same approach as 07.10-Ne_global):
#   - currentNe2 (default) : current Ne, panmictic assumption, 50%/90% CI
#   - currentNe2 (-x)      : sub-structure within each cluster (FST, s)
#                            only if N >= MIN_N_STRUCT individuals
#   - GONE2 (default)      : temporal Ne trajectory, panmictic assumption
#   - GONE2 (-x)           : structure-corrected trajectory
#                            only if N >= MIN_N_STRUCT individuals
#
# Parameters:
#   -s 200000 : currentNe2 internal subsampling to 200k SNP pairs
#   -r 1      : uniform 1 cM/Mb recombination rate (no Dicorynia map)
#   -t 8      : threads (= SLURM cpus-per-task)
#   -k 0      : random-mating (currentNe2 default only; incompatible with -x)
#   MAF >= 0.05 computed within each cluster subset
#   PED pre-thinned to ~1.5M SNPs (< MAXLOCI=2M compile limit)
#
# Site exclusions: St_georges and Petit_Croissant (confirmed 2026-07-07)
#
# Output: 07b-population_structure_pur_inds/07b.4-ne/{TAG}/{POP}/
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/1.90b7
module load bioinfo/currentNe2/1.0.1
module load bioinfo/GONE2/1.0.2
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

CPUS="${SLURM_CPUS_PER_TASK:-8}"
N_CHR=13

# -------------------------------------------------------------------
# STEP 1: Paths and parameters
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
WD="${BASE}/07b-population_structure_pur_inds"

INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
SAMPLE_SHEET="${BASE}/metadata/sample_sheet_complete_filtered.csv"
SITES_BY_CLUSTERS="${BASE}/metadata/sites_by_clusters.csv"
LISTS_DIR="${WD}/07b.0-purity_filtering/individual_lists"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07b.4-ne.R"

# Reuse chr13 VCF from 07.10 (already filtered + chromosomes renamed 1-13)
VCF_CHR13="${BASE}/07-population_structure/workflow_full/07.10-Ne_global/tmp/vcf_chr13.vcf.gz"

OUT_BASE="${WD}/07b.4-ne"
LOG_FILE="${WD}/logs/07b.4-ne.log"
TMP_DIR="${OUT_BASE}/tmp"

# Minimum individuals to run Ne tools; -x mode requires more
MIN_N_NE=10
MIN_N_STRUCT=20

for f in "$INPUT_VCF" "$SAMPLE_SHEET" "$SITES_BY_CLUSTERS" "$VCF_CHR13" "$R_SCRIPT"; do
    [ -s "$f" ] || { echo "ERROR: missing required input: $f"; exit 1; }
done

mkdir -p "$OUT_BASE" "$TMP_DIR" "${WD}/logs"
echo "START $(date -Iseconds) 07b.4 Ne estimation per cluster (purity-filtered)" > "$LOG_FILE"
echo "Parameters: MIN_N_NE=${MIN_N_NE} | MIN_N_STRUCT=${MIN_N_STRUCT} | CPUs=${CPUS}" | tee -a "$LOG_FILE"
echo "Chr13 VCF (reused): $VCF_CHR13" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Build site-exclusion list (St_georges + Petit_Croissant)
# -------------------------------------------------------------------
EXCLUDE_FILE="${TMP_DIR}/exclude_sites.txt"
awk -F',' 'NR>1 && ($4 == "St_georges" || $4 == "Petit_Croissant") {print $1}' \
    "$SAMPLE_SHEET" > "$EXCLUDE_FILE"
N_EXCL=$(wc -l < "$EXCLUDE_FILE")
echo "Site exclusion: $N_EXCL individuals from St_georges + Petit_Croissant" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Build site → cluster (K=3) map from sites_by_clusters.csv
#   Col 1 = site name, col 3 = K=3 assignment
#   Exclude St_georges and Petit_Croissant
# -------------------------------------------------------------------
SITE_CLUSTER="${TMP_DIR}/site_cluster_k3.tsv"
awk -F',' 'NR>1 && $1 != "St_georges" && $1 != "Petit_Croissant" {print $1"\t"$3}' \
    "$SITES_BY_CLUSTERS" > "$SITE_CLUSTER"
echo "Site→cluster map (K=3): $(wc -l < "$SITE_CLUSTER") sites" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: Loop over purity thresholds
# -------------------------------------------------------------------
for THRESH in 70 80 90 95; do
    TAG="T${THRESH}"
    KEEP_RAW="${LISTS_DIR}/inds_${TAG}_all.txt"

    if [ ! -s "$KEEP_RAW" ]; then
        echo "WARN: keep file not found for $TAG — skipping" | tee -a "$LOG_FILE"
        continue
    fi

    OUT_T="${OUT_BASE}/${TAG}"
    mkdir -p "${OUT_T}"
    echo "=== $TAG ===" | tee -a "$LOG_FILE"

    # STEP 4a: Extract purity-filtered IIDs (col 2 of keep file = IID, FID=0)
    PURITY_IDS="${OUT_T}/purity_ids.txt"
    awk '{print $2}' "$KEEP_RAW" > "$PURITY_IDS"
    echo "  Purity-filtered: $(wc -l < "$PURITY_IDS") individuals" | tee -a "$LOG_FILE"

    # STEP 4b: Per-cluster keep files
    #   Intersect: purity-filtered × K=3 cluster × site-excluded
    #   Output format: IID<tab>IID (PLINK --double-id convention)
    for POP in Pop_1 Pop_2 Pop_3; do
        OUT_POP="${OUT_T}/${POP}"
        mkdir -p "${OUT_POP}/plink" \
                 "${OUT_POP}/currentNe2" \
                 "${OUT_POP}/GONE2/default" \
                 "${OUT_POP}/GONE2/structure"

        POP_KEEP="${OUT_POP}/keep.txt"

        # Build keep: purity IDs × K=3 cluster assignment × not excluded
        # Output: IID<tab>IID for PLINK --double-id
        awk -v pop="$POP" \
        'FILENAME==ARGV[1] { pur[$1]=1; next }
         FILENAME==ARGV[2] { clust[$1]=$2; next }
         FILENAME==ARGV[3] { excl[$1]=1; next }
         FNR==1 { next }
         {
             n=split($0, f, ",")
             id=f[1]; site=f[4]
             gsub(/\r/, "", id); gsub(/\r/, "", site)
             if (id in excl) next
             if (!(site in clust)) next
             if (clust[site] != pop) next
             if (!(id in pur)) next
             print id"\t"id
         }' \
         "$PURITY_IDS" "$SITE_CLUSTER" "$EXCLUDE_FILE" "$SAMPLE_SHEET" \
         > "$POP_KEEP"

        N_IND=$(wc -l < "$POP_KEEP")
        echo "  ${POP}: $N_IND individuals" | tee -a "$LOG_FILE"

        if [ "$N_IND" -lt "$MIN_N_NE" ]; then
            echo "  WARN: ${POP} has $N_IND individuals (< ${MIN_N_NE}) — skipping Ne tools" | tee -a "$LOG_FILE"
            continue
        fi

        LABEL="${POP}_${TAG}"
        PLINK_PREFIX="${OUT_POP}/plink/${LABEL}"
        NE2_PED="${OUT_POP}/plink/${LABEL}_ne2.ped"
        NE2_MAP="${OUT_POP}/plink/${LABEL}_ne2.map"

        # ----------------------------------------------------------
        # STEP 5: PLINK VCF → BED (MAF computed within cluster subset)
        # ----------------------------------------------------------
        if [ ! -s "${PLINK_PREFIX}.bed" ]; then
            echo "  STEP 5: PLINK VCF→BED — ${POP} ${TAG}" | tee -a "$LOG_FILE"
            plink \
                --vcf "$VCF_CHR13" \
                --double-id \
                --keep "$POP_KEEP" \
                --chr-set ${N_CHR} \
                --maf 0.05 \
                --geno 0.1 \
                --make-bed \
                --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1
        fi

        if [ ! -s "${PLINK_PREFIX}.bed" ]; then
            echo "  ERROR: PLINK failed for ${POP} ${TAG}" | tee -a "$LOG_FILE"
            continue
        fi

        N_SNP=$(wc -l < "${PLINK_PREFIX}.bim")
        echo "  ${PLINK_PREFIX}: $N_IND ind, $N_SNP SNPs (MAF>=0.05)" | tee -a "$LOG_FILE"

        # ----------------------------------------------------------
        # STEP 5b: Pre-thin to PED (target ~1.5M SNPs, < MAXLOCI=2M)
        # ----------------------------------------------------------
        if [ ! -s "$NE2_PED" ]; then
            THIN_RATIO=$(awk -v n="$N_SNP" 'BEGIN {t=1500000/n; print (t>0.95?0.95:t)}')
            echo "  Thinning: --thin ${THIN_RATIO} from ${N_SNP} SNPs (target 1.5M)" | tee -a "$LOG_FILE"
            plink \
                --bfile "$PLINK_PREFIX" \
                --chr-set ${N_CHR} \
                --thin "$THIN_RATIO" \
                --recode \
                --out "${OUT_POP}/plink/${LABEL}_ne2" >> "$LOG_FILE" 2>&1
        fi

        N_SNP_NE2=$(wc -l < "$NE2_MAP")
        echo "  Ne2 PED: $N_SNP_NE2 SNPs" | tee -a "$LOG_FILE"

        # ----------------------------------------------------------
        # STEP 6: currentNe2 default (panmictic)
        # ----------------------------------------------------------
        CNE2_DEFAULT="${OUT_POP}/currentNe2/${LABEL}_default.txt"
        if [ ! -s "$CNE2_DEFAULT" ]; then
            echo "  STEP 6: currentNe2 default — ${POP} ${TAG}" | tee -a "$LOG_FILE"
            (
                cd "${OUT_POP}/plink"
                currentne2 \
                    -s 200000 \
                    -r 1 \
                    -t "${CPUS}" \
                    -k 0 \
                    -o "$CNE2_DEFAULT" \
                    "${LABEL}_ne2.ped" >> "$LOG_FILE" 2>&1
            ) || echo "  WARN: currentNe2 default returned non-zero — ${POP} ${TAG}" | tee -a "$LOG_FILE"
            echo "  currentNe2 default done — ${POP} ${TAG}" | tee -a "$LOG_FILE"
        else
            echo "  currentNe2 default: reusing existing output" | tee -a "$LOG_FILE"
        fi

        # ----------------------------------------------------------
        # STEP 7: currentNe2 -x (sub-structure within cluster)
        #         Only if N >= MIN_N_STRUCT
        # ----------------------------------------------------------
        CNE2_MIX="${OUT_POP}/currentNe2/${LABEL}_mix.txt"
        if [ "$N_IND" -ge "$MIN_N_STRUCT" ]; then
            if [ ! -s "$CNE2_MIX" ]; then
                echo "  STEP 7: currentNe2 -x — ${POP} ${TAG}" | tee -a "$LOG_FILE"
                (
                    cd "${OUT_POP}/plink"
                    currentne2 \
                        -x \
                        -s 200000 \
                        -r 1 \
                        -t "${CPUS}" \
                        -o "$CNE2_MIX" \
                        "${LABEL}_ne2.ped" >> "$LOG_FILE" 2>&1
                ) || echo "  WARN: currentNe2 -x returned non-zero — ${POP} ${TAG}" | tee -a "$LOG_FILE"
                echo "  currentNe2 -x done — ${POP} ${TAG}" | tee -a "$LOG_FILE"
            else
                echo "  currentNe2 -x: reusing existing output" | tee -a "$LOG_FILE"
            fi
        else
            echo "  currentNe2 -x: skipped (N=${N_IND} < MIN_N_STRUCT=${MIN_N_STRUCT})" | tee -a "$LOG_FILE"
        fi

        # ----------------------------------------------------------
        # STEP 8: GONE2 default (panmictic temporal trajectory)
        #         Must run from output dir; PED/MAP symlinked as ${LABEL}.*
        # ----------------------------------------------------------
        GONE2_DEF_NE="${OUT_POP}/GONE2/default/${LABEL}.ped_GONE2_Ne"
        if [ ! -s "$GONE2_DEF_NE" ]; then
            echo "  STEP 8: GONE2 default — ${POP} ${TAG}" | tee -a "$LOG_FILE"
            ln -sf "$NE2_MAP" "${OUT_POP}/GONE2/default/${LABEL}.map"
            ln -sf "$NE2_PED" "${OUT_POP}/GONE2/default/${LABEL}.ped"
            (
                cd "${OUT_POP}/GONE2/default"
                gone2 -r 1 -t "${CPUS}" "${LABEL}.ped" >> "$LOG_FILE" 2>&1
            ) || echo "  WARN: GONE2 default returned non-zero — ${POP} ${TAG}" | tee -a "$LOG_FILE"
            echo "  GONE2 default done — ${POP} ${TAG}" | tee -a "$LOG_FILE"
        else
            echo "  GONE2 default: reusing existing output" | tee -a "$LOG_FILE"
        fi

        # ----------------------------------------------------------
        # STEP 9: GONE2 -x (structure-corrected trajectory)
        #         Only if N >= MIN_N_STRUCT
        # ----------------------------------------------------------
        GONE2_STR_NE="${OUT_POP}/GONE2/structure/${LABEL}.ped_GONE2_Ne_mix"
        if [ "$N_IND" -ge "$MIN_N_STRUCT" ]; then
            if [ ! -s "$GONE2_STR_NE" ]; then
                echo "  STEP 9: GONE2 -x — ${POP} ${TAG}" | tee -a "$LOG_FILE"
                ln -sf "$NE2_MAP" "${OUT_POP}/GONE2/structure/${LABEL}.map"
                ln -sf "$NE2_PED" "${OUT_POP}/GONE2/structure/${LABEL}.ped"
                (
                    cd "${OUT_POP}/GONE2/structure"
                    gone2 -x -r 1 -t "${CPUS}" "${LABEL}.ped" >> "$LOG_FILE" 2>&1
                ) || echo "  WARN: GONE2 -x returned non-zero — ${POP} ${TAG}" | tee -a "$LOG_FILE"
                echo "  GONE2 -x done — ${POP} ${TAG}" | tee -a "$LOG_FILE"
            else
                echo "  GONE2 -x: reusing existing output" | tee -a "$LOG_FILE"
            fi
        else
            echo "  GONE2 -x: skipped (N=${N_IND} < MIN_N_STRUCT=${MIN_N_STRUCT})" | tee -a "$LOG_FILE"
        fi

        echo "  DONE ${POP} ${TAG} $(date -Iseconds)" | tee -a "$LOG_FILE"
    done

    # ----------------------------------------------------------
    # STEP 10: R compilation per threshold
    # ----------------------------------------------------------
    echo "  R compilation — $TAG" | tee -a "$LOG_FILE"
    Rscript --vanilla "$R_SCRIPT" \
        "$OUT_BASE" \
        "$TAG" >> "$LOG_FILE" 2>&1 \
        || echo "  WARN: R script failed for $TAG — check log" | tee -a "$LOG_FILE"

    echo "  DONE $TAG — plots: ${OUT_T}/plots/" | tee -a "$LOG_FILE"
done

# -------------------------------------------------------------------
# STEP 11: R cross-threshold summary (all T combined)
# -------------------------------------------------------------------
echo "STEP 11: R cross-threshold summary" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$OUT_BASE" \
    "ALL" >> "$LOG_FILE" 2>&1 \
    || echo "WARN: R cross-threshold summary failed — check log" | tee -a "$LOG_FILE"

echo "DONE $(date -Iseconds) 07b.4 Ne estimation completed" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_BASE}/"

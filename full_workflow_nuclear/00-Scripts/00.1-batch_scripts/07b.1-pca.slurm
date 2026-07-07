#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J 07b.1-pca
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment
# -------------------------------------------------------------------
module purge
module load bioinfo/PLINK/1.90b7
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

# -------------------------------------------------------------------
# STEP 1: Paths
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
WD="${BASE}/07b-population_structure_pur_inds"

VCF_LD="${WD}/input/snps_TE_filtered.ld_pruned.vcf.gz"
Q_FILE="${WD}/input/admixture_input.3.named.Q"
SAMPLE_SHEET="${BASE}/metadata/sample_sheet_complete_filtered.csv"
LISTS_DIR="${WD}/07b.0-purity_filtering/individual_lists"
OUT_BASE="${WD}/07b.1-pca"
LOG_FILE="${WD}/logs/07b.1-pca.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07b.1-pca_plot.R"

# Q-column to cluster mapping confirmed by 07b.0:
#   Q1 = Pop_1 | Q2 = Pop_2 | Q3 = Pop_3
# Used in awk cluster map generation below.

for f in "$VCF_LD" "$Q_FILE" "$SAMPLE_SHEET" "$R_SCRIPT"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: required input missing: $f"
        exit 1
    fi
done

echo "START $(date -Iseconds) 07b.1 PCA on purity-filtered individuals" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Loop over purity thresholds
# -------------------------------------------------------------------
for THRESH in 70 80 90 95; do
    TAG="T${THRESH}"
    KEEP="${LISTS_DIR}/inds_${TAG}_all.txt"

    if [ ! -s "$KEEP" ]; then
        echo "WARN: keep file not found: $KEEP — skipping $TAG" | tee -a "$LOG_FILE"
        continue
    fi

    N_IND=$(wc -l < "$KEEP")
    N_PC=$(( N_IND < 21 ? N_IND - 1 : 20 ))
    OUT_T="${OUT_BASE}/${TAG}"
    mkdir -p "$OUT_T"

    echo "--- $TAG: $N_IND individuals, $N_PC PCs ---" | tee -a "$LOG_FILE"

    # STEP 2a: Extract individual IDs (col 2 of keep file: IND_ID)
    IDS_FILE="${OUT_T}/keep_ids.txt"
    cut -f2 "$KEEP" > "$IDS_FILE"

    # STEP 2b: Convert keep file to FID=IID format required by PLINK --double-id
    KEEP_PLINK="${OUT_T}/keep_plink_${TAG}.txt"
    awk '{print $2"\t"$2}' "$KEEP" > "$KEEP_PLINK"

    # STEP 2c: PLINK — filter VCF and convert to BED
    PLINK_PREFIX="${OUT_T}/plink_input_${TAG}"
    plink \
        --vcf "$VCF_LD" \
        --allow-extra-chr \
        --keep "$KEEP_PLINK" \
        --double-id \
        --set-missing-var-ids @:#\$1,\$2 \
        --make-bed \
        --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

    if [ ! -s "${PLINK_PREFIX}.bed" ]; then
        echo "ERROR: PLINK BED conversion failed for $TAG" | tee -a "$LOG_FILE"
        continue
    fi

    # STEP 2d: PLINK — PCA
    PCA_PREFIX="${OUT_T}/pca_${TAG}"
    plink \
        --bfile "$PLINK_PREFIX" \
        --allow-extra-chr \
        --pca ${N_PC} header tabs \
        --out "$PCA_PREFIX" >> "$LOG_FILE" 2>&1

    echo "  PCA done → ${PCA_PREFIX}.eigenvec" | tee -a "$LOG_FILE"

    # STEP 2e: Build cluster map (Pop_1/2/3 from Q file)
    # Q column mapping: Q2=Pop_1, Q3=Pop_2, Q4=Pop_3 (Q file: $1=ID $2=Q1 $3=Q2 $4=Q3)
    CLUSTER_MAP="${OUT_T}/cluster_map_${TAG}.tsv"
    awk 'NR==FNR {keep[$1]=1; next}
         $1 in keep {
             if ($2 >= $3 && $2 >= $4) cl="Pop_1"
             else if ($3 >= $2 && $3 >= $4) cl="Pop_2"
             else cl="Pop_3"
             print $1"\t"cl
         }' "$IDS_FILE" "$Q_FILE" > "$CLUSTER_MAP"

    # STEP 2f: Build site map (sample_id → site, from sample sheet)
    SITE_MAP="${OUT_T}/site_map_${TAG}.tsv"
    awk -F',' 'NR==FNR {keep[$1]=1; next}
               FNR>1 && $1 in keep && $4!="" {print $1"\t"$4}' \
        "$IDS_FILE" "$SAMPLE_SHEET" > "$SITE_MAP"

    # STEP 2g: R plots (cluster + site coloring in one call)
    echo "  R plots — $TAG" | tee -a "$LOG_FILE"
    Rscript --vanilla "$R_SCRIPT" \
        "$OUT_T" \
        "$CLUSTER_MAP" \
        "$SITE_MAP" \
        "${PCA_PREFIX}.eigenvec" \
        "${PCA_PREFIX}.eigenval" \
        "$TAG" >> "$LOG_FILE" 2>&1

    echo "  Plots written to ${OUT_T}/plots/" | tee -a "$LOG_FILE"
done

echo "DONE $(date -Iseconds) 07b.1 PCA completed" | tee -a "$LOG_FILE"
echo "Outputs: $OUT_BASE"

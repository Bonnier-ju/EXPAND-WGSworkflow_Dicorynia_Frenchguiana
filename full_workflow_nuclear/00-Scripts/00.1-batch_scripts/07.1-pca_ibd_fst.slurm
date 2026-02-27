#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J 07.1-pca_ibd_fst
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/1.90b7
# Optional for automatic PCA plotting (Genotoul FAQ: prefer statistics/R/4.3.0)
module load statistics/R/4.3.0 >/dev/null 2>&1 || module load devel/R >/dev/null 2>&1 || true

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
VCF_66="${BASE}/06-variant_filtering/workflow_test/06.6-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
VCF_65="${BASE}/06-variant_filtering/workflow_test/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"

SAMPLE_SHEET="${SAMPLE_SHEET:-${BASE}/metadata/sample_sheet_test.csv}"
POP_COLUMN="${POP_COLUMN:-site}"   # ex: site or project

OUT_DIR="${BASE}/07-population_structure/workflow_test/07.1-pca_ibd_fst"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"
R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.1-pca_plot.R"

LOG_FILE="${LOG_DIR}/07.1-pca_ibd_fst.log"
SUMMARY_TXT="${OUT_DIR}/07.1-summary.txt"

if [ -s "$VCF_66" ]; then
    INPUT_VCF="$VCF_66"
    INPUT_STEP="06.6"
else
    INPUT_VCF="$VCF_65"
    INPUT_STEP="06.5"
fi

echo "START $(date -Iseconds) 07.1 PCA/IBD/FST" > "$LOG_FILE"
echo "Input VCF: $INPUT_VCF (from $INPUT_STEP)" | tee -a "$LOG_FILE"
echo "Sample sheet: $SAMPLE_SHEET (POP_COLUMN=$POP_COLUMN)" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found: $SAMPLE_SHEET" | tee -a "$LOG_FILE"
    exit 1
fi

VCF_SAMPLES="${TMP_DIR}/vcf_samples.txt"
bcftools query -l "$INPUT_VCF" | sort -u > "$VCF_SAMPLES"

MAP_RAW="${TMP_DIR}/sample_population_raw.tsv"
MAP_FINAL="${OUT_DIR}/sample_population_map.tsv"
KEEP_FILE="${TMP_DIR}/keep.samples.txt"
CLUSTER_FILE="${TMP_DIR}/plink.cluster.txt"
POPS_LIST="${TMP_DIR}/populations.list"

awk -F',' -v popcol="$POP_COLUMN" '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="sample_id") sid=i
        if ($i==popcol) pop=i
    }
    if (!sid || !pop) {
        print "ERROR: sample_id or " popcol " column not found in sample sheet" > "/dev/stderr"
        exit 2
    }
    next
}
{
    gsub(/\r/, "", $sid)
    gsub(/\r/, "", $pop)
    if ($sid != "" && $pop != "" && $pop != "NA" && $pop != "na") {
        print $sid "\t" $pop
    }
}
' "$SAMPLE_SHEET" > "$MAP_RAW"

join -t $'\t' -1 1 -2 1 \
    <(sort -k1,1 "$MAP_RAW") \
    <(awk '{print $1"\t"$1}' "$VCF_SAMPLES" | sort -k1,1) \
    | awk -F'\t' '{print $1"\t"$2}' > "$MAP_FINAL"

if [ ! -s "$MAP_FINAL" ]; then
    echo "ERROR: no sample overlap between VCF and metadata for POP_COLUMN=$POP_COLUMN" | tee -a "$LOG_FILE"
    exit 1
fi

awk -F'\t' '{print $1"\t"$1}' "$MAP_FINAL" > "$KEEP_FILE"
awk -F'\t' '{print $1"\t"$1"\t"$2}' "$MAP_FINAL" > "$CLUSTER_FILE"
awk -F'\t' '{print $2}' "$MAP_FINAL" | sort -u > "$POPS_LIST"

N_SAMPLES=$(wc -l < "$MAP_FINAL")
N_POPS=$(wc -l < "$POPS_LIST")
echo "Samples kept: $N_SAMPLES" | tee -a "$LOG_FILE"
echo "Populations detected: $N_POPS" | tee -a "$LOG_FILE"

PLINK_PREFIX="${OUT_DIR}/genotypes"
echo "STEP: build PLINK bed/bim/fam" | tee -a "$LOG_FILE"
plink \
    --vcf "$INPUT_VCF" \
    --allow-extra-chr \
    --keep "$KEEP_FILE" \
    --double-id \
    --set-missing-var-ids @:#\$1,\$2 \
    --make-bed \
    --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "${PLINK_PREFIX}.bed" ] || [ ! -s "${PLINK_PREFIX}.bim" ] || [ ! -s "${PLINK_PREFIX}.fam" ]; then
    echo "ERROR: PLINK conversion failed" | tee -a "$LOG_FILE"
    exit 1
fi

PCA_PREFIX="${OUT_DIR}/pca"
echo "STEP: PCA" | tee -a "$LOG_FILE"
plink \
    --bfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --pca 20 header tabs \
    --out "$PCA_PREFIX" >> "$LOG_FILE" 2>&1

IBD_PREFIX="${OUT_DIR}/ibd"
echo "STEP: IBD (pairwise genome)" | tee -a "$LOG_FILE"
plink \
    --bfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --genome full \
    --out "$IBD_PREFIX" >> "$LOG_FILE" 2>&1

FST_PREFIX="${OUT_DIR}/fst"
FST_STATUS="SKIPPED"
if [ "$N_POPS" -ge 2 ]; then
    echo "STEP: FST" | tee -a "$LOG_FILE"
    plink \
        --bfile "$PLINK_PREFIX" \
        --allow-extra-chr \
        --within "$CLUSTER_FILE" \
        --fst \
        --out "$FST_PREFIX" >> "$LOG_FILE" 2>&1
    FST_STATUS="DONE"
else
    echo "WARN: FST skipped (<2 populations)" | tee -a "$LOG_FILE"
fi

# Optional: automatic PCA plotting with R
if [ "${RUN_R_PLOT:-1}" = "1" ]; then
    if [ -s "$R_PLOT_SCRIPT" ] && command -v Rscript >/dev/null 2>&1; then
        echo "STEP: PCA plotting in R" | tee -a "$LOG_FILE"
        Rscript --vanilla "$R_PLOT_SCRIPT" "$OUT_DIR" >> "$LOG_FILE" 2>&1 || \
            echo "WARN: R plotting step failed (analysis outputs remain valid)" | tee -a "$LOG_FILE"
    else
        echo "WARN: R plotting skipped (Rscript or script file missing)" | tee -a "$LOG_FILE"
    fi
fi

{
    echo "07.1 PCA / IBD / FST summary"
    echo "Date: $(date -Iseconds)"
    echo "Input step: $INPUT_STEP"
    echo "Input VCF: $INPUT_VCF"
    echo "Samples used: $N_SAMPLES"
    echo "Populations used: $N_POPS"
    echo "Population map: $MAP_FINAL"
    echo
    echo "Outputs:"
    echo "  PCA: ${PCA_PREFIX}.eigenvec ; ${PCA_PREFIX}.eigenval"
    echo "  IBD: ${IBD_PREFIX}.genome"
    echo "  FST: ${FST_PREFIX}.* (${FST_STATUS})"
    echo "  PLINK: ${PLINK_PREFIX}.bed/.bim/.fam"
    echo "  PCA plots (R): ${OUT_DIR}/plots/* (if RUN_R_PLOT=1 and R available)"
} > "$SUMMARY_TXT"

echo "DONE $(date -Iseconds) 07.1 completed" | tee -a "$LOG_FILE"
cat "$SUMMARY_TXT"

#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J 07.1-pca
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
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
VCF_66="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
VCF_65="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"

SAMPLE_SHEET="${SAMPLE_SHEET:-${BASE}/metadata/sample_sheet_complete_filtered.csv}"
POP_COLUMN="${POP_COLUMN:-site}"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.1-pca_ibd_fst"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"
R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.1-pca_plot.R"

LOG_FILE="${LOG_DIR}/07.1-pca_ibd_fst.log"
SUMMARY_TXT="${OUT_DIR}/07.1-summary.txt"

# Individuals excluded from specific PCA groups.
# These lists are used to build per-group keep files.
# Group 2 (no_outgroup): exclude outgroup individuals (Cameroun/Benin + Venezuela)
OUTGROUP_IDS="DBR_AGN DBR_AGO DBR_AGP DBR_AGQ DBR_AGR"
# Group 3 (guiana_only): additionally exclude Treemutation project individuals
TREEMUTATION_IDS="DBR_AGJ DBR_AGK DBR_AGL DBR_AGM"

if [ -s "$VCF_66" ]; then
    INPUT_VCF="$VCF_66"
    INPUT_STEP="06.6b"
else
    INPUT_VCF="$VCF_65"
    INPUT_STEP="06.5"
fi

echo "START $(date -Iseconds) 07.1 PCA (3 independent analyses)" > "$LOG_FILE"
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

# -------------------------------------------------------------------
# Build full sample → population map
# -------------------------------------------------------------------
VCF_SAMPLES="${TMP_DIR}/vcf_samples.txt"
bcftools query -l "$INPUT_VCF" | sort -u > "$VCF_SAMPLES"

MAP_RAW="${TMP_DIR}/sample_population_raw.tsv"
MAP_FINAL="${OUT_DIR}/sample_population_map.tsv"

awk -F',' -v popcol="$POP_COLUMN" '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="sample_id") sid=i
        if ($i==popcol) pop=i
    }
    if (!sid || !pop) {
        print "ERROR: sample_id or " popcol " column not found" > "/dev/stderr"
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
    echo "ERROR: no sample overlap between VCF and metadata" | tee -a "$LOG_FILE"
    exit 1
fi

N_TOTAL=$(wc -l < "$MAP_FINAL")
echo "Total samples in map: $N_TOTAL" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# Build keep files for each PCA group
# -------------------------------------------------------------------
# Group 1: all samples
KEEP_ALL="${TMP_DIR}/keep_all.txt"
awk -F'\t' '{print $1"\t"$1}' "$MAP_FINAL" > "$KEEP_ALL"

# Group 2: exclude outgroups
KEEP_NO_OG="${TMP_DIR}/keep_no_outgroup.txt"
EXCLUDE_NO_OG=$(echo "$OUTGROUP_IDS" | tr ' ' '\n' | awk '{print "^"$1"$"}' | paste -sd'|')
awk -F'\t' -v pat="$EXCLUDE_NO_OG" '$1 !~ pat {print $1"\t"$1}' "$MAP_FINAL" > "$KEEP_NO_OG"

# Group 3: exclude outgroups + Treemutation
KEEP_GUIANA="${TMP_DIR}/keep_guiana_only.txt"
ALL_EXCLUDE="$OUTGROUP_IDS $TREEMUTATION_IDS"
EXCLUDE_GUIANA=$(echo "$ALL_EXCLUDE" | tr ' ' '\n' | awk '{print "^"$1"$"}' | paste -sd'|')
awk -F'\t' -v pat="$EXCLUDE_GUIANA" '$1 !~ pat {print $1"\t"$1}' "$MAP_FINAL" > "$KEEP_GUIANA"

N_ALL=$(wc -l < "$KEEP_ALL")
N_NO_OG=$(wc -l < "$KEEP_NO_OG")
N_GUIANA=$(wc -l < "$KEEP_GUIANA")
echo "Group all_samples  : ${N_ALL} individuals" | tee -a "$LOG_FILE"
echo "Group no_outgroup  : ${N_NO_OG} individuals" | tee -a "$LOG_FILE"
echo "Group guiana_only  : ${N_GUIANA} individuals" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# Helper function: run one PCA
# run_pca <group_label> <keep_file> <n_individuals>
# -------------------------------------------------------------------
run_pca() {
    local group="$1"
    local keep="$2"
    local n_ind="$3"
    local n_pc
    n_pc=$(( n_ind < 20 ? n_ind - 1 : 20 ))

    local plink_prefix="${OUT_DIR}/genotypes_${group}"
    local pca_prefix="${OUT_DIR}/pca_${group}"

    echo "STEP: PLINK conversion — ${group} (n=${n_ind})" | tee -a "$LOG_FILE"
    plink \
        --vcf "$INPUT_VCF" \
        --allow-extra-chr \
        --keep "$keep" \
        --double-id \
        --set-missing-var-ids @:#\$1,\$2 \
        --make-bed \
        --out "$plink_prefix" >> "$LOG_FILE" 2>&1

    if [ ! -s "${plink_prefix}.bed" ]; then
        echo "ERROR: PLINK conversion failed for group ${group}" | tee -a "$LOG_FILE"
        return 1
    fi

    echo "STEP: PCA ${n_pc} axes — ${group}" | tee -a "$LOG_FILE"
    plink \
        --bfile "$plink_prefix" \
        --allow-extra-chr \
        --pca ${n_pc} header tabs \
        --out "$pca_prefix" >> "$LOG_FILE" 2>&1

    echo "INFO: PCA ${group} done → ${pca_prefix}.eigenvec" | tee -a "$LOG_FILE"
}

# -------------------------------------------------------------------
# Run the three independent PCAs
# -------------------------------------------------------------------
run_pca "all_samples"  "$KEEP_ALL"    "$N_ALL"
run_pca "no_outgroup"  "$KEEP_NO_OG"  "$N_NO_OG"
run_pca "guiana_only"  "$KEEP_GUIANA" "$N_GUIANA"

# -------------------------------------------------------------------
# R plotting — one call per group, each with its own eigenvec/eigenval
# -------------------------------------------------------------------
if [ "${RUN_R_PLOT:-1}" = "1" ]; then
    if [ -s "$R_PLOT_SCRIPT" ] && command -v Rscript >/dev/null 2>&1; then
        for group in all_samples no_outgroup guiana_only; do
            echo "STEP: PCA plotting in R — ${group}" | tee -a "$LOG_FILE"
            Rscript --vanilla "$R_PLOT_SCRIPT" \
                "$OUT_DIR" \
                "$MAP_FINAL" \
                "${OUT_DIR}/pca_${group}.eigenvec" \
                "${OUT_DIR}/pca_${group}.eigenval" \
                "$group" >> "$LOG_FILE" 2>&1 || \
                echo "WARN: R plotting failed for ${group}" | tee -a "$LOG_FILE"
        done
    else
        echo "WARN: R plotting skipped (Rscript or script file missing)" | tee -a "$LOG_FILE"
    fi
fi

{
    echo "07.1 PCA summary — 3 independent analyses"
    echo "Date: $(date -Iseconds)"
    echo "Input step: $INPUT_STEP"
    echo "Input VCF: $INPUT_VCF"
    echo
    echo "Groups:"
    echo "  all_samples : ${N_ALL} individuals"
    echo "  no_outgroup : ${N_NO_OG} individuals (excl. ${OUTGROUP_IDS})"
    echo "  guiana_only : ${N_GUIANA} individuals (excl. outgroups + ${TREEMUTATION_IDS})"
    echo
    echo "PCA outputs:"
    for group in all_samples no_outgroup guiana_only; do
        echo "  ${group}: ${OUT_DIR}/pca_${group}.eigenvec"
    done
    echo "Plots: ${OUT_DIR}/plots/<group>/"
} > "$SUMMARY_TXT"

echo "DONE $(date -Iseconds) 07.1 completed" | tee -a "$LOG_FILE"
cat "$SUMMARY_TXT"

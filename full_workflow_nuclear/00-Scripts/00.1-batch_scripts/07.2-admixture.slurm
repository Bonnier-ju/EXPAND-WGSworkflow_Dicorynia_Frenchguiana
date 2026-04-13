#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=120:00:00
#SBATCH -J 07.2-admixture
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

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
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"

VCF_66="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
VCF_65="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"

if [ -s "$VCF_66" ]; then
    INPUT_VCF="$VCF_66"
    INPUT_STEP="06.6b"
else
    INPUT_VCF="$VCF_65"
    INPUT_STEP="06.5"
fi

SAMPLE_SHEET="${SAMPLE_SHEET:-${BASE}/metadata/sample_sheet_complete_filtered.csv}"
POP_COLUMN="${POP_COLUMN:-site}"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.2-admixture"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/07.2-admixture.log"

R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.2-admixture_plot.R"

# K range and cross-validation folds
K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-10}"
CV_FOLDS="${CV_FOLDS:-10}"
ADMIXTURE_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Individuals excluded from specific groups (same definitions as 07.1-pca)
# Group no_outgroup: exclude Cameroun/Benin + Venezuela outgroups
OUTGROUP_IDS="DBR_AGN DBR_AGO DBR_AGP DBR_AGQ DBR_AGR"
# Group guiana_only: additionally exclude Treemutation project individuals
TREEMUTATION_IDS="DBR_AGJ DBR_AGK DBR_AGL DBR_AGM"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 07.2 ADMIXTURE (3 groups)" > "$LOG_FILE"
echo "Input VCF: $INPUT_VCF (from $INPUT_STEP)" | tee -a "$LOG_FILE"
echo "K range: ${K_MIN}..${K_MAX} | CV folds: ${CV_FOLDS} | threads: ${ADMIXTURE_THREADS}" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found: $SAMPLE_SHEET" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build full sample -> population map
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
    gsub(/\r/, "", $sid); gsub(/\r/, "", $pop)
    if ($sid != "" && $pop != "" && $pop != "NA" && $pop != "na")
        print $sid "\t" $pop
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
# STEP 3: Build keep files for each group
# -------------------------------------------------------------------
# Keep files use FID=0 (--const-fid convention in PLINK)
KEEP_ALL="${TMP_DIR}/keep_all.txt"
awk -F'\t' '{print "0\t"$1}' "$MAP_FINAL" > "$KEEP_ALL"

KEEP_NO_OG="${TMP_DIR}/keep_no_outgroup.txt"
EXCLUDE_NO_OG=$(echo "$OUTGROUP_IDS" | tr ' ' '\n' | awk '{print "^"$1"$"}' | paste -sd'|')
awk -F'\t' -v pat="$EXCLUDE_NO_OG" '$1 !~ pat {print "0\t"$1}' "$MAP_FINAL" > "$KEEP_NO_OG"

KEEP_GUIANA="${TMP_DIR}/keep_guiana_only.txt"
ALL_EXCLUDE="$OUTGROUP_IDS $TREEMUTATION_IDS"
EXCLUDE_GUIANA=$(echo "$ALL_EXCLUDE" | tr ' ' '\n' | awk '{print "^"$1"$"}' | paste -sd'|')
awk -F'\t' -v pat="$EXCLUDE_GUIANA" '$1 !~ pat {print "0\t"$1}' "$MAP_FINAL" > "$KEEP_GUIANA"

N_ALL=$(wc -l < "$KEEP_ALL")
N_NO_OG=$(wc -l < "$KEEP_NO_OG")
N_GUIANA=$(wc -l < "$KEEP_GUIANA")
echo "Group all_samples : ${N_ALL} individuals" | tee -a "$LOG_FILE"
echo "Group no_outgroup : ${N_NO_OG} individuals" | tee -a "$LOG_FILE"
echo "Group guiana_only : ${N_GUIANA} individuals" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# Helper function: run one complete ADMIXTURE analysis
# run_admixture_group <group_label> <keep_file>
# -------------------------------------------------------------------
run_admixture_group() {
    local group="$1"
    local keep_file="$2"

    local group_dir="${OUT_DIR}/${group}"
    local plink_prefix="${group_dir}/admixture_input"
    mkdir -p "$group_dir"

    echo "--- GROUP: ${group} ---" | tee -a "$LOG_FILE"

    # PLINK VCF -> BED with sample subsetting
    echo "STEP: PLINK conversion — ${group}" | tee -a "$LOG_FILE"
    plink \
        --vcf "$INPUT_VCF" \
        --make-bed \
        --allow-extra-chr \
        --const-fid \
        --keep "$keep_file" \
        --out "$plink_prefix" >> "$LOG_FILE" 2>&1

    if [ ! -s "${plink_prefix}.bed" ]; then
        echo "ERROR: PLINK conversion failed for group ${group}" | tee -a "$LOG_FILE"
        return 1
    fi

    local n_ind n_snp
    n_ind=$(wc -l < "${plink_prefix}.fam")
    n_snp=$(wc -l < "${plink_prefix}.bim")
    echo "INFO: ${n_ind} individuals, ${n_snp} SNPs" | tee -a "$LOG_FILE"

    # Rename chromosomes to 0 (ADMIXTURE requires numeric chromosome IDs)
    awk '{$1=0; print $0}' "${plink_prefix}.bim" > "${plink_prefix}_tmp.bim"
    mv "${plink_prefix}_tmp.bim" "${plink_prefix}.bim"

    # Run ADMIXTURE for K = K_MIN..K_MAX
    # ADMIXTURE writes Q/P files to the current working directory
    echo "STEP: ADMIXTURE K=${K_MIN}..${K_MAX} — ${group}" | tee -a "$LOG_FILE"
    cd "$group_dir"
    for K in $(seq "$K_MIN" "$K_MAX"); do
        echo "  K=${K}..." | tee -a "$LOG_FILE"
        admixture \
            --cv="$CV_FOLDS" \
            -j"$ADMIXTURE_THREADS" \
            "${plink_prefix}.bed" "$K" \
            | tee "${group_dir}/log_K${K}.out" >> "$LOG_FILE"
    done

    # Collect cross-validation errors
    grep "CV error" "${group_dir}"/log_K*.out \
        > "${group_dir}/cross_validation.txt" || true
    echo "CV errors (${group}):" | tee -a "$LOG_FILE"
    cat "${group_dir}/cross_validation.txt" | tee -a "$LOG_FILE"

    # Associate individual names with Q matrices
    local ids_file="${group_dir}/individual_ids.txt"
    awk '{print $2}' "${plink_prefix}.fam" > "$ids_file"

    for K in $(seq "$K_MIN" "$K_MAX"); do
        local q_file="${group_dir}/admixture_input.${K}.Q"
        local named_q="${group_dir}/admixture_input.${K}.named.Q"
        if [ -s "$q_file" ]; then
            paste "$ids_file" "$q_file" > "$named_q"
        else
            echo "  WARN: Q file missing for K=${K}" | tee -a "$LOG_FILE"
        fi
    done

    # R plotting
    if [ -s "$R_PLOT_SCRIPT" ]; then
        echo "STEP: R plotting — ${group}" | tee -a "$LOG_FILE"
        Rscript --vanilla "$R_PLOT_SCRIPT" \
            "$group_dir" \
            "$MAP_FINAL" \
            "$K_MIN" \
            "$K_MAX" \
            "admixture_input" >> "$LOG_FILE" 2>&1 \
            || echo "WARN: R plotting failed for group ${group}" | tee -a "$LOG_FILE"
    fi

    echo "DONE group ${group}" | tee -a "$LOG_FILE"
}

# -------------------------------------------------------------------
# STEP 4: Run the three independent ADMIXTURE analyses
# -------------------------------------------------------------------
run_admixture_group "all_samples"  "$KEEP_ALL"
run_admixture_group "no_outgroup"  "$KEEP_NO_OG"
run_admixture_group "guiana_only"  "$KEEP_GUIANA"

echo "DONE $(date -Iseconds) 07.2 ADMIXTURE completed — 3 groups" | tee -a "$LOG_FILE"
echo "Outputs:"
for group in all_samples no_outgroup guiana_only; do
    echo "  ${group}: ${OUT_DIR}/${group}/plots/"
done

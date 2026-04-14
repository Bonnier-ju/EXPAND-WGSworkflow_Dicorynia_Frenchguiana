#!/bin/bash
#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J 07.4-ibd
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/VCFtools/0.1.16
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Inputs / outputs
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"

# VCF: LD-pruned SNPs (independent markers for FST estimation)
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"

SAMPLE_SHEET="${SAMPLE_SHEET:-${BASE}/metadata/sample_sheet_complete_filtered.csv}"
GEOLOC_FILE="${GEOLOC_FILE:-${BASE}/metadata/geoloc_site.csv}"
POP_COLUMN="${POP_COLUMN:-site}"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.4-ibd"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
FST_DIR="${OUT_DIR}/fst_vcftools"
LOG_FILE="${LOG_DIR}/07.4-ibd.log"

R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.4-ibd_mantel.R"

# Mantel test permutations
N_PERM="${N_PERM:-999}"

# Outgroup IDs to exclude (same as 07.1 / 07.2 — Cameroun/Benin + Venezuela)
OUTGROUP_IDS="DBR_AGN DBR_AGO DBR_AGP DBR_AGQ DBR_AGR"
# Treemutation project IDs (excluded in guiana_only group)
TREEMUTATION_IDS="DBR_AGJ DBR_AGK DBR_AGL DBR_AGM"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR" "$FST_DIR"
echo "START $(date -Iseconds) 07.4 IBD analysis" > "$LOG_FILE"
echo "Input VCF: $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Geoloc file: $GEOLOC_FILE" | tee -a "$LOG_FILE"
echo "Mantel permutations: $N_PERM" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"; exit 1
fi
if [ ! -s "$GEOLOC_FILE" ]; then
    echo "ERROR: geoloc file not found: $GEOLOC_FILE" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Patch VCFv4.3 → v4.2 if needed (vcftools compatibility)
# -------------------------------------------------------------------
VCF_VERSION=$(bcftools view -h "$INPUT_VCF" | grep "^##fileformat=" | sed 's/##fileformat=//')
VCFTOOLS_VCF="$INPUT_VCF"
if [ "$VCF_VERSION" = "VCFv4.3" ]; then
    VCFTOOLS_VCF="${TMP_DIR}/input_v42.vcf.gz"
    if [ ! -s "$VCFTOOLS_VCF" ]; then
        echo "INFO: patching VCFv4.3 -> v4.2 for vcftools" | tee -a "$LOG_FILE"
        zcat "$INPUT_VCF" | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
            | bgzip -c > "$VCFTOOLS_VCF"
        bcftools index -t "$VCFTOOLS_VCF"
    fi
fi

# -------------------------------------------------------------------
# STEP 3: Build population map and keep files for each group
# -------------------------------------------------------------------
VCF_SAMPLES="${TMP_DIR}/vcf_samples.txt"
bcftools query -l "$VCFTOOLS_VCF" | sort -u > "$VCF_SAMPLES"

MAP_FILE="${OUT_DIR}/sample_population_map.tsv"
awk -F',' -v popcol="$POP_COLUMN" '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="sample_id") sid=i
        if ($i==popcol) pop=i
    }
    next
}
{
    gsub(/\r/, "", $sid); gsub(/\r/, "", $pop)
    if ($sid != "" && $pop != "" && $pop != "NA" && $pop != "na")
        print $sid "\t" $pop
}
' "$SAMPLE_SHEET" \
| sort -k1,1 \
| join -t $'\t' -1 1 -2 1 \
    - \
    <(awk '{print $1"\t"$1}' "$VCF_SAMPLES" | sort -k1,1) \
| awk -F'\t' '{print $1"\t"$2}' > "$MAP_FILE"

# Build keep file: guiana_only only (outgroups + treemutation excluded)
KEEP_GUIANA="${TMP_DIR}/keep_guiana_only.txt"
ALL_EXCL="$OUTGROUP_IDS $TREEMUTATION_IDS"
EXCLUDE_G=$(echo "$ALL_EXCL" | tr ' ' '\n' | awk '{print "^"$1"$"}' | paste -sd'|')
awk -v pat="$EXCLUDE_G" '$1 !~ pat {print $1}' "$MAP_FILE" > "$KEEP_GUIANA"

echo "INFO: guiana_only=$(wc -l < "$KEEP_GUIANA") individuals" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: Pairwise FST per population pair (vcftools --weir-fst-pop)
# -------------------------------------------------------------------
run_fst_group() {
    local group="$1"
    local keep_file="$2"

    local group_dir="${FST_DIR}/${group}"
    mkdir -p "$group_dir"
    echo "--- FST group: ${group} ---" | tee -a "$LOG_FILE"

    # Build population lists restricted to the keep file
    local pop_map="${TMP_DIR}/map_${group}.tsv"
    grep -Ff "$keep_file" "$MAP_FILE" > "$pop_map" || true

    local pops
    pops=$(awk '{print $2}' "$pop_map" | sort -u)
    local pop_array=()
    while IFS= read -r p; do pop_array+=("$p"); done <<< "$pops"
    local n_pops=${#pop_array[@]}

    if [ "$n_pops" -lt 2 ]; then
        echo "WARN: fewer than 2 populations in group ${group}, skipping FST" | tee -a "$LOG_FILE"
        return
    fi

    # Write per-population sample files
    local pop_list_dir="${group_dir}/pop_lists"
    mkdir -p "$pop_list_dir"
    for pop in "${pop_array[@]}"; do
        local safe
        safe=$(echo "$pop" | sed 's/[^A-Za-z0-9._-]/_/g')
        awk -v p="$pop" '$2==p {print $1"\t"$1}' "$pop_map" \
            > "${pop_list_dir}/${safe}.txt"
    done

    # Pairwise FST: iterate over all pairs
    local fst_summary="${group_dir}/pairwise_fst_summary.tsv"
    echo -e "pop1\tpop2\tmean_fst\tweighted_fst" > "$fst_summary"

    for (( i=0; i<n_pops; i++ )); do
        for (( j=i+1; j<n_pops; j++ )); do
            local p1="${pop_array[$i]}"
            local p2="${pop_array[$j]}"
            local s1 s2
            s1=$(echo "$p1" | sed 's/[^A-Za-z0-9._-]/_/g')
            s2=$(echo "$p2" | sed 's/[^A-Za-z0-9._-]/_/g')
            local pair_prefix="${group_dir}/${s1}_vs_${s2}"

            echo "  FST: ${p1} vs ${p2}" | tee -a "$LOG_FILE"
            vcftools \
                --gzvcf "$VCFTOOLS_VCF" \
                --keep "${pop_list_dir}/${s1}.txt" \
                --keep "${pop_list_dir}/${s2}.txt" \
                --weir-fst-pop "${pop_list_dir}/${s1}.txt" \
                --weir-fst-pop "${pop_list_dir}/${s2}.txt" \
                --out "$pair_prefix" >> "$LOG_FILE" 2>&1

            # Extract mean and weighted FST from log
            local mean_fst weighted_fst
            mean_fst=$(grep "^Mean Fst" "${pair_prefix}.log" 2>/dev/null \
                | awk '{print $NF}' || echo "NA")
            weighted_fst=$(grep "^Weighted Fst" "${pair_prefix}.log" 2>/dev/null \
                | awk '{print $NF}' || echo "NA")

            echo -e "${p1}\t${p2}\t${mean_fst}\t${weighted_fst}" >> "$fst_summary"
        done
    done

    echo "FST summary written: $fst_summary" | tee -a "$LOG_FILE"
}

run_fst_group "guiana_only" "$KEEP_GUIANA"

# -------------------------------------------------------------------
# STEP 5: R — Mantel tests and plots
# -------------------------------------------------------------------
echo "STEP 5: R Mantel tests and IBD plots" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$OUT_DIR" \
    "$MAP_FILE" \
    "$GEOLOC_FILE" \
    "$N_PERM" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 07.4 IBD analysis completed" | tee -a "$LOG_FILE"
echo "Outputs: $OUT_DIR"

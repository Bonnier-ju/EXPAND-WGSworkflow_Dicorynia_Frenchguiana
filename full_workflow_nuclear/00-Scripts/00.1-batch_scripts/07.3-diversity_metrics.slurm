#!/bin/bash
#SBATCH -p workq
#SBATCH --time=12:00:00
#SBATCH -J 07.3-diversity_metrics
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/VCFtools/0.1.16
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
# Diversity metrics are computed on the TE-filtered VCF (06.5) —
# NOT the LD-pruned set (06.6b). Rare variants (low MAF) must be
# retained for accurate estimation of π, Tajima's D, and private
# alleles. LD pruning would remove rare variants and bias these metrics.
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
# VCF for vcftools (pi, Tajima's D, private alleles) — all SNPs, rare variants retained
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
# VCF for hierfstat (Ho, He, FIS, Ar) — LD-pruned independent markers
VCF_66="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
SAMPLE_SHEET="${SAMPLE_SHEET:-${BASE}/metadata/sample_sheet_complete_filtered.csv}"
POP_COLUMN="${POP_COLUMN:-site}"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.3-diversity_metrics/07.3.1-metrics_by_sites"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
POP_DIR="${OUT_DIR}/per_population"
LOG_FILE="${LOG_DIR}/07.3-diversity_metrics.log"

R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.3-diversity_metrics.R"

# Window size (bp) for Tajima's D calculation.
# 10 kb is a standard window for nuclear SNP data in trees.
TAJIMA_WINDOW="${TAJIMA_WINDOW:-10000}"

# Minimum number of individuals required per site to include in
# allelic richness rarefaction (sites with N < MIN_N_AR are computed
# for all other metrics but excluded from the Ar rarefaction).
MIN_N_AR="${MIN_N_AR:-4}"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR" "$POP_DIR"
echo "START $(date -Iseconds) 07.3 diversity metrics" > "$LOG_FILE"
echo "Input VCF: $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Sample sheet: $SAMPLE_SHEET (POP_COLUMN=$POP_COLUMN)" | tee -a "$LOG_FILE"
echo "Tajima D window: ${TAJIMA_WINDOW} bp | Min N for Ar: ${MIN_N_AR}" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: VCFv4.3 → v4.2 patch for vcftools compatibility
# -------------------------------------------------------------------
# vcftools 0.1.16 only supports up to VCFv4.2. Patch the header once
# and use the patched copy for all subsequent vcftools calls.
VCF_VERSION=$(bcftools view -h "$INPUT_VCF" | grep "^##fileformat=" | sed 's/##fileformat=//')
VCFTOOLS_VCF="$INPUT_VCF"
if [ "$VCF_VERSION" = "VCFv4.3" ]; then
    VCFTOOLS_VCF="${TMP_DIR}/input_v42.vcf.gz"
    if [ ! -s "$VCFTOOLS_VCF" ]; then
        echo "INFO: patching VCFv4.3 → v4.2 for vcftools compatibility" | tee -a "$LOG_FILE"
        zcat "$INPUT_VCF" \
            | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
            | bgzip -c > "$VCFTOOLS_VCF"
        bcftools index -t "$VCFTOOLS_VCF" >> "$LOG_FILE" 2>&1
    fi
    echo "INFO: vcftools will use patched VCF: $VCFTOOLS_VCF" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 3: Build population sample lists from metadata
# -------------------------------------------------------------------
echo "STEP 3: building population sample lists" | tee -a "$LOG_FILE"

VCF_SAMPLES="${TMP_DIR}/vcf_samples.txt"
bcftools query -l "$INPUT_VCF" | sort -u > "$VCF_SAMPLES"

MAP_FILE="${OUT_DIR}/sample_population_map.tsv"
awk -F',' -v popcol="$POP_COLUMN" '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="sample_id") sid=i
        if ($i==popcol) pop=i
    }
    if (!sid || !pop) { print "ERROR: columns not found" > "/dev/stderr"; exit 2 }
    next
}
{
    gsub(/\r/, "", $sid); gsub(/\r/, "", $pop)
    if ($sid != "" && $pop != "" && $pop != "NA") print $sid "\t" $pop
}
' "$SAMPLE_SHEET" | sort -k1,1 | join -t $'\t' -1 1 -2 1 \
    - <(awk '{print $1"\t"$1}' "$VCF_SAMPLES" | sort -k1,1) \
    | awk -F'\t' '{print $1"\t"$2}' | sort -k1,1 > "$MAP_FILE"

# Exclude outgroup populations and treemutation project samples
# Outgroup populations: Cameroun_Benin, Venezuela
# Treemutation individual IDs: DBR_AGJ DBR_AGK DBR_AGL DBR_AGM
EXCLUDE_POPS="Cameroun_Benin Venezuela"
TREEMUTATION_IDS="DBR_AGJ DBR_AGK DBR_AGL DBR_AGM"
POP_EXCL_PAT=$(echo "$EXCLUDE_POPS"    | tr ' ' '\n' | awk '{print "^"$0"$"}' | paste -sd'|')
ID_EXCL_PAT=$(echo  "$TREEMUTATION_IDS" | tr ' ' '\n' | awk '{print "^"$0"$"}' | paste -sd'|')
awk -v epop="$POP_EXCL_PAT" -v eid="$ID_EXCL_PAT" \
    '$2 !~ epop && $1 !~ eid' "$MAP_FILE" > "${TMP_DIR}/map_guiana_only.tsv"
mv "${TMP_DIR}/map_guiana_only.tsv" "$MAP_FILE"

POPS_LIST="${TMP_DIR}/populations.list"
awk -F'\t' '{print $2}' "$MAP_FILE" | sort -u > "$POPS_LIST"
N_POPS=$(wc -l < "$POPS_LIST")
echo "INFO: ${N_POPS} populations detected (outgroups and treemutation excluded)" | tee -a "$LOG_FILE"

# Per-population keep files (sample_id  sample_id for vcftools)
while IFS= read -r pop; do
    safe=$(echo "$pop" | sed 's/[^A-Za-z0-9._-]/_/g')
    awk -F'\t' -v p="$pop" '$2==p{print $1}' "$MAP_FILE" \
        > "${POP_DIR}/${safe}.samples.txt"
done < "$POPS_LIST"

# -------------------------------------------------------------------
# STEP 4: Per-population vcftools statistics
# -------------------------------------------------------------------
# For each population: π (per site), individual heterozygosity (--het),
# allele frequencies (--freq2), and Tajima's D (windowed).
# vcftools outputs are written to POP_DIR/<safe_pop>/ subdirectories.
echo "STEP 4: per-population vcftools stats" | tee -a "$LOG_FILE"

while IFS= read -r pop; do
    safe=$(echo "$pop" | sed 's/[^A-Za-z0-9._-]/_/g')
    keep_file="${POP_DIR}/${safe}.samples.txt"
    n_ind=$(wc -l < "$keep_file")
    pop_out_dir="${POP_DIR}/${safe}"
    mkdir -p "$pop_out_dir"

    echo "  Population: ${pop} (n=${n_ind})" | tee -a "$LOG_FILE"

    # Nucleotide diversity π per site
    vcftools \
        --gzvcf "$VCFTOOLS_VCF" \
        --keep "$keep_file" \
        --site-pi \
        --out "${pop_out_dir}/pi" >> "$LOG_FILE" 2>&1

    # Per-individual observed heterozygosity (F statistics)
    vcftools \
        --gzvcf "$VCFTOOLS_VCF" \
        --keep "$keep_file" \
        --het \
        --out "${pop_out_dir}/het" >> "$LOG_FILE" 2>&1

    # Allele frequencies (for He, private alleles, Ar)
    vcftools \
        --gzvcf "$VCFTOOLS_VCF" \
        --keep "$keep_file" \
        --freq2 \
        --out "${pop_out_dir}/freq" >> "$LOG_FILE" 2>&1

    # Tajima's D (windowed)
    vcftools \
        --gzvcf "$VCFTOOLS_VCF" \
        --keep "$keep_file" \
        --TajimaD "$TAJIMA_WINDOW" \
        --out "${pop_out_dir}/tajima" >> "$LOG_FILE" 2>&1

done < "$POPS_LIST"

# -------------------------------------------------------------------
# STEP 4: R — compile metrics and produce summary table + plots
# -------------------------------------------------------------------
echo "STEP 5: R compilation of diversity metrics" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "$OUT_DIR" \
    "$MAP_FILE" \
    "$POPS_LIST" \
    "$MIN_N_AR" \
    "$VCF_66" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 07.3 diversity metrics completed" | tee -a "$LOG_FILE"
echo "Outputs: $OUT_DIR"

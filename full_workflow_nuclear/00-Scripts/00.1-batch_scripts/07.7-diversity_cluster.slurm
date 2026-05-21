#!/bin/bash
#SBATCH -p workq
#SBATCH --time=12:00:00
#SBATCH -J 07.7-diversity_cluster
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   K_GLOBAL=2 sbatch 07.7-diversity_cluster.slurm
#   K_GLOBAL=3 sbatch 07.7-diversity_cluster.slurm
#
# Parameters:
#   K_GLOBAL      — global ADMIXTURE K defining the populations (2, 3 or 4)
#   TAJIMA_WINDOW — window size in bp for Tajima's D (default 10000)
#   MIN_N_AR      — minimum N for allelic richness rarefaction (default 4)
#
# Computes diversity metrics treating each genetic cluster as one population:
#   - nucleotide diversity pi (vcftools --site-pi)
#   - individual heterozygosity (vcftools --het)
#   - allele frequencies (vcftools --freq2)
#   - Tajima's D windowed (vcftools --TajimaD)
#   - Ho, He, Fis, allelic richness, private alleles (R / hierfstat)
#   - pairwise Fst between clusters (R / hierfstat)
#
# Output: 07.3-diversity_metrics/07.3.${K_GLOBAL}-K=${K_GLOBAL}/
# -------------------------------------------------------------------

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
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL="${K_GLOBAL:-2}"
TAJIMA_WINDOW="${TAJIMA_WINDOW:-10000}"
MIN_N_AR="${MIN_N_AR:-4}"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
VCF_66="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"

R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.3-diversity_metrics.R"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.3-diversity_metrics/07.3.${K_GLOBAL}-K=${K_GLOBAL}"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
POP_DIR="${OUT_DIR}/per_cluster"
MAP_FILE="${OUT_DIR}/sample_cluster_map.tsv"
POPS_LIST="${TMP_DIR}/clusters.list"
LOG_FILE="${LOG_DIR}/07.7-diversity_K${K_GLOBAL}.log"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR" "$POP_DIR"
echo "START $(date -Iseconds) 07.7 diversity by cluster K=${K_GLOBAL}" > "$LOG_FILE"
echo "Input VCF : $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Cluster assignments: $CLUSTER_ASSIGNMENTS" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"; exit 1
fi
if [ ! -s "$CLUSTER_ASSIGNMENTS" ]; then
    echo "ERROR: cluster_assignments.tsv not found — run 07.0-prep_cluster_assignments.R first" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build sample -> cluster map from cluster_assignments.tsv
# Map format: sample_id <tab> Pop_X  (no header, two columns)
# -------------------------------------------------------------------
K_COL="K${K_GLOBAL}"
echo "STEP 2: building sample -> cluster map (column ${K_COL})" | tee -a "$LOG_FILE"

awk -v kcol="$K_COL" '
NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == kcol) col = i
    if (!col) { print "ERROR: column " kcol " not found" > "/dev/stderr"; exit 2 }
    next
}
$col != "" && $col != "NA" { print $1 "\t" $col }
' "$CLUSTER_ASSIGNMENTS" > "$MAP_FILE"

if [ ! -s "$MAP_FILE" ]; then
    echo "ERROR: empty map — check cluster_assignments.tsv" | tee -a "$LOG_FILE"; exit 1
fi

awk '{print $2}' "$MAP_FILE" | sort -u > "$POPS_LIST"
N_POPS=$(wc -l < "$POPS_LIST")
echo "INFO: ${N_POPS} clusters detected for K=${K_GLOBAL}" | tee -a "$LOG_FILE"
cat "$POPS_LIST" | tee -a "$LOG_FILE"

# Per-cluster sample files
while IFS= read -r pop; do
    awk -v p="$pop" '$2 == p { print $1 }' "$MAP_FILE" > "${POP_DIR}/${pop}.samples.txt"
    N=$(wc -l < "${POP_DIR}/${pop}.samples.txt")
    echo "  ${pop}: ${N} individuals" | tee -a "$LOG_FILE"
done < "$POPS_LIST"

# -------------------------------------------------------------------
# STEP 3: VCFv4.3 -> v4.2 patch for vcftools compatibility
# -------------------------------------------------------------------
VCF_VERSION=$(bcftools view -h "$INPUT_VCF" | grep "^##fileformat=" | sed 's/##fileformat=//')
VCFTOOLS_VCF="$INPUT_VCF"
if [ "$VCF_VERSION" = "VCFv4.3" ]; then
    VCFTOOLS_VCF="${TMP_DIR}/input_v42.vcf.gz"
    if [ ! -s "$VCFTOOLS_VCF" ]; then
        echo "INFO: patching VCFv4.3 → v4.2 for vcftools" | tee -a "$LOG_FILE"
        zcat "$INPUT_VCF" \
            | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
            | bgzip -c > "$VCFTOOLS_VCF"
        bcftools index -t "$VCFTOOLS_VCF" >> "$LOG_FILE" 2>&1
    fi
fi

# -------------------------------------------------------------------
# STEP 4: Per-cluster vcftools statistics
# -------------------------------------------------------------------
echo "STEP 4: per-cluster vcftools stats" | tee -a "$LOG_FILE"

while IFS= read -r pop; do
    keep_file="${POP_DIR}/${pop}.samples.txt"
    n_ind=$(wc -l < "$keep_file")
    pop_out="${POP_DIR}/${pop}"
    mkdir -p "$pop_out"
    echo "  ${pop} (n=${n_ind})" | tee -a "$LOG_FILE"

    vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$keep_file" \
        --site-pi --out "${pop_out}/pi" >> "$LOG_FILE" 2>&1

    vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$keep_file" \
        --het --out "${pop_out}/het" >> "$LOG_FILE" 2>&1

    vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$keep_file" \
        --freq2 --out "${pop_out}/freq" >> "$LOG_FILE" 2>&1

    vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$keep_file" \
        --TajimaD "$TAJIMA_WINDOW" --out "${pop_out}/tajima" >> "$LOG_FILE" 2>&1

done < "$POPS_LIST"

# -------------------------------------------------------------------
# STEP 5: R — compile metrics + plots
# Reuses the existing 07.3-diversity_metrics.R:
#   arg1 = out_dir       (07.3.K-K=K directory)
#   arg2 = map_file      (sample_id -> cluster_name)
#   arg3 = pops_list     (list of cluster names)
#   arg4 = min_n_ar      (min N for Ar rarefaction)
#   arg5 = vcf_ld_pruned (LD-pruned VCF for hierfstat)
# The R script reads per_population/ subdirectories — we pass POP_DIR
# via a symlink so the script finds them under the expected path.
# -------------------------------------------------------------------
echo "STEP 5: R diversity metrics compilation" | tee -a "$LOG_FILE"

# Expose per-cluster stats under the path the R script expects
POP_SYMLINK="${OUT_DIR}/per_population"
[ -L "$POP_SYMLINK" ] && rm "$POP_SYMLINK"
ln -s "$POP_DIR" "$POP_SYMLINK"

Rscript --vanilla "$R_SCRIPT" \
    "$OUT_DIR" \
    "$MAP_FILE" \
    "$POPS_LIST" \
    "$MIN_N_AR" \
    "$VCF_66" >> "$LOG_FILE" 2>&1 \
    || echo "WARN: R script returned non-zero — check log" | tee -a "$LOG_FILE"

echo "DONE $(date -Iseconds) 07.7 diversity K=${K_GLOBAL} completed" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}/"

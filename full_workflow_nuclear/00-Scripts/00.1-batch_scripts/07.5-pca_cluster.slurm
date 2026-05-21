#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J 07.5-pca_cluster
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   K_GLOBAL=2 POP_ID=1 sbatch 07.5-pca_cluster.slurm
#   K_GLOBAL=3 POP_ID=2 sbatch 07.5-pca_cluster.slurm
#
# Parameters:
#   K_GLOBAL  — global ADMIXTURE K used to define populations (2, 3 or 4)
#   POP_ID    — population index to analyse (1 .. n_pop for this K)
#
# Runs PCA on the subset of guiana_only individuals belonging to
# Pop_${POP_ID} under the K=${K_GLOBAL} grouping.
# Output: 07.1-pca_ibd_fst/07.1.${K_GLOBAL}-K=${K_GLOBAL}/Pop_${POP_ID}/
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/1.90b7
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL="${K_GLOBAL:-4}"
POP_ID="${POP_ID:-4}"
POP_LABEL="Pop_${POP_ID}"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
VCF_66="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
VCF_65="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"
SITE_MAP="${BASE}/07-population_structure/workflow_full/07.1-pca_ibd_fst/07.1.1-K=1/sample_population_map.tsv"

R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.1-pca_plot.R"

if [ -s "$VCF_66" ]; then
    INPUT_VCF="$VCF_66"; INPUT_STEP="06.6b"
else
    INPUT_VCF="$VCF_65"; INPUT_STEP="06.5"
fi

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.1-pca_ibd_fst/07.1.${K_GLOBAL}-K=${K_GLOBAL}/${POP_LABEL}"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/07.5-pca_K${K_GLOBAL}_${POP_LABEL}.log"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 07.5 PCA cluster K=${K_GLOBAL} ${POP_LABEL}" > "$LOG_FILE"
echo "Input VCF : $INPUT_VCF ($INPUT_STEP)" | tee -a "$LOG_FILE"
echo "Cluster assignments: $CLUSTER_ASSIGNMENTS" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"; exit 1
fi
if [ ! -s "$CLUSTER_ASSIGNMENTS" ]; then
    echo "ERROR: cluster_assignments.tsv not found — run 07.0-prep_cluster_assignments.R first" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build keep file for this K_GLOBAL / POP_LABEL
# -------------------------------------------------------------------
# cluster_assignments.tsv columns: id_genoscope site K2 K3 K4
# Select column K${K_GLOBAL} and filter rows matching Pop_${POP_ID}
K_COL="K${K_GLOBAL}"
KEEP_FILE="${TMP_DIR}/keep_${POP_LABEL}.txt"

awk -v kcol="$K_COL" -v pop="$POP_LABEL" '
NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == kcol) col = i
    if (!col) { print "ERROR: column " kcol " not found" > "/dev/stderr"; exit 2 }
    next
}
$col == pop { print $1 "\t" $1 }
' "$CLUSTER_ASSIGNMENTS" > "$KEEP_FILE"

N_IND=$(wc -l < "$KEEP_FILE")
echo "Individuals in K=${K_GLOBAL} ${POP_LABEL}: ${N_IND}" | tee -a "$LOG_FILE"

if [ "$N_IND" -lt 3 ]; then
    echo "ERROR: fewer than 3 individuals — cannot run PCA" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 3: PLINK VCF -> BED (subset to this population)
# -------------------------------------------------------------------
PLINK_PREFIX="${OUT_DIR}/genotypes_${POP_LABEL}"
echo "STEP 3: PLINK conversion" | tee -a "$LOG_FILE"
plink \
    --vcf "$INPUT_VCF" \
    --allow-extra-chr \
    --keep "$KEEP_FILE" \
    --double-id \
    --set-missing-var-ids @:#\$1,\$2 \
    --make-bed \
    --out "$PLINK_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "${PLINK_PREFIX}.bed" ]; then
    echo "ERROR: PLINK conversion failed" | tee -a "$LOG_FILE"; exit 1
fi
N_SNP=$(wc -l < "${PLINK_PREFIX}.bim")
echo "INFO: ${N_IND} individuals, ${N_SNP} SNPs" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: PCA
# -------------------------------------------------------------------
N_PC=$(( N_IND < 20 ? N_IND - 1 : 20 ))
PCA_PREFIX="${OUT_DIR}/pca_${POP_LABEL}"
echo "STEP 4: PCA — ${N_PC} axes" | tee -a "$LOG_FILE"
plink \
    --bfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --pca ${N_PC} header tabs \
    --out "$PCA_PREFIX" >> "$LOG_FILE" 2>&1

echo "INFO: PCA done → ${PCA_PREFIX}.eigenvec" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 5: R plotting
# The existing 07.1-pca_plot.R is reused:
#   arg1 = output base dir (plots written to out_dir/plots/<group_label>/)
#   arg2 = pop map file (sample_id -> site, for coloring by site)
#   arg3 = eigenvec file
#   arg4 = eigenval file
#   arg5 = group_label (used as title and subfolder name)
# -------------------------------------------------------------------
if [ -s "$R_PLOT_SCRIPT" ]; then
    GROUP_LABEL="K${K_GLOBAL}_${POP_LABEL}"
    echo "STEP 5: R plotting — ${GROUP_LABEL}" | tee -a "$LOG_FILE"
    Rscript --vanilla "$R_PLOT_SCRIPT" \
        "$OUT_DIR" \
        "$SITE_MAP" \
        "${PCA_PREFIX}.eigenvec" \
        "${PCA_PREFIX}.eigenval" \
        "$GROUP_LABEL" >> "$LOG_FILE" 2>&1 \
        || echo "WARN: R plotting failed" | tee -a "$LOG_FILE"
else
    echo "WARN: R plot script not found: $R_PLOT_SCRIPT" | tee -a "$LOG_FILE"
fi

echo "DONE $(date -Iseconds) 07.5 PCA K=${K_GLOBAL} ${POP_LABEL}" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}/plots/K${K_GLOBAL}_${POP_LABEL}/"

#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J 07b.3-diversity_metrics
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
# STEP 1: Paths and parameters
# -------------------------------------------------------------------
BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
WD="${BASE}/07b-population_structure_pur_inds"

# Full-SNP VCF (06.5, all SNPs, no LD pruning) — for vcftools metrics
VCF_ALL="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
# LD-pruned VCF (06.6b, independent markers) — for hierfstat Ho/He/FIS/Ar
VCF_LD="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"

SAMPLE_SHEET="${BASE}/metadata/sample_sheet_complete_filtered.csv"
SITES_BY_CLUSTERS="${BASE}/metadata/sites_by_clusters.csv"
LISTS_DIR="${WD}/07b.0-purity_filtering/individual_lists"
OUT_BASE="${WD}/07b.3-diversity_metrics"
LOG_FILE="${WD}/logs/07b.3-diversity_metrics.log"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07b.3-diversity_metrics.R"

# Rarefaction sample size for Ar and PAR
MIN_N_AR=10

# Tajima's D sliding window (bp) — standard for nuclear data in trees
TAJIMA_WINDOW=10000

for f in "$VCF_ALL" "$VCF_LD" "$SAMPLE_SHEET" "$SITES_BY_CLUSTERS" "$R_SCRIPT"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: missing required input: $f"
        exit 1
    fi
done

TMP_DIR="${OUT_BASE}/tmp"
mkdir -p "$OUT_BASE" "$TMP_DIR" "${WD}/logs"
echo "START $(date -Iseconds) 07b.3 diversity metrics on purity-filtered individuals" > "$LOG_FILE"
echo "Parameters: MIN_N_AR=${MIN_N_AR} | TAJIMA_WINDOW=${TAJIMA_WINDOW}" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Build site → cluster (K=3) mapping
#   Source: sites_by_clusters.csv, column K=3 (col 3)
#   Excluded: St_georges, Petit_Croissant
# -------------------------------------------------------------------
SITE_CLUSTER="${TMP_DIR}/site_cluster_k3.tsv"
awk -F',' 'NR>1 && $1 != "St_georges" && $1 != "Petit_Croissant" {print $1"\t"$3}' \
    "$SITES_BY_CLUSTERS" > "$SITE_CLUSTER"
echo "Site→cluster map (K=3): $(wc -l < "$SITE_CLUSTER") sites" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Build full sample → cluster map
#   sample_sheet ($1=sample_id, $4=site) × site_cluster (TSV)
#   Excluded individuals: DBR_ADT + outgroups + Angela/treemutation
# -------------------------------------------------------------------
FULL_SAMPLE_MAP="${TMP_DIR}/full_sample_cluster_map.tsv"
awk '
NR==FNR {
    clust[$1] = $2; next
}
FNR==1 { next }
{
    n = split($0, f, ",")
    id = f[1]; site = f[4]
    gsub(/\r/, "", id); gsub(/\r/, "", site)
    if (id in excl) next
    if (site in clust) print id"\t"clust[site]
}
BEGIN {
    split("DBR_ADT DBR_AGJ DBR_AGK DBR_AGL DBR_AGM DBR_AGN DBR_AGP DBR_AGQ DBR_AGR",
          arr, " ")
    for (i in arr) excl[arr[i]] = 1
}
' "$SITE_CLUSTER" "$SAMPLE_SHEET" > "$FULL_SAMPLE_MAP"
echo "Full sample→cluster map: $(wc -l < "$FULL_SAMPLE_MAP") individuals" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: VCF version patch for vcftools (v4.3 → v4.2, done once)
# -------------------------------------------------------------------
VCF_VERSION=$(bcftools view -h "$VCF_ALL" | grep "^##fileformat=" | sed 's/##fileformat=//')
VCFTOOLS_VCF="$VCF_ALL"
if [ "$VCF_VERSION" = "VCFv4.3" ]; then
    VCFTOOLS_VCF="${TMP_DIR}/input_v42.vcf.gz"
    if [ ! -s "$VCFTOOLS_VCF" ]; then
        echo "INFO: patching VCFv4.3 → v4.2 for vcftools" | tee -a "$LOG_FILE"
        zcat "$VCF_ALL" \
            | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
            | bgzip -c > "$VCFTOOLS_VCF"
        bcftools index -t "$VCFTOOLS_VCF" >> "$LOG_FILE" 2>&1
    fi
fi

# -------------------------------------------------------------------
# STEP 5: Loop over purity thresholds
# -------------------------------------------------------------------
for THRESH in 70 80 90 95; do
    TAG="T${THRESH}"
    KEEP_RAW="${LISTS_DIR}/inds_${TAG}_all.txt"

    if [ ! -s "$KEEP_RAW" ]; then
        echo "WARN: keep file not found for $TAG — skipping" | tee -a "$LOG_FILE"
        continue
    fi

    OUT_T="${OUT_BASE}/${TAG}"
    CLUSTER_DIR="${OUT_T}/per_cluster"
    mkdir -p "$CLUSTER_DIR"

    # Skip threshold if R plots already exist (threshold fully complete)
    if ls "${OUT_T}/plots/"*.png 2>/dev/null | grep -q .; then
        echo "=== $TAG: already completed — skipping ===" | tee -a "$LOG_FILE"
        continue
    fi

    echo "=== $TAG ===" | tee -a "$LOG_FILE"

    # STEP 5a: Extract purity-filtered IDs (col 2 of keep file)
    KEEP_IDS="${OUT_T}/keep_ids.txt"
    cut -f2 "$KEEP_RAW" > "$KEEP_IDS"

    # STEP 5b: Build per-threshold sample → cluster map
    #   Intersect purity-filtered IDs with full sample→cluster map
    SAMPLE_MAP="${OUT_T}/sample_cluster_map.tsv"
    awk 'NR==FNR {keep[$1]=1; next} $1 in keep' \
        "$KEEP_IDS" "$FULL_SAMPLE_MAP" > "$SAMPLE_MAP"
    N_MAPPED=$(wc -l < "$SAMPLE_MAP")
    echo "  Mapped individuals: $N_MAPPED" | tee -a "$LOG_FILE"

    # STEP 5c: Build per-cluster sample lists (IDs only, for vcftools --keep)
    for CLUSTER in Pop_1 Pop_2 Pop_3; do
        awk -v c="$CLUSTER" '$2==c {print $1}' "$SAMPLE_MAP" \
            > "${CLUSTER_DIR}/${CLUSTER}.samples.txt"
        N_C=$(wc -l < "${CLUSTER_DIR}/${CLUSTER}.samples.txt")
        echo "  $CLUSTER: $N_C individuals" | tee -a "$LOG_FILE"
    done

    # STEP 5d: Per-cluster vcftools statistics
    echo "  Running vcftools per cluster..." | tee -a "$LOG_FILE"
    for CLUSTER in Pop_1 Pop_2 Pop_3; do
        KEEP_C="${CLUSTER_DIR}/${CLUSTER}.samples.txt"
        N_C=$(wc -l < "$KEEP_C")
        if [ "$N_C" -eq 0 ]; then
            echo "  WARN: $CLUSTER has 0 individuals — skipping vcftools" | tee -a "$LOG_FILE"
            continue
        fi
        POP_OUT="${CLUSTER_DIR}/${CLUSTER}"
        mkdir -p "$POP_OUT"

        # Skip vcftools if all output files already exist
        if [ -s "${POP_OUT}/pi.sites.pi" ] && [ -s "${POP_OUT}/het.het" ] && \
           [ -s "${POP_OUT}/freq.frq" ]    && [ -s "${POP_OUT}/tajima.Tajima.D" ]; then
            echo "  vcftools already done: $CLUSTER — skipping" | tee -a "$LOG_FILE"
        else
            # Nucleotide diversity π per site
            vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$KEEP_C" \
                --site-pi --out "${POP_OUT}/pi" >> "$LOG_FILE" 2>&1

            # Per-individual observed heterozygosity (F statistics)
            vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$KEEP_C" \
                --het --out "${POP_OUT}/het" >> "$LOG_FILE" 2>&1

            # Allele frequencies — N_CHR and freq retained for private allele rarefaction
            vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$KEEP_C" \
                --freq2 --out "${POP_OUT}/freq" >> "$LOG_FILE" 2>&1

            # Tajima's D (windowed)
            vcftools --gzvcf "$VCFTOOLS_VCF" --keep "$KEEP_C" \
                --TajimaD "$TAJIMA_WINDOW" --out "${POP_OUT}/tajima" >> "$LOG_FILE" 2>&1

            echo "  vcftools done: $CLUSTER ($N_C ind)" | tee -a "$LOG_FILE"
        fi
    done

    # STEP 5e: R — compile metrics and produce summary table + plots
    echo "  R compilation — $TAG" | tee -a "$LOG_FILE"
    Rscript --vanilla "$R_SCRIPT" \
        "$OUT_T" \
        "$SAMPLE_MAP" \
        "$MIN_N_AR" \
        "$VCF_LD" \
        "$CLUSTER_DIR" \
        "$TAG" >> "$LOG_FILE" 2>&1 \
        || echo "  WARN: R script failed for $TAG — check log" | tee -a "$LOG_FILE"

    echo "  DONE $TAG — plots: ${OUT_T}/plots/" | tee -a "$LOG_FILE"
done

echo "DONE $(date -Iseconds) 07b.3 diversity metrics completed" | tee -a "$LOG_FILE"
echo "Outputs: $OUT_BASE"

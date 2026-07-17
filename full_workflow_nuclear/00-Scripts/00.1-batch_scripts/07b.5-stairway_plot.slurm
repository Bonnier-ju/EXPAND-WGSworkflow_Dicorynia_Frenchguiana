#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=180:00:00
#SBATCH -J 07b.5-stairway_plot
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   sbatch 07b.5-stairway_plot.slurm
#
# Demographic inference via Stairway Plot 2 for K=3 genetic clusters
# (Pop_1 / Pop_2 / Pop_3) on purity-filtered individuals at four
# thresholds: T70 / T80 / T90 / T95.
# Sites St_georges and Petit_Croissant excluded from all analyses.
#
# Pipeline per threshold × population:
#   STEP 4a : Build per-cluster sample list (purity filter + cluster)
#   STEP 4b : Extract per-cluster VCF (within-cluster polymorphic)
#   STEP 4c : easySFS preview (logged for reference)
#   STEP 4d : easySFS folded SFS (full 2n projection)
#   STEP 4e : Parse SFS → write Stairway Plot 2 blueprint
#   STEP 4f : Stairway Plot phase 1 (Stairbuilder → training script)
#   STEP 4g : Stairway Plot phase 2 (training replicates via sarray)
#   STEP 4h : Stairway Plot phase 3 (summary + built-in plot)
#
# L values (callable non-TE sites, GQ >= 20) reused from 07.11b:
#   K=3 populations are identical in site composition to the
#   corresponding K=4 populations (Pop_4 = St_georges + Petit_Croissant,
#   which are excluded here, so no overlap):
#     Pop_1: 215,357,483 bp  |  Pop_2: 221,273,354 bp
#     Pop_3: 193,750,488 bp
#
# Mutation rate: mu=5e-9 per site per generation (central estimate;
#   R script shows sensitivity for mu = 2..8 x 10^-9).
#
# Site independence: LD-pruned VCF (06.6b) used (same as 07.11).
#
# Populations run SEQUENTIALLY per threshold (memory safety for easySFS).
# Training replicates run in PARALLEL via sarray (200 per population).
#
# Output: 07b-population_structure_pur_inds/07b.5-stairway_plot/
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/easySFS/0.0.1
module load bioinfo/stairway_plot/2.1.2
module load devel/java/17.0.6
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL=3

# L per population (non-TE callable sites; same groupings as K=4 Pop_1/2/3)
declare -A L_POP=([Pop_1]=215357483 [Pop_2]=221273354 [Pop_3]=193750488)
MU_REF="5e-9"
YEAR_PER_GEN=30
NINPUT=200
PCT_TRAIN=0.67

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
WD="${BASE}/07b-population_structure_pur_inds"
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
SAMPLE_SHEET="${BASE}/metadata/sample_sheet_complete_filtered.csv"
SITES_BY_CLUSTERS="${BASE}/metadata/sites_by_clusters.csv"
LISTS_DIR="${WD}/07b.0-purity_filtering/individual_lists"
STAIRWAY_ES_DIR="/usr/local/bioinfo/src/stairway_plot/stairway_plot_v2.1.2/stairway_plot_es"
R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07b.5-stairway_plot.R"
OUT_BASE="${WD}/07b.5-stairway_plot"
TMP_DIR="${OUT_BASE}/tmp"

mkdir -p "$OUT_BASE" "$TMP_DIR" "${WD}/logs"
LOG_FILE="${WD}/logs/07b.5-stairway_plot.log"
echo "START $(date -Iseconds) 07b.5 Stairway Plot K=${K_GLOBAL}" > "$LOG_FILE"
echo "L callable: Pop_1=${L_POP[Pop_1]} Pop_2=${L_POP[Pop_2]} Pop_3=${L_POP[Pop_3]}" | tee -a "$LOG_FILE"
echo "mu_ref: ${MU_REF} | year_per_gen: ${YEAR_PER_GEN}" | tee -a "$LOG_FILE"

for f in "$INPUT_VCF" "$SAMPLE_SHEET" "$SITES_BY_CLUSTERS"; do
    [ -s "$f" ] || { echo "ERROR: not found: $f" | tee -a "$LOG_FILE"; exit 1; }
done

# -------------------------------------------------------------------
# STEP 2: Site → cluster (K=3) mapping
#   Column 3 of sites_by_clusters.csv holds K=3 assignments.
#   St_georges and Petit_Croissant excluded from all 07b analyses.
# -------------------------------------------------------------------
SITE_CLUSTER="${TMP_DIR}/site_cluster_k3.tsv"
awk -F',' 'NR>1 && $1!="St_georges" && $1!="Petit_Croissant" {print $1"\t"$3}' \
    "$SITES_BY_CLUSTERS" > "$SITE_CLUSTER"
echo "Site→cluster map (K=3): $(wc -l < "$SITE_CLUSTER") sites" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Full sample → cluster map (excluding outlier individuals)
# -------------------------------------------------------------------
FULL_SAMPLE_MAP="${TMP_DIR}/full_sample_cluster_map.tsv"
awk '
NR==FNR { clust[$1]=$2; next }
FNR==1  { next }
{
    n=split($0, f, ","); id=f[1]; site=f[4]
    gsub(/\r/,"",id); gsub(/\r/,"",site)
    if (id in excl) next
    if (site in clust) print id"\t"clust[site]
}
BEGIN {
    split("DBR_ADT DBR_AGJ DBR_AGK DBR_AGL DBR_AGM DBR_AGN DBR_AGP DBR_AGQ DBR_AGR",
          arr, " ")
    for (i in arr) excl[arr[i]]=1
}
' "$SITE_CLUSTER" "$SAMPLE_SHEET" > "$FULL_SAMPLE_MAP"
echo "Full sample→cluster map: $(wc -l < "$FULL_SAMPLE_MAP") individuals" | tee -a "$LOG_FILE"

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
    mkdir -p "${OUT_T}/logs"

    # Threshold-level skip: R plots already written → fully done
    if ls "${OUT_T}/plots/"*.png 2>/dev/null | grep -q .; then
        echo "=== ${TAG}: already completed — skipping ===" | tee -a "$LOG_FILE"
        continue
    fi

    echo "=== ${TAG} ===" | tee -a "$LOG_FILE"

    # Per-threshold purity-filtered IDs (col 2 = sample ID)
    KEEP_IDS="${OUT_T}/keep_ids.txt"
    cut -f2 "$KEEP_RAW" > "$KEEP_IDS"

    # Intersect purity-filtered IDs with full cluster map
    SAMPLE_MAP="${OUT_T}/sample_cluster_map.tsv"
    awk 'NR==FNR {keep[$1]=1; next} $1 in keep' \
        "$KEEP_IDS" "$FULL_SAMPLE_MAP" > "$SAMPLE_MAP"
    echo "  Mapped individuals: $(wc -l < "$SAMPLE_MAP")" | tee -a "$LOG_FILE"

    # ------------------------------------------------------------------
    # Per-population pipeline (run sequentially)
    # ------------------------------------------------------------------
    for POP in Pop_1 Pop_2 Pop_3; do
        OUT_POP="${OUT_T}/${POP}"
        mkdir -p "${OUT_POP}/easySFS" "${OUT_POP}/logs"

        POP_LOG="${OUT_POP}/logs/stairway_${POP}_${TAG}.log"
        SUMMARY_FILE="${OUT_POP}/stairway/${POP} (${TAG}).final.summary"

        echo "=== START ${POP} ${TAG} $(date -Iseconds) ===" | tee -a "$POP_LOG" "$LOG_FILE"

        # STEP 4a: Build per-cluster sample list
        CLUSTER_SAMPLES="${OUT_POP}/logs/samples_${POP}_${TAG}.txt"
        awk -v pop="$POP" '$2==pop {print $1}' "$SAMPLE_MAP" > "$CLUSTER_SAMPLES"
        N_IND=$(wc -l < "$CLUSTER_SAMPLES")
        NPROJ=$((2 * N_IND))
        echo "  ${N_IND} individuals (projection: ${NPROJ})" | tee -a "$POP_LOG"

        if [ "$N_IND" -lt 4 ]; then
            echo "WARN: fewer than 4 individuals in ${POP} ${TAG} — skipping" | tee -a "$POP_LOG" "$LOG_FILE"
            continue
        fi

        # Pop-level skip: summary already exists
        if [ -s "$SUMMARY_FILE" ]; then
            echo "  Summary already exists — skipping ${POP} ${TAG}" | tee -a "$POP_LOG" "$LOG_FILE"
            continue
        fi

        # easySFS population map (sample_id<TAB>population)
        POP_MAP="${OUT_POP}/logs/pop_map_${POP}_${TAG}.txt"
        awk -v pop="$POP" '{ print $1 "\t" pop }' "$CLUSTER_SAMPLES" > "$POP_MAP"

        # STEP 4b: Extract per-cluster VCF (within-cluster polymorphic biallelic)
        POP_VCF="${OUT_POP}/${POP}_${TAG}.vcf.gz"
        if [ -s "$POP_VCF" ]; then
            echo "  VCF already exists, skipping extraction" | tee -a "$POP_LOG"
        else
            echo "  Extracting within-cluster VCF..." | tee -a "$POP_LOG"
            bcftools view \
                --samples-file "$CLUSTER_SAMPLES" \
                --min-ac 1:minor \
                --max-alleles 2 \
                -O z -o "$POP_VCF" \
                "$INPUT_VCF" >> "$POP_LOG" 2>&1
            bcftools index "$POP_VCF" >> "$POP_LOG" 2>&1
            N_SNPS=$(bcftools index -n "$POP_VCF")
            echo "  INFO: ${N_SNPS} within-cluster polymorphic SNPs" | tee -a "$POP_LOG"
        fi

        # STEP 4c: easySFS preview (logged for reference)
        PREVIEW_LOG="${OUT_POP}/logs/easySFS_preview_${POP}_${TAG}.log"
        echo "  Running easySFS preview..." | tee -a "$POP_LOG"
        easySFS.py \
            -i "$POP_VCF" \
            -p "$POP_MAP" \
            -a \
            --ploidy 2 \
            --preview \
            > "$PREVIEW_LOG" 2>&1 || true
        cat "$PREVIEW_LOG" | tee -a "$POP_LOG"

        # STEP 4d: easySFS folded SFS (full 2n projection)
        # -a: keep all SNPs; folded by default (no --unfolded).
        EASYSFS_OUT="${OUT_POP}/easySFS"
        SFS_FILE="${EASYSFS_OUT}/dadi/${POP}-${NPROJ}.sfs"

        if [ -s "$SFS_FILE" ]; then
            echo "  easySFS output already exists, skipping" | tee -a "$POP_LOG"
        else
            echo "  Running easySFS (proj=${NPROJ})..." | tee -a "$POP_LOG"
            easySFS.py \
                -i "$POP_VCF" \
                -p "$POP_MAP" \
                -a \
                --ploidy 2 \
                --proj "${NPROJ}" \
                -f \
                -o "$EASYSFS_OUT" \
                >> "$POP_LOG" 2>&1
        fi

        if [ ! -s "$SFS_FILE" ]; then
            echo "ERROR: easySFS output not found: $SFS_FILE" | tee -a "$POP_LOG" "$LOG_FILE"
            continue
        fi

        # STEP 4e: Parse SFS and write Stairway Plot 2 blueprint
        #
        # easySFS dadi format:
        #   Line 1: "NPROJ+1 folded POPNAME"
        #   Line 2: SFS values index 0 (monomorphic) to NPROJ
        #   Line 3: mask
        # For Stairway Plot: extract indices 1 .. NPROJ/2 (folded half).
        #
        # nrand: 4 breakpoint counts capped at NPROJ/2 (avoids overfitting).
        # Standard values (7,15,22,28) used when HALF >= 28; otherwise
        # proportionally scaled to maintain strict monotonic increase.
        HALF=$((NPROJ / 2))

        SFS_STR=$(awk -v half="$HALF" 'NR==2 {
            for (i=2; i<=half+1; i++) printf "%s ", $i
        }' "$SFS_FILE" | sed 's/ $//')

        echo "  SFS parsed: ${HALF} bins (nseq=${NPROJ})" | tee -a "$POP_LOG"

        if [ "$HALF" -ge 28 ]; then
            NR1=7 ; NR2=15 ; NR3=22 ; NR4=28
        else
            NR1=$(( HALF / 4 ))
            [ "$NR1" -lt 2 ] && NR1=2
            NR2=$(( HALF / 2 ))
            NR3=$(( HALF * 3 / 4 ))
            NR4="$HALF"
            [ "$NR2" -le "$NR1" ] && NR2=$(( NR1 + 1 ))
            [ "$NR3" -le "$NR2" ] && NR3=$(( NR2 + 1 ))
            [ "$NR4" -le "$NR3" ] && NR4=$(( NR3 + 1 ))
        fi
        NRAND=$(printf "%d\t%d\t%d\t%d" "$NR1" "$NR2" "$NR3" "$NR4")

        BLUEPRINT="${OUT_POP}/${POP}_${TAG}.blueprint"
        PLOT_TITLE="${POP} (${TAG})"

        cat > "$BLUEPRINT" << BPEOF
popid: ${POP}
nseq: ${NPROJ}
L: ${L_POP[${POP}]}
whether_folded: true
SFS: ${SFS_STR}
pct_training: ${PCT_TRAIN}
nrand: ${NRAND}
project_dir: stairway
stairway_plot_dir: ${STAIRWAY_ES_DIR}
ninput: ${NINPUT}
mu: ${MU_REF}
year_per_generation: ${YEAR_PER_GEN}
plot_title: ${PLOT_TITLE}
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12
BPEOF

        echo "  Blueprint written: $BLUEPRINT" | tee -a "$POP_LOG"
        echo "  nrand: ${NR1} ${NR2} ${NR3} ${NR4}" | tee -a "$POP_LOG"

        # STEP 4f: Stairway Plot Phase 1 — Stairbuilder
        # Generates ${POP}_${TAG}.blueprint.sh (training + summary commands).
        # Must run from OUT_POP so project_dir (stairway/) is created there.
        STEP_SH="${OUT_POP}/${POP}_${TAG}.blueprint.sh"
        if [ -s "$STEP_SH" ]; then
            echo "  Stairbuilder output already exists, skipping" | tee -a "$POP_LOG"
        else
            echo "  Running Stairbuilder..." | tee -a "$POP_LOG"
            (
                cd "$OUT_POP"
                java -cp "$CLASSPATH" Stairbuilder "${POP}_${TAG}.blueprint" >> "$POP_LOG" 2>&1
            )
        fi

        if [ ! -s "$STEP_SH" ]; then
            echo "ERROR: Stairbuilder failed — blueprint.sh not generated" | tee -a "$POP_LOG" "$LOG_FILE"
            continue
        fi

        # STEP 4g: Stairway Plot Phase 2 — training replicates via sarray
        # sarray inherits CWD → must cd to OUT_POP before calling it.
        # sarray may exit non-zero for edge-case numerical failures (normal);
        # continue to summary regardless.
        TRAIN_SH="${OUT_POP}/step1_train_${POP}_${TAG}.sh"
        grep '^java.*Stairway_fold_training' "$STEP_SH" > "$TRAIN_SH"
        N_TRAIN=$(wc -l < "$TRAIN_SH")
        echo "  Submitting ${N_TRAIN} training jobs via sarray..." | tee -a "$POP_LOG"
        (
            cd "$OUT_POP"
            sarray --wait "$TRAIN_SH" || \
                echo "  WARN: sarray completed with some failed jobs" | tee -a "$POP_LOG"
            rm -f *.out
        )
        echo "  Training done for ${POP} ${TAG}" | tee -a "$POP_LOG"

        # STEP 4h: Stairway Plot Phase 3 — summary + built-in plot
        # xvfb-run provides a virtual X11 display required for PDF generation.
        SUMMARY_SH="${OUT_POP}/step2_summary_${POP}_${TAG}.sh"
        grep -v '^java.*Stairway_fold_training' "$STEP_SH" > "$SUMMARY_SH"
        echo "  Running summary..." | tee -a "$POP_LOG"
        (
            cd "$OUT_POP"
            xvfb-run bash "$SUMMARY_SH" >> "$POP_LOG" 2>&1
        ) || echo "  WARN: summary step failed for ${POP} ${TAG}" | tee -a "$POP_LOG"

        if [ -s "$SUMMARY_FILE" ]; then
            N_ROWS=$(wc -l < "$SUMMARY_FILE")
            echo "  Summary: ${SUMMARY_FILE} (${N_ROWS} time points)" | tee -a "$POP_LOG"
        else
            echo "WARN: summary not found: ${SUMMARY_FILE}" | tee -a "$POP_LOG" "$LOG_FILE"
        fi

        echo "=== DONE ${POP} ${TAG} $(date -Iseconds) ===" | tee -a "$POP_LOG" "$LOG_FILE"
    done

    # R visualization for this threshold
    echo "  R visualization — ${TAG}" | tee -a "$LOG_FILE"
    if [ -s "$R_PLOT_SCRIPT" ]; then
        Rscript --vanilla "$R_PLOT_SCRIPT" \
            "$OUT_BASE" "$TAG" "$MU_REF" "$YEAR_PER_GEN" \
            >> "$LOG_FILE" 2>&1 \
            || echo "WARN: R visualization failed for ${TAG}" | tee -a "$LOG_FILE"
    fi

done

# -------------------------------------------------------------------
# STEP 5: Cross-threshold R summary (all 4 thresholds combined)
# -------------------------------------------------------------------
echo "STEP 5: Cross-threshold R summary" | tee -a "$LOG_FILE"
if [ -s "$R_PLOT_SCRIPT" ]; then
    Rscript --vanilla "$R_PLOT_SCRIPT" \
        "$OUT_BASE" "ALL" "$MU_REF" "$YEAR_PER_GEN" \
        >> "$LOG_FILE" 2>&1 \
        || echo "WARN: cross-threshold R summary failed" | tee -a "$LOG_FILE"
fi

echo "DONE $(date -Iseconds) 07b.5 Stairway Plot K=${K_GLOBAL}" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_BASE}/"

#!/bin/bash
#SBATCH -p unlimitq
#SBATCH --time=96:00:00
#SBATCH -J 07.11-stairway_plot
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   sbatch 07.11-stairway_plot.slurm
#   Optional: K_GLOBAL=4 (default)
#
# Demographic inference via Stairway Plot 2 on the 4 ADMIXTURE clusters
# (K=4 partition, 155 French Guiana individuals only).
#
# Pipeline per cluster:
#   STEP 3a : Build per-cluster sample list + easySFS population map
#   STEP 3b : Extract per-cluster VCF (within-cluster polymorphic biallelic)
#   STEP 3c : easySFS preview (logged for reference)
#   STEP 3d : easySFS folded SFS (full 2n projection)
#   STEP 3e : Parse SFS → write Stairway Plot 2 blueprint
#   STEP 3f : Stairway Plot phase 1 (Stairbuilder → training script)
#   STEP 3g : Stairway Plot phase 2 (training replicates via sarray)
#   STEP 3h : Stairway Plot phase 3 (summary + built-in plot)
#
# L estimation (callable sites):
#   Reference genome: Dgu_HS1_HYBRID_SCAFFOLD.fa → 550,240,880 bp
#   TE-masked bases (Dgu_HS1_HYBRID_SCAFFOLD.fa.out.gff.gz): 296,419,896 bp
#   L = 550,240,880 - 296,419,896 = 253,820,984 callable sites
#
# Mutation rate: mu=5e-9 per site per generation (central estimate;
#   R script shows sensitivity for mu = 2..8 × 10^-9 by rescaling axes).
#
# Note on site independence: Stairway Plot 2 assumes independent sites.
#   The TE-filtered VCF contains ~5-10M SNPs per cluster without LD pruning.
#   LD decay in tropical trees is typically short-range (<1 kb), so the
#   impact on SFS shape is expected to be minor, but should be acknowledged
#   in the methods. Using all SNPs (vs LD-pruned) maximises SFS accuracy.
#
# Clusters run SEQUENTIALLY (memory safety for easySFS with ~5-10M SNPs).
# Training replicates run in PARALLEL via sarray (200 per population).
#
# Output: 07-population_structure/workflow_full/07.11-stairway_plot/
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
K_GLOBAL="${K_GLOBAL:-4}"

# L = non-TE callable genome size (bp)
L_CALLABLE=253820984
# Reference mutation rate (see R script for sensitivity to other values)
MU_REF="5e-9"
YEAR_PER_GEN=30
NINPUT=200     # Stairway Plot bootstrap replicates
PCT_TRAIN=0.67

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"
STAIRWAY_ES_DIR="/usr/local/bioinfo/src/stairway_plot/stairway_plot_v2.1.2/stairway_plot_es"
R_PLOT_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.11-stairway_plot.R"
OUT_BASE="${BASE}/07-population_structure/workflow_full/07.11-stairway_plot"

mkdir -p "${OUT_BASE}/logs"
LOG_FILE="${OUT_BASE}/logs/07.11-stairway_plot.log"
echo "START $(date -Iseconds) 07.11 Stairway Plot K=${K_GLOBAL}" > "$LOG_FILE"
echo "Input VCF  : $INPUT_VCF" | tee -a "$LOG_FILE"
echo "L callable : ${L_CALLABLE} bp" | tee -a "$LOG_FILE"
echo "mu_ref     : ${MU_REF}" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found" | tee -a "$LOG_FILE"; exit 1
fi
if [ ! -s "$CLUSTER_ASSIGNMENTS" ]; then
    echo "ERROR: cluster_assignments.tsv not found" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build list of 155 FG individuals (from cluster_assignments.tsv)
# -------------------------------------------------------------------
K_COL="K${K_GLOBAL}"
FG_SAMPLES="${OUT_BASE}/logs/fg_samples.txt"
awk 'NR > 1 { print $1 }' "$CLUSTER_ASSIGNMENTS" > "$FG_SAMPLES"
N_FG=$(wc -l < "$FG_SAMPLES")
echo "FG individuals: ${N_FG}" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Per-population pipeline (clusters run sequentially)
# -------------------------------------------------------------------
for POP in Pop_1 Pop_2 Pop_3 Pop_4; do
    OUT_POP="${OUT_BASE}/${POP}_K${K_GLOBAL}"
    mkdir -p "${OUT_POP}/easySFS" "${OUT_POP}/logs"

    POP_LOG="${OUT_POP}/logs/stairway_${POP}.log"
    echo "=== START ${POP} $(date -Iseconds) ===" | tee -a "$POP_LOG" "$LOG_FILE"

    # ----------------------------------------------------------------
    # STEP 3a: Per-cluster sample list and easySFS population map
    # ----------------------------------------------------------------
    CLUSTER_SAMPLES="${OUT_POP}/logs/samples_${POP}.txt"
    POP_MAP="${OUT_POP}/logs/pop_map_${POP}.txt"

    awk -v kcol="$K_COL" -v pop="$POP" '
    NR == 1 {
        for (i = 1; i <= NF; i++) if ($i == kcol) col = i
        if (!col) { print "ERROR: column " kcol " not found" > "/dev/stderr"; exit 2 }
        next
    }
    $col == pop { print $1 }
    ' "$CLUSTER_ASSIGNMENTS" > "$CLUSTER_SAMPLES"

    N_IND=$(wc -l < "$CLUSTER_SAMPLES")
    NPROJ=$((2 * N_IND))
    echo "  ${POP}: ${N_IND} individuals (projection: ${NPROJ})" | tee -a "$POP_LOG"

    if [ "$N_IND" -lt 4 ]; then
        echo "ERROR: fewer than 4 individuals in ${POP} — skipping" | tee -a "$POP_LOG" "$LOG_FILE"
        continue
    fi

    # easySFS population map: sample_id<TAB>population
    awk -v pop="$POP" '{ print $1 "\t" pop }' "$CLUSTER_SAMPLES" > "$POP_MAP"

    # ----------------------------------------------------------------
    # STEP 3b: Extract per-cluster VCF (within-cluster polymorphic)
    # Keeps biallelic SNPs polymorphic within the cluster (MAC >= 1:minor).
    # ----------------------------------------------------------------
    POP_VCF="${OUT_POP}/${POP}_K${K_GLOBAL}.vcf.gz"
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
        echo "  INFO: ${N_SNPS} within-cluster polymorphic SNPs for ${POP}" | tee -a "$POP_LOG"
    fi

    # ----------------------------------------------------------------
    # STEP 3c: easySFS preview (logged for reference)
    # ----------------------------------------------------------------
    PREVIEW_LOG="${OUT_POP}/logs/easySFS_preview_${POP}.log"
    echo "  Running easySFS preview..." | tee -a "$POP_LOG"
    easySFS.py \
        -i "$POP_VCF" \
        -p "$POP_MAP" \
        -a \
        --ploidy 2 \
        --preview \
        > "$PREVIEW_LOG" 2>&1 || true
    cat "$PREVIEW_LOG" | tee -a "$POP_LOG"

    # ----------------------------------------------------------------
    # STEP 3d: easySFS folded SFS (full 2n projection)
    # Folded (minor allele) SFS is the default in easySFS (no --unfolded flag).
    # -a: keep all SNPs (appropriate for WGS data, not RAD loci).
    # ----------------------------------------------------------------
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

    # ----------------------------------------------------------------
    # STEP 3e: Parse SFS and write Stairway Plot 2 blueprint
    #
    # easySFS dadi format:
    #   Line 1: "NPROJ+1 folded POPNAME"
    #   Line 2: space-separated SFS values (index 0=monomorphic to NPROJ)
    #   Line 3: mask (1=masked, 0=unmasked)
    #
    # For Stairway Plot (folded, nseq=NPROJ):
    #   SFS = values at indices 1 .. NPROJ/2 (skip index 0 = monomorphic;
    #   indices NPROJ/2+1 .. NPROJ are 0 due to folding).
    #
    # nrand: breakpoint complexities to try (tab-separated).
    #   Capped at NPROJ/2 to avoid overfitting for small populations.
    # ----------------------------------------------------------------
    HALF=$((NPROJ / 2))

    SFS_STR=$(awk -v half="$HALF" 'NR==2 {
        for (i=2; i<=half+1; i++) printf "%s ", $i
    }' "$SFS_FILE" | sed 's/ $//')

    echo "  SFS parsed: ${HALF} bins for ${POP} (nseq=${NPROJ})" | tee -a "$POP_LOG"

    # nrand: 4 breakpoint counts, evenly spaced across [floor(HALF/4), HALF].
    # For large populations (HALF >= 28) use the standard values 7,15,22,28.
    # For smaller populations (Pop_3 HALF=27, Pop_4 HALF=16) values are scaled
    # proportionally to HALF so they remain ≤ HALF (avoids overfitting) and
    # are strictly monotonically increasing (required by Stairway Plot 2).
    if [ "$HALF" -ge 28 ]; then
        NR1=7 ; NR2=15 ; NR3=22 ; NR4=28
    else
        NR1=$(( HALF / 4 ))
        [ "$NR1" -lt 2 ] && NR1=2
        NR2=$(( HALF / 2 ))
        NR3=$(( HALF * 3 / 4 ))
        NR4="$HALF"
        # Ensure strict monotonicity after integer division
        [ "$NR2" -le "$NR1" ] && NR2=$(( NR1 + 1 ))
        [ "$NR3" -le "$NR2" ] && NR3=$(( NR2 + 1 ))
        [ "$NR4" -le "$NR3" ] && NR4=$(( NR3 + 1 ))
    fi
    NRAND=$(printf "%d\t%d\t%d\t%d" "$NR1" "$NR2" "$NR3" "$NR4")

    BLUEPRINT="${OUT_POP}/${POP}_K${K_GLOBAL}.blueprint"

    cat > "$BLUEPRINT" << BPEOF
popid: ${POP}
nseq: ${NPROJ}
L: ${L_CALLABLE}
whether_folded: true
SFS: ${SFS_STR}
pct_training: ${PCT_TRAIN}
nrand: ${NRAND}
project_dir: stairway
stairway_plot_dir: ${STAIRWAY_ES_DIR}
ninput: ${NINPUT}
mu: ${MU_REF}
year_per_generation: ${YEAR_PER_GEN}
plot_title: ${POP} (K=${K_GLOBAL})
xrange: 0.1,10000
yrange: 0,0
xspacing: 2
yspacing: 2
fontsize: 12
BPEOF

    echo "  Blueprint written: $BLUEPRINT" | tee -a "$POP_LOG"
    echo "  nrand: ${NR1} ${NR2} ${NR3} ${NR4}" | tee -a "$POP_LOG"

    # ----------------------------------------------------------------
    # STEP 3f: Stairway Plot Phase 1 — Stairbuilder
    # Generates ${POP}_K${K_GLOBAL}.blueprint.sh (training + summary)
    # Stairbuilder is run from OUT_POP so that project_dir (stairway/)
    # is created relative to OUT_POP.
    # ----------------------------------------------------------------
    STEP_SH="${OUT_POP}/${POP}_K${K_GLOBAL}.blueprint.sh"

    if [ -s "$STEP_SH" ]; then
        echo "  Stairbuilder output already exists, skipping" | tee -a "$POP_LOG"
    else
        echo "  Running Stairbuilder..." | tee -a "$POP_LOG"
        (
            cd "$OUT_POP"
            java -cp "$CLASSPATH" Stairbuilder "${POP}_K${K_GLOBAL}.blueprint" >> "$POP_LOG" 2>&1
        )
    fi

    if [ ! -s "$STEP_SH" ]; then
        echo "ERROR: Stairbuilder failed — blueprint.sh not generated" | tee -a "$POP_LOG" "$LOG_FILE"
        continue
    fi

    # ----------------------------------------------------------------
    # STEP 3g: Stairway Plot Phase 2 — training replicates via sarray
    # Each replicate is one java Stairway_fold_training_testing7 call.
    # sarray parallelises all replicates as separate SLURM jobs.
    # sarray inherits the CWD → must cd to OUT_POP before calling it.
    # ----------------------------------------------------------------
    SUMMARY_FILE="${OUT_POP}/stairway/${POP}.final.summary"
    if [ -s "$SUMMARY_FILE" ]; then
        echo "  Summary already exists, skipping training" | tee -a "$POP_LOG"
    else
        TRAIN_SH="${OUT_POP}/step1_train_${POP}.sh"
        grep '^java.*Stairway_fold_training' "$STEP_SH" > "$TRAIN_SH"
        N_TRAIN=$(wc -l < "$TRAIN_SH")
        echo "  Submitting ${N_TRAIN} training jobs via sarray..." | tee -a "$POP_LOG"
        (
            cd "$OUT_POP"
            # sarray may exit non-zero if a minority of replicates fail (normal
            # for numerical edge cases); continue to summary step regardless.
            sarray --wait "$TRAIN_SH" || \
                echo "  WARN: sarray completed with some failed jobs for ${POP}" | tee -a "$POP_LOG"
            # Remove per-replicate SLURM output files created by sarray
            rm -f *.out
        )
        echo "  Training done for ${POP}" | tee -a "$POP_LOG"

        # ------------------------------------------------------------
        # STEP 3h: Stairway Plot Phase 3 — summary + built-in plot
        # All non-training commands from blueprint.sh (summary + R plot).
        # xvfb-run provides a virtual X11 display for PDF generation.
        # ------------------------------------------------------------
        SUMMARY_SH="${OUT_POP}/step2_summary_${POP}.sh"
        grep -v '^java.*Stairway_fold_training' "$STEP_SH" > "$SUMMARY_SH"
        echo "  Running summary..." | tee -a "$POP_LOG"
        (
            cd "$OUT_POP"
            xvfb-run bash "$SUMMARY_SH" >> "$POP_LOG" 2>&1
        ) || echo "  WARN: summary step failed for ${POP} (check POP_LOG)" | tee -a "$POP_LOG"
    fi

    if [ -s "$SUMMARY_FILE" ]; then
        N_ROWS=$(wc -l < "$SUMMARY_FILE")
        echo "  Summary: ${SUMMARY_FILE} (${N_ROWS} time points)" | tee -a "$POP_LOG"
    else
        echo "WARN: summary not found: ${SUMMARY_FILE}" | tee -a "$POP_LOG" "$LOG_FILE"
    fi

    echo "=== DONE ${POP} $(date -Iseconds) ===" | tee -a "$POP_LOG" "$LOG_FILE"
done

# -------------------------------------------------------------------
# STEP 4: R visualization (combined Ne(t) + mu sensitivity)
# -------------------------------------------------------------------
echo "STEP 4: R visualization" | tee -a "$LOG_FILE"
if [ -s "$R_PLOT_SCRIPT" ]; then
    Rscript --vanilla "$R_PLOT_SCRIPT" \
        "$OUT_BASE" \
        "$K_GLOBAL" \
        "$MU_REF" \
        "$YEAR_PER_GEN" \
        >> "$LOG_FILE" 2>&1 \
        || echo "WARN: R visualization failed" | tee -a "$LOG_FILE"
else
    echo "WARN: R script not found: $R_PLOT_SCRIPT" | tee -a "$LOG_FILE"
fi

echo "DONE $(date -Iseconds) 07.11 Stairway Plot" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_BASE}/"

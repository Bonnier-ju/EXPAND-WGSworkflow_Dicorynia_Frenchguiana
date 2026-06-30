#!/bin/bash
#SBATCH -p workq
#SBATCH --time=48:00:00
#SBATCH -J 07.8-Ne_cluster
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   K_GLOBAL=2 sbatch 07.8-Ne_cluster.slurm
#   K_GLOBAL=3 sbatch 07.8-Ne_cluster.slurm
#   K_GLOBAL=4 sbatch 07.8-Ne_cluster.slurm
#
# Parameters:
#   K_GLOBAL  — global ADMIXTURE K defining the populations (2, 3 or 4)
#   N_CHR     — number of chromosome-level scaffolds used (default 13)
#
# Estimates effective population size (Ne) per genetic cluster:
#   - currentNe  v1: current Ne with 50%/90% CI (neural-network LD)
#   - currentNe2 v2: current Ne, panmictic + structure-corrected (-x)
#   - GONE   v1    : temporal Ne trajectory over ~2000 generations
#   - GONE2  v2    : temporal Ne trajectory, panmictic + structure (-x)
#
# Both analyses use only the 13 chromosome-level Super-Scaffolds (≥30 Mb)
# renamed to integers 1-13 (see metadata/scaffold_info.txt).
# Input: TE-filtered SNP VCF (06.5, NOT LD-pruned — LD signal required).
#
# Output: 07-population_structure/workflow_full/07.5-Ne_estimation/07.5.K-K=${K_GLOBAL}/
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/1.90b7
module load bioinfo/currentNe/37daed5
module load bioinfo/GONE/134996d
module load bioinfo/GONE2/1.0.2
module load bioinfo/currentNe2/1.0.1
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
K_GLOBAL="${K_GLOBAL:-2}"
N_CHR=13
CPUS="${SLURM_CPUS_PER_TASK:-8}"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.8-Ne_cluster.R"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.5-Ne_estimation/07.5.${K_GLOBAL}-K=${K_GLOBAL}"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
POP_DIR="${OUT_DIR}/per_cluster"
MAP_FILE="${OUT_DIR}/sample_cluster_map.tsv"
POPS_LIST="${TMP_DIR}/clusters.list"
LOG_FILE="${LOG_DIR}/07.8-Ne_K${K_GLOBAL}.log"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR" "$POP_DIR"
echo "START $(date -Iseconds) 07.8 Ne estimation by cluster K=${K_GLOBAL}" > "$LOG_FILE"
echo "Input VCF          : $INPUT_VCF"          | tee -a "$LOG_FILE"
echo "Cluster assignments: $CLUSTER_ASSIGNMENTS" | tee -a "$LOG_FILE"
echo "CPUs               : $CPUS"               | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"; exit 1
fi
if [ ! -s "$CLUSTER_ASSIGNMENTS" ]; then
    echo "ERROR: cluster_assignments.tsv not found — run 07.0-prep_cluster_assignments.R first" \
        | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build sample -> cluster map from cluster_assignments.tsv
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

while IFS= read -r pop; do
    # Two-column format (FID IID) required by PLINK --keep
    awk -v p="$pop" '$2 == p { print $1 "\t" $1 }' "$MAP_FILE" \
        > "${POP_DIR}/${pop}.samples.txt"
    N=$(wc -l < "${POP_DIR}/${pop}.samples.txt")
    echo "  ${pop}: ${N} individuals" | tee -a "$LOG_FILE"
done < "$POPS_LIST"

# -------------------------------------------------------------------
# STEP 3: Shared VCF preparation — filter to 13 chr-level scaffolds
#         and rename Super-Scaffold_X -> integer (done once, reused)
# -------------------------------------------------------------------
echo "STEP 3: preparing chromosome-filtered VCF (13 scaffolds)" | tee -a "$LOG_FILE"

CHR_REGIONS="Super-Scaffold_1,Super-Scaffold_2,Super-Scaffold_4,Super-Scaffold_5,\
Super-Scaffold_6,Super-Scaffold_7,Super-Scaffold_8,Super-Scaffold_12,\
Super-Scaffold_13,Super-Scaffold_14,Super-Scaffold_16,Super-Scaffold_17,\
Super-Scaffold_21"

# Rename file: Super-Scaffold_X -> integer (tab-separated, bcftools format)
CHR_RENAME="${TMP_DIR}/chr_rename.txt"
cat > "$CHR_RENAME" << 'EOF'
Super-Scaffold_1	1
Super-Scaffold_2	2
Super-Scaffold_4	3
Super-Scaffold_5	4
Super-Scaffold_6	5
Super-Scaffold_7	6
Super-Scaffold_8	7
Super-Scaffold_12	8
Super-Scaffold_13	9
Super-Scaffold_14	10
Super-Scaffold_16	11
Super-Scaffold_17	12
Super-Scaffold_21	13
EOF

VCF_CHR13="${TMP_DIR}/vcf_chr13.vcf.gz"
if [ ! -s "$VCF_CHR13" ]; then
    echo "INFO: filtering + renaming chromosomes" | tee -a "$LOG_FILE"
    bcftools view -r "$CHR_REGIONS" "$INPUT_VCF" \
        | bcftools annotate --rename-chrs "$CHR_RENAME" \
        | bgzip -c > "$VCF_CHR13"
    bcftools index -t "$VCF_CHR13" >> "$LOG_FILE" 2>&1
    N_SNPS=$(bcftools view -H "$VCF_CHR13" | wc -l)
    echo "INFO: ${N_SNPS} SNPs retained on ${N_CHR} chromosomes" | tee -a "$LOG_FILE"
else
    echo "INFO: reusing existing $VCF_CHR13" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 4: Per-cluster analyses
# -------------------------------------------------------------------
echo "STEP 4: per-cluster PLINK conversion + currentNe + GONE" | tee -a "$LOG_FILE"

while IFS= read -r pop; do
    keep_file="${POP_DIR}/${pop}.samples.txt"
    n_ind=$(wc -l < "$keep_file")
    plink_dir="${POP_DIR}/${pop}/plink"
    cne_dir="${POP_DIR}/${pop}/currentNe"
    gone_dir="${POP_DIR}/${pop}/GONE"
    mkdir -p "$plink_dir" "$cne_dir" "$gone_dir"

    echo "--- ${pop} (n=${n_ind}) ---" | tee -a "$LOG_FILE"

    # ------------------------------------------------------------------
    # 4a: PLINK — VCF -> BED
    #     then BED -> MAP/PED for GONE (full SNP set; GONE subsamples
    #     to maxNSNP=2000/chr internally)
    # ------------------------------------------------------------------
    # Sample size warning: currentNe and GONE require >= 25 individuals for reliable estimates
    if [ "$n_ind" -lt 25 ]; then
        echo "  WARN: [${pop}] n=${n_ind} < 25 — currentNe and GONE estimates will be unstable (wide CI)" \
            | tee -a "$LOG_FILE"
    fi

    # ------------------------------------------------------------------
    # 4a: PLINK — VCF -> BED (skip if BED already exists)
    # ------------------------------------------------------------------
    if [ -s "${plink_dir}/${pop}.bed" ]; then
        N_SNP=$(wc -l < "${plink_dir}/${pop}.bim")
        echo "  [${pop}] reusing existing BED: ${n_ind} individuals, ${N_SNP} SNPs" | tee -a "$LOG_FILE"
    else
        echo "  [${pop}] PLINK conversion (BED + full MAP/PED)" | tee -a "$LOG_FILE"
        plink --vcf "$VCF_CHR13" \
              --double-id \
              --keep "$keep_file" \
              --chr-set ${N_CHR} \
              --maf 0.05 \
              --geno 0.1 \
              --make-bed \
              --out "${plink_dir}/${pop}" \
              >> "$LOG_FILE" 2>&1
        N_SNP=$(wc -l < "${plink_dir}/${pop}.bim")
        echo "  [${pop}] BED: ${n_ind} individuals, ${N_SNP} SNPs (MAF>=0.05, geno<=10%)" | tee -a "$LOG_FILE"

        plink --bfile "${plink_dir}/${pop}" \
              --chr-set ${N_CHR} \
              --recode \
              --out "${plink_dir}/${pop}" \
              >> "$LOG_FILE" 2>&1

        # Subsampled MAP/PED for currentNe (hard limit: MAXLOCI=2,000,000)
        # ~16.4M SNPs x 0.06 ≈ 980k SNPs — sufficient and well within the limit
        echo "  [${pop}] PLINK thinning for currentNe (~1M SNPs)" | tee -a "$LOG_FILE"
        plink --bfile "${plink_dir}/${pop}" \
              --chr-set ${N_CHR} \
              --thin 0.06 \
              --recode \
              --out "${plink_dir}/${pop}_cne" \
              >> "$LOG_FILE" 2>&1
        N_SNP_CNE=$(wc -l < "${plink_dir}/${pop}_cne.map")
        echo "  [${pop}] subsampled: ${N_SNP_CNE} SNPs for currentNe" | tee -a "$LOG_FILE"
    fi

    # ------------------------------------------------------------------
    # 4b: currentNe (skip if output already exists)
    #     Must run from plink_dir so the MAP file is co-located with PED
    #     -k 0 : each offspring from a new random pairing
    #     -t    : number of threads (= SLURM cpus-per-task)
    # ------------------------------------------------------------------
    if [ -s "${cne_dir}/${pop}_currentNe" ]; then
        echo "  [${pop}] reusing existing currentNe output" | tee -a "$LOG_FILE"
    else
        echo "  [${pop}] currentNe" | tee -a "$LOG_FILE"
        (
            cd "$plink_dir"
            currentNe \
                -t "${CPUS}" \
                -k 0 \
                -o "${cne_dir}/${pop}_currentNe" \
                "${pop}_cne.ped" ${N_CHR} >> "$LOG_FILE" 2>&1
        )
        echo "  [${pop}] currentNe done" | tee -a "$LOG_FILE"
    fi

    # ------------------------------------------------------------------
    # 4c: GONE (skip if output already exists)
    #     Must run from gone_dir (working directory for temp files +
    #     output). Symlink MAP/PED so GONE can find them by basename.
    #     INPUT_PARAMETERS_FILE is generated with threads=${CPUS}.
    #     maxNSNP=2000: beyond 2000 SNPs/chr no gain in precision.
    # ------------------------------------------------------------------
    if [ -s "${gone_dir}/Output_Ne_${pop}" ]; then
        echo "  [${pop}] reusing existing GONE output" | tee -a "$LOG_FILE"
    else
    echo "  [${pop}] GONE" | tee -a "$LOG_FILE"

    ln -sf "${plink_dir}/${pop}.map" "${gone_dir}/${pop}.map"
    ln -sf "${plink_dir}/${pop}.ped" "${gone_dir}/${pop}.ped"

    cat > "${gone_dir}/INPUT_PARAMETERS_FILE" << EOF
PHASE=2
cMMb=1
DIST=1
NGEN=2000
NBIN=400
MAF=0.05
ZERO=1
maxNCHROM=-99
maxNSNP=2000
hc=0.05
REPS=40
threads=${CPUS}
EOF

    (
        cd "$gone_dir"
        /usr/local/bioinfo/src/GONE/GONE-134996d/Linux/script_GONE.sh "${pop}" >> "$LOG_FILE" 2>&1
    ) || echo "WARN: GONE returned non-zero for ${pop} — check log" | tee -a "$LOG_FILE"

    echo "  [${pop}] GONE done" | tee -a "$LOG_FILE"
    fi  # end GONE checkpoint

    # ------------------------------------------------------------------
    # 4d: Pre-thin PED for GONE2 (~1.5M SNPs)
    #     MAXLOCI=2,000,000 is a compile-time array limit: both tools read
    #     the entire PED sequentially. Files with 5-6M SNPs crash on the
    #     2M boundary. Fix: pre-thin to ~1.5M SNPs with PLINK --thin so
    #     both tools can read the file (1.5M < MAXLOCI).
    #     This file is used AS-IS by GONE2 (all 1.5M SNPs — its binned
    #     algorithm is fast regardless of SNP count).
    #     currentNe2 receives the same file but subsamples internally to
    #     200k SNPs via -s 200000 (~1.6h/run vs ~89h with 1.5M effective).
    #     The thinned MAP/PED are symlinked as ${pop}.ped in GONE2 subdirs
    #     so the output name (${pop}.ped_GONE2_Ne) is preserved for the R
    #     parser.
    # ------------------------------------------------------------------
    NE2_PED="${plink_dir}/${pop}_ne2.ped"
    NE2_MAP="${plink_dir}/${pop}_ne2.map"

    if [ -s "$NE2_PED" ]; then
        N_SNP_NE2=$(wc -l < "$NE2_MAP")
        echo "  [${pop}] reusing existing ne2 PED: ${N_SNP_NE2} SNPs" | tee -a "$LOG_FILE"
    else
        N_SNP_FULL=$(wc -l < "${plink_dir}/${pop}.bim")
        THIN_NE2=$(awk -v n="$N_SNP_FULL" 'BEGIN {t=1500000/n; print (t>0.95?0.95:t)}')
        echo "  [${pop}] thinning: --thin ${THIN_NE2} from ${N_SNP_FULL} SNPs (target 1.5M for GONE2)" | tee -a "$LOG_FILE"
        plink \
            --bfile "${plink_dir}/${pop}" \
            --chr-set ${N_CHR} \
            --thin "$THIN_NE2" \
            --recode \
            --out "${plink_dir}/${pop}_ne2" \
            >> "$LOG_FILE" 2>&1
        N_SNP_NE2=$(wc -l < "$NE2_MAP")
        echo "  [${pop}] ne2 PED: ${N_SNP_NE2} SNPs retained" | tee -a "$LOG_FILE"
    fi

    # ------------------------------------------------------------------
    # 4e: currentNe2 — default (panmictic) + structure (-x)
    #     Input: ne2-thinned PED (~1.5M SNPs, safely below MAXLOCI=2M).
    #     -s 200000: internal subsampling to 200k SNPs for computation;
    #       200k gives ~20B pairs → ~1.6h/run at 8 threads (vs ~89h at
    #       1.5M effective). Precision for Ne estimation is excellent at
    #       200k (>10× the baddoch validation example of 157k SNPs).
    #     -k 0 : random pairing (default only; incompatible with -x).
    # ------------------------------------------------------------------
    cne2_dir="${POP_DIR}/${pop}/currentNe2"
    mkdir -p "$cne2_dir"

    CNE2_DEFAULT="${cne2_dir}/${pop}_currentNe2_default.txt"
    CNE2_STRUCT="${cne2_dir}/${pop}_currentNe2_mix.txt"

    if [ -s "$CNE2_DEFAULT" ]; then
        echo "  [${pop}] reusing existing currentNe2 default output" | tee -a "$LOG_FILE"
    else
        echo "  [${pop}] currentNe2 default" | tee -a "$LOG_FILE"
        (
            cd "$plink_dir"
            currentne2 \
                -s 200000 \
                -r 1 \
                -t "${CPUS}" \
                -k 0 \
                -o "$CNE2_DEFAULT" \
                "${pop}_ne2.ped" >> "$LOG_FILE" 2>&1
        ) || echo "WARN: currentNe2 default returned non-zero for ${pop} — check log" | tee -a "$LOG_FILE"
        echo "  [${pop}] currentNe2 default done" | tee -a "$LOG_FILE"
    fi

    if [ -s "$CNE2_STRUCT" ]; then
        echo "  [${pop}] reusing existing currentNe2 -x output" | tee -a "$LOG_FILE"
    else
        echo "  [${pop}] currentNe2 -x (structure correction)" | tee -a "$LOG_FILE"
        (
            cd "$plink_dir"
            currentne2 \
                -x \
                -s 200000 \
                -r 1 \
                -t "${CPUS}" \
                -o "$CNE2_STRUCT" \
                "${pop}_ne2.ped" >> "$LOG_FILE" 2>&1
        ) || echo "WARN: currentNe2 -x returned non-zero for ${pop} — check log" | tee -a "$LOG_FILE"
        echo "  [${pop}] currentNe2 -x done" | tee -a "$LOG_FILE"
    fi

    # ------------------------------------------------------------------
    # 4f: GONE2 — default (panmictic) + structure (-x)
    #     Input: ne2-thinned PED (~1.5M SNPs), symlinked as ${pop}.ped in
    #     each subdir so output (${pop}.ped_GONE2_Ne) matches R parser.
    #     Run in dedicated subdirectories: both options produce the same
    #     filename — separate dirs avoid overwrite.
    # ------------------------------------------------------------------
    gone2_dir="${POP_DIR}/${pop}/GONE_2"
    g2_def_dir="${gone2_dir}/default"
    g2_str_dir="${gone2_dir}/structure"
    mkdir -p "$g2_def_dir" "$g2_str_dir"

    GONE2_DEFAULT_NE="${g2_def_dir}/${pop}.ped_GONE2_Ne"
    GONE2_STRUCT_NE="${g2_str_dir}/${pop}.ped_GONE2_Ne"

    if [ -s "$GONE2_DEFAULT_NE" ]; then
        echo "  [${pop}] reusing existing GONE2 default output" | tee -a "$LOG_FILE"
    else
        echo "  [${pop}] GONE2 default" | tee -a "$LOG_FILE"
        ln -sf "$NE2_MAP" "${g2_def_dir}/${pop}.map"
        ln -sf "$NE2_PED" "${g2_def_dir}/${pop}.ped"
        (
            cd "$g2_def_dir"
            gone2 -r 1 -t "${CPUS}" "${pop}.ped" >> "$LOG_FILE" 2>&1
        ) || echo "WARN: GONE2 default returned non-zero for ${pop} — check log" | tee -a "$LOG_FILE"
        echo "  [${pop}] GONE2 default done" | tee -a "$LOG_FILE"
    fi

    if [ -s "$GONE2_STRUCT_NE" ]; then
        echo "  [${pop}] reusing existing GONE2 -x output" | tee -a "$LOG_FILE"
    else
        echo "  [${pop}] GONE2 -x (structure correction)" | tee -a "$LOG_FILE"
        ln -sf "$NE2_MAP" "${g2_str_dir}/${pop}.map"
        ln -sf "$NE2_PED" "${g2_str_dir}/${pop}.ped"
        (
            cd "$g2_str_dir"
            gone2 -x -r 1 -t "${CPUS}" "${pop}.ped" >> "$LOG_FILE" 2>&1
        ) || echo "WARN: GONE2 -x returned non-zero for ${pop} — check log" | tee -a "$LOG_FILE"
        echo "  [${pop}] GONE2 -x done" | tee -a "$LOG_FILE"
    fi

done < "$POPS_LIST"

# -------------------------------------------------------------------
# STEP 5: R — compile results and produce plots
# -------------------------------------------------------------------
echo "STEP 5: R compilation and plots" | tee -a "$LOG_FILE"

Rscript --vanilla "$R_SCRIPT" \
    "$OUT_DIR" \
    "$POPS_LIST" \
    "$K_GLOBAL" >> "$LOG_FILE" 2>&1 \
    || echo "WARN: R script returned non-zero — check log" | tee -a "$LOG_FILE"

echo "DONE $(date -Iseconds) 07.8 Ne K=${K_GLOBAL} completed" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}/"

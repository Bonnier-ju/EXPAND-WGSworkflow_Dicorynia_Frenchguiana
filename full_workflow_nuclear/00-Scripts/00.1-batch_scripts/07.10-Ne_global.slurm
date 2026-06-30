#!/bin/bash
#SBATCH -p workq
#SBATCH --time=48:00:00
#SBATCH -J 07.10-Ne_global
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# USAGE
#   sbatch 07.10-Ne_global.slurm
#
# No parameters required. Runs Ne analyses on the pooled set of all
# 155 French Guiana individuals (outgroups absent from
# cluster_assignments.tsv are automatically excluded).
#
# Tools:
#   - currentNe2 (default)  : current Ne, panmictic assumption, 50%/90% CI
#   - currentNe2 (-x)       : metapopulation Ne + FST + migration rate +
#                             inferred number of subpopulations (s)
#   - GONE2 (default)       : temporal Ne trajectory, panmictic assumption
#   - GONE2 (-x)            : structure-corrected temporal Ne trajectory
#
# Scientific rationale:
#   Running both tools on the unsubdivided dataset provides an LD-based
#   inference of population structure entirely independent of ADMIXTURE.
#   The estimated s (number of subpopulations) from currentNe2 -x can be
#   compared to the ADMIXTURE K_best to cross-validate the cluster signal.
#   The comparison of GONE2 default vs -x trajectories reveals whether the
#   global signal is driven by structure or genuine demographic change.
#
# Parameters:
#   Input PED  : pre-thinned to ~1.5M SNPs (< MAXLOCI=2M compile limit).
#                Tools read the entire file sequentially — 5-6M SNP files
#                crash on the 2M boundary even when -s requests fewer SNPs.
#   -s 200000  : currentNe2 only — internal subsampling to 200k effective
#                SNPs (~20B pairs, ~1.6h/run). Not used for GONE2 (binned
#                algorithm; fast regardless of SNP count).
#   -r 1       : uniform 1 cM/Mb recombination rate (no genetic map for
#                Dicorynia guianensis; same assumption as GONE v1 cMMb=1)
#   -t 8       : threads (= SLURM cpus-per-task)
#   -k 0       : random-mating assumption (currentNe2 default only;
#                incompatible with -x per Santiago et al. 2025 README)
#
# Output: 07-population_structure/workflow_full/07.10-Ne_global/
# -------------------------------------------------------------------

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/PLINK/1.90b7
module load bioinfo/currentNe2/1.0.1
module load bioinfo/GONE2/1.0.2
module load statistics/R/4.3.0
export R_LIBS_USER="${HOME}/R/R"

# -------------------------------------------------------------------
# STEP 1: Parameters and paths
# -------------------------------------------------------------------
N_CHR=13
CPUS="${SLURM_CPUS_PER_TASK:-8}"
LABEL="guiana_all"

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
INPUT_VCF="${BASE}/06-variant_filtering/workflow_full/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
CLUSTER_ASSIGNMENTS="${BASE}/metadata/cluster_assignments.tsv"
R_SCRIPT="${BASE}/00-Scripts/00.2-R_scripts/07.10-Ne_global.R"

OUT_DIR="${BASE}/07-population_structure/workflow_full/07.10-Ne_global"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
PLINK_DIR="${OUT_DIR}/plink"
CNE2_DIR="${OUT_DIR}/currentNe2"
GONE2_DIR="${OUT_DIR}/GONE2"
LOG_FILE="${LOG_DIR}/07.10-Ne_global.log"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR" "$PLINK_DIR" "$CNE2_DIR" \
         "${GONE2_DIR}/default" "${GONE2_DIR}/structure"

echo "START $(date -Iseconds) 07.10 Ne global — all 155 French Guiana individuals" > "$LOG_FILE"
echo "Input VCF          : $INPUT_VCF"          | tee -a "$LOG_FILE"
echo "Cluster assignments: $CLUSTER_ASSIGNMENTS" | tee -a "$LOG_FILE"
echo "Label              : $LABEL"               | tee -a "$LOG_FILE"
echo "CPUs               : $CPUS"               | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found: $INPUT_VCF" | tee -a "$LOG_FILE"; exit 1
fi
if [ ! -s "$CLUSTER_ASSIGNMENTS" ]; then
    echo "ERROR: cluster_assignments.tsv not found" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build keep file — all French Guiana individuals
# cluster_assignments.tsv contains only the 155 FG individuals;
# the 4 Cameroun/Benin outgroups are absent from this file.
# PLINK --keep expects: FID IID (tab-separated, --double-id convention).
# -------------------------------------------------------------------
echo "STEP 2: building French Guiana keep file" | tee -a "$LOG_FILE"

KEEP_FILE="${TMP_DIR}/${LABEL}.keep.txt"
awk 'NR > 1 && NF >= 1 { print $1 "\t" $1 }' "$CLUSTER_ASSIGNMENTS" > "$KEEP_FILE"

N_IND=$(wc -l < "$KEEP_FILE")
echo "INFO: ${N_IND} French Guiana individuals in keep file" | tee -a "$LOG_FILE"

if [ "$N_IND" -lt 10 ]; then
    echo "ERROR: fewer than 10 individuals — check cluster_assignments.tsv" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 3: Chromosome-filtered VCF
# Filter to the 13 chromosome-level Super-Scaffolds (>= 30 Mb) and
# rename Super-Scaffold_X -> integers 1-13 for PLINK compatibility.
# -------------------------------------------------------------------
echo "STEP 3: preparing chromosome-filtered VCF (13 scaffolds)" | tee -a "$LOG_FILE"

CHR_REGIONS="Super-Scaffold_1,Super-Scaffold_2,Super-Scaffold_4,Super-Scaffold_5,\
Super-Scaffold_6,Super-Scaffold_7,Super-Scaffold_8,Super-Scaffold_12,\
Super-Scaffold_13,Super-Scaffold_14,Super-Scaffold_16,Super-Scaffold_17,\
Super-Scaffold_21"

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
    echo "INFO: ${N_SNPS} SNPs on ${N_CHR} chromosomes" | tee -a "$LOG_FILE"
else
    echo "INFO: reusing existing $VCF_CHR13" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 4: PLINK — VCF -> BED, then BED -> PED/MAP
# MAF >= 0.05 computed across all 155 individuals.
# Both formats are produced: BED for QC stats, PED/MAP for Ne tools.
# -------------------------------------------------------------------
echo "STEP 4: PLINK VCF -> BED + PED (${N_IND} individuals)" | tee -a "$LOG_FILE"

if [ -s "${PLINK_DIR}/${LABEL}.bed" ]; then
    N_SNP=$(wc -l < "${PLINK_DIR}/${LABEL}.bim")
    echo "INFO: reusing existing BED: ${N_IND} individuals, ${N_SNP} SNPs" | tee -a "$LOG_FILE"
else
    plink \
        --vcf "$VCF_CHR13" \
        --double-id \
        --keep "$KEEP_FILE" \
        --chr-set ${N_CHR} \
        --maf 0.05 \
        --geno 0.1 \
        --make-bed \
        --out "${PLINK_DIR}/${LABEL}" \
        >> "$LOG_FILE" 2>&1
    N_SNP=$(wc -l < "${PLINK_DIR}/${LABEL}.bim")
    echo "INFO: BED — ${N_IND} individuals, ${N_SNP} SNPs (MAF>=0.05, geno<=10%)" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 4b: Pre-thin PED for currentNe2 and GONE2
# MAXLOCI=2,000,000 is a compile-time array limit: both tools read
# the entire PED sequentially before any internal subsampling, so
# files with 5-6M SNPs crash on the 2M boundary. Fix: pre-thin to
# ~1.5M SNPs via PLINK --thin (25% margin below 2M). The thin ratio
# is computed dynamically from the BIM line count.
# The thinned PED is symlinked as ${LABEL}.ped in each GONE2 subdir
# so the output name (${LABEL}.ped_GONE2_Ne) is preserved for the R
# parser.
# -------------------------------------------------------------------
NE2_PED="${PLINK_DIR}/${LABEL}_ne2.ped"
NE2_MAP="${PLINK_DIR}/${LABEL}_ne2.map"

if [ -s "$NE2_PED" ]; then
    N_SNP_NE2=$(wc -l < "$NE2_MAP")
    echo "INFO: reusing existing ne2 PED: ${N_SNP_NE2} SNPs" | tee -a "$LOG_FILE"
else
    N_SNP_FULL=$(wc -l < "${PLINK_DIR}/${LABEL}.bim")
    THIN_NE2=$(awk -v n="$N_SNP_FULL" 'BEGIN {t=1500000/n; print (t>0.95?0.95:t)}')
    echo "INFO: thinning for currentNe2/GONE2: --thin ${THIN_NE2} from ${N_SNP_FULL} SNPs (target 1.5M)" | tee -a "$LOG_FILE"
    plink \
        --bfile "${PLINK_DIR}/${LABEL}" \
        --chr-set ${N_CHR} \
        --thin "$THIN_NE2" \
        --recode \
        --out "${PLINK_DIR}/${LABEL}_ne2" \
        >> "$LOG_FILE" 2>&1
    N_SNP_NE2=$(wc -l < "$NE2_MAP")
    echo "INFO: ne2 PED: ${N_SNP_NE2} SNPs retained" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 5: currentNe2 — default (panmictic)
# Input: ne2-thinned PED (~1.5M SNPs, safely below MAXLOCI=2M).
# -s 200000: internal subsampling to 200k SNPs for computation;
#   200k → ~20B pairs → ~1.6h/run at 8 threads (vs ~89h at 1.5M
#   effective). Precision for Ne estimation is excellent at 200k.
# -k 0: random-mating assumption (incompatible with -x).
# Must run from PLINK_DIR so the MAP file is co-located with PED.
# -------------------------------------------------------------------
echo "STEP 5: currentNe2 default (panmictic assumption)" | tee -a "$LOG_FILE"

CNE2_DEFAULT="${CNE2_DIR}/${LABEL}_currentNe2_default.txt"
if [ -s "$CNE2_DEFAULT" ]; then
    echo "INFO: reusing existing currentNe2 default output" | tee -a "$LOG_FILE"
else
    (
        cd "$PLINK_DIR"
        currentne2 \
            -s 200000 \
            -r 1 \
            -t "${CPUS}" \
            -k 0 \
            -o "$CNE2_DEFAULT" \
            "${LABEL}_ne2.ped" >> "$LOG_FILE" 2>&1
    ) || echo "WARN: currentNe2 default returned non-zero — check log" | tee -a "$LOG_FILE"
    echo "INFO: currentNe2 default done" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 6: currentNe2 -x (population structure)
# Estimates: Ne_metapopulation, NT (total Ne across subpopulations),
# FST (pairwise average), m (migration rate), s (number of
# inferred subpopulations). The estimated s is directly comparable
# to ADMIXTURE K_best as an independent cross-validation of structure.
# -s 200000: same rationale as Step 5.
# -k incompatible with -x (omitted per Santiago et al. 2025 README).
# -------------------------------------------------------------------
echo "STEP 6: currentNe2 -x (structure inference)" | tee -a "$LOG_FILE"

CNE2_STRUCT="${CNE2_DIR}/${LABEL}_currentNe2_mix.txt"
if [ -s "$CNE2_STRUCT" ]; then
    echo "INFO: reusing existing currentNe2 -x output" | tee -a "$LOG_FILE"
else
    (
        cd "$PLINK_DIR"
        currentne2 \
            -x \
            -s 200000 \
            -r 1 \
            -t "${CPUS}" \
            -o "$CNE2_STRUCT" \
            "${LABEL}_ne2.ped" >> "$LOG_FILE" 2>&1
    ) || echo "WARN: currentNe2 -x returned non-zero — check log" | tee -a "$LOG_FILE"
    echo "INFO: currentNe2 -x done" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 7: GONE2 — default (panmictic temporal trajectory)
# Input: ne2-thinned PED (~1.5M SNPs), symlinked as ${LABEL}.ped so
# output name (${LABEL}.ped_GONE2_Ne) is preserved for the R parser.
# Run in GONE2/default/ to avoid filename collision with -x run.
# -------------------------------------------------------------------
echo "STEP 7: GONE2 default (panmictic trajectory)" | tee -a "$LOG_FILE"

GONE2_DEF_NE="${GONE2_DIR}/default/${LABEL}.ped_GONE2_Ne"
if [ -s "$GONE2_DEF_NE" ]; then
    echo "INFO: reusing existing GONE2 default output" | tee -a "$LOG_FILE"
else
    ln -sf "$NE2_MAP" "${GONE2_DIR}/default/${LABEL}.map"
    ln -sf "$NE2_PED" "${GONE2_DIR}/default/${LABEL}.ped"
    (
        cd "${GONE2_DIR}/default"
        gone2 -r 1 -t "${CPUS}" "${LABEL}.ped" >> "$LOG_FILE" 2>&1
    ) || echo "WARN: GONE2 default returned non-zero — check log" | tee -a "$LOG_FILE"
    echo "INFO: GONE2 default done" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 8: GONE2 -x (structure-corrected temporal trajectory)
# Run in GONE2/structure/ (separate subdir to avoid filename collision).
# -------------------------------------------------------------------
echo "STEP 8: GONE2 -x (structure-corrected trajectory)" | tee -a "$LOG_FILE"

GONE2_STR_NE="${GONE2_DIR}/structure/${LABEL}.ped_GONE2_Ne"
if [ -s "$GONE2_STR_NE" ]; then
    echo "INFO: reusing existing GONE2 -x output" | tee -a "$LOG_FILE"
else
    ln -sf "$NE2_MAP" "${GONE2_DIR}/structure/${LABEL}.map"
    ln -sf "$NE2_PED" "${GONE2_DIR}/structure/${LABEL}.ped"
    (
        cd "${GONE2_DIR}/structure"
        gone2 -x -r 1 -t "${CPUS}" "${LABEL}.ped" >> "$LOG_FILE" 2>&1
    ) || echo "WARN: GONE2 -x returned non-zero — check log" | tee -a "$LOG_FILE"
    echo "INFO: GONE2 -x done" | tee -a "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 9: R — compile results and produce summary tables + plots
# -------------------------------------------------------------------
echo "STEP 9: R compilation and plots" | tee -a "$LOG_FILE"

Rscript --vanilla "$R_SCRIPT" \
    "$OUT_DIR" \
    "$LABEL" >> "$LOG_FILE" 2>&1 \
    || echo "WARN: R script returned non-zero — check log" | tee -a "$LOG_FILE"

echo "DONE $(date -Iseconds) 07.10 Ne global completed" | tee -a "$LOG_FILE"
echo "Outputs: ${OUT_DIR}/"

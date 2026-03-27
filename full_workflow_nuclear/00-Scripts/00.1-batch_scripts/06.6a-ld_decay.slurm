#!/bin/bash
#SBATCH -p workq
#SBATCH --time=06:00:00
#SBATCH -J 06.6a-ld_decay
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/PLINK/2.00a4
module load bioinfo/Bcftools/1.9
module load bioinfo/VCFtools/0.1.16
module load bioinfo/samtools/1.23
module load statistics/R/4.3.0

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
# Input: TE-filtered VCF from 06.5 (17,436,796 SNPs, 164 samples)
# Computing LD on all 17M SNPs would be prohibitive in time and memory.
# We therefore work on a random subsample of SNPs to estimate the
# genome-wide LD decay curve reliably. 20,000 SNPs provides sufficient
# pairs within a 1 Mb window to fit a smooth decay curve.
IN_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_full/06.5-transposable_element_filtering"
IN_VCF="${IN_DIR}/snps_TE_filtered.vcf.gz"
IN_TBI="${IN_VCF}.tbi"

OUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_full/06.6a-ld_decay"
LOG_DIR="${OUT_DIR}/logs"
PLOT_DIR="${OUT_DIR}/plots"
LOG_FILE="${LOG_DIR}/06.6a-ld_decay.log"

R_SCRIPT="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.2-R_scripts/06.6a-ld_decay.R"

SUBSET_PREFIX="${OUT_DIR}/snps_subset_for_ld"
SUBSET_VCF="${SUBSET_PREFIX}.vcf.gz"
LD_PREFIX="${OUT_DIR}/ld_decay"

# Number of SNPs randomly sampled for LD computation.
# 20,000 SNPs gives ~250,000 pairs within a 1 Mb window — sufficient
# for a reliable decay curve while keeping runtime under 1 hour.
N_SNPS_SUBSET="${N_SNPS_SUBSET:-20000}"

# Maximum distance (bp) over which pairwise r² is computed.
# 1 Mb is large enough to observe full LD decay in most tree species,
# where LD typically decays within 1–100 kb.
LD_WINDOW_BP="${LD_WINDOW_BP:-1000000}"

PLINK_THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$PLOT_DIR"
echo "START $(date -Iseconds) 06.6a LD decay estimation" > "$LOG_FILE"

if [ ! -s "$IN_VCF" ] || [ ! -s "$IN_TBI" ]; then
    echo "ERROR: input VCF missing: $IN_VCF" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Random subsample of SNPs
# -------------------------------------------------------------------
# Why subsample?
# Computing all pairwise r² for 17M SNPs is computationally intractable.
# A random subsample of 20,000 SNPs evenly represents the full SNP
# landscape and is sufficient to estimate LD decay across the genome.
# --thin-count draws exactly N_SNPS_SUBSET SNPs randomly across all contigs.
echo "STEP 2: subsampling ${N_SNPS_SUBSET} random SNPs" | tee -a "$LOG_FILE"
plink2 \
    --vcf "$IN_VCF" \
    --allow-extra-chr \
    --threads "$PLINK_THREADS" \
    --thin-count "$N_SNPS_SUBSET" \
    --recode vcf bgz \
    --out "$SUBSET_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "${SUBSET_VCF}" ]; then
    echo "ERROR: subset VCF missing" | tee -a "$LOG_FILE"
    exit 1
fi

# plink2 exports VCF v4.3 but vcftools 0.1.16 only supports up to v4.2.
# Patch the header to v4.2 — the format differences are minor and
# do not affect the genotype data used by vcftools --geno-r2.
SUBSET_VCF_V42="${SUBSET_PREFIX}.v42.vcf.gz"
zcat "$SUBSET_VCF" \
    | sed 's/##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
    | bgzip -c > "$SUBSET_VCF_V42"
bcftools index --threads "$PLINK_THREADS" -t "$SUBSET_VCF_V42" >> "$LOG_FILE" 2>&1

echo "INFO: subset VCF created: $(bcftools view -H "$SUBSET_VCF_V42" | wc -l) SNPs" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Compute pairwise r² with vcftools --geno-r2
# -------------------------------------------------------------------
# vcftools --geno-r2 computes the squared correlation coefficient (r²)
# between all pairs of sites within a given distance window.
# Output format: CHR  POS1  POS2  N_INDV  R^2
#
# Why r² and not D'?
# r² is the standard metric for LD decay analyses in population genomics:
#   - Directly interpretable as proportion of variance explained
#   - Decreases monotonically with distance (D' does not)
#   - Independent of allele frequency (unlike D)
#   - Threshold of r² = 0.1–0.2 is typically used to define effective
#     LD block boundaries and to set the pruning window
#
# --min-r2 0: keep all pairs including those with r²=0 to fit
# the full decay curve including the asymptote.
echo "STEP 3: computing pairwise r² within ${LD_WINDOW_BP} bp" | tee -a "$LOG_FILE"
vcftools \
    --gzvcf "$SUBSET_VCF_V42" \
    --geno-r2 \
    --ld-window-bp "$LD_WINDOW_BP" \
    --min-r2 0.0 \
    --out "$LD_PREFIX" >> "$LOG_FILE" 2>&1

if [ ! -s "${LD_PREFIX}.geno.ld" ]; then
    echo "ERROR: LD output file missing: ${LD_PREFIX}.geno.ld" | tee -a "$LOG_FILE"
    exit 1
fi

N_PAIRS=$(tail -n +2 "${LD_PREFIX}.geno.ld" | wc -l)
echo "INFO: ${N_PAIRS} SNP pairs computed" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 4: Plot LD decay curve and estimate window recommendation
# -------------------------------------------------------------------
# The R script bins pairwise r² values by distance, fits a smooth
# decay curve, and identifies the distance at which mean r² drops
# below 0.2 (standard threshold used to set the LD pruning window).
# Output:
#   - ld_decay_curve.png : r² vs distance plot
#   - ld_window_recommendation.txt : recommended window for 06.6
echo "STEP 4: plotting LD decay curve" | tee -a "$LOG_FILE"
Rscript --vanilla "$R_SCRIPT" \
    "${LD_PREFIX}.geno.ld" \
    "$PLOT_DIR" >> "$LOG_FILE" 2>&1

echo "DONE $(date -Iseconds) 06.6a completed" | tee -a "$LOG_FILE"
echo "LD table  : ${LD_PREFIX}.geno.ld"
echo "Plots     : $PLOT_DIR"
echo "See ld_window_recommendation.txt for suggested window size for 06.6"

# Summary of steps:
# Step 1: Input TE-filtered VCF (17,436,796 SNPs) from 06.5
# Step 2: Random subsample of 20,000 SNPs (plink2 --thin-count)
# Step 3: Pairwise r² within 1 Mb (vcftools --geno-r2)
# Step 4: R script — bin r² by distance, smooth curve, window recommendation

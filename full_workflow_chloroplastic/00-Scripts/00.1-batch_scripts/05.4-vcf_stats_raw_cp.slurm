#!/bin/bash
#SBATCH -p workq
#SBATCH --time=00:15:00
#SBATCH -J cp-05.4-vcf_stats_raw
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load bioinfo/Bcftools/1.9

WORKFLOW="workflow_full"

VCF="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/05-variant_calling/${WORKFLOW}/genotypegvcfs/jointgenotypes.raw.vcf.gz"
OUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/05-variant_calling/${WORKFLOW}/genotypegvcfs"

echo "=== TOTAL SITES ==="
bcftools view -H "$VCF" | wc -l

echo "=== SNPs ==="
bcftools view -H -v snps "$VCF" | wc -l

echo "=== INDELs ==="
bcftools view -H -v indels "$VCF" | wc -l

echo "=== SITES MULTIALLELIQUES ==="
bcftools view -H -m3 "$VCF" | wc -l

echo "=== STATS GENERALES (bcftools stats) ==="
bcftools stats "$VCF" > "${OUT_DIR}/jointgenotypes.raw.stats.txt"
grep "^SN" "${OUT_DIR}/jointgenotypes.raw.stats.txt"

echo "DONE stats saved to ${OUT_DIR}/jointgenotypes.raw.stats.txt"

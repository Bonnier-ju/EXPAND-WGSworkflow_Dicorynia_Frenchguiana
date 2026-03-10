#!/bin/bash
#SBATCH -p workq
#SBATCH --time=00:30:00
#SBATCH -J cp-06.7-filtration_summary
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load bioinfo/Bcftools/1.9

# -------------------------------------------------------------------
# STEP 0: Paramètres
# -------------------------------------------------------------------
WORKFLOW="workflow_full"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
OUT_DIR="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/filtration_summary"
mkdir -p "$OUT_DIR"

SUMMARY_TSV="${OUT_DIR}/filtration_summary.tsv"
SUMMARY_TXT="${OUT_DIR}/filtration_summary.txt"

# VCFs à chaque étape
VCF_RAW="${BASE_DIR}/05-variant_calling/${WORKFLOW}/genotypegvcfs/jointgenotypes.raw.vcf.gz"
VCF_HARD="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/hard_filter/variants.hard_filtered.pass.vcf.gz"
VCF_NORM="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/normalization/variants.normalized.vcf.gz"
VCF_FILT="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/missing_data_filter/variants.filtered.vcf.gz"

# -------------------------------------------------------------------
# Fonction : compte SNPs / INDELs / total depuis un VCF
# -------------------------------------------------------------------
count_variants() {
    local vcf="$1"
    if [ ! -s "$vcf" ]; then
        echo "NA NA NA NA NA"
        return
    fi
    local total snps indels multi missing_pct
    total=$(bcftools view -H "$vcf" | wc -l)
    snps=$(bcftools view -H -v snps "$vcf" | wc -l)
    indels=$(bcftools view -H -v indels "$vcf" | wc -l)
    multi=$(bcftools view -H -m3 "$vcf" | wc -l)
    # taux de missing moyen sur les sites (INFO/AN vs INFO/AN_max)
    # bcftools stats pour missing
    missing_pct=$(bcftools stats "$vcf" | awk '/^SN.*number of no-ALTs/{print $4}')
    echo "$total $snps $indels $multi"
}

# -------------------------------------------------------------------
# STEP 1: Compter les variants à chaque étape
# -------------------------------------------------------------------
echo "INFO: counting variants at each step..."

read total_raw   snps_raw   indels_raw   multi_raw   <<< $(count_variants "$VCF_RAW")
read total_hard  snps_hard  indels_hard  multi_hard  <<< $(count_variants "$VCF_HARD")
read total_norm  snps_norm  indels_norm  multi_norm  <<< $(count_variants "$VCF_NORM")
read total_filt  snps_filt  indels_filt  multi_filt  <<< $(count_variants "$VCF_FILT")

# -------------------------------------------------------------------
# STEP 2: Calcul des pertes relatives
# -------------------------------------------------------------------
pct_lost() {
    local before="$1" after="$2"
    if [ "$before" = "NA" ] || [ "$after" = "NA" ] || [ "$before" -eq 0 ]; then
        echo "NA"
    else
        awk -v b="$before" -v a="$after" 'BEGIN{printf "%.1f", (b-a)/b*100}'
    fi
}

lost_hard=$(pct_lost "$total_raw"  "$total_hard")
lost_norm=$(pct_lost "$total_hard" "$total_norm")
lost_filt=$(pct_lost "$total_norm" "$total_filt")
lost_total=$(pct_lost "$total_raw" "$total_filt")

snp_lost_hard=$(pct_lost "$snps_raw"  "$snps_hard")
snp_lost_norm=$(pct_lost "$snps_hard" "$snps_norm")
snp_lost_filt=$(pct_lost "$snps_norm" "$snps_filt")
snp_lost_total=$(pct_lost "$snps_raw" "$snps_filt")

# -------------------------------------------------------------------
# STEP 3: Écriture du tableau TSV
# -------------------------------------------------------------------
{
echo -e "step\tdescription\ttotal_sites\tSNPs\tINDELs\tmultiallelic\tlost_vs_prev_pct"
echo -e "05.3_raw\tGenotypeGVCFs raw VCF\t${total_raw}\t${snps_raw}\t${indels_raw}\t${multi_raw}\t-"
echo -e "06.1_hard_filter\tGATK VariantFiltration (PASS only)\t${total_hard}\t${snps_hard}\t${indels_hard}\t${multi_hard}\t${lost_hard}%"
echo -e "06.2_normalization\tbcftools norm (split multiallelic)\t${total_norm}\t${snps_norm}\t${indels_norm}\t${multi_norm}\t${lost_norm}%"
echo -e "06.3_missing_filter\tF_MISSING<0.2 + MAC>=1\t${total_filt}\t${snps_filt}\t${indels_filt}\t${multi_filt}\t${lost_filt}%"
} > "$SUMMARY_TSV"

# -------------------------------------------------------------------
# STEP 4: Rapport lisible en texte
# -------------------------------------------------------------------
{
echo "================================================================"
echo " FILTRATION SUMMARY — Chloroplast Workflow — ${WORKFLOW}"
echo " $(date)"
echo "================================================================"
echo ""
echo "VARIANT COUNTS AT EACH FILTERING STEP"
echo "----------------------------------------------------------------"
printf "%-30s %10s %10s %10s %10s %10s\n" "Step" "Total" "SNPs" "INDELs" "Multi." "Loss(%)"
echo "----------------------------------------------------------------"
printf "%-30s %10s %10s %10s %10s %10s\n" "05.3 raw VCF"            "$total_raw"  "$snps_raw"  "$indels_raw"  "$multi_raw"  "-"
printf "%-30s %10s %10s %10s %10s %10s\n" "06.1 hard filter (PASS)" "$total_hard" "$snps_hard" "$indels_hard" "$multi_hard" "${lost_hard}%"
printf "%-30s %10s %10s %10s %10s %10s\n" "06.2 normalization"      "$total_norm" "$snps_norm" "$indels_norm" "$multi_norm" "${lost_norm}%"
printf "%-30s %10s %10s %10s %10s %10s\n" "06.3 missing+MAC filter" "$total_filt" "$snps_filt" "$indels_filt" "$multi_filt" "${lost_filt}%"
echo "----------------------------------------------------------------"
printf "%-30s %10s %10s %10s %10s %10s\n" "TOTAL LOSS (raw → final)" \
    "$([ "$total_raw" != "NA" ] && [ "$total_filt" != "NA" ] && echo "$((total_raw - total_filt))" || echo "NA")" \
    "$([ "$snps_raw"  != "NA" ] && [ "$snps_filt"  != "NA" ] && echo "$((snps_raw  - snps_filt))"  || echo "NA")" \
    "-" "-" "${lost_total}%"
echo ""
echo "SNP RETENTION"
echo "----------------------------------------------------------------"
printf "  Raw SNPs          : %s\n"  "$snps_raw"
printf "  After hard filter : %s  (-%s%%)\n" "$snps_hard" "$snp_lost_hard"
printf "  After norm        : %s  (-%s%%)\n" "$snps_norm" "$snp_lost_norm"
printf "  After missing/MAC : %s  (-%s%%)\n" "$snps_filt" "$snp_lost_filt"
printf "  Overall SNP loss  : %s%%\n" "$snp_lost_total"
echo ""
echo "FILTER THRESHOLDS APPLIED"
echo "----------------------------------------------------------------"
echo "  06.1 SNP  : QD<2 | FS>60 | MQ<40 | SOR>3 | MQRankSum<-12.5 | ReadPosRankSum<-8"
echo "  06.1 INDEL: QD<2 | FS>200 | SOR>10 | ReadPosRankSum<-20"
echo "  06.2      : bcftools norm -m -any (split multiallelic) + norm -f REF"
echo "  06.3      : F_MISSING > 0.2 removed | MAC < 1 removed"
echo "================================================================"
} > "$SUMMARY_TXT"

# Affichage dans le log SLURM
cat "$SUMMARY_TXT"

echo ""
echo "DONE outputs:"
echo "  TSV : $SUMMARY_TSV"
echo "  TXT : $SUMMARY_TXT"

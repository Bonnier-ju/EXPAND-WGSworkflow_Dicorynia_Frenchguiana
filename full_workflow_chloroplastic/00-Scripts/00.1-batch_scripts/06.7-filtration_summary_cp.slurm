#!/bin/bash
#SBATCH -p workq
#SBATCH --time=01:00:00
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
# STEP 0: Parameters
# -------------------------------------------------------------------
WORKFLOW="workflow_full"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
OUT_DIR="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/filtration_summary"
mkdir -p "$OUT_DIR"

SUMMARY_TSV="${OUT_DIR}/filtration_summary.tsv"
SUMMARY_TXT="${OUT_DIR}/filtration_summary.txt"

# VCFs at each filtering step
VCF_RAW="${BASE_DIR}/05-variant_calling/${WORKFLOW}/genotypegvcfs/jointgenotypes.raw.vcf.gz"
VCF_HARD="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/hard_filter/variants.hard_filtered.pass.vcf.gz"
VCF_NORM="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/normalization/variants.normalized.vcf.gz"
VCF_FILT="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/missing_data_filter/variants.filtered.vcf.gz"

# Hard filter log (contains per-filter counts from GATK)
HARD_FILTER_LOG="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/hard_filter/hard_filter_cp.log"
# Missing data filter log (contains excluded samples list)
MISSING_FILTER_LOG="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/missing_data_filter/missing_data_filter_cp.log"
SAMPLES_EXCLUDED="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/missing_data_filter/samples_high_missing.txt"

# -------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------

# Count SNPs / INDELs / multiallelic / singletons from a VCF
# Returns: total snps indels multi singletons tstv
vcf_stats() {
    local vcf="$1"
    if [ ! -s "$vcf" ]; then
        echo "NA NA NA NA NA NA"
        return
    fi
    local total snps indels multi singletons ts tv tstv
    total=$(bcftools view -H "$vcf" | wc -l)
    snps=$(bcftools view -H -v snps "$vcf" | wc -l)
    indels=$(bcftools view -H -v indels "$vcf" | wc -l)
    multi=$(bcftools view -H -m3 "$vcf" | wc -l)
    # Singletons: SNP sites where AC=1 (seen in exactly 1 sample)
    singletons=$(bcftools view -v snps "$vcf" | bcftools view -H -c 1:minor -C 1:minor 2>/dev/null | wc -l || echo "NA")
    # Ts/Tv from bcftools stats
    tstv=$(bcftools stats -r Contig1 "$vcf" 2>/dev/null | awk '/^TSTV/{printf "%.2f", $5}')
    [ -z "$tstv" ] && tstv="NA"
    echo "$total $snps $indels $multi $singletons $tstv"
}

# Percentage lost between two values
pct_lost() {
    local before="$1" after="$2"
    if [ "$before" = "NA" ] || [ "$after" = "NA" ] || [ "$before" -eq 0 ] 2>/dev/null; then
        echo "NA"
    else
        awk -v b="$before" -v a="$after" 'BEGIN{printf "%.1f", (b-a)/b*100}'
    fi
}

# -------------------------------------------------------------------
# STEP 1: Collect stats at each step
# -------------------------------------------------------------------
echo "INFO: collecting stats for each VCF..."

read total_raw   snps_raw   indels_raw   multi_raw   sing_raw   tstv_raw   <<< $(vcf_stats "$VCF_RAW")
read total_hard  snps_hard  indels_hard  multi_hard  sing_hard  tstv_hard  <<< $(vcf_stats "$VCF_HARD")
read total_norm  snps_norm  indels_norm  multi_norm  sing_norm  tstv_norm  <<< $(vcf_stats "$VCF_NORM")
read total_filt  snps_filt  indels_filt  multi_filt  sing_filt  tstv_filt  <<< $(vcf_stats "$VCF_FILT")

# -------------------------------------------------------------------
# STEP 2: Compute losses
# -------------------------------------------------------------------
lost_hard=$(pct_lost "$total_raw"  "$total_hard")
lost_norm=$(pct_lost "$total_hard" "$total_norm")
lost_filt=$(pct_lost "$total_norm" "$total_filt")
lost_total=$(pct_lost "$total_raw" "$total_filt")

snp_lost_hard=$(pct_lost "$snps_raw"  "$snps_hard")
snp_lost_norm=$(pct_lost "$snps_hard" "$snps_norm")
snp_lost_filt=$(pct_lost "$snps_norm" "$snps_filt")
snp_lost_total=$(pct_lost "$snps_raw" "$snps_filt")

# -------------------------------------------------------------------
# STEP 3: 06.1 — per-filter breakdown
# Count how many SNPs failed each individual filter (from merged VCF)
# -------------------------------------------------------------------
echo "INFO: computing per-filter breakdown for 06.1..."
VCF_SNP_FILTERED="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/hard_filter/snps.filtered.vcf.gz"
VCF_INDEL_FILTERED="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/hard_filter/indels.filtered.vcf.gz"

snp_fail_QD="NA"; snp_fail_FS="NA"; snp_fail_MQ="NA"
snp_fail_SOR="NA"; snp_fail_MQR="NA"; snp_fail_RPR="NA"
indel_fail_QD="NA"; indel_fail_FS="NA"; indel_fail_SOR="NA"; indel_fail_RPR="NA"

if [ -s "$VCF_SNP_FILTERED" ]; then
    snp_fail_QD=$(bcftools view -H "$VCF_SNP_FILTERED" | awk '$7~/QD2/'        | wc -l)
    snp_fail_FS=$(bcftools view -H "$VCF_SNP_FILTERED" | awk '$7~/FS60/'       | wc -l)
    snp_fail_MQ=$(bcftools view -H "$VCF_SNP_FILTERED" | awk '$7~/MQ40/'       | wc -l)
    snp_fail_SOR=$(bcftools view -H "$VCF_SNP_FILTERED" | awk '$7~/SOR3/'      | wc -l)
    snp_fail_MQR=$(bcftools view -H "$VCF_SNP_FILTERED" | awk '$7~/MQRankSum/' | wc -l)
    snp_fail_RPR=$(bcftools view -H "$VCF_SNP_FILTERED" | awk '$7~/ReadPos/'   | wc -l)
fi
if [ -s "$VCF_INDEL_FILTERED" ]; then
    indel_fail_QD=$(bcftools view -H "$VCF_INDEL_FILTERED" | awk '$7~/QD2/'    | wc -l)
    indel_fail_FS=$(bcftools view -H "$VCF_INDEL_FILTERED" | awk '$7~/FS200/'  | wc -l)
    indel_fail_SOR=$(bcftools view -H "$VCF_INDEL_FILTERED" | awk '$7~/SOR10/' | wc -l)
    indel_fail_RPR=$(bcftools view -H "$VCF_INDEL_FILTERED" | awk '$7~/ReadPos/' | wc -l)
fi

# -------------------------------------------------------------------
# STEP 4: 06.2 — normalization details
# -------------------------------------------------------------------
snps_gained_norm="NA"
if [ "$snps_norm" != "NA" ] && [ "$snps_hard" != "NA" ]; then
    snps_gained_norm=$((snps_norm - snps_hard))
fi
multi_split="NA"
if [ "$multi_hard" != "NA" ]; then
    multi_split=$multi_hard
fi

# -------------------------------------------------------------------
# STEP 5: 06.3 — excluded samples and missing data stats
# -------------------------------------------------------------------
n_samples_raw=$(bcftools query -l "$VCF_FILT" 2>/dev/null | wc -l || echo "NA")
n_excluded=0
excluded_list="none"
if [ -s "$SAMPLES_EXCLUDED" ]; then
    n_excluded=$(wc -l < "$SAMPLES_EXCLUDED")
    excluded_list=$(paste -sd ',' "$SAMPLES_EXCLUDED")
fi

# Average missing rate in final VCF
avg_missing_filt=$(bcftools stats -s - "$VCF_FILT" 2>/dev/null | \
    awk '/^PSC/{total+=$4+$5+$6; miss+=$6} END{if(total>0) printf "%.1f%%", miss/total*100; else print "NA"}' || echo "NA")

# -------------------------------------------------------------------
# STEP 6: Write TSV summary
# -------------------------------------------------------------------
{
echo -e "step\tdescription\ttotal_sites\tSNPs\tINDELs\tmultiallelic\tsingletons\tTs_Tv\tlost_vs_prev_pct"
echo -e "05.3_raw\tGenotypeGVCFs raw VCF\t${total_raw}\t${snps_raw}\t${indels_raw}\t${multi_raw}\t${sing_raw}\t${tstv_raw}\t-"
echo -e "06.1_hard_filter\tGATK VariantFiltration PASS only\t${total_hard}\t${snps_hard}\t${indels_hard}\t${multi_hard}\t${sing_hard}\t${tstv_hard}\t-${lost_hard}%"
echo -e "06.2_normalization\tbcftools norm split+left-align\t${total_norm}\t${snps_norm}\t${indels_norm}\t${multi_norm}\t${sing_norm}\t${tstv_norm}\t-${lost_norm}%"
echo -e "06.3_missing_MAC\tF_MISSING<0.2 + MAC>=1\t${total_filt}\t${snps_filt}\t${indels_filt}\t${multi_filt}\t${sing_filt}\t${tstv_filt}\t-${lost_filt}%"
} > "$SUMMARY_TSV"

# -------------------------------------------------------------------
# STEP 7: Write detailed text report
# -------------------------------------------------------------------
{
echo "=================================================================="
echo " FILTRATION SUMMARY — Chloroplast Workflow — ${WORKFLOW}"
echo " $(date)"
echo "=================================================================="
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " OVERVIEW — VARIANT COUNTS AT EACH STEP"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
printf "%-28s %8s %8s %8s %8s %10s %7s %8s\n" \
    "Step" "Total" "SNPs" "INDELs" "Multi." "Singletons" "Ts/Tv" "Loss(%)"
echo "------------------------------------------------------------------"
printf "%-28s %8s %8s %8s %8s %10s %7s %8s\n" \
    "05.3 raw VCF"            "$total_raw"  "$snps_raw"  "$indels_raw"  "$multi_raw"  "$sing_raw"  "$tstv_raw"  "-"
printf "%-28s %8s %8s %8s %8s %10s %7s %8s\n" \
    "06.1 hard filter (PASS)" "$total_hard" "$snps_hard" "$indels_hard" "$multi_hard" "$sing_hard" "$tstv_hard" "-${lost_hard}%"
printf "%-28s %8s %8s %8s %8s %10s %7s %8s\n" \
    "06.2 normalization"      "$total_norm" "$snps_norm" "$indels_norm" "$multi_norm" "$sing_norm" "$tstv_norm" "-${lost_norm}%"
printf "%-28s %8s %8s %8s %8s %10s %7s %8s\n" \
    "06.3 missing+MAC filter" "$total_filt" "$snps_filt" "$indels_filt" "$multi_filt" "$sing_filt" "$tstv_filt" "-${lost_filt}%"
echo "------------------------------------------------------------------"
printf "%-28s %8s %8s\n" "TOTAL LOSS (raw → final)" \
    "$([ "$total_raw" != "NA" ] && [ "$total_filt" != "NA" ] && echo "$((total_raw - total_filt))" || echo "NA")" \
    "-${lost_total}%"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " STEP 06.1 — GATK HARD FILTER DETAILS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " SNP filters applied (GATK BEST PRACTICES — haploid chloroplast):"
printf "   %-30s %8s  %s\n" "Filter" "Rejected" "Threshold"
printf "   %-30s %8s  %s\n" "QD (Quality by Depth)"        "$snp_fail_QD"  "QD < 2.0"
printf "   %-30s %8s  %s\n" "FS (Strand Bias Fisher)"      "$snp_fail_FS"  "FS > 60.0"
printf "   %-30s %8s  %s\n" "MQ (Mapping Quality RMS)"     "$snp_fail_MQ"  "MQ < 40.0"
printf "   %-30s %8s  %s\n" "SOR (Strand Odds Ratio)"      "$snp_fail_SOR" "SOR > 3.0"
printf "   %-30s %8s  %s\n" "MQRankSum"                    "$snp_fail_MQR" "MQRankSum < -12.5"
printf "   %-30s %8s  %s\n" "ReadPosRankSum"               "$snp_fail_RPR" "ReadPosRankSum < -8.0"
echo ""
echo " INDEL filters applied:"
printf "   %-30s %8s  %s\n" "QD (Quality by Depth)"        "$indel_fail_QD"  "QD < 2.0"
printf "   %-30s %8s  %s\n" "FS (Strand Bias Fisher)"      "$indel_fail_FS"  "FS > 200.0"
printf "   %-30s %8s  %s\n" "SOR (Strand Odds Ratio)"      "$indel_fail_SOR" "SOR > 10.0"
printf "   %-30s %8s  %s\n" "ReadPosRankSum"               "$indel_fail_RPR" "ReadPosRankSum < -20.0"
echo ""
printf "   SNPs retained after hard filter : %s / %s (-%s%%)\n" "$snps_hard" "$snps_raw" "$snp_lost_hard"
printf "   Ts/Tv ratio improvement         : %s → %s\n" "$tstv_raw" "$tstv_hard"
echo "   Note: sites may fail multiple filters simultaneously"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " STEP 06.2 — NORMALIZATION DETAILS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
printf "   Multiallelic sites split        : %s\n"     "$multi_split"
printf "   SNP records before split        : %s\n"     "$snps_hard"
printf "   SNP records after split         : %s\n"     "$snps_norm"
printf "   Net SNP records gained          : %s\n"     "$snps_gained_norm"
printf "   Multiallelic sites remaining    : %s\n"     "$multi_norm"
printf "   Ts/Tv after normalization       : %s\n"     "$tstv_norm"
echo "   Operations: bcftools norm -m -any (split) + bcftools norm -f REF (left-align indels)"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " STEP 06.3 — MISSING DATA + MAC FILTER DETAILS"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
printf "   Thresholds applied:\n"
printf "     Site missing rate threshold   : 20%% (F_MISSING > 0.2 removed)\n"
printf "     Sample missing rate threshold : 50%% (samples > 50%% missing excluded)\n"
printf "     Min minor allele count (MAC)  : 1   (monomorphic sites removed)\n"
echo ""
printf "   Samples before filter  : %s\n"  "$((n_samples_raw + n_excluded))"
printf "   Samples excluded       : %s\n"  "$n_excluded"
printf "   Excluded sample(s)     : %s\n"  "$excluded_list"
printf "   Samples retained       : %s\n"  "$n_samples_raw"
printf "   Avg missing rate final : %s\n"  "$avg_missing_filt"
echo ""
printf "   SNPs before filter     : %s\n"  "$snps_norm"
printf "   SNPs after filter      : %s  (-%s%%)\n" "$snps_filt" "$snp_lost_filt"
printf "   Singletons retained    : %s  (MAC=1 kept)\n" "$sing_filt"
printf "   Ts/Tv after filter     : %s\n"  "$tstv_filt"
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo " SNP RETENTION SUMMARY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
printf "   Raw SNPs (05.3)        : %s\n"        "$snps_raw"
printf "   After 06.1 hard filter : %s  (-%s%%)\n" "$snps_hard" "$snp_lost_hard"
printf "   After 06.2 norm        : %s  (+%s records from split)\n" "$snps_norm" "$snps_gained_norm"
printf "   After 06.3 missing/MAC : %s  (-%s%%)\n" "$snps_filt" "$snp_lost_filt"
printf "   Overall SNP loss       : -%s%%\n"      "$snp_lost_total"
echo ""

echo "=================================================================="
echo " FINAL VCF READY FOR HAPLOTYPE ANALYSIS (06.4) AND CONSENSUS (06.5)"
printf "   File : %s\n" "$VCF_FILT"
printf "   SNPs : %s | INDELs : %s | Samples : %s | Ts/Tv : %s\n" \
    "$snps_filt" "$indels_filt" "$n_samples_raw" "$tstv_filt"
echo "=================================================================="

} > "$SUMMARY_TXT"

# Print to SLURM log
cat "$SUMMARY_TXT"

echo ""
echo "DONE outputs:"
echo "  TSV : $SUMMARY_TSV"
echo "  TXT : $SUMMARY_TXT"

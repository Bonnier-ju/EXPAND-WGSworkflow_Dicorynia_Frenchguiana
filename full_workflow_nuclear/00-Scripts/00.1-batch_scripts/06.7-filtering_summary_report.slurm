#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J 06.7-filter_summary
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load bioinfo/Bcftools/1.9

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test"
OUT_DIR="${BASE}/06.7-filtering_summary"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

SUMMARY_TSV="${OUT_DIR}/06.7-filtering_summary.tsv"
REASONS_TSV="${OUT_DIR}/06.7-filtering_reasons.tsv"
SUMMARY_TXT="${OUT_DIR}/06.7-filtering_summary.txt"
LOG_FILE="${LOG_DIR}/06.7-filtering_summary.log"

echo "START $(date -Iseconds) 06.7 filtering summary" > "$LOG_FILE"

count_variants() {
    local vcf="$1"
    if [ -s "$vcf" ]; then
        bcftools view -H "$vcf" | wc -l
    else
        echo "NA"
    fi
}

count_samples_vcf() {
    local vcf="$1"
    if [ -s "$vcf" ]; then
        bcftools query -l "$vcf" | wc -l
    else
        echo "NA"
    fi
}

to_int_or_zero() {
    case "$1" in
        ''|NA) echo 0 ;;
        *) echo "$1" ;;
    esac
}

pct_removed() {
    local in_n="$1"
    local out_n="$2"
    if [ "$in_n" = "NA" ] || [ "$out_n" = "NA" ] || [ "$in_n" = "0" ]; then
        echo "NA"
    else
        awk -v i="$in_n" -v o="$out_n" 'BEGIN{printf "%.2f", ((i-o)/i)*100}'
    fi
}

write_step_row() {
    local step_id="$1"
    local input_file="$2"
    local output_file="$3"
    local n_samples="$4"
    local n_in="$5"
    local n_out="$6"
    local notes="$7"

    local in_i out_i removed pct
    in_i="$(to_int_or_zero "$n_in")"
    out_i="$(to_int_or_zero "$n_out")"
    if [ "$n_in" = "NA" ] || [ "$n_out" = "NA" ]; then
        removed="NA"
        pct="NA"
    else
        removed=$((in_i - out_i))
        pct="$(pct_removed "$n_in" "$n_out")"
    fi

    echo -e "${step_id}\t${input_file}\t${output_file}\t${n_samples}\t${n_in}\t${n_out}\t${removed}\t${pct}\t${notes}" >> "$SUMMARY_TSV"
}

echo -e "step_id\tinput_file\toutput_file\tn_samples\tn_variants_input\tn_variants_output\tn_removed\tpct_removed\tnotes" > "$SUMMARY_TSV"
echo -e "step_id\treason\tcount\tnotes" > "$REASONS_TSV"

# -----------------------------
# 06.1
# -----------------------------
S61_IN="${BASE}/05-variant_calling/workflow_test/genotypegvcfs/jointgenotypes.raw.vcf.gz"
S61_OUT="${BASE}/06.1-normalization_selection/joint.norm.biallelic.snps.vcf.gz"
S61_IN_N="$(count_variants "$S61_IN")"
S61_OUT_N="$(count_variants "$S61_OUT")"
S61_SAMPLES="$(count_samples_vcf "$S61_OUT")"
write_step_row "06.1" "$S61_IN" "$S61_OUT" "$S61_SAMPLES" "$S61_IN_N" "$S61_OUT_N" "Normalization + biallelic SNP selection"

# -----------------------------
# 06.2
# -----------------------------
S62_IN="${BASE}/06.1-normalization_selection/joint.norm.biallelic.snps.vcf.gz"
S62_OUT="${BASE}/06.2-technical_hard_filter/snps.hardfiltered.vcf.gz"
S62_TAGGED="${BASE}/06.2-technical_hard_filter/snps.filtered.tagged.vcf.gz"
S62_IN_N="$(count_variants "$S62_IN")"
S62_OUT_N="$(count_variants "$S62_OUT")"
S62_SAMPLES="$(count_samples_vcf "$S62_OUT")"
write_step_row "06.2" "$S62_IN" "$S62_OUT" "$S62_SAMPLES" "$S62_IN_N" "$S62_OUT_N" "Hard filtering (QD/QUAL/FS/SOR/MQ/MQRankSum/ReadPosRankSum/DP)"

if [ -s "$S62_TAGGED" ]; then
    bcftools query -f '%FILTER\n' "$S62_TAGGED" | tr ';' '\n' | sort | uniq -c | while read -r c filt; do
        [ -z "${filt:-}" ] && continue
        echo -e "06.2\t${filt}\t${c}\tFILTER tags from tagged VCF (non-exclusive counts)" >> "$REASONS_TSV"
    done
fi

# -----------------------------
# 06.3
# -----------------------------
S63_IN="${BASE}/06.2-technical_hard_filter/snps.hardfiltered.vcf.gz"
S63_OUT="${BASE}/06.3-missing_data_sample_filtering/snps.missing_filtered.vcf.gz"
S63_IMISS="${BASE}/06.3-missing_data_sample_filtering/snps.site_missing_filtered.imiss"
S63_RM="${BASE}/06.3-missing_data_sample_filtering/samples_to_remove_missing.txt"
S63_KEPT="${BASE}/06.3-missing_data_sample_filtering/samples_kept_missing.txt"
S63_IN_N="$(count_variants "$S63_IN")"
S63_OUT_N="$(count_variants "$S63_OUT")"
S63_SAMPLES="$(count_samples_vcf "$S63_OUT")"
write_step_row "06.3" "$S63_IN" "$S63_OUT" "$S63_SAMPLES" "$S63_IN_N" "$S63_OUT_N" "Missing data filtering: sites + sample missingness"

if [ -s "$S63_RM" ]; then
    echo -e "06.3\tsamples_removed_missingness\t$(wc -l < "$S63_RM")\tF_MISS above threshold" >> "$REASONS_TSV"
fi
if [ -s "$S63_KEPT" ]; then
    echo -e "06.3\tsamples_kept_after_missingness\t$(wc -l < "$S63_KEPT")\tF_MISS within threshold" >> "$REASONS_TSV"
fi
if [ -s "$S63_IMISS" ]; then
    echo -e "06.3\timiss_rows\t$(( $(wc -l < "$S63_IMISS") - 1 ))\tNumber of individuals evaluated in .imiss" >> "$REASONS_TSV"
fi

# -----------------------------
# 06.4
# -----------------------------
S64_IN="${BASE}/06.3-missing_data_sample_filtering/snps.missing_filtered.vcf.gz"
S64_OUT="${BASE}/06.4-allele_frequency_population_filter/snps_non_informative_filtered.vcf.gz"
S64_IN_N="$(count_variants "$S64_IN")"
S64_OUT_N="$(count_variants "$S64_OUT")"
S64_SAMPLES="$(count_samples_vcf "$S64_OUT")"
write_step_row "06.4" "$S64_IN" "$S64_OUT" "$S64_SAMPLES" "$S64_IN_N" "$S64_OUT_N" "MAC filtering (default MIN_MAC=3), remove non-informative/singleton/doubleton"
if [ "$S64_IN_N" != "NA" ] && [ "$S64_OUT_N" != "NA" ]; then
    echo -e "06.4\tmac_below_threshold_removed\t$((S64_IN_N - S64_OUT_N))\tRemoved by MAC threshold (and non-informative variants)" >> "$REASONS_TSV"
fi

# -----------------------------
# 06.5
# -----------------------------
S65_IN="${BASE}/06.4-allele_frequency_population_filter/snps_non_informative_filtered.vcf.gz"
S65_OUT="${BASE}/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
S65_TE_BED="${BASE}/06.5-transposable_element_filtering/te_regions.merged.bed"
S65_IN_N="$(count_variants "$S65_IN")"
S65_OUT_N="$(count_variants "$S65_OUT")"
S65_SAMPLES="$(count_samples_vcf "$S65_OUT")"
write_step_row "06.5" "$S65_IN" "$S65_OUT" "$S65_SAMPLES" "$S65_IN_N" "$S65_OUT_N" "TE overlap filtering using GFF-derived BED"
if [ -s "$S65_TE_BED" ]; then
    echo -e "06.5\tte_intervals_merged\t$(wc -l < "$S65_TE_BED")\tMerged TE intervals used for exclusion" >> "$REASONS_TSV"
fi
if [ "$S65_IN_N" != "NA" ] && [ "$S65_OUT_N" != "NA" ]; then
    echo -e "06.5\tvariants_overlapping_te_removed\t$((S65_IN_N - S65_OUT_N))\tVariants removed by TE overlap" >> "$REASONS_TSV"
fi

# -----------------------------
# 06.6
# -----------------------------
S66_IN="${BASE}/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
S66_OUT="${BASE}/06.6-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
S66_LOG="${BASE}/06.6-ld_pruning/logs/06.6-ld_pruning.log"
S66_IN_N="$(count_variants "$S66_IN")"
S66_OUT_N="$(count_variants "$S66_OUT")"
S66_SAMPLES="$(count_samples_vcf "$S66_IN")"

if [ -s "$S66_OUT" ]; then
    write_step_row "06.6" "$S66_IN" "$S66_OUT" "$S66_SAMPLES" "$S66_IN_N" "$S66_OUT_N" "LD pruning (--indep-pairwise)"
    if [ "$S66_IN_N" != "NA" ] && [ "$S66_OUT_N" != "NA" ]; then
        echo -e "06.6\tld_pruned_variants_removed\t$((S66_IN_N - S66_OUT_N))\tVariants removed by LD pruning" >> "$REASONS_TSV"
    fi
else
    write_step_row "06.6" "$S66_IN" "NA" "$S66_SAMPLES" "$S66_IN_N" "NA" "LD pruning not completed"
    if [ -s "$S66_LOG" ] && rg -q "less than 50 samples" "$S66_LOG"; then
        echo -e "06.6\tld_pruning_skipped_n_lt_50\t1\tPLINK2 refused LD estimation with <50 samples" >> "$REASONS_TSV"
    else
        echo -e "06.6\tld_pruning_missing_output\t1\tOutput VCF not found" >> "$REASONS_TSV"
    fi
fi

{
    echo "06.7 Filtering Summary Report"
    echo "Generated: $(date -Iseconds)"
    echo
    echo "Main summary table: $SUMMARY_TSV"
    echo "Reasons table     : $REASONS_TSV"
    echo
    awk -F'\t' 'NR==1 || $1=="06.1" || $1=="06.2" || $1=="06.3" || $1=="06.4" || $1=="06.5" || $1=="06.6" {print}' "$SUMMARY_TSV"
} > "$SUMMARY_TXT"

cat "$SUMMARY_TXT"
echo "DONE $(date -Iseconds) 06.7 filtering summary" | tee -a "$LOG_FILE"

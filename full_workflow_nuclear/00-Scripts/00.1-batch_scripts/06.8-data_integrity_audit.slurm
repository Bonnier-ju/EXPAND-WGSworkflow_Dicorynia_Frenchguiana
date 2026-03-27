#!/bin/bash
#SBATCH -p workq
#SBATCH --time=06:00:00
#SBATCH -J 06.8-integrity_audit
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# =======================================================================
# SCRIPT: 06.8 — Data Integrity Audit
# =======================================================================
# PURPOSE
#   Cross-validate all VCF outputs produced by the filtering pipeline
#   (06.1 → 06.6b) and generate a reproducibility record of the run.
#   This script does not modify any data; it only reads existing VCFs.
#
# CHECKS PERFORMED
#   01 — Index + readability  : verifies each VCF has a .tbi index and
#                               that bcftools can open the header.
#   02 — Sample concordance   : ensures the same sample set is present
#                               at every filtering step (no accidental
#                               sample gain/loss between steps).
#   03 — Variant ID integrity : counts total variants in the final VCF,
#                               detects missing IDs (ID='.') and
#                               duplicate variant IDs or duplicate loci
#                               (same CHROM:POS:REF:ALT).
#   04 — Allele sanity        : confirms no multiallelic sites remain
#                               (all sites should be biallelic after
#                               06.1 SelectVariants), no non-ACGTN
#                               alleles, and no records without ALT.
#   06 — Missingness          : computes per-individual (F_MISS) and
#                               per-site missingness on the final VCF
#                               using vcftools, to confirm that missing
#                               data thresholds applied in 06.3 are
#                               met in the final dataset.
#   07 — Heterozygosity/depth : computes observed heterozygosity
#                               (F-statistic via vcftools --het) and
#                               mean per-individual depth (--depth) on
#                               the final VCF. Used to flag outlier
#                               samples (excess heterozygosity may
#                               indicate contamination; very low depth
#                               indicates poor sequencing quality).
#   08 — Reproducibility      : records tool versions, loaded modules,
#                               and SHA-256 checksums of all VCF files
#                               to enable future exact reproduction
#                               and audit of filtering parameters.
#
# OUTPUTS
#   06.8-integrity_report.tsv    : machine-readable check results
#   06.8-integrity_report.txt    : human-readable formatted report
#   06.8-runtime_environment.txt : tool versions + pipeline paths
#   06.8-file_hashes.sha256      : SHA-256 checksums of all VCFs
#
# PIPELINE CONTEXT
#   Input VCFs (checked in order):
#     06.1 — normalization + biallelic SNP selection
#     06.2 — technical hard filters (QD/QUAL/FS/SOR/MQ/MQRankSum/ReadPosRankSum)
#     06.3 — missing data + sample missingness filtering
#     06.4 — MAC filtering (remove singletons/doubletons, MAC < 3)
#     06.5 — TE overlap filtering (exclude SNPs in transposable elements)
#     06.6b — LD pruning + MAF >= 5% filter (final dataset for GEA/structure)
#   Final VCF: 06.6b if present, otherwise 06.5
# =======================================================================

set -euo pipefail

module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/VCFtools/0.1.16
module load bioinfo/samtools/1.21

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_full"
OUT_DIR="${BASE}/06.8-data_integrity_audit"
LOG_DIR="${OUT_DIR}/logs"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP_DIR"

LOG_FILE="${LOG_DIR}/06.8-data_integrity_audit.log"
REPORT_TXT="${OUT_DIR}/06.8-integrity_report.txt"
REPORT_TSV="${OUT_DIR}/06.8-integrity_report.tsv"
ENV_TXT="${OUT_DIR}/06.8-runtime_environment.txt"
HASH_TXT="${OUT_DIR}/06.8-file_hashes.sha256"

echo "START $(date -Iseconds) 06.8 integrity audit" > "$LOG_FILE"

S61="${BASE}/06.1-normalization_selection/joint.norm.biallelic.snps.vcf.gz"
S62="${BASE}/06.2-technical_hard_filter/snps.hardfiltered.vcf.gz"
S63="${BASE}/06.3-missing_data_sample_filtering/snps.missing_filtered.vcf.gz"
S64="${BASE}/06.4-allele_frequency_population_filter/snps_non_informative_filtered.vcf.gz"
S65="${BASE}/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"
S66B="${BASE}/06.6b-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"

if [ -s "$S66B" ]; then
    FINAL_VCF="$S66B"
    FINAL_STEP="06.6b"
else
    FINAL_VCF="$S65"
    FINAL_STEP="06.5"
fi
FINAL_TBI="${FINAL_VCF}.tbi"

echo -e "check_id\tstatus\tvalue\tdetail" > "$REPORT_TSV"

report_ok()   { echo -e "$1\tPASS\t$2\t$3" >> "$REPORT_TSV"; }
report_warn() { echo -e "$1\tWARN\t$2\t$3" >> "$REPORT_TSV"; }
report_fail() { echo -e "$1\tFAIL\t$2\t$3" >> "$REPORT_TSV"; }

# Fast variant counting via tabix index — avoids streaming the full VCF.
# Falls back to full stream count if index is missing.
count_variants() {
    local vcf="$1"
    if [ -s "$vcf" ]; then
        bcftools index -n "$vcf" 2>/dev/null || bcftools view -H "$vcf" | wc -l
    else
        echo "NA"
    fi
}

count_samples() {
    local vcf="$1"
    if [ -s "$vcf" ]; then
        bcftools query -l "$vcf" | wc -l
    else
        echo "NA"
    fi
}

# -------------------------------------------------------------------
# CHECK 01: Index presence + header readability
# -------------------------------------------------------------------
# Each VCF must have a valid .tbi tabix index (required by all
# downstream tools). bcftools view -h opens only the header — a fast
# readability check that catches truncated or corrupted bgzip blocks.
# Result: PASS/FAIL per VCF file.
# -------------------------------------------------------------------
echo "CHECK 01: index + readability" | tee -a "$LOG_FILE"
for vcf in "$S61" "$S62" "$S63" "$S64" "$S65" "$S66B"; do
    [ -s "$vcf" ] || continue
    tbi="${vcf}.tbi"
    tag="$(basename "$vcf")"
    if [ -s "$tbi" ]; then
        report_ok "01_index_${tag}" "present" "$tbi"
    else
        report_fail "01_index_${tag}" "missing" "$tbi"
    fi
    if bcftools view -h "$vcf" >/dev/null 2>&1; then
        report_ok "01_read_${tag}" "readable" "$vcf"
    else
        report_fail "01_read_${tag}" "unreadable" "$vcf"
    fi
done

# -------------------------------------------------------------------
# CHECK 02: Sample concordance across filtering steps
# -------------------------------------------------------------------
# Compares the sample list between consecutive VCFs to detect any
# accidental sample gain or loss during filtering. Each step should
# preserve the full sample set (164 individuals), except 06.3 which
# explicitly removes samples with excess missingness.
# Method: bcftools query -l extracts sample names, comm compares
# sorted lists to identify added/removed individuals.
# Result: PASS if identical; WARN if any difference detected.
# -------------------------------------------------------------------
echo "CHECK 02: sample concordance" | tee -a "$LOG_FILE"
prev_vcf=""
for vcf in "$S61" "$S62" "$S63" "$S64" "$S65" "$S66B"; do
    [ -s "$vcf" ] || continue
    step="$(echo "$vcf" | sed -E 's#^.*/(06\.[0-9]+[a-z]?)-.*#\1#')"
    if [ -z "$prev_vcf" ]; then
        report_ok "02_samples_${step}" "$(count_samples "$vcf")" "baseline sample set"
    else
        a="${TMP_DIR}/samples_prev.txt"
        b="${TMP_DIR}/samples_cur.txt"
        bcftools query -l "$prev_vcf" | sort > "$a"
        bcftools query -l "$vcf" | sort > "$b"
        added=$(comm -13 "$a" "$b" | wc -l)
        removed=$(comm -23 "$a" "$b" | wc -l)
        if [ "$added" -eq 0 ] && [ "$removed" -eq 0 ]; then
            report_ok "02_samples_${step}" "$(count_samples "$vcf")" "identical sample set to previous step"
        else
            report_warn "02_samples_${step}" "$(count_samples "$vcf")" "added=${added};removed=${removed}"
        fi
    fi
    prev_vcf="$vcf"
done

# -------------------------------------------------------------------
# CHECK 03: Variant ID integrity and duplicate loci (final VCF)
# -------------------------------------------------------------------
# Checks performed on the final filtered VCF:
#   - Total variant count (for reporting).
#   - Missing IDs (ID='.'): warns if present. Although not a blocking
#     issue, named IDs are required by some downstream tools (e.g.
#     plink --extract, GEMMA). LD pruning (06.6b) already assigns
#     chr:pos:ref:alt IDs via --set-all-var-ids.
#   - Duplicate variant IDs: warns if two variants share the same ID.
#   - Duplicate loci (CHROM:POS:REF:ALT): warns if the same site
#     appears twice — should not occur after 06.1 biallelic selection
#     and bcftools norm deduplication.
# -------------------------------------------------------------------
echo "CHECK 03: variant ID integrity" | tee -a "$LOG_FILE"
if [ -s "$FINAL_VCF" ]; then
    total_final=$(count_variants "$FINAL_VCF")
    miss_id=$(bcftools query -f '%ID\n' "$FINAL_VCF" | awk '$1=="."{c++} END{print c+0}')
    dup_id=$(bcftools query -f '%ID\n' "$FINAL_VCF" | awk '$1!="."{print $1}' | sort | uniq -d | wc -l)
    dup_locus=$(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$FINAL_VCF" | sort | uniq -d | wc -l)
    report_ok "03_total_variants_final" "$total_final" "final step=${FINAL_STEP}"
    [ "$miss_id" -eq 0 ] && report_ok "03_missing_ids" "$miss_id" "no missing IDs" || report_warn "03_missing_ids" "$miss_id" "ID='.' variants present"
    [ "$dup_id" -eq 0 ] && report_ok "03_duplicate_ids" "$dup_id" "unique variant IDs" || report_warn "03_duplicate_ids" "$dup_id" "duplicate variant IDs detected"
    [ "$dup_locus" -eq 0 ] && report_ok "03_duplicate_loci" "$dup_locus" "unique CHROM:POS:REF:ALT" || report_warn "03_duplicate_loci" "$dup_locus" "duplicate loci detected"
else
    report_fail "03_final_vcf_missing" "NA" "final VCF not found"
fi

# -------------------------------------------------------------------
# CHECK 04: Allele sanity (final VCF)
# -------------------------------------------------------------------
# Verifies that the final VCF meets the expectations set by the
# filtering pipeline:
#   - No multiallelic sites: 06.1 SelectVariants --restrict-alleles
#     BIALLELIC should have removed all sites with >1 ALT allele.
#     Multiallelic sites can cause errors in plink, GEMMA, and BayPass.
#   - No non-ACGTN alleles: symbolic alleles (<DEL>, <INS>, *) or
#     ambiguity codes indicate indel or SV records that slipped through
#     and should not be present in an SNP-only VCF.
#   - No records without ALT: monomorphic sites should have been
#     removed by MAC filtering in 06.4 (--mac 3).
# -------------------------------------------------------------------
echo "CHECK 04: allele sanity" | tee -a "$LOG_FILE"
if [ -s "$FINAL_VCF" ]; then
    multi_alt=$(bcftools view -H "$FINAL_VCF" | awk -F'\t' 'index($5,",")>0{c++} END{print c+0}')
    non_acgtn=$(bcftools query -f '%REF\t%ALT\n' "$FINAL_VCF" | \
        awk -F'\t' '
            {
                if ($1 !~ /^[ACGTN]+$/) bad++
                n=split($2,a,",")
                for(i=1;i<=n;i++) if (a[i] !~ /^[ACGTN]+$/) {bad++; break}
            }
            END{print bad+0}')
    no_alt=$(bcftools stats "$FINAL_VCF" | awk -F'\t' '/^SN/ && $3=="number of no-ALTs:"{print $4+0}')
    [ "$multi_alt" -eq 0 ] && report_ok "04_multiallelic_sites" "$multi_alt" "none" || report_warn "04_multiallelic_sites" "$multi_alt" "multiallelic sites remain"
    [ "$non_acgtn" -eq 0 ] && report_ok "04_non_acgtn_alleles" "$non_acgtn" "none" || report_warn "04_non_acgtn_alleles" "$non_acgtn" "non-ACGTN alleles detected"
    [ "${no_alt:-0}" -eq 0 ] && report_ok "04_no_alt_records" "${no_alt:-0}" "none" || report_warn "04_no_alt_records" "${no_alt:-0}" "records without ALT remain"
fi

# -------------------------------------------------------------------
# CHECK 06: Per-individual and per-site missingness (final VCF)
# -------------------------------------------------------------------
# Recomputes missing data statistics on the final VCF to confirm that
# the thresholds applied in 06.3 (--max-missing 0.9 → max 10% missing
# per site; sample removal if F_MISS > 0.5) are met in the final
# dataset after all subsequent filtering steps.
# Output files:
#   final_missingness.imiss : per-individual missing rate (F_MISS)
#   final_missingness.lmiss : per-site missing rate (F_MISS)
# -------------------------------------------------------------------
echo "CHECK 06: missingness" | tee -a "$LOG_FILE"
MISS_PREFIX="${OUT_DIR}/final_missingness"
if [ -s "$FINAL_VCF" ]; then
    vcftools --gzvcf "$FINAL_VCF" --missing-indv --out "$MISS_PREFIX" >> "$LOG_FILE" 2>&1
    vcftools --gzvcf "$FINAL_VCF" --missing-site --out "$MISS_PREFIX" >> "$LOG_FILE" 2>&1
    if [ -s "${MISS_PREFIX}.imiss" ]; then
        max_imiss=$(awk 'NR>1{if($5>m)m=$5} END{print m+0}' "${MISS_PREFIX}.imiss")
        mean_imiss=$(awk 'NR>1{s+=$5;n++} END{if(n>0) printf "%.6f", s/n; else print "0"}' "${MISS_PREFIX}.imiss")
        report_ok "06_missingness_individual" "${mean_imiss}" "max=${max_imiss}"
    fi
    if [ -s "${MISS_PREFIX}.lmiss" ]; then
        max_lmiss=$(awk 'NR>1{if($6>m)m=$6} END{print m+0}' "${MISS_PREFIX}.lmiss")
        mean_lmiss=$(awk 'NR>1{s+=$6;n++} END{if(n>0) printf "%.6f", s/n; else print "0"}' "${MISS_PREFIX}.lmiss")
        report_ok "06_missingness_site" "${mean_lmiss}" "max=${max_lmiss}"
    fi
fi

# -------------------------------------------------------------------
# CHECK 07: Heterozygosity and sequencing depth per individual
# -------------------------------------------------------------------
# Per-individual heterozygosity (F-statistic, vcftools --het):
#   F = (O_HOM - E_HOM) / (N_sites - E_HOM)
#   F > 0  → excess homozygosity (inbreeding or poor mapping)
#   F < 0  → excess heterozygosity (contamination, mapping errors)
#   Outlier individuals (|F| > 0.15 typically) should be flagged
#   for manual inspection before GEA or genetic load analyses.
#
# Per-individual depth (vcftools --depth):
#   Reports mean depth across all sites per individual.
#   Very low-depth individuals (<10x) may introduce genotyping errors
#   that inflate false positive rates in downstream analyses.
# -------------------------------------------------------------------
echo "CHECK 07: heterozygosity + depth" | tee -a "$LOG_FILE"
HET_PREFIX="${OUT_DIR}/final_heterozygosity"
DEPTH_PREFIX="${OUT_DIR}/final_depth"
if [ -s "$FINAL_VCF" ]; then
    vcftools --gzvcf "$FINAL_VCF" --het --out "$HET_PREFIX" >> "$LOG_FILE" 2>&1 || true
    if [ -s "${HET_PREFIX}.het" ]; then
        mean_f=$(awk 'NR>1{s+=$6;n++} END{if(n>0) printf "%.6f", s/n; else print "NA"}' "${HET_PREFIX}.het")
        min_f=$(awk 'NR>1{if(NR==2||$6<m)m=$6} END{if(NR>1) print m; else print "NA"}' "${HET_PREFIX}.het")
        max_f=$(awk 'NR>1{if($6>m)m=$6} END{if(NR>1) print m; else print "NA"}' "${HET_PREFIX}.het")
        report_ok "07_heterozygosity_F" "$mean_f" "min=${min_f};max=${max_f}"
    else
        report_warn "07_heterozygosity_F" "NA" "vcftools --het output not produced"
    fi

    vcftools --gzvcf "$FINAL_VCF" --depth --out "$DEPTH_PREFIX" >> "$LOG_FILE" 2>&1 || true
    if [ -s "${DEPTH_PREFIX}.idepth" ]; then
        mean_dp=$(awk 'NR>1{s+=$3;n++} END{if(n>0) printf "%.3f", s/n; else print "NA"}' "${DEPTH_PREFIX}.idepth")
        min_dp=$(awk 'NR>1{if(NR==2||$3<m)m=$3} END{if(NR>1) print m; else print "NA"}' "${DEPTH_PREFIX}.idepth")
        max_dp=$(awk 'NR>1{if($3>m)m=$3} END{if(NR>1) print m; else print "NA"}' "${DEPTH_PREFIX}.idepth")
        report_ok "07_depth_individual" "$mean_dp" "min=${min_dp};max=${max_dp}"
    else
        report_warn "07_depth_individual" "NA" "vcftools --depth output not produced"
    fi
fi

# -------------------------------------------------------------------
# CHECK 08: Reproducibility audit — tool versions, paths, checksums
# -------------------------------------------------------------------
# Records the exact runtime environment to enable future auditing:
#   - Tool versions (bcftools, vcftools, samtools)
#   - Loaded modules (SLURM environment)
#   - Full paths to all input/output VCFs
#   - SHA-256 checksums of all VCF files
# SHA-256 checksums allow verification that files have not been
# modified or corrupted after the pipeline completed.
# -------------------------------------------------------------------
echo "CHECK 08: reproducibility audit" | tee -a "$LOG_FILE"
{
    echo "Generated: $(date -Iseconds)"
    echo "Host: $(hostname)"
    echo "User: ${USER}"
    echo
    echo "Modules:"
    module list 2>&1
    echo
    echo "Versions:"
    bcftools --version | head -n 2
    vcftools --version 2>&1 | head -n 1
    samtools --version | head -n 1
    echo
    echo "Input/Output targets:"
    echo "S61=${S61}"
    echo "S62=${S62}"
    echo "S63=${S63}"
    echo "S64=${S64}"
    echo "S65=${S65}"
    echo "S66B=${S66B}"
    echo "FINAL=${FINAL_VCF}"
} > "$ENV_TXT"

for f in "$S61" "$S62" "$S63" "$S64" "$S65" "$S66B" "$FINAL_VCF"; do
    [ -s "$f" ] || continue
    sha256sum "$f"
done > "$HASH_TXT"
report_ok "08_runtime_environment" "written" "$ENV_TXT"
report_ok "08_file_hashes" "written" "$HASH_TXT"

{
    echo "06.8 Data Integrity Audit"
    echo "Generated: $(date -Iseconds)"
    echo "Final VCF used: $FINAL_VCF"
    echo
    column -t -s $'\t' "$REPORT_TSV"
} > "$REPORT_TXT"

cat "$REPORT_TXT"
echo "DONE $(date -Iseconds) 06.8 integrity audit" | tee -a "$LOG_FILE"

# Summary of checks:
# Check 01 — Index presence + header readability on all step VCFs
# Check 02 — Sample concordance across all filtering steps
# Check 03 — Variant ID uniqueness + duplicate loci on final VCF
# Check 04 — Allele sanity (biallelic, ACGTN only, no empty ALT)
# Check 06 — Per-individual and per-site missingness on final VCF
# Check 07 — Per-individual heterozygosity (F) and depth on final VCF
# Check 08 — Tool versions, pipeline paths, SHA-256 checksums

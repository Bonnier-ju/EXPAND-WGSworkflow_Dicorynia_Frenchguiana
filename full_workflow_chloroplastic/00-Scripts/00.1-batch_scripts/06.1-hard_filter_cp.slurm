#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J cp-06.1-hard_filter
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load devel/java/17.0.6
module load bioinfo/GATK/4.2.6.1
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete.csv"
REF_GENOME="$(awk -F, 'NR==2{print $7}' "$SAMPLE_SHEET")"
INPUT_VCF="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/05-variant_calling/workflow_full/genotypegvcfs/jointgenotypes.raw.vcf.gz"
OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/06-haplotype_filtering/workflow_full/hard_filter"
TMP_DIR="${OUTPUT_DIR}/tmp"
LOG_FILE="${OUTPUT_DIR}/hard_filter_cp.log"

JAVA_MEM_GB=14

SNP_RAW="${OUTPUT_DIR}/snps.raw.vcf.gz"
INDEL_RAW="${OUTPUT_DIR}/indels.raw.vcf.gz"
SNP_FILTERED="${OUTPUT_DIR}/snps.filtered.vcf.gz"
INDEL_FILTERED="${OUTPUT_DIR}/indels.filtered.vcf.gz"
MERGED_VCF="${OUTPUT_DIR}/variants.filtered.merged.vcf.gz"
OUTPUT_VCF="${OUTPUT_DIR}/variants.hard_filtered.pass.vcf.gz"

mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found or empty: $INPUT_VCF" >&2; exit 1
fi
if [ -z "$REF_GENOME" ] || [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found: $REF_GENOME" >&2; exit 1
fi
if [ -s "$OUTPUT_VCF" ]; then
    echo "SKIP: output already exists: $OUTPUT_VCF"; exit 0
fi

echo "START $(date -Iseconds) hard_filter cp" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Separate SNPs and INDELs
# -------------------------------------------------------------------
echo "INFO: extracting SNPs" >> "$LOG_FILE"
gatk --java-options "-Xmx${JAVA_MEM_GB}G" SelectVariants \
    -R "$REF_GENOME" -V "$INPUT_VCF" \
    --select-type-to-include SNP \
    --tmp-dir "$TMP_DIR" \
    -O "$SNP_RAW" >> "$LOG_FILE" 2>&1

echo "INFO: extracting INDELs" >> "$LOG_FILE"
gatk --java-options "-Xmx${JAVA_MEM_GB}G" SelectVariants \
    -R "$REF_GENOME" -V "$INPUT_VCF" \
    --select-type-to-include INDEL \
    --tmp-dir "$TMP_DIR" \
    -O "$INDEL_RAW" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 3: Hard filter SNPs (seuils adaptés chloroplaste haploïde)
# -------------------------------------------------------------------
echo "INFO: hard filtering SNPs" >> "$LOG_FILE"
gatk --java-options "-Xmx${JAVA_MEM_GB}G" VariantFiltration \
    -R "$REF_GENOME" -V "$SNP_RAW" \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "FS > 60.0"             --filter-name "FS60" \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --tmp-dir "$TMP_DIR" \
    -O "$SNP_FILTERED" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 4: Hard filter INDELs
# -------------------------------------------------------------------
echo "INFO: hard filtering INDELs" >> "$LOG_FILE"
gatk --java-options "-Xmx${JAVA_MEM_GB}G" VariantFiltration \
    -R "$REF_GENOME" -V "$INDEL_RAW" \
    --filter-expression "QD < 2.0"               --filter-name "QD2" \
    --filter-expression "FS > 200.0"             --filter-name "FS200" \
    --filter-expression "SOR > 10.0"             --filter-name "SOR10" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --tmp-dir "$TMP_DIR" \
    -O "$INDEL_FILTERED" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 5: Merge SNPs + INDELs filtrés
# -------------------------------------------------------------------
echo "INFO: merging filtered SNPs and INDELs" >> "$LOG_FILE"
gatk --java-options "-Xmx${JAVA_MEM_GB}G" MergeVcfs \
    -I "$SNP_FILTERED" -I "$INDEL_FILTERED" \
    --TMP_DIR "$TMP_DIR" \
    -O "$MERGED_VCF" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 6: Garder uniquement les variants PASS
# -------------------------------------------------------------------
echo "INFO: selecting PASS variants only" >> "$LOG_FILE"
gatk --java-options "-Xmx${JAVA_MEM_GB}G" SelectVariants \
    -R "$REF_GENOME" -V "$MERGED_VCF" \
    --exclude-filtered \
    --tmp-dir "$TMP_DIR" \
    -O "$OUTPUT_VCF" >> "$LOG_FILE" 2>&1

if [ -s "$OUTPUT_VCF" ]; then
    n_snp=$(gatk --java-options "-Xmx4G" CountVariants -V "$SNP_FILTERED" --exclude-filtered 2>/dev/null | tail -1)
    n_indel=$(gatk --java-options "-Xmx4G" CountVariants -V "$INDEL_FILTERED" --exclude-filtered 2>/dev/null | tail -1)
    echo "DONE $(date -Iseconds) snps_pass=${n_snp} indels_pass=${n_indel} output=${OUTPUT_VCF}" | tee -a "$LOG_FILE"
else
    echo "ERROR: hard filter output missing" | tee -a "$LOG_FILE"; exit 1
fi

rm -rf "$TMP_DIR"

# -------------------------------------------------------------------
# Description
# Step 1: SelectVariants SNPs / INDELs depuis le VCF brut (05.3)
# Step 2: VariantFiltration SNPs — seuils GATK best practices (chloroplaste haploïde)
# Step 3: VariantFiltration INDELs
# Step 4: MergeVcfs SNPs + INDELs filtrés
# Step 5: SelectVariants --exclude-filtered → VCF final PASS uniquement
# -------------------------------------------------------------------

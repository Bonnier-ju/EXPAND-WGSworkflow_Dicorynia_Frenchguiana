#!/bin/bash
#SBATCH -p workq
#SBATCH --time=12:00:00
#SBATCH -J 06.5-te_filter
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load devel/python/Python-3.11.1
module load bioinfo/Bcftools/1.9
module load bioinfo/samtools/1.21
module load bioinfo/bedtools/2.30.0
module load bioinfo/BEDOPS/2.4.41


# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# -------------------------------------------------------------------
INPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.4-allele_frequency_population_filter"
INPUT_VCF="${INPUT_DIR}/snps_non_informative_filtered.vcf.gz"
INPUT_TBI="${INPUT_VCF}.tbi"

# Repeat/TE annotation (GFF/GFF.gz from reference annotation)
TE_GFF="${TE_GFF:-/home/jbonnier/work/chapitre_5/raw_data/reference_genome/nuclear/Dgu_HS1_HYBRID_SCAFFOLD.fa.out.gff.gz}"

OUTPUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/06-variant_filtering/workflow_test/06.5-transposable_element_filtering"
LOG_DIR="${OUTPUT_DIR}/logs"
STATS_DIR="${OUTPUT_DIR}/stats"
TMP_DIR="${OUTPUT_DIR}/tmp"
LOG_FILE="${LOG_DIR}/06.5-te_filtering.log"

TE_BED_RAW="${OUTPUT_DIR}/te_regions.raw.bed"
TE_BED="${OUTPUT_DIR}/te_regions.merged.bed"
OUTPUT_VCF="${OUTPUT_DIR}/snps_TE_filtered.vcf.gz"
OUTPUT_TBI="${OUTPUT_VCF}.tbi"

RAW_STATS="${STATS_DIR}/snps_non_informative_filtered.stats.txt"
TE_STATS="${STATS_DIR}/snps_TE_filtered.stats.txt"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$STATS_DIR" "$TMP_DIR"
echo "START $(date -Iseconds) 06.5 TE filtering" > "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ] || [ ! -s "$INPUT_TBI" ]; then
    echo "ERROR: input VCF missing: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$TE_GFF" ]; then
    echo "ERROR: TE GFF file not found: $TE_GFF" | tee -a "$LOG_FILE"
    exit 1
fi
if ! command -v gff2bed >/dev/null 2>&1; then
    echo "ERROR: gff2bed not found in environment (BEDOPS module)." | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 2: Build TE BED from GFF annotation
# -------------------------------------------------------------------
if [[ "$TE_GFF" =~ \.gz$ ]]; then
    zcat "$TE_GFF" | gff2bed | awk -F'\t' '{print $1"\t"$2"\t"$3}' > "$TE_BED_RAW"
else
    gff2bed < "$TE_GFF" | awk -F'\t' '{print $1"\t"$2"\t"$3}' > "$TE_BED_RAW"
fi

if [ ! -s "$TE_BED_RAW" ]; then
    echo "WARNING: no TE regions detected in GFF. Output VCF will match input." | tee -a "$LOG_FILE"
    cp "$INPUT_VCF" "$OUTPUT_VCF"
    cp "$INPUT_TBI" "$OUTPUT_TBI"
    bcftools stats "$INPUT_VCF" > "$RAW_STATS"
    cp "$RAW_STATS" "$TE_STATS"
    echo "DONE $(date -Iseconds) 06.5 completed (no TE region found)" | tee -a "$LOG_FILE"
    exit 0
fi

sort -k1,1 -k2,2n "$TE_BED_RAW" | bedtools merge -i - > "$TE_BED"
echo "INFO: TE intervals raw=$(wc -l < "$TE_BED_RAW") merged=$(wc -l < "$TE_BED")" | tee -a "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Remove TE-overlapping variants
# -------------------------------------------------------------------
bedtools intersect -header -v -a "$INPUT_VCF" -b "$TE_BED" | bgzip -c > "$OUTPUT_VCF"
bcftools index -t "$OUTPUT_VCF"

if [ ! -s "$OUTPUT_VCF" ] || [ ! -s "$OUTPUT_TBI" ]; then
    echo "ERROR: TE-filtered VCF output missing" | tee -a "$LOG_FILE"
    exit 1
fi

# -------------------------------------------------------------------
# STEP 4: Stats before/after TE filtering
# -------------------------------------------------------------------
bcftools stats "$INPUT_VCF" > "$RAW_STATS"
bcftools stats "$OUTPUT_VCF" > "$TE_STATS"

echo "DONE $(date -Iseconds) 06.5 completed" | tee -a "$LOG_FILE"
echo "Input VCF : $INPUT_VCF"
echo "TE GFF    : $TE_GFF"
echo "TE BED    : $TE_BED"
echo "Output VCF: $OUTPUT_VCF"
echo "Stats dir : $STATS_DIR"

# Step 1: Input VCF from 06.4 (snps_non_informative_filtered.vcf.gz)
# Step 2: Extract and merge TE regions from annotation GFF
# Step 3: Remove variants overlapping TE regions (bedtools intersect -v)
# Step 4: Produce bcftools stats before/after filtering

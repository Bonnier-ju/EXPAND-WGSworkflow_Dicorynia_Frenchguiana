#!/bin/bash
#SBATCH -p workq
#SBATCH --time=02:00:00
#SBATCH -J cp-06.4-haplotype_tables
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# To run on the test dataset, change WORKFLOW:
#   WORKFLOW="workflow_test_9ind"
# -------------------------------------------------------------------
WORKFLOW="workflow_full"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
INPUT_VCF="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/missing_data_filter/variants.filtered.vcf.gz"
OUTPUT_DIR="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/haplotype_tables"
LOG_FILE="${OUTPUT_DIR}/haplotype_tables_cp.log"

GENOTYPE_MATRIX="${OUTPUT_DIR}/genotype_matrix.tsv"
HAPLOTYPE_TABLE="${OUTPUT_DIR}/haplotype_table.tsv"
HAPLOTYPE_FREQ="${OUTPUT_DIR}/haplotype_frequencies.tsv"

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found or empty: $INPUT_VCF" >&2; exit 1
fi
if [ -s "$HAPLOTYPE_TABLE" ]; then
    echo "SKIP: output already exists: $HAPLOTYPE_TABLE"; exit 0
fi

echo "START $(date -Iseconds) haplotype_tables cp" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Extract genotype matrix
#   Columns: CHROM, POS, REF, ALT, then one GT per sample (0 or 1)
#   Haploid chloroplast: GT = 0 (ref) or 1 (alt), never heterozygous
# -------------------------------------------------------------------
echo "INFO: extracting genotype matrix from VCF" >> "$LOG_FILE"

# Retrieve header with sample names
SAMPLES=$(bcftools query -l "$INPUT_VCF" | tr '\n' '\t' | sed 's/\t$//')
echo -e "CHROM\tPOS\tREF\tALT\t${SAMPLES}" > "$GENOTYPE_MATRIX"

# Extract genotypes (haploid: 0 or 1, replace missing . with NA)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$INPUT_VCF" 2>> "$LOG_FILE" | \
    sed 's/\.\//NA/g; s/\//NA/g; s/|/\t/g' >> "$GENOTYPE_MATRIX"

n_sites=$(tail -n +2 "$GENOTYPE_MATRIX" | wc -l)
echo "INFO: genotype matrix: ${n_sites} sites" >> "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Build the haplotype table
#   Haplotype = concatenation of alleles across all sites for one sample
#   Assign a unique ID (H1, H2, ...) to each distinct haplotype
# -------------------------------------------------------------------
echo "INFO: building haplotype table" >> "$LOG_FILE"

awk 'BEGIN{OFS="\t"}
NR==1 {
    # Get sample names (columns 5+)
    for(i=5; i<=NF; i++) samples[i]=$i
    next
}
{
    # Build haplotype string for each sample
    for(i=5; i<=NF; i++) {
        haplo[i] = haplo[i] $i
    }
}
END {
    # Identify unique haplotypes and assign IDs
    h_id = 0
    for(i=5; i<=length(samples)+4; i++) {
        h = haplo[i]
        if(!(h in haplo_id)) {
            h_id++
            haplo_id[h] = "H" h_id
        }
    }
    # Write haplotype table
    print "sample_id", "haplotype_id", "haplotype_string"
    for(i=5; i<=length(samples)+4; i++) {
        h = haplo[i]
        print samples[i], haplo_id[h], h
    }
}' "$GENOTYPE_MATRIX" > "$HAPLOTYPE_TABLE"

# -------------------------------------------------------------------
# STEP 4: Compute haplotype frequencies
# -------------------------------------------------------------------
echo "INFO: computing haplotype frequencies" >> "$LOG_FILE"

n_total=$(tail -n +2 "$HAPLOTYPE_TABLE" | wc -l)

awk -v total="$n_total" 'BEGIN{OFS="\t"}
NR==1 {next}
{
    count[$2]++
    # Keep one representative string per haplotype
    if(!($2 in haplo_str)) haplo_str[$2]=$3
}
END {
    print "haplotype_id", "count", "frequency", "haplotype_string"
    for(h in count) {
        printf "%s\t%d\t%.4f\t%s\n", h, count[h], count[h]/total, haplo_str[h]
    }
}' "$HAPLOTYPE_TABLE" | sort -t$'\t' -k1,1V >> "$HAPLOTYPE_FREQ"

# Add header to frequency file
head -1 <(echo -e "haplotype_id\tcount\tfrequency\thaplotype_string") | cat - "$HAPLOTYPE_FREQ" > "${HAPLOTYPE_FREQ}.tmp"
mv "${HAPLOTYPE_FREQ}.tmp" "$HAPLOTYPE_FREQ"

if [ -s "$HAPLOTYPE_TABLE" ] && [ -s "$HAPLOTYPE_FREQ" ]; then
    n_haplo=$(tail -n +2 "$HAPLOTYPE_FREQ" | wc -l)
    echo "DONE $(date -Iseconds) sites=${n_sites} haplotypes_uniques=${n_haplo} output=${OUTPUT_DIR}" | tee -a "$LOG_FILE"
else
    echo "ERROR: haplotype table output missing" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# Description
# Step 1: bcftools query → genotype matrix (sites × samples)
# Step 2: Concatenate alleles per sample → haplotype string
# Step 3: Assign unique IDs (H1, H2, ...) to each distinct haplotype
# Step 4: Compute haplotype frequencies
# Outputs:
#   genotype_matrix.tsv      : full matrix CHROM/POS/REF/ALT + GTs
#   haplotype_table.tsv      : sample_id | haplotype_id | haplotype_string
#   haplotype_frequencies.tsv: haplotype_id | count | frequency | string
# Haploid chloroplast: GT = 0/1 only, no heterozygous genotypes expected
# -------------------------------------------------------------------

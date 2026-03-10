#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J cp-06.5-consensus_sequences
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# -------------------------------------------------------------------
# STEP 0: Environment setup
# -------------------------------------------------------------------
module purge
module load bioinfo/Bcftools/1.9
module load bioinfo/samtools/1.21

# -------------------------------------------------------------------
# STEP 1: Define inputs/outputs
# To run on the test dataset, change WORKFLOW and SAMPLE_SHEET:
#   WORKFLOW="workflow_test_9ind"
#   SAMPLE_SHEET="...sample_sheet_chloroplastic_test.csv"
# -------------------------------------------------------------------
WORKFLOW="workflow_full"
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete.csv"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
REF_GENOME="$(awk -F, 'NR==2{print $7}' "$SAMPLE_SHEET")"
INPUT_VCF="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/missing_data_filter/variants.filtered.vcf.gz"
OUTPUT_DIR="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/consensus_sequences"
LOG_FILE="${OUTPUT_DIR}/consensus_sequences_cp.log"

PER_SAMPLE_DIR="${OUTPUT_DIR}/per_sample"
MULTIFASTA="${OUTPUT_DIR}/all_samples_consensus.fasta"
SNP_VCF="${OUTPUT_DIR}/variants.snps_only.vcf.gz"

mkdir -p "$OUTPUT_DIR" "$PER_SAMPLE_DIR"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found or empty: $INPUT_VCF" >&2; exit 1
fi
if [ -z "$REF_GENOME" ] || [ ! -s "$REF_GENOME" ]; then
    echo "ERROR: reference genome not found: $REF_GENOME" >&2; exit 1
fi
if [ ! -s "${REF_GENOME}.fai" ]; then
    echo "INFO: creating reference fai index" >> "$LOG_FILE"
    samtools faidx "$REF_GENOME"
fi
if [ -s "$MULTIFASTA" ]; then
    echo "SKIP: output already exists: $MULTIFASTA"; exit 0
fi

echo "START $(date -Iseconds) consensus_sequences cp" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Extract SNPs only for consensus generation
#   INDELs are excluded because they alter sequence length, which
#   breaks coordinate-alignment across samples (trimAl requires
#   all sequences to be the same length).
#   SNP substitutions preserve the reference length (160,511 bp),
#   guaranteeing true pseudo-alignment.
# -------------------------------------------------------------------
echo "INFO: extracting SNP-only VCF for consensus generation" >> "$LOG_FILE"
bcftools view -v snps -O z -o "$SNP_VCF" "$INPUT_VCF" >> "$LOG_FILE" 2>&1
bcftools index --tbi "$SNP_VCF" >> "$LOG_FILE" 2>&1
n_snps=$(bcftools stats "$SNP_VCF" 2>/dev/null | grep "^SN" | grep "number of SNPs" | awk '{print $NF}')
echo "INFO: SNP-only VCF: ${n_snps} SNPs retained for consensus" >> "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 3: Generate one consensus FASTA sequence per sample
#   -s SAMPLE  : select the sample
#   -H 1       : use the first (and only) haplotype (haploid chloroplast)
#   -f REF     : apply SNP variants onto the reference genome
#   Missing sites (./.) : retained as N in the consensus
# -------------------------------------------------------------------
echo "INFO: generating per-sample consensus sequences" >> "$LOG_FILE"

> "$MULTIFASTA"
done_count=0
fail_count=0

while IFS= read -r sample; do
    [ -z "$sample" ] && continue

    OUT_FA="${PER_SAMPLE_DIR}/${sample}.consensus.fasta"

    if [ -s "$OUT_FA" ]; then
        echo "INFO: SKIP ${sample} (already exists)" >> "$LOG_FILE"
        cat "$OUT_FA" >> "$MULTIFASTA"
        done_count=$((done_count + 1))
        continue
    fi

    if bcftools consensus \
            -s "$sample" \
            -H 1 \
            -f "$REF_GENOME" \
            "$SNP_VCF" 2>> "$LOG_FILE" | \
        sed "s/^>.*$/>${sample}/" > "$OUT_FA"; then

        if [ -s "$OUT_FA" ]; then
            cat "$OUT_FA" >> "$MULTIFASTA"
            echo "DONE sample=${sample}" >> "$LOG_FILE"
            done_count=$((done_count + 1))
        else
            echo "ERROR: empty consensus for ${sample}" | tee -a "$LOG_FILE"
            fail_count=$((fail_count + 1))
        fi
    else
        echo "ERROR: bcftools consensus failed for ${sample}" | tee -a "$LOG_FILE"
        fail_count=$((fail_count + 1))
    fi

done < <(bcftools query -l "$SNP_VCF")

# -------------------------------------------------------------------
# STEP 4: Validate the final multi-FASTA
# -------------------------------------------------------------------
if [ -s "$MULTIFASTA" ]; then
    n_seq=$(grep -c "^>" "$MULTIFASTA")
    echo "DONE $(date -Iseconds) sequences=${n_seq} success=${done_count} fail=${fail_count} output=${MULTIFASTA}" | tee -a "$LOG_FILE"
    if [ "$fail_count" -gt 0 ]; then
        echo "WARNING: ${fail_count} samples failed — check log for details" | tee -a "$LOG_FILE"
    fi
else
    echo "ERROR: multi-FASTA output is empty" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# Description
# Step 1: bcftools view -v snps → SNP-only VCF
#         INDELs excluded: applying indels per sample yields sequences
#         of different lengths, breaking coordinate-alignment required
#         by trimAl (06.6). SNPs preserve reference length (160,511 bp).
# Step 2: For each sample in the SNP-only VCF:
#   bcftools consensus -s SAMPLE -H 1 -f REF
#   → apply SNP alleles onto the reference sequence
#   → missing sites retained as N
#   → FASTA header renamed to sample_id
# Step 3: Concatenate all individual FASTAs into a multi-FASTA
# Outputs:
#   variants.snps_only.vcf.gz           : SNP-only VCF (intermediate)
#   per_sample/{sample}.consensus.fasta : one sequence per sample
#   all_samples_consensus.fasta         : multi-FASTA ready for trimAl (06.6)
# Haploid chloroplast: -H 1 (single haplotype, no phasing required)
# -------------------------------------------------------------------

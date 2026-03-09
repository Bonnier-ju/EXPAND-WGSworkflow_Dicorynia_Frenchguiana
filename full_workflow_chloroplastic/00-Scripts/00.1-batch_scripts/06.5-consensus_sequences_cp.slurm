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
# Pour tester sur le jeu de données test, changer WORKFLOW et SAMPLE_SHEET :
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
# STEP 2: Générer une séquence consensus FASTA par individu
#   -s SAMPLE  : sélectionner l'individu
#   -H 1       : utiliser le premier (et unique) haplotype (chloroplaste haploïde)
#   -f REF     : appliquer les variants sur le génome de référence
#   Sites manquants (./.) : conservés comme N dans le consensus
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
        ((done_count++))
        continue
    fi

    if bcftools consensus \
            -s "$sample" \
            -H 1 \
            -f "$REF_GENOME" \
            "$INPUT_VCF" 2>> "$LOG_FILE" | \
        sed "s/^>.*$/>${sample}/" > "$OUT_FA"; then

        if [ -s "$OUT_FA" ]; then
            cat "$OUT_FA" >> "$MULTIFASTA"
            echo "DONE sample=${sample}" >> "$LOG_FILE"
            ((done_count++))
        else
            echo "ERROR: empty consensus for ${sample}" | tee -a "$LOG_FILE"
            ((fail_count++))
        fi
    else
        echo "ERROR: bcftools consensus failed for ${sample}" | tee -a "$LOG_FILE"
        ((fail_count++))
    fi

done < <(bcftools query -l "$INPUT_VCF")

# -------------------------------------------------------------------
# STEP 3: Vérification du multi-FASTA final
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
# Step 1: Pour chaque individu du VCF filtré (06.3) :
#   bcftools consensus -s SAMPLE -H 1 -f REF
#   → applique les allèles de l'individu sur la séquence de référence
#   → sites manquants conservés comme N
#   → en-tête FASTA renommé au sample_id
# Step 2: Concaténation de tous les FASTA individuels en un multi-FASTA
# Sorties :
#   per_sample/{sample}.consensus.fasta : une séquence par individu
#   all_samples_consensus.fasta         : multi-FASTA prêt pour alignement phylo (07.1)
# Chloroplaste haploïde : -H 1 (un seul haplotype, pas de phase à résoudre)
# -------------------------------------------------------------------

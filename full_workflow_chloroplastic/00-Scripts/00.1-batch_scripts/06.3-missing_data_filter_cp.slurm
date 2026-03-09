#!/bin/bash
#SBATCH -p workq
#SBATCH --time=02:00:00
#SBATCH -J cp-06.3-missing_data_filter
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
# Pour tester sur le jeu de données test, changer WORKFLOW et SAMPLE_SHEET :
#   WORKFLOW="workflow_test_9ind"
#   SAMPLE_SHEET="...sample_sheet_chloroplastic_test.csv"
# Seuils (chloroplaste haploïde) :
#   MAX_SITE_MISSING  : taux de données manquantes max par site (0.2 = 20%)
#   MAX_SAMPLE_MISSING: taux de données manquantes max par individu (0.5 = 50%)
#   MIN_MAC           : nombre minimal d'allèles mineurs (1 = au moins 1 copie)
# -------------------------------------------------------------------
WORKFLOW="workflow_full"
SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete.csv"

BASE_DIR="/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic"
INPUT_VCF="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/normalization/variants.normalized.vcf.gz"
OUTPUT_DIR="${BASE_DIR}/06-haplotype_filtering/${WORKFLOW}/missing_data_filter"
LOG_FILE="${OUTPUT_DIR}/missing_data_filter_cp.log"

MAX_SITE_MISSING=0.2
MAX_SAMPLE_MISSING=0.5
MIN_MAC=1

SAMPLES_TO_EXCLUDE="${OUTPUT_DIR}/samples_high_missing.txt"
OUTPUT_VCF="${OUTPUT_DIR}/variants.filtered.vcf.gz"

mkdir -p "$OUTPUT_DIR"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: input VCF not found or empty: $INPUT_VCF" >&2; exit 1
fi
if [ -s "$OUTPUT_VCF" ]; then
    echo "SKIP: output already exists: $OUTPUT_VCF"; exit 0
fi

echo "START $(date -Iseconds) missing_data_filter cp" > "$LOG_FILE"

# -------------------------------------------------------------------
# STEP 2: Identifier les individus avec trop de données manquantes
#   bcftools stats -s- : stats par échantillon
#   Extraire le taux de génotypes manquants par individu
# -------------------------------------------------------------------
echo "INFO: computing per-sample missing rate" >> "$LOG_FILE"
bcftools stats -s - "$INPUT_VCF" 2>> "$LOG_FILE" | \
    awk -v thresh="$MAX_SAMPLE_MISSING" '
        /^PSC/ {
            sample=$3; total=$4+$5+$6; missing=$6
            if (total > 0 && (missing/total) > thresh) print sample
        }
    ' > "$SAMPLES_TO_EXCLUDE"

n_excluded=$(wc -l < "$SAMPLES_TO_EXCLUDE")
echo "INFO: ${n_excluded} samples excluded (missing > ${MAX_SAMPLE_MISSING})" >> "$LOG_FILE"
if [ "$n_excluded" -gt 0 ]; then
    echo "INFO: excluded samples:" >> "$LOG_FILE"
    cat "$SAMPLES_TO_EXCLUDE" >> "$LOG_FILE"
fi

# -------------------------------------------------------------------
# STEP 3: Filtrer sites et individus
#   -e 'F_MISSING > MAX_SITE_MISSING' : exclure sites trop incomplets
#   -c MIN_MAC:minor                  : garder sites avec ≥1 allèle mineur
#   --samples-file avec ^ prefix      : exclure les individus à exclure
# -------------------------------------------------------------------
echo "INFO: filtering sites (missing > ${MAX_SITE_MISSING}) and MAC < ${MIN_MAC}" >> "$LOG_FILE"

EXCLUDE_OPT=""
if [ "$n_excluded" -gt 0 ]; then
    EXCLUDE_OPT="--samples-file ^${SAMPLES_TO_EXCLUDE}"
fi

bcftools view \
    -e "F_MISSING > ${MAX_SITE_MISSING}" \
    -c "${MIN_MAC}:minor" \
    ${EXCLUDE_OPT} \
    -O z -o "$OUTPUT_VCF" \
    "$INPUT_VCF" >> "$LOG_FILE" 2>&1

# -------------------------------------------------------------------
# STEP 4: Indexer le VCF filtré
# -------------------------------------------------------------------
echo "INFO: indexing filtered VCF" >> "$LOG_FILE"
bcftools index --tbi "$OUTPUT_VCF" >> "$LOG_FILE" 2>&1

if [ -s "$OUTPUT_VCF" ]; then
    n_var=$(bcftools stats "$OUTPUT_VCF" 2>/dev/null | grep "^SN" | grep "number of records" | awk '{print $NF}')
    n_samples=$(bcftools query -l "$OUTPUT_VCF" 2>/dev/null | wc -l)
    echo "DONE $(date -Iseconds) variants=${n_var} samples=${n_samples} output=${OUTPUT_VCF}" | tee -a "$LOG_FILE"
else
    echo "ERROR: missing data filter output missing" | tee -a "$LOG_FILE"; exit 1
fi

# -------------------------------------------------------------------
# Description
# Step 1: bcftools stats -s- → taux de génotypes manquants par individu
# Step 2: Exclusion individus avec > MAX_SAMPLE_MISSING (défaut 50%)
# Step 3: bcftools view -e 'F_MISSING > 0.2' → filtre sites incomplets
#         -c 1:minor → garde sites polymorphes (MAC ≥ 1)
# Step 4: Indexation TBI
# Chloroplaste haploïde : pas de filtre HWE, pas de filtre MAF minimum strict
# -------------------------------------------------------------------

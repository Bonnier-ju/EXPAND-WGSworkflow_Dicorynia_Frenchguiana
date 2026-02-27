#!/bin/bash
#SBATCH -p workq
#SBATCH --time=04:00:00
#SBATCH -J 07.3-pop_inputs
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load bioinfo/Bcftools/1.9

BASE="/home/jbonnier/work/chapitre_5/full_workflow_nuclear"
VCF_66="${BASE}/06-variant_filtering/workflow_test/06.6-ld_pruning/snps_TE_filtered.ld_pruned.vcf.gz"
VCF_65="${BASE}/06-variant_filtering/workflow_test/06.5-transposable_element_filtering/snps_TE_filtered.vcf.gz"

SAMPLE_SHEET="${SAMPLE_SHEET:-${BASE}/metadata/sample_sheet_test.csv}"
POP_COLUMN="${POP_COLUMN:-site}"   # ex: site, project

OUT_DIR="${BASE}/07-population_structure/workflow_test/07.3-population_inputs"
LOG_DIR="${OUT_DIR}/logs"
SAMPLES_DIR="${OUT_DIR}/population_samples"
mkdir -p "$OUT_DIR" "$LOG_DIR" "$SAMPLES_DIR"

LOG_FILE="${LOG_DIR}/07.3-population_inputs.log"
SUMMARY_FILE="${OUT_DIR}/07.3-summary.txt"

if [ -s "$VCF_66" ]; then
    INPUT_VCF="$VCF_66"
    INPUT_STEP="06.6"
else
    INPUT_VCF="$VCF_65"
    INPUT_STEP="06.5"
fi

echo "START $(date -Iseconds) 07.3 population input preparation" > "$LOG_FILE"
echo "Input VCF: $INPUT_VCF (from $INPUT_STEP)" | tee -a "$LOG_FILE"
echo "Sample sheet: $SAMPLE_SHEET (POP_COLUMN=$POP_COLUMN)" | tee -a "$LOG_FILE"

if [ ! -s "$INPUT_VCF" ]; then
    echo "ERROR: missing input VCF: $INPUT_VCF" | tee -a "$LOG_FILE"
    exit 1
fi
if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: missing sample sheet: $SAMPLE_SHEET" | tee -a "$LOG_FILE"
    exit 1
fi

VCF_SAMPLES="${OUT_DIR}/vcf_samples.list"
MAP_RAW="${OUT_DIR}/sample_population_raw.tsv"
MAP_FINAL="${OUT_DIR}/sample_population_map.tsv"
KEEP_FILE="${OUT_DIR}/keep.samples.txt"
CLUSTER_FILE="${OUT_DIR}/plink.cluster.txt"
POPS_LIST="${OUT_DIR}/populations.list"

bcftools query -l "$INPUT_VCF" | sort -u > "$VCF_SAMPLES"

awk -F',' -v popcol="$POP_COLUMN" '
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="sample_id") sid=i
    if ($i==popcol) pop=i
  }
  if (!sid || !pop) {
    print "ERROR: sample_id or " popcol " not found in header" > "/dev/stderr"
    exit 2
  }
  next
}
{
  gsub(/\r/, "", $sid)
  gsub(/\r/, "", $pop)
  if ($sid != "") {
    if ($pop == "" || $pop == "NA" || $pop == "na") $pop = "UNKNOWN"
    print $sid "\t" $pop
  }
}
' "$SAMPLE_SHEET" > "$MAP_RAW"

join -t $'\t' -1 1 -2 1 \
  <(sort -k1,1 "$MAP_RAW") \
  <(awk '{print $1"\t"$1}' "$VCF_SAMPLES" | sort -k1,1) \
  | awk -F'\t' '{print $1"\t"$2}' > "$MAP_FINAL"

if [ ! -s "$MAP_FINAL" ]; then
    echo "ERROR: no overlap between sample sheet and VCF samples" | tee -a "$LOG_FILE"
    exit 1
fi

awk -F'\t' '{print $1"\t"$1}' "$MAP_FINAL" > "$KEEP_FILE"
awk -F'\t' '{print $1"\t"$1"\t"$2}' "$MAP_FINAL" > "$CLUSTER_FILE"
awk -F'\t' '{print $2}' "$MAP_FINAL" | sort -u > "$POPS_LIST"

while read -r pop; do
    [ -z "$pop" ] && continue
    safe=$(echo "$pop" | sed 's/[^A-Za-z0-9._-]/_/g')
    awk -F'\t' -v p="$pop" '$2==p{print $1}' "$MAP_FINAL" > "${SAMPLES_DIR}/${safe}.samples.txt"
    awk -F'\t' -v p="$pop" '$2==p{print $1"\t"$1}' "$MAP_FINAL" > "${SAMPLES_DIR}/${safe}.keep.txt"
done < "$POPS_LIST"

N_SAMPLES=$(wc -l < "$MAP_FINAL")
N_POPS=$(wc -l < "$POPS_LIST")
{
    echo "07.3 Population input preparation summary"
    echo "Date: $(date -Iseconds)"
    echo "Input step: $INPUT_STEP"
    echo "Input VCF: $INPUT_VCF"
    echo "Sample sheet: $SAMPLE_SHEET"
    echo "Population column: $POP_COLUMN"
    echo "Samples retained: $N_SAMPLES"
    echo "Populations retained: $N_POPS"
    echo
    echo "Outputs:"
    echo "  - $MAP_FINAL"
    echo "  - $KEEP_FILE"
    echo "  - $CLUSTER_FILE"
    echo "  - $POPS_LIST"
    echo "  - ${SAMPLES_DIR}/*.samples.txt and *.keep.txt"
} > "$SUMMARY_FILE"

echo "DONE $(date -Iseconds) 07.3 completed" | tee -a "$LOG_FILE"
cat "$SUMMARY_FILE"

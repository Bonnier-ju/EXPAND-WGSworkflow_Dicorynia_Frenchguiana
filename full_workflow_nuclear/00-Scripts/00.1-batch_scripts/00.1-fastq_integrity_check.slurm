#!/bin/bash
#SBATCH -p workq
#SBATCH --time=72:00:00
#SBATCH -J 00-fastq_integrity_check
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.out
#SBATCH --array=1-8
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u

# Preflight checks for paired-end FASTQ files from a sample sheet.
# Expected CSV format: sample_id in column 1, read1/read2 in last 2 columns.

# Defaults allow running directly with sbatch without extra arguments.
DEFAULT_SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_invalid_trimmed_73.csv"
DEFAULT_OUT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/01-quality_control/workflow_retrim_73/fastq_integrity_report"
SAMPLES_PER_TASK=10

if [ "$#" -gt 2 ]; then
    echo "Usage: $0 [sample_sheet.csv] [output_dir]" >&2
    exit 1
fi

SAMPLE_SHEET="${1:-$DEFAULT_SAMPLE_SHEET}"
OUT_DIR="${2:-$DEFAULT_OUT_DIR}"
mkdir -p "$OUT_DIR"

TS="$(date +%Y%m%d_%H%M%S)"
TASK_TAG="${SLURM_ARRAY_TASK_ID:-1}"
REPORT="${OUT_DIR}/preflight_report_${TS}.tsv"
SUMMARY="${OUT_DIR}/preflight_summary_${TS}.txt"
REPORT="${OUT_DIR}/preflight_report_${TS}_task${TASK_TAG}.tsv"
SUMMARY="${OUT_DIR}/preflight_summary_${TS}_task${TASK_TAG}.txt"

if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found: $SAMPLE_SHEET" >&2
    exit 1
fi

echo -e "sample_id\tread1\tread2\texists\tgzip_ok\tpaired_count_ok\tpaired_id_ok\tstatus\tnotes" > "$REPORT"

TOTAL=0
OK_COUNT=0
FAIL_COUNT=0

check_pair_names() {
    # Compare first 20k read headers after normalizing /1 and /2 suffixes.
    # Return 0 if headers are synchronized, 1 otherwise.
    local r1="$1"
    local r2="$2"
    local t1 t2
    t1="$(mktemp)"
    t2="$(mktemp)"

    zcat "$r1" | awk 'NR%4==1 {h=$1; sub(/^@/,"",h); sub(/\/[12]$/,"",h); print h; if (++n==20000) exit}' > "$t1"
    zcat "$r2" | awk 'NR%4==1 {h=$1; sub(/^@/,"",h); sub(/\/[12]$/,"",h); print h; if (++n==20000) exit}' > "$t2"

    if paste "$t1" "$t2" | awk '$1 != $2 {bad=1; exit} END {exit bad ? 1 : 0}'; then
        rm -f "$t1" "$t2"
        return 0
    else
        rm -f "$t1" "$t2"
        return 1
    fi
}

TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1 "\t" $(NF-1) "\t" $NF}' | sed -n "${START_LINE},${END_LINE}p" | while IFS=$'\t' read -r sample_id read1 read2; do
    TOTAL=$((TOTAL + 1))
    notes=""

    exists="OK"
    gzip_ok="OK"
    paired_count_ok="OK"
    paired_id_ok="OK"
    status="PASS"

    if [ -z "$sample_id" ]; then
        status="FAIL"
        exists="FAIL"
        gzip_ok="NA"
        paired_count_ok="NA"
        paired_id_ok="NA"
        notes="empty_sample_id"
        echo -e "${sample_id}\t${read1}\t${read2}\t${exists}\t${gzip_ok}\t${paired_count_ok}\t${paired_id_ok}\t${status}\t${notes}" >> "$REPORT"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    if [ ! -s "$read1" ] || [ ! -s "$read2" ]; then
        status="FAIL"
        exists="FAIL"
        gzip_ok="NA"
        paired_count_ok="NA"
        paired_id_ok="NA"
        notes="missing_or_empty_fastq"
        echo -e "${sample_id}\t${read1}\t${read2}\t${exists}\t${gzip_ok}\t${paired_count_ok}\t${paired_id_ok}\t${status}\t${notes}" >> "$REPORT"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    if ! gzip -t "$read1" 2>/dev/null || ! gzip -t "$read2" 2>/dev/null; then
        status="FAIL"
        gzip_ok="FAIL"
        paired_count_ok="NA"
        paired_id_ok="NA"
        notes="gzip_corrupted_or_truncated"
        echo -e "${sample_id}\t${read1}\t${read2}\t${exists}\t${gzip_ok}\t${paired_count_ok}\t${paired_id_ok}\t${status}\t${notes}" >> "$REPORT"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    c1="$(zcat "$read1" | awk 'END {print NR/4}')"
    c2="$(zcat "$read2" | awk 'END {print NR/4}')"
    if [ "$c1" != "$c2" ]; then
        status="FAIL"
        paired_count_ok="FAIL"
        notes="read_count_mismatch:${c1}_vs_${c2}"
    fi

    if ! check_pair_names "$read1" "$read2"; then
        status="FAIL"
        paired_id_ok="FAIL"
        if [ -n "$notes" ]; then
            notes="${notes};pair_header_mismatch"
        else
            notes="pair_header_mismatch"
        fi
    fi

    if [ "$status" = "PASS" ]; then
        OK_COUNT=$((OK_COUNT + 1))
    else
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi

    echo -e "${sample_id}\t${read1}\t${read2}\t${exists}\t${gzip_ok}\t${paired_count_ok}\t${paired_id_ok}\t${status}\t${notes}" >> "$REPORT"
done

# Recompute totals from report to avoid subshell side effects from the while-loop pipe.
TOTAL="$(tail -n +2 "$REPORT" | wc -l | awk '{print $1}')"
OK_COUNT="$(awk -F'\t' 'NR>1 && $8=="PASS" {n++} END {print n+0}' "$REPORT")"
FAIL_COUNT="$(awk -F'\t' 'NR>1 && $8=="FAIL" {n++} END {print n+0}' "$REPORT")"

{
    echo "Preflight FASTQ validation"
    echo "Sample sheet: $SAMPLE_SHEET"
    echo "Report: $REPORT"
    echo "Total samples: $TOTAL"
    echo "PASS: $OK_COUNT"
    echo "FAIL: $FAIL_COUNT"
    echo ""
    echo "Notes:"
    echo "- This script validates file integrity and PE coherence."
    echo "- Contamination checks require dedicated tools (FastQC/MultiQC, Kraken2, FastQ Screen)."
} > "$SUMMARY"

cat "$SUMMARY"

if [ "$FAIL_COUNT" -gt 0 ]; then
    exit 2
fi

exit 0

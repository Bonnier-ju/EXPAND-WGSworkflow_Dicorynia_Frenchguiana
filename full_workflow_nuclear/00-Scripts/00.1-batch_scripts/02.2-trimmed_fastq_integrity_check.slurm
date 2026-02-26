#!/bin/bash
#SBATCH -p workq
#SBATCH --time=24:00:00
#SBATCH -J 02.2-trimmed_fastq_check
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.err
#SBATCH --mem=12G
#SBATCH --cpus-per-task=3
#SBATCH --array=1-18
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u
set -o pipefail

SAMPLE_SHEET="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_complete.csv"
TRIMMED_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/02-filtering_indexing/workflow_full/cleaned_reads"
REPORT_DIR="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/02-filtering_indexing/workflow_full/trimmed_integrity_report"
SAMPLES_PER_TASK=10
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

mkdir -p "$REPORT_DIR"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found or empty: $SAMPLE_SHEET" >&2
    exit 1
fi

# Columns used from sample sheet:
# sample_id (1), raw_read1 (last-1), raw_read2 (last)
ALL_FILES_LIST="${REPORT_DIR}/trimmed_check_samples.all.list"
TASK_FILES_LIST="${REPORT_DIR}/trimmed_check_samples.task${TASK_ID}.list"
tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1 "\t" $(NF-1) "\t" $NF}' > "$ALL_FILES_LIST"
sed -n "${START_LINE},${END_LINE}p" "$ALL_FILES_LIST" > "$TASK_FILES_LIST"

REPORT_TSV="${REPORT_DIR}/trimmed_integrity_report_task${TASK_ID}.tsv"
SUMMARY_TXT="${REPORT_DIR}/trimmed_integrity_summary_task${TASK_ID}.txt"

DUP_SAMPLE_IDS="${REPORT_DIR}/duplicate_sample_ids_task${TASK_ID}.txt"
DUP_READ_PAIRS="${REPORT_DIR}/duplicate_read_pairs_task${TASK_ID}.txt"

cut -f1 "$ALL_FILES_LIST" | sort | uniq -d > "$DUP_SAMPLE_IDS"
awk -F'\t' '{k=$2"\t"$3; c[k]++} END {for (k in c) if (c[k] > 1) print k}' "$ALL_FILES_LIST" | sort > "$DUP_READ_PAIRS"

echo -e "sample_id\tstatus\tgzip_r1\tgzip_r2\tfastq_blocks_r1\tfastq_blocks_r2\tread_count_r1\tread_count_r2\tpair_sync\tdetail\tr1_path\tr2_path" > "$REPORT_TSV"

add_reason() {
    local current="$1"
    local new_reason="$2"
    if [ -z "$current" ]; then
        echo "$new_reason"
    else
        echo "${current};${new_reason}"
    fi
}

check_header_format() {
    local fq="$1"
    zcat "$fq" | awk '
        NR%4==1 {
            if ($0 !~ /^@/) {
                print "header_missing_at_line_" NR
                exit 1
            }
            split(substr($0,2), a, " ")
            id=a[1]
            if (id !~ /^[A-Za-z0-9_.:-]+(\/[12])?$/) {
                print "abnormal_header_id_at_line_" NR " id=" id
                exit 1
            }
        }
        END {
            print "ok"
            exit 0
        }'
}

count_reads_and_blocks() {
    local fq="$1"
    zcat "$fq" | awk '
        END {
            if (NR % 4 != 0) {
                print "INVALID_BLOCKS\tNA"
                exit 1
            }
            print "OK\t" NR/4
            exit 0
        }'
}

check_pair_sync() {
    local r1="$1"
    local r2="$2"
    paste \
        <(zcat "$r1" | awk 'NR%4==1 {h=$0; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/\/1$/,"",h); print h}') \
        <(zcat "$r2" | awk 'NR%4==1 {h=$0; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/\/2$/,"",h); print h}') \
    | awk '
        BEGIN {line=0}
        {
            line++
            if ($1 != $2) {
                print "name_mismatch_at_pair_" line " r1=" $1 " r2=" $2
                exit 1
            }
        }
        END {
            print "ok"
            exit 0
        }'
}

while IFS=$'\t' read -r sample_id raw_r1 raw_r2; do
    [ -z "$sample_id" ] && continue

    r1="${TRIMMED_DIR}/${sample_id}_R1.trimmed.fastq.gz"
    r2="${TRIMMED_DIR}/${sample_id}_R2.trimmed.fastq.gz"

    status="PASS"
    reasons=""
    gzip_r1="OK"
    gzip_r2="OK"
    blocks_r1="OK"
    blocks_r2="OK"
    count_r1="NA"
    count_r2="NA"
    pair_sync="SKIP"

    # Metadata coherence checks
    if [ -z "$raw_r1" ] || [ -z "$raw_r2" ]; then
        reasons="$(add_reason "$reasons" "metadata_missing_raw_paths")"
    fi

    if grep -Fxq "$sample_id" "$DUP_SAMPLE_IDS"; then
        reasons="$(add_reason "$reasons" "metadata_duplicate_sample_id")"
    fi

    if grep -Fqx "$raw_r1	$raw_r2" "$DUP_READ_PAIRS"; then
        reasons="$(add_reason "$reasons" "metadata_duplicate_raw_read_pair")"
    fi

    # File presence
    if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then
        reasons="$(add_reason "$reasons" "missing_trimmed_fastq")"
        status="FAIL"
        echo -e "${sample_id}\t${status}\tNA\tNA\tNA\tNA\tNA\tNA\tSKIP\t${reasons}\t${r1}\t${r2}" >> "$REPORT_TSV"
        continue
    fi

    # Gzip integrity
    if ! gzip -t "$r1" 2>/dev/null; then
        gzip_r1="FAIL"
        reasons="$(add_reason "$reasons" "gzip_corruption_r1")"
    fi
    if ! gzip -t "$r2" 2>/dev/null; then
        gzip_r2="FAIL"
        reasons="$(add_reason "$reasons" "gzip_corruption_r2")"
    fi

    # Header format check (only if gzip OK)
    if [ "$gzip_r1" = "OK" ]; then
        hdr1_msg="$(check_header_format "$r1" 2>/dev/null)"
        if [ "$hdr1_msg" != "ok" ]; then
            reasons="$(add_reason "$reasons" "abnormal_header_r1:${hdr1_msg}")"
        fi
    fi
    if [ "$gzip_r2" = "OK" ]; then
        hdr2_msg="$(check_header_format "$r2" 2>/dev/null)"
        if [ "$hdr2_msg" != "ok" ]; then
            reasons="$(add_reason "$reasons" "abnormal_header_r2:${hdr2_msg}")"
        fi
    fi

    # FASTQ blocks + read counts (only if gzip OK)
    if [ "$gzip_r1" = "OK" ]; then
        rc1="$(count_reads_and_blocks "$r1" 2>/dev/null)"
        blocks_r1="$(echo "$rc1" | cut -f1)"
        count_r1="$(echo "$rc1" | cut -f2)"
        if [ "$blocks_r1" != "OK" ]; then
            reasons="$(add_reason "$reasons" "invalid_fastq_blocks_r1")"
        fi
    fi
    if [ "$gzip_r2" = "OK" ]; then
        rc2="$(count_reads_and_blocks "$r2" 2>/dev/null)"
        blocks_r2="$(echo "$rc2" | cut -f1)"
        count_r2="$(echo "$rc2" | cut -f2)"
        if [ "$blocks_r2" != "OK" ]; then
            reasons="$(add_reason "$reasons" "invalid_fastq_blocks_r2")"
        fi
    fi

    # Same read count check
    if [ "$count_r1" != "NA" ] && [ "$count_r2" != "NA" ] && [ "$count_r1" != "$count_r2" ]; then
        reasons="$(add_reason "$reasons" "read_count_mismatch_r1_r2")"
    fi

    # Full R1/R2 name synchronization
    if [ "$gzip_r1" = "OK" ] && [ "$gzip_r2" = "OK" ] && [ "$blocks_r1" = "OK" ] && [ "$blocks_r2" = "OK" ] && [ "$count_r1" = "$count_r2" ]; then
        pair_msg="$(check_pair_sync "$r1" "$r2" 2>/dev/null)"
        if [ "$pair_msg" = "ok" ]; then
            pair_sync="OK"
        else
            pair_sync="FAIL"
            reasons="$(add_reason "$reasons" "pair_sync_fail:${pair_msg}")"
        fi
    fi

    if [ -n "$reasons" ]; then
        status="FAIL"
    fi

    if [ -z "$reasons" ]; then
        reasons="ok"
    fi

    echo -e "${sample_id}\t${status}\t${gzip_r1}\t${gzip_r2}\t${blocks_r1}\t${blocks_r2}\t${count_r1}\t${count_r2}\t${pair_sync}\t${reasons}\t${r1}\t${r2}" >> "$REPORT_TSV"
done < "$TASK_FILES_LIST"

total_checked=$(awk 'NR>1 {n++} END {print n+0}' "$REPORT_TSV")
total_pass=$(awk -F'\t' 'NR>1 && $2=="PASS" {n++} END {print n+0}' "$REPORT_TSV")
total_fail=$(awk -F'\t' 'NR>1 && $2=="FAIL" {n++} END {print n+0}' "$REPORT_TSV")

{
    echo "Trimmed FASTQ integrity check"
    echo "Task ID           : $TASK_ID"
    echo "Sample lines      : $START_LINE-$END_LINE"
    echo "Checked samples   : $total_checked"
    echo "PASS              : $total_pass"
    echo "FAIL              : $total_fail"
    echo "Report            : $REPORT_TSV"
    echo "Duplicate sample IDs file : $DUP_SAMPLE_IDS"
    echo "Duplicate raw pairs file  : $DUP_READ_PAIRS"
} > "$SUMMARY_TXT"

cat "$SUMMARY_TXT"

# Step 1: Verify metadata coherence (duplicates and missing raw paths)
# Step 2: Verify trimmed FASTQ presence for each sample
# Step 3: Test gzip integrity with "gzip -t"
# Step 4: Check FASTQ header format and block structure
# Step 5: Verify R1/R2 read count equality and full pair synchronization
# Step 6: Write one global TSV + summary with explicit failure reasons per sample

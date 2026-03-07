#!/bin/bash
#SBATCH -p workq
#SBATCH --time=08:00:00
#SBATCH -J cp-04-mapping_qc
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/output/output-%x-%A_%a.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/00-Scripts/00.1-batch_scripts/error/error-%x-%A_%a.err
#SBATCH --array=1-18
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

# -------------------------------------------------------------------
# Post-alignment QC (chloroplastic workflow)
# Input : *.sorted.bam + *.sorted.bam.bai from step 03.1
# Output: per-sample QC reports + per-task summary table
# -------------------------------------------------------------------

set -euo pipefail

module purge
module load bioinfo/samtools/1.21

SAMPLE_SHEET="${SAMPLE_SHEET:-/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/metadata/sample_sheet_chloroplastic_complete.csv}"
ALIGN_DIR="${ALIGN_DIR:-/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/03-alignment/workflow_full/alignments}"
QC_BASE="${QC_BASE:-/home/jbonnier/work/chapitre_5/full_workflow_chloroplastic/04-mapping_qc/workflow_full}"

PER_SAMPLE_DIR="${QC_BASE}/per_sample"
SUMMARY_DIR="${QC_BASE}/summaries"
LOG_DIR="${QC_BASE}/logs"

mkdir -p "$PER_SAMPLE_DIR" "$SUMMARY_DIR" "$LOG_DIR"

if [ ! -s "$SAMPLE_SHEET" ]; then
    echo "ERROR: sample sheet not found or empty: $SAMPLE_SHEET" >&2
    exit 1
fi

FILES_LIST="${SUMMARY_DIR}/samples.list"
tail -n +2 "$SAMPLE_SHEET" | awk -F, '{print $1}' > "$FILES_LIST"

if [ ! -s "$FILES_LIST" ]; then
    echo "ERROR: empty samples list generated from $SAMPLE_SHEET" >&2
    exit 1
fi

SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-10}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-1}"
START_LINE=$(( (TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_TASK - 1 ))

GLOBAL_SUMMARY="${SUMMARY_DIR}/mapping_qc_summary.tsv"
SUMMARY_LOCK="${SUMMARY_DIR}/mapping_qc_summary.lock"

echo "=== Step 1/4: Select samples for current array task ==="
echo "TASK_ID=${TASK_ID} START_LINE=${START_LINE} END_LINE=${END_LINE} SAMPLES_PER_TASK=${SAMPLES_PER_TASK}"

TASK_FILES_LIST="${TMPDIR:-/tmp}/cp_mapping_qc_samples_${SLURM_JOB_ID:-manual}_${TASK_ID}.list"
sed -n "${START_LINE},${END_LINE}p" "$FILES_LIST" > "$TASK_FILES_LIST"

if [ ! -s "$TASK_FILES_LIST" ]; then
    echo "WARNING: no samples selected for task ${TASK_ID} (lines ${START_LINE}-${END_LINE})."
fi

while IFS= read -r sample_id; do
    [ -z "$sample_id" ] && continue

    BAM="${ALIGN_DIR}/${sample_id}.sorted.bam"
    BAI="${BAM}.bai"
    SAMPLE_LOG="${LOG_DIR}/${sample_id}.mapping_qc.log"
    SAMPLE_TSV="${PER_SAMPLE_DIR}/${sample_id}.mapping_qc.tsv"
    FLAGSTAT_OUT="${PER_SAMPLE_DIR}/${sample_id}.flagstat.txt"
    IDXSTATS_OUT="${PER_SAMPLE_DIR}/${sample_id}.idxstats.txt"
    STATS_OUT="${PER_SAMPLE_DIR}/${sample_id}.samtools_stats.txt"
    DEPTH_OUT="${PER_SAMPLE_DIR}/${sample_id}.depth_q20.tsv"

    echo "=== Step 2/4: QC for sample ${sample_id} ===" | tee "$SAMPLE_LOG"

    if [ ! -s "$BAM" ] || [ ! -s "$BAI" ]; then
        echo "ERROR: BAM/BAI missing for ${sample_id}: $BAM $BAI" | tee -a "$SAMPLE_LOG"
        {
            flock -x 200
            if [ ! -s "$GLOBAL_SUMMARY" ]; then
                printf "sample_id\tbam_path\ttotal_reads_in_bam\tmapq20_reads\tmapq20_pct\tproper_pair_reads\tproper_pair_pct\tavg_depth_q20\tbreadth_1x_q20\tbreadth_5x_q20\tbreadth_10x_q20\tnote\n" > "$GLOBAL_SUMMARY"
            fi
            printf "%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tMISSING_BAM_OR_BAI\n" "$sample_id" "$BAM" >> "$GLOBAL_SUMMARY"
        } 200>"$SUMMARY_LOCK"
        continue
    fi

    echo "Running samtools flagstat / idxstats / stats..." | tee -a "$SAMPLE_LOG"
    samtools flagstat -@ 2 "$BAM" > "$FLAGSTAT_OUT"
    samtools idxstats "$BAM" > "$IDXSTATS_OUT"
    samtools stats -@ 2 "$BAM" > "$STATS_OUT"

    echo "Computing depth and coverage breadth (Q20)..." | tee -a "$SAMPLE_LOG"
    samtools depth -a -Q 20 "$BAM" > "$DEPTH_OUT"

    TOTAL_READS=$(samtools view -c "$BAM")
    MAPQ20_READS=$(samtools view -c -q 20 "$BAM")
    PROPER_PAIR_READS=$(samtools view -c -f 2 "$BAM")

    MAPQ20_PCT=$(awk -v a="$MAPQ20_READS" -v b="$TOTAL_READS" 'BEGIN{if(b>0) printf "%.4f", (a/b)*100; else print "0.0000"}')
    PROPER_PAIR_PCT=$(awk -v a="$PROPER_PAIR_READS" -v b="$TOTAL_READS" 'BEGIN{if(b>0) printf "%.4f", (a/b)*100; else print "0.0000"}')

    AVG_DEPTH_Q20=$(awk '{s+=$3;n++} END{if(n>0) printf "%.6f", s/n; else print "0.000000"}' "$DEPTH_OUT")
    BREADTH_1X=$(awk '{if($3>=1) c++; n++} END{if(n>0) printf "%.6f", c/n; else print "0.000000"}' "$DEPTH_OUT")
    BREADTH_5X=$(awk '{if($3>=5) c++; n++} END{if(n>0) printf "%.6f", c/n; else print "0.000000"}' "$DEPTH_OUT")
    BREADTH_10X=$(awk '{if($3>=10) c++; n++} END{if(n>0) printf "%.6f", c/n; else print "0.000000"}' "$DEPTH_OUT")

    # Important: BAM from cp-03 uses '-F 4', so unmapped reads are already removed.
    NOTE="MAPPED_ONLY_BAM"

    printf "sample_id\tbam_path\ttotal_reads_in_bam\tmapq20_reads\tmapq20_pct\tproper_pair_reads\tproper_pair_pct\tavg_depth_q20\tbreadth_1x_q20\tbreadth_5x_q20\tbreadth_10x_q20\tnote\n" > "$SAMPLE_TSV"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample_id" "$BAM" "$TOTAL_READS" "$MAPQ20_READS" "$MAPQ20_PCT" "$PROPER_PAIR_READS" "$PROPER_PAIR_PCT" \
        "$AVG_DEPTH_Q20" "$BREADTH_1X" "$BREADTH_5X" "$BREADTH_10X" "$NOTE" >> "$SAMPLE_TSV"

    {
        flock -x 200
        if [ ! -s "$GLOBAL_SUMMARY" ]; then
            printf "sample_id\tbam_path\ttotal_reads_in_bam\tmapq20_reads\tmapq20_pct\tproper_pair_reads\tproper_pair_pct\tavg_depth_q20\tbreadth_1x_q20\tbreadth_5x_q20\tbreadth_10x_q20\tnote\n" > "$GLOBAL_SUMMARY"
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$sample_id" "$BAM" "$TOTAL_READS" "$MAPQ20_READS" "$MAPQ20_PCT" "$PROPER_PAIR_READS" "$PROPER_PAIR_PCT" \
            "$AVG_DEPTH_Q20" "$BREADTH_1X" "$BREADTH_5X" "$BREADTH_10X" "$NOTE" >> "$GLOBAL_SUMMARY"
    } 200>"$SUMMARY_LOCK"

    echo "DONE sample=${sample_id}" | tee -a "$SAMPLE_LOG"
done < "$TASK_FILES_LIST"

rm -f "$TASK_FILES_LIST"

echo "=== Step 3/4: Task summary ==="
echo "Global summary updated: $GLOBAL_SUMMARY"

echo "=== Step 4/4: End of cp-04 mapping QC task ${TASK_ID} ==="

cat <<'EOF'
------------------------------------------------------------
Script description (cp-04.1):
1) Reads sample IDs from chloroplastic sample sheet.
2) For each sample BAM from cp-03:
   - samtools flagstat
   - samtools idxstats
   - samtools stats
   - samtools depth -Q 20 (per-position depth)
3) Computes key QC metrics:
   - total reads in BAM
   - MAPQ>=20 reads and percentage
   - proper-pair reads and percentage
   - average depth at Q20
   - breadth of coverage at 1x / 5x / 10x (Q20)
4) Writes:
   - per-sample reports in 04-mapping_qc/.../per_sample/
   - per-task summary table in 04-mapping_qc/.../summaries/

Note:
Current cp-03 alignment keeps only mapped reads (-F 4).
So BAM-level percentages are interpreted on mapped reads only.
------------------------------------------------------------
EOF

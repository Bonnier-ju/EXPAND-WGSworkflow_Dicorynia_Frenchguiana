#!/bin/bash
#SBATCH -p workq
#SBATCH --time=01:00:00
#SBATCH -J 04-filter_samples
#SBATCH -o /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/output/output-%x-%j.out
#SBATCH -e /home/jbonnier/work/chapitre_5/full_workflow_nuclear/00-Scripts/00.1-batch_scripts/error/error-%x-%j.out
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=julien.bonnier@inrae.fr
#SBATCH --mail-type=END,FAIL

set -u
set -o pipefail

# -------------------------------------------------------------------
# STEP 0: Define thresholds and paths
# -------------------------------------------------------------------
MIN_MAPPED_PERCENT="${MIN_MAPPED_PERCENT:-70}"
MIN_AVG_DEPTH="${MIN_AVG_DEPTH:-5}"

SAMPLE_SHEET_IN="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_test.csv"
MAPPING_SUMMARY="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/04-coverage_mapping/workflow_test/mapping_coverage_summary.tsv"

SAMPLE_SHEET_KEEP="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_test_filtered.csv"
SAMPLE_SHEET_REMOVED="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/metadata/sample_sheet_test.removed_low_quality.csv"
FILTER_REPORT="/home/jbonnier/work/chapitre_5/full_workflow_nuclear/04-coverage_mapping/workflow_test/sample_filtering_report.txt"

if [ ! -s "$SAMPLE_SHEET_IN" ]; then
    echo "ERROR: input sample sheet not found or empty: $SAMPLE_SHEET_IN" >&2
    exit 1
fi

if [ ! -s "$MAPPING_SUMMARY" ]; then
    echo "ERROR: mapping summary not found or empty: $MAPPING_SUMMARY" >&2
    exit 1
fi

# -------------------------------------------------------------------
# STEP 1: Filter samples using mapped_percent + avg_depth_q20
# -------------------------------------------------------------------
awk -F',' \
    -v OFS=',' \
    -v summary="$MAPPING_SUMMARY" \
    -v keep="$SAMPLE_SHEET_KEEP" \
    -v removed="$SAMPLE_SHEET_REMOVED" \
    -v report="$FILTER_REPORT" \
    -v min_map="$MIN_MAPPED_PERCENT" \
    -v min_depth="$MIN_AVG_DEPTH" '
BEGIN {
    while ((getline line < summary) > 0) {
        if (line ~ /^sample_id\t/) {
            continue
        }
        split(line, f, "\t")
        sid = f[1]
        total_reads[sid] = f[2]
        mapped_reads[sid] = f[3]
        mapped_percent[sid] = f[4] + 0
        avg_depth[sid] = f[5] + 0
        has_summary[sid] = 1
    }
    close(summary)
}
NR == 1 {
    header = $0
    print header > keep
    print header ",total_reads,mapped_reads,mapped_percent,avg_depth_q20,removal_reason" > removed
    next
}
{
    sid = $1
    total++

    if (!(sid in has_summary)) {
        reason = "missing_mapping_coverage_summary"
        print $0 ",NA,NA,NA,NA," reason >> removed
        removed_count++
        next
    }

    fail_map = (mapped_percent[sid] < min_map)
    fail_depth = (avg_depth[sid] < min_depth)

    if (fail_map || fail_depth) {
        reason = ""
        if (fail_map) {
            reason = reason "mapped_percent<" min_map
        }
        if (fail_depth) {
            if (reason != "") {
                reason = reason ";"
            }
            reason = reason "avg_depth_q20<" min_depth
        }
        print $0 "," total_reads[sid] "," mapped_reads[sid] "," mapped_percent[sid] "," avg_depth[sid] "," reason >> removed
        removed_count++
    } else {
        print $0 >> keep
        kept_count++
    }
}
END {
    print "Sample filtering report" > report
    print "Date: " strftime("%Y-%m-%d %H:%M:%S") >> report
    print "Input sample sheet: " ENVIRON["SAMPLE_SHEET_IN"] >> report
    print "Mapping summary: " summary >> report
    print "Thresholds: mapped_percent >= " min_map ", avg_depth_q20 >= " min_depth >> report
    print "Total samples: " total >> report
    print "Kept samples: " kept_count >> report
    print "Removed samples: " removed_count >> report
    print "Kept sample sheet: " keep >> report
    print "Removed sample sheet: " removed >> report
}
' "$SAMPLE_SHEET_IN"

# -------------------------------------------------------------------
# STEP 2: Print output locations
# -------------------------------------------------------------------
echo "Filtering finished."
echo "Kept sample sheet: $SAMPLE_SHEET_KEEP"
echo "Removed sample sheet: $SAMPLE_SHEET_REMOVED"
echo "Report: $FILTER_REPORT"

# Step 1: Read per-sample QC metrics from mapping_coverage_summary.tsv
# - Uses mapped_percent and avg_depth_q20
# - Loads metrics by sample_id
#
# Step 2: Evaluate each sample from sample_sheet_test.csv
# - Keep sample if mapped_percent >= MIN_MAPPED_PERCENT and avg_depth_q20 >= MIN_AVG_DEPTH
# - Otherwise write sample to removed table with explicit removal reason
#
# Step 3: Write three outputs
# - sample_sheet_test.filtered.csv: samples kept for downstream analysis
# - sample_sheet_test.removed_low_quality.csv: excluded samples + QC metrics + reason
# - sample_filtering_report.txt: thresholds and final counts
